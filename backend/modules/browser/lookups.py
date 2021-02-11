"""Lookup functions for the variant browser."""

import logging
import re

import db

from . import error

SEARCH_LIMIT = 10000

REGION_REGEX = re.compile(r'^\s*(\d+|X|Y|M|MT)\s*([-:]?)\s*(\d*)-?([\dACTG]*)-?([ACTG]*)')


def autocomplete(dataset: str, query: str, ds_version: str = None) -> list:
    """
    Provide autocomplete suggestions based on the query.

    Args:
        dataset (str): short name of dataset
        query (str): the query to compare to the available gene names
        ds_version (str): the dataset version

    Returns:
        list: A list of genes names whose beginning matches the query

    """
    try:
        ref_set = db.get_dataset_version(dataset, ds_version).reference_set
    except AttributeError as err:
        raise error.NotFoundError(f'Reference set not found for dataset {dataset}.') from err
    query = (db.Gene.select(db.Gene.name)
             .where(((db.Gene.name.startswith(query)) &
                     (db.Gene.reference_set == ref_set))))
    gene_names = [str(gene.name) for gene in query]
    return gene_names


def get_awesomebar_result(dataset: str, query: str, ds_version: str = None) -> tuple:
    """
    Parse the search input.

    Datatype is one of:

    * `gene`
    * `transcript`
    * `variant`
    * `dbsnp_variant_set`
    * `region`

    Identifier is one of:

    * ensembl ID for gene
    * variant ID string for variant (eg. 1-1000-A-T)
    * region ID string for region (eg. 1-1000-2000)

    Follow these steps:

    * if query is an ensembl ID, return it
    * if a gene symbol, return that gene's ensembl ID
    * if an RSID, return that variant's string

    Args:
        dataset (str): short name of dataset
        query (str): the search query
        ds_version (str): the dataset version

    Returns:
        tuple: (datatype, identifier)

    """
    # pylint: disable=too-many-return-statements,too-many-branches
    query = query.strip()

    # Parse Variant types
    try:
        variant = get_variants_by_rsid(dataset, query.lower(), ds_version=ds_version)
    except (error.NotFoundError, error.ParsingError):
        pass
    else:
        if len(variant) == 1:
            return ('variant', variant[0]['variant_id'])
        return ('dbsnp_variant_set', variant[0]['rsid'])

    # Gene
    try:
        gene = get_gene_by_name(dataset, query)
        return 'gene', gene['gene_id']
    except error.NotFoundError:
        pass

    # Capital letters for all other queries
    query = query.upper()
    try:
        gene = get_gene_by_name(dataset, query)
        return 'gene', gene['gene_id']
    except error.NotFoundError:
        pass

    # Ensembl formatted queries
    if query.startswith('ENS'):
        # Gene
        try:
            gene = get_gene(dataset, query)
            return 'gene', gene['gene_id']
        except error.NotFoundError:
            pass
        # Transcript
        try:
            transcript = get_transcript(dataset, query)
            return 'transcript', transcript['transcript_id']
        except error.NotFoundError:
            pass

    # Region and variant queries
    query = query[3:] if query.startswith('CHR') else query

    match = REGION_REGEX.match(query)
    if match:
        target = match.group(0)
        target_type = 'region'
        if match.group(2) == ":":
            target = target.replace(":", "-")

        if match.group(5) and set(match.group(4)).issubset(set("ACGT")):
            target_type = 'variant'
            try:
                get_raw_variant(dataset, int(match.group(3)), match.group(1),
                                match.group(4), match.group(5), ds_version)
            except error.NotFoundError:
                target_type = 'not_found'

        return target_type, target

    return 'not_found', query


def get_coverage_for_bases(dataset: str, chrom: str, start_pos: int,
                           end_pos: int = None, ds_version: str = None) -> list:
    """
    Get the coverage for the list of bases given by start_pos->end_pos, inclusive.

    Args:
        dataset (str): short name for the dataset
        chrom (str): chromosome
        start_pos (int): first position of interest
        end_pos (int): last position of interest; if None it will be set to start_pos
        ds_version (str): version of the dataset

    Returns:
        list: coverage dicts for the region of interest. None if failed

    """
    dataset_version = db.get_dataset_version(dataset, ds_version)
    if not dataset_version:
        raise error.NotFoundError(f'Unable to find the dataset version in the database')

    if end_pos is None:
        end_pos = start_pos
    coverage = list(db.Coverage.select()
                    .where((db.Coverage.pos >= start_pos) &
                           (db.Coverage.pos <= end_pos) &
                           (db.Coverage.chrom == chrom) &
                           (db.Coverage.dataset_version == dataset_version.id))
                    .dicts())
    return coverage


def get_coverage_for_transcript(dataset: str, chrom: str, start_pos: int,
                                end_pos: int = None, ds_version: str = None) -> list:
    """
    Get the coverage for the list of bases given by start_pos->end_pos, inclusive.

    Args:
        dataset (str): short name for the dataset
        chrom (str): chromosome
        start_pos (int): first position of interest
        end_pos (int): last position of interest; if None it will be set to start_pos
        ds_version (str): version of the dataset

    Returns:
        list: coverage dicts for the region of interest

    """
    # Is this function still relevant with postgres?
    # Only entries with reported cov are in database
    coverage_array = get_coverage_for_bases(dataset, chrom, start_pos, end_pos, ds_version)
    # only return coverages that have coverage (if that makes any sense?)
    # return coverage_array
    if not coverage_array:
        return []
    covered = [c for c in coverage_array if c['mean']]
    return covered


def get_exons_in_transcript(dataset: str, transcript_id: str, ds_version: str = None) -> list:
    """
    Retrieve exons associated with the given transcript id.

    Args:
        dataset (str): short name of the dataset
        transcript_id (str): the id of the transcript
        ds_version (str): dataset version

    Returns:
        list: dicts with values for each exon sorted by start position

    """
    try:
        ref_set = db.get_dataset_version(dataset, ds_version).reference_set
    except AttributeError:
        logging.info(f'get_exons_in_transcript({dataset}, ' +
                     f'{transcript_id}): unable to find dataset dbid')
        raise error.NotFoundError(f'Reference set not found for dataset {dataset}.')

    try:
        transcript = (db.Transcript
                      .select()
                      .join(db.Gene)
                      .where((db.Transcript.transcript_id == transcript_id) &
                             (db.Gene.reference_set == ref_set))
                      .get())
    except db.Transcript.DoesNotExist as err:
        logging.info('get_exons_in_transcript({dataset}, {transcript_id}): ' +
                     'unable to retrieve transcript')
        raise error.NotFoundError(f'Transcript {transcript_id} not found in reference.') from err
    wanted_types = ('CDS', 'UTR', 'exon')
    features = sorted(list(db.Feature.select()
                           .where((db.Feature.transcript == transcript) &
                                  (db.Feature.feature_type in wanted_types))
                           .dicts()),
                      key=lambda k: k['start'])
    if not features:
        raise error.NotFoundError(f'No features found for transcript {transcript_id} in reference.')
    return features


def get_gene(dataset: str, gene_id: str, ds_version: str = None) -> dict:
    """
    Retrieve gene by gene id.

    Args:
        dataset (str): short name of the dataset
        gene_id (str): the id of the gene
        ds_version (str): dataset version

    Returns:
        dict: values for the gene; None if not found

    """
    try:
        ref_set = db.get_dataset_version(dataset, ds_version).reference_set
    except AttributeError as err:
        raise error.NotFoundError(f'Reference set not found for dataset {dataset}.') from err

    try:
        return db.Gene.select().where((db.Gene.gene_id == gene_id) &
                                      (db.Gene.reference_set == ref_set)).dicts().get()
    except db.Gene.DoesNotExist as err:
        raise error.NotFoundError(f'Gene {gene_id} not found in reference data.') from err


def get_gene_by_dbid(gene_dbid: str) -> dict:
    """
    Retrieve gene by gene database id.

    Args:
        gene_dbid (str): the database id of the gene

    Returns:
        dict: values for the gene; empty if not found

    """
    try:
        return db.Gene.select().where(db.Gene.id == gene_dbid).dicts().get()
    except db.Gene.DoesNotExist:
        return {}
    except ValueError:
        return {}


def get_gene_by_name(dataset: str, gene_name: str, ds_version: str = None) -> dict:
    """
    Retrieve gene by gene_name.

    Args:
        dataset (str): short name of the dataset
        gene_name (str): the id of the gene
        ds_version (str): dataset version

    Returns:
        dict: values for the gene; empty if not found

    """
    try:
        ref_set = db.get_dataset_version(dataset, ds_version).reference_set
    except AttributeError as err:
        raise error.NotFoundError(f'Reference set not found for dataset {dataset}.') from err

    try:
        return (db.Gene.select()
                .where((db.Gene.reference_set == ref_set) &
                       (db.Gene.name == gene_name))
                .dicts()
                .get())
    except db.Gene.DoesNotExist:
        try:
            return (db.GeneOtherNames.select(db.Gene)
                    .join(db.Gene)
                    .where((db.GeneOtherNames.name == gene_name) &
                           (db.Gene.reference_set == ref_set))
                    .dicts()
                    .get())
        except db.GeneOtherNames.DoesNotExist as err:
            logging.info(f'get_gene_by_name({dataset}, {gene_name}): unable to retrieve gene')
            raise error.NotFoundError(f'Gene {gene_name} not found in reference data') from err


def get_genes_in_region(dataset: str, chrom: str, start_pos: int,
                        stop_pos: int, ds_version: str = None) -> dict:
    """
    Retrieve genes located within a region.

    Args:
        dataset (str): short name of the dataset
        chrom (str): chromosome name
        start_pos (int): start of region
        stop_pos (int): end of region
        ds_version (str): dataset version

    Returns:
        dict: values for the gene; empty if not found

    """
    try:
        ref_set = db.get_dataset_version(dataset, ds_version).reference_set
    except AttributeError as err:
        raise error.NotFoundError(f'Reference set not found for dataset {dataset}.') from err

    genes = db.Gene.select().where((db.Gene.reference_set == ref_set) &
                                   (db.Gene.start <= stop_pos) &
                                   (db.Gene.stop >= start_pos) &
                                   (db.Gene.chrom == chrom)).dicts()
    return genes


def get_raw_variant(dataset: str, pos: int, chrom: str, ref: str,  # pylint: disable=too-many-arguments
                    alt: str, ds_version: str = None) -> dict:
    """
    Retrieve variant by position and change.

    Args:
        dataset (str): short name of the reference set
        pos (int): position of the variant
        chrom (str): name of the chromosome
        ref (str): reference sequence
        alt (str): variant sequence
        ds_version (str): dataset version

    Returns:
        dict: values for the variant; None if not found

    """
    dataset_version = db.get_dataset_version(dataset, ds_version)
    if not dataset_version:
        raise error.NotFoundError(f'Unable to find the dataset version in the database')

    try:
        variant = (db.Variant
                   .select()
                   .where((db.Variant.pos == pos) &
                          (db.Variant.ref == ref) &
                          (db.Variant.alt == alt) &
                          (db.Variant.chrom == chrom) &
                          (db.Variant.dataset_version == dataset_version))
                   .dicts()
                   .get())
        variant['genes'] = [gene['gene_id'] for gene in
                            db.VariantGenes.select(db.Gene.gene_id)
                            .join(db.Gene)
                            .where(db.VariantGenes.variant == variant['id'])
                            .dicts()]
        variant['transcripts'] = [transcript['transcript_id'] for transcript in
                                  db.VariantTranscripts.select(db.Transcript.transcript_id)
                                  .join(db.Transcript)
                                  .where(db.VariantTranscripts.variant == variant['id'])
                                  .dicts()]
        return variant
    except db.Variant.DoesNotExist as err:
        logging.info('get_raw_variant({}, {}, {}, {}, {}, {}): unable to retrieve variant'
                     .format(dataset, pos, chrom, ref, alt, dataset_version.id))
        raise error.NotFoundError(f'Variant {chrom}-{pos}-{ref}-{alt} not found') from err


def get_transcript(dataset: str, transcript_id: str, ds_version: str = None) -> dict:
    """
    Retrieve transcript by transcript id.

    Also includes exons as ['exons']

    Args:
        dataset (str): short name of the dataset
        transcript_id (str): the id of the transcript
        ds_version (str): dataset version

    Returns:
        dict: values for the transcript, including exons; None if not found

    """
    try:
        ref_set = db.get_dataset_version(dataset, ds_version).reference_set
    except AttributeError as err:
        raise error.NotFoundError(f'Reference set not found for dataset {dataset}.') from err
    try:
        transcript = (db.Transcript
                      .select(db.Transcript, db.Gene.gene_id)
                      .join(db.Gene)
                      .where((db.Transcript.transcript_id == transcript_id) &
                             (db.Gene.reference_set == ref_set))
                      .dicts()
                      .get())
        transcript['exons'] = get_exons_in_transcript(dataset, transcript_id)
        return transcript
    except db.Transcript.DoesNotExist as err:
        logging.info(f'get_transcript({dataset}, {transcript_id}): unable to retrieve transcript')
        raise error.NotFoundError(f'Transcript {transcript_id} not found in reference') from err


def get_transcripts_in_gene(dataset: str, gene_id: str, ds_version: str = None) -> list:
    """
    Get the transcripts associated with a gene.

    Args:
        dataset (str): short name of the reference set
        gene_id (str): id of the gene
        ds_version (str): dataset version

    Returns:
        list: transcripts (dict) associated with the gene; empty if no hits

    """
    try:
        ref_set = db.get_dataset_version(dataset, ds_version).reference_set
    except AttributeError as err:
        logging.info(f'get_transcripts_in_gene({dataset}, {gene_id}): unable to get ref dbid')
        raise error.NotFoundError(f'Reference set not found for dataset {dataset}.') from err

    try:
        gene = (db.Gene.select()
                .where((db.Gene.reference_set == ref_set) &
                       (db.Gene.gene_id == gene_id))
                .dicts()
                .get())
    except db.Gene.DoesNotExist as err:
        logging.info('get_transcripts_in_gene({dataset}, {gene_id}): unable to retrieve gene')
        raise error.NotFoundError(f'Gene {gene_id} not found in reference data') from err

    return list(db.Transcript.select()
                .where(db.Transcript.gene == gene['id'])
                .dicts())


def get_transcripts_in_gene_by_dbid(gene_dbid: int) -> list:
    """
    Get the transcripts associated with a gene.

    Args:
        gene_dbid (int): database id of the gene

    Returns:
        list: transcripts (dict) associated with the gene; empty if no hits

    """
    return list(db.Transcript.select()
                .where(db.Transcript.gene == gene_dbid)
                .dicts())


def get_variant(dataset: str, pos: int, chrom: str, ref: str,  # pylint: disable=too-many-arguments
                alt: str, ds_version: str = None) -> dict:
    """
    Retrieve variant by position and change.

    Args:
        dataset (str): short name of the dataset
        pos (int): position of the variant
        chrom (str): name of the chromosome
        ref (str): reference sequence
        alt (str): variant sequence
        ds_version (str): version of the dataset

    Returns:
        dict: values for the variant; None if not found

    """
    variant = get_raw_variant(dataset, pos, chrom, ref, alt, ds_version)
    variant = get_raw_variant(dataset, pos, chrom, ref, alt, ds_version)
    if variant and variant.get('rsid') and not str(variant['rsid']).startswith('rs'):
        variant['rsid'] = 'rs{}'.format(variant['rsid'])
    return variant


def get_variants_by_rsid(dataset: str, rsid: str, ds_version: str = None) -> list:
    """
    Retrieve variants by their associated rsid.

    Args:
        dataset (str): short name of dataset
        rsid (str): rsid of the variant (starting with rs)
        ds_version (str): version of the dataset

    Returns:
        list: variants as dict; no hits returns None

    """
    dataset_version = db.get_dataset_version(dataset, ds_version)
    if not dataset_version:
        raise error.NotFoundError(f'Unable to find the dataset version in the database')

    if not rsid.startswith('rs'):
        logging.warning(f'get_variants_by_rsid({dataset}, {rsid}): rsid not starting with rs')
        raise error.ParsingError('rsid not starting with rs')

    try:
        int_rsid = int(rsid.lstrip('rs'))
    except ValueError as err:
        logging.warning(f'get_variants_by_rsid({dataset}, {rsid}): not an integer after rs')
        raise error.ParsingError('Not an integer after rs') from err

    variants = (db.Variant
                .select()
                .where((db.Variant.rsid == int_rsid) &
                       (db.Variant.dataset_version == dataset_version))
                .dicts())

    if not variants:
        raise error.NotFoundError('No variants found for rsid {rsid}')
    return variants


def get_variants_in_gene(dataset: str, gene_id: str, ds_version: str = None) -> list:
    """
    Retrieve variants present inside a gene.

    Args:
        dataset (str): short name of the dataset
        gene_id (str): id of the gene
        ds_version (str): version of the dataset

    Returns:
        list: values for the variants

    """
    dataset_version = db.get_dataset_version(dataset, ds_version)
    if not dataset_version:
        raise error.NotFoundError(f'Unable to find the dataset version in the database')
    gene = get_gene(dataset, gene_id, ds_version)
    if not gene:
        raise error.NotFoundError(f'Gene {gene_id} not found in reference data')

    variants = list(db.Variant.select()
                    .join(db.VariantGenes)
                    .where((db.VariantGenes.gene == gene['id']) &
                           (db.Variant.dataset_version == dataset_version))
                    .dicts())
    for variant in variants:
        if not variant['hom_count']:
            variant['hom_count'] = 0
        variant['filter'] = variant['filter_string']

    for variant in variants:
        if variant['rsid']:
            variant['rsid'] = 'rs{}'.format(variant['rsid'])
    return variants


def get_variants_in_region(dataset: str, chrom: str, start_pos: int,
                           end_pos: int, ds_version: str = None) -> list:
    """
    Variants that overlap a region.

    Args:
        dataset (str): short name of the dataset
        chrom (str): name of the chromosom
        start_pos (int): start of the region
        end_pos (int): start of the region
        ds_version (str): version of the dataset

    Returns:
        list: variant dicts, None if no hits

    """
    dataset_version = db.get_dataset_version(dataset, ds_version)
    if not dataset_version:
        raise error.NotFoundError(f'Unable to find the dataset version in the database')
    query = (db.Variant
             .select()
             .where((db.Variant.pos >= start_pos) &
                    (db.Variant.pos <= end_pos) &
                    (db.Variant.chrom == chrom) &
                    (db.Variant.dataset_version == dataset_version))
             .dicts())

    variants = list(query)

    for variant in variants:
        if not variant['hom_count']:
            variant['hom_count'] = 0
        variant['filter'] = variant['filter_string']

    for variant in variants:
        if variant['rsid']:
            variant['rsid'] = 'rs{}'.format(variant['rsid'])
    return variants


def get_variants_in_transcript(dataset: str, transcript_id: str, ds_version: str = None) -> list:
    """
    Retrieve variants inside a transcript.

    Args:
        dataset (str): short name of the dataset
        transcript_id (str): id of the transcript (ENST)
        ds_version (str): version of the dataset

    Returns:
        list: values for the variant; None if not found

    """
    dataset_version = db.get_dataset_version(dataset, ds_version)
    if not dataset_version:
        raise error.NotFoundError(f'Unable to find the dataset version in the database')

    transcript = get_transcript(dataset, transcript_id, ds_version)
    if not transcript:
        raise error.NotFoundError(f'Transcript {transcript_id} not found in reference data')

    variants = list(db.Variant.select()
                    .join(db.VariantTranscripts)
                    .where((db.VariantTranscripts.transcript == transcript['id']) &
                           (db.Variant.dataset_version == dataset_version))
                    .dicts())

    for variant in variants:
        if not variant['hom_count']:
            variant['hom_count'] = 0
        variant['filter'] = variant['filter_string']

    for variant in variants:
        variant['vep_annotations'] = [anno for anno in variant['vep_annotations']
                                      if anno['Feature'] == transcript_id]
        if variant['rsid']:
            variant['rsid'] = 'rs{}'.format(variant['rsid'])
    return variants
