"""Request handlers for the variant browser."""

import logging

import db
import handlers

from . import error
from . import lookups
from . import utils


class Autocomplete(handlers.UnsafeHandler):
    """Provide autocompletion for protein names based on current query."""

    def get(self, dataset: str, query: str, ds_version: str = None):
        """
        Provide autocompletion for protein names based on current query.

        Args:
            dataset (str): dataset short name
            query (str): query
            ds_version (str): dataset version

        """
        dataset, ds_version = utils.parse_dataset(dataset, ds_version)
        ret = {}

        results = lookups.autocomplete(dataset, query, ds_version)
        ret = {'values': sorted(list(set(results)))[:20]}

        self.finish(ret)


class Download(handlers.UnsafeHandler):
    """Download variants in CSV format."""

    def get(self, dataset: str, datatype: str, item: str,  # pylint: disable=too-many-arguments
            ds_version: str = None, filter_type: str = None):
        """
        Download variants in CSV format.

        Will filter the variants if filter_type is provided.

        Args:
            dataset (str): dataset short name
            datatype (str): type of data
            item (str): query item
            ds_version (str): dataset version
            filter_type (str): type of filter to apply

        """
        # ctrl.filterVariantsBy~ctrl.filterIncludeNonPass
        dataset, ds_version = utils.parse_dataset(dataset, ds_version)
        filename = f'{dataset}_{datatype}_{item}.csv'
        self.set_header('Content-Type', 'text/csv')
        self.set_header(f'content-Disposition',
                        f'attachment; filename={filename}')

        data = utils.get_variant_list(dataset, datatype, item, ds_version)
        # filter variants based on what is shown
        if filter_type:
            filters = filter_type.split('~')
            if filters[1] == 'false':
                data['variants'] = [variant for variant in data['variants']
                                    if variant['filter_string'] == 'PASS']
            if filters[0] == 'mislof':
                data['variants'] = [variant for variant in data['variants']
                                    if variant['major_consequence'] == 'missense'
                                    or 'LoF' in variant['flags']]
            elif 'lof' in filters[0]:
                data['variants'] = [variant for variant in data['variants']
                                    if 'LoF' in variant['flags']]

        # Write header
        self.write(','.join([h[1] for h in data['headers']]) + '\n')

        for variant in data['variants']:
            headers = [h[0] for h in data['headers']]
            self.write(','.join(map(str, [variant[h] for h in headers])) + '\n')


class GetCoverage(handlers.UnsafeHandler):
    """Retrieve coverage."""

    def get(self, dataset: str, datatype: str, item: str, ds_version: str = None):
        """
        Retrieve coverage.

        Args:
            dataset (str): dataset short name
            datatype (str): type of data
            item (str): query item
            ds_version (str): dataset version

        """
        dataset, ds_version = utils.parse_dataset(dataset, ds_version)
        try:
            ret = utils.get_coverage(dataset, datatype, item, ds_version)
        except error.NotFoundError as err:
            self.send_error(status_code=404, reason=str(err))
            return
        except (error.ParsingError, error.MalformedRequest) as err:
            self.send_error(status_code=400, reason=str(err))
            return
        self.finish(ret)


class GetCoveragePos(handlers.UnsafeHandler):
    """Retrieve coverage range."""

    def get(self, dataset: str, datatype: str, item: str, ds_version: str = None):
        """
        Retrieve coverage range.

        Args:
            dataset (str): dataset short name
            datatype (str): type of data
            item (str): query item
            ds_version (str): dataset version

        """
        dataset, ds_version = utils.parse_dataset(dataset, ds_version)
        try:
            ret = utils.get_coverage_pos(dataset, datatype, item, ds_version)
        except error.NotFoundError as err:
            self.send_error(status_code=404, reason=str(err))
            return
        except (error.ParsingError, error.MalformedRequest) as err:
            self.send_error(status_code=400, reason=str(err))
            return

        self.finish(ret)


class GetGene(handlers.UnsafeHandler):
    """Request information about a gene."""

    def get(self, dataset: str, gene: str, ds_version: str = None):
        """
        Request information about a gene.

        Args:
            dataset (str): short name of the dataset
            gene (str): the gene id
            ds_version (str): dataset version

        """
        dataset, ds_version = utils.parse_dataset(dataset, ds_version)
        gene_id = gene

        ret = {'gene': {'gene_id': gene_id}}

        # Gene
        try:
            gene = lookups.get_gene(dataset, gene_id, ds_version)
        except error.NotFoundError as err:
            self.send_error(status_code=404, reason=str(err))
            return

        ret['gene'] = gene

        # Add exons from transcript
        transcript = lookups.get_transcript(dataset, gene['canonical_transcript'], ds_version)
        ret['exons'] = []
        for exon in sorted(transcript['exons'], key=lambda k: k['start']):
            ret['exons'] += [{'start': exon['start'],
                              'stop': exon['stop'],
                              'type': exon['feature_type']}]

        # Transcripts
        transcripts_in_gene = lookups.get_transcripts_in_gene(dataset, gene_id, ds_version)
        if transcripts_in_gene:
            ret['transcripts'] = []
            for transcript in transcripts_in_gene:
                ret['transcripts'] += [{'transcript_id': transcript['transcript_id']}]

        # temporary fix for names
        gene['gene_name'] = gene['name']
        gene['full_gene_name'] = gene['full_name']

        self.finish(ret)


class GetRegion(handlers.UnsafeHandler):
    """Request information about genes in a region."""

    def get(self, dataset: str, region: str, ds_version: str = None):
        """
        Request information about genes in a region.

        Args:
            dataset (str): short name of the dataset
            region (str): the region in the format chr-startpos-endpos
            ds_version (str): dataset version

        """
        dataset, ds_version = utils.parse_dataset(dataset, ds_version)

        try:
            chrom, start, stop = utils.parse_region(region)
        except error.ParsingError as err:
            self.send_error(status_code=400, reason=str(err))
            logging.error('GetRegion: unable to parse region ({})'.format(region))
            return

        ret = {'region': {'chrom': chrom,
                          'start': start,
                          'stop':  stop}}

        if utils.is_region_too_large(start, stop):
            self.send_error(status_code=400, reason='Region too large')
            return

        genes_in_region = lookups.get_genes_in_region(dataset, chrom, start, stop, ds_version)
        if genes_in_region:
            ret['region']['genes'] = []
            for gene in genes_in_region:
                ret['region']['genes'] += [{'gene_id': gene['gene_id'],
                                            'gene_name': gene['name'],
                                            'full_gene_name': gene['full_name']}]
        self.finish(ret)


class GetTranscript(handlers.UnsafeHandler):
    """Request information about a transcript."""

    def get(self, dataset: str, transcript: str, ds_version: str = None):
        """
        Request information about a transcript.

        Args:
            dataset (str): short name of the dataset
            transcript (str): the transcript id

        """
        dataset, ds_version = utils.parse_dataset(dataset, ds_version)
        transcript_id = transcript
        ret = {'transcript': {},
               'gene': {}}

        # Add transcript information
        try:
            transcript = lookups.get_transcript(dataset, transcript_id, ds_version)
        except error.NotFoundError as err:
            self.send_error(status_code=404, reason=str(err))
            return

        ret['transcript']['id'] = transcript['transcript_id']
        ret['transcript']['number_of_CDS'] = len([t for t in transcript['exons']
                                                  if t['feature_type'] == 'CDS'])

        # Add exon information
        ret['exons'] = []
        for exon in sorted(transcript['exons'], key=lambda k: k['start']):
            ret['exons'] += [{'start': exon['start'],
                              'stop': exon['stop'],
                              'type': exon['feature_type']}]

        # Add gene information
        gene = lookups.get_gene_by_dbid(transcript['gene'])
        ret['gene']['id'] = gene['gene_id']
        ret['gene']['name'] = gene['name']
        ret['gene']['full_name'] = gene['full_name']
        ret['gene']['canonical_transcript'] = gene['canonical_transcript']

        gene_transcripts = lookups.get_transcripts_in_gene_by_dbid(transcript['gene'])
        ret['gene']['transcripts'] = [g['transcript_id'] for g in gene_transcripts]

        self.finish(ret)


class GetVariant(handlers.UnsafeHandler):
    """Request information about a gene."""

    def get(self, dataset: str, variant: str, ds_version: str = None):
        """
        Request information about a gene.

        Args:
            dataset (str): short name of the dataset
            variant (str): variant in the format chrom-pos-ref-alt

        """
        # pylint: disable=too-many-locals,too-many-branches,too-many-statements
        dataset, ds_version = utils.parse_dataset(dataset, ds_version)

        ret = {'variant': {}}
        # Variant
        split_var = variant.split('-')
        if len(split_var) != 4:
            logging.error('GetVariant: unable to parse variant ({})'.format(variant))
            self.send_error(status_code=400, reason=f'Unable to parse variant {variant}')
        try:
            split_var[1] = int(split_var[1])
        except ValueError:
            logging.error('GetVariant: position not an integer ({})'.format(variant))
            self.send_error(status_code=400,
                            reason=f'Position is not an integer in variant {variant}')
            return
        orig_variant = variant
        try:
            variant = lookups.get_variant(dataset, split_var[1], split_var[0],
                                          split_var[2], split_var[3], ds_version)
        except error.NotFoundError as err:
            logging.info('Variant not found ({})'.format(orig_variant))
            self.send_error(status_code=404, reason=str(err))
            return

        # Just get the information we need
        for item in ["variant_id", "chrom", "pos", "ref", "alt", "rsid", "allele_num",
                     "allele_freq", "allele_count", "orig_alt_alleles", "site_quality",
                     "quality_metrics", "transcripts", "genes"]:
            ret['variant'][item] = variant[item]
        ret['variant']['filter'] = variant['filter_string']

        # Variant Effect Predictor (VEP) annotations
        # https://www.ensembl.org/info/docs/tools/vep/vep_formats.html
        ret['variant']['consequences'] = []
        if 'vep_annotations' in variant:
            utils.add_consequence_to_variant(variant)
            variant['vep_annotations'] = \
                utils.remove_extraneous_vep_annotations(variant['vep_annotations'])
            # Adds major_consequence
            variant['vep_annotations'] = utils.order_vep_by_csq(variant['vep_annotations'])
            ret['variant']['annotations'] = {}
            for annotation in variant['vep_annotations']:
                annotation['HGVS'] = utils.get_proper_hgvs(annotation)

                # Add consequence type to the annotations if it doesn't exist
                consequence_type = (annotation['Consequence'].split('&')[0]
                                    .replace("_variant", "")
                                    .replace('_prime_', '\'')
                                    .replace('_', ' '))
                if consequence_type not in ret['variant']['annotations']:
                    ret['variant']['annotations'][consequence_type] = \
                        {'gene': {'name': annotation['SYMBOL'],
                                  'id': annotation['Gene']},
                         'transcripts': []}

                ret['variant']['annotations'][consequence_type]['transcripts'] += \
                    [{'id': annotation['Feature'],
                      'sift': annotation['SIFT'].rstrip("()0123456789"),
                      'polyphen': annotation['PolyPhen'].rstrip("()0123456789"),
                      'canonical': annotation['CANONICAL'],
                      'modification': annotation['HGVSp'].split(":")[1] if ':' in annotation['HGVSp'] else None}]  # pylint: disable=line-too-long

        # Dataset frequencies.
        # This is reported per variable in the database data, with dataset
        # information inside the variables, so here we reorder to make the
        # data easier to use in the template

        # get the variant for other datasets with the same reference_set
        curr_dsv = db.get_dataset_version(dataset, ds_version)
        dsvs = [db.get_dataset_version(dset.short_name) for dset in db.Dataset.select()
                if dset.short_name != dataset]
        # if the only available version is not released yet
        dsvs = list(filter(lambda dsv: dsv, dsvs))
        dsvs = [dsv for dsv in dsvs if dsv.reference_set == curr_dsv.reference_set]
        dsv_groups = [(curr_dsv, variant)]
        for dsv in dsvs:
            try:
                hit = lookups.get_variant(dsv.dataset.short_name, split_var[1], split_var[0],
                                          split_var[2], split_var[3], dsv.version)
            except error.NotFoundError:
                continue
            dsv_groups.append((dsv, hit))

        frequencies = {'headers': [['Dataset', 'pop'],
                                   ['Allele Count', 'acs'],
                                   ['Allele Number', 'ans'],
                                   ['Number of Homozygotes', 'homs'],
                                   ['Allele Frequency', 'freq']],
                       'datasets': {},
                       'total': {}}
        term_map = {'allele_num': 'ans',
                    'allele_count': 'acs',
                    'allele_freq': 'freq',
                    'hom_count': 'homs'}

        for dsv_group in dsv_groups:
            ds_name = dsv_group[0].dataset.short_name

            if ds_name not in frequencies['datasets']:
                frequencies['datasets'][ds_name] = {'pop': ds_name}
            for item in term_map:
                if term_map[item] not in frequencies['total']:
                    frequencies['total'][term_map[item]] = 0
                if dsv_group[1][item] is None:
                    frequencies['datasets'][ds_name][term_map[item]] = 0
                else:
                    frequencies['datasets'][ds_name][term_map[item]] = dsv_group[1][item]
                    frequencies['total'][term_map[item]] += dsv_group[1][item]

            if 'freq' in frequencies['total']:
                frequencies['total']['freq'] = \
                    frequencies['total']['acs']/frequencies['total']['ans']
        ret['variant']['pop_freq'] = frequencies

        self.finish(ret)


class GetVariants(handlers.UnsafeHandler):
    """Retrieve variants."""

    def get(self, dataset: str, datatype: str, item: str, ds_version: str = None):
        """
        Retrieve variants.

        Args:
            dataset (str): short name of the dataset
            datatype (str): gene, region, or transcript
            item (str): item to query

        """
        dataset, ds_version = utils.parse_dataset(dataset, ds_version)
        try:
            ret = utils.get_variant_list(dataset, datatype, item, ds_version)
        except error.NotFoundError as err:
            self.send_error(status_code=404, reason=str(err))
            return
        except (error.ParsingError, error.MalformedRequest) as err:
            self.send_error(status_code=400, reason=str(err))
            return

        # inconvenient way of doing humpBack-conversion
        headers = []
        for a, h in ret['headers']:
            n = a[0] + "".join([b[0].upper() + b[1:] for b in a.split("_")])[1:]
            headers += [[n, h]]
        ret['headers'] = headers
        self.finish(ret)


class Search(handlers.UnsafeHandler):
    """Perform a search for the wanted object."""

    def get(self, dataset: str, query: str, ds_version: str = None):
        """
        Perform a search for the wanted object.

        Args:
            dataset (str): short name of the dataset
            query (str): search query

        """
        dataset, ds_version = utils.parse_dataset(dataset, ds_version)
        ret = {"dataset": dataset, "value": None, "type": None}

        datatype, identifier = lookups.get_awesomebar_result(dataset, query, ds_version)

        if datatype == "dbsnp_variant_set":
            datatype = "dbsnp"

        ret["type"] = datatype
        ret["value"] = identifier

        self.finish(ret)
