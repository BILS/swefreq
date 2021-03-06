"""
Tests for utils.py
"""

import pytest

from .. import error
from .. import lookups
from .. import utils


def test_add_consequence_to_variants():
    """
    Test add_consequence_to_variants()
    """
    variants = []
    variants.append(lookups.get_variant('SweGen', 16252949, '22', 'T', 'G'))
    variants.append(lookups.get_variant('SweGen', 16269934, '22', 'A', 'G'))

    utils.add_consequence_to_variants(variants)
    assert variants[0]['major_consequence'] == 'downstream_gene_variant'
    assert variants[1]['major_consequence'] == 'missense_variant'


def test_add_consequence_to_variant():
    """
    Test add_consequence_to_variant()
    """
    variant = lookups.get_variant('SweGen', 16252949, '22', 'T', 'G')
    utils.add_consequence_to_variant(variant)
    assert variant['major_consequence'] == 'downstream_gene_variant'

    variant['vep_annotations'][0]['Consequence'] = "stop_gained"
    utils.add_consequence_to_variant(variant)
    assert variant['category'] == 'lof_variant'
    assert variant['major_consequence'] == 'stop_gained'

    variant = lookups.get_variant('SweGen', 16269985, '22', 'C', 'G')
    utils.add_consequence_to_variant(variant)
    assert variant['category'] == 'other_variant'
    assert variant['major_consequence'] == 'intron_variant'

    variant = lookups.get_variant('SweGen', 16277852, '22', 'C', 'T')
    utils.add_consequence_to_variant(variant)
    assert variant['major_consequence'] == 'synonymous_variant'
    assert variant['category'] == 'synonymous_variant'

    variant = lookups.get_variant('SweGen', 16269934, '22', 'A', 'G')
    utils.add_consequence_to_variant(variant)
    assert variant['category'] == 'missense_variant'
    assert variant['major_consequence'] == 'missense_variant'

    variant['vep_annotations'] = []
    utils.add_consequence_to_variant(variant)
    assert variant['major_consequence'] == ''

    # bad variant
    variant = {}
    utils.add_consequence_to_variant(variant)
    assert not variant


def test_annotation_severity():
    """
    Test annotation_severity()
    """
    variant = lookups.get_variant('SweGen', 16269934, '22', 'A', 'G')
    res = utils.annotation_severity(variant['vep_annotations'][0])
    assert res == -26.9


def test_get_coverage():
    """
    Test get_coverage()
    """
    res = utils.get_coverage('SweGen', 'gene', 'ENSG00000231565')
    assert len(res['coverage']) == 144
    res = utils.get_coverage('SweGen', 'region', '22-46615715-46615880')
    assert len(res['coverage']) == 17
    res = utils.get_coverage('SweGen', 'transcript', 'ENST00000438441')
    assert len(res['coverage']) == 144

    # bad regions
    with pytest.raises(error.ParsingError):
        res = utils.get_coverage('SweGen', 'region', '22-46615715asd-46615880')
    # is seen as 22:46615715-46615880-46615880
    with pytest.raises(error.NotFoundError):
        utils.get_coverage('SweGen', 'region', '22:46615715-46615880')

    # no coverage found
    with pytest.raises(error.NotFoundError):
        res = utils.get_coverage('BAD_SET', 'transcript', 'ENST00000438441')['coverage']

    with pytest.raises(error.MalformedRequest):
        res = utils.get_coverage('SweGen', 'region', '22-1-1000000')


def test_get_coverage_pos():
    """
    Test get_coverage_pos()
    """
    res = utils.get_coverage_pos('SweGen', 'gene', 'ENSG00000231565')
    assert res['chrom'] == '22'
    assert res['start'] == 16364817
    assert res['stop'] == 16366254
    res = utils.get_coverage_pos('SweGen', 'region', '22-46615715-46615880')
    assert res['chrom'] == '22'
    assert res['start'] == 46615715
    assert res['stop'] == 46615880
    res = utils.get_coverage_pos('SweGen', 'transcript', 'ENST00000438441')
    assert res['chrom'] == '22'
    assert res['start'] == 16364817
    assert res['stop'] == 16366254

    # bad requests
    with pytest.raises(error.NotFoundError):
        utils.get_coverage_pos('BAD_SET', 'transcript', 'ENST00000438441')
    with pytest.raises(error.NotFoundError):
        utils.get_coverage_pos('SweGen', 'transcript', 'ENST1234321')
    with pytest.raises(error.NotFoundError):
        utils.get_coverage_pos('SweGen', 'gene', 'ENSG1234321')
    with pytest.raises(error.ParsingError):
        utils.get_coverage_pos('BAD_SET', 'region', '1:1:1:1')

    # too large request
    with pytest.raises(error.MalformedRequest):
        utils.get_coverage_pos('SweGen', 'region', '1-1-10000000')


def test_data_structures():
    """
    Test the constants
    """
    assert len(utils.CSQ_ORDER) == len(set(utils.CSQ_ORDER))  # No duplicates
    assert all(csq == utils.REV_CSQ_ORDER_DICT[utils.CSQ_ORDER_DICT[csq]]
               for csq in utils.CSQ_ORDER)


def test_get_flags_from_variant():
    """
    Test get_flags_from_variant()
    """
    fake_variant = {'vep_annotations': [{'LoF': 'LC', 'LoF_flags': 'something'},
                                        {'LoF': '', 'LoF_flags': ''},
                                        {'LoF': 'LC', 'LoF_flags': 'something'}]}
    flags = utils.get_flags_from_variant(fake_variant)
    assert flags == ['LC LoF', 'LoF flag']

    fake_variant = {'vep_annotations': [{'LoF': 'LC', 'LoF_flags': 'something'},
                                        {'LoF': 'HC', 'LoF_flags': 'something'}]}
    flags = utils.get_flags_from_variant(fake_variant)
    assert flags == ['LoF', 'LoF flag']

    fake_variant = {'mnps': 'no idea', 'vep_annotations': []}
    flags = utils.get_flags_from_variant(fake_variant)
    assert flags == ['MNP']


def test_get_proper_hgvs():
    """
    Test get_proper_hgvs()
    """
    annotation = {'HGVSc': 'ENST00000343518.6:c.35C>T',
                  'HGVSp': 'ENSP00000340610.6:p.Ser12Phe',
                  'major_consequence': 'splice_donor_variant'}
    assert utils.get_proper_hgvs(annotation) == 'c.35C>T'
    annotation['major_consequence'] = 'coding_sequence_variant'
    assert utils.get_proper_hgvs(annotation) == 'p.Ser12Phe'
    assert not utils.get_proper_hgvs(dict())


def test_get_protein_hgvs():
    """
    Test get_protein_hgvs()
    """
    annotation = {'HGVSc': 'ENST00000343518.6:c.35C>T',
                  'HGVSp': 'ENSP00000340610.6:p.Ser12Phe'}
    result = utils.get_protein_hgvs(annotation)
    assert result == 'p.Ser12Phe'
    annotation = {'HGVSc': 'ENST00000343518.6:c.27G>A',
                  'HGVSp': 'ENST00000343518.6:c.27G>A(p.%3D)',
                  'Protein_position': '9',
                  'Amino_acids': 'P'}
    result = utils.get_protein_hgvs(annotation)
    assert result == 'p.Pro9Pro'
    annotation['Amino_acids'] = 'Z'
    assert not utils.get_protein_hgvs(annotation)
    assert not utils.get_protein_hgvs(dict())


def test_get_transcript_hgvs():
    """
    Test get_transcript_hgvs()

    """
    annotation = {'HGVSc': 'ENST00000343518.6:c.35C>T',
                  'HGVSp': 'ENSP00000340610.6:p.Ser12Phe'}
    assert utils.get_transcript_hgvs(annotation) == 'c.35C>T'
    assert not utils.get_transcript_hgvs(dict())


def test_get_variant_list():
    """
    Test get_variant_list()
    """
    res = utils.get_variant_list('SweGen', 'gene', 'ENSG00000231565')
    assert len(res['variants']) == 178
    res = utils.get_variant_list('SweGen', 'region', '22-16360000-16361200')
    assert len(res['variants']) == 13
    res = utils.get_variant_list('SweGen', 'transcript', 'ENST00000438441')
    assert len(res['variants']) == 178
    res = utils.get_variant_list('SweGen', 'region', '22-16272587')
    assert len(res['variants']) == 4

    # bad requests
    with pytest.raises(error.NotFoundError):
        utils.get_variant_list('SweGen', 'transcript', 'ENSTWEIRD')
    with pytest.raises(error.NotFoundError):
        utils.get_variant_list('Bad_dataset', 'transcript', 'ENSTWEIRD')
    with pytest.raises(error.NotFoundError):
        utils.get_variant_list('SweGen', 'gene', 'ENSG1234321')
    with pytest.raises(error.ParsingError):
        utils.get_variant_list('SweGen', 'region', '1-1-1-1-1')

    # too large region
    with pytest.raises(error.MalformedRequest):
        utils.get_variant_list('SweGen', 'region', '22-1-1000000')


def test_order_vep_by_csq():
    """
    Test order_vep_by_csq()
    """
    annotation = [{'Consequence': 'frameshift_variant'},
                  {'Consequence': 'transcript_ablation'},
                  {'Consequence': 'mature_miRNA_variant'}]
    expected = [{'Consequence': 'transcript_ablation',
                 'major_consequence': 'transcript_ablation'},
                {'Consequence': 'frameshift_variant',
                 'major_consequence': 'frameshift_variant'},
                {'Consequence': 'mature_miRNA_variant',
                 'major_consequence': 'mature_miRNA_variant'}]
    result = utils.order_vep_by_csq(annotation)
    assert result == expected
    assert utils.order_vep_by_csq([dict()]) == [{'major_consequence': ''}]


def test_parse_dataset():
    assert utils.parse_dataset('SweGen') == ('SweGen', None)
    assert utils.parse_dataset('SweGen', '180101') == ('SweGen', '180101')
    assert utils.parse_dataset('hg19:SweGen:180101') == ('SweGen', '180101')


def test_parse_region():
    assert utils.parse_region('1-2-3') == ('1', 2, 3)
    assert utils.parse_region('X-15-30') == ('X', 15, 30)
    assert utils.parse_region('1-2') == ('1', 2, 2)

    # bad regions
    with pytest.raises(error.ParsingError):
        print(utils.parse_region('1:2:2'))
    with pytest.raises(error.ParsingError):
        utils.parse_region('1-2-2-2')
    with pytest.raises(error.ParsingError):
        utils.parse_region('asdfgh')
    with pytest.raises(error.ParsingError):
        utils.parse_region('X-15-z')
    with pytest.raises(error.ParsingError):
        utils.parse_region('X-y-15')


def test_remove_extraneous_vep_annotations():
    """
    Test remove_extraneous_vep_annotations()
    """
    annotation = [{'Consequence': 'frameshift_variant'},
                  {'Consequence': 'feature_elongation&TF_binding_site_variant'}]
    assert utils.remove_extraneous_vep_annotations(annotation) == \
        [{'Consequence': 'frameshift_variant'}]


def test_worst_csq_from_csq():
    """
    Test worst_csq_from_csq()
    """
    variant = lookups.get_variant('SweGen', 16269941, '22', 'G', 'C')
    res = utils.worst_csq_from_csq(variant['vep_annotations'][0]['Consequence'])
    assert res == 'splice_region_variant'
    res = utils.worst_csq_from_csq('non_coding_exon_variant&nc_transcript_variant')
    assert res == 'non_coding_exon_variant'


def test_worst_csq_from_list():
    """
    Test worst_csq_from_list()
    """
    csqs = ['frameshift_variant', 'missense_variant']
    assert utils.worst_csq_from_list(csqs) == 'frameshift_variant'


def test_worst_csq_index():
    """
    Test worst_csq_index()
    """
    csqs = ['frameshift_variant', 'missense_variant']
    assert utils.worst_csq_index(csqs) == 4


def test_worst_csq_with_vep():
    """
    Test worst_csq_from_vep()
    """
    veps = [{'SYMBOL': '1', 'Consequence': 'intergenic_variant', 'CANONICAL': ''},
            {'SYMBOL': '2', 'Consequence': 'frameshift_variant', 'CANONICAL': ''},
            {'SYMBOL': '3', 'Consequence': 'intron_variant', 'CANONICAL': ''},
            {'SYMBOL': '4', 'Consequence': 'stop_lost', 'CANONICAL': ''}]
    res = utils.worst_csq_with_vep(veps)
    assert res == {'SYMBOL': '2', 'Consequence': 'frameshift_variant',
                   'CANONICAL': '', 'major_consequence': 'frameshift_variant'}

    veps = [{'SYMBOL': '1', 'Consequence': 'frameshift_variant', 'CANONICAL': 'YES'},
            {'SYMBOL': '2', 'Consequence': 'frameshift_variant', 'CANONICAL': ''},
            {'SYMBOL': '3', 'Consequence': 'intron_variant', 'CANONICAL': ''},
            {'SYMBOL': '4', 'Consequence': 'stop_lost', 'CANONICAL': ''}]
    res = utils.worst_csq_with_vep(veps)
    assert res == {'SYMBOL': '1', 'Consequence': 'frameshift_variant',
                   'CANONICAL': 'YES', 'major_consequence': 'frameshift_variant'}
    assert not utils.worst_csq_with_vep([])
