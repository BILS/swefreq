"""
Test the browser handlers
"""

import requests
import json

BASE_URL="http://localhost:4000"

def test_get_autocomplete():
    """
    Test GetAutocomplete.get()
    """
    dataset = 'SweGen'

    query = 'PA'
    response = requests.get('{}/api/datasets/{}/browser/autocomplete/{}'.format(BASE_URL, dataset, query))
    data = json.loads(response.text)
    assert set(data["values"]) == set(["PABPC1P9", "PACSIN2", "PANX2", "PARP4P3",
                                       "PARVB", "PARVG", "PATZ1", "PAXBP1", "PAXBP1-AS1"])


def test_download():
    """
    Test Download.get()
    """
    dataset = 'SweGen'

    data_type = 'transcript'
    data_item = 'ENST00000438441'
    response = requests.get('{}/api/datasets/{}/browser/download/{}/{}'.format(BASE_URL, dataset, data_type, data_item))
    assert len(response.text.split('\n')) == 407


def test_get_coverage():
    """
    Test GetCoverage.get()
    """
    dataset = 'SweGen'

    data_type = 'transcript'
    data_item = 'ENST00000438441'
    response = requests.get('{}/api/datasets/{}/browser/coverage/{}/{}'.format(BASE_URL, dataset, data_type, data_item))
    data = json.loads(response.text)
    assert len(data['coverage']) == 144


def test_get_coverage_pos():
    """
    Test GetCoveragePos.get()
    """
    dataset = 'SweGen'
    data_type = 'region'
    data_item = '22-100001-100101'
    response = requests.get('{}/api/datasets/{}/browser/coverage_pos/{}/{}'.format(BASE_URL, dataset, data_type, data_item))
    cov_pos = json.loads(response.text)
    assert cov_pos['start'] == 100001
    assert cov_pos['stop'] == 100101
    assert cov_pos['chrom'] == '22'


def test_get_gene():
    """
    Test GetGene.get()
    """
    dataset = 'SweGen'
    gene_id = 'ENSG00000015475'
    response = requests.get('{}/api/datasets/{}/browser/gene/{}'.format(BASE_URL, dataset, gene_id))
    expected = {"name": "BID", "canonicalTranscript": "ENST00000317361", "chrom": "22", "strand": "-", "geneName": "BID"}
    gene = json.loads(response.text)

    for value in expected:
        assert gene['gene'][value] == expected[value]
    assert len(gene['exons']) == 14
    assert len(gene['transcripts']) == 10

    dataset = 'SweGen'
    gene_id = 'ENSG00000015475'
    response = requests.get('{}/api/datasets/{}/browser/gene/{}'.format(BASE_URL, dataset, gene_id))
    expected = {"name": "BID", "canonicalTranscript": "ENST00000317361", "chrom": "22", "strand": "-", "geneName": "BID"}
    gene = json.loads(response.text)
    

def test_get_region():
    """
    Test GetRegion.get()
    """
    dataset = 'SweGen'
    region_def = '22-46615715-46615880'
    response = requests.get('{}/api/datasets/{}/browser/region/{}'.format(BASE_URL, dataset, region_def))
    region = json.loads(response.text)
    assert region == {'region': {'chrom': '22', 'start': 46615715, 'stop': 46615880, 'limit': 100000}}

    region_def = '22-46A1615715-46615880'
    response = requests.get('{}/api/datasets/{}/browser/region/{}'.format(BASE_URL, dataset, region_def))
    assert response.status_code == 400

    region_def = '22-46A1615715-46615880'
    response = requests.get('{}/api/datasets/{}/browser/region/{}'.format(BASE_URL, dataset, region_def))
    assert response.status_code == 400


def test_get_transcript():
    """
    Test GetTranscript.get()
    """
    dataset = 'SweGen'
    transcript_id = 'ENST00000317361'
    response = requests.get('{}/api/datasets/{}/browser/transcript/{}'.format(BASE_URL, dataset, transcript_id))
    transcript = json.loads(response.text)

    assert transcript['gene']['id'] == 'ENSG00000015475'
    assert len(transcript['exons']) == 14


def test_get_variant():
    """
    Test GetVariant.get()
    """
    dataset = 'SweGen'
    variant_id = '22-16057464-G-A'
    response = requests.get('{}/api/datasets/{}/browser/variant/{}'.format(BASE_URL, dataset, variant_id))
    variant = json.loads(response.text)

    assert variant['variant']['variantId'] == '22-16057464-G-A'
    assert variant['variant']['genes'] == ['ENSG00000233866']
    assert variant['variant']['transcripts'] == ['ENST00000424770']

    variant_id = '21-9435852-T-C'
    version = '20161223'
    response = requests.get('{}/api/datasets/{}/browser/variant/{}'.format(BASE_URL, dataset, variant_id))
    assert response.status_code == 404
    response = requests.get('{}/api/datasets/{}/version/{}/browser/variant/{}'.format(BASE_URL, dataset, version, variant_id))
    variant = json.loads(response.text)
    assert variant['variant']['variantId'] == '21-9435852-T-C'


def test_get_variants():
    """
    Test GetVariants.get()
    """
    dataset = 'SweGen'

    data_type = 'gene'
    data_item = 'ENSG00000231565'
    response = requests.get('{}/api/datasets/{}/browser/variants/{}/{}'.format(BASE_URL, dataset, data_type, data_item))
    data = json.loads(response.text)
    assert len(data['variants']) == 405

    data_type = 'region'
    data_item = '22-46615715-46615880'
    response = requests.get('{}/api/datasets/{}/browser/variants/{}/{}'.format(BASE_URL, dataset, data_type, data_item))
    data = json.loads(response.text)
    assert len(data['variants']) == 3

    data_type = 'transcript'
    data_item = 'ENST00000438441'
    response = requests.get('{}/api/datasets/{}/browser/variants/{}/{}'.format(BASE_URL, dataset, data_type, data_item))
    data = json.loads(response.text)
    assert len(data['variants']) == 405


def test_search():
    """
    Test Search.get()
    """
    dataset = 'SweGen'

    query = 'NF1P3'
    response = requests.get('{}/api/datasets/{}/browser/search/{}'.format(BASE_URL, dataset, query))
    data = json.loads(response.text)
    assert data['type'] == 'gene'
    assert data['value'] == 'ENSG00000183249'

    query = '21-9411281-T-C'
    version = '20161223'
    response = requests.get('{}/api/datasets/{}/version/{}/browser/search/{}'.format(BASE_URL, dataset, version, query))
    data = json.loads(response.text)
    assert data['type'] == 'variant'
    assert data['value'] == '21-9411281-T-C'
