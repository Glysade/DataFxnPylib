import re
from typing import NamedTuple, List

from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

from ruse.bio.sequence_align import copy_features_to_gapped_sequence


class IgRegion(NamedTuple):
    name: str
    start: int
    end: int


class IgHit(NamedTuple):
    chain: str
    query: str
    target: str
    evalue: float
    bits: float
    query_start: int
    query_end: int
    target_id: str
    target_def: str = ""

    def identifier(self) -> str:
        return self.target_id


class IgResult(NamedTuple):
    regions: List[IgRegion]
    hits: List[IgHit]


def parse_igblastp_results(file: str) -> IgResult:
    return _parse_igblast_results(file, 'igblastp')


def parse_igblastn_results(file: str) -> IgResult:
    return _parse_igblast_results(file, 'igblastn')


def _parse_igblast_results(file: str, method: str) -> IgResult:
    with open(file, 'r') as fh:
        lines: List[str] = fh.read().splitlines()

    lines.reverse()
    line: str = lines.pop()
    assert method.upper() in line

    regions: List[IgRegion] = []
    while 'Alignment summary' not in line:
        if '# 0 hits found' in line:
            raise ValueError("No hits found")
        line = lines.pop()
    line = lines.pop()
    while line and not line.startswith('Total'):
        (name, start_str, end_str, l_str, _, _, _, _) = line.split('\t')
        start = int(start_str)
        end = int(end_str)
        assert end - start + 1 == int(l_str)
        regions.append(IgRegion(name, start - 1, end - 1))
        line = lines.pop()

    while '# Fields: query' not in line:
        line = lines.pop()
    line = lines.pop()
    m = re.search(r'^# (\d+) hits? found', line)
    assert m
    n_hits: int = int(m.group(1))

    line = lines.pop()
    hits: List[IgHit] = []
    while line and not line.startswith('# '):
        (chain, query_id, subject_id, _, l_str, _, _, _, query_start_str, query_end_str, subject_start_str,
         subject_end_str,
         evalue_str, bits_str, query,
         subject, _) = line.split()
        query_start = int(query_start_str) - 1
        query_end = int(query_end_str) - 1
        assert query_end - query_start + 1 == int(l_str)
        evalue = float(evalue_str)
        bits = float(bits_str)
        hit = IgHit(chain, query, subject, evalue, bits, query_start, query_end, subject_id)
        hits.append(hit)
        line = lines.pop()

    assert n_hits == len(hits)
    return IgResult(regions, hits)


def annotate_igblast_regions(hit: IgResult, query: SeqRecord, ungap_query: bool = True) -> SeqRecord:
    if ungap_query and '-' in query.seq:
        ungapped_query = SeqRecord(query.seq.ungap(), query.id, query.name, query.description)
        copy_features_to_gapped_sequence(query, ungapped_query)
        query = ungapped_query

    for region in hit.regions:
        feature = SeqFeature(FeatureLocation(region.start, region.end + 1), type='region', qualifiers={
            'note': [region.name]
        })
        query.features.append(feature)

    return query
