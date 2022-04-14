import os
import uuid
from enum import Enum
from typing import NamedTuple, List, Dict

from Bio import AlignIO
from Bio import SeqIO
from Bio.Emboss.Applications import NeedleCommandline, WaterCommandline, StretcherCommandline
from Bio.SeqRecord import SeqRecord

from ruse.bio.bio_data_table_helper import add_sequences_to_data_table, add_sequences_to_data_table_as_genbank_column
from ruse.bio.sequence_align import copy_features_to_gapped_sequence
from ruse.util.data_table import DataTable
from ruse.util.util import get_memory_in_gb


class PairwiseAlignmentMethod(Enum):
    NEEDLE = 1
    WATER = 2


class PairwiseAlignment(NamedTuple):
    score: float
    query: SeqRecord
    target: SeqRecord


def needle_pairwise_alignment(query: SeqRecord, targets: List[SeqRecord]) -> List[PairwiseAlignment]:
    # substitute stretcher for needle, if we don't have much memory
    if get_memory_in_gb() > 3.0:
        return _pairwise_alignment(PairwiseAlignmentMethod.NEEDLE, query, targets)
    else:
        return stretcher_pairwise_alignment(query, targets)


def water_pairwise_alignment(query: SeqRecord, targets: List[SeqRecord]) -> List[PairwiseAlignment]:
    return _pairwise_alignment(PairwiseAlignmentMethod.WATER, query, targets)


def _to_pairwise_alignment(query, align, target, local):
    aligned_sequences = [s for s in align]
    aligned_query = aligned_sequences[0]
    aligned_target = aligned_sequences[1]
    if query.features:
        copy_features_to_gapped_sequence(query, aligned_query, local)
    if target.features:
        copy_features_to_gapped_sequence(target, aligned_target, local)
    score = align.annotations['score']
    return PairwiseAlignment(score, aligned_query, aligned_target)


def _write_in_file(basename: str, seqs: List[SeqRecord], ind: str) -> str:
    seqs = [SeqRecord(s.seq, id=s.id, name='', description='') for s in seqs]
    in_file = '{}_{}in.fasta'.format(basename, ind)
    with open(in_file, 'w') as fh:
        SeqIO.write(seqs, fh, 'fasta')
    return in_file


def _cleanup(basename: str) -> None:
    for path in ["{}_ain.fasta", "{}_bin.fasta", "{}_out.align"]:
        path = path.format(basename)
        if os.path.exists(path):
            os.remove(path)


def stretcher_pairwise_alignment(query: SeqRecord, targets: List[SeqRecord]) -> List[
    PairwiseAlignment]:
    # TODO determine if sequence is protein or dna and adjust opts
    basename = str(uuid.uuid4())
    a_in_file = _write_in_file(basename, [query], 'a')

    def align_target(target):
        b_in_file = _write_in_file(basename, [target], 'b')
        out_file = '{}_out.align'.format(basename)
        # substitute stretcher for needle, due to memory requirements
        # Can only do single target sequence
        args = {'asequence': a_in_file, 'bsequence': b_in_file, 'outfile': out_file,
                'gapopen': 12, 'gapextend': 2, 'aformat': 'srspair'}
        cline = StretcherCommandline(**args)
        print("emboss command line is {}".format(cline))
        stdout, stderr = cline()
        alignments = list(AlignIO.parse(out_file, "emboss"))
        assert len(alignments) == 1
        return _to_pairwise_alignment(query, alignments[0], target, False)

    pairwise_alignments = [align_target(t) for t in targets]
    _cleanup(basename)
    return pairwise_alignments


def _pairwise_alignment(method: PairwiseAlignmentMethod, query: SeqRecord, targets: List[SeqRecord]) -> List[
    PairwiseAlignment]:
    # TODO determine if sequence is protein or dna and adjust opts

    if not targets:
        return []

    basename = str(uuid.uuid4())
    a_in_file = _write_in_file(basename, [query], 'a')
    b_in_file = _write_in_file(basename, targets, 'b')
    out_file = '{}_out.align'.format(basename)

    # default gap values from documentation- they shouldn't really need to be specified
    args = {'asequence': a_in_file, 'bsequence': b_in_file, 'outfile': out_file,
            'gapopen': 10.0, 'gapextend': 0.5}
    if method == PairwiseAlignmentMethod.NEEDLE:
        cline = NeedleCommandline(**args)
    elif method == PairwiseAlignmentMethod.WATER:
        cline = WaterCommandline(**args)
    else:
        raise ValueError()
    print("emboss command line is {}".format(cline))
    stdout, stderr = cline()

    alignments = list(AlignIO.parse(out_file, "emboss"))
    assert len(alignments) <= len(targets)
    pairwise_alignments = [_to_pairwise_alignment(query, a, t, method == PairwiseAlignmentMethod.WATER) for a, t in
                           zip(alignments, targets)]

    _cleanup(basename)

    return pairwise_alignments


def add_alignment_to_data_table(data_table: DataTable, alignments: List[PairwiseAlignment],
                                id_to_index_map: Dict[str, int],
                                genbank_targets: bool = False) -> None:
    alignment_map = {alignment.target.id: alignment for alignment in alignments}
    index_to_id_map = {v: k for k, v in id_to_index_map.items()}
    alignment_rows = [alignment_map.get(index_to_id_map[row_num]) for row_num in range(len(data_table.data))]
    aligned_sequences = ['{}|{}'.format(alignment.query.seq, alignment.target.seq) if alignment else None
                         for alignment in alignment_rows]
    scores = [alignment.score if alignment else None for alignment in alignment_rows]
    queries = [alignment.query if alignment else None for alignment in alignment_rows]
    targets = [alignment.target if alignment else None for alignment in alignment_rows]
    for target, query in zip(targets, queries):
        query.id = target.id

    data_table.add_column('Pairwise Alignment', 'string', content_type='chemical/x-sequence-pair',
                          data=aligned_sequences)
    data_table.add_column('Pairwise Alignment Score', 'float', data=scores)
    add_sequences_to_data_table(data_table, queries, [], 'Pairwise Alignment Query')
    if genbank_targets:
        add_sequences_to_data_table_as_genbank_column(data_table, targets, [], 'Pairwise Alignment Target')
    else:
        add_sequences_to_data_table(data_table, targets, [], 'Pairwise Alignment Target')
