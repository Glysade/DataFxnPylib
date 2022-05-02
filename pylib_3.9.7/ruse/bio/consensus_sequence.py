"""
Copyright (C) 2017-2022 Glysade, LLC
"""

from collections import Counter
from typing import List, Optional

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ruse.bio.bio_data_table_helper import data_table_column_to_sequence, sequence_to_genbank_base64_str
from ruse.bio.bio_util import is_genbank_column
from ruse.util.data_table import DataTable


def add_consensus_sequence_to_data_table(data_table: DataTable, sequence_column: int,
                                         min_frequency: Optional[float] = .0,
                                         id_column: Optional[int] = None) -> None:
    sequences = data_table_column_to_sequence(data_table, sequence_column)
    consensus_sequence = determine_consensus_sequence(sequences, min_frequency=min_frequency)
    column_def = data_table.columns[sequence_column]
    genbank_format = is_genbank_column(column_def)
    encoded_sequence = sequence_to_genbank_base64_str(consensus_sequence) if genbank_format \
        else str(consensus_sequence.seq)

    new_row = [None] * len(data_table.columns)
    new_row[sequence_column] = encoded_sequence
    if id_column:
        new_row[id_column] = consensus_sequence.name
    data_table.data.append(new_row)


def determine_consensus_sequence(sequences: List[SeqRecord], min_frequency: float = .0) -> SeqRecord:
    size = max(len(s) for s in sequences)
    consensus_items = [_most_common_item_at_position(sequences, p, min_frequency) for p in range(size)]

    name = "consensus"
    if min_frequency == 1.0:
        name = "conserved"
    consensus_seq = SeqRecord(Seq(''.join(consensus_items)), id=name, name=name)
    return consensus_seq


def _most_common_item_at_position(sequences: List[SeqRecord], position: int, min_frequency: float) -> str:
    def value_at_position(seq):
        if len(seq) > position:
            return seq[position]
        else:
            return '-'

    items = [value_at_position(s) for s in sequences]
    items = [i for i in items if i != '-']
    counts = Counter(items).most_common(1)
    if not counts:
        return '-'
    (item, count) = counts[0]

    n_sequences = len(sequences)
    if min_frequency == .0:
        return item
    elif min_frequency == 1.0:
        if count == n_sequences:
            return item
        else:
            return '-'
    else:
        freq = float(count) / float(n_sequences)
        if freq >= min_frequency:
            return item
        else:
            return '-'
