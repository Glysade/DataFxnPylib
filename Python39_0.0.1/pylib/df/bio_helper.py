from typing import List

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from df.data_transfer import ColumnData, DataType
from ruse.bio.bio_data_table_helper import genbank_base64_str_to_sequence, sequence_to_genbank_base64_str, \
    string_to_sequence


def values_to_sequences(column: ColumnData) -> List[SeqRecord]:
    content_type = column.contentType
    if content_type == 'chemical/x-sequence':
        sequences = [string_to_sequence(s, index) if s else None for (index, s) in
                     enumerate(column.values)]
    elif content_type == 'chemical/x-genbank':
        sequences = [genbank_base64_str_to_sequence(s, index) if s else None for (index, s) in
                     enumerate(column.values)]
    else:
        raise ValueError(f'Unable to process content type {content_type}')
    return sequences


def sequences_to_column(sequences: List[SeqRecord], column_name: str, genbank_output=True) -> ColumnData:
    def encode_seq(seq):
        if not seq:
            return None
        if genbank_output:
            return sequence_to_genbank_base64_str(seq)
        return str(seq.seq)

    values = [encode_seq(s) for s in sequences]
    if genbank_output:
        return ColumnData(name=column_name, dataType=DataType.BINARY, contentType='chemical/x-genbank', values=values)
    else:
        return ColumnData(name=column_name, dataType=DataType.STRING, contentType='chemical/x-sequence', values=values)