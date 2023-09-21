"""
=============
bio_helper.py
=============

Helper functions for data functions that manipulate biological sequences.

`Biopython <https://biopython.org/>`_ SeqRecord objects are used to represent sequences

"""

from io import StringIO
from typing import List, Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from df.data_transfer import ColumnData, DataType, DataFunctionRequest, string_input_field
from ruse.bio.bio_data_table_helper import genbank_base64_str_to_sequence, sequence_to_genbank_base64_str, \
    string_to_sequence


def column_to_sequences(column: ColumnData, id_column: Optional[ColumnData] = None) -> List[Optional[SeqRecord]]:
    """
    Converts a Spotfire column into a list of sequence records

    :param column:  the Spotfire column
    :param id_column:  if set, row values from this column are used to set the sequence identifier
    :return: sequence records
    """
    content_type = column.contentType
    if content_type == 'chemical/x-sequence':
        sequences = [string_to_sequence(s.strip().replace(' ', ''), index) if s else None for (index, s) in
                     enumerate(column.values)]
    elif content_type == 'chemical/x-genbank':
        sequences = [genbank_base64_str_to_sequence(s.strip().replace(' ', ''), index) if s else None for (index, s) in
                     enumerate(column.values)]
    else:
        raise ValueError(f'Unable to process content type {content_type}')
    if id_column:
        for seq, seq_id in zip(sequences, id_column.values):
            if seq and id:
                seq.id = seq_id
    return sequences


def sequences_to_column(sequences: List[Optional[SeqRecord]], column_name: str, genbank_output=True) -> ColumnData:
    """
    Converts a list of sequence records to a Spotfire column object

    :param sequences: The list of sequences
    :param column_name: The name of the new column
    :param genbank_output: if set true (the default) the returned column will contain binary encoded Genbank records (otherwise string sequences will be used)
    :return: Spotfire column
    """
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


def query_from_request(request: DataFunctionRequest, input_field_name: str = 'query') -> SeqRecord:
    """
    Converts a string input field into sequence request.  The function will attempt to parse Genbank, then fasta
    formats.
    If neither works the input field data will be used a raw sequence.

    :param request: the request
    :param input_field_name: the input field name
    :return: the extracted sequence
    """
    query = string_input_field(request, input_field_name)
    if not query:
        raise ValueError()
    sequence = None
    for fmt in ['gb', 'fasta']:
        try:
            with StringIO(query) as fh:
                sequence = SeqIO.read(fh, fmt)
            if sequence:
                sequence.seq = sequence.seq.strip().replace(' ', '')
                break
        except ValueError:
            pass
    if not sequence:
        query = ''.join(query.split()).upper().strip()
        sequence = SeqRecord(Seq(query), 'Query')

    return sequence
