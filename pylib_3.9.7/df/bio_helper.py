"""
=============
bio_helper.py
=============

Helper functions for data functions that manipulate biological sequences.

`Biopython <https://biopython.org/>`_ SeqRecord objects are used to represent sequences

"""
import base64
import gzip
from io import StringIO
import traceback
from typing import List, Optional, Tuple

from Bio import SeqIO
from Bio.PDB import PDBParser
from Bio.PDB.Structure import Structure
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from df.data_transfer import ColumnData, DataType, DataFunctionRequest, string_input_field, \
                             Notification, NotificationLevel
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
        # "".join(s.split()) removes all whitespace from a string
        sequences = [string_to_sequence("".join(s.split()), index) if s else None for (index, s) in
                     enumerate(column.values)]
    elif content_type == 'chemical/x-genbank':
        sequences = [genbank_base64_str_to_sequence("".join(s.split()), index) if s else None for (index, s) in
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

def column_to_structures(column: ColumnData, id_column: Optional[ColumnData] = None)  \
        -> Tuple[List[Optional[Structure]], List[Optional[Notification]]]:
    """
    Converts a Spotfire column into a list of PDB files

    :param column:  the Spotfire column
    :param id_column:  if set, row values from this column are used to set the structure identifier

    :return: Tuple of List of Structure objects and List of Notification objects
    """

    content_type = column.contentType
    if content_type != 'chemical/x-pdb':
        raise ValueError(f'Unable to process content type {content_type} as Structure Column.')

    parser = PDBParser(QUIET=True)
    structures = []
    notifications = []

    try:
        for index, data in enumerate(column.values):
            if id_column:
                identifier = id_column.values[index]
            else:
                identifier = f'Row {index}'

            pdb_zip = base64.b64decode(data)
            pdb_data = gzip.decompress(pdb_zip).decode()

            with StringIO(pdb_data) as structure_fh:
                structure = parser.get_structure(identifier, structure_fh)
                structure.id = identifier

            structures.append(structure)
    except Exception as ex:
        structures.append(None)
        notifications.append(Notification(level = NotificationLevel.WARNING,
                                          title = 'column_to_structures',
                                          summary = f'Error for structure {identifier}./n{ex.__class__} - {ex}',
                                          details = f'An error occurred during parsing of the compressed structure.\n' +
                                                    f'{traceback.format_exc()}'))

    return structures, notifications

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

def generate_color_gradient(start_color, end_color, num_steps):
    """
    Extract RGB components from hex values

    :param start_color: hex color for start of gradient
    :param end_color: hex color for end of gradient
    :param num_steps: number of color steps in gradient
    :return: a list of hex colors constituting the gradient
    """

    start_rgb = [int(start_color[i:i+2], 16) / 255.0 for i in (1, 3, 5)]
    end_rgb = [int(end_color[i:i+2], 16) / 255.0 for i in (1, 3, 5)]

    # Calculate step size for each RGB component
    r_step = (end_rgb[0] - start_rgb[0]) / (num_steps - 1)
    g_step = (end_rgb[1] - start_rgb[1]) / (num_steps - 1)
    b_step = (end_rgb[2] - start_rgb[2]) / (num_steps - 1)

    # Generate a linear gradient in RGB space
    color_gradient = [
        '#{0:02x}{1:02x}{2:02x}'.format(
            int((start_rgb[0] + i * r_step) * 255),
            int((start_rgb[1] + i * g_step) * 255),
            int((start_rgb[2] + i * b_step) * 255)
        ) for i in range(num_steps)
    ]

    return color_gradient