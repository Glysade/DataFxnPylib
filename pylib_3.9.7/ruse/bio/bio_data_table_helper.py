"""
========================
bio_data_table_helper.py
========================

Copyright (C) 2017 Anodyne Informatics, LLC

Functions for adding sequence data to and extracting sequence data from :class:`ruse.util.data_table.DataTable` objects

- Encoding and decoding data table cells containing sequence data
- Creating data tables and data table columns from Blast search and sequnce alignment results

"""

import base64
import gzip
import re
import uuid

from typing import List, Dict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import BytesIO, StringIO

from ruse.bio.bio_util import is_dna, is_defined_sequence
from ruse.bio.blast_parse import BlastResult, BlastResults, MultipleBlastResults
from ruse.util.data_table import DataTable
from ruse.util import log


def genbank_cell_to_sequence(data_table, column_no: int, row_index: int) -> SeqRecord:
    """
    Given a data table and cell decodes the cell to to a Biopython SeqRecord.

    The cell must have been created using a GenBank formatted sequence record that has been gzipped then encoded as a Base64(utf8)
    string.

    :param data_table: The data table
    :param column_no: column index of cell
    :param row_index: row index of cell
    :return: Sequence as :class:`Bio.SeqRecord.SeqRecord`
    """
    """Converts base64 gzipped genbank entries to a sequence record"""
    cell = data_table.data[row_index][column_no]
    record = genbank_base64_str_to_sequence(cell, row_index)
    return record


def string_to_sequence(data: str, row_index: int) -> SeqRecord:
    record = SeqRecord(Seq(data), id=f'row_{row_index}', name='', description='')
    record.annotations['__row_index__'] = str(row_index)
    return record


def genbank_base64_str_to_sequence(data: str, row_index: int) -> SeqRecord:
    """
    Decodes a string to a Biopython SeqRecord

    The string  must have been created using a GenBank formatted sequence record that has been gzipped then
    encoded in Base64(utf8)

    :param data: the string
    :return: Sequence as :class:`Bio.SeqRecord.SeqRecord`
    """
    binary = base64.b64decode(data)
    value = gzip.decompress(binary)
    for charset in ['ascii', 'utf-8', 'latin-1']:
        try:
            genbank_str = value.decode(charset)
            with StringIO(genbank_str) as fh:
                record = SeqIO.read(fh, 'gb')
            record.annotations['__row_index__'] = str(row_index)
            return record
        except UnicodeDecodeError:
            pass
        except ValueError as err:
            raise ValueError(f"Genbank entry appears to be malformed: {err}")

    raise UnicodeDecodeError("Unable to decode genbank cell")


def sequence_to_genbank_base64_str(record: SeqRecord) -> str:
    """
    Converts a Biopython sequence record to an encoded string.  The sequenced is converted to a Genbank format string,
    then gzipped and encoded in Base64(utf8)

    :param record: Sequence as :class:`Bio.SeqRecord.SeqRecord`
    :return: Encoded string
    """

    seq = record.seq
    # these 2 checks are to make sure we can write records read from fasta into genbank format
    if 'molecule_type' not in record.annotations:
        if is_dna(str(seq)):
            record.annotations['molecule_type'] = 'DNA'
        else:
            record.annotations['molecule_type'] = 'protein'
    if len(record.name) > 16:
        record.name = record.name[:16]
    genbank_str = record.format('gb')
    data = gzip.compress(genbank_str.encode('utf8'))
    return base64.b64encode(data).decode('utf-8')


def ids_to_row_map(sequences: List[SeqRecord], missing: List[int]) -> Dict[str, int]:
    all_sequences = sequences.copy()
    for index in missing:
        all_sequences.insert(index, None)

    return {s.id: i for i, s in enumerate(all_sequences) if s}


def uniq_sequence_ids(sequences: List[SeqRecord]) -> None:
    id_counts = {}
    for record in sequences:
        id = record.id
        if id in id_counts:
            count = id_counts[id]
            count += 1
            id_counts[id] = count
            record.id = f'{id}_{count}'
        else:
            id_counts[id] = 1


def data_table_column_to_sequence(data_table: DataTable, column_no: int, missing_indices: List[int] = []) -> \
        [SeqRecord]:
    """
    Decodes a data table sequence column to a list of Biopython SeqRecords.  The content-type definition for the column
    is used to determine the column encoding

    :param data_table: The data table
    :param column_no: column index to decode
    :return: Tuple containing list of SeqRecords contained in column and list of missing ids
    """
    """Convert a column containing sequence data to a list of SeqRecords"""

    column_def = data_table.columns[column_no]

    if column_def['dataType'] == 'string' and column_def['properties']['ContentType'] == 'chemical/x-sequence':
        sequences = [string_to_sequence(row[column_no], row_index) for row_index, row in enumerate(data_table.data) if
                     row[column_no]]
    elif column_def['dataType'] == 'binary' and column_def['properties']['ContentType'] == 'chemical/x-genbank':
        sequences = [genbank_cell_to_sequence(data_table, column_no, row_no) for
                     row_no, row in enumerate(data_table.data) if row[column_no]]
    else:
        raise ValueError("column {} is not a sequence column".format(column_no))

    missing_rows = [row_no for row_no, row in enumerate(data_table.data) if row[column_no] is None]
    missing_indices.clear()
    missing_indices.extend(missing_rows)
    uniq_sequence_ids(sequences)
    return sequences


def data_table_cell_to_sequence(data_table: DataTable, row_no: int, column_no: int) -> SeqRecord:
    """
    Given a data table and cell decodes the cell to to a Biopython SeqRecord.

    The content-type definition for the column is used to determine the cell encoding

    :param data_table:  The data table
    :param row_no: Row index
    :param column_no: Column index
    :return: Sequence as :class:`Bio.SeqRecord.SeqRecord`
    """
    if not data_table.data[row_no][column_no]:
        return None

    column_def = data_table.columns[column_no]
    if column_def['dataType'] == 'string' and column_def['properties']['ContentType'] == 'chemical/x-sequence':
        return string_to_sequence(data_table.data[row_no][column_no], row_no)
    if column_def['dataType'] == 'binary' and column_def['properties']['ContentType'] == 'chemical/x-genbank':
        return genbank_cell_to_sequence(data_table, column_no, row_no)

    raise ValueError("column {} is not a sequence column".format(column_no))


# TODO merge 2 add sequences to data table methods
def add_sequences_to_data_table(data_table: DataTable, sequences: List[SeqRecord], missing_indices: List[int],
                                column_name: str = 'Aligned Sequence', add_dendrogram_names: bool = False) -> None:
    """
    Add matching/aligned sequences to a new column as sequence strings.  Sequences are assigned to rows based on their id
    which should map to an existing row name. The sequence is added as a text field

    :param data_table: The data table
    :param sequences: List of new SeqRecords to add in a new column.
    :param missing_indices: Ordered list of null input sequence indices
    :param column_name: New column name
    :param add_dendrogram_names: Add a second column containing sequence names that is used for construction of a dendrogram
    """

    assert len(missing_indices) + len(sequences) == len(data_table.data)
    new_values = [None] * len(data_table.data)
    new_ids = [None] * len(data_table.data)

    # can add data by array index or sequence name
    for index, seq in enumerate(sequences):
        if not seq:
            continue

        new_values[index] = str(seq.seq)
        new_ids[index] = seq.id

    for index in missing_indices:
        new_values.insert(index, None)
        new_ids.insert(index, None)

    for row_no, row in enumerate(data_table.data):
        row.append(new_values[row_no])
        if add_dendrogram_names:
            row.append(new_ids[row_no])

    column_name = data_table.unique_column_name(column_name)
    dendrogram_column = data_table.unique_column_name('Dendrogram')
    data_table.columns.append(
        {'name': column_name, 'dataType': 'string', 'properties': {'ContentType': 'chemical/x-sequence'},
         'id': str(uuid.uuid4())})
    if add_dendrogram_names:
        data_table.columns.append({'name': dendrogram_column, 'dataType': 'string', 'properties': {}})


def add_sequences_to_data_table_as_genbank_column(data_table: DataTable, sequences: List[SeqRecord],
                                                  missing_indices: List[int],
                                                  column_name: str = 'Aligned Sequence',
                                                  add_dendrogram_names: bool = False) -> None:
    """
    Add matching/aligned sequences to a new column as sequence strings.  Sequences are assigned to rows based on their id
    which should map to an existing row name. The sequence is added as an encoded Genbank record

    :param data_table: The data table
    :param sequences: List of new SeqRecords to add in a new column.
    :param missing_indices: Ordered list of missing null sequences
    :param column_name: New column name
    :param add_dendrogram_names: Add a second column containing sequence names that is used for construction of a dendrogram
    """

    assert len(data_table.data) == len(sequences) + len(missing_indices)
    new_values = [None] * len(data_table.data)
    new_ids = [None] * len(data_table.data)

    # can index sequence through name mapping, or just array position
    for index, seq in enumerate(sequences):
        if not seq:
            continue
        new_values[index] = sequence_to_genbank_base64_str(seq)
        new_ids[index] = seq.id

    for index in missing_indices:
        new_values.insert(index, None)
        new_ids.insert(index, None)

    for row_no, row in enumerate(data_table.data):
        row.append(new_values[row_no])
        if add_dendrogram_names:
            row.append(new_ids[row_no])

    column_name = data_table.unique_column_name(column_name)
    dendrogram_column = data_table.unique_column_name('Dendrogram')
    data_table.columns.append(
        {'name': column_name, 'dataType': 'binary', 'properties': {'ContentType': 'chemical/x-genbank'},
         'id': str(uuid.uuid4())})
    if add_dendrogram_names:
        data_table.columns.append({'name': dendrogram_column, 'dataType': 'string', 'properties': {}})


def add_blast_results_to_data_table(data_table: DataTable, blast_results: BlastResult,
                                    id_to_index_map: Dict[str, int]) -> None:
    """
    Adds columns from a :class:`ruse.bio.blast_parse.BlastResult` object into a :class:`ruse.util.data_table.DataTable`

    This is used when a sequence column from a data table is used to create a blast database.  From a search against that
    database hit information is appended as columns to the data table.  Each hit is appended to the row that contains the
    target sequence.  The columns include the aligned sequence pair, blast expect value, blast score and number of bits.

    :param data_table: The data table
    :param blast_results: The blast results
    """

    # 4 values per row: alignment, expect, score, bits
    for row in data_table.data:
        row.extend([None] * 4)

    for hit in blast_results.hits:
        name = hit.data_table_name()
        if name not in id_to_index_map:
            raise ValueError('Cannot find row for sequence name {}'.format(name))
        row_no = id_to_index_map[name]
        align_str = '{}|{}'.format(hit.query, hit.target)
        data_table.data[row_no][-4:] = [align_str, hit.evalue, hit.score, hit.bits]

    data_table.add_column('Aligned Sequence', 'string', content_type='chemical/x-sequence-pair', fill_rows=False)
    data_table.add_column('Evalue', 'float', fill_rows=False)
    data_table.add_column('Score', 'integer', fill_rows=False)
    data_table.add_column('Bits', 'float', fill_rows=False)


def create_data_table_from_sequences(sequences: List[SeqRecord]) -> DataTable:
    """
    Create a :class:`ruse.util.data_table.DataTable` from a list of Biopython :class:`Bio.SeqRecord.SeqRecord`.
    The sequences are added as strings. Two columns are created: an id column and a sequence column

    :param sequences: Input sequences as a list of :class:`Bio.SeqRecord.SeqRecord`
    :return: The data table
    """

    data_table = DataTable()

    data_table.data = [[seq.id, str(seq.seq)] for seq in sequences]
    for index, seq in enumerate(sequences):
        seq.annotations['__row_index__'] = index

    data_table.columns = [{'name': 'id', 'dataType': 'string',
                           'id': str(uuid.uuid4())},
                          {'name': 'Sequence', 'dataType': 'string',
                           'properties': {'ContentType': 'chemical/x-sequence'},
                           'id': str(uuid.uuid4())}]
    return data_table


def create_data_table_from_genbank_sequences(sequences: List[SeqRecord]) -> DataTable:
    """
    Create a :class:`ruse.util.data_table.DataTable` from a list of Biopython :class:`Bio.SeqReord.SeqRecord`.
    A single column is created containing the sequences as gzipped Base64 encoded Genbank records.

    :param sequences: Input sequences as a list of :class:`Bio.SeqRecord.SeqRecord`
    :return: The :class:`ruse.util.data_table.DataTable` data table
    """

    data_table = DataTable()

    data_table.data = [[seq.id, sequence_to_genbank_base64_str(seq)] for seq in sequences]
    data_table.columns = [
        {'name': 'id', 'dataType': 'string',
         'id': str(uuid.uuid4())},
        {'name': 'Sequence', 'dataType': 'binary', 'properties': {'ContentType': 'chemical/x-genbank'},
         'id': str(uuid.uuid4())}]
    return data_table


def create_data_table_from_blast_results(query: SeqRecord, blast_results: BlastResults, database_name: str,
                                         query_column_type: str = None) -> DataTable:
    """
    Creates a new :class:`ruse.util.data_table.DataTable` from :class:`ruse.bio.blast_parse.BlastResult`

    The new table contains one row for each target sequence that matches the query:

    - the aligned sequences as a text pair (query and target matches separated by a '|')
    - the id of the target
    - the definition of the target
    - blast expect value
    - blast score
    - blast number of matching bits
    - the target sequence record as a gzipped Base64 encoded genbank record

    :param query: :class:`Bio.SeqRecord.SeqRecordSeqRecord` object containing query sequence
    :param blast_results: Blast results from a web search
    :return: The :class:`ruse.util.data_table.DataTable` data table
    """

    table = DataTable()

    include_pdb_uris = database_name == 'pdb'

    log.info("database_name = " + database_name)
    log.info("include_pdb_uris = " + str(include_pdb_uris))

    genbank_query = query_column_type == 'genbank'
    query_seq_str = sequence_to_genbank_base64_str(query) if genbank_query else str(query.seq)
    for hit in blast_results.hits:
        align_str = '{}|{}'.format(hit.query, hit.target)
        genbank_str = sequence_to_genbank_base64_str(hit.target_record) if hit.target_record else None

        if include_pdb_uris:
            parts = hit.target_id.split('|')
            if len(parts) == 5 and parts[2] == 'pdb':
                uri = 'https://files.rcsb.org/download/%s.pdb' % parts[3]
            elif len(parts) == 3 and parts[0] == 'pdb':
                uri = 'https://files.rcsb.org/download/%s.pdb' % parts[1]
            else:
                uri = None
                log.warn('Unable to extract PDB entry from Sequence id = {}'.format(hit.target_id))
            log.info('Sequence id = {}, uri = {}'.format(hit.target_id, uri))
            table.data.append(
                [align_str, hit.target_id, hit.target_def, hit.evalue, hit.score, hit.bits, query_seq_str, genbank_str,
                 uri])
        else:
            table.data.append(
                [align_str, hit.target_id, hit.target_def, hit.evalue, hit.score, hit.bits, query_seq_str, genbank_str])

    if genbank_query:
        query_seq_data_type = 'binary'
        query_seq_column_type = 'chemical/x-genbank'
    else:
        query_seq_data_type = 'string'
        query_seq_column_type = 'chemical/x-sequence'
    table.columns = [{'name': 'Aligned Sequence', 'dataType': 'string',
                      'properties': {'ContentType': 'chemical/x-sequence-pair'},
                      'id': str(uuid.uuid4())},
                     {'name': 'Target Id', 'dataType': 'string'},
                     {'name': 'Target Definition', 'dataType': 'string'},
                     {'name': 'Evalue', 'dataType': 'float'},
                     {'name': 'Score', 'dataType': 'integer'},
                     {'name': 'Bits', 'dataType': 'float'},
                     {'name': 'Query Sequence', 'dataType': query_seq_data_type,
                      'properties': {'ContentType': query_seq_column_type},
                      'id': str(uuid.uuid4())},
                     {'name': 'Target Sequence', 'dataType': 'binary',
                      'properties': {'ContentType': 'chemical/x-genbank'},
                      'id': str(uuid.uuid4())}]

    if include_pdb_uris:
        table.columns.append(
            {'name': 'URL', 'dataType': 'string', 'properties': {'ContentType': 'chemical/x-uri', "Dimension": "3"}})

    return table


def create_data_table_from_multiple_blast_searches(query_sequences: List[SeqRecord], results: MultipleBlastResults,
                                                   sequence_column_type: str = None) -> DataTable:
    """
    Creates a new :class:`DataTable` from :class:`MultipleBlastResults` obtained by performing blast alignment searched
    using the query sequences

    The new table contains one row for each target sequence that matches the query:

    - the aligned sequences as a text pair (query and target matches separated by a '|')
    - the id of the query
    - the definition of the query
    - the id of the target
    - the definition of the target
    - blast expect value
    - blast score
    - blast number of matching bits
    - the query  sequence record either as a gzipped Base64 encoded genbank record or sequence string
    - the target sequence record either as a gzipped Base64 encoded genbank record or sequence string

    Note that target results from NT databases may include full length chromosomes and Genbank entries for those are
    large, in which case the string encoding should be preferred

    :param query_sequences:  query sequences
    :param results: multiple blast search results
    :param sequence_column_type: should be either 'genbank' or 'string', sets how query and target sequences are encoded in the table
    :return: The :class:`ruse.util.data_table.DataTable` data table
    """

    table = DataTable()
    assert len(query_sequences) == len(results.query_hits)

    if sequence_column_type == 'genbank':
        genbank_sequence_columns = True
    elif sequence_column_type == 'string':
        genbank_sequence_columns = False
    else:
        # genbank_sequence_columns = True if results.number_hits() < 1000 else False
        genbank_sequence_columns = False

    for query_sequence, query_result in zip(query_sequences, results.query_hits):
        query_seq_str = sequence_to_genbank_base64_str(query_sequence) if genbank_sequence_columns else str(
            query_sequence.seq)
        for hit in query_result.hits:
            align_str = '{}|{}'.format(hit.query, hit.target)
            if hit.target_record and is_defined_sequence(hit.target_record.seq):
                hit_seq_str = sequence_to_genbank_base64_str(
                    hit.target_record) if genbank_sequence_columns else str(hit.target_record.seq)
            else:
                hit_seq_str = None
            table.data.append(
                [align_str, query_result.query_id, query_result.query_def, hit.target_id, hit.target_def,
                 hit.evalue, hit.score, hit.bits, query_seq_str, hit_seq_str])

    if genbank_sequence_columns:
        seq_data_type = 'binary'
        seq_column_type = 'chemical/x-genbank'
    else:
        seq_data_type = 'string'
        seq_column_type = 'chemical/x-sequence'
    table.columns = [{'name': 'Aligned Sequence', 'dataType': 'string',
                      'properties': {'ContentType': 'chemical/x-sequence-pair'},
                      'id': str(uuid.uuid4())},
                     {'name': 'Query Id', 'dataType': 'string'},
                     {'name': 'Query Definition', 'dataType': 'string'},
                     {'name': 'Target Id', 'dataType': 'string'},
                     {'name': 'Target Definition', 'dataType': 'string'},
                     {'name': 'Evalue', 'dataType': 'float'},
                     {'name': 'Score', 'dataType': 'integer'},
                     {'name': 'Bits', 'dataType': 'float'},
                     {'name': 'Query Sequence', 'dataType': seq_data_type,
                      'properties': {'ContentType': seq_column_type},
                      'id': str(uuid.uuid4())},
                     {'name': 'Target Sequence', 'dataType': seq_data_type,
                      'properties': {'ContentType': seq_column_type},
                      'id': str(uuid.uuid4())}]

    return table
