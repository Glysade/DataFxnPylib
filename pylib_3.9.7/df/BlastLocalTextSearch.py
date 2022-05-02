import os
from typing import List, Optional

from Bio.SeqRecord import SeqRecord

from df.bio_helper import sequences_to_column, query_from_request
from df.data_transfer import DataFunction, DataFunctionRequest, DataFunctionResponse, string_input_field, \
    integer_input_field, ColumnData, DataType, TableData, boolean_input_field
from ruse.bio.antibody import numbering_and_regions_from_sequence, ANTIBODY_NUMBERING_COLUMN_PROPERTY
from ruse.bio.bio_data_table_helper import sequence_to_genbank_base64_str
from ruse.bio.bio_util import is_defined_sequence
from ruse.bio.blast_parse import MultipleBlastResults, build_common_alignments
from ruse.bio.blast_search import BlastSearch, BlastSearchType
from ruse.bio.blast_utils import blast_databases_list
from ruse.util.log import error


def run_blast_search(sequences: List[SeqRecord], request: DataFunctionRequest, output_table_name: str,
                     sequence_column_type: Optional[str] = None, show_query_id: bool = False):
    database_name = string_input_field(request, 'databaseName')
    method = string_input_field(request, 'method')
    max_hits = integer_input_field(request, 'maxHits')
    show_multiple_alignments = boolean_input_field(request, 'showMultipleAlignments', False)
    if not sequences or len(sequences) > 1:
        show_multiple_alignments = False
    if method in ['TBLASTX', 'BLASTX']:
        show_multiple_alignments = False
    search = BlastSearch()
    options = {}
    if max_hits:
        options['max_target_seqs'] = max_hits

    available_database = blast_databases_list()
    if database_name not in available_database:
        msg = 'Blast database "{}" not present'.format(database_name)
        error(msg)
        raise ValueError(msg)

    search_type = BlastSearchType.from_string(method)

    results: Optional[MultipleBlastResults] = None
    if sequences:
        search.multiple_query_search_blast_database(sequences, database_name, search_type, None, options)
        if search.error is not None:
            raise ValueError(f'Blast search failed with error {search.error}')
        results = MultipleBlastResults()
        results.parse(search.output_file())

        results.retrieve_local_targets()
        results.complement_target_sequence()
        results.truncate_target_sequences(10_000)
        assert len(sequences) == len(results.query_hits)

    if not sequence_column_type:
        sequence_column_type = 'chemical/x-sequence'
    genbank = sequence_column_type == 'chemical/x-genbank'

    alignments = []
    target_ids = []
    target_definitions = []
    query_ids = []
    query_definitions = []
    e_values = []
    scores = []
    bits = []
    query_sequences = []
    target_sequences = []
    multiple_alignments = []
    if show_multiple_alignments:
        alignments.append(None)
        target_ids.append(None)
        target_definitions.append(None)
        query_ids.append(None)
        query_definitions.append(None)
        e_values.append(None)
        scores.append(None)
        bits.append(None)
        query_sequences.append(None)
        target_sequences.append(None)

    if sequences:
        for query_sequence, query_result in zip(sequences, results.query_hits):
            query_seq_str = sequence_to_genbank_base64_str(query_sequence) if genbank else str(
                query_sequence.seq)
            for hit in query_result.hits:
                align_str = '{}|{}'.format(hit.query, hit.target)
                if hit.target_record and is_defined_sequence(hit.target_record.seq):
                    hit_seq_str = sequence_to_genbank_base64_str(
                        hit.target_record) if genbank else str(hit.target_record.seq)
                else:
                    hit_seq_str = None

                alignments.append(align_str)
                target_ids.append(hit.target_id)
                target_definitions.append(hit.target_def)
                query_ids.append(query_sequence.id)
                query_definitions.append(query_sequence.description)
                e_values.append(hit.evalue)
                scores.append(hit.score)
                bits.append(hit.bits)
                query_sequences.append(query_seq_str)
                target_sequences.append(hit_seq_str)

            if show_multiple_alignments:
                multiple_alignments = build_common_alignments(query_sequence, query_result.hits)

    sequence_data_type = DataType.BINARY if genbank else DataType.STRING
    aligned_sequences_column = ColumnData(name='Aligned Sequence Pairs', dataType=DataType.STRING,
                                          contentType='chemical/x-sequence-pair', values=alignments)
    target_id_column = ColumnData(name='Target Id', dataType=DataType.STRING, values=target_ids)
    target_definition_column = ColumnData(name='Target Definition', dataType=DataType.STRING,
                                          values=target_definitions)
    query_id_column = ColumnData(name='Query Id', dataType=DataType.STRING, values=query_ids)
    query_definition_column = ColumnData(name='Query Definition', dataType=DataType.STRING,
                                         values=query_definitions)
    e_value_column = ColumnData(name='EValue', dataType=DataType.DOUBLE, values=e_values)
    score_column = ColumnData(name='Score', dataType=DataType.LONG, values=scores)
    bit_column = ColumnData(name='Bits', dataType=DataType.DOUBLE, values=bits)
    query_sequence_column = ColumnData(name='Query Sequence', dataType=sequence_data_type,
                                       contentType=sequence_column_type, values=query_sequences)
    target_sequence_column = ColumnData(name='Target Sequence', dataType=sequence_data_type,
                                        contentType=sequence_column_type, values=target_sequences)
    columns = [aligned_sequences_column,
               target_id_column,
               target_definition_column,
               e_value_column, score_column, bit_column,
               query_sequence_column,
               target_sequence_column]
    if show_query_id:
        columns.insert(3, query_definition_column)
        columns.insert(3, query_id_column)
    if show_multiple_alignments:
        multiple_alignment_column = sequences_to_column(multiple_alignments, 'Aligned Sequence', True)
        antibody_numbering = numbering_and_regions_from_sequence(multiple_alignments[0])
        if antibody_numbering:
            multiple_alignment_column.properties[
                ANTIBODY_NUMBERING_COLUMN_PROPERTY] = antibody_numbering.to_column_json()
        columns.insert(0, multiple_alignment_column)
    output_table = TableData(tableName=output_table_name,
                             columns=columns)
    response = DataFunctionResponse(outputTables=[output_table])
    return response


class BlastLocalTextSearch(DataFunction):
    """
    Performs a BLAST search of user selected databases using a query entered as text
    """

    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:
        query_sequence = query_from_request(request)
        blastdb = string_input_field(request, 'blastDbPath')
        os.environ['BLASTDB'] = blastdb
        return run_blast_search([query_sequence], request, 'Blast text search results')
