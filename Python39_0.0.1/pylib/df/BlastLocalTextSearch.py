from typing import List, Optional

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from df.data_transfer import DataFunction, DataFunctionRequest, DataFunctionResponse, string_input_field, \
    integer_input_field, ColumnData, DataType, TableData
from ruse.bio.bio_data_table_helper import sequence_to_genbank_base64_str
from ruse.bio.bio_util import is_defined_sequence
from ruse.bio.blast_parse import MultipleBlastResults
from ruse.bio.blast_search import BlastSearch, BlastSearchType
from ruse.bio.blast_utils import blast_databases_list
from ruse.util.log import error


def run_blast_search(sequences: List[SeqRecord], request: DataFunctionRequest,
                     sequence_column_type: Optional[str] = None):
    database_name = string_input_field(request, 'databaseName')
    method = string_input_field(request, 'method')
    max_hits = integer_input_field(request, 'maxHits')

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

    search.multiple_query_search_blast_database(sequences, database_name, search_type, None, options)
    if search.error is not None:
        return {'error': search.error}
    results = MultipleBlastResults()
    results.parse(search.output_file())

    results.retrieve_local_targets()
    results.complement_target_sequence()
    results.truncate_target_sequences(10_000)

    if not sequence_column_type:
        sequence_column_type = 'chemical/x-sequence'
    genbank = sequence_column_type == 'chemical/x-genbank'
    assert len(sequences) == len(results.query_hits)

    alignments = []
    query_ids = []
    query_definitions = []
    target_ids = []
    target_definitions = []
    e_values = []
    scores = []
    bits = []
    query_sequences = []
    target_sequences = []

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
            query_ids.append(query_result.query_id)
            query_definitions.append(query_result.query_def)
            target_ids.append(hit.target_id)
            target_definitions.append(hit.target_def)
            e_values.append(hit.evalue)
            scores.append(hit.score)
            bits.append(hit.bits)
            query_sequences.append(query_seq_str)
            target_sequences.append(hit_seq_str)

    sequence_data_type = DataType.BINARY if genbank else DataType.STRING
    aligned_sequences_column = ColumnData(name='Aligned Sequence', dataType=DataType.STRING,
                                          contentType='chemical/x-sequence-pair', values=alignments)
    query_id_column = ColumnData(name='Query Id', dataType=DataType.STRING, values=query_ids)
    query_definition_column = ColumnData(name='Query Definition', dataType=DataType.STRING,
                                         values=query_definitions)
    target_id_column = ColumnData(name='Target Id', dataType=DataType.STRING, values=target_ids)
    target_definition_column = ColumnData(name='Target Definition', dataType=DataType.STRING,
                                          values=target_definitions)
    e_value_column = ColumnData(name='EValue', dataType=DataType.FLOAT, values=e_values)
    score_column = ColumnData(name='Score', dataType=DataType.INTEGER, values=scores)
    bit_column = ColumnData(name='Bits', dataType=DataType.FLOAT, values=bits)
    query_sequence_column = ColumnData(name='Query Sequence', dataType=sequence_data_type,
                                       contentType=sequence_column_type, values=query_sequences)
    target_sequence_column = ColumnData(name='Target Sequence', dataType=sequence_data_type,
                                        contentType=sequence_column_type, values=target_sequences)
    output_table = TableData(tableName='Blast search results',
                             columns=[aligned_sequences_column, query_id_column, query_definition_column,
                                      target_id_column,
                                      target_definition_column, e_value_column, score_column, bit_column,
                                      query_sequence_column,
                                      target_sequence_column])
    response = DataFunctionResponse(outputTables=[output_table])
    return response


class BlastLocalTextSearch(DataFunction):

    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:
        query = string_input_field(request, 'query')
        if not query:
            raise ValueError()
        sequence = SeqRecord(Seq(query), 'Query')
        return run_blast_search([sequence], request)
