import uuid

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from df.data_transfer import DataFunction, DataFunctionRequest, DataFunctionResponse, string_input_field, ColumnData, \
    DataType, integer_input_field, TableData
from ruse.bio.bio_data_table_helper import sequence_to_genbank_base64_str
from ruse.bio.blast_parse import BlastResults
import ruse.bio.blast_search


class BlastWebSearch(DataFunction):
    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:
        args = {}
        query_input_field = string_input_field(request, 'query')
        if not query_input_field:
            raise ValueError()
        query = SeqRecord(Seq(query_input_field), id=str(uuid.uuid4()), name='', description='')

        database_name = string_input_field(request, 'databaseName')
        max_hits = integer_input_field(request, 'maxHits')

        args['hitlist_size'] = max_hits
        args['database_name'] = database_name

        search = ruse.bio.blast_search.BlastWebSearch()
        search.search_blast_database(query, **args)
        results = BlastResults()
        results.parse(search.output_file())
        results.retrieve_targets()
        results.complement_target_sequence()

        include_pdb_uris = database_name == 'pdb'
        # query_seq_str = str(query.seq)
        alignments = []
        target_ids = []
        target_definitions = []
        e_values = []
        scores = []
        bits = []
        query_sequences = []
        target_sequences = []
        uris = []

        for hit in results.hits:
            align_str = '{}|{}'.format(hit.query, hit.target)
            genbank_str = sequence_to_genbank_base64_str(hit.target_record) if hit.target_record else None

            alignments.append(align_str)
            target_ids.append(hit.target_id)
            target_definitions.append(hit.target_def)
            scores.append(hit.score)
            e_values.append(hit.evalue)
            bits.append(hit.bits)
            # query_sequences.append(query_seq_str)
            query_sequences.append(hit.query)
            target_sequences.append(genbank_str)
            if include_pdb_uris:
                parts = hit.target_id.split('|')
                if len(parts) == 5 and parts[2] == 'pdb':
                    uri = 'https://files.rcsb.org/download/%s.pdb' % parts[3]
                elif len(parts) == 3 and parts[0] == 'pdb':
                    uri = 'https://files.rcsb.org/download/%s.pdb' % parts[1]
                else:
                    uri = None
                uris.append(uri)

        aligned_sequences_column = ColumnData(name='Aligned Sequence', dataType=DataType.STRING,
                                              contentType='chemical/x-sequence-pair', values=alignments)
        target_id_column = ColumnData(name='Target Id', dataType=DataType.STRING, values=target_ids)
        target_definition_column = ColumnData(name='Target Definition', dataType=DataType.STRING,
                                              values=target_definitions)
        e_value_column = ColumnData(name='EValue', dataType=DataType.FLOAT, values=e_values)
        score_column = ColumnData(name='Score', dataType=DataType.INTEGER, values=scores)
        bit_column = ColumnData(name='Bits', dataType=DataType.FLOAT, values=bits)
        query_sequence_column = ColumnData(name='Query Sequence', dataType=DataType.STRING,
                                           contentType='chemical/x-sequence', values=query_sequences)
        target_sequence_column = ColumnData(name='Target Sequence', dataType=DataType.BINARY,
                                            contentType='chemical/x-genbank', values=target_sequences)
        columns = [aligned_sequences_column, target_id_column, target_definition_column, e_value_column, score_column,
                   bit_column, query_sequence_column, target_sequence_column]
        if include_pdb_uris:
            uri_column = ColumnData(name='URL', dataType=DataType.STRING, contentType='chemical/x-uri',
                                    properties={'Dimension': '3'}, values=uris)
            columns.append(uri_column)
        output_table = TableData(tableName='Blast search results',
                                 columns=columns)
        response = DataFunctionResponse(outputTables=[output_table])
        return response
