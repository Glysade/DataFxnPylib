import ruse.bio.blast_search
from df.bio_helper import sequences_to_column, query_from_request
from df.data_transfer import DataFunction, DataFunctionRequest, DataFunctionResponse, string_input_field, ColumnData, \
    DataType, integer_input_field, TableData, boolean_input_field
from ruse.bio.bio_data_table_helper import sequence_to_genbank_base64_str
from ruse.bio.blast_parse import BlastResults, build_common_alignments


class BlastWebSearch(DataFunction):
    """
    Perform a remote Blast sequence search against an Entrez database
    """

    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:
        args = {}
        query = query_from_request(request)

        database_name = string_input_field(request, 'databaseName')
        max_hits = integer_input_field(request, 'maxHits')
        show_multiple_alignments = boolean_input_field(request, 'showMultipleAlignments', False)

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
        multiple_alignments = []
        if show_multiple_alignments:
            alignments.append(None)
            target_ids.append(None)
            target_definitions.append(None)
            e_values.append(None)
            scores.append(None)
            bits.append(None)
            query_sequences.append(None)
            target_sequences.append(None)
            uris.append(None)

        for hit in results.hits:
            align_str = '{}|{}'.format(hit.query, hit.target)
            genbank_str = sequence_to_genbank_base64_str(hit.target_record) if hit.target_record else None

            alignments.append(align_str)
            target_ids.append(hit.target_id)
            target_definitions.append(hit.target_def)
            scores.append(hit.score)
            e_values.append(hit.evalue)
            bits.append(hit.bits)
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

        if show_multiple_alignments:
            multiple_alignments = build_common_alignments(query, results.hits)
        aligned_sequences_column = ColumnData(name='Sequence Pairs', dataType=DataType.STRING,
                                              contentType='chemical/x-sequence-pair', values=alignments)
        target_id_column = ColumnData(name='Target Id', dataType=DataType.STRING, values=target_ids)
        target_definition_column = ColumnData(name='Target Definition', dataType=DataType.STRING,
                                              values=target_definitions)
        e_value_column = ColumnData(name='EValue', dataType=DataType.DOUBLE, values=e_values)
        score_column = ColumnData(name='Score', dataType=DataType.LONG, values=scores)
        bit_column = ColumnData(name='Bits', dataType=DataType.DOUBLE, values=bits)
        query_sequence_column = ColumnData(name='Query Sequence', dataType=DataType.STRING,
                                           contentType='chemical/x-sequence', values=query_sequences)
        target_sequence_column = ColumnData(name='Target Sequence', dataType=DataType.BINARY,
                                            contentType='chemical/x-genbank', values=target_sequences)
        columns = [aligned_sequences_column, target_id_column, target_definition_column, e_value_column, score_column,
                   bit_column, query_sequence_column, target_sequence_column]
        if show_multiple_alignments:
            multiple_alignment_column = sequences_to_column(multiple_alignments, 'Aligned Sequence', True)
            columns.insert(0, multiple_alignment_column)
        if include_pdb_uris:
            uri_column = ColumnData(name='URL', dataType=DataType.STRING, contentType='chemical/x-uri',
                                    properties={'Dimension': '3'}, values=uris)
            columns.append(uri_column)
        output_table = TableData(tableName='Blast search results',
                                 columns=columns)
        response = DataFunctionResponse(outputTables=[output_table])
        return response
