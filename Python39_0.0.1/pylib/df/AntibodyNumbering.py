import json

from df.bio_helper import values_to_sequences
from df.data_transfer import ColumnData, DataFunctionRequest, DataFunctionResponse, DataFunction, string_input_field, \
    DataType
from ruse.bio.antibody import align_antibody_sequences
from ruse.bio.bio_data_table_helper import sequence_to_genbank_base64_str


class AntibodyNumbering(DataFunction):
    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:
        input_column = next(iter(request.inputColumns.values()))
        input_column.remove_nulls()
        input_sequences = values_to_sequences(input_column)
        scheme = string_input_field(request, 'scheme', 'chothia')
        align_information = align_antibody_sequences(input_sequences, scheme)
        output_sequences = align_information.aligned_sequences
        numbering = align_information.numbering
        regions = align_information.regions

        numbering_data = [{'domain': n.domain, 'position': n.query_position, 'label': n.label()} for n in numbering]
        region_data = [r.to_data() for r in regions]
        antibody_numbering = {
            'scheme': scheme,
            'numbering': numbering_data,
            'regions': region_data
        }
        numbering_json = json.dumps(antibody_numbering)

        rows = [sequence_to_genbank_base64_str(s) for s in output_sequences]
        properties = {'antibodyNumbering': numbering_json}
        output_column = ColumnData(name='Aligned Sequence', dataType=DataType.BINARY, properties=properties,
                                   contentType='chemical/x-genbank', values=rows)
        output_column.insert_nulls(input_column.missing_null_positions)
        response = DataFunctionResponse(outputColumns=[output_column])
        return response
