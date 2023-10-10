import json
import os

from df.bio_helper import column_to_sequences
from df.data_transfer import ColumnData, DataFunctionRequest, DataFunctionResponse, DataFunction, string_input_field, \
    DataType, input_field_to_column

from ImmuneBuilder import ABodyBuilder2

class AntibodyStructurePrediction(DataFunction):
    """
    Predicts antibody structure from heavy and light chanin sequences
    """

    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:

        predictor = ABodyBuilder2()

        light_chains = column_to_sequences(input_field_to_column(request, 'lightChain'))
        heavy_chains = column_to_sequences(input_field_to_column(request, 'heavyChain'))
        ids = input_field_to_column(request, 'idColumn').values
        output_dir = string_input_field(request, 'outputDirectory')

        filenames = []

        for light, heavy, id in zip(light_chains, heavy_chains, ids):
            if light is None or heavy is None:
                filenames.append(None)
            else:
                sequences = {'H': str(heavy.seq).upper(), 'L': str(light.seq).upper()}
                antibody = predictor.predict(sequences)

                filename = os.path.join(output_dir, '_'.join([id, 'predicted.pdb']))
                filenames.append(filename)
                antibody.save(filename)

        output_column = ColumnData(name = 'Structure Files', dataType = DataType.STRING,
                                   contentType = 'chemical/x-uri', values = filenames)

        return DataFunctionResponse(outputColumns = [output_column])