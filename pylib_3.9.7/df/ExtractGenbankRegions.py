from typing import Optional, Dict

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from df.bio_helper import column_to_sequences, sequences_to_column
from df.data_transfer import DataFunction, DataFunctionRequest, DataFunctionResponse, integer_input_field, \
    string_input_field


def _feature_to_name(feature: SeqFeature, feature_key: str, feature_qualifier: str) -> Optional[str]:
    if feature.type != feature_key:
        return None
    if feature_qualifier not in feature.qualifiers:
        return None
    values = feature.qualifiers[feature_qualifier]
    if not values or len(values) != 1:
        return None
    return values[0]


def _extract_region(sequence: Optional[SeqRecord], feature_key: str, feature_qualifier: str, name: str) -> Optional[
    Seq]:
    if not sequence:
        return None
    feature: SeqFeature
    for feature in sequence.features:
        feature_name = _feature_to_name(feature, feature_key, feature_qualifier)
        if not feature_name:
            continue
        if feature_name == name:
            seq = feature.extract(sequence.seq)
            return SeqRecord(seq, name)
    return None


class ExtractGenbankRegions(DataFunction):

    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:
        input_column = next(iter(request.inputColumns.values()))
        input_column.remove_nulls()
        input_sequences = column_to_sequences(input_column)
        max_number_of_regions = integer_input_field(request, 'maximumNumberRegions')
        feature_key = string_input_field(request, 'featureKey')
        feature_qualifier = string_input_field(request, 'featureQualifier')

        # first find all regions - dict preserves order
        region_names_dict: Dict[str] = {}
        finished = False
        for sequence in input_sequences:
            if finished:
                break
            if not sequence:
                continue
            feature: SeqFeature
            for feature in sequence.features:
                if finished:
                    break
                name = _feature_to_name(feature, feature_key, feature_qualifier)
                if name and name not in region_names_dict:
                    region_names_dict[name] = None
                    if len(region_names_dict) >= max_number_of_regions:
                        break

        region_names = list(region_names_dict.keys())

        output_columns = []
        for region_name in region_names:
            feature_sequences = [_extract_region(s, feature_key, feature_qualifier, region_name) for s in
                                 input_sequences]
            output_column = sequences_to_column(feature_sequences, f'{input_column.name} {region_name}', False)
            output_columns.append(output_column)
        response = DataFunctionResponse(outputColumns=output_columns)
        return response
