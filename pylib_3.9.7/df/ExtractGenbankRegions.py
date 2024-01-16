from typing import Optional, Dict

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from df.bio_helper import column_to_sequences, sequences_to_column
from df.data_transfer import DataFunction, DataFunctionRequest, DataFunctionResponse, DataType, \
                             ColumnData, TableData, integer_input_field, string_input_field, input_field_to_column


def _feature_to_name(feature: SeqFeature, feature_key: str, feature_qualifier: str) -> Optional[str]:
    if feature.type != feature_key:
        return None
    if feature_qualifier:
        if feature_qualifier not in feature.qualifiers:
            return None

        values = feature.qualifiers[feature_qualifier]
        if not values or len(values) != 1:
            return None
    else:
        values = [feature_key]

    return values[0]


def _extract_region(sequence: Optional[SeqRecord], feature_key: str, feature_qualifier: str, name: str) \
        -> [Optional[Seq]]:
    if not sequence:
        return None

    feature_sequences = []
    feature: SeqFeature
    for feature in sequence.features:
        feature_name = _feature_to_name(feature, feature_key, feature_qualifier)

        if not feature_name:
            continue

        if feature_name == name:
            try:
                seq = feature.extract(sequence.seq)
            except Exception as ex:
                # need to create a notification here
                continue

            # return SeqRecord(seq, name)
            feature_sequences.append(SeqRecord(seq, name, features = [feature]))

    return feature_sequences
    # return None


class ExtractGenbankRegions(DataFunction):
    """
    Extracts all regions matching a feature key and qualifier from a column of genbank sequences into new columns
    """

    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:
        input_column = next(iter(request.inputColumns.values()))
        input_column.remove_nulls()
        input_sequences = column_to_sequences(input_column)

        sequence_ids = input_field_to_column(request, 'uiIDColumn')
        if sequence_ids:
            sequence_id_column_name = sequence_ids.name
            sequence_ids_no_nulls = [seq_id for idx, seq_id in enumerate(sequence_ids.values)
                                     if idx not in input_column.missing_null_positions]
            # create fake ids for any missing ids
            sequence_ids_no_nulls = [seq_id if seq_id else f'Row {i + 1}' for i, seq_id in enumerate(sequence_ids_no_nulls)]
        else:
            # create fake ids if no ID column was selected
            sequence_id_column_name = 'ID'
            sequence_ids_no_nulls = [f'Row {i + 1}' for i in range(len(input_sequences))]

        max_number_of_regions = integer_input_field(request, 'maximumNumberRegions')

        feature_key = string_input_field(request, 'featureKey')
        feature_qualifier = string_input_field(request, 'featureQualifier')

        # first find all regions - dict preserves order
        region_names_dict: Dict[str] = {}
        for sequence, sequence_id in zip(input_sequences, sequence_ids_no_nulls):
            if not sequence:
                continue
            feature: SeqFeature
            for feature in sequence.features:
                name = _feature_to_name(feature, feature_key, feature_qualifier)
                if name and name not in region_names_dict:
                    region_names_dict[name] = None
                    if len(region_names_dict) >= max_number_of_regions:
                        break

        region_names = list(region_names_dict.keys())
        output_columns = []
        feature_sequence_ids = []

        for region_name in region_names:
            feature_sequences = [_extract_region(s, feature_key, feature_qualifier, region_name)
                                 for s in input_sequences]
            feature_sequence_ids = [i for i, seq in zip(sequence_ids_no_nulls, input_sequences) if seq]


            if feature_qualifier:
                # verify handling of list of lists
                output_column = sequences_to_column([feature_sequence[0] if len(feature_sequence) != 0 else None
                                                     for feature_sequence in feature_sequences],
                                        f'{input_column.name} {region_name}', False)
                output_column.insert_nulls(input_column.missing_null_positions)
                output_columns.append(output_column)

                response = DataFunctionResponse(outputColumns = output_columns)
            else:
                output_feature_sequences = []
                output_feature_ids = []
                output_feature_descr = []
                output_feature_start = []
                output_feature_end = []
                for feature_seqs, sequence_id in zip(feature_sequences, feature_sequence_ids):
                    output_feature_ids.extend([sequence_id] * len(feature_seqs))
                    output_feature_sequences.extend(feature_seqs)

                    for feature_seq in feature_seqs:
                        qualifiers = [v for k, v in feature_seq.features[0].qualifiers.items()]
                        output_feature_descr.append(qualifiers[0][0])
                        output_feature_start.append(feature_seq.features[0].location.start)
                        output_feature_end.append(feature_seq.features[0].location.end)
                        # output_feature_start.extend([feature_seq.features[0].location.start for feature_seq in feature_seqs])
                        # output_feature_end.extend([feature_seq.features[0].location.end for feature_seq in feature_seqs])

                output_columns = [
                                  ColumnData(name = sequence_id_column_name, dataType = DataType.STRING,
                                             values = output_feature_ids),
                                  sequences_to_column(output_feature_sequences, f'{input_column.name} {region_name}', False),
                                  ColumnData(name = 'Description', dataType = DataType.STRING,
                                             values = output_feature_descr),
                                  ColumnData(name = 'Start', dataType = DataType.INTEGER,
                                             values = output_feature_start),
                                  ColumnData(name = 'end', dataType = DataType.INTEGER,
                                             values = output_feature_end)
                ]

                output_table = TableData(tableName = f'{input_column.name} {feature_key}', columns = output_columns)

                response = DataFunctionResponse(outputTables = [output_table])

        return response
