"""
===============
data_table.py
===============

Copyright (C) 2017-2022 Glysade, LLC

The data table class
"""
import re
import uuid

from typing import List, Dict, Union, Tuple, Optional

import ruse.util.ruse_resource
from ruse.util.frozen import Frozen
# Types for Data Tables- need to adhere to this for DataTableType and ColumnInfo structures
# or Spotfire won't be able to read.
from ruse.util.util import get_task_id

DataTypes = Union[str, int, float]
DataTableType = List[List[DataTypes]]
ColumnDataType = Union[str, Dict[str, Union[str, Dict[str, str]]]]
ColumnInfo = List[Dict[str, ColumnDataType]]

# Types for Json passed to/from Spotfire- Json should comprise

# input fields should be simple dictionary with settings and a single options dict which can be passed
# through to programs such as blast

# InputFieldValues = Union[DataTypes, Dict[str, DataTypes], Dict[str, Dict[str, DataTypes]]]
InputFieldValues = Dict[str, Dict[str, DataTypes]]
InputFields = Dict[str, InputFieldValues]

# Json should have three keys 'columns', 'data' and 'inputFields' with appropriate Json as values
JsonValues = Union[DataTableType, ColumnInfo, InputFields]
JsonType = Dict[str, JsonValues]


class DataTable(Frozen):
    """
    A class for representing tabular data


    Attributes:

        - data: A 2D array of all data in the table.  Each cell can have an float, int or string value- with columns having the same type.  Array is row major- i.e. implemented as a list of rows

        - columns: Column definition structure.  The columns attribute is a list of dictionaries, with the following entries:

             - name: column name

             - dataType: data type of the column (e.g. float, int or string)

             - properties: Dictionary of properties, including ContentType which is the mime type of the column (e.g. chemical/x-sequence-type)

             - id: column id.  Normally a generated UUID

        - output_columns: a list of output columns for use in task chaining

    """

    def __init__(self, columns: ColumnInfo = None, data: DataTableType = None, output_columns=None):
        """
        Constructor.  Call with no arguments to create an empty data table

        :param columns: defines any columns in the data table
        :param data: a 2D array to set cell data
        :param output_columns: used to specify any output columns
        """
        if output_columns is None:
            output_columns = {}
        self.output_columns = output_columns
        if columns is not None and data is not None:
            self.data = data
            self.columns = columns
        else:
            self.data = []
            self.columns = []

    def to_data(self) -> JsonType:
        """
        Convert to a raw data structure, suitable for encoding into JSON

        :return: Raw data representation
        """

        return {'columns': self.columns, 'data': self.data, 'outputColumns': self.output_columns, 'format': 'dataTable',
                'version': 1.0}

    def column_values(self, column_no: int) -> List[DataTypes]:
        """
        Return all values for a column

        :param column_no: the column index
        :return:
        """

        return [row[column_no] for row in self.data]

    @classmethod
    def column_definition(cls, name: str, data_type: str, content_type: str = None,
                          properties: Dict[str, str] = None) -> \
            ColumnDataType:
        """
        Creates a column definition for a column.  This should only be used when creating a new table, otherwise
        column names may not be unique

        :param name: the column name
        :param data_type: the column data type: one of float, int, binary or string
        :param content_type: the mime type or content type for the column
        :param properties: Optional additional properties
        :return: column definition as dictionary
        """
        if not properties:
            properties = {}
        column = {'name': name, 'dataType': data_type, 'id': str(uuid.uuid4()), 'properties': properties}
        if content_type and 'ContentType' not in column['properties']:
            column['properties']['ContentType'] = content_type
        return column

    def unique_column_name(self, name: str) -> str:
        """
        Creates a unique name for a column, assuming that it is the next column added.  If there is no column name in
        the data table that has the base column name, then the base name will be the column name.  Otherwise,
        the column number will be appended to the base column name.  If that name is in use then a integer will be
        found to append that creates a unique column name

        :param name: base name for the column
        :return: unique name for the column
        """

        def column_base_name(col_name: str) -> Tuple[str, int]:
            col_name = col_name.strip()
            match = re.match(r'^(.*?)\s*\((\d+)\)$', col_name)
            return (match.group(1), int(match.group(2))) if match else (col_name, 0)

        column_names = [column['name'] for column in self.columns if column['name']]
        if name in column_names:
            column_no = 2
            while "{} ({})".format(name, column_no) in column_names:
                column_no += 1
            name = "{} ({})".format(name, column_no)
        return name

    def add_column(self, base_name: str, data_type: str, props: Dict[str, str] = None, content_type: str = None,
                   fill_rows: bool = True, data: Optional[List[DataTypes]] = None) -> int:
        """
        Adds a column to the data table.  Does not change the class data attribute unless fill_rows is True

        :param base_name: a base name for the column
        :param data_type: the data type (float, integer or string) for the column
        :param props: dictionary of column properties
        :param content_type: the content type or mime type for the column
        :param fill_rows: set True (the default) to resize data to include empty cells (value None) for the new column
        :param data: optional data for the new column
        :return: the index of the new column
        """
        if not props:
            props = {}
        name = self.unique_column_name(base_name)

        # Many of the tests are run in an environment without a task id directory- in this case just
        # set the task_id to  1
        try:
            task_id = str(get_task_id())
        except ValueError:
            task_id = 1

        # add new column to output_columns dictionary
        if self.output_columns:
            if str(task_id) in self.output_columns:
                self.output_columns[str(task_id)].append({"name": name, "baseName": base_name})
            else:
                self.output_columns[str(task_id)] = [{"name": name, "baseName": base_name}]
        else:
            self.output_columns = {str(task_id): [{"name": name, "baseName": base_name}]}

        if content_type:
            props['ContentType'] = content_type
        column = {'name': name, 'dataType': data_type, 'properties': props.copy()}
        self.columns.append(column)
        if data:
            assert len(data) == len(self.data)
            for row, new_data in zip(self.data, data):
                row.append(new_data)
        elif fill_rows:
            for row in self.data:
                row.append(None)
        return len(self.columns) - 1

    def set_value(self, value: Union[str, int, float], row_idx: int, column_idx: int) -> None:
        """
        Set a cell value

        :param value: cell value
        :param column_idx: column index
        :param row_idx: row index
        """

        self.data[row_idx][column_idx] = value

    def content_type(self, column_idx: int):
        """
        Returns the content type for the column at an index

        :param column_idx: the column index
        :return: the content type for the column
        """
        return self.columns[column_idx]['properties']['ContentType']

    # importing RuseResource will create a circular reference
    def create_from_resource(self, resource: 'ruse.util.ruse_resource.RuseResource') -> None:
        """
        Initializes data table contents for a resource.  The data table will comprise one cell containing the resource id

        :param resource: ruse resource as :class:`ruse.util.ruse_resource.RuseResource`
        """

        assert self.data == []
        assert self.columns == []

        self.data = [[resource.resource_id]]
        self.columns.append(self.column_definition('Resource ID', 'int', 'ruse/x-resource'))

    def to_json(self) -> JsonType:
        """
        Converts the datatype to a simple python data structure that can be converted to JSON

        :return: Json compatible data
        """
        result = {'format': 'dataTable', 'version': 1.0}
        if self.table:
            result['table'] = self.table
        if self.columns:
            result['columns'] = self.columns
        if self.data:
            result['data'] = self.data
        if self.output_columns:
            result['outputColumns'] = self.output_columns
        return result

    def __iadd__(self, other):
        """
        Adds/appends another data table to this one, using row axis.  Columns are mapped from the other data table to
        this one using names and rows appended

        :param other: the data table to append to this one
        :return: this data table
        """

        dst_columns = [column['name'] for column in self.columns]
        column_idx_map = {}
        for column in other.columns:
            if column['name'] in dst_columns:
                column_idx_map[column['name']] = dst_columns.index(column['name'])
            else:
                column_idx_map[column['name']] = self.add_column(column['name'], column['dataType'],
                                                                 column['properties'])

        ncols = len(self.columns)
        for srcrow in other.data:
            dstrow = [None] * ncols
            for idx, column in enumerate(other.columns):
                dstrow[column_idx_map[column['name']]] = srcrow[idx]
            self.data.append(dstrow)

        if other.output_columns:
            if not self.output_columns:
                self.output_columns = other.output_columns
            else:
                for key, val in other.output_columns.items():
                    if key not in self.output_columns:
                        self.output_columns[key] = val

        return self

    def append_columns(self, other) -> None:
        """
        Appends the columns in another table to this one.  Assumes each table has the same number of rows

        :param other: the other data table
        """

        assert (len(self.data) == len(other.data))

        # TODO make sure column names are unique
        self.columns.extend(other.columns)
        for row, other_row in zip(self.data, other.data):
            row.extend(other_row)

    def get_column_index(self, name: str):
        """
        Return the index of the column with the given name, or -1 if no column is found

        :param name: column name
        :return: index of column in data table
        """
        for idx, column in enumerate(self.columns):
            if column['name'] == name:
                return idx
        return -1
