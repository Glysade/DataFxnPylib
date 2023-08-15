"""
================
data_transfer.py
================

Copyright (C) 2022 Glysade, LLC

Data transfer model objects and utility functions for Glysade Python data functions

We use classes derived from `Pydantic Models <https://pydantic-docs.helpmanual.io/>`_.  This allows request and
response validation and serialization to and from JSON.
"""

from abc import ABC, abstractmethod
from enum import Enum
from typing import Optional, Any, Dict, List

from pydantic import BaseModel


# I wanted to use this type instead of any for column/inputField values, but it coerces values to strings when used
# with pydantic
# DataValue = Union[str, int, float, bool]


class DataType(Enum):
    """
    A string enum to specify Spotfire column data types and input field data types.  The Spotfire column type and
    equivalent Python type are listed below.

    Members:
       * BOOLEAN: Spotfire boolean data type, associated value should be Python bool
       * STRING: Spotfire string data type, Python str
       * BINARY: Spotfire binary data type, Python str value is the result of gzipping then Base64 encoding the
         original raw value.
       * INTEGER: integer data type, equivalent to Spotfire 32 bit integer column type (Integer). Python int.
         Since Python ints may be 64 bit numbers use this as a return type with care.
       * LONG: long data type, equivalent to Spotfire 64 bit integer column type (LongInteger).  Python int.
         Since Python will automatically use
         64 bit integers when required, use this type for return columns to avoid data truncation
       * FLOAT: 64 bit double data type.  Python float
       * DOUBLE: 64 bit double data type.  Python float. There is no difference between the double and float data type
       * STRING_LIST: a JSON list of strings.  Used in input fields only
       * INTEGER_LIST: a JSON list of integers.  Used in input fields only
       * DOUBLE_lIST: a JSON list of double.  Used in input fields only
    """
    BOOLEAN = 'boolean'
    STRING = 'string'
    BINARY = 'binary'
    INTEGER = 'integer'
    LONG = 'long'
    FLOAT = 'float'
    DOUBLE = 'double'
    STRING_LIST = 'list(string)'
    INTEGER_LIST = 'list(integer)'
    DOUBLE_LIST = 'list(double)'


class InputFieldSelectorType(Enum):
    """
    A string enum to specify if an input field is used to select input data from Spotfire

    Members:
        * COLUMN: the input field is used to select data in a column
        * TABLE: the input fields is used to select a table (not currently implemented)
    """
    COLUMN = 'column'
    TABLE = 'table'


class InputField(BaseModel):
    """
    A class to represent a data function input field.  This is not normally constructed by the user but is
    created by parsing the data function input JSON.
    """
    id: str
    """
    The id name for the input field (as defined in the data function definition)`
    """
    contentType: Optional[str]
    """
    Optional content type for the input field (e.g. chemical/x-smiles)
    """
    dataType: DataType
    """
    The data type of the input field
    """
    data: Any
    """
    The value for the input field.  The type used is based on the dataType (:class:`DataType`) attribute
    """
    selectorType: Optional[InputFieldSelectorType]
    """
    Optional selector type.  For example, the input field may be used to select a column, in which case the 
    data value is the column id.
    """
    limitBy: Optional[str]
    """
    Any Limit By set for an input column.  Available only in Charts version 2.3.0 or later
    """
    limitByMarking: Optional[str]
    """
    The marking (if any) used by limitBy. Available only in Charts version 2.3.0 or later
    """


class ColumnData(BaseModel):
    """
    A data model class containing a Spotfire column.
    Used when the data function contains column inputs and should be instantiated when the data function is
    returning column outputs
    """
    name: str
    """
    The column name
    """
    dataType: DataType
    """
    The data type stored in the column
    """
    properties: Dict[str, str] = {}
    """
    Spotfire column properties
    """
    contentType: Optional[str]
    """
    Optional content type for the column (e.g. chemical/s-smiles)
    """
    values: List[Optional[Any]]
    """
    A list of row values for the column. The type used is based on the dataType (:class:`DataType`) attribute or None if
    the cell is empty
    """
    missing_null_positions: List[int] = []
    """
    An array to track null values if they have been removed
    """

    columnId: Optional[str]
    """
    The Spotfire input Column Id.  Set for input columns only in Charts version 2.3.0 or later
    """

    limitBy: Optional[str]
    """
    Any Limit By set for an input column.  Available only in Charts version 2.3.0 or later
    """

    limitByMarking: Optional[str]
    """
    The marking (if any) used by limitBy. Available only in Charts version 2.3.0 or later
    """

    # TODO: maybe just use column data as transfer object move this functionality to another class or utility functions
    def remove_nulls(self) -> None:
        """
        Remove null data values and store their position in the :attr:`missing_null_positions` attribute
        """
        positions = [i for i, v in enumerate(self.values) if v is None]
        positions.reverse()
        for p in positions:
            self.values.pop(p)
        positions.reverse()
        self.missing_null_positions = positions

    def insert_nulls(self, positions: Optional[List[int]] = None, offset: int = 0) -> None:
        """
        Inserts nulls (or None) into the values array at the specified positions.
        The attribute :attr:`missing_null_positions` will be used if the position argument is not supplied.

        For example, you may remove null values from an input column, create an output column then use the
        :attr:`missing_null_positions` attribute from the input column as an argument to correctly insert nulls into
        the output column.

        :param positions: optional list of null positions
        :param offset: optional offset for null positions
        """
        if not positions:
            positions = self.missing_null_positions
        for p in positions:
            self.values.insert(p + offset, None)
        self.missing_null_positions = []


class DataFunctionDocumentProperty(BaseModel):
    """
    A Spotfire document property.  Either an existing property in the request, or a property in the reponse to be
    added or updated.
    """

    name: str
    """
    The name of the document property
    """

    dataType: DataType
    """
    The data type stored in the property
    """

    description: Optional[str]
    """
    An optional description for the property
    """

    value: Any
    """
    The value of the document property
    """


class TableData(BaseModel):
    """
    Model class for Spotfire data table.  Instantiate this to return a new table from a data function
    """
    tableName: str
    """
    A suggested name for the data table
    """
    columns: List[ColumnData]
    """
    A list of columns in the data table
    """


class DataFunctionRequest(BaseModel):
    """
    Model class for data function requests
    """
    inputColumns: Dict[str, ColumnData] = []
    """
    Input columns in the request.  This dictionary is keyed by column id and limit by.
    """
    inputFields: Dict[str, InputField] = {}
    """
    Input fields keyed by input field :attr:`InputField.id`.  If an input field value has an associated input column then
    :attr:`InputField.selectorType` will be COLUMN and :attr:`InputField.id` will be the column id (that can be used to 
    retrieve column data from :attr:`inputColumns`).
    """
    serviceName: str
    """
    Data function service name.  If we are running a Python script this will be Script. Otherwise it is the name of
    the module and class that will run the data function.
    """
    id: str
    """
    Unique UUID4 identifier for the request
    """
    script: Optional[str]
    """
    An optional Python script that contains code to execute the data function.
    """
    documentProperties: List[DataFunctionDocumentProperty] = []
    """
    An array of the current document properties
    """


class DataFunctionResponse(BaseModel):
    """
    Model class for data function responses
    """
    outputColumns: List[ColumnData] = []
    """
    An array of output columns (to be added to the data function's input table)
    """

    outputTables: List[TableData] = []
    """
    An array of output tables
    """

    documentProperties: List[DataFunctionDocumentProperty] = []
    """
    An array of document properties (to add or update the current document properties)
    """


class DataFunction(ABC):
    """
    Base class to be used when defining builtin data function classes that are run from their module files
    """

    @abstractmethod
    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:
        """
        Abstract method that runs a data function.

        :param request:
        """
        pass


def _input_field_value(request: DataFunctionRequest, field: str) -> Optional[Any]:
    input_fields = request.inputFields
    if field not in input_fields:
        return None
    return input_fields[field].data


def _check_data_type(request: DataFunctionRequest, field: str, expected_type: DataType) -> None:
    data_type = request.inputFields[field].dataType
    if data_type != expected_type:
        raise ValueError(f'Input field {field} has data type {data_type} expected {expected_type}')


def _check_data_types(request: DataFunctionRequest, field: str, expected_types: List[DataType]) -> None:
    data_type = request.inputFields[field].dataType
    if data_type not in expected_types:
        msg = ','.join(str(e) for e in expected_types)
        raise ValueError(f'Input field {field} has data type {data_type} expected one of {msg}')


def string_input_field(request: DataFunctionRequest, field: str, default_value: Optional[str] = None) -> Optional[str]:
    """
    A helper method for retrieving a string input field from a request.

    :param request: the request
    :param field: the input field
    :param default_value: a default value to return if the input field has no data
    :return:
    """
    data = _input_field_value(request, field)
    if data is None:
        return default_value
    _check_data_type(request, field, DataType.STRING)
    return data


def integer_input_field(request: DataFunctionRequest, field: str, default_value: Optional[int] = None) -> Optional[int]:
    """
    A helper method for retrieving an integer input field from a request.

    :param request: the request
    :param field: the input field
    :param default_value: a default value to return if the input field has no data
    :return:
    """
    data = _input_field_value(request, field)
    if data is None:
        return default_value
    _check_data_types(request, field, [DataType.LONG, DataType.INTEGER])
    return data


def boolean_input_field(request: DataFunctionRequest, field: str, default_value: Optional[bool] = None) -> Optional[
    bool]:
    """
    A helper method for retrieving a boolean input field from a request.

    :param request: the request
    :param field: the input field
    :param default_value: a default value to return if the input field has no data
    :return:
    """
    data = _input_field_value(request, field)
    if data is None:
        return default_value
    _check_data_type(request, field, DataType.BOOLEAN)
    return data


def double_input_field(request: DataFunctionRequest, field: str, default_value: Optional[float] = None) -> Optional[
    float]:
    """
    A helper method for retrieving a double input field from a request.

    :param request: the request
    :param field: the input field
    :param default_value: a default value to return if the input field has no data
    :return:
    """
    data = _input_field_value(request, field)
    if data is None:
        return default_value
    _check_data_types(request, field, [DataType.DOUBLE, DataType.FLOAT])
    return data


def binary_input_field(request: DataFunctionRequest, field: str, default_value: Optional[str] = None) -> Optional[str]:
    """
    A helper method for retrieving a binary input field from a request.

    :param request: the request
    :param field: the input field
    :param default_value: a default value to return if the input field has no data
    :return:
    """
    data = _input_field_value(request, field)
    if data is None:
        return default_value
    _check_data_type(request, field, DataType.BINARY)
    return data


def string_list_input_field(request: DataFunctionRequest, field: str, default_value: Optional[List[str]] = None) -> \
        Optional[List[str]]:
    """
    A helper method for retrieving a string list input field from a request.

    :param request: the request
    :param field: the input field
    :param default_value: a default value to return if the input field has no data
    :return:
    """
    data = _input_field_value(request, field)
    if data is None:
        return default_value
    _check_data_type(request, field, DataType.STRING_LIST)
    return data


def _columnLookupKey(input_field: InputField, column_id: str) -> str:
    if not input_field.limitBy:
        return column_id
    key = column_id + '_' + input_field.limitBy
    if input_field.limitByMarking:
        key = key + "_" + input_field.limitByMarking
    return key


def input_field_to_column_key(request: DataFunctionRequest, field: str) -> Optional[str]:
    """
    A helper method for retrieving the data column lookup key corresponding to a string input field that is a single column selector

    @param request: the request
    @param field: the input field
    @return: the lookup key
    """
    column_id = string_input_field(request, field)
    if not column_id:
        return None
    key = _columnLookupKey(request.inputFields[field], column_id)
    return key


def input_field_to_column_keys(request: DataFunctionRequest, field: str) -> Optional[List[str]]:
    """
    A helper method for retrieving the data columns lookup keys corresponding to a string input field that selects multiple columns

    @param request: the request
    @param field: the input field
    @return: the lookup keys
    """
    column_ids = string_list_input_field(request, field)
    if not column_ids:
        return None
    keys = [_columnLookupKey(request.inputFields[field], column_id) for column_id in column_ids]
    return keys


def input_field_to_column(request: DataFunctionRequest, field: str) -> Optional[ColumnData]:
    """
    A helper method for retrieving the data column corresponding to a string input field that is a single column selector

    @param request: the request
    @param field: the input field
    @return: the column
    """
    key = input_field_to_column_key(request, field)
    if not key:
        return None
    column = request.inputColumns[key]
    return column


def input_field_to_columns(request: DataFunctionRequest, field: str) -> Optional[List[ColumnData]]:
    """
    A helper method for retrieving the data columns corresponding to a string input field that selects multiple columns

    @param request: the request
    @param field: the input field
    @return: the columns
    """
    keys = input_field_to_column_keys(request, field)
    if not keys:
        return None
    columns = [request.inputColumns[key] for key in keys]
    return columns
