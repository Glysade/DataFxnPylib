"""
==========================
chem_data_table_helper.py
==========================

Copyright (C) 2017-2022 Glysade, LLC

Functions for adding molecules to and extracting molecules from :class:`ruse.util.data_table.DataTable` objects

- Encoding and decoding data table cells containing molecules
- Creating data tables and data table columns from lists of molecules

Uses RDKit for chemical structure handling
"""

import base64
import gzip
import uuid

from rdkit import Chem

from rdkit.Chem.rdchem import Mol
from typing import Iterable, List

from ruse.rdkit.rdkit_utils import type_to_format, mol_to_string, string_to_mol, is_three_dimensional
from ruse.rdkit.rgroup import Rgroup, RgroupDecomposer
from ruse.util.data_table import DataTable


def _compress_type(type: str):
    """
    Should this type be compressed and Base64 encoded

    :param type: chemical mime type
    :return: True if this type should be compressed
    """

    return type not in ['chemical/x-smiles', 'chemical/x-daylight-smiles', 'chemical/x-smarts',
                        'chemical/x-daylight-smarts']


def _encode_mol_cell(type: str, value: str) -> str:
    """
    Encodes structure information for inclusion in a :class:`ruse.util.data_table.DataTable` cell.  Smiles strings are
     added as is, while MDL formats will be gzipped and BASE64 encoded

    :param type: chemical mime type
    :param value: structure
    :return: encoded structure
    """

    if _compress_type(type):
        data = gzip.compress(value.encode('utf-8'))
        return base64.b64encode(data).decode('utf-8')
    else:
        return value


def _decode_mol_cell(type: str, data_type: str, contents: str) -> str:
    """
    Decodes a :class:`ruse.util.data_table.DataTable` cell containing an encoded chemical structure of the
    supplied mime type.

    :param type:  chemical mime type of cell contents
    :param contents:  cell contents
    :return: decoded chemical structure
    """
    if data_type == 'binary':
        binary_data = base64.b64decode(contents)
        value = gzip.decompress(binary_data)
        for charset in ['ascii', 'utf-8', 'latin-1']:
            try:
                value = value.decode(charset)
                return value
            except UnicodeDecodeError:
                pass
        raise UnicodeDecodeError("Unable to decode mol cell")
    else:
        return contents


# This signature that is compatible with Typing breaks Sphinx.  We can't just import Rocs
# as we get circular import errors
# def create_data_table_from_rocs_search(rocs: 'ruse.services.rocs_service.Rocs') -> DataTable:
# so remove argument typing for now
def create_data_table_from_rocs_search(rocs) -> DataTable:
    """
    Create a new :class:`ruse.util.data_table.DataTable` from the results of an OpenEye ROCS search.  The new table has
    columns for the query and target structures, rank, combo score, shape and color tanimoto.

    :param rocs: :class:`ruse.services.rocs_service.Rocs` object containing ROCS search result
    :return: the new data table
    """

    data_table = DataTable()
    content_type = 'chemical/x-mdl-molfile'
    format = type_to_format(content_type)

    with rocs.queries() as queries:
        for idx, query in enumerate(queries):
            query_cell = _encode_mol_cell(content_type, mol_to_string(format, query))
            for hit, hit_information in rocs.hits_and_hit_information_for_query(idx):
                hit_cell = _encode_mol_cell(content_type, mol_to_string(format, hit))
                data_table.data.append([
                    query_cell, hit_cell, hit_information.rank, hit_information.combo_tanimoto,
                    hit_information.shape_tanimoto, hit_information.color_tanimoto
                ])

    data_table.columns = [
        data_table.column_definition("Query", "binary", content_type, {'Dimension': '3'}),
        data_table.column_definition("Hit", "binary", content_type, {'Dimension': '3'}),
        data_table.column_definition("Rank", "integer"),
        data_table.column_definition("Tanimoto Combo", "float"),
        data_table.column_definition("Shape Tanimoto", "float"),
        data_table.column_definition("Color Tanimoto", "float")
    ]

    return data_table


def data_table_column_to_mols(data_table: DataTable, column_index: int, include_empty: bool = False) -> Iterable[Mol]:
    """
    Convert a column in a data table containing chemical structures to a generator of RDKit molecules

    :param data_table: the data table as :class:`ruse.util.data_table.DataTable`
    :param column_index: index of the column containing chemical structures
    :param include_empty: set True to include empty cells as None in the iterator
    :return: A generator of :class:`rdkit.Chem.rdchemRdKit.Mol` molecules
    """

    column_info = data_table.columns[column_index]
    content_type = column_info['properties']['ContentType']
    data_type = column_info['dataType']
    for row_idx, row in enumerate(data_table.data):
        if row[column_index]:
            mol_string = _decode_mol_cell(content_type, data_type, row[column_index])
            mol = string_to_mol(type_to_format(content_type), mol_string)
            if mol:
                mol.SetIntProp("RowIdx", row_idx)
                yield mol
            else:
                print("Failed to convert {} content type {} to mol!".format(mol_string, content_type))
                if include_empty:
                    yield None
        elif include_empty:
            yield None


def add_mols_as_data_table_column(data_table: DataTable, mols: Iterable[Mol],
                                  content_type: str = "chemical/x-mdl-molfile",
                                  title: str = "Molecular conformation", three_dimensional: bool = None) -> None:
    """
    Add the molecules as a new column to a data table.  They should have an SD tag containing the row index

    :param data_table: The data table
    :param mols: A iterator of :class:`rdkit.Chem.rdchemRdKit.Mol` molecules
    :param content_type: The content type to store the molecules as
    :param title: Optional column title
    :param three_dimensional: set true if the molecules are 3D.  If None :func:`ruse.rdkit.rdkit_util.is_three_dimensional` is used to determine 3D nature
    """

    new_values = [None] * len(data_table.data)
    for mol in mols:
        if not mol.HasProp("RowIdx"):
            raise ValueError("No row index tag for input mol")
        row_idx = mol.GetIntProp("RowIdx")
        if three_dimensional is None:
            three_dimensional = is_three_dimensional(mol)
        mol_string = mol_to_string(type_to_format(content_type), mol)
        new_values[row_idx] = _encode_mol_cell(content_type, mol_string)

    for row_no, row in enumerate(data_table.data):
        row.append(new_values[row_no])

    properties = {'Dimension': '3'} if three_dimensional else {}
    data_type = "binary" if _compress_type(content_type) else "string"
    data_table.add_column(title, data_type, properties, content_type=content_type, fill_rows=False)


def create_data_table_from_mols(mols: Iterable[Mol], content_type: str = "chemical/x-smiles",
                                title: str = "Molecule") -> DataTable:
    """
    Create a new table with a single column containing the input molecules

    :param mols: A iterator of :class:`rdkit.Chem.rdchemRdKit.Mol` molecules
    :param content_type: The chemical mime type specifying the format to be used to store the molecules
    :param title: title to use for new column
    :return: The new data table as :class:`ruse.util.data_table.DataTable`
    """

    data_table = DataTable()
    data_table.data = [[_encode_mol_cell(content_type, mol_to_string(type_to_format(content_type), mol))] for mol in
                       mols]
    data_table.columns = [
        DataTable.column_definition(title, "binary" if _compress_type(content_type) else "string",
                                    content_type)]
    return data_table


def data_table_column_to_local_files(data_table: DataTable, column_index: int, prefix: str = None,
                                     num_mols_per_file: int = 1) -> List[str]:
    """
    Writes all structures in a column in batches to local MOL files.  The files will be named <prefix>_<starting_row_index>.mol

    :param data_table:  The data table to extract the molecules from
    :param column_index:  Column number of molecular column
    :param prefix: File name prefix (UUID will be used if this is not set)
    :param num_mols_per_file: Batch size, defaults to 1
    :return: A list of the file names created
    """

    if not prefix:
        prefix = str(uuid.uuid4())
    suffix = 'mol' if num_mols_per_file == 1 else 'sdf'
    files = []
    fh = None
    for mol in data_table_column_to_mols(data_table, column_index=column_index, include_empty=False):
        row_idx = mol.GetIntProp("RowIdx")
        if not fh or row_idx % num_mols_per_file == 0:
            if fh:
                fh.close()
            file_out = "{}_{}.{}".format(prefix, row_idx, suffix)
            files.append(file_out)
            fh = Chem.SDWriter(file_out)
        fh.write(mol)

    if fh:
        fh.close()

    return files


def create_data_table_from_rgroup_decomposition(decomp: RgroupDecomposer,
                                                core_format: str = 'chemical/x-mdl-molfile') -> DataTable:
    """
    Creates a data table from the results of an R Group decomposition.  Multiple structure columns are created in the
    output table: one for the structures, one for matching cores and a column for each rgroup

    :param decomp: results of core decomposition of a set of molecules as :class:`ruse.rdkit.rgroup.RgroupDecomposer`
    :param core_format: the chemical mime type to encode the cores in, should be chemical/x-mdl-molfile or chemical/x-smarts.  Defaults to chemical/x-mdl-molfile
    :return: the data table
    """
    mol_columns = decomp.to_molecule_grid(False, True)
    for column_no, mols in enumerate(mol_columns):
        if column_no == 0:
            data_table = create_data_table_from_mols(mols, 'chemical/x-smiles', 'Molecule')
        elif column_no == 1:
            core_table = create_data_table_from_mols(mols, core_format, 'Core')
            data_table.append_columns(core_table)
        else:
            for row_idx, mol in enumerate(mols):
                if mol:
                    mol.SetIntProp("RowIdx", row_idx)
            mols_to_add = [m for m in mols if m]
            title = "Molecule" if column_no == 1 else 'R{}'.format(column_no - 1)
            add_mols_as_data_table_column(data_table, mols_to_add, 'chemical/x-smiles', title, three_dimensional=False)
    return data_table


def merge_data_table_from_rgroup_decomposition(table: DataTable, decomp: RgroupDecomposer,
                                               core_format: str = 'chemical/x-mdl-molfile') -> DataTable:
    """
    Merges the results of an R Group decomposition with an existing data table.  Multiple structure columns are created in the
    output table: one for matching cores and a column for each rgroup

    :param table: original input data table
    :param decomp: results of core decomposition of a set of molecules as :class:`ruse.rdkit.rgroup.RgroupDecomposer`
    :param core_format: the chemical mime type to encode the cores in, should be chemical/x-mdl-molfile or chemical/x-smarts.  Defaults to chemical/x-mdl-molfile
    :return: the data table
    """

    assert (len(table.data) == len(decomp.decomposition))
    core_column_no = len(table.columns)
    smi_props = {'ContentType': 'chemical/x-smiles'}
    # for now force core format to smiles, as that is what SAR map expects
    core_format = 'chemical/x-smiles'
    table.add_column('Core', 'string', {'ContentType': core_format})
    for group_no in decomp.r_group_numbers:
        table.add_column('R{}'.format(group_no), 'string', smi_props)

    new_data = list()
    for table_row, mol_decomps in zip(table.data, decomp.decomposition):
        if not mol_decomps:
            new_data.append(table_row)
        else:
            for mol_decomp in mol_decomps:
                new_row = list(table_row)
                new_row[core_column_no] = _encode_mol_cell(core_format,
                                                           mol_to_string(type_to_format(core_format),
                                                                         mol_decomp.core.display_mol()))
                new_row[core_column_no + 1:] = [mol_to_string(type_to_format('chemical/x-smiles'), s.mol) if s else None
                                                for s in mol_decomp.sidechains]
                new_data.append(new_row)
    table.data = new_data
    return table


def merge_data_table_from_rdkit_rgroup_decomposition(table: DataTable, decomp: Rgroup,
                                               core_format: str = 'chemical/x-mdl-molfile') -> DataTable:
    """
    Merges the results of an R Group decomposition with an existing data table.  Multiple structure columns are created in the
    output table: one for matching cores and a column for each rgroup

    :param table: original input data table
    :param decomp: results of core decomposition of a set of molecules as :class:`ruse.rdkit.rgroup.RgroupDecomposer`
    :param core_format: the chemical mime type to encode the cores in, should be chemical/x-mdl-molfile or chemical/x-smarts.  Defaults to chemical/x-mdl-molfile
    :return: the data table
    """

    decomposition_rows = decomp.to_molecule_grid()
    assert (len(table.data) == len(decomposition_rows))
    core_column_no = len(table.columns)
    smi_props = {'ContentType': 'chemical/x-smiles'}
    # for now force core format to smiles, as that is what SAR map expects
    core_format = 'chemical/x-smiles'
    table.add_column('Core', 'string', {'ContentType': core_format})
    for group_no in decomp.rgroup_labels()[1:]:
        table.add_column('{}'.format(group_no), 'string', smi_props)

    new_data = list()
    for table_row, decomp_row in zip(table.data, decomposition_rows):
        new_row = list(table_row)
        new_row[core_column_no] = _encode_mol_cell(core_format,
                                                   mol_to_string(type_to_format(core_format),
                                                                 decomp_row[0])) if decomp_row[0] else None;
        new_row[core_column_no + 1:] = [mol_to_string(type_to_format('chemical/x-smiles'), s.mol) if s else None
                                        for s in decomp_row[1:]]
        new_data.append(new_row)
    table.data = new_data
    return table
