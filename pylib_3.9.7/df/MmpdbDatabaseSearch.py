import base64
import gzip
import sqlite3
from typing import List, Tuple

from mmpdblib.analysis_algorithms import TransformResult
from rdkit import Chem
from rdkit.Chem import rdFMCS, rdDepictor
from rdkit.Chem.rdFMCS import MCSResult
from rdkit.Chem.rdchem import Mol

from df.chem_helper import input_field_to_molecule, molecules_to_column
from df.data_transfer import DataFunction, DataFunctionRequest, DataFunctionResponse, string_input_field, \
    DataType, ColumnData, TableData
from ruse.rdkit.mmp import result_to_mmp_transform, transform, RadiusAndRuleGroup, build_fingerprint, MmpTransform
from ruse.rdkit.rdkit_utils import sanitize_mol, smiles_to_mol


def database_properties(database_path: str) -> List[str]:
    connection = sqlite3.connect(database_path)
    cursor = connection.cursor()
    properties = []
    for row in cursor.execute('select name from property_name'):
        properties.append(row[0])
    cursor.close()
    connection.close()
    return properties


def mmpdb_query(request: DataFunctionRequest) -> Tuple[str, Mol, str]:
    input_query_mol: Mol = input_field_to_molecule(request, 'queryData')
    # Sanitize to ensure conversion to whole structure from MOL block
    sanitize_mol(input_query_mol)
    query_smiles: Mol = Chem.MolToSmiles(input_query_mol, True)
    # MMPDB queries use smiles, but mol to smiles may be lossy, so rebuild query mol from query smiles
    query_mol = smiles_to_mol(query_smiles)
    if input_query_mol.GetNumConformers() > 0:
        mcs: MCSResult = rdFMCS.FindMCS([input_query_mol, query_mol])
        if mcs.numAtoms > 2 and mcs.numBonds > 2:
            rdDepictor.GenerateDepictionMatching2DStructure(query_mol, input_query_mol, refPatt=mcs.queryMol)
    # input_query = string_input_field(request, 'queryData')
    input_query = Chem.MolToMolBlock(query_mol)
    mol_gzip = gzip.compress(bytes(input_query, 'utf-8'))
    query_property = base64.b64encode(mol_gzip).decode('utf-8')

    return query_smiles, query_mol, query_property


def search_mmpdb(request: DataFunctionRequest, property_names: list[str],
                 mmpdb_database_file: str, output_table_name: str) -> DataFunctionResponse:
    query_smiles, query_mol, query_property = mmpdb_query(request)
    mmp_result: TransformResult = transform(query_smiles, property_names, mmpdb_database_file)
    mmp_transforms: MmpTransform = result_to_mmp_transform(query_mol, mmp_result)

    groups: List[RadiusAndRuleGroup] = [g for p in mmp_transforms.products for g in p.group_by_rule_and_radius()]

    queries: List[Mol] = []
    products: List[Mol] = []
    froms: List[Mol] = []
    tos: List[Mol] = []
    radii: List[int] = []
    rules: List[int] = []
    transform_fingerprints: List[str] = []
    fingerprints: List[str] = []
    properties: List[List[float]] = []
    for _ in range(len(property_names)):
        properties.append([])

    for g in groups:
        queries.append(g.query)
        products.append(g.product)
        froms.append(Chem.MolFromSmiles(g.from_smiles))
        tos.append(Chem.MolFromSmiles(g.to_smiles))
        radii.append(g.radius)
        rules.append(g.rule_env_id)
        transform_fingerprints.append(build_fingerprint(query_mol, g.query, 'variable'))
        fingerprints.append(build_fingerprint(query_mol, g.query))
        for p, p_list in zip(g.properties, properties):
            p_list.append(p)

    query_column = molecules_to_column(queries, 'Query', DataType.BINARY)
    query_column.properties['mmpInputQuery'] = query_property
    product_column = molecules_to_column(products, 'Product', DataType.BINARY)
    from_column = molecules_to_column(froms, 'From', DataType.BINARY)
    to_column = molecules_to_column(tos, 'To', DataType.BINARY)
    radius_column = ColumnData(name='RADIUS', dataType=DataType.LONG, values=radii)
    rule_column = ColumnData(name='Env_rule', dataType=DataType.LONG, values=rules,
                             properties={'mmpdbPath': mmpdb_database_file})
    transform_fingerprint_column = ColumnData(name='TransformFingerprint', dataType=DataType.BINARY,
                                              values=transform_fingerprints, contentType='application/octet-stream')
    fingerprint_column = ColumnData(name='Fingerprint', dataType=DataType.BINARY,
                                    values=fingerprints, contentType='application/octet-stream')
    property_columns: List[ColumnData] = []
    for name, values in zip(property_names, properties):
        column = ColumnData(name=f'Delta {name}', dataType=DataType.DOUBLE, values=values)
        property_columns.append(column)

    columns = [query_column, product_column, from_column, to_column, radius_column, rule_column, fingerprint_column,
               transform_fingerprint_column]
    columns.extend(property_columns)
    output_table = TableData(tableName=output_table_name, columns=columns)
    response = DataFunctionResponse(outputTables=[output_table])
    return response


class MmpdbDatabaseSearch(DataFunction):
    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:
        mmpdb_database_file = string_input_field(request, 'mmpdbDatabasePath')

        property_names = database_properties(mmpdb_database_file)
        return search_mmpdb(request, property_names, mmpdb_database_file, 'MMPDB database search')
