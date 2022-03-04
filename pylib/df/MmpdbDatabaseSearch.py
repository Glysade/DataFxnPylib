import hashlib
import json
import multiprocessing
import os.path
import re
import tempfile
from datetime import datetime
from typing import List, Final

try:
    from mmpdblib import commandline as mmp
except ImportError:
    from mmpdblib import cli as mmp
from mmpdblib.analysis_algorithms import TransformResult
from pydantic import BaseModel
from rdkit import Chem
from rdkit.Chem import rdFMCS, rdDepictor
from rdkit.Chem.rdFMCS import MCSResult
from rdkit.Chem.rdchem import Mol

from df.chem_helper import column_to_molecules, input_field_to_molecule, molecules_to_column
from df.data_transfer import DataFunction, DataFunctionRequest, DataFunctionResponse, string_input_field, \
    string_list_input_field, DataType, ColumnData, TableData
from ruse.rdkit.mmp import result_to_mmp_transform, transform, RadiusAndRuleGroup, build_fingerprint, MmpTransform
from ruse.rdkit.rdkit_utils import sanitize_mol, smiles_to_mol

VERSION: Final[str] = '0.0.2'


class MmpdbProperties(BaseModel):
    date_created: datetime
    checksum: str
    property_names: List[str]
    database_file: str
    version: str = VERSION


def build_mmp_database(mmpdb_dir: str, request: DataFunctionRequest) -> MmpdbProperties:
    fragments_file = os.path.join(mmpdb_dir, 'mmpdb.fragments.gz')
    database_file = os.path.join(mmpdb_dir, 'mmpdb.mmpdb')
    property_file = os.path.join(mmpdb_dir, 'mmpdb.props')
    smiles_file = os.path.join(mmpdb_dir, 'mmpdb.smiles')
    settings_file = os.path.join(mmpdb_dir, 'mmpdb.settings')

    structure_column_id = string_input_field(request, 'structureColumn')
    structures = column_to_molecules(request.inputColumns[structure_column_id])
    smiles = [Chem.MolToSmiles(m) if m else None for m in structures]
    property_column_ids = string_list_input_field(request, 'propertyFields')
    properties = [request.inputColumns[column_id].values for column_id in property_column_ids]
    property_names = ['_'.join(request.inputColumns[column_id].name.split()) for column_id in property_column_ids]

    input_data = {'smiles': smiles, 'properties': properties}
    input_data_str = json.dumps(input_data)
    input_data_checksum = hashlib.md5(input_data_str.encode()).hexdigest()

    if not os.path.exists(mmpdb_dir):
        os.makedirs(mmpdb_dir)

    if os.path.exists(settings_file):
        with open(settings_file, 'r') as fh:
            settings_json = fh.read()

        settings = MmpdbProperties.parse_raw(settings_json)
        if settings.checksum == input_data_checksum and settings.version == VERSION:
            return settings
        else:
            os.remove(settings_file)

    if os.path.exists(database_file):
        os.remove(database_file)

    with open(property_file, 'w') as prop_fh, open(smiles_file, 'w') as smi_fh:
        property_name_str = ' '.join(property_names)
        prop_fh.write(f'ID {property_name_str}\n')
        for compound_num in range(0, len(smiles)):
            smi = smiles[compound_num]
            if not smi:
                continue
            props = [p[compound_num] for p in properties]
            if all(p is None for p in props):
                continue
            prop_str = ' '.join('*' if p is None else str(p) for p in props)
            smi_fh.write(f'{smi} {compound_num}\n')
            prop_fh.write(f'{compound_num} {prop_str}\n')

    # build fragments file
    n_proc = multiprocessing.cpu_count()
    # args = ['fragment', '--in', 'smi', '--out', 'fragments.gz', '--num-jobs', str(n_proc), '--output',
    args = ['fragment', '--in', 'smi', '--num-jobs', str(n_proc), '--output',
            fragments_file, '--delimiter',
            'space', smiles_file]
    try:
        mmp.main(args=args, standalone_mode=False)
    except TypeError:
        mmp.main(args)

    # index fragments
    args = ['index', fragments_file, '-o', database_file]
    try:
        mmp.main(args=args, standalone_mode=False)
    except TypeError:
        mmp.main(args)

    # add properties to database
    args = ['loadprops', '-p', property_file, database_file]
    try:
        mmp.main(args=args, standalone_mode=False)
    except TypeError:
        mmp.main(args)

    os.remove(smiles_file)
    os.remove(property_file)
    os.remove(fragments_file)

    settings = MmpdbProperties(checksum=input_data_checksum, date_created=datetime.now(), property_names=property_names,
                               database_file=database_file)
    with open(settings_file, 'w') as fh:
        fh.write(settings.json())

    return settings


class MmpdbDatabaseSearch(DataFunction):
    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:
        structure_column_id = string_input_field(request, 'structureColumn')
        mmpdb_dir = re.sub('[^A-Za-z0-9 ]+', '', structure_column_id)
        mmpdb_dir = re.sub(' ', '_', mmpdb_dir)
        mmpdb_dir = os.path.join(tempfile.gettempdir(), 'ChemChartsData', f'mmpdb_{mmpdb_dir}')

        settings = build_mmp_database(mmpdb_dir, request)
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
        mmp_result: TransformResult = transform(query_smiles, settings.property_names, settings.database_file)
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
        for _ in range(len(settings.property_names)):
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
        product_column = molecules_to_column(products, 'Product', DataType.BINARY)
        from_column = molecules_to_column(froms, 'From', DataType.BINARY)
        to_column = molecules_to_column(tos, 'To', DataType.BINARY)
        radius_column = ColumnData(name='RADIUS', dataType=DataType.INTEGER, values=radii)
        rule_column = ColumnData(name='Env_rule', dataType=DataType.INTEGER, values=rules,
                                 properties={'mmpdbPath': settings.database_file})
        transform_fingerprint_column = ColumnData(name='TransformFingerprint', dataType=DataType.BINARY,
                                                  values=transform_fingerprints, contentType='application/octet-stream')
        fingerprint_column = ColumnData(name='Fingerprint', dataType=DataType.BINARY,
                                        values=fingerprints, contentType='application/octet-stream')
        property_columns: List[ColumnData] = []
        for name, values in zip(settings.property_names, properties):
            column = ColumnData(name=f'Delta {name}', dataType=DataType.FLOAT, values=values)
            property_columns.append(column)

        columns = [query_column, product_column, from_column, to_column, radius_column, rule_column, fingerprint_column,
                   transform_fingerprint_column]
        columns.extend(property_columns)
        output_table = TableData(tableName='MMPDB table search', columns=columns)
        response = DataFunctionResponse(outputTables=[output_table])
        return response
