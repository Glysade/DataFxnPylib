from typing import List, Dict, Optional, Union, Tuple, Any

from rdkit import Chem
from rdkit.Chem.rdChemReactions import ChemicalReaction
from rdkit.Chem.rdchem import Mol

from df.chem_helper import column_to_molecules, molecules_to_column
from df.data_transfer import DataFunction, DataFunctionRequest, DataFunctionResponse, string_input_field, \
    string_list_input_field, DataType, ColumnData, TableData
from ruse.rdkit.rdkit_utils import sanitize_mol, string_to_reaction


class ReactionTableSearch(DataFunction):
    """
    Finds reactant and product pairs in the input columns using a reaction. The reaction must include atom maps.
    """

    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:

        def delta(from_row: int, to_row: int, props: List[List[Optional[float]]]) -> Tuple[
            List[Optional[float]], List[Optional[float]], List[Optional[float]]]:
            ds: list[Optional[float]] = []
            froms: list[Optional[float]] = []
            tos: list[Optional[float]] = []
            for p in props:
                p1 = p[from_row]
                p2 = p[to_row]
                if p1 is not None and p2 is not None:
                    diff = p2 - p1
                else:
                    diff = None
                froms.append(p1)
                tos.append(p2)
                ds.append(diff)
            return froms, tos, ds

        def carry_over(from_row: int, to_row: int, values: List[List[Optional[Any]]]) -> Tuple[
            List[Optional[Any]], List[Optional[Any]]]:
            froms: list[Optional[Any]] = []
            tos: list[Optional[Any]] = []
            for p in values:
                p1 = p[from_row]
                p2 = p[to_row]
                froms.append(p1)
                tos.append(p2)
            return froms, tos

        def mol_to_smiles(mol: Mol) -> Optional[str]:
            try:
                mol = Chem.RemoveHs(mol)
                smiles = Chem.MolToSmiles(mol)
                return smiles
            except:
                return None

        rxn_field = request.inputFields['reactionQuery']
        rxn_text = str(rxn_field.data)
        rxn_content_type = rxn_field.contentType
        rxn: ChemicalReaction = string_to_reaction(rxn_content_type, rxn_text)
        if rxn.GetNumProductTemplates() != 1:
            raise ValueError('Number of products in reaction is not one')
        if rxn.GetNumReactantTemplates() != 1:
            raise ValueError('Number of reactants in reaction is not one')
        reactant = rxn.GetReactants()[0]
        chiral_reactant = False
        for atom in reactant.GetAtoms():
            if atom.GetChiralTag() in [Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW]:
                chiral_reactant = True
                break

        structure_column_id = string_input_field(request, 'structureColumn')
        structures = column_to_molecules(request.inputColumns[structure_column_id])
        structure_to_row: Dict[str, int] = {mol_to_smiles(mol): r for r, mol in enumerate(structures) if mol}
        id_column_id = string_input_field(request, 'idColumn')
        id_column = request.inputColumns[id_column_id]
        ids: List[Union[str, int]] = id_column.values
        property_column_ids = string_list_input_field(request, 'propertyColumns')
        properties = [request.inputColumns[column_id].values for column_id in property_column_ids]
        property_names = ['_'.join(request.inputColumns[column_id].name.split()) for column_id in property_column_ids]
        additional_column_ids = string_list_input_field(request, 'additionalColumns')
        if additional_column_ids:
            additional_values = [request.inputColumns[column_id].values for column_id in additional_column_ids]
            additional_names = ['_'.join(request.inputColumns[column_id].name.split()) for column_id in
                                additional_column_ids]
            additional_columns = [request.inputColumns[column_id] for column_id in additional_column_ids]

        lhs: List[Mol] = []
        rhs: List[Mol] = []
        from_ids: List[Union[str, int]] = []
        to_ids: List[Union[str, int]] = []
        deltas: List[List[Optional[float]]] = []
        from_properties: List[List[Optional[float]]] = []
        to_properties: List[List[Optional[float]]] = []
        from_values: List[List[Optional[Any]]] = []
        to_values: List[List[Optional[Any]]] = []
        for _ in range(len(property_names)):
            deltas.append([])
            from_properties.append([])
            to_properties.append([])
        if additional_column_ids:
            for _ in range(len(additional_columns)):
                from_values.append([])
                to_values.append([])

        for row, structure in enumerate(structures):
            if not structure:
                continue
            from_id = ids[row]
            if from_id is None:
                continue
            product_list = None
            try:
                skip_run = False
                if chiral_reactant:
                    if not structure.HasSubstructMatch(reactant, useChirality=True):
                        skip_run = True
                if not skip_run:
                    product_list = rxn.RunReactants([structure])
                else:
                    product_list = []
            except:
                if not product_list:
                    continue
            for products in product_list:
                assert len(products) == 1
                product: Mol = products[0]
                sanitize_mol(product)
                try:
                    product = Chem.RemoveHs(product)
                    product_smiles = Chem.MolToSmiles(product)
                except:
                    break
                if product_smiles in structure_to_row:
                    other_row = structure_to_row[product_smiles]
                    to_id = ids[other_row]
                    if to_id is None:
                        continue
                    lhs.append(structure)
                    rhs.append(structures[other_row])
                    from_ids.append(from_id)
                    to_ids.append(to_id)
                    f, t, d = delta(row, other_row, properties)
                    for i, v in enumerate(d):
                        deltas[i].append(v)
                    for i, v in enumerate(f):
                        from_properties[i].append(v)
                    for i, v in enumerate(t):
                        to_properties[i].append(v)
                    if additional_column_ids:
                        f, t = carry_over(row, other_row, additional_values)
                        for i, v in enumerate(f):
                            from_values[i].append(v)
                        for i, v in enumerate(t):
                            to_values[i].append(v)

        lhs_column = molecules_to_column(lhs, 'LHS', DataType.BINARY)
        rhs_column = molecules_to_column(rhs, 'RHS', DataType.BINARY)
        lhs_id_column = ColumnData(name=f'LHS {id_column.name}', dataType=id_column.dataType, values=from_ids)
        rhs_id_column = ColumnData(name=f'RHS {id_column.name}', dataType=id_column.dataType, values=to_ids)
        columns = [lhs_column, rhs_column, lhs_id_column, rhs_id_column]
        for index, name in enumerate(property_names):
            column = ColumnData(name=f'LHS {name}', dataType=DataType.DOUBLE, values=from_properties[index])
            columns.append(column)
            column = ColumnData(name=f'RHS {name}', dataType=DataType.DOUBLE, values=to_properties[index])
            columns.append(column)
            column = ColumnData(name=f'Delta {name}', dataType=DataType.DOUBLE, values=deltas[index])
            columns.append(column)
        if additional_column_ids:
            for index, name in enumerate(additional_names):
                additional_column = additional_columns[index]
                column = ColumnData(name=f'LHS {name}', dataType=additional_column.dataType,
                                    values=from_values[index], contentType=additional_column.contentType)
                columns.append(column)
                column = ColumnData(name=f'RHS {name}', dataType=additional_column.dataType,
                                    values=to_values[index], contentType=additional_column.contentType)
                columns.append(column)

        output_table = TableData(tableName='Pair Finder Results', columns=columns)
        response = DataFunctionResponse(outputTables=[output_table])
        return response
