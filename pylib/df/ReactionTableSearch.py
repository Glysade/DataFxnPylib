from typing import List, Dict, Optional

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import Mol

from df.chem_helper import column_to_molecules, molecules_to_column
from df.data_transfer import DataFunction, DataFunctionRequest, DataFunctionResponse, string_input_field, \
    string_list_input_field, DataType, ColumnData, TableData
from ruse.rdkit.rdkit_utils import sanitize_mol


class ReactionTableSearch(DataFunction):
    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:

        def delta(from_row: int, to_row: int, props: List[List[Optional[float]]]) -> List[Optional[float]]:
            ds: list[Optional[float]] = []
            for p in props:
                p1 = p[from_row]
                p2 = p[to_row]
                if p1 is not None and p2 is not None:
                    diff = p2 - p1
                else:
                    diff = None
                ds.append(diff)
            return ds

        rxn_field = request.inputFields['reactionQuery']
        rxn_text = str(rxn_field.data)
        rxn_content_type = rxn_field.contentType
        if rxn_content_type in ['chemical/x-smiles', 'chemical/x-smarts']:
            if '>>' not in rxn_text:
                raise ValueError(f'Input {rxn_text} is not a reaction smarts')
            rxn = AllChem.ReactionFromSmarts(rxn_text)
        elif rxn_content_type in ['chemical/x-mdl-molfile', 'chemical/x-mdl-molfile-v3000']:
            if '$RXN' not in rxn_text:
                raise ValueError(f'Input is not a RXN block: {rxn_text}')
            rxn = AllChem.ReactionFromRxnBlock(rxn_text)
        else:
            raise ValueError(f'Unable to convert content type {rxn_content_type} to reaction')
        if rxn.GetNumProductTemplates() != 1:
            raise ValueError('Number of products in reaction is not one')
        if rxn.GetNumReactantTemplates() != 1:
            raise ValueError('Number of reactants in reaction is not one')

        structure_column_id = string_input_field(request, 'structureColumn')
        structures = column_to_molecules(request.inputColumns[structure_column_id])
        property_column_ids = string_list_input_field(request, 'propertyFields')
        structure_to_row: Dict[str, int] = {Chem.MolToSmiles(mol): r for r, mol in enumerate(structures) if mol}
        properties = [request.inputColumns[column_id].values for column_id in property_column_ids]
        property_names = ['_'.join(request.inputColumns[column_id].name.split()) for column_id in property_column_ids]

        lhs: List[Mol] = []
        rhs: List[Mol] = []
        from_rows: List[int] = []
        to_rows: List[int] = []
        deltas: List[List[Optional[float]]] = []
        for _ in range(len(property_names)):
            deltas.append([])

        for row, structure in enumerate(structures):
            if not structure:
                continue
            product_list = rxn.RunReactants([structure])
            for products in product_list:
                assert len(products) == 1
                product: Mol = products[0]
                sanitize_mol(product)
                product_smiles = Chem.MolToSmiles(product)
                if product_smiles in structure_to_row:
                    other_row = structure_to_row[product_smiles]
                    lhs.append(structure)
                    rhs.append(structures[other_row])
                    from_rows.append(row)
                    to_rows.append(other_row)
                    d = delta(row, other_row, properties)
                    for i, v in enumerate(d):
                        deltas[i].append(v)

        lhs_column = molecules_to_column(lhs, 'LHS', DataType.BINARY)
        rhs_column = molecules_to_column(rhs, 'RHS', DataType.BINARY)
        lhs_row_column = ColumnData(name='LHS row', dataType=DataType.INTEGER, values=from_rows)
        rhs_row_column = ColumnData(name='RHS row', dataType=DataType.INTEGER, values=to_rows)
        columns = [lhs_column, rhs_column, lhs_row_column, rhs_row_column]
        for name, d in zip(property_names, deltas):
            column = ColumnData(name=f'Delta {name}', dataType=DataType.FLOAT, values=d)
            columns.append(column)
        output_table = TableData(tableName='Reaction table search', columns=columns)
        response = DataFunctionResponse(outputTables=[output_table])
        return response
