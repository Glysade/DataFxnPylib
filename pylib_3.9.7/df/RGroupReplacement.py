import re
from collections import defaultdict
from itertools import product
from json import load
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdDepictor, rdMolAlign

from df.chem_helper import column_to_molecules, input_field_to_molecule, \
    molecules_to_column
from df.data_transfer import DataFunction, DataFunctionRequest, \
    DataFunctionResponse, DataType, \
    TableData, boolean_input_field, string_input_field


def load_replacements_data() -> dict[str: tuple[list[str], list[str]]]:
    """
    Load the replacements data from the JSON file.  Returns a dict
    keyed on the SMILES of the thing to be replaced, containing a tuple
    with the lists of SMILES of first and second layer replacements.
    :return:
    """
    im_here = Path(__file__)
    substs_file = im_here.parent.parent.parent / 'Data' / 'top500_R_replacements.json'
    with open(substs_file, 'r') as f:
        raw_substs_data = load(f)

    # it's a lot more convenient in use if the outer list is turned
    # into a dict.
    substs_data = {}
    for raw_subst in raw_substs_data:
        smi = raw_subst['@SMILES']
        layer1 = []
        layer2 = []
        # a vagary of the conversion from the original XML
        # to the JSON used here is that if the layer
        # only had 1 entry, it's not put in a list, just
        # into a dict
        first_layers = raw_subst['first_layer']
        if isinstance(first_layers, dict):
            first_layers = [first_layers]
        for first_layer in first_layers:
            layer1.append(first_layer['@SMILES'])
            try:
                second_layers = first_layer['second_layer']
                if isinstance(second_layers, dict):
                    second_layers = [second_layers]
                for second_layer in second_layers:
                    layer2.append(second_layer['@SMILES'])
            except KeyError:
                # not all substituents have second layers
                pass
        substs_data[smi] = (layer1, layer2)

    return substs_data


def make_core_and_rgroups(mol: Chem.Mol, core_query: Chem.Mol) -> list[tuple[Chem.Mol, list[str]]]:
    """
    Use the core_query to do an R Group strip on mol.  That involves
    matching the core_query onto the molecule, and removing any
    substituents off it, using Chem.FragmentOnBonds.  If the
    core_query matches more than once, all possibilities are returned.
    For these purposes, we can discount symmetrical matches of the same
    atoms.
    Returns a list of the cores as a molecule, the R Groups as SMILES
    strings, with atom mappings set up so molzip can be used
    conveniently.

    :param mol:
    :param core_query:
    :return:
    """
    ret_mols = []
    matches = mol.GetSubstructMatches(core_query)
    for match in matches:
        # take a copy so temporary markings aren't propogated
        mol_cp = Chem.Mol(mol)
        bonds_to_go = []
        dummy_labels = []
        for match_at in match:
            at = mol_cp.GetAtomWithIdx(match_at)
            at.SetProp('_GL_CORE_', f'{match_at}')
            for bond in mol_cp.GetAtomWithIdx(match_at).GetBonds():
                other_atom = bond.GetOtherAtomIdx(match_at)
                if other_atom not in match:
                    btg_num = len(bonds_to_go)
                    dummy_labels.append((btg_num + 1001, btg_num + 1001))
                    bonds_to_go.append(bond.GetIdx())
        frag_mol = Chem.FragmentOnBonds(mol_cp, bonds_to_go,
                                        dummyLabels=dummy_labels)
        # the fragmentation labels the dummy atoms at the break points
        # with isotope numbers.  Convert these to atom map numbers so
        # that molzip can be used
        for atom in frag_mol.GetAtoms():
            if atom.GetIsotope() > 1000:
                atom.SetAtomMapNum(atom.GetIsotope())
                atom.SetIsotope(0)
        frags = Chem.GetMolFrags(frag_mol, asMols=True)
        frag_core = None
        r_group_smis = []
        bidentates = []
        stars = re.compile(r'(\[\*:\d+\])')
        for frag in frags:
            for atom in frag.GetAtoms():
                if atom.HasProp('_GL_CORE_'):
                    frag_core = frag
                    break
            if frag != frag_core:
                smi = Chem.MolToSmiles(frag)
                num_stars = len(stars.findall(smi))
                if num_stars > 1:
                    bidentates.append(frag)
                else:
                    r_group_smis.append(smi)

        # if we have bidentate substituents, add them back to the core
        if bidentates:
            for bid in bidentates:
                frag_core = Chem.CombineMols(frag_core, bid)
            frag_core = Chem.molzip(frag_core)
        ret_mols.append((frag_core, r_group_smis))
    return ret_mols


def make_rgroup_lookup_smi(rgroup: str) -> tuple[str, int]:
    """
    Take the R Group SMILES that is expected to have an atom mapping
    on the dummy atom, and return a new SMILES with the atom mapping
    removed and the atom mapping number.  There should be only 1
    atom mapping (we're not dealing with substituents that attached
    twice).
    :param rgroup:
    :return:
    """
    mol = Chem.MolFromSmiles(rgroup)
    atom_map_num = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            atom_map_num = atom.GetAtomMapNum()
            atom.SetAtomMapNum(0)
            break

    return Chem.MolToSmiles(mol), atom_map_num


def make_mapped_rgroups(smis: list[str], atom_map_num: int) -> list[Chem.Mol]:
    """
    Take the SMIlES strings, create a Mol for each one and add the atom
    map num to the dummy atom.
    :param smis:
    :param atom_map_num:
    :return:
    """
    ret_mols = []
    for smi in smis:
        smi = smi.replace('*', f'[*:{atom_map_num}]')
        ret_mols.append(Chem.MolFromSmiles(smi))

    return ret_mols


def make_rgroups_for_substs(rgroup_smi: str, atom_map_num: int,
                            substs_data: dict[str: tuple[list[str], list[str]]],
                            use_layer1: bool, use_layer2: bool) -> tuple[list[Chem.Mol], list[Chem.Mol]]:
    """
    Lookup the replacements for rgroup_smi in substs_data and build
    molecules from them with the atom_map_num added to the dummy atom
    so they are ready to go straight into molzip.  Returns lists of
    the layer1 and layer2 replacements as requested.  If rgroup_smi
    isn't in substs_data, just returns rgroup_smi suitably mackled so
    that no replacement is performed in the end.  That's a bit
    inefficient, but probably good enough.

    :param rgroup_smi:
    :param atom_map_num:
    :param substs_data:
    :param use_layer1:
    :param use_layer2:
    :return:
    """
    layer1_mols = []
    layer2_mols = []
    if use_layer1:
        try:
            layer1_smis = substs_data[rgroup_smi][0]
        except KeyError:
            layer1_smis = [rgroup_smi]
        layer1_mols = make_mapped_rgroups(layer1_smis, atom_map_num)
    if use_layer2:
        try:
            layer2_smis = substs_data[rgroup_smi][1]
        except KeyError:
            layer2_smis = [rgroup_smi]
        layer2_mols = make_mapped_rgroups(layer2_smis, atom_map_num)

    return layer1_mols, layer2_mols


def make_analogues(core_and_rgroups: list[tuple[Chem.Mol, list[str]]],
                   substs_data: dict[str: tuple[list[str], list[str]]],
                   use_layer1: bool, use_layer2: bool) -> list[Chem.Mol]:
    """
    Take the list of core and r groups for the molecule, and make
    analogues of each using the relevant replacement in substs_data.
    :param core_and_rgroups:
    :param substs_data:
    :param use_layer1:
    :param usr_layer2:
    :return:
    """
    analogues = []
    for core, rgroups in core_and_rgroups:
        rgroup_repls = []
        for rgroup in rgroups:
            rgroup_lookup, atom_map_num = make_rgroup_lookup_smi(rgroup)
            layer1_mols, layer2_mols = \
                make_rgroups_for_substs(rgroup_lookup, atom_map_num,
                                        substs_data, use_layer1, use_layer2)
            rgroup_repls.append(layer1_mols + layer2_mols)
        for substs in product(*rgroup_repls):
            analogue = Chem.Mol(core)
            for s in substs:
                analogue = Chem.CombineMols(analogue, s)
            analogue = Chem.molzip(analogue)
            analogues.append(analogue)
    return analogues


def align_analogue_to_parent(analogue: Chem.Mol, parent: Chem.Mol) -> None:
    """
    Use the information in atom props _GL_CORE_ to align the analogue
    so that corresponding atoms are on top of the parent.
    :param analogue:
    :param parent:
    :return:
    """
    atom_map = []
    for at in analogue.GetAtoms():
        try:
            atom_map.append((at.GetIdx(), int(at.GetProp('_GL_CORE_'))))
        except KeyError:
            pass
    rdMolAlign.AlignMol(analogue, parent, atomMap=atom_map)


def replace_rgroups(mols: list[Chem.Mol], core_query: Chem.Mol,
                    use_layer1: bool, use_layer2: bool,
                    input_column_name: str) -> TableData:
    """
    Take the list of molecules, use the core_query to define R Groups
    of it, and then create a table of analogues where the R Groups are
    replaced using Bajorath's table.  Returns a table where each row is
    a new analogue, with the first column the parent molecule.
    :param mols: list of molecules from which analogues are to be
                 derived ([Chem.mol, ]).
    :param core_query: query molecule defining the core.  Not all
                       molecules must match, but clearly only those
                       that do will produce analogues (Chem.Mol).
    :param use_layer1: Use the layer 1 replacements in the table (bool)
    :param use_layer2: Use the layer 1 replacements in the table (bool)
    :param input_column_name: (str)
    :return: The table of analogues (TableData)
    """
    substs_data = load_replacements_data()

    all_analogues = []
    analogue_parents = []
    for mol in mols:
        if mol is not None and mol and mol.HasSubstructMatch(core_query):
            rdDepictor.Compute2DCoords(mol)
            core_and_rgroups = make_core_and_rgroups(mol, core_query)
            analogues = make_analogues(core_and_rgroups, substs_data,
                                       use_layer1, use_layer2)
            all_analogues.extend(analogues)
            analogue_parents.extend([mol] * len(analogues))

    analogue_smiles = defaultdict(int)
    final_analogues = []
    final_parents = []
    for analogue, parent in zip(all_analogues, analogue_parents):
        smi = Chem.MolToSmiles(analogue)
        if not analogue_smiles[smi]:
            align_analogue_to_parent(analogue, parent)
            final_analogues.append(analogue)
            rdDepictor.Compute2DCoords(analogue)
            final_parents.append(parent)
        analogue_smiles[smi] += 1

    table_name = f'Analogues of {input_column_name}'
    parent_col = molecules_to_column(final_parents, f'Parents {input_column_name}', DataType.BINARY)
    analogue_col = molecules_to_column(final_analogues, 'Analogues', DataType.BINARY)
    return TableData(tableName=table_name, columns=[parent_col, analogue_col])


class RGroupReplacement(DataFunction):

    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:
        raise(ValueError(request))
        column_id = string_input_field(request, 'structureColumn')
        input_column = request.inputColumns[column_id]
        mols = column_to_molecules(input_column)
        core_smarts = string_input_field(request, 'coreSMARTS')
        if core_smarts:
            core_query = Chem.MolFromSmarts(core_smarts)
        else:
            core_query = input_field_to_molecule(request, 'coreSketcher')
        use_layer1 = boolean_input_field(request, 'useLayer1')
        use_layer2 = boolean_input_field(request, 'useLayer2')
        analogues_table = replace_rgroups(mols, core_query, use_layer1,
                                          use_layer2, input_column.name)
        response = DataFunctionResponse(outputTables=[analogues_table])
        return response