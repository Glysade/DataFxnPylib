"""
============
rgroup.py
============

Contains classes for performing r group decomposition
"""

# Note: RDKit and attachment labels
# using SetIsotope to set r-group number you get labels like 1* and smiles like '[1*]*1CCC([2*])C1[3*]'
# using SetAtomMapNumber you get labels like *:1 and smiles like 'C1CC([*:2])C([*:3])*1[*:1]'
# using both you get labels like 1*:1 and smiles like 'C1CC([2*:2])C([3*:3])*1[1*:1]'
# RDkit uses Isotope for both rgroup numbers and core atom numbers for dummy atoms in fragment decomposition
# RDKit uses AtomMapNumber for reaction mapping
# In this code I hijack AtomMapNumber to track rgroup numbers in sidechains during fragment decomposition

from typing import List, NamedTuple, Optional, Tuple, Set, Dict, Any

import functools
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, Mol, RWMol, Atom, SanitizeFlags
from rdkit.Chem.Descriptors import HeavyAtomMolWt

from ruse.rdkit.rdkit_utils import remove_atom_mappings, sanitize_mol
from ruse.util.frozen import Frozen

try:
    from rdkit.Chem.rdRGroupDecomposition import RGroupDecomposition, RGroupDecompositionParameters, RGroupMatching, \
    RGroupCoreAlignment, RGroupScore, RGroupLabelling
except (NameError, ImportError):
    pass


class AttachmentPoint(NamedTuple):
    """
    A named tuple that stores information about core R group attachment points

    Attributes:

     - core_index: the atom index of the attachment point on the core

     - group_num: the associated R group number

    """
    core_index: int
    group_num: int
    required: bool


class Attachment(NamedTuple):
    """
    A named tuple that describes an R-group and associated attachment point

    Attributes:

        - attachment_point: the location of this R group on the core as an :class:`AttachmentPoint` object

        - mol: the R-group as a :class:`rdkit.Chem.Mol` object

    """
    attachment_point: AttachmentPoint
    mol: Mol

    def __eq__(self, other: 'Attachment'):
        """
        Returns true if two attachment points are identical (same structures and attachment points)
        :param other:
        :return:
        """
        if other is None:
            return False
        if self.attachment_point != other.attachment_point:
            return False
        return self.same_sidechain(other)

    def __ne__(self, other) -> bool:
        return not self.__eq__(other)

    def __repr__(self):
        return 'R{} index {}, {}'.format(self.attachment_point.group_num, self.attachment_point.core_index,
                                         Chem.MolToSmiles(self.mol))

    def heavy_sidechain(self) -> None:
        """
        :return: True if this sidechain is not a hydrogen
        """
        return self.mol and HeavyAtomMolWt(self.mol) > 0

    def relabel_sidechain(self) -> None:
        """
        Relabel attachment atom with R group number- should only be done once all matching is completed
        """
        # this should only be done after all matching is finished, otherwise this attachment could incorrectly
        # match an additional core
        match = False
        for atom in self.mol.GetAtoms():
            if atom.GetAtomicNum() == 0 and atom.GetIsotope() == self.attachment_point.core_index:
                # the rdkit model uses isotope for smarts labels, but this does, not work for smiles
                # use smiles reaction class labels instead
                atom.SetAtomMapNum(self.attachment_point.group_num)
                match = True
        assert (match)

    def same_sidechain(self, other: 'Attachment', strip_atom_map_nums: bool = False) -> bool:
        if strip_atom_map_nums:
            strip_mol = remove_atom_mappings(self.mol)
            other_strip_mol = remove_atom_mappings(other.mol)
            return Chem.MolToSmiles(strip_mol, True) == Chem.MolToSmiles(other_strip_mol, True)
        else:
            return Chem.MolToSmiles(self.mol, True) == Chem.MolToSmiles(other.mol, True)


def is_r_group_atom(atom: Atom) -> Optional[int]:
    """
    Returns true if this is an R group atom.

    :param atom: An rdkit atom
    :return: The rgroup number, or none if the atom does not define an rgroup
    """
    # in an rdkit query mol a wildcard atom has atomic number 0 and the r group number is stored in the isotope property
    if atom.GetAtomicNum() == 0 and atom.GetIsotope() > 0:
        return atom.GetIsotope()
    # alternatively for SMARTS queries we can use reaction mapping indices on wildcard atoms
    elif atom.GetAtomicNum() == 0 and atom.GetAtomMapNum() > 0:
        return atom.GetAtomMapNum()
    elif atom.HasProp("rgroup") and atom.GetIntProp("rgroup") > 0:
        return atom.GetIntProp("rgroup")
    return None


class Core(Frozen):
    """
    A class to represent a core query for decomposition

    If the core contains explicit R groups they are used to define attachment points- otherwise any atom in the core
    may be used as an attachment point

    Attributes:

        - index: the index of the core in a multi-core decomposition

        - mol: a class:`rdkit.Chem.Mol` molecule that contains the original core definition

        - attachment_points: a list of R-group attachment points or :class:`AttachmentPoint` objects

        - matching_mol: a class:`rdkit.Chem.Mol` molecule that contains the core matcher.  This differs from the mol attribute as it will have attachment points removed
    """

    def __init__(self, index: int, mol: Mol, allow_rgroups_at_any_free_valance=False):
        """
        Constructor. Finds attachment points and constructs core matching molecule

        :param index: core index in multi-core decomposition
        :param mol: a class:`rdkit.Chem.Mol` molecule that defined the core
        """

        self.index = index
        self.mol = mol
        self.attachment_points = None  # type: List[AttachmentPoint]

        # find attachment points:

        mol = self.mol
        matching_mol = RWMol(mol)
        self.matching_mol = matching_mol  # type: Mol
        r_group_atoms = [a for a in matching_mol.GetAtoms() if is_r_group_atom(a)]
        group_numbers = []
        self.attachment_points = []

        if r_group_atoms:
            # explict attachment points in the core query- allow to attach only to these
            atoms_and_groups = []

            for r_group_atom in r_group_atoms:
                group_num = is_r_group_atom(r_group_atom)
                if group_num in group_numbers:
                    raise ValueError('Group number {} appears more than once in core'.format(group_num))
                group_numbers.append(group_num)
                assert (group_num > 0)

                neighbours = r_group_atom.GetNeighbors()
                if r_group_atom.GetAtomicNum() == 0 and len(neighbours) == 1:
                    # if it looks like the labelled R-Group is not part of the core remove it
                    neighbour = neighbours[0]
                    matching_mol.RemoveAtom(r_group_atom.GetIdx())
                    atoms_and_groups.append((neighbour, group_num))
                else:
                    r_group_atom.ClearProp('rgroup')
                    r_group_atom.SetAtomMapNum(0)
                    atoms_and_groups.append((r_group_atom, group_num))

            # add reaction class labels to original core for depiction
            for atom in [a for a in mol.GetAtoms() if is_r_group_atom(a)]:
                group_num = is_r_group_atom(atom)
                atom.SetAtomMapNum(group_num)
                atom.SetIsotope(0)

            self.attachment_points = [AttachmentPoint(a.GetIdx(), r, True) for a, r in atoms_and_groups]

        # once we have looked at symmetry we can merge explicit and implicit attachment points
        if not r_group_atoms or allow_rgroups_at_any_free_valance:
            # find implicit attachment points
            # sanitize_mol(matching_mol)
            matching_mol.UpdatePropertyCache(False)
            group_num = 0
            required_attachment_points = self.attachment_points.copy()

            periodic_table = Chem.GetPeriodicTable()
            for index, atom in enumerate(matching_mol.GetAtoms()):
                if any(a for a in required_attachment_points if a.core_index == index):
                    continue
                if atom.GetSymbol() == '*' or atom.HasQuery():
                    n_attachments = 1
                else:
                    bonds = atom.GetBonds()
                    has_query_bond = any(b.HasQuery() for b in bonds)
                    if has_query_bond:
                        # implicit valence is wrong if atom has query bonds
                        # assign a valence of 1 to query bonds- this is wrong if the query is something like
                        # double or aromatic- but in that case we will have too many attachment points and they
                        # should just be filtered out
                        # At the moment we restrict the number of attachments per atom to 1 anyway
                        bond_valances = [1 if b.HasQuery() else b.GetValenceContrib(atom) for b in bonds]
                        bond_valance = int(round(sum(bond_valances)))
                        n_attachments = periodic_table.GetDefaultValence(atom.GetSymbol()) - bond_valance
                    else:
                        n_attachments = atom.GetImplicitValence() + atom.GetNumExplicitHs()
                if n_attachments > 0:
                    # For now only allow one attachment per core atom- may wish to revisit this or make an option.
                    # would need to do more with the sidechain matching to make this happen
                    n_attachments = 1
                    atom.SetNumExplicitHs(0)
                    for _ in range(n_attachments):
                        group_num += 1
                        while group_num in group_numbers:
                            group_num += 1
                        self.attachment_points.append(AttachmentPoint(index, group_num, False))

        self._display_mol = None
        self.attachment_points.sort(key=lambda a: a.group_num)
        self.user_group_numbers = group_numbers

    def display_mol(self):
        if self._display_mol:
            return self._display_mol
        display_mol = Chem.RWMol(self.matching_mol)
        self._display_mol = display_mol

        for attachment_point in self.attachment_points:
            if attachment_point:
                atom = display_mol.GetAtomWithIdx(attachment_point.core_index)
                r_group_atom = Atom(0)
                group_num = attachment_point.group_num
                r_group_atom.SetAtomMapNum(0)
                r_group_atom.SetIsotope(group_num)
                new_idx = display_mol.AddAtom(r_group_atom)
                display_mol.AddBond(atom.GetIdx(), new_idx, Chem.BondType.SINGLE)

        return display_mol

    def order_attachment_points(self, group_numbers: List[int]):
        group_map = {a.group_num: a for a in self.attachment_points}
        self.attachment_points = []
        for group_num in group_numbers:
            self.attachment_points.append(group_map.setdefault(group_num, None))

    def remove_groups(self, group_indices_to_remove):
        group_num_to_remove = []
        for group_idx in group_indices_to_remove:
            if self.attachment_points[group_idx]:
                group_num_to_remove.append(self.attachment_points[group_idx].group_num)
            del self.attachment_points[group_idx]
        for atom in self.matching_mol.GetAtoms():
            group_num = atom.GetAtomMapNum()
            if group_num > 0 and group_num in group_num_to_remove:
                atom.SetAtomMapNum(0)
                atom.ClearProp('label')
                atom.ClearProp('rgroup')

    def match_attachment_points_to_sidechain(self, sidechain: Mol) -> List[AttachmentPoint]:
        """
        Returns a list of attachment points if this side chain matches the core.  If a match is found labels the wildcard atom
        in the sidechain with the R group number of the attachment point.

        Typically there will only be one unique attachment point- though a cyclic sidechain may match 2 or more core
        attachment points and an attachment point may appear more than once

        :param sidechain: sidechain fragment as class:`rdkit.Chem.Mol` molecule
        :return: Matching attachment points or None
        """
        attachment_points = list()
        for atom in sidechain.GetAtoms():
            # after Chem.ReplaceCore(core, replaceIndex=True) the fragment wildcard atoms (attachment points)
            # will have the matching core atom index set in the isotope property
            if atom.GetAtomicNum() == 0:
                match = False
                for ap in self.attachment_points:
                    if ap and atom.GetIsotope() == ap.core_index:
                        atom.SetAtomMapNum(ap.core_index)
                        attachment_points.append(ap)
                        match = True
                if not match:
                    return list()
        return attachment_points

    def fix_sidechain(self, attachment_point: AttachmentPoint, sidechain: Mol) -> Mol:
        """
        Restores cyclicity in sidechains which have multiple mappings to the same attachment point. After
        Chem.ReplaceCore cycles at the R-group positions will be broken and replaced with multiple wildcard atoms

        :param attachment_point: attachment point as :class:`AttachmentPoint` object
        :param sidechain: sidechain fragment as :class:`rdkit.Chem.Mol` molecule
        :return: fixed structure
        """
        mol = RWMol(sidechain)
        matching_atoms = [atom for atom in sidechain.GetAtoms() if
                          atom.GetAtomicNum() == 0 and atom.GetIsotope() == attachment_point.core_index]
        assert (len(matching_atoms)) > 0
        if len(matching_atoms) == 1:
            return sidechain
        keep_atom = matching_atoms[0]
        remove_atoms = matching_atoms[1:]
        for remove_atom in remove_atoms:
            neighbour = remove_atom.GetNeighbors()[0]
            bond_type = mol.GetBondBetweenAtoms(neighbour.GetIdx(), remove_atom.GetIdx()).GetBondType()
            mol.RemoveAtom(remove_atom.GetIdx())
            mol.AddBond(keep_atom.GetIdx(), neighbour.GetIdx(), bond_type)
        return mol

    def match_sidechains(self, mol: Mol, match: List[int], sidechains: List[Mol]) -> List[Optional[Attachment]]:
        """
        Match all side chains in a decomposition to the r group positions in this core.

        The list returned is the same length as self.attachment_points (the list of available R-group positions).
        If there is an side chain attached at that position the output list will contain an :class:`Attachment`
        object, otherwise it will be None

        :param sidechains: a list of side chains as :class:`rdkit.Chem.Mol` molecules
        :return: A list of optional attachments
        """

        def sidechain_mol_for_ap(attachment_point, sidechain_matches):
            if not attachment_point:
                return None
            sidechain_mols = [self.fix_sidechain(attachment_point, sidechains[i]) for i, sidechain_list in
                              enumerate(sidechain_matches) if attachment_point in sidechain_list]
            if not sidechain_mols:
                matching_atom = mol.GetAtomWithIdx(match[attachment_point.core_index])
                if matching_atom.GetNumImplicitHs() > 0:
                    h_mol = Chem.MolFromSmiles("[*][H]")
                    h_mol.GetAtomWithIdx(0).SetIsotope(attachment_point.core_index)
                    return Attachment(attachment_point, h_mol)
                return None
            elif len(sidechain_mols) == 1:
                return Attachment(attachment_point, sidechain_mols[0])
            else:
                combineMol = functools.reduce(Chem.CombineMols, sidechain_mols)
                return Attachment(attachment_point, combineMol)

        sidechain_matches = [self.match_attachment_points_to_sidechain(sc) for sc in sidechains]
        # check all sidechains match
        if any(not apl for apl in sidechain_matches):
            return None

        # convert sidechain lists to attachment list
        attachments = [sidechain_mol_for_ap(ap, sidechain_matches) for ap in self.attachment_points]
        # check all (defined/required) attachment points are used
        # if bond orders are not defined a attachment point may not exist (e.g. the attachment exists in hexane,
        # but not benzene for a ring template)
        if any(ap and ap.required and not a for a, ap in zip(attachments, self.attachment_points)):
            return None

        return attachments


class DecompositionSummary(NamedTuple):
    """
    A class that is use to store molecular decompositions.

    The list of side chains is the same length as core.attachment_points (the list of available R-group positions).
    If there is an side chain attached at that position the output list will contain an :class:`Attachment`
    object, otherwise it will be None

    Attributes:

        - molecule: the molecule that is decomposed (as :class:`rdkit.Chem.Mol`)

        - core: the decomposition core (class :class:`Core`)

        - sidechains: A list of optional attachments

        - mapping: mapping of core to molecule

    """
    molecule: Mol
    core: Core
    sidechains: List[Optional[Attachment]]
    mapping: Tuple[int, ...]

    def __str__(self):
        mol_str = Chem.MolToSmiles(self.molecule, True)
        core_str = Chem.MolToSmarts(self.core.mol, True)
        sidechain_strs = [Chem.MolToSmiles(s.mol, True) if s else "NC" for s in self.sidechains]
        return 'MOL {}\nCORE #{} {}\nGROUPS\n{}\n'.format(mol_str, self.core.index, core_str, '\n'.join(sidechain_strs))

    def __eq__(self, other: 'DecompositionSummary') -> bool:
        """
        :param other: another decomposition
        :return: True if the two decompositions are the same
        """
        if other is None:
            return False
        if self.core.index != other.core.index:
            return False
        for sidechain, other_sidechain in zip(self.sidechains, other.sidechains):
            if sidechain != other_sidechain:
                return False
        return True

    def __ne__(self, other: 'DecompositionSummary') -> bool:
        return not self.__eq__(other)

    def contains(self, other: 'DecompositionSummary') -> bool:
        """
        :param other: another decomposition
        :return: True if this decomposition is a superset of the other decomposition
        """
        if self.core.index != other.core.index:
            return False
        if len([s for s in self.sidechains if self.heavy_sidechain(s)]) < len(
                [s for s in self.sidechains if self.heavy_sidechain(s)]):
            return False
        for sidechain, other_sidechain in zip(self.sidechains, other.sidechains):
            if not other_sidechain:
                continue
            if not sidechain:
                return False
            if not other_sidechain.heavy_sidechain():
                continue
            if not sidechain.heavy_sidechain():
                return False
            elif sidechain != other_sidechain:
                return False
        return True

    def symmetric(self, other: 'DecompositionSummary') -> bool:
        if self.core.index != other.core.index:
            return False
        if len(self.sidechains) != len(other.sidechains):
            return False
        if set(self.mapping) != set(other.mapping):
            return False
        for ap, sc in zip(self.core.attachment_points, self.sidechains):
            if not ap or not sc:
                continue
            if not sc.heavy_sidechain():
                continue
            mol_index = self.mapping[ap.core_index]
            assert (sc.attachment_point == ap)
            other_core_index = other.mapping.index(mol_index)
            other_sc = next(s for s in other.sidechains if s and s.attachment_point.core_index == other_core_index)
            if not sc.same_sidechain(other_sc, strip_atom_map_nums=True):
                return False
        return True

    @classmethod
    def heavy_sidechain(cls, sidechain: AttachmentPoint):
        return sidechain and sidechain.heavy_sidechain()

    def to_molecule_list(self, generate_coordinates: bool = False) -> List[Optional[Mol]]:
        """
        Converts this decomposition to a list of rdkit molecules. The list contains molecule, core and R groups.
        If an R group is not present in the molecule it will be None in the list

        :param generate_coordinates: set True to create 2D coordinates for depiction
        :return: the list of optional molecules
        """
        mol_list = []
        mol_list.append(self.molecule)
        mol_list.append(self.core.display_mol())
        mol_list.extend([s.mol if s else None for s in self.sidechains])
        if generate_coordinates:
            for mol in [s.mol for s in self.sidechains if s]:
                AllChem.Compute2DCoords(mol)
            AllChem.Compute2DCoords(self.molecule)
        if hasattr(self.molecule, '__sssAtoms'):
            delattr(self.molecule, '__sssAtoms')
        return mol_list

    def relabel_attachments(self) -> None:
        """
        Relabel attachment atoms with R group number- should only be done once all matching is completed
        """
        # sidechains molecules may be shared between attachments (two sidechains may have the same molecule)
        sidechain_mols = {s.mol for s in self.sidechains if s}
        for sidechain in sidechain_mols:
            if sidechain:
                for atom in sidechain.GetAtoms():
                    for attachment_point in self.core.attachment_points:
                        if attachment_point and atom.GetAtomicNum() == 0 and atom.GetIsotope() == attachment_point.core_index:
                            # the rdkit model uses isotope for smarts labels, but this does, not work for smiles
                            # use smiles reaction class labels instead
                            atom.SetIntProp('molAtomMapNumber', attachment_point.group_num)
                    if atom.GetAtomicNum() == 0 and atom.GetIsotope() > 0:
                        atom.SetIsotope(0)


class SymmetryDecompositionSummary:

    def __init__(self, molecule: Mol, core: Core, sidechains: List[Optional[Attachment]],
                 mapping: Tuple[int, ...]):
        self.decomposition_summary = DecompositionSummary(molecule, core, sidechains, mapping)
        self.symmetric_decompositions = []  # type: List[DecompositionSummary]

    def add_symmetric_decomposition(self, decomp: 'SymmetryDecompositionSummary'):
        if all(decomp.decomposition_summary != d for d in self.symmetric_decompositions):
            self.symmetric_decompositions.append(decomp.decomposition_summary)

    def __getattr__(self, item: object):
        if hasattr(self.decomposition_summary, item):
            return getattr(self.decomposition_summary, item)
        raise AttributeError

    def relabel_attachments(self):
        self.decomposition_summary.relabel_attachments()
        for sd in self.symmetric_decompositions:
            sd.relabel_attachments()


class RgroupDecomposer(Frozen):
    """
    A class to perform multi-core decomposition of a set of molecules

    Attributes:

        - cores: a list of cores to decompose molecules

        - molecules: the input molecules

        - decomposition: the results of decomposition.  This a list with one entry for each input molecule.
          If there are any decompositions for a molecule its entry will be a list of those decompositions
          (each of class :class:`DecompositionSummary`), if there are no decompositions for a molecule
          the entry will be None

        - match_first_core_only: set true if only the results for the first matching core should be included.
          By default, matches for all cores are included. If this is set and the first matching core has multiple
          matches the one will the largest number of R groups will be retained
    """

    def __init__(self):
        """
        Empty constructor
        """
        self.cores = None  # type: List[Core]
        self.molecules = None  # type: List[Mol]
        self.decomposition = None  # type: List[Optional[List[DecompositionSummary]]]
        self.match_first_core_only = False
        # Normally either the R group is specified in the query or we allow rgroups at any free valance
        # set allow_rgroups_at_any_free_valance True to allow both
        self.allow_rgroups_at_any_free_valance = False
        self.r_group_numbers = []

    def decompose(self, core_smarts: List[str], structure_smiles: List[str]) -> None:
        """
        Performs decomposition

        :param core_smarts: A list of core smarts patterns
        :param structure_smiles: A list of molecule smiles
        """
        molecules = [Chem.MolFromSmiles(m) for m in structure_smiles]
        core_molecules = [Chem.MolFromSmarts(c) for c in core_smarts]

        self.decompose_molecules(core_molecules, molecules)

    def decompose_molecules(self, core_molecules: List[Mol], molecules: List[Mol],
                            process_symmetric_groups=False) -> None:
        """
        Performs decomposition

        :param core_molecules: a list of query core molecules (each of :class:`rdkit.Chem.Mol`)
        :param molecules: a list of target molecules to decompose (each of :class:`rdkit.Chem.Mol`)
        :return:
        """
        cores = [Core(index, c, self.allow_rgroups_at_any_free_valance) for index, c in enumerate(core_molecules)]
        self.cores = cores
        self.molecules = molecules
        group_numbers = set()
        for core in cores:
            for ap in core.attachment_points:
                group_numbers.add(ap.group_num)
        group_numbers = list(group_numbers)
        group_numbers.sort()
        self.r_group_numbers = group_numbers
        for core in cores:
            core.order_attachment_points(group_numbers)

        decomposition = []
        match_no = 0
        for mol_no, mol in enumerate(molecules):
            if not mol:
                decomposition.append(None)
                continue
            if any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
                print('Mol No {} smiles {} contains dummy atoms, skipping'.format(mol_no, Chem.MolToSmiles(mol)))
                decomposition.append(None)
                continue
            mol_decomposition = []
            for core in cores:
                # Need to get all symmetric matches as, while the core may symmetric, the R group labels on the core may not be
                for match in mol.GetSubstructMatches(core.matching_mol, uniquify=False):
                    rest = Chem.ReplaceCore(mol, core.matching_mol, match, labelByIndex=True, replaceDummies=True)
                    sidechains = Chem.GetMolFrags(rest, asMols=True, sanitizeFrags=False)
                    # Isotope in sidechain dummy atoms has core index
                    sidechain_matches = core.match_sidechains(mol, match, sidechains)
                    # AtomMapNum is sidechain dummy atom has our rgroup index
                    if sidechain_matches:
                        match_no += 1
                        args = [mol, core, sidechain_matches, match]
                        core_decomp = SymmetryDecompositionSummary(
                            *args) if process_symmetric_groups else DecompositionSummary(*args)
                        present = False
                        for index, other_core_decomp in enumerate(mol_decomposition):
                            if core_decomp.core.index == other_core_decomp.core.index:
                                # check to see if this decomp is present
                                if core_decomp == other_core_decomp:
                                    present = True
                                    break
                                # check to see if a superset of this decomp is already present
                                elif other_core_decomp.contains(core_decomp):
                                    present = True
                                    break
                                # check to see if this decomp is a superset of another decomp- replace if it is
                                elif core_decomp.contains(other_core_decomp):
                                    mol_decomposition[index] = core_decomp
                                    present = True
                                    break
                                # check to see if this decomp is a symmetric match of another decomp
                                elif core_decomp.symmetric(other_core_decomp):
                                    present = True
                                    if process_symmetric_groups:
                                        other_core_decomp.add_symmetric_decomposition(core_decomp)
                                    break
                                    pass
                        if not present:
                            mol_decomposition.append(core_decomp)
                if self.match_first_core_only and len(mol_decomposition) > 0:
                    # if first matching core has multiple mappings, pick one with largest number of heavy-atom r -groups
                    mol_decomposition.sort(key=lambda x: len([s for s in x.sidechains if s and s.heavy_sidechain()]),
                                           reverse=True)
                    mol_decomposition = [mol_decomposition[0]]
                    break
            # relabel attachment points to r groups only once all matching is done
            for decomp in mol_decomposition:
                decomp.relabel_attachments()
            decomposition.append(mol_decomposition if len(mol_decomposition) else None)

        self.decomposition = decomposition
        self.remove_empty_groups()

    def remove_empty_groups(self):
        def first_free_group_num() -> int:
            group_num = 1
            while group_num in self.r_group_numbers:
                group_num += 1
            return group_num

        def rename_groups(mol: Mol) -> None:
            for atom in mol.GetAtoms():
                group_num = atom.GetAtomMapNum()
                if group_num > 0:
                    new_group_num = r_group_mapping[group_num]
                    atom.SetAtomMapNum(0)
                    atom.SetIsotope(new_group_num)

        user_groups = set()
        for core in self.cores:
            for group_num in core.user_group_numbers:
                user_groups.add(group_num)
        group_idxs_to_remove = []
        for group_idx, group_num in enumerate(self.r_group_numbers):
            if group_num in user_groups:
                continue
            if self.is_empty_group(group_idx):
                group_idxs_to_remove.append(group_idx)

        group_idxs_to_remove.sort(reverse=True)
        for core in self.cores:
            core.remove_groups(group_idxs_to_remove)
        for group_idx in group_idxs_to_remove:
            del self.r_group_numbers[group_idx]

        for decomp in self.decomposition:
            if decomp:
                for s in decomp:
                    for group_idx in group_idxs_to_remove:
                        if len(s.sidechains) > group_idx:
                            del s.sidechains[group_idx]

        r_group_mapping = {}
        for index, group_num in enumerate(self.r_group_numbers):
            if group_num in user_groups:
                r_group_mapping[group_num] = group_num
                continue
            free_group_num = first_free_group_num()
            new_group_num = free_group_num if free_group_num < group_num else group_num
            r_group_mapping[group_num] = new_group_num
            self.r_group_numbers[index] = new_group_num

        for core in self.cores:
            core.attachment_points = [
                AttachmentPoint(a.core_index, r_group_mapping[a.group_num], a.required) if a else None for a in
                core.attachment_points]
        for decomp in self.decomposition:
            if decomp:
                for s in decomp:
                    mols = {s.mol for s in s.sidechains if s}
                    for mol in mols:
                        rename_groups(mol)

    def is_empty_group(self, group_num):
        core_decomps = [c for d in self.decomposition if d for c in d if c]
        column = [c.sidechains[group_num].mol for c in core_decomps if
                  c and len(c.sidechains) > group_num and c.sidechains[group_num]]
        if len(column) == 0:
            return True
        mol_weights = [HeavyAtomMolWt(m) for m in column if m]
        return all(w == 0 for w in mol_weights)

    def to_molecule_grid(self, generate_rgroup_coords=False, column_major=False, include_missing=False)\
            -> List[List[Optional[Mol]]]:
        """
        Converts the decompositions to a grid of molecules

        Each row of the grid contains molecule, core and R groups.
        If an R group is not present in the molecule it will be None in the list

        :param generate_rgroup_coords: set True to generate 2D molecular co-ordinates for depiction
        :param column_major: set True to return the transpose of the grid
        :return: the grid of molecules
        """
        if include_missing:
            mol_lists = []
            for mol, decomp in zip(self.molecules, self.decomposition):
                if decomp:
                    for d in decomp:
                        mol_lists.append(d.to_molecule_list(generate_rgroup_coords))
                else:
                    mol_lists.append([mol])
        else:
            decompositions = [d for d in self.decomposition if d]
            mol_lists = [d.to_molecule_list(generate_rgroup_coords) for ds in decompositions for d in ds]
        if not mol_lists:
            return []
        max_len = max([len(m) for m in mol_lists])
        if not mol_lists:
            return mol_lists
        for m in mol_lists:
            m.extend([None] * (max_len - len(m)))
        return mol_lists if not column_major else [list(c) for c in zip(*mol_lists)]

    def number_r_groups(self):
        """
        Returns the number of available R-groups in the decomposition

        :return: number of R-groups
        """
        return len(self.r_group_numbers)


def mol_to_cores(mol_block: str, print_information: bool = False) -> List[Mol]:
    """
    Converts an MDL mol block containing multiple query cores into a list of Rdkit query molecules

    :param mol_block:
    :param print_information: set True to print out summary information
    :return: a list of core molecules of :class:`rdkit.Chem.Mol`
    """
    mol = Chem.MolFromMolBlock(mol_block)
    if not mol:
        mol = Chem.MolFromMolBlock(mol_block, sanitize=False)
    try:
        cores = list(Chem.GetMolFrags(mol, asMols=True))
    except:
        cores = list(Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False))
    # this is an attempted fix for CDJS molfiles that specify R groups in the ATOMS block
    for core in cores:
        for atom in core.GetAtoms():
            if atom.GetAtomicNum() == 0 and atom.GetIsotope() != 0:
                if "R{}".format(atom.GetIsotope()) != atom.GetSymbol():
                    atom.SetProp('_MolFileRLabel', str(atom.GetIsotope()))
                    atom.SetProp('dummyLabel', "R{}".format(atom.GetIsotope()))

    smarts = [Chem.MolToSmarts(c, True) for c in cores]
    if print_information:
        print("query smarts is {}".format(' '.join(smarts)))
        for c in cores:
            print('Core')
            _atom_labels(c)
    return cores


def _atom_labels(mol: Mol) -> None:
    for atom in mol.GetAtoms():
        t = atom.GetIsotope()
        s = atom.GetSymbol()
        ind = atom.GetIdx()
        no = atom.GetAtomicNum()

        print("Atom {} istope {} sym {} atom no {}".format(ind, t, s, no))


class Rgroup(Frozen):
    """
    A class to perform r group decomposition using new Rdkit methods.
    """

    def __init__(self):
        self.cores = None  # type: List[Mol]
        self.molecules = None  # type: List[Mol]
        self.failed = set()  # type: Set[int]
        self.columns = {}  # type: Dict[str, List[Mol]]
        self.rg = None  # type: RGroupDecomposition
        self.always_match_at_any_position = False

    def decompose(self, core_smarts: List[str], structure_smiles: List[str]):

        molecules = [Chem.MolFromSmiles(m) for m in structure_smiles]
        cores = [Chem.MolFromSmarts(c) for c in core_smarts]
        self.decompose_molecules(cores, molecules)

    def decompose_molecules(self, cores: List[Mol], molecules: List[Mol], extra_options: Dict[str, Any] = {}) -> None:
        self.cores = cores
        self.molecules = molecules

        def has_labelled_rgroups(mol):
            for atom in mol.GetAtoms():
                if is_r_group_atom(atom):
                    return True
            return False

        rdBase.DisableLog("rdApp.debug")
        has_labels = all(has_labelled_rgroups(core) for core in cores)

        options = RGroupDecompositionParameters()
        # options.matchingStrategy = RGroupMatching.GreedyChunks
        options.scoreMethod = RGroupScore.FingerprintVariance
        options.onlyMatchAtRGroups = False if self.always_match_at_any_position else has_labels
        options.alignment = RGroupCoreAlignment.MCS
        options.alignment = RGroupCoreAlignment.NoAlignment
        options.rgroupLabelling = RGroupLabelling.AtomMap
        for opt in extra_options:
            setattr(options, opt, extra_options[opt])
        rg = RGroupDecomposition(cores, options)

        for i, m in enumerate(molecules):
            if not m:
                self.failed.add(i)
            elif rg.Add(m) == -1:
                self.failed.add(i)

        rg.Process()
        self.rg = rg
        self.columns = rg.GetRGroupsAsColumns()

    def to_molecule_grid(self, generate_rgroup_coords=False, column_major=False, include_missing=True) -> \
            List[List[Optional[Mol]]]:
        """
        Converts the decompositions to a grid of molecules

        Each row of the grid contains molecule, core and R groups.
        If an R group is not present in the molecule it will be None in the list

        :param generate_rgroup_coords: set True to generate 2D molecular co-ordinates for depiction
        :param column_major: set True to return the transpose of the grid
        :return: the grid of molecules
        """

        # TODO generate coordinates option

        columns = self.columns
        missing_indices = self.failed
        molecules = self.molecules
        rows = []
        # The labels should be properly ordered
        labels = self.rgroup_labels()
        columns_index = 0

        for idx, mol in enumerate(molecules):
            row = [mol]
            if idx in missing_indices:
                if include_missing:
                    row.extend([None] * len(columns))
            else:
                for label in labels:
                    row.append(columns[label][columns_index])
                columns_index += 1
            rows.append(row)

        return rows if not column_major else [list(c) for c in zip(*rows)]

    def rgroup_labels(self):
        def label_index(label):
            if label == 'Core':
                return 0
            return int(label[1:])
        labels = list(self.columns.keys())
        labels.sort(key=label_index)
        return labels
