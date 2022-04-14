import itertools
from itertools import islice
from typing import Tuple, List, NamedTuple, Dict, FrozenSet, Set

from rdkit import Chem
from rdkit.Chem.rdchem import Mol

from ruse.util.frozen import Frozen


class CorrespondenceNode(NamedTuple):
    a_index: int
    b_index: int


class FoundCliques(Frozen):

    def __init__(self):
        self.max_size = 0  # type: int
        self.cliques = set()  # type: Set[FrozenSet[int]]

    def add_clique(self, clique: FrozenSet[int]) -> None:
        clique_size = len(clique)
        if clique_size < self.max_size:
            return
        elif clique_size == self.max_size:
            self.cliques.add(clique)
        elif clique_size > self.max_size:
            self.cliques = {clique}
            self.max_size = clique_size
        else:
            assert False


class CorrespondenceGraph(Frozen):

    def __init__(self, nodes: List[CorrespondenceNode], edges: List[Tuple[int, int]]):
        self.nodes = nodes
        self.edges = edges
        self.neighbours = [self._node_neighbours(i) for i in range(len(nodes))]

    def find_cliques(self) -> List[FrozenSet[int]]:
        vertices = set([v for v in range(len(self.nodes))])
        cliques = FoundCliques()
        self._bronkerbosch(set(), vertices, set(), cliques)
        return list(cliques.cliques)

    def _node_neighbours(self, node_index: int) -> FrozenSet[int]:
        neighbours1 = [e[1] for e in self.edges if e[0] == node_index]
        neighbours2 = [e[0] for e in self.edges if e[1] == node_index]
        return frozenset(neighbours1 + neighbours2)

    def _bronkerbosch(self, r: Set[int], p: Set[int], x: Set[int], cliques: FoundCliques) -> None:
        """
        See https://en.wikipedia.org/wiki/Bron%E2%80%93Kerbosch_algorithm

        :param r: All vertices in current clique
        :param p: Some vertices that could be added to the clique
        :param x: Vertices that do not belong to the clique
        :param cliques: List of found cliques
        :return:
        """
        if len(p) + len(r) < cliques.max_size:
            return
        if len(p) == 0 and len(x) == 0:
            cliques.add_clique(frozenset(r))
            return
        # don't think pivoting makes much of a difference
        # for v in list(p):
        pivot = max(p, key=lambda n: len(self.neighbours[n]))
        for v in list(p.difference(self.neighbours[pivot])):
            neighbours = self.neighbours[v]
            self._bronkerbosch(r.union([v]), p.intersection(neighbours), x.intersection(neighbours), cliques)
            p.remove(v)
            x.add(v)


class MCS(Frozen):

    def __init__(self, mol_a: Mol, mol_b: Mol, mapping: List[Tuple[int, int]] = [], check_degree: bool = True):
        self.mol_a = mol_a
        self.mol_b = mol_b
        self.mapping = mapping
        self.check_degree = check_degree
        self.fixed_a = frozenset([t[0] for t in mapping])
        self.fixed_b = frozenset([t[1] for t in mapping])
        self.a_to_b = {t[0]: t[1] for t in mapping}

    def build_correspondence_graph(self) -> CorrespondenceGraph:
        def mapped_neighbours(atom, mapped):
            neighbours = [n.GetIdx() for n in atom.GetNeighbors()]
            return [n for n in neighbours if n in mapped]

        unmapped_a = [a for a in range(self.mol_a.GetNumAtoms()) if a not in self.fixed_a]
        unmapped_b = [a for a in range(self.mol_b.GetNumAtoms()) if a not in self.fixed_b]
        fixed_neighbours_a = [m for a in unmapped_a for m in
                              mapped_neighbours(self.mol_a.GetAtomWithIdx(a), self.fixed_a)]
        fixed_neighbours_b = [m for a in unmapped_b for m in
                              mapped_neighbours(self.mol_b.GetAtomWithIdx(a), self.fixed_b)]
        for n in set(fixed_neighbours_a):
            unmapped_a.append(n)
        for n in set(fixed_neighbours_b):
            unmapped_b.append(n)

        group_a = self.group_by_atom_type(self.mol_a, unmapped_a)
        group_b = self.group_by_atom_type(self.mol_b, unmapped_b)
        atom_types = set(group_a.keys()).intersection(set(group_b.keys()))

        nodes = []
        for atom_type in atom_types:
            for a_index in group_a[atom_type]:
                for b_index in group_b[atom_type]:
                    atom_a = self.mol_a.GetAtomWithIdx(a_index)
                    atom_b = self.mol_b.GetAtomWithIdx(b_index)
                    assert (atom_a.GetSymbol() == atom_b.GetSymbol())
                    if atom_a.GetAtomicNum() == 0:
                        # match wildcards with same mapping number only
                        if atom_a.GetAtomMapNum() == atom_b.GetAtomMapNum():
                            nodes.append(CorrespondenceNode(a_index, b_index))
                        continue
                    if a_index in self.a_to_b and b_index != self.a_to_b[a_index]:
                        continue
                    if atom_a.IsInRing() and not atom_b.IsInRing():
                        continue
                    if not atom_a.IsInRing() and atom_b.IsInRing():
                        continue
                    if atom_a.IsInRing() and atom_b.IsInRing():
                        ring_sizes = [s for s in range(3, 20) if atom_a.IsInRingSize(s)]
                        if not any(atom_b.IsInRingSize(s) for s in ring_sizes):
                            continue
                    # Not matching degree result in a large MCS, but sometimes the results are hard to interpret
                    if self.check_degree:
                        if atom_a.GetDegree() == atom_b.GetDegree():
                            nodes.append(CorrespondenceNode(a_index, b_index))
                    else:
                        nodes.append(CorrespondenceNode(a_index, b_index))

        edges = []
        for index1, node1 in enumerate(nodes):
            for index2, node2 in islice(enumerate(nodes), index1 + 1, None):
                if node1.a_index == node2.a_index:
                    continue
                if node1.b_index == node2.b_index:
                    continue

                bond_a = self.mol_a.GetBondBetweenAtoms(node1.a_index, node2.a_index)
                bond_b = self.mol_b.GetBondBetweenAtoms(node1.b_index, node2.b_index)
                if bond_a is None and bond_b is None:
                    edges.append((index1, index2))
                elif bond_a and bond_b:
                    if bond_a.GetBondType() == bond_b.GetBondType():
                        edges.append((index1, index2))

        return CorrespondenceGraph(nodes, edges)

    @staticmethod
    def check_clique(graph: CorrespondenceGraph, clique: List[int]) -> None:
        for n1, n2 in itertools.permutations(clique, 2):
            present = (n1, n2) in graph.edges or (n2, n1) in graph.edges
            if not present:
                raise RuntimeError("Unable to find edge between {} and {}".format(n1, n2))

    @staticmethod
    def clique_to_mapping(graph: CorrespondenceGraph, clique: List[int]) -> List[Tuple[int, int]]:
        nodes = [graph.nodes[c] for c in clique]
        return [(n.a_index, n.b_index) for n in nodes]

    @staticmethod
    def mapping_to_clique(graph: CorrespondenceGraph, mapping: List[Tuple[int, int]]) -> List[int]:
        def node(a, b):
            return next(i for i, n in enumerate(graph.nodes) if n.a_index == a and n.b_index == b)

        clique = [node(a, b) for a, b in mapping]
        return clique

    def print_mapping(self, mapping: List[Tuple[int, int]]) -> None:
        for index_a, index_b in mapping:
            atom_a = self.mol_a.GetAtomWithIdx(index_a)
            atom_b = self.mol_b.GetAtomWithIdx(index_b)

            print("Atom {} {} -> {} {}"
                  .format(atom_a.GetIdx(), atom_a.GetSymbol(), atom_b.GetIdx(), atom_b.GetSymbol()))

    @staticmethod
    def group_by_atom_type(mol: Mol, atom_indices: List[int]) -> Dict[int, List[int]]:
        grouping = {}
        atoms = [mol.GetAtomWithIdx(i) for i in atom_indices]
        for atom in atoms:
            grouping.setdefault(atom.GetAtomicNum(), []).append(atom.GetIdx())
        return grouping

    def find_mcs(self):
        correspondence_graph = self.build_correspondence_graph()
        cliques = correspondence_graph.find_cliques()

        mappings = [self.clique_to_mapping(correspondence_graph, c) for c in cliques]
        return mappings


def _combination_ok(fragments: List[Mol], mappings: Tuple[Tuple[int, ...], ...]) -> bool:
    """
    Evaluates the atom mappings from a set of substructure matches and returns false if any atom is mapped more than
    once (i.e if the substructure matches overlap)

    :param fragments:
    :param mappings:
    :return:
    """
    mapped_atoms = set()
    for fragment, map in zip(fragments, mappings):
        for frag_idx, mol_idx in enumerate(map):
            if fragment.GetAtomWithIdx(frag_idx).GetAtomicNum() == 0:
                continue
            if mol_idx in mapped_atoms:
                return False
            mapped_atoms.add(mol_idx)
    return True


def _remove_wildcard(fragments: List[Mol], mappings: Tuple[Tuple[int, ...], ...]) -> List[List[int]]:
    """
    Removes attachment points from a set of substructure matches

    :param fragments:
    :param mappings:
    :return:
    """
    new_mappings = []
    for fragment, map in zip(fragments, mappings):
        new_map = [mol_idx for frag_idx, mol_idx in enumerate(map)
                   if fragment.GetAtomWithIdx(frag_idx).GetAtomicNum() != 0]
        new_mappings.append(new_map)
    return new_mappings


def _use_clique(map: List[Tuple[int, int]], clique: List[Tuple[int, int]]) -> bool:
    """
    Return true if we want to keep this clique.  We discard cliques with small disconnected mappings.

    :param map:
    :param clique:
    :return:
    """
    a_indices = [t[0] for t in clique]
    a_mapping = [t[0] for t in map]
    in_a_mapping = [a for a in a_indices if a in a_mapping]
    if len(a_indices) < 4 and len(in_a_mapping) == 0:
        return False
    return True


def _fix_missing_product_bonds(mol_a: Mol, mol_b: Mol, map: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """
    If we have two mapped bonded atoms in the query which are not bonded in the product, we want to remove one of
    them from the mapping, as there is a discontinuity in the mapping here.

    The query is split on the bond and the atom in the smaller fragment is chosen for removal.

    :param mol_a:
    :param mol_b:
    :param map:
    :return:
    """
    a_to_b = {m[0]: m[1] for m in map}
    a_indices, _ = zip(*map)
    missing_bond_atoms = []

    for bond in mol_a.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        if a1 in a_indices and a2 in a_indices:
            b1 = a_to_b[a1]
            b2 = a_to_b[a2]
            if not mol_b.GetBondBetweenAtoms(b1, b2):
                missing_bond_atoms.append((a1, a2, bond))

    removed_atoms = []
    for a1, a2, bond in missing_bond_atoms:
        if a1 in removed_atoms or a2 in removed_atoms:
            continue
        if bond.GetBeginAtom().GetAtomicNum() == 0:
            removed_atoms.append(a2)
            continue
        if bond.GetEndAtom().GetAtomicNum() == 0:
            removed_atoms.append(a1)
            continue
        frag_mol = Chem.FragmentOnBonds(mol_a, [bond.GetIdx()], dummyLabels=[(a1, a2)])
        frags = Chem.GetMolFrags(frag_mol)
        if len(frags) == 2:
            smallFrag = frags[0] if len(frags[0]) < len(frags[1]) else frags[1]
            smallFragAtoms = [frag_mol.GetAtomWithIdx(a) for a in smallFrag]
            link_atoms = [a for a in smallFragAtoms if a.GetAtomicNum() == 0 and a.GetIsotope() in [a1, a2]]
            assert len(link_atoms) == 1
            index_to_remove = link_atoms[0].GetIsotope()
            removed_atoms.append(index_to_remove)
        else:
            removed_atoms.append(a1)

    new_mapping = [(a, b) for (a, b) in map if a not in removed_atoms]
    return new_mapping


def refine_clique(mol_a: Mol, mol_b: Mol, map: List[Tuple[int, int]], clique: List[Tuple[int, int]]) -> List[
    Tuple[int, int]]:
    """
    Refines the clique by removing atoms which have no path back to an attachment point. This routine will only work
    for MMP combined molecules.

    :param mol_a:
    :param mol_b:
    :param map:
    :param clique:
    :return:
    """
    a_indices = [t[0] for t in clique]
    a_mapping = {t[0] for t in map}.union(a_indices)

    def has_path_attachment(a, attach):
        if a == attach:
            return True
        path = set(Chem.GetShortestPath(mol_a, a, attach))
        return all(p in a_mapping for p in path)

    attachment_indices = [atom.GetIdx() for atom in mol_a.GetAtoms() if atom.HasProp('molAtomMapNumber')]
    for a in list(a_indices):
        if not any(has_path_attachment(a, attach) for attach in attachment_indices):
            a_indices.remove(a)
    return [(a, b) for a, b in clique if a in a_indices]


def mcs_with_constant_smiles(mol_a: Mol, mol_b: Mol, constant_smiles: str, check_degree: bool = True) \
        -> List[Tuple[int, int]]:
    """
    Perform MCS between two molecules, with the constraint that the constant pattern maps to both molecules

    :param mol_a:
    :param mol_b:
    :param constant_smiles:
    :return:
    """

    # split constant smiles into fragment molecules
    constant_fragments = [Chem.MolFromSmarts(s) for s in constant_smiles.split('.')]
    # and find matches in both molecules (must enumerate all matches- set uniquify false)
    mol_a_matches = [mol_a.GetSubstructMatches(fragment, uniquify=False) for fragment in constant_fragments]
    mol_b_matches = [mol_b.GetSubstructMatches(fragment, uniquify=False) for fragment in constant_fragments]

    # build no-overlapping fragment conbinations
    mol_a_matches_product = [_remove_wildcard(constant_fragments, m) for m in itertools.product(*mol_a_matches)
                             if _combination_ok(constant_fragments, m)]
    mol_b_matches_product = [_remove_wildcard(constant_fragments, m) for m in itertools.product(*mol_b_matches)
                             if _combination_ok(constant_fragments, m)]
    # iterate through all matches
    cliques = []
    mapping = []
    for mol_a_mappings in mol_a_matches_product:
        for mol_b_mappings in mol_b_matches_product:
            a_indices = [i for f in mol_a_mappings for i in f]
            b_indices = [i for f in mol_b_mappings for i in f]
            # mapping is the A -> B correspondence  from the constant smiles
            mapping = [(a, b) for a, b in zip(a_indices, b_indices)]
            # use MCS algorithm to find the remaining matches
            mcs = MCS(mol_a, mol_b, mapping, check_degree=check_degree)
            cliques.extend([(c, mapping) for c in mcs.find_mcs()])

    cliques = [(_filter_partial_rings(mol_a, c), mapping) for c, mapping in cliques if _use_clique(mapping, c)]
    cliques = [(refine_clique(mol_a, mol_b, mapping, c), mapping) for c, mapping in cliques]
    if not cliques:
        if not mapping:
            raise RuntimeError('Unable to map constant fragments to structure!')
        return mapping
    # TODO - handle symmetry by storing alternative maximum cliques (may not be required now we expect full
    # ring matches including degree)
    clique, mapping = max(cliques, key=lambda c: len(c[0]))
    # remove atoms in query which have fewer bonds than in product
    clique = [(a, b) for a, b in clique if mol_a.GetAtomWithIdx(a).GetDegree() >= mol_b.GetAtomWithIdx(b).GetDegree()]
    complete_mapping = list(set(mapping).union(clique))
    # For mapped atoms if there is a bond in the product that is not present in the query, remove that query atom.
    # In most cases will not be required after the degree test.
    complete_mapping = _fix_missing_product_bonds(mol_a, mol_b, complete_mapping)
    return complete_mapping


def _filter_partial_rings(mol: Mol, mapping: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """
    Remove partial ring matches from the MCS
    :param mol:
    :param mapping:
    :return:
    """
    a_indices, _ = zip(*mapping)
    ring_info = mol.GetRingInfo()
    ring_atoms = set()
    for ring in ring_info.AtomRings():
        if all(r in a_indices for r in ring):
            ring_atoms.update(ring)

    a_indices = [a for a in a_indices if not mol.GetAtomWithIdx(a).IsInRing() or a in ring_atoms]
    return [(a, b) for a, b in mapping if a in a_indices]


def _refine_fragment_clique(mol_a: Mol, clique: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """
    Refines the clique by removing atoms which have no path back to an attachment point. This routine will only work
    for MMP combined molecules.

    :param mol_a:
    :param mol_b:
    :param map:
    :param clique:
    :return:
    """
    a_indices = [t[0] for t in clique]

    def has_path_attachment(a, attach):
        if a == attach:
            return True
        path = set(Chem.GetShortestPath(mol_a, a, attach))
        return all(p in a_indices for p in path)

    attachment_indices = [atom.GetIdx() for atom in mol_a.GetAtoms() if atom.GetAtomicNum() == 0]
    for a in list(a_indices):
        if not any(has_path_attachment(a, attach) for attach in attachment_indices):
            a_indices.remove(a)
    return [(a, b) for a, b in clique if a in a_indices]


def mcs_match_fragments(mol_a: Mol, mol_b: Mol, check_degree: bool = True) -> List[Tuple[int, int]]:
    """
    Assuming mol_a and mol_b represent a transform create an MCS mapping between the two

    :param mol_a:
    :param mol_b:
    :param check_degree:
    :return:
    """
    mcs = MCS(mol_a, mol_b, check_degree=check_degree)
    cliques = mcs.find_mcs()
    cliques = [_refine_fragment_clique(mol_a, c) for c in cliques if c]
    cliques = [_filter_partial_rings(mol_a, c) for c in cliques if c]
    cliques = [_fix_missing_product_bonds(mol_a, mol_b, c) for c in cliques if c]
    clique = max(cliques, key=lambda c: len(c)) if cliques else []

    return clique


def mcs_match_fragment_smiles(smi_a: str, smi_b: str, check_degree: bool = True) -> List[Tuple[int, int]]:
    """
    Assuming smi_a and smi_b represent a transform create an MCS mapping between the two

    :param smi_a:
    :param smi_b:
    :param check_degree:
    :return:
    """
    mol_a = Chem.MolFromSmarts(smi_a)
    mol_b = Chem.MolFromSmarts(smi_b)
    return mcs_match_fragments(mol_a, mol_b, check_degree)


def smiles_to_mcs_with_constant_smiles(smiles_a: str, smiles_b: str, constant_smiles: str, check_degree: bool = True) \
        -> List[Tuple[int, int]]:
    return mcs_with_constant_smiles(Chem.MolFromSmiles(smiles_a), Chem.MolFromSmiles(smiles_b), constant_smiles,
                                    check_degree)
