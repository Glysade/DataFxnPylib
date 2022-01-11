"""
=================
phylo_tree.py
=================

Copyright (C) 2017 Anodyne Informatics, LLC

Classes to handle the generation and processing of phylogenetic trees from multiple sequence alignments.

"""

import os
import uuid
from glob import glob
from typing import Union, List, Dict, Tuple, Optional

from Bio import AlignIO, Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.BaseTree import Tree
from Bio.Phylo.Applications import RaxmlCommandline, FastTreeCommandline
from Bio.Phylo.Newick import Clade

from ruse.bio.bio_util import is_dna
from ruse.util.frozen import Frozen
from ruse.util.util import is_exe

OptionDataTypes = Union[str, int, float]


class PhyloTreeBuilder(Frozen):
    """
    A class to build a phylogenetic tree from a multiple sequence alignment using an external application

    Attributes:
        * input_aln_file: Clustal format multiple sequence alignment file used to construct tree
        * basename: UUID string used to construct output files
        * alignment: a :class:`Bio.AlignIO` object of input alignment
    """

    def __init__(self):
        """
        Empty constructor
        """
        self.input_aln_file = None  # type: str
        self.basename = str(uuid.uuid4())
        self.alignment = None  # type: MultipleSeqAlignment

    def _transform_align_output(self, new_format: str = 'phylip') -> str:
        """
        Transforms the input alignment file into a new format using :class:`Bio.AlignIO`

        :param new_format: format to transform the file to (defaults to 'phylip')
        :return: file name of transformed file
        """
        assert self.alignment
        transformed_file = "{}.{}".format(self.basename, new_format)
        AlignIO.write([self.alignment], transformed_file, new_format)
        return transformed_file

    def _read_alignment_file(self) -> MultipleSeqAlignment:
        """
        Reads a ClustalO alignment file

        :return: a :class:`Bio.Align.MultipleSeqAlignment` object containing the alignment
        """
        alignment = AlignIO.read(self.input_aln_file, "clustal")
        self.alignment = alignment
        return alignment

    def build_raxml_tree(self, input_aln_file: str) -> Tree:
        """
        Creates a Biopython Tree object (of class :class:`Bio.Phylo.BaseTree.Tree`) by running raxml
        (see `https://sco.h-its.org/exelixis/web/software/raxml/index.html`).  The raxml executable must be on the PATH.

        This is not a recommended method to create trees.  It is very time consuming and deprecated in favor of ExaML.

        :param input_aln_file: path of clustal alignment file
        :return: tree created from running RaxML
        """

        self.input_aln_file = input_aln_file
        self._read_alignment_file()
        phylip_file = self._transform_align_output()

        example_seq = self.alignment[0]
        # see https://sco.h-its.org/exelixis/web/software/raxml/hands_on.html for explanation of model types
        model = 'PROTGAMMAWAG' if not is_dna(example_seq) else 'GTRGAMMA'
        cmdline = RaxmlCommandline(sequences=phylip_file, model=model, name=self.basename)
        _, _ = cmdline()
        tree_file = 'RAxML_bestTree.{}'.format(self.basename)
        tree = Phylo.read(tree_file, 'newick')
        return tree

    @classmethod
    def find_fasttree_exe(cls) -> Optional[str]:
        """
        Returns the full path to an FastTree (`http://www.microbesonline.org/fasttree/`) executable.
        Windows and Linux FastTree binaries are included in this distribution.

        :return: full path to the program, of None if no executable can be found
        """

        base = os.path.normpath("../../bin")
        if os.name == 'nt':
            return 'FastTree.exe'
            # TODO remove once everything seems to be working
            exe = os.path.abspath(
                os.path.join(os.path.dirname(__file__), base, "win_bin", "FastTree.exe"))
        elif os.name == 'posix':
            exe = os.path.abspath(
                os.path.join(os.path.dirname(__file__), base, "linux_bin", "FastTreeDbl"))
        else:
            return None

        if os.path.isfile(exe) and is_exe(exe):
            return exe

        print("Not using fasttree: file {} is not present or not executable!".format(exe))
        return None

    def build_fasttree_tree(self, input_aln_file: str) -> Tree:
        """
        Creates a Biopython Tree object (of class :class:`Bio.Phylo.BaseTree.Tree`) by running FastTree
        (see `http://www.microbesonline.org/fasttree/`).

        Assumes that FastTree binaries are available in this distribution.
        For nucleotide sequences the -gtr and -nt arguments to FastTree are set.
        No arguments are set for protein alignment.

        :param input_aln_file: path of clustal alignment file
        :return: tree created from running FastTree
        """
        self.input_aln_file = input_aln_file
        self._read_alignment_file()
        fasta_file = self._transform_align_output('fasta')

        example_seq = self.alignment[0]
        out_file = '{}.tree'.format(self.basename)

        args = {'input': fasta_file, 'out': out_file}
        if is_dna(str(example_seq.seq)):
            args['gtr'] = True
            args['nt'] = True
        exe = self.find_fasttree_exe()
        assert exe
        cmdline = FastTreeCommandline(exe, **args)
        #print(cmdline)
        _, _ = cmdline()
        tree = Phylo.read(out_file, 'newick')
        tree.root_at_midpoint()
        tree.ladderize()
        return tree

    def cleanup(self):
        """
        Removes all output and intermediate files
        :return:
        """
        for file in glob('*{}*'.format(self.basename)):
            os.remove(file)


class PhyloTree(Frozen):
    """
    A class to facilitate operations on pylogenetic trees.  Handles conversion of Biopython trees to JSON and
    calculation of leaf distances

    Attributes:
        * tree: Bioython pyhlogenetic tree (class :class:`Bio.Phylo.BaseTree.Tree`)
        * data_tree: Python dictionary representation of tree, that is JSON compatible
    """

    def __init__(self, tree: Tree):
        """
        Constructor.  Converts Biopython tree object of :class:`Bio.Phylo.BaseTree.Tree` to a JSON compatible
        dictionary.

        :param tree: Bioython pyhlogenetic tree
        """
        self.tree = tree  # type: Tree
        self.data_tree = self._tree_to_data()  # Type: Dict

    # could use Biopython's search and traversal methods in preference to these

    def _tree_to_data(self) -> Dict:
        """
        Converts the guide tree to a nested dictionary data structure. Uses a non-recursive method

        :return: tree as python data structure
        """

        root = {}
        tree = self.tree
        if tree is not None:
            nodes_to_process = [(tree.clade, root)]
            while len(nodes_to_process):
                (clade, node_data) = nodes_to_process.pop()
                self._add_node_to_parent(clade, node_data, nodes_to_process)
        return root

    @classmethod
    def _add_node_to_parent(cls, clade: Clade, node: Dict, nodes_to_process: List[Tuple[Clade, Dict]]) -> None:
        """
        Helper function to create guide tree data.  Set node properties from clade and adds all children to list of
        nodes to process.

        :param clade: Biopython class for node as :class:`Bio.Phylo.Newick.Clade`
        :param node: Empty dictionary to store node properties in
        :param nodes_to_process: List of current nodes to to process
        :return:
        """
        if clade.branch_length is not None:
            node['distance'] = clade.branch_length
        if clade.name is not None:
            node['name'] = clade.name
        if len(clade.clades):
            node['children'] = [dict() for _ in range(len(clade.clades))]
            for child_clade, child in zip(clade.clades, node['children']):
                nodes_to_process.append((child_clade, child))

    def tree_to_data_recursive(self) -> Dict:
        """
        Converts the guide tree to a nested dictionary data structure. Uses a recursive method used to
        validate :func:`_tree_to_data`

        :return: tree as python data structure
        """

        return self._clade_tree_to_data(self.tree.clade) if self.tree is not None else {}

    @classmethod
    def _clade_tree_to_data(cls, clade: Clade) -> Dict:
        node = {}
        if clade.branch_length is not None:
            node['distance'] = clade.branch_length
        if clade.name is not None:
            node['name'] = clade.name
        if len(clade.clades):
            node['children'] = [cls._clade_tree_to_data(c) for c in clade.clades]
        return node

    def leaf_distances(self, last_distance_only: bool = False) -> Dict[str, float]:
        """
        Non-recursive function to find the distance to the leaves

        :param last_distance_only: set true to return the distance from leaf to parent rather than cumulative distance to the root.  This is appropriate for UPGMA trees (such as those created by ClustalO) which have constant root to leaf distances.
        :return: A dictionary mapping leaf names to distances
        """

        data_tree = self.data_tree
        leaves = dict()
        nodes_to_process = list()
        nodes_to_process.append((0, data_tree))

        while len(nodes_to_process):
            (distance, node) = nodes_to_process.pop()
            if 'children' in node:
                for child in node['children']:
                    nodes_to_process.append((distance + child['distance'], child))
            elif 'name' in node:
                leaves[node['name']] = node['distance'] if last_distance_only else distance
            else:
                raise ValueError

        return leaves
