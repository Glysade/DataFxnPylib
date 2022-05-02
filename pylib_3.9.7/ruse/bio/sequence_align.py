"""
=================
sequence_align.py
=================

Copyright (C) 2017-2022 Glysade, LLC

A classs to handle multiple sequence alignment (MSA)

"""

import os
import uuid
from enum import Enum
from typing import Union, List, Dict, Tuple, Optional

from Bio import AlignIO, Phylo
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import ClustalOmegaCommandline, ClustalwCommandline, MuscleCommandline
from Bio.Nexus.Trees import Tree
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation, AbstractPosition, AfterPosition, \
    BeforePosition, BetweenPosition, ExactPosition, OneOfPosition, UnknownPosition, WithinPosition
from Bio.SeqRecord import SeqRecord

from ruse.bio.bio_util import is_dna, is_protein, ok_sequence
from ruse.bio.phylo_tree import PhyloTreeBuilder, PhyloTree
from ruse.util.frozen import Frozen

OptionDataTypes = Union[str, int, float]


class SequenceAlignmentMethod(Enum):
    CLUSTALO = 1
    MUSCLE = 2
    CLUSTALW = 3


class MultipleSequenceAlignment(Frozen):
    """
    A class to perform MSA using ClustalO on input sequences

    Attributes:
        * options: a key-value dictionary of options to be passed to ClustalO
        * input_sequences: a List of :class:`Bio.SeqRecord.SeqRecord` sequence objects to be aligned
        * aligned_sequences: the input sequences aligned as a List of :class:`Bio.SeqRecord.SeqRecord`
        * alignment: Biopython parsing of alignment as :class:`Bio.Align.Bio.MultipleSeqAlignment`
        * basename: UUID string for common resource naming
        * tree: guide tree for alignment as :class:`Bio.Nexus.Trees.Tree`
        * distances: distance of each sequence from root in guide tree
    """

    def __init__(self):
        """
        Empty constructor
        """
        self.options = {}  # type: Dict[str, OptionDataTypes]
        self.input_sequences = []  # type: List[SeqRecord]
        self.alignment = None  # type: MultipleSeqAlignment
        self.aligned_sequences = []  # type: List[SeqRecord]
        self.basename = str(uuid.uuid4())  # type: str
        self.tree = None  # type: Tree
        self.distances = None  # type: List[float]
        self.distance_lookup = None  # type Dict[str, float]
        self.alignment_method = None  # type: SequenceAlignmentMethod

    def align_sequences(self, options: Dict[str, OptionDataTypes], sequences: List[SeqRecord],
                        alignment_method: SequenceAlignmentMethod = SequenceAlignmentMethod.CLUSTALO) -> None:
        """
        Performs MSA using ClustalO on the input sequences and creates alignment attributes. A guide tree with also be
        created if there are 3 or more input sequences

        :param options: a key-value dictionary of options to be passed to ClustalO
        :param sequences: a List of :class:`Bio.SeqRecord.SeqRecord` sequence objects to be aligned
        """

        self.input_sequences = sequences
        self.options = options
        self.alignment_method = alignment_method

        # create input file
        in_file = "{}_in.fasta".format(self.basename)
        out_file = "{}_out.aln".format(self.basename)
        guide_file = "{}_in.dnd".format(self.basename)  # _in to be compatable with clustal!

        # make sure we only have id in input labels
        in_sequences = [SeqRecord(s.seq, id=s.id, name='', description='') for s in sequences]
        # purge gap only or empty sequences:
        in_sequences = [s for s in in_sequences if ok_sequence(s)]

        if not in_sequences:
            aligned_sequences = []
            self.alignment = None
        elif len(in_sequences) == 1:
            aligned_sequences = in_sequences
        else:
            with open(in_file, 'w') as f:
                SeqIO.write(in_sequences, f, "fasta")

            if alignment_method == SequenceAlignmentMethod.CLUSTALO:
                # create and execute clustalo command line
                # TODO - need to pass options
                if os.name == 'nt':
                    clustalo_exe = 'clustalo.exe'
                elif os.name == 'posix':
                    clustalo_exe = 'clustalo'
                else:
                    raise RuntimeError("Unable to determine clustal exe for os {}".format(os.name))
                clustalo_cline = ClustalOmegaCommandline(clustalo_exe, infile=in_file, outfile=out_file, distmat_full=False,
                                                         # outputorder='tree-order',
                                                         outfmt="clustal", verbose=True, auto=True,
                                                         guidetree_out=guide_file)
                # print("Clustalo command line {}".format(clustalo_cline))
                stdout, stderr = clustalo_cline()

            elif alignment_method == SequenceAlignmentMethod.CLUSTALW:
                # clustalw.  Note that aligned sequence names may be truncated in clustalw output
                # clustalw is currently not supported and not working!
                # clustalw_cline = ClustalwCommandline('clustalw', infile=in_file, outfile=out_file, outorder='INPUT')
                clustalw_cline = ClustalwCommandline('clustalw', infile=in_file, outfile=out_file)
                stdout, stderr = clustalw_cline()

            elif alignment_method == SequenceAlignmentMethod.MUSCLE:
                if os.name == 'nt':
                    muscle_exe = 'muscle3.8.31_i86win32.exe'
                elif os.name == 'posix':
                    muscle_exe = os.path.abspath(
                        os.path.join(os.path.dirname(__file__), "../../bin/linux_bin/muscle"))
                else:
                    raise RuntimeError("Unable to determine muscle exe for os {}".format(os.name))
                # tree2 to save newick guide tree after 2nd iteration, stable to have output order match input order
                # clw is clustal flag, stable is no longer supported, so rearrange sequences after alignment
                muscle_cline = MuscleCommandline(muscle_exe, input=in_file, out=out_file, clw=True,
                                                 tree2=guide_file)
                stdout, stderr = muscle_cline()

            else:
                raise AttributeError('Unknown alignment method {}'.format(str(alignment_method)))

            align = AlignIO.read(out_file, "clustal")
            aligned_sequences = [s for s in align]

            # muscle sequence output order is not the same as input order- reorder here assuming that ids are unique
            # if alignment_method == SequenceAlignmentMethod.MUSCLE:
            # Do this anyway as
            # not all input sequences are present in the output- the output is ordered
            # as the same as the input, but some sequences may not be included in the alignment
            # (for example if they are gap only)
            aligned_sequences_lookup = {s.id: s for s in aligned_sequences}
            aligned_sequences = [aligned_sequences_lookup.get(s.id) for s in sequences]
            self.alignment = align

        self.aligned_sequences = aligned_sequences
        if len(sequences) > 1:
            if PhyloTreeBuilder.find_fasttree_exe():
                phylo_tree_builder = PhyloTreeBuilder()
                self.tree = phylo_tree_builder.build_fasttree_tree(out_file)
                phylo_tree_builder.cleanup()
            else:
                self.tree = Phylo.read(guide_file, "newick")
                self.tree.ladderize()

    def add_distances(self) -> None:
        """
        Add an evolutionary distance for each sequence using the guide tree or fasttree
        """
        if not self.input_sequences:
            distance_lookup = {}
        elif len(self.input_sequences) == 1:
            distance_lookup = {self.input_sequences[0].id: 1.0}
        else:
            tree = PhyloTree(self.tree)
            distance_lookup = tree.leaf_distances(last_distance_only=not PhyloTreeBuilder.find_fasttree_exe())
        assert (len(distance_lookup) == len([s for s in self.aligned_sequences if s]))
        self.distance_lookup = distance_lookup
        self.distances = [distance_lookup.get(s.id) for s in self.input_sequences]

    def copy_features_to_aligned_sequences(self) -> None:
        """
        Copies sequence features from input sequences to aligned sequences
        """

        id_to_input = {s.id: s for s in self.input_sequences}
        for output_sequence in self.aligned_sequences:
            input_sequence = id_to_input[output_sequence.id]
            copy_features_to_gapped_sequence(input_sequence, output_sequence)

    def clean_files(self) -> None:
        """
        Delete MSA ClustalO files
        """

        for path in ["{}_in.fasta", "{}_out.aln", "{}_out.dnd", "{}_in.dnd"]:
            path = path.format(self.basename)
            if os.path.exists(path):
                os.remove(path)


def copy_features_to_gapped_sequence(in_seq: SeqRecord, out_seq: SeqRecord, local: bool = False) -> None:
    """
    Copies sequence features from an input sequence to another sequence, adjusting for inserted or deleted gaps

    :param in_seq: The input sequence with no gaps as a :class:`Bio.SeqRecord.SeqRecord` object
    :param out_seq: The gapped aligned output sequence as a :class:`Bio.SeqRecord.SeqRecord` object
    :param local: Set true if the alignment is local (partial continuous match)
    """

    if local:
        ungapped = str(out_seq.seq.ungap('-'))
        start = str(in_seq.seq).index(ungapped)
        end = start + len(ungapped)
        in_seq = in_seq[start:end]
    # when mapping from in_seq to out_seq account for gaps in both sequences
    in_position = 0
    new_positions = [None] * (len(in_seq) + 1)
    for out_position, out_nt in enumerate(out_seq.seq):
        if in_position == len(in_seq):
            break
        in_nt = in_seq[in_position]
        while in_nt == '-':
            new_positions[in_position] = -1
            in_position += 1
            if in_position == len(in_seq):
                break
            in_nt = in_seq[in_position]
        if out_nt != '-':
            if out_nt != in_nt:
                raise ValueError('copy_features_to_aligned_sequence: sequence does not match')
            new_positions[in_position] = out_position
            in_position += 1
            if in_position == len(in_seq):
                break
    while in_position < len(in_seq):
        in_nt = in_seq[in_position]
        if in_nt != '-':
            raise ValueError('copy_features_to_aligned_sequence: sequence does not match')
        new_positions[in_position] = -1
        in_position += 1

    new_positions[len(in_seq)] = new_positions[in_position - 1] + 1
    _build_new_features(in_seq, out_seq, new_positions)


def _build_new_features(in_seq: SeqRecord, out_seq: SeqRecord, new_positions: List[int]) -> None:
    """
    Copies sequence features from an input sequence to an aligned sequence, adjusting for inserted gaps

    :param in_seq: The input sequence with no gaps as a :class:`Bio.SeqRecord.SeqRecord` object
    :param out_seq: The gapped aligned output sequence as a :class:`Bio.SeqRecord.SeqRecord` object
    :param new_positions: A list that maps positions in the input sequence to the aligned sequence, used to map features onto the output sequence
    """

    for feature in in_seq.features:
        old_location = feature.location
        new_location = None
        if isinstance(old_location, CompoundLocation):
            new_parts = [_adjust_location(l, new_positions) for l in old_location.parts]
            if all(new_parts):
                new_location = CompoundLocation([_adjust_location(l, new_positions) for l in old_location.parts],
                                                old_location.operator)
        else:
            new_location = _adjust_location(old_location, new_positions)
        if new_location is not None:
            new_feature = SeqFeature(location=new_location, type=feature.type,
                                     location_operator=feature.location_operator,
                                     id=feature.id, qualifiers=dict(feature.qualifiers.items()))
            out_seq.features.append(new_feature)

    out_seq.dbxrefs = in_seq.dbxrefs[:]
    out_seq.annotations = in_seq.annotations.copy()


def _adjust_location(location: FeatureLocation, new_positions: List[int]) -> Optional[FeatureLocation]:
    """
    Determines a feature location in the aligned sequence

    :param location: A location in the input sequence
    :param new_positions: A list that maps positions in the input sequence to the aligned sequence
    :return: The equivalent location in the aligned sequence
    """

    new_start = _adjust_position(location.start, new_positions)
    if new_start is None:
        return None
    # biopython end is exclusive, but I'm not sure what to do if the end position is not exact
    end = location.end
    if isinstance(end, ExactPosition):
        end = ExactPosition(end - 1)
        new_end = _adjust_position(end, new_positions)
        new_end = ExactPosition(new_end + 1)
    else:
        new_end = _adjust_position(location.end, new_positions)
    if not new_end:
        return None
    return FeatureLocation(new_start, new_end, location.strand)


def _adjust_position(position: AbstractPosition, new_positions: List[int]) -> Optional[AbstractPosition]:
    """
    Maps a sequence position from the input sequence to the  aligned sequence

    :param position: A position in the input sequence
    :param new_positions: A list that maps positions in the input sequence to the aligned sequence
    :return:  The equivalent position in the aligned sequence
    """

    new_position = new_positions[position]
    if new_position == -1:
        return None
    if isinstance(position, AfterPosition):
        return AfterPosition(new_position)

    elif isinstance(position, BeforePosition):
        return BeforePosition(new_position)

    elif isinstance(position, BetweenPosition):
        # this class may be obsolete
        # this code path has not been tested
        new_left = new_positions[position._left]
        new_right = new_positions[position._right]
        if new_left == -1 or new_right == -1:
            return None;
        return BetweenPosition(new_position, new_left, new_right)

    elif isinstance(position, ExactPosition):
        return ExactPosition(new_position)

    elif isinstance(position, OneOfPosition):
        # this code path has not been tested
        old_choices = position.position_choices
        new_choices = [_adjust_position(old_choice, new_positions) for old_choice in old_choices]
        if not all(new_choices):
            return None
        return OneOfPosition(new_position, new_choices)

    elif isinstance(position, UnknownPosition):
        # afaik this should not occur for genbank files
        raise ValueError

    elif isinstance(position, WithinPosition):
        # this code path has not been tested
        new_left = new_positions[position._left]
        new_right = new_positions[position._right]
        if new_left == -1 or new_right == -1:
            return None;
        return WithinPosition(new_position, new_left, new_right)

    raise ValueError
