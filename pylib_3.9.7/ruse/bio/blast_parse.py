"""
==============
blast_parse.py
==============

Copyright (C) 2017-2022 Glysade, LLC

Classes for parsing the results from Blast searches
"""

from __future__ import annotations
import re
from io import StringIO
from typing import List, Dict, Union, Optional, TYPE_CHECKING

from Bio import Entrez
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import nt_search

from ruse.bio.bio_util import entrez_email, is_defined_sequence
from ruse.bio.blast_search import BlastSearchType
from ruse.bio.blast_utils import BlastRecords, retrieve_entrez_records, SequenceMatcher
if TYPE_CHECKING:
    from ruse.bio.igblast_parse import IgHit
from ruse.bio.sequence_align import copy_features_to_gapped_sequence
from ruse.util.frozen import Frozen

Entrez.email = entrez_email()

# allowed blast information data types
DataTypes = Union[str, int, float]


class BlastResult(Frozen):
    """
    Class to hold an NCBI Blast result, comprising a single alignment between query and target sub sequences.

    Attributes:
        * All constructor parameters
        * target_record (:class:`SeqRecord`) the target sequence downloaded from Entrez
    """

    # genbank and refseq accession and version numbers, see https://www.ncbi.nlm.nih.gov/Sequin/acc.html
    accession_version_pattern = re.compile(r"^[A-Z]+_?\d+.\d+$")
    # general id pattern
    accession_pattern = re.compile(r"^[A-Z_.\d]+$")

    # database tags see https://blast.advbiocomp.com/doc/FAQ-Indexing.html
    ncbi_tags = 'bbm bbs dbj emb gb gim gnl gp lcl pat pdb pir prf ref sp tpd tpe tpg'.split(' ')

    def __init__(self, target_id: str, target_def: str, target_length: int, align_length: int, query: str, target: str,
                 evalue: float, score: int, bits: float, accession: str, query_start: int, query_end: int,
                 search_type: BlastSearchType):
        """
        Constructs the class from NCBI blast search output fields

        :param target_id: the id or the target/database sequence
        :param target_def: the definition for the target sequence
        :param target_length: the full length of the target seqence
        :param align_length: the length of the matching alignment
        :param query: The matching query sequence
        :param target: The matching target sequence
        :param evalue: The Blast expect value for the alignment
        :param score: The Blast score
        :param bits: The Blast bit count
        :param accession: Target accession
        """

        self.target_id = target_id
        self.target_def = target_def
        self.target_length = target_length
        self.align_length = align_length
        self.query = query
        self.target = target
        self.evalue = evalue
        self.score = score
        self.bits = bits
        self.accession = accession
        self.query_start = query_start
        self.query_end = query_end
        self.search_type = search_type
        self.target_record = None  # type: Optional[SeqRecord]

    @classmethod
    def fasta_header_identifier(cls, header: str) -> Optional[str]:
        """
        Attempt to extract an identifier that can be used as an Entrez query from a header string.  The header string
        is a list of identifiers separated by  the '|' character- sit the first line of FASTA files and can be extracted
        from definitions or ids in Blast results

        :param header: string to extract identifier from
        """

        # BL_ORD_ID is an identifier used by makeblastdb for internal record tracking
        if 'BL_ORD_ID' not in header:
            if ' ' in header:
                header = header.split(' ', 1)[0]

            if '|' in header:
                terms = re.split("\\|+", header)
                if terms[0] == 'gi':
                    terms = terms[2:]
                if len(terms) == 3 and terms[0] == 'pdb' and re.match(r'[A-Z]$', terms[2]):
                    # account for PDB strand matches
                    return "{}_{}".format(terms[1], terms[2])
                if len(terms) > 1:
                    if terms[0] in cls.ncbi_tags and cls.accession_pattern.match(terms[1]):
                        return terms[1]

            # check for new style header
            # (see https://www.ncbi.nlm.nih.gov/news/09-17-2014-simple-FASTA-headers-genomes-FTP)
            if cls.accession_pattern.match(header):
                return header

        return None

    def identifier(self) -> str:
        """
        Parse target ids and target definition for an identifier that can be used for Entrez search

        :return: an identifier for Entrez search
        """

        acc = self.fasta_header_identifier(self.target_id)
        if acc is not None:
            return acc
        acc = self.fasta_header_identifier(self.target_def)
        if acc is not None:
            return acc

        # this should not happen for public data or databases derived from public data
        print(
            "Unable to find Entrez identifier in hit target id {} or definition {}".format(self.target_id,
                                                                                           self.target_def))

        return None

    def data_table_name(self) -> str:
        """
        Used when the blast hit is the result of searching a blast database created from sequences extracted from a
        data table column

        :return: Sequence name
        """
        name = self.target_def
        match = re.match(r'(\S+) ', name)
        if match is not None:
            name = match.group(1)
        return name

    def complement_target_sequence(self) -> None:
        """
        Reverse complements the target record sequence if the match is on the other strand (blastn only)

        Determination of reverse complement match handles ambiguous DNA in the query, but not on the target.
        :class:`SequenceMatcher` has an exhaustive matcher, but this is too slow for use on larger target sequences
        """

        if self.search_type != BlastSearchType.BLASTN:
            return

        if self.target_record is None:
            return

        seq = self.target_record.seq
        if len(seq) == 0:
            return
        if not is_defined_sequence(seq):
            return

        # hide molecule type while creating reverse complement
        # Entrez records may label sequences as mRNA, whereas the sequence is DNA
        # on parsing Biopython sets the molecule type to RNA and the reverse complement barfs when is gets DNA
        molecule_type = None
        if 'molecule_type' in self.target_record.annotations:
            molecule_type = self.target_record.annotations['molecule_type']
            del self.target_record.annotations['molecule_type']

        # should really use the SequenceMatcher class from blast_utils, but it's very slow
        # nt_search will not catch ambiguous dna matches on the target reference sequence
        if len(nt_search(str(seq.reverse_complement()), self.target.replace('-', ''))) > 1:
            self.target_record = self.target_record.reverse_complement(id=True, name=True, description=True,
                                                                       features=True, annotations=True,
                                                                       letter_annotations=True, dbxrefs=True)
        if molecule_type:
            self.target_record.annotations['molecule_type'] = molecule_type;

    def truncate_target_record(self, maximum_length: int) -> None:
        assert maximum_length % 2 == 0
        if self.target_record is None:
            return
        if len(self.target_record) <= maximum_length:
            return
        target_match = self.target.replace('-', '')
        if self.search_type == BlastSearchType.BLASTN:
            match = SequenceMatcher.match_sequences(str(self.target_record.seq), target_match)
        else:
            match = str(self.target_record.seq).find(target_match)
        if match == -1:
            self.target_record = self.target_record[0: maximum_length]
        else:
            mid = match + (int)(len(target_match) / 2)
            half = (int)(maximum_length / 2)
            start = mid - half
            end = mid + half
            if start < 0:
                start = 0
                end = maximum_length
            elif end > len(self.target_record):
                end = len(self.target_record)
                start = len(self.target_record) - maximum_length
            self.target_record = self.target_record[start: end]
            assert len(self.target_record) == maximum_length


class MultipleBlastResults(Frozen):
    """
    A class to store the results of multiple query searches against a single Blast database

    Attributes:
        * query_hits: a list of :class:`BlastResults`, with one result class for each query
        * database: the name of the database searched
        * search_type: a :class:`BlastSearchType` describing the blast search parameters
    """

    def __init__(self):
        """
        Empty constructor
        """
        self.query_hits = []  # type: List[BlastResults]
        self.database = None  # type: str
        self.search_type = None  # type: BlastSearchType

    def parse(self, file: str) -> None:
        """
        Uses Biopython's :class:`Bio.Blast.NCBIXML` parser to extract relevant information from a Blast search's XML output.  Each
        query's search is added as a :class:`BlastResults`.

        :param file: inout handle to blast results in XML format
        """
        with open(file, 'r') as fh:
            for blast_record in NCBIXML.parse(fh):
                results = BlastResults()
                results.read_record(blast_record)
                self.query_hits.append(results)
                self.database = results.database
                self.search_type = results.search_type

    def retrieve_targets(self, retrieve_missing_locally: bool = False) -> None:
        """
        Retrieves target records from NCBI Entrez non-redundant databases

        :param retrieve_missing_locally: if True extract missing full length target sequences from the local database (assumes search is local)
        """
        """Attempt to download all target sequences from Entrez"""
        target_ids = [hit.identifier() for results in self.query_hits for hit in results.hits if
                      hit.target_record is None]
        db = 'Nucleotide' if self.search_type.database_type == 'nucl' else 'Protein'
        records = retrieve_entrez_records(db, target_ids)
        record_map = {record.id: record for record in records}
        missing = False
        for results in self.query_hits:
            for hit in results.hits:
                if hit.target_record is None:
                    target_id = hit.identifier()
                    if target_id in record_map and is_defined_sequence(record_map[target_id].seq):
                        hit.target_record = record_map[target_id]
                    else:
                        print("Failed to find matching record for hit {}".format(hit.target_id))
                        missing = True
        if missing and retrieve_missing_locally:
            self.retrieve_local_targets()

    def retrieve_local_targets(self):
        """
        Retrieves all full-length target sequences from the local blast or
        emboss database (assumes search is local)
        :return:
        """

        local_records = BlastRecords()
        target_ids = [hit.target_id for results in self.query_hits for hit in results.hits if
                      hit.target_record is None]
        local_records.retrieve_from_local_database(self.database, target_ids)
        for results in self.query_hits:
            for hit in results.hits:
                if hit.target_record is None:
                    if hit.target_id in local_records.record_map:
                        hit.target_record = local_records.record_map[hit.target_id].record
                    else:
                        print("Failed to find matching local record for hit {}".format(hit.target_id))

    def complement_target_sequence(self) -> None:
        """
        For blastn searches, if it looks like the match is on the other strand, reverse complements the matching
        target record, see :func:`BlastResult.complement_target_sequence`
        """

        if self.search_type != BlastSearchType.BLASTN:
            return
        for results in self.query_hits:
            for hit in results.hits:
                if hit.target_record is not None:
                    hit.complement_target_sequence()

    def truncate_target_sequences(self, maximum_length: int = 10_000) -> None:
        """
        Truncate full length target sequences

        :param maximum_length: Maxumum target length
        """
        for results in self.query_hits:
            for hit in results.hits:
                if hit.target_record is not None and len(hit.target_record) > maximum_length:
                    hit.truncate_target_record(maximum_length)

    def number_hits(self):
        return len([hit for results in self.query_hits for hit in results.hits])


class BlastResults(Frozen):
    """
    A class to hold results for a query.  A container for a list of :class:`BlastResult`

    Attributes:
        * query_id: query identifier
        * query_def: query definition
        * database: the name of the database searched
        * search_type: a :class:`ruse.bio.blast_search.BlastSearchType` describing the blast search parameters
        * hits: List of associated :class:`BlastResult` results

    """

    def __init__(self):
        """
        No argument constructor
        """
        self.hits = []  # type: List[BlastResult]
        self.database = None  # type: str
        self.query_id = None  # type: str
        self.query_def = None  # type: str
        self.search_type = None  # type: BlastSearchType

    def parse(self, file: str) -> None:
        """
        Parses query search results from a blast result file.  There should only  be one result in the file

        :param file: file-like handle to blast search results in XML format
        """
        with open(file, 'r') as fh:
            blast_record = NCBIXML.read(fh)
            self.read_record(blast_record)

    def read_record(self, blast_record) -> None:
        """
        Processes a blast search record from Biopython extracting target hits into a list of  class:`BlastResult` results

        :param blast_record: record returned from Biopython :func:`Bio.Blast.NCBIXML.read` or :func:`Bio.Blast.NCBIXML.parse`

        """
        self.database = blast_record.database
        self.hits = []
        self.query_id = blast_record.query_id
        self.query_def = blast_record.query
        self.search_type = BlastSearchType.from_string(blast_record.application)
        for alignment in blast_record.alignments:
            hsp = alignment.hsps[0]
            attr = {'target_id': alignment.hit_id, 'target_def': alignment.hit_def,
                    'target_length': alignment.length, 'align_length': hsp.align_length,
                    'query': hsp.query, 'target': hsp.sbjct, 'evalue': hsp.expect,
                    'score': int(round(hsp.score)), 'bits': hsp.bits, 'accession': alignment.accession,
                    'query_start': hsp.query_start-1,
                    'query_end': hsp.query_end-1,
                    'search_type': self.search_type}
            result = BlastResult(**attr)
            self.hits.append(result)

    def to_mapping(self) -> Dict[str, BlastResult]:
        """
        Maps target_ids to results

        :return: a map of target identifiers to  Blast class:`BlastResult` results
        """
        return {hit['target_id']: hit for hit in self.hits}

    def to_data(self) -> Dict[str, List[Dict[str, DataTypes]]]:
        """
        Maps target ids to dictionaries of results

        :return: a map of target identifiers to  Blast class:`BlastResult` results
        """

        results_data = [r.__dict__ for r in self.hits]
        return {'hits': results_data}

    def retrieve_targets(self) -> None:
        """
        Retrieve from NCBI entrez genbank records for all targets
        """
        db = 'Nucleotide' if self.search_type.database_type == 'nucl' else 'Protein'
        target_ids = [h.identifier() for h in self.hits]
        target_ids = [t for t in target_ids if t is not None]
        records = retrieve_entrez_records(db, target_ids)
        record_map = {record.id: record for record in records}
        for hit in self.hits:
            target_id = hit.identifier()
            if target_id is not None:
                if target_id in record_map:
                    hit.target_record = record_map[target_id]
                else:
                    print("Failed to find matching record for hit {}".format(hit.target_id))

    def retrieve_local_targets(self) -> None:
        """
        Retrieve from a local database full length sequences for all targets without full length records
        """
        local_records = BlastRecords()
        target_ids = [h.target_id for h in self.hits if h.target_record is None]
        local_records.retrieve_from_local_database(self.database, target_ids)
        for hit in self.hits:
            if hit.target_record:
                pass
            elif hit.target_id in local_records.record_map:
                hit.target_record = local_records.record_map[hit.target_id].record
            else:
                print("Failed to find matching record for hit {}".format(hit.target_id))

    def complement_target_sequence(self) -> None:
        """
        For BLASTN searches replace the target sequence by the complement if it looks like the match is on the
        other strand
        :return:
        """
        if self.search_type != BlastSearchType.BLASTN:
            return
        for hit in self.hits:
            if hit.target_record:
                hit.complement_target_sequence()


def build_common_alignments(query: SeqRecord, hits: List[Optional[Union[BlastResult, IgHit]]]) -> List[SeqRecord]:

    # create alignment from hit
    def build_hit_sequence(hit: Optional[Union[BlastResult, IgHit]]) -> Optional[SeqRecord]:
        if not hit:
            return None
        buf = StringIO()
        hit_pos = 0
        q_seq = hit.query
        q_len = len(q_seq)
        t_seq = hit.target
        for i, res in enumerate(str(query_seq)):
            query_gaps = gaps[i] if i in gaps else 0
            if i < hit.query_start or i > hit.query_end:
                buf.write('-')
                for _ in range(query_gaps):
                    buf.write('-')
            else:
                query_gaps_processed = 0
                assert res == q_seq[hit_pos]
                buf.write(t_seq[hit_pos])
                hit_pos += 1
                while hit_pos < q_len and q_seq[hit_pos] == '-':
                    query_gaps_processed += 1
                    buf.write(t_seq[hit_pos])
                    hit_pos += 1
                assert query_gaps_processed <= query_gaps
                for _ in range(query_gaps-query_gaps_processed):
                    buf.write('-')
        hit_align = buf.getvalue()
        assert len(hit_align) == len(query_align)
        id = hit.identifier()
        if not id:
            id = hit.target_id
        hit_record = SeqRecord(Seq(hit_align), id, hit.target_id, hit.target_def)
        return hit_record

    # first find all gaps in query matches
    query_seq = query.seq.ungap()
    gaps = {}
    for hit in hits:
        if not hit:
            continue
        ungapped_pos = hit.query_start-1

        in_gap = False
        for i, base in enumerate(hit.query):
            if base == '-':
                if not in_gap:
                    gap_length = 1
                    in_gap = True
                else:
                    gap_length += 1
            else:
                if in_gap:
                    if ungapped_pos in gaps:
                        if gap_length > gaps[ungapped_pos]:
                            gaps[ungapped_pos] = gap_length
                    else:
                        gaps[ungapped_pos] = gap_length
                    in_gap = False
                ungapped_pos += 1
                assert query_seq[ungapped_pos] == base

    buf = StringIO()
    for i, res in enumerate(str(query_seq)):
        buf.write(res)
        query_gaps = gaps[i] if i in gaps else 0
        for _ in range(query_gaps):
            buf.write('-')

    # now build alignments
    query_align_seq = Seq(buf.getvalue())
    query_align = SeqRecord(query_align_seq, query.id, query.name, query.description, query.dbxrefs)
    copy_features_to_gapped_sequence(query, query_align)
    aligned_hits = [build_hit_sequence(h) for h in hits]
    aligned_hits.insert(0, query_align)

    return aligned_hits




