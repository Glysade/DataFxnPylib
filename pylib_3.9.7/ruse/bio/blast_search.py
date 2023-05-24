"""
===============
blast_search.py
===============

Copyright (C) 2017-2022 Glysade, LLC

Classes and functions for performing Blast searches against local and NCBI-hosted databases
"""

import os
import subprocess
import uuid
from enum import Enum
from typing import List, Dict, Union, NamedTuple, Optional

from Bio import Entrez
from Bio import SeqIO
from Bio.Application import ApplicationError
from Bio.Blast import NCBIWWW
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.SeqRecord import SeqRecord

from ruse.bio.bio_util import entrez_email, is_dna, ok_sequence
from ruse.util.frozen import Frozen

Entrez.email = entrez_email()
ncbi_blast_version = '2.10.1+'


def full_path_to_blast_exe(exe: str) -> str:
    """
    Finds the full path to an NCBI executable

    :param exe: the name of the program
    :return: full path to the program
    """
    base = os.path.normpath("../../bin")
    if os.name == 'nt':
        if not exe.lower().endswith('.exe'):
            exe += '.exe'
        return exe
    elif os.name == 'posix':
        if exe.lower().endswith('.exe'):
            exe = exe[:-4]
        return exe
    else:
        raise ValueError(f'Unknown OS {os.name}')


class BlastDatabase(Frozen):
    """
    A class to handle creation of blast databases

    Attributes:
        * name: the name of the database
        * sequence_count: the number of sequences in the database

    """

    def __init__(self, name: str = None):
        """
        Constructor

        :param name: the name of the database

        """

        self.name = name
        self.sequence_count = None  # type: int

    def build_database(self, sequences: List[SeqRecord], database_name: str = None, type: str = None) -> str:
        """
        Constructs a local database from the input sequences

        :param sequences: a list of :class:`Bio.SeqRecord.SeqRecord` sequence records
        :param database_name: the name of the database.  If None it will be created from a UUID
        :param type: The database type (prot or ncul).  If None it will be set from the sequences
        :return: the name of the database
        """

        sequences = [SeqRecord(s.seq.ungap(), s.id, s.name, s.description) for s in sequences if ok_sequence(s)]
        fasta_file = "{}.fasta".format(str(uuid.uuid4()))
        with open(fasta_file, 'w') as f:
            SeqIO.write(sequences, f, "fasta")

        if type is None:
            type = 'nucl' if is_dna(str(sequences[0].seq)) else 'prot'
        if type not in ['nucl', 'prot']: raise ValueError

        if database_name is None:
            database_name = self.name
        if database_name is None:
            database_name = "{}_{}.db".format(str(uuid.uuid4()), type)

        # print("type {} database name {}".format(type, database_name))

        formatdb_exe = full_path_to_blast_exe('makeblastdb')
        args = [formatdb_exe, '-dbtype', type, '-in', fasta_file, '-out', database_name]
        # print("creating blast database {} args {}".format(database_name, args))
        subprocess.run(args, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        os.remove(fasta_file)

        self.name = database_name
        self.sequence_count = len(sequences)

        return database_name

    def clean_database(self) -> None:
        """
        Removes all database files
        """

        for path in ["{}.phr", "{}.pin", "{}.psq", "{}.nhr", "{}.nin", "{}.nsq", "{}.pdb", "{}.pot", "{}.ptf",
                     "{}.pto", "{}.ndb", "{}.not", "{}.ntf", "{}.nto"]:
            path = path.format(self.name)
            if os.path.exists(path):
                os.remove(path)


class BlastSearchType(Enum):
    """
    A Enum class for the different types of Blast searches

    Members:
        * BLASTN: compare nucl query with nucl
        * TBLASTX: compare nucl query with nucl (using translations/both strands to protein sequences)
        * TBLASTN: compare protein query (using translations/both strands) with nucl
        * BLASTP: compare protein with protein
        * BLASTX: compare nucl query (using translations/both strands) with protein
    """

    # compare nucl query with nucl
    BLASTN = 1
    # compare nucl query with nucl (using translations/both strands to protein sequences)
    TBLASTX = 2
    # compare protein query (using translations/both strands) with nucl
    TBLASTN = 3
    # compare protein with protein
    BLASTP = 4
    # compare nucl query (using translations/both strands) with protein
    BLASTX = 5

    def query_type(self) -> str:
        """
        The expected query type for this search

        :return: "nucl" or "prot"
        """

        if self in [self.BLASTN, self.BLASTX, self.TBLASTX]:
            return 'nucl'
        else:
            return 'prot'

    def database_type(self) -> str:
        """
        The expected database type for this search

        :return: "nucl" or "prot"
        """
        if self in [self.BLASTN, self.TBLASTX, self.TBLASTN]:
            return 'nucl'
        else:
            return 'prot'

    def exe(self) -> str:
        """
        Determine full patch to command line search program

        :return: full path to executable
        """
        return full_path_to_blast_exe(self.name.lower())

    @classmethod
    def from_string(cls, value) -> 'BlastSearchType':
        """
        Generates a search type enum from the string search name

        :param value: The string name fo a search type
        :return: The equivalent Enum class
        """
        return cls.__members__[value.upper()]

    @classmethod
    def default_search_type(cls, query_type: str, target_type: str) -> 'BlastSearchType':
        """
        Returns default search type from query and target sequence type

        :param query_type: one of 'nucl' or 'prot'
        :param target_type: one of 'nucl' or 'prot'
        :return: Blast search type
        """

        assert target_type in ['nucl', 'prot']
        assert query_type in ['nucl', 'prot']

        if query_type == 'prot':
            if target_type == 'prot':
                return BlastSearchType.BLASTP
            else:
                return BlastSearchType.TBLASTN
        if query_type == 'nucl':
            if target_type == 'nucl':
                return BlastSearchType.BLASTN
            else:
                return BlastSearchType.BLASTX

        assert False


OptionDataTypes = Union[str, int, float]
SearchType = Union[BlastSearchType, str]


class BlastSearch(Frozen):
    """
    A class to manage a blast search of multiple queries against a local database

    Attributes:
        * queries: A list of :class:`Bio.SeqRecord.SeqRecord` query sequences
        * database_name: The name of the target database
        * search_type: Search type of :class:`BlastSearchType` value
        * options: Dict of allowable options (see :func:`BlastSearch.multiple_query_search_blast_database`)
    """

    def __init__(self, ):
        """
        Empty constructor
        """

        self.queries = None  # type: List[SeqRecord]
        self.database_name = None  # type: str
        self.search_name = None  # type: str
        self.search_type = None  # type: BlastSearchType
        self.options = {}  # type: Dict[str, OptionDataTypes]
        self.error = None  # type: str

    def search_blast_database(self, query: SeqRecord, database_name: str, search_type: SearchType,
                              search_name: str = None,
                              options: Dict[str, OptionDataTypes] = {}) -> str:
        """
        Search a single query against a blast database.  Delegates to :func:`multiple_query_search_blast_database`


        :param query:  A :class:`Bio.SeqRecord.SeqRecord` query sequence
        :param database_name: The name of the database to search
        :param search_type: Search type of :class:`BlastSearchType` value
        :param search_name: The name of the search (if none a UUID will be used to assign ont)
        :param options: Optional dict of allowable options (see :func:`BlastSearch.multiple_query_search_blast_database`).
        :return: The search name
        """

        # make some reasonable defaults
        if 'max_target_seqs' not in options and 'num_alignments' not in options and 'num_descriptions' not in options:
            options['max_target_seqs'] = 250
        return self.multiple_query_search_blast_database([query], database_name, search_type, search_name, options)

    def multiple_query_search_blast_database(self, queries: List[SeqRecord], database_name: str,
                                             search_type: SearchType,
                                             search_name: str = None,
                                             options: Dict[str, OptionDataTypes] = {}):
        """
        Searches multiple queries against a local Blast database, uses :func"`Bio.Blast.Applications.NcbiblastxCommandline`

        Options can be any option available to the particular blast program.  Default options are:
            num_threads: number of processors
            max_target_seqs: 10
            max_hsps: 1 (multiple hsps are not processed by the parse classes)


        :param queries:  A list of :class:`Bio.SeqRecord.SeqRecord` query sequences
        :param database_name: The name of the database to search
        :param search_type: Search type of :class:`BlastSearchType` value
        :param search_name: The name of the search (if none a UUID will be used to assign ont)
        :param options: Optional dict of allowable options as option name, value pairs
        :return: The search name
        """

        if search_name is None:
            search_name = uuid.uuid4()
        self.search_name = search_name

        fasta_file = "{}.fasta".format(search_name)
        blast_queries = [SeqRecord(q.seq.ungap(), q.id, q.name, q.description) for q in queries]
        with open(fasta_file, 'w') as f:
            SeqIO.write(blast_queries, f, "fasta")
        output_file = self.output_file()

        # local import to prevent circular errors
        from ruse.util.util import num_processors
        args = {'cmd': search_type.exe(), 'outfmt': 5, 'db': database_name, 'query': fasta_file, 'out': output_file,
                'num_threads': num_processors()}
        # in options max_target_seqs cannot be used with num_alignments or num_descriptions
        if 'max_target_seqs' in options:
            if 'num_alignments' in options or 'num_descriptions' in options:
                raise ValueError("Blast option max_target_seqs cannot be used with num_alignments or num_descriptions ")

        # make some reasonable defaults, so as to keep down number of hits
        if 'max_target_seqs' not in options and 'num_alignments' not in options and 'num_descriptions' not in options:
            options['max_target_seqs'] = 10
        if 'max_hsps' not in options:
            options['max_hsps'] = 1
        args.update(options)

        cmd = NcbiblastxCommandline(**args)
        # print("NCBI command line is {}".format(cmd))
        try:
            cmd()
        except ApplicationError as ex:
            self.error = str(ex)
            print("Exception {} on blast search".format(ex))

        self.queries = queries
        self.database_name = database_name
        self.search_type = search_type
        self.options = options

        return search_name

    def output_file(self) -> str:
        """

        :return: the name of the xml output file
        """

        return "{}.xml".format(self.search_name)

    def clean_search(self) -> None:
        """
        Remove all search results
        """

        if os.path.exists(self.output_file()):
            os.remove(self.output_file())
        fasta_file = "{}.fasta".format(self.search_name)
        if os.path.exists(fasta_file):
            os.remove(fasta_file)


class BlastCreateAndSearch(Frozen):
    """
    A class that creates a blast database on the fly from supplied sequences and performs a single search against that
    database

    Attributes:
        * query: A :class:`Bio.SeqRecord.SeqRecord` query sequence
        * search_name: the name for the search
        * search_type:  Search type of :class:`BlastSearchType` value
        * options: Blast search options as detailed in :func:`BlastSearch.multiple_query_search_blast_database`
        * search: A :class:`BlastSearch` object that runs the local search
        * database: The name for the database
        * target_sequences: A list of :class:`Bio.SeqRecord.SeqRecord` target sequences
        * query_type: The query sequence type. Either 'prot' or 'nucl'
        * target_type: The target sequence type. Either 'prot' or 'nucl'

    """

    def __init__(self):
        """
        Empty constructor
        """
        self.query = None  # type: SeqRecord
        self.search_name = None  # type: str
        self.search_type = None  # type: BlastSearchType
        self.options = {}  # type: Dict[str, OptionDataTypes]
        self.search = None  # type: BlastSearch
        self.database = None  # type: BlastDatabase
        self.target_sequences = None  # type: List[SeqRecord]
        self.query_type = None  # type: str
        self.target_type = None  # type: str

    def search_blast_sequences(self, query: SeqRecord, sequences: List[SeqRecord],
                               search_type: Optional[SearchType] = None,
                               search_name: Optional[str] = None, query_type: Optional[str] = None,
                               options: Dict[str, OptionDataTypes] = {}) -> str:
        """
        Creates a blast database from the target sequences, then performs blast search using the query sequence

        :param query: The :class:`Bio.SeqRecord.SeqRecord` query sequence
        :param sequences: A list of :class:`Bio.SeqRecord.SeqRecord` target sequences
        :param search_type: Optional search type of :class:`BlastSearchType` value.  Will guess from sequences if not present
        :param search_name: Optional search name.  Assigned using UUID if not present
        :param options: Blast search options as detailed in :func:`BlastSearch.multiple_query_search_blast_database`
        :return:
        """
        if search_type is not None and type(search_type) is not BlastSearchType:
            search_type = BlastSearchType.from_string(search_type)

        if search_type is None:
            target_type = 'nucl' if is_dna(str(sequences[0].seq)) else 'prot'
            if query_type is None:
                query_type = 'nucl' if is_dna(str(query.seq)) else 'prot'
            search_type = BlastSearchType.default_search_type(query_type, target_type)
        else:
            target_type = search_type.database_type()
            query_type = search_type.query_type()

        if search_name is None:
            search_name = str(uuid.uuid4())

        self.search_name = search_name
        self.search_type = search_type
        self.query = query
        self.target_sequences = sequences
        self.query_type = query_type
        self.target_type = target_type

        # print("Building blast database of {} {} sequences and performing {} search using {} query".
        #      format(len(sequences), target_type, search_type, query_type))

        database = BlastDatabase()
        self.database = database
        database.build_database(sequences, search_name, target_type)

        search = BlastSearch()
        self.search = search
        search.search_blast_database(query, search_name, search_type, search_name, options)

        return search_name

    def clean_search(self) -> None:
        """
        Remove all search and database files
        """
        if self.database:
            self.database.clean_database()
        if self.search:
            self.search.clean_search();


class BlastWebDatabase(NamedTuple):
    name: str
    protein: bool


class BlastWebSearch(Frozen):
    """
    A class to do a BLAST search using the QBLAST server at NCBI or a cloud service provider, using :func:`Bio.Blast.NCBIWWW`

    Attributes:
        * query: The :class:`Bio.SeqRecord.SeqRecord` query sequence
        * database_name: The remote database to search
        * search_name: Search name.  Assigned using UUID if not present
        * search_type: Search type of :class:`BlastSearchType` value.
        * options: Blast search options
    """

    known_databases: List[BlastWebDatabase] = {
        BlastWebDatabase('nr', True),
        BlastWebDatabase('refseq_select_prot', True),
        BlastWebDatabase('refseq_protein', True),
        BlastWebDatabase('swissprot', True),
        BlastWebDatabase('pataa', True),
        BlastWebDatabase('pdb', True),
        BlastWebDatabase('env_nr', True),
        BlastWebDatabase('tsa_nr', True),
        BlastWebDatabase('nt', False),
        BlastWebDatabase('refseq_select_nucl', False),
        BlastWebDatabase('refseq_rna', False),
        BlastWebDatabase('est', False),
        BlastWebDatabase('pat', False),
        BlastWebDatabase('gss', False),
        BlastWebDatabase('dbsts', False),
        BlastWebDatabase('env_nt', False)
    }

    def __init__(self):
        """
        Empty constructor
        """
        self.query = None  # type: SeqRecord
        self.database_name = None  # type: str
        self.search_name = None  # type: str
        self.search_type = None  # type: BlastSearchType
        self.options = {}  # type: Dict[str, OptionDataTypes]

    def search_blast_database(self, query: SeqRecord, database_name: str = None, search_type: SearchType = None,
                              search_name: str = None, hitlist_size: int = 100, query_type: str = None,
                              expect: Optional[float] = None,
                              word_size: Optional[int] = None,
                              options: Dict[str, OptionDataTypes] = {}) -> str:
        """
        Performs a remote blast search against the NCBI QBLAST server

        :param query: The :class:`Bio.SeqRecord.SeqRecord` query sequence
        :param database_name: The remote database to search. If not present use 'nt' for nucleotide query sequence and 'nr' for a protein query
        :param search_type: Optional search type of :class:`BlastSearchType` value.  Will guess from sequences and database if not present
        :param query_type: Optional query type- should be one of 'nucl' or 'prot', if not present guess from sequence
        :param search_name: Optional search name.  Assigned using UUID if not present
        :param options: Dictionary of Qblast search options, passed to :func:`Bio.Blast.NCBIWWW`
        :return: search name
        """

        if query_type is None:
            query_type = 'nucl' if is_dna(str(query.seq)) else 'prot'
        if database_name is None:
            database_name = 'nt' if query_type == 'nucl' else 'nr'
        databases = [d for d in self.known_databases if d.name == database_name]
        if len(databases) != 1:
            raise ValueError("Unknown web blast database {}".format(database_name))
        database = databases[0]

        # refseq_select can be specified for either protein or nucleotide searches
        # suffix was added above and in .yaml to differentiate between selections
        # remove suffix and create appropriate BlastWebDatabase object
        if 'refseq_select' in database.name:
            database = BlastWebDatabase('refseq_select', database.protein)

        if search_type is None:
            target_type = 'prot' if database.protein else 'nucl'
            search_type = BlastSearchType.default_search_type(query_type, target_type)
        if search_name is None:
            search_name = str(uuid.uuid4())
        self.search_name = search_name

        args = {'program': search_type.name.lower(), 'database': database_name, 'sequence': query.seq,
                'hitlist_size': hitlist_size, 'expect': expect}
        if expect is not None:
            args['expect'] = expect
        if word_size is not None:
            args['word_size'] = word_size
        args.update(options)

        print("Performing web blast search {}".format(search_name))
        with NCBIWWW.qblast(**args) as rh:
            with open(self.output_file(), 'w') as fh:
                record = rh.read()
                fh.write(record)

        self.query = query
        self.database_name = database_name
        self.search_type = search_type
        self.search_name = search_name
        self.options = options

        return search_name

    def output_file(self) -> str:
        """

        :return: the xml result file
        """

        return "{}.xml".format(self.search_name)

    def clean_search(self) -> None:
        """
        Cleans search results from hard drive
        """
        os.remove(self.output_file())
