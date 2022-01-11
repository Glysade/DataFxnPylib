from Bio.Application import _Option, _Switch, AbstractCommandline

import numbers

from Bio.Emboss.Applications import _EmbossCommandLine, _EmbossMinimalCommandLine, SeqretCommandline


class AnarciCommandLine(AbstractCommandline):

    def __init__(self, cmd: str = 'ANARCI', **kwargs: str) -> None:
        numbering_schemes = ['kabat', 'aho', 'wolfguy', 'imgt', 'a', 'c', 'chothia', 'i', 'k', 'm', 'w', 'martin']
        restrict_chains = ['ig', 'tr', 'heavy', 'light', 'H', 'K', 'L', 'A', 'B']
        species = ['alpaca', 'rabbit', 'rhesus', 'pig', 'rat', 'human', 'mouse']

        self.parameters = [
            _Option(['--sequence', '-i', 'sequence'],
                    'A sequence or an input fasta file',
                    filename=True,
                    equate=False),
            _Option(['--outfile', '-o', 'outfile'],
                    'The output file to use. Default is stdout',
                    filename=True,
                    equate=False),
            _Option(['--scheme', '-s', 'scheme'],
                    'Which numbering scheme should be used',
                    checker_function=lambda n: n in numbering_schemes,
                    equate=False
                    ),
            _Option(['--restrict', '-r', 'restrict'],
                    'Restrict ANARCI to only recognise certain types of receptor chains',
                    checker_function=lambda c: c in restrict_chains,
                    equate=False),
            _Switch(['--csv', 'csv'],
                    ' Write the output in csv format'),
            _Option(['--outfile_hits', '-ht', 'outfile_hits'],
                    'Output file for domain hit tables for each sequence',
                    filename=True,
                    equate=False),
            _Option(['--hmmerpath', '-hp', 'hmmerpath'],
                    'The path to the directory containing hmmer programs',
                    filename=True,
                    equate=False),
            _Option(['--ncpu', '-p', 'ncpu'],
                    'Number of parallel processes to use',
                    equate=False),
            _Switch(['--assign_germline', 'assign_germline'],
                    'Assign the v and j germlines to the sequence.'),
            _Option(['--use_species', 'use_species'],
                    'Use a specific species in the germline assignment',
                    checker_function=lambda s: s in species,
                    equate=False),
            _Option(['--bit_score_threshold', 'bit_score_threshold'],
                    'Change the bit score threshold',
                    checker_function=lambda b: isinstance(b, numbers.Integral),
                    equate=False)]
        super().__init__(cmd, **kwargs)


class RestrictCommandLine(_EmbossCommandLine):

    def __init__(self, cmd='restrict', **kwargs: str) -> None:
        self.parameters = [
            _Option(['-sequence', 'sequence'],
                    'Nucleotide sequence(s) filename',
                    filename=True,
                    is_required=True),
            _Option(['-enzymes', 'enzymes'],
                    'Comma separated list of enzymes or all'),
            _Option(['-sitelen', 'sitelen'],
                    'This sets the minimum length of the restriction enzyme recognition site (default 4)'),
            _Option(["-rformat", "rformat"],
                    "Specify the report format to output in.")
        ]
        super().__init__(cmd, **kwargs)


# Extends Biopython Seqret command line to use the feature option
class SeqretFeatureCommandline(SeqretCommandline):

    def __init__(self, cmd="seqret", **kwargs):
        super().__init__(cmd, **kwargs)

    def __getattribute__(self, item):
        if item == 'parameters' and 'feature_switch_enabled' not in self.__dict__:
            self.__dict__['feature_switch_enabled'] = True
            feature_switch = _Switch(["-feature", "feature"], "Process features"),
            self.__dict__['parameters'].append(feature_switch[0])
        return super().__getattribute__(item)


# Simple wrapper round showdb that only handles a few options
class ShowDbCommandLine(_EmbossCommandLine):

      def __init__(self, cmd='showdb', **kwargs: str) -> None:
        self.parameters = [
            _Switch(['-protein', 'protein'],
                    "Display protein databases"),
            _Switch(['-nucleic', 'nucleic'],
                    "Display nucleotide databases"),
            _Switch(['-type', 'type'],
                    "Display 'type' column"),
        ]
        super().__init__(cmd, **kwargs)
