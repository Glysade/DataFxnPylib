import collections
import json
import os
import re
import uuid
from typing import List, NamedTuple, Optional, Dict, Final

from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from ruse.bio.applications import AnarciCommandLine
from ruse.bio.bio_util import sequences_to_file
from ruse.bio.sequence_align import copy_features_to_gapped_sequence

ANTIBODY_NUMBERING_COLUMN_PROPERTY: Final[str] = "antibodyNumbering"


class AntibodyNumber(NamedTuple):
    domain: int
    chain: str
    position: int
    insertion: Optional[str]
    query_position: Optional[int]

    def label(self):
        insertion = self.insertion if self.insertion else ''
        return '{}{}{}'.format(self.chain, self.position, insertion)

    def __lt__(self, other: 'AntibodyNumber') -> bool:
        if self.domain != other.domain:
            return self.domain < other.domain
        if self.chain != other.chain:
            return self.chain == 'L'
        if self.position != other.position:
            return self.position < other.position
        if self.insertion is None and other.insertion:
            return True
        if self.insertion and other.insertion:
            return self.insertion < other.insertion
        return False

    def __eq__(self, other: 'AntibodyNumber') -> bool:
        return self.domain == other.domain and self.chain == other.chain and self.position == other.position \
               and self.insertion == other.insertion


class AntibodyNumberMapping(NamedTuple):
    domain: int
    chain: str
    position: int
    insertion: Optional[str]
    residue: str
    query_position: Optional[int]

    def label(self):
        insertion = self.insertion if self.insertion else ''
        return '{}{}{}'.format(self.chain, self.position, insertion)

    def to_antibody_number(self) -> AntibodyNumber:
        return AntibodyNumber(self.domain, self.chain, self.position, self.insertion, self.query_position)

    def matches(self, number: AntibodyNumber) -> bool:
        return self.chain == number.chain and self.position == number.position and self.insertion and number.insertion


class ChainRegions(NamedTuple):
    domain: int
    chain: str
    start: AntibodyNumber
    end: AntibodyNumber
    cdr1_start: AntibodyNumber
    cdr1_end: AntibodyNumber
    cdr2_start: AntibodyNumber
    cdr2_end: AntibodyNumber
    cdr3_start: AntibodyNumber
    cdr3_end: AntibodyNumber

    def to_data(self) -> Dict:
        start = self.start.query_position
        end = self.end.query_position
        start1 = self.cdr1_start.query_position
        end1 = self.cdr1_end.query_position
        start2 = self.cdr2_start.query_position
        end2 = self.cdr2_end.query_position
        start3 = self.cdr3_start.query_position
        end3 = self.cdr3_end.query_position
        chain = self.chain
        data = list()

        data.append({'name': '{}FR1'.format(chain), 'start': start, 'end': start1})
        data.append({'name': 'CDR-{}1'.format(chain), 'start': start1, 'end': end1 + 1})
        data.append({'name': '{}FR2'.format(chain), 'start': end1 + 1, 'end': start2})
        data.append({'name': 'CDR-{}2'.format(chain), 'start': start2, 'end': end2 + 1})
        data.append({'name': '{}FR3'.format(chain), 'start': end2 + 1, 'end': start3})
        data.append({'name': 'CDR-{}3'.format(chain), 'start': start3, 'end': end3 + 1})
        data.append({'name': '{}FR4'.format(chain), 'start': end3 + 1, 'end': end + 1})

        return {'domain': self.domain, 'chain': self.chain, 'regions': data}


class AntibodyAlignmentResult(NamedTuple):
    aligned_sequences: List[SeqRecord]
    numbering: List[AntibodyNumber]
    regions: List[ChainRegions]
    numbering_scheme: str
    cdr_definition: str

    def to_column_json(self) -> str:
        numbering_data = [{'domain': n.domain, 'position': n.query_position, 'label': n.label()} for n in
                          self.numbering]
        region_data = [r.to_data() for r in self.regions]
        antibody_numbering = {
            'scheme': self.cdr_definition,
            'numbering_scheme': self.numbering_scheme,
            'numbering': numbering_data,
            'regions': region_data
        }
        numbering_json = json.dumps(antibody_numbering)
        return numbering_json


def label_antibody_sequences(sequences: List[SeqRecord], numbering_scheme: str = 'kabat',
                             cdr_definition: str = 'kabat') -> List[List[AntibodyNumberMapping]]:
    mappings = _create_antibody_mappings(sequences, numbering_scheme)
    for mapping, sequence in zip(mappings, sequences):
        _annotate_sequence(sequence, mapping, numbering_scheme, cdr_definition)
    return mappings


def align_antibody_sequences(sequences: List[SeqRecord], numbering_scheme: str = 'kabat',
                             cdr_definition: str = 'kabat') -> AntibodyAlignmentResult:
    mappings = label_antibody_sequences(sequences, numbering_scheme, cdr_definition)
    alignments = _do_align_antibody_sequences(sequences, mappings, numbering_scheme, cdr_definition)
    return alignments


def _do_align_antibody_sequences(sequences: List[SeqRecord],
                                 mapping: List[List[AntibodyNumberMapping]], numbering_scheme: str,
                                 cdr_definition: str) -> AntibodyAlignmentResult:
    labellings = set()
    for seq_map in mapping:
        for pos in seq_map:
            labellings.add(AntibodyNumber(pos.domain, pos.chain, pos.position, pos.insertion, None))

    labellings = sorted(labellings)
    input_positions = [0] * len(sequences)
    aligned_seqs = [""] * len(sequences)
    labelling_positions = [0] * len(labellings)
    start = True
    position = 0
    seqs = [s.seq for s in sequences]

    for label_no, label in enumerate(labellings):

        sequence_labels = [_find_in_mapping(map, label, input_position) for map, input_position in
                           zip(mapping, input_positions)]
        if all(sl is None or sl.query_position is None for sl in sequence_labels):
            diff = 0
            label_position = position
        else:
            if start:
                start = False
                n_residues_before = [
                    next((pos.query_position for pos in seq_map if pos.query_position and pos.matches(label)), 0) for
                    seq_map in mapping]
                n_prefix = max(n_residues_before)
                position += n_prefix
                for seq_no, sequence in enumerate(seqs):
                    padding = n_prefix - n_residues_before[seq_no]
                    input_positions[seq_no] = n_residues_before[seq_no]
                    aligned_seqs[seq_no] = "-" * padding + sequence[0:n_residues_before[seq_no]] + aligned_seqs[seq_no]

            diff = max(l.query_position - o for l, o in zip(sequence_labels, input_positions) if
                       l and l.query_position is not None)
            label_position = position + diff

        for seq_no, sequence in enumerate(seqs):

            sequence_label = sequence_labels[seq_no]
            last = input_positions[seq_no]
            if sequence_label and sequence_label.query_position is not None:
                input_positions[seq_no] = sequence_label.query_position + 1
                n_insert = sequence_label.query_position - last
                for x in range(n_insert):
                    aligned_seqs[seq_no] += sequence[last + x]
                aligned_seqs[seq_no] += '-' * (diff - n_insert)
                aligned_seqs[seq_no] += sequence[sequence_label.query_position]
            else:
                aligned_seqs[seq_no] += '-' * (diff + 1)

            assert (len(aligned_seqs[seq_no]) == label_position + 1)

        position = label_position + 1
        labelling_positions[label_no] = label_position

    n_residues_after = max(len(s) - p for s, p in zip(seqs, input_positions)) if seqs else 0;
    for seq_no, sequence in enumerate(seqs):
        pos = input_positions[seq_no]
        padding = n_residues_after - len(sequence) + pos
        aligned_seqs[seq_no] += sequence[pos:]
        aligned_seqs[seq_no] += '-' * padding

    aligned_records = [SeqRecord(id=s.id, name=s.name, seq=seq, annotations={'molecule_type': 'protein'}) for s, seq in
                       zip(sequences, aligned_seqs)]
    for init, align in zip(sequences, aligned_records):
        copy_features_to_gapped_sequence(init, align)

    numbering = [AntibodyNumber(n.domain, n.chain, n.position, n.insertion, p) for p, n in
                 zip(labelling_positions, labellings)]
    regions = _find_all_regions(numbering, numbering_scheme, cdr_definition)

    return AntibodyAlignmentResult(aligned_records, numbering, regions, numbering_scheme, cdr_definition)


def _find_in_mapping(mapping: List[AntibodyNumberMapping], position: AntibodyNumber, input_position: int) -> Optional[
    AntibodyNumberMapping]:
    for m in mapping:
        if m.query_position is None:
            continue
        if m.query_position < input_position:
            continue
        antibody_number = m.to_antibody_number()
        if position < antibody_number:
            return None
        if antibody_number == position:
            return m
    return None


class AnarciDomain(NamedTuple):
    sequence_start: int
    sequence_end: int
    numbers: List[AntibodyNumberMapping]


def _create_antibody_mappings(sequences: List[SeqRecord], numbering_scheme: str) -> List[List[AntibodyNumberMapping]]:
    base = str(uuid.uuid4())
    in_file = 'seq_in_{}.fasta'.format(base)
    out_file = 'anarci_numbering_{}.txt'.format(base)
    sequences_to_file(in_file, sequences)

    lines: List[str] = []
    if sequences:
        command = AnarciCommandLine(cmd='ANARCI', sequence=in_file, scheme=numbering_scheme, outfile=out_file,
                                    restrict='ig')
        stdout, stderr = command()

        with open(out_file) as fh:
            lines = fh.readlines()

        for file in [in_file, out_file]:
            if os.path.exists(file):
                os.remove(file)

    mappings = []
    current_mapping = []
    current_domain = None
    domain_no = 0
    in_domain = False

    for line in lines:
        if line.startswith('//'):
            # print("start")
            mappings.append(current_mapping)
            current_domain = None
            current_mapping = []
        elif line.startswith('# Domain '):
            # print("Domain line {}"+line)
            terms = line.split(' ')
            domain_no = int(terms[2])
            in_domain = False
        elif line.startswith('#|') and not line.startswith('#|species'):
            # print("Info line {}"+line)
            items = line.split('|')
            current_domain = AnarciDomain(int(items[-3]), int(items[-2]), list())
            current_mapping.append(current_domain)
            assert not in_domain
            in_domain = True
        elif not line.startswith('#'):
            chain = line[0:1]
            assert chain == 'L' or chain == 'H'
            position = int(line[2:7])
            insertion = line[8:9].strip()
            if not insertion:
                insertion = None
            residue = line[10:11]
            current_domain.numbers.append(AntibodyNumberMapping(domain_no, chain, position, insertion, residue, None))

    assert len(mappings) == len(sequences)
    for mapping, record in zip(mappings, sequences):
        for domain in mapping:
            mapping_seq_arr = [m.residue for m in domain.numbers if m.residue != '-']
            mapping_seq = ''.join(mapping_seq_arr)
            assert mapping_seq in record.seq
            assert mapping_seq in record.seq[domain.sequence_start: domain.sequence_end + 1]

    mappings = [_match_to_sequence(s, m) for m, s in zip(mappings, sequences)]
    return mappings


def _match_to_sequence(record: SeqRecord, domains: List[AnarciDomain]) -> List[AntibodyNumberMapping]:
    domain_mappings = [_match_to_domain(record, d) for d in domains]
    return [m for dm in domain_mappings for m in dm]


def _match_to_domain(record: SeqRecord, domain: AnarciDomain) -> List[AntibodyNumberMapping]:
    mapping = domain.numbers
    match = ''.join([m.residue for m in mapping if m.residue != '-'])
    record_position = str(record.seq).find(match)
    assert record_position >= 0

    queue = collections.deque(mapping)
    mapping_to_record = []
    while queue:
        m = queue.popleft()
        if m.residue == '-':
            mapping_to_record.append(m)
        else:
            assert record[record_position] == m.residue
            mapping_to_record.append(
                AntibodyNumberMapping(m.domain, m.chain, m.position, m.insertion, m.residue, record_position))
            record_position += 1
    return mapping_to_record


def _find_all_regions(mapping: List[AntibodyNumber], mapping_scheme: str,
                      cdr_definition: str) -> List[ChainRegions]:
    def regions_in_domain(domain):
        filtered_mapping = [m for m in mapping if m.domain == domain]
        return _find_regions(filtered_mapping, mapping_scheme, cdr_definition)

    domains = set((m.domain for m in mapping))
    regions = [i for d in domains for i in regions_in_domain(d) if i]
    return regions


def _find_regions(mapping: List[AntibodyNumber], numbering_scheme: str,
                  cdr_definition: str) -> List[Optional[ChainRegions]]:
    # see CDR definitions in http://www.bioinf.org.uk/abs/info.html#cdrdef
    # This table is a little unclear as the end of H1 is not specified if the numbering is neither Kabat or Chothia
    # I have assumed that we use Chothia for everything except Kabat

    # numbering_scheme numbers residues
    # cdr_definition assigns regions

    # the Wolfguy numbering is not supported
    assert numbering_scheme in ['chothia', 'kabat', 'imgt', 'martin', 'aho']
    assert cdr_definition in ['chothia', 'kabat', 'imgt']
    regions = list()

    if cdr_definition == 'kabat':

        r = _find_chain_regions(mapping, 'L', 'L24', 'L34', 'L50', 'L56', 'L89', 'L97')
        regions.append(r)

        if numbering_scheme == 'kabat':
            h1_end = 'H35B'
        else:
            h1_end = 'H35'

        r = _find_chain_regions(mapping, 'H', 'H31', h1_end, 'H50', 'H65', 'H95', 'H102')
        regions.append(r)

    elif cdr_definition == 'chothia':

        r = _find_chain_regions(mapping, 'L', 'L24', 'L34', 'L50', 'L56', 'L89', 'L97')
        regions.append(r)

        if numbering_scheme == 'kabat':
            h35a = _find_antibody_number(mapping, 'H35A')
            h35a_present = h35a and h35a.position == 35 and h35a.insertion == 'A'
            h35b = _find_antibody_number(mapping, 'H35B')
            h35b_present = h35b and h35b.position == 35 and h35b.insertion == 'B'
            if not h35a_present and not h35b_present:
                h1_end = 'H32'
            elif h35a_present and h35b_present:
                h1_end = 'H34'
            elif h35a_present:
                h1_end = 'H33'
            else:
                raise ValueError()
        else:
            h1_end = 'H32'

        r = _find_chain_regions(mapping, 'H', 'H26', h1_end, 'H52', 'H56', 'H95', 'H102')
        regions.append(r)

    elif cdr_definition == 'imgt':

        r = _find_chain_regions(mapping, 'L', 'L27', 'L32', 'L50', 'L51', 'L89', 'L97')
        regions.append(r)

        if numbering_scheme == 'kabat':
            h1_end = 'H35B'
        else:
            h1_end = 'H33'
        r = _find_chain_regions(mapping, 'H', 'H26', h1_end, 'H51', 'H56', 'H93', 'H102')
        regions.append(r)

    else:
        raise ValueError()

    return regions


def _annotate_sequence(record: SeqRecord, number_mapping: List[AntibodyNumberMapping], numbering_scheme: str,
                       cdr_definition: str = None):
    def region_note(region_, name_):
        return ['antibody_label: {}'.format(name_),
                'antibody_domain: {}'.format(region_.domain),
                'antibody_chain: {}'.format(region_.chain),
                'cdr_definition: {}'.format(cdr_definition),
                'antibody_scheme: {}'.format(numbering_scheme)]

    mapping = [n.to_antibody_number() for n in number_mapping]
    all_regions = _find_all_regions(mapping, numbering_scheme, cdr_definition)
    for region in all_regions:

        start = region.start.query_position
        end = region.end.query_position
        start1 = region.cdr1_start.query_position
        end1 = region.cdr1_end.query_position
        start2 = region.cdr2_start.query_position
        end2 = region.cdr2_end.query_position
        start3 = region.cdr3_start.query_position
        end3 = region.cdr3_end.query_position
        chain = region.chain

        if not start1 or not end1 or not start2 or not end2 or not start3 or not end3:
            return

        if start <= start1:
            name = '{}FR1'.format(chain)
            fr1_feature = SeqFeature(FeatureLocation(start, start1), type='region',
                                     qualifiers={'note': region_note(region, name)})
            record.features.append(fr1_feature)

        if start1 <= end1 + 1:
            name = 'CDR-{}1'.format(chain)
            h1_feature = SeqFeature(FeatureLocation(start1, end1 + 1), type='region',
                                    qualifiers={'note': region_note(region, name)})
            record.features.append(h1_feature)

        if end1 + 1 <= start2:
            name = '{}FR2'.format(chain)
            fr1_feature = SeqFeature(FeatureLocation(end1 + 1, start2), type='region',
                                     qualifiers={'note': region_note(region, name)})
            record.features.append(fr1_feature)

        if start2 <= end2 + 1:
            name = 'CDR-{}2'.format(chain)
            h1_feature = SeqFeature(FeatureLocation(start2, end2 + 1), type='region',
                                    qualifiers={'note': region_note(region, name)})
            record.features.append(h1_feature)

        if end2 + 1 <= start3:
            name = '{}FR3'.format(chain)
            fr1_feature = SeqFeature(FeatureLocation(end2 + 1, start3), type='region',
                                     qualifiers={'note': region_note(region, name)})
            record.features.append(fr1_feature)

        if start3 <= end3 + 1:
            name = 'CDR-{}3'.format(chain)
            h1_feature = SeqFeature(FeatureLocation(start3, end3 + 1), type='region',
                                    qualifiers={'note': region_note(region, name)})
            record.features.append(h1_feature)

        if end3 + 1 <= end + 1:
            name = '{}FR4'.format(chain)
            fr1_feature = SeqFeature(FeatureLocation(end3 + 1, end + 1), type='region',
                                     qualifiers={'note': region_note(region, name)})
            record.features.append(fr1_feature)

    for num in mapping:
        if num.query_position is not None:
            num_feature = SeqFeature(FeatureLocation(num.query_position, num.query_position + 1), type='misc_feature',
                                     qualifiers={'note': ['antibody_number: {}'.format(num.label()),
                                                          'antibody_domain: {}'.format(num.domain),
                                                          'antibody_chain: {}'.format(num.chain),
                                                          'cdr_definition: {}'.format(cdr_definition),
                                                          'antibody_scheme: {}'.format(numbering_scheme)]})
            record.features.append(num_feature)


def numbering_and_regions_from_sequence(seq: SeqRecord) -> Optional[AntibodyAlignmentResult]:
    feature: SeqFeature
    domain_pattern = re.compile(r'^antibody_domain: (\d+)$')
    chain_pattern = re.compile(r'^antibody_chain: (\w+)$')
    cdr_pattern = re.compile(r'^cdr_definition: (\w+)$')
    scheme_pattern = re.compile(r'^antibody_scheme: (\w+)$')
    number_pattern = re.compile(r'^antibody_number: (\w+)$')
    label_pattern = re.compile(r'^([A-Z])(\d+)([A-Z]*)$')

    cdr_definition: Optional[str] = None
    numbering_scheme: Optional[str] = None
    numbering: List[AntibodyNumber] = []
    for feature in seq.features:
        if feature.type != 'misc_feature':
            continue
        if 'note' in feature.qualifiers:
            notes: List[str] = feature.qualifiers['note']
            domain: Optional[int] = None
            chain: Optional[str] = None
            numbering_label: Optional[str] = None

            if all((n.startswith('antibody_') or n.startswith('cdr_') for n in notes)):
                for note in notes:
                    if not cdr_definition:
                        match = cdr_pattern.match(note)
                        if match:
                            cdr_definition = match.group(1)
                            continue
                    if not numbering_scheme:
                        match = scheme_pattern.match(note)
                        if match:
                            numbering_scheme = match.group(1)
                            continue

                    match = domain_pattern.match(note)
                    if match:
                        domain = int(match.group(1))
                        continue
                    match = chain_pattern.match(note)
                    if match:
                        chain = match.group(1)
                        continue
                    match = number_pattern.match(note)
                    if match:
                        numbering_label = match.group(1)
                        continue

            if domain and chain and numbering_label and numbering_scheme:
                match = label_pattern.match(numbering_label)
                if match:
                    assert chain == match.group(1)
                    position = int(match.group(2))
                    insertion = match.group(3)
                    if not insertion:
                        insertion = None
                    number = AntibodyNumber(domain, chain, position, insertion, feature.location.start)
                    numbering.append(number)
    if not numbering:
        return None

    regions = _find_all_regions(numbering, numbering_scheme, cdr_definition)
    return AntibodyAlignmentResult([seq], numbering, regions, numbering_scheme, cdr_definition)


def _find_chain_regions(mapping: List[AntibodyNumber], chain: str, start1: str, end1: str,
                        start2: str, end2: str, start3: str, end3: str) -> Optional[ChainRegions]:
    start = _find_antibody_number(mapping, chain)
    end = _find_antibody_number(mapping, chain, last=True)

    if not start or not end:
        return None

    start1 = _find_antibody_number_from_str(mapping, start1)
    end1 = _find_antibody_number_from_str(mapping, end1)
    start2 = _find_antibody_number_from_str(mapping, start2)
    end2 = _find_antibody_number_from_str(mapping, end2)
    start3 = _find_antibody_number_from_str(mapping, start3)
    end3 = _find_antibody_number_from_str(mapping, end3)

    if not start1 or not end1 or not start2 or not end2 or not start3 or not end3:
        return None

    return ChainRegions(mapping[0].domain, chain, start, end, start1, end1, start2, end2, start3, end3)


def _find_antibody_number_from_str(mapping: List[AntibodyNumber], label: str) -> Optional[AntibodyNumber]:
    if "matcher" not in _find_antibody_number_from_str.__dict__:
        _find_antibody_number_from_str.matcher = re.compile(r'^([HL])(\d+)([A-Z]?)$')
    match = _find_antibody_number_from_str.matcher.match(label)
    assert match
    chain, position, insertion = match.groups()
    if not insertion:
        insertion = None
    return _find_antibody_number(mapping, chain, int(position), insertion)


def _find_antibody_number(mapping: List[AntibodyNumber], chain: str, position: Optional[int] = None,
                          insertion: Optional[str] = None, last: bool = None) -> Optional[AntibodyNumber]:
    last_match = None
    last_insert_match = None
    for num in mapping:
        if num.query_position is None:
            continue
        if num.chain == chain and num.position == position and insertion:
            last_insert_match = num
        if (num.chain == chain and num.position == position and num.insertion == insertion) or (
                num.chain == chain and position and num.position > position) or (
                num.chain == chain and position is None):
            match = num
            if num.chain == chain and position and num.position > position and last_insert_match:
                match = last_insert_match
            if not last:
                return match
            else:
                last_match = match
    return last_match
