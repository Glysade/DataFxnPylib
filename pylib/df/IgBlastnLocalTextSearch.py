import os.path
import uuid

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from df.IgBlastpLocalTextSearch import process_ig_search_results
from df.bio_helper import query_from_request
from df.data_transfer import DataFunction, DataFunctionRequest, DataFunctionResponse, string_input_field
from ruse.bio.applications import IgblastnCommandLine
from ruse.bio.igblast_parse import parse_igblastn_results


class IgBlastnLocalTextSearch(DataFunction):

    def execute(self, request: DataFunctionRequest) -> DataFunctionResponse:
        query_sequence = query_from_request(request)

        id = str(uuid.uuid4())
        v_database = string_input_field(request, 'germlineDbV')
        v_database = os.path.join(os.environ['IGDATA'], 'database', v_database)
        j_database = string_input_field(request, 'germlineDbJ')
        j_database = os.path.join(os.environ['IGDATA'], 'database', j_database)
        d_database = string_input_field(request, 'germlineDbD')
        d_database = os.path.join(os.environ['IGDATA'], 'database', d_database)
        organism = string_input_field(request, 'organism')
        domain_system = string_input_field(request, 'domainSystem')
        format = '"7 std qseq sseq btop"'

        query_file = f'query_{id}.fasta'
        with open(query_file, 'w') as fh:
            ungapped = SeqRecord(query_sequence.seq.ungap(), query_sequence.id)
            SeqIO.write([ungapped], fh, 'fasta')

        out_file = f'igblastn_{id}.out'
        command = IgblastnCommandLine(germline_db_V=v_database, germline_db_D=d_database, germline_db_J=j_database,
                                      organism=organism, domain_system=domain_system,
                                      outfmt=format, query=query_file, out=out_file, show_translation='True')
        stdout, stderr = command()

        # currently, get the 'auxilary data file could not be found' error
        #if stderr or stdout:
        #    raise RuntimeError(f'IGBLASTP failed stderr: {stderr} stdout: {stdout}')

        results = parse_igblastn_results(out_file)
        response = process_ig_search_results(query_sequence, results, 'Igblastn Results')
        return response
