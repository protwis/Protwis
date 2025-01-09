from build.management.commands.base_build import Command as BaseBuild
from django.conf import settings
from protein.models import Protein, ProteinSegment, ProteinConformation, ProteinState, ProteinFamily
from structure.models import Structure, Rotamer
from structure.functions import BlastSearch
from Bio.Blast import NCBIXML, NCBIWWW

import subprocess, shlex, os


class Command(BaseBuild):
    help = 'Blastp search custom dbs'

    def add_arguments(self, parser):
        super(Command, self).add_arguments(parser=parser)
        parser.add_argument('-q', help='Query sequence(s) in FASTA format', default=False, type=str, nargs='+')
        parser.add_argument('-d', help='Query database', default=False, type=str)
        parser.add_argument('--make_db', help='''Create and use custom database. Single argument: (Available presets) 1. xtal - only proteins with structure.
                                                                                                                      2. all - all GPCRs
                                                                                                                      3. fasta text input
                                                                                 Multiple arguments: specific protein entry names''', default=False, type=str, nargs='+')

    def handle(self, *args, **options):
        blastdb = None
        if options['d']:
            blastdb = options['d'] ### FIXME import/parse blast db
        else:
            blastdb = 'blastp_out.fasta'
            if options['make_db']:
                if len(options['make_db'])>1:
                    prots = Protein.objects.filter(entry_name__in=options['make_db'])
                    ### FIXME
                elif len(options['make_db'])==1:
                    prots = []
                    fasta = ''
                    ### xtal preset
                    if options['make_db']==['xtal']:
                        structs = Structure.objects.all().exclude(structure_type__slug__startswith='af-')
                        for i in structs:
                            if i.protein_conformation.protein.parent not in prots:
                                prots.append(i.protein_conformation.protein.parent)
                                fasta+='>{}\n{}\n'.format(i.protein_conformation.protein.parent.entry_name, i.protein_conformation.protein.parent.sequence)
                    elif options['make_db']==['all']:
                        receptor_fams = ProteinFamily.objects.filter(name__startswith='Class')
                        prots = Protein.objects.filter(accession__isnull=False, family__parent__parent__parent__in=receptor_fams)
                        for i in prots:
                            fasta+='>{}\n{}\n'.format(i.entry_name, i.sequence)
                    else:
                        fasta+='>{}\n{}\n'.format('single input', options['make_db'][0])
                    with open('./blastp_out.fasta','w') as f:
                        f.write(fasta)
                make_db_command = shlex.split('makeblastdb -in blastp_out.fasta -dbtype prot -parse_seqids')
                subprocess.call(make_db_command)

        if options['q']:
            for q in options['q']:
                if blastdb:
                    bs = BlastSearch(blastdb=blastdb, top_results=1)
                    out = bs.run(q)
                    for o in out:
                        print(o[0])
                        print(o[1])
                else:
                    bs = BlastSearch()
                    out = bs.run(q)
                    for o in out:
                        for i in o:
                            print(i)

        # if blastdb=='blastp_out.fasta':
        #   files = os.listdir()
        #   for f in files:
        #       if 'blastp_out.fasta' in f:
        #           os.remove(f)

class CustomBlast():
    help = 'Run query on blastdb'

    def __init__(self, fasta_file):
        self.blastdb = fasta_file
        self.blast_search = BlastSearch(blastdb=self.blastdb, top_results=20)

    def run(self, query):
        out = self.blast_search.run(query)
        return out

