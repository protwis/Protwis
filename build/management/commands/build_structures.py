from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection
from django.utils.text import slugify
from django.db import IntegrityError

from build.management.commands.base_build import Command as BaseBuild
from build.management.commands.build_ligand_functions import get_or_create_ligand, match_id_via_unichem
from protein.models import (Protein, ProteinConformation, ProteinState, ProteinAnomaly, ProteinAnomalyType,
    ProteinSegment)
from residue.models import ResidueGenericNumber, ResidueNumberingScheme, Residue, ResidueGenericNumberEquivalent
from common.models import WebLink, WebResource, Publication
from common.tools import test_model_updates
from structure.models import Structure, StructureType, StructureStabilizingAgent,PdbData, Rotamer, Fragment
from construct.functions import *

from contactnetwork.models import *
import contactnetwork.interaction as ci
from contactnetwork.cube import compute_interactions

from Bio.PDB import PDBParser, PPBuilder, Polypeptide
from Bio import pairwise2

from structure.assign_generic_numbers_gpcr import GenericNumbering
from structure.functions import StructureBuildCheck, ParseStructureCSV
from ligand.models import Ligand, LigandType, LigandRole, LigandPeptideStructure
from interaction.models import *
from interaction.views import runcalculation_2022, regexaa, check_residue, extract_fragment_rotamer
from residue.functions import dgn

import django.apps
import logging
import os
import re
import yaml
import time
import gc
from collections import OrderedDict
import json
from urllib.request import urlopen
from Bio.PDB import parse_pdb_header
from Bio.PDB.Selection import *

# import traceback


## FOR VIGNIR ORDERED DICT YAML IMPORT/DUMP
_mapping_tag = yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG
def dict_constructor(loader, node):
    return OrderedDict(loader.construct_pairs(node))

def represent_ordereddict(dumper, data):
    value = []

    for item_key, item_value in data.items():
        node_key = dumper.represent_data(item_key)
        node_value = dumper.represent_data(item_value)

        value.append((node_key, node_value))

    return yaml.nodes.MappingNode(u'tag:yaml.org,2002:map', value)

yaml.add_representer(OrderedDict, represent_ordereddict)
yaml.add_constructor(_mapping_tag, dict_constructor)

class Command(BaseBuild):
    help = 'Reads source data and creates pdb structure records'

    def add_arguments(self, parser):
        parser.add_argument('-p', '--proc',
            type=int,
            action='store',
            dest='proc',
            default=1,
            help='Number of processes to run')
        parser.add_argument('-s', '--structure',
            dest='structure',
            help='Structure to import (PDB ID)',
            nargs='+')
        parser.add_argument('-u', '--purge',
            action='store_true',
            dest='purge',
            default=False,
            help='Purge existing records')
        parser.add_argument('--skip_cn',
            action='store_true',
            default=False,
            help='Skip building contact network for test build')
        parser.add_argument('-i', '--incremental',
            action='store_true',
            dest='incremental',
            default=False,
            help='Incremental update to structures for small live update')
        parser.add_argument('--debug',
            action='store_true',
            dest='debug',
            default=False,
            help='Print info for debugging')

    tracker = {}
    all_models = django.apps.apps.get_models()[6:]
    test_model_updates(all_models, tracker, initialize=True)

    # source file directory
    pdb_data_dir = os.sep.join([settings.DATA_DIR, 'structure_data', 'pdbs'])

    ### USE below to fix seg ends
    xtal_seg_end_file = os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'mod_xtal_segends.yaml'])
    with open(xtal_seg_end_file, 'r') as f:
        xtal_seg_ends = yaml.load(f, Loader=yaml.FullLoader)

    xtal_anomalies_file = os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'all_anomalies.yaml'])
    with open(xtal_anomalies_file, 'r') as f2:
        xtal_anomalies = yaml.load(f2, Loader=yaml.FullLoader)

    xtal_representatives = os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'xtal_representatives.yaml'])
    with open(xtal_representatives, 'r') as f3:
        xtal_representatives = yaml.load(f3, Loader=yaml.FullLoader)

    s = ProteinSegment.objects.all()
    segments = {}
    for segment in s:
        segments[segment.slug] = segment

    parsed_pdb = None

    with open(os.sep.join([settings.DATA_DIR, 'residue_data', 'unnatural_amino_acids.yaml']), 'r') as f_yaml:
        raw_uaa = yaml.safe_load(f_yaml)
        unnatural_amino_acids = {}
        for i, j in raw_uaa.items():
            unnatural_amino_acids[i] = j

    exp_method_dict = {'X-ray': 'X-ray diffraction', 'cryo-EM': 'Electron microscopy', 'Electron crystallography': 'Electron crystallography'}


    def handle(self, *args, **options):
        # delete any existing structure data
        if options['purge']:
            try:
                self.purge_structures()
                self.tracker = {}
                test_model_updates(self.all_models, self.tracker, initialize=True)
            except Exception as msg:
                print(msg)
                self.logger.error(msg)

        if options['skip_cn']:
            self.run_contactnetwork=False
        else:
            self.run_contactnetwork=True

        self.construct_errors, self.rotamer_errors, self.contactnetwork_errors, self.interaction_errors = [],[],[],[]

        self.parsed_structures = ParseStructureCSV()
        self.parsed_structures.parse_ligands()
        self.parsed_structures.parse_nanobodies()
        self.parsed_structures.parse_fusion_proteins()
        self.parsed_structures.parse_ramp()
        self.parsed_structures.parse_grk()

        if options['structure']:
            self.parsed_structures.pdb_ids = [i for i in self.parsed_structures.pdb_ids if i in options['structure'] or i.lower() in options['structure']]

        if options['incremental']:
            self.incremental_mode = True
        else:
            self.incremental_mode = False

        self.debug = options['debug']

        try:
            self.logger.info('CREATING STRUCTURES')
            # run the function twice (once for representative structures, once for non-representative)
            # iterations = 1
            # for i in range(1,iterations+1):
            self.prepare_input(options['proc'], self.parsed_structures.pdb_ids)
            test_model_updates(self.all_models, self.tracker, check=True)
            self.logger.info('COMPLETED CREATING STRUCTURES')
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

        print('Construct errors:')
        print(self.construct_errors)
        print('Rotamer erros:')
        print(self.rotamer_errors)
        print('Contact network errors')
        print(self.contactnetwork_errors)
        print('Interaction errors')
        print(self.interaction_errors)

    @staticmethod
    def queryset_iterator(qs, batchsize = 5000, gc_collect = True):
        # See https://www.guguweb.com/2020/03/27/optimize-django-memory-usage/
        iterator = qs.values_list('pk', flat=True).order_by('pk').distinct().iterator()
        eof = False
        while not eof:
            primary_key_buffer = []
            try:
                while len(primary_key_buffer) < batchsize:
                    primary_key_buffer.append(iterator.__next__())
            except StopIteration:
                eof = True
            for obj in qs.filter(pk__in=primary_key_buffer).order_by('pk').iterator():
                yield obj
            if gc_collect:
                gc.collect()

    def purge_structures(self):
        Structure.objects.all().delete()
        ResidueFragmentInteraction.objects.all().delete()
        ResidueFragmentInteractionType.objects.all().delete()
        StructureLigandInteraction.objects.all().delete()
        #Remove previous Rotamers/Residues to prepare repopulate
        Fragment.objects.all().delete()
        Rotamer.objects.all().delete()

        # FOR DEBUGGING - using queryset_iterator as otherwise the
        # memory was exceeded when deleting >300K PDBs simultaneously
        #for obj in self.queryset_iterator(PdbData.objects.all()):
        #    obj.delete()
        PdbData.objects.all().delete()

    def create_rotamers(self, structure, pdb_path,d):
        wt_lookup = {} #used to match WT seq_number to WT residue record
        pdbseq = {} #used to keep track of pdbseq residue positions vs index in seq
        ref_positions = {} #WT postions in alignment
        mapped_seq = {} # index in contruct, tuple of AA and WT [position,AA]
        debug = False
        preferred_chain = structure.preferred_chain
        wt_pdb_lookup = []

        if len(preferred_chain.split(','))>1: #if A,B
            preferred_chain = preferred_chain.split(',')[0]


        AA = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D',
             'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G',
             'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K',
             'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S',
             'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V',
             'YCM':'C', 'CSD':'C', 'TYS':'Y', 'SEP':'S', 'TPO':'T', 'U8S':'U'} #non-standard AAs

        atom_num_dict = {'E':9, 'S':6, 'Y':12, 'G':4, 'A':5, 'V':7, 'M':8, 'L':8, 'I':8, 'T':7, 'F':11, 'H':10, 'K':9,
                         'D':8, 'C':6, 'R':11, 'P':7, 'Q':9, 'N':8, 'W':14}


        entry_name = d['construct_crystal']['uniprot']

        # print(d['xml_segments'])
        # print(d['deletions'])
        deletions = []
        ## Remove WT ranges that arent in xtal (no need to try to map them)
        # print(d)
        if 'deletions' in d:
            for del_range in d['deletions']:
                if del_range['start']==146 and structure.pdb_code.index=='4K5Y':
                    #Manual fix for faulty 4K5Y annotation
                    continue
                for i in range(del_range['start'],del_range['end']+1):
                    deletions.append(i)
            #print("Annotation missing WT residues",d['deletions'])
        removed = []
        ## Remove segments that arent receptor (tags, fusion etc)
        if 'xml_segments' in d:
            for seg in d['xml_segments']:
                if seg[1]:
                    # Odd rules to fit everything..
                    # print(seg[1][0], entry_name)
                    if seg[1][0]!=entry_name and seg[-1]!=True and seg[1][0]!='Uncharacterized protein' and 'receptor' not in seg[1][0]:
                        if seg[0].split("_")[1]==preferred_chain:
                            #print(seg[2],seg[3]+1)
                            #for i in range(seg[2],seg[3]+1):
                            # print(seg)
                            for i in seg[6]:
                                removed.append(i)
        # Reset removed, since it causes more problems than not

        removed, deletions = construct_structure_annotation_override(structure.pdb_code.index, removed, deletions)

        if self.debug:
            print('Deletions: ', deletions)
            print('Removed: ', removed)
        if len(deletions)>len(d['wt_seq'])*0.9:
            #if too many deletions
            removed = []
            deletions = []

        s = self.parsed_pdb
        chain = s[preferred_chain] #select only one chain (avoid n-mer receptors)

        ppb=PPBuilder()
        seq = ''
        i = 1

        # Checking Na+ atom in xtal
        parent_prot_conf = ProteinConformation.objects.get(protein=structure.protein_conformation.protein.parent)
        try:
            wt_2x50 = Residue.objects.get(protein_conformation=parent_prot_conf, display_generic_number__label=dgn('2x50',parent_prot_conf))
            wt_3x39 = Residue.objects.get(protein_conformation=parent_prot_conf, display_generic_number__label=dgn('3x39',parent_prot_conf))
            if wt_2x50.amino_acid=='D' and wt_3x39.amino_acid=='S':
                if chain[wt_2x50.sequence_number].get_resname()=='ASP' and chain[wt_3x39.sequence_number].get_resname()=='SER':
                    v_2x50 = chain[wt_2x50.sequence_number]['OD1'].get_vector()
                    v_3x39 = chain[wt_3x39.sequence_number]['OG'].get_vector()
                    all_resis = uniqueify(chain)
                    for r in all_resis:
                        id_ = r.get_id()
                        if id_[0]=='H_ NA':
                            v_na = r['NA'].get_vector()
                            d_2x50 = (v_na-v_2x50).norm()
                            d_3x39 = (v_na-v_3x39).norm()
                            if d_2x50<3 and d_3x39<3:
                                structure.sodium = True
                                structure.save()
                                break
        except:
            pass

        check_1000 = 0
        prev_id = 0
        bigjump = False
        all_pdb_residues_in_chain = 0
        for pp in ppb.build_peptides(chain, aa_only=False): #remove >1000 pos (fusion protein / gprotein)
            for i,res in enumerate(pp,1 ):
                all_pdb_residues_in_chain += 1
                residue_id = res.get_full_id()

        if len(removed)+100>all_pdb_residues_in_chain:
            print(structure,'More (or almost) sequence set to be removed from sequence',len(removed),' than exists',all_pdb_residues_in_chain,' removing removed[]')
            #print(removed)
            removed = []

        for pp in ppb.build_peptides(chain, aa_only=False): #remove >1000 pos (fusion protein / gprotein)
            for i,res in enumerate(pp,1 ):
                id = res.id
                residue_id = res.get_full_id()
                if id[1] in removed:
                    chain.detach_child(id)
                    continue
                # if id[1]<600:
                #     check_1000 += 1
                #     #need check_1000 to catch structures where they lie in 1000s (4LDE, 4LDL, 4LDO, 4N4W, 4QKX)
                # if structure.pdb_code.index in ["4RWD","3SN6","4L6R"] and id[1]>1000:
                #     last_valid = 0
                #     bigjump = True
                #     removed.append(id[1])
                # if (id[1]-prev_id)>100 and check_1000>150:
                #     last_valid = prev_id
                #     bigjump = True
                # if bigjump:
                #     if (id[1]-last_valid)<100 or (id[1]<1000 and (id[1]-last_valid)<300 ):
                #         bigjump = False
                # if (id[1]>1000 or bigjump) and check_1000>150 and not (structure.pdb_code.index=="4PHU" and id[1]>2000):
                #     chain.detach_child(id)
                #     #print("removing",id)
                #     removed.append(id[1])
                prev_id = id[1]
        ranges = []
        for k, g in groupby(enumerate(removed), lambda x:x[0]-x[1]):
            group = list(map(itemgetter(1), g))
            ranges.append((group[0], group[-1]))
        if debug: print("Removed XTAL positions due to not being WT receptor",ranges)
        i = 1
        for pp in ppb.build_peptides(chain, aa_only=False):
            seq += str(pp.get_sequence()) #get seq from fasta (only chain A)
            for residue in pp:
                residue_id = residue.get_full_id()
                chain = residue_id[2]
                if chain not in pdbseq:
                    pdbseq[chain] = {}
                pos = residue_id[3][1]

                if residue.resname != "NH2": # skip amidation of peptide
                    pdbseq[chain][pos] = [i, AA[residue.resname]]
                    i += 1

        parent_seq_protein = str(structure.protein_conformation.protein.parent.sequence)
        # print(structure.protein_conformation.protein.parent.entry_name)
        rs = Residue.objects.filter(protein_conformation__protein=structure.protein_conformation.protein.parent).prefetch_related('display_generic_number','generic_number','protein_segment')

        if structure.pdb_code.index.upper() in self.xtal_seg_ends:
            seg_ends = self.xtal_seg_ends[structure.pdb_code.index.upper()]
        else:
            seg_ends = []
            # print('No SEG ENDS info for {}'.format(structure.pdb_code.index))
            self.logger.info('No SEG ENDS info for {}'.format(structure.pdb_code.index))

        parent_seq = ""
        for i,r in enumerate(rs.exclude(sequence_number__in = deletions),1): #required to match WT position to a record (for duplication of GN values)
            # wt_lookup[r.sequence_number] = r
            wt_lookup[i] = r
            parent_seq += r.amino_acid
        # if parent_seq != parent_seq_protein:
        #     print('Residues sequence differ from sequence in protein',structure.protein_conformation.protein.parent.entry_name,structure.pdb_code.index)

        if len(wt_lookup)==0:
            print("No residues for",structure.protein_conformation.protein.parent.entry_name)
            return None

        if self.debug:
            print('parent_seq-pdb_seq length',len(parent_seq),'pdb_seq',len(seq))
            print(parent_seq)
            print(seq)
        #align WT with structure seq -- make gaps penalties big, so to avoid too much overfitting

        if structure.pdb_code.index=='6U1N':
            seq = seq[:265]
        elif structure.pdb_code.index in ['1GZM', '3C9L']:
            seq = seq[:-3]
        elif structure.pdb_code.index=='8WU1':
            seq = seq[:-13]
        if structure.pdb_code.index in ['6NBI','6NBF','6NBH','6U1N','6M1H','6PWC','7JVR','7SHF','7EJ0','7EJ8','7EJA','7EJK','7VVJ','7TS0','7W6P','7W7E','8IRS',
                                        '8FLQ','8FLR','8FLS','8FLU','8FU6','8IRU','7Y35','7Y36','8TB7','8SZI','8TZQ','8U02','8W8Q','8W8R','8W8S','8YN2']:
            pw2 = pairwise2.align.localms(parent_seq, seq, 3, -4, -3, -1)
        elif structure.pdb_code.index in ['6KUX','6KUY','6KUW','7SRS']:
            pw2 = pairwise2.align.localms(parent_seq, seq, 3, -4, -4, -1.5)
        elif structure.pdb_code.index in ['7YMJ']:
            pw2 = pairwise2.align.localms(parent_seq, seq, 3, -5, -4, -4)
        else:
            pw2 = pairwise2.align.localms(parent_seq, seq, 3, -4, -5, -2)

        gaps = 0
        unmapped_ref = {}
        ref_seq, temp_seq = str(pw2[0][0]), str(pw2[0][1])
        # if structure.pdb_code.index in ['5WIU','5WIV']:
        #     temp_seq = temp_seq[:144]+'D'+temp_seq[145:]
        #     temp_seq = temp_seq[:149]+'-'+temp_seq[150:]
        # Custom fixes for alignment issues
        if structure.pdb_code.index=='5ZKP':
            ref_seq = ref_seq[:197]+'-'+ref_seq[198:]
            ref_seq = ref_seq[:198]+'A'+ref_seq[199:]
        elif structure.pdb_code.index in ['5VEW','5VEX']:
            ref_seq = ref_seq[:164]+'IG'+ref_seq[167:]
            temp_seq = temp_seq[:166]+temp_seq[167:]
        elif structure.pdb_code.index in ['3V2W']:
            ref_seq = ref_seq[:201]+ref_seq[202:]
            temp_seq = temp_seq[:207]+temp_seq[208:]
        elif structure.pdb_code.index in ['3V2Y']:
            ref_seq = ref_seq[:209]+ref_seq[210:]
            temp_seq = temp_seq[:215]+temp_seq[216:]
        elif structure.pdb_code.index in ['6KUX','6KUY']:
            ref_seq = ref_seq[:416]+('-'*(416-233))+ref_seq[416:]
            temp_seq = temp_seq[:233]+('-'*(416-233))+temp_seq[233:]
        elif structure.pdb_code.index in ['6KJV','6KK1','6KK7']:
            ref_seq = ref_seq[:176]+ref_seq[177:]
            temp_seq = temp_seq[:178]+temp_seq[179:]
        elif structure.pdb_code.index in ['6LN2']:
            ref_seq = ref_seq[:292]+ref_seq[293:]
            temp_seq = temp_seq[:294]+temp_seq[295:]
        elif structure.pdb_code.index in ['6WHA']:
            temp_seq = temp_seq[:76]+'P--'+temp_seq[79:]
        elif structure.pdb_code.index in ['6X18']:
            temp_seq = temp_seq[:105]+'S------'+temp_seq[112:]
        elif structure.pdb_code.index in ['6K41']:
            temp_seq = temp_seq[:71]+'W------'+temp_seq[78:]
        elif structure.pdb_code.index in ['6LPB']:
            temp_seq = temp_seq[:316]+'S--------'+temp_seq[325:]
        elif structure.pdb_code.index in ['6RZ5']:
            temp_seq = temp_seq[:221]+'K-'+temp_seq[223:]
        elif structure.pdb_code.index in ['6TPK']:
            temp_seq = temp_seq[:30]+'H-----'+temp_seq[36:]
        elif structure.pdb_code.index in ['6W25']:
            temp_seq = temp_seq[:94]+'T--'+temp_seq[97:]
        elif structure.pdb_code.index in ['7C61']:
            temp_seq = temp_seq[:207]+'-'+temp_seq[207:225]+'-'+temp_seq[225:238]+'-'+temp_seq[238:251]+'-'+temp_seq[251:]
            ref_seq = ref_seq[:304]+'----'+ref_seq[304:]
        elif structure.pdb_code.index in ['6WHC']:
            temp_seq = temp_seq[:202]+'S----------'+temp_seq[213:]
        elif structure.pdb_code.index=='5T1A':
            temp_seq = temp_seq[:229]+temp_seq[231:238]+temp_seq[239:242]+temp_seq[245:]
            ref_seq = ref_seq[:224]+ref_seq[227:233]+ref_seq[236:]
        elif structure.pdb_code.index=='5UEN':
            temp_seq = temp_seq[:218]+temp_seq[222:]
            ref_seq = ref_seq[:210]+ref_seq[214:]
        elif structure.pdb_code.index=='6DO1':
            temp_seq = temp_seq[:228]+temp_seq[230:]
            ref_seq = ref_seq[:225]+ref_seq[227:]
        elif structure.pdb_code.index=='7KH0':
            temp_seq = temp_seq[:240]+'S'+temp_seq[240:262]+temp_seq[263:]
        elif structure.pdb_code.index=='7BB6':
            temp_seq = temp_seq[:231]+'A'+temp_seq[231:257]+temp_seq[258:]
        elif structure.pdb_code.index=='7C4S':
            temp_seq = temp_seq[:224]+'S'+temp_seq[224:234]+temp_seq[235:]
        elif structure.pdb_code.index=='7M3J':
            temp_seq = temp_seq[:100]+'I'+temp_seq[100:115]+temp_seq[116:337]+'N'+temp_seq[337:380]+temp_seq[381:]
        elif structure.pdb_code.index=='7DUQ':
            temp_seq = temp_seq[:105]+'S------'+temp_seq[112:]
        elif structure.pdb_code.index in ['7KI0','7KI1']:
            temp_seq = temp_seq[:105]+'S-------'+temp_seq[113:]
        elif structure.pdb_code.index=='7MTQ':
            temp_seq = temp_seq[:694]+'E---'+temp_seq[698:]
        elif structure.pdb_code.index=='7FD9':
            temp_seq = temp_seq[:99]+'S'+temp_seq[99:118]+temp_seq[119:]
        elif structure.pdb_code.index in ['6ZFZ', '6ZG4', '6ZG9']:
            temp_seq = temp_seq[:28]+temp_seq[29:207]+temp_seq[208:]
            ref_seq = ref_seq[:23]+ref_seq[24:211]+ref_seq[212:]
        elif structure.pdb_code.index in ['7NA7', '7NA8']:
            temp_seq = temp_seq[:242]+'R'+temp_seq[242:253]+temp_seq[254:]
        elif structure.pdb_code.index in ['7F8V', '7F8W']:
            temp_seq = temp_seq[:247]+'L'+temp_seq[247:323]+temp_seq[324:]
        elif structure.pdb_code.index=='7RTB':
            temp_seq = temp_seq[:107]+'S'+temp_seq[107:113]+temp_seq[114:]
        elif structure.pdb_code.index=='6Z4Q':
            temp_seq = temp_seq[:131]+'H-'+temp_seq[133:]
        elif structure.pdb_code.index=='6ZA8':
            temp_seq = temp_seq[:212]+'L'+temp_seq[212:219]+temp_seq[220:]
        elif structure.pdb_code.index=='7EO4':
            temp_seq = temp_seq[:36]+'I'+temp_seq[36:44]+temp_seq[45:]
        elif structure.pdb_code.index in ['7EWP','7EWR']:
            temp_seq = temp_seq[:667]+'S'+temp_seq[667:701]+temp_seq[702:]
        elif structure.pdb_code.index in ['7T10','7T11']:
            temp_seq = temp_seq[:237]+temp_seq[239:246]+'R'+temp_seq[251:]
            ref_seq = ref_seq[:244]+ref_seq[245:253]+ref_seq[258:]
        elif structure.pdb_code.index in ['7FIY']:
            temp_seq = temp_seq[:306]+'R---'+temp_seq[310:]
        elif structure.pdb_code.index in ['7B6W']:
            temp_seq = temp_seq[:318]+'L----'+temp_seq[323:]
            temp_seq = temp_seq[:245]+temp_seq[246:282]+'-S'+temp_seq[283:]
        elif structure.pdb_code.index in ['7SIL','7SIM','7SIN']:
            temp_seq = temp_seq[:702]+'L'+temp_seq[702:720]+temp_seq[721:]
        elif structure.pdb_code.index in ['7WIH']:
            temp_seq = temp_seq[:26]+'E'+temp_seq[26:33]+temp_seq[34:94]+'L'+temp_seq[94:117]+temp_seq[118:]
        elif structure.pdb_code.index=='7SBF':
            temp_seq = temp_seq[:8]+temp_seq[10:]
            ref_seq = ref_seq[2:]
        elif structure.pdb_code.index=='7RA3':
            temp_seq = temp_seq[:306]+'R'+temp_seq[306:311]+temp_seq[312:]
        elif structure.pdb_code.index in ['7EJ8']:
            temp_seq = temp_seq[:182]+'P'+temp_seq[182:200]+temp_seq[201:]
        elif structure.pdb_code.index in ['7EJ0']:
            temp_seq = temp_seq[:242]+'R'+temp_seq[242:377]+temp_seq[378:]
        elif structure.pdb_code.index in ['7EJK']:
            temp_seq = temp_seq[:182]+'P'+temp_seq[182:200]+temp_seq[201:242]+'R'+temp_seq[242:379]+temp_seq[380:]
        elif structure.pdb_code.index=='7EZC':
            temp_seq = temp_seq[:146]+'Q'+temp_seq[146:155]+temp_seq[156:]
        elif structure.pdb_code.index=='7RBT':
            temp_seq = temp_seq[:306]+'R'+temp_seq[306:311]+temp_seq[312:]
        elif structure.pdb_code.index=='7WUJ':
            temp_seq = temp_seq[:158]+'T'+temp_seq[158:163]+temp_seq[164:]
        elif structure.pdb_code.index=='5JQH':
            temp_seq = temp_seq[:211]+'--D'+temp_seq[214:]
        elif structure.pdb_code.index=='2YCW':
            temp_seq = temp_seq[:242]+'R'+temp_seq[242:270]+temp_seq[271:]
        elif structure.pdb_code.index=='7EPT':
            temp_seq = temp_seq[:197]+'SA'+temp_seq[197:208]+temp_seq[210:]
        elif structure.pdb_code.index=='7SK5':
            temp_seq = temp_seq[:186]+'S--'+temp_seq[189:]
        elif structure.pdb_code.index=='7WU9':
            temp_seq = temp_seq[:261]+'Q'+temp_seq[261:275]+temp_seq[276:]
        elif structure.pdb_code.index=='7SRS':
            temp_seq = temp_seq[:246]+'VRLLS'+61*'-'+'R'+temp_seq[313:]
        elif structure.pdb_code.index=='7UL2':
            ref_seq = ref_seq[:257]+ref_seq[259:305]+ref_seq[306:]
            temp_seq = temp_seq[:264]+'SVRL'+19*'-'+'LSGS'+temp_seq[291:296]+temp_seq[298:308]+temp_seq[309:]
        elif structure.pdb_code.index=='7UL3':
            ref_seq = ref_seq[:202]+ref_seq[206:211]+ref_seq[212:236]+ref_seq[242:244]+ref_seq[246:255]+ref_seq[256:]
            temp_seq = temp_seq[:208]+'LKSVRLLS'+5*'-'+'SRE'+temp_seq[235:247]+'L'+temp_seq[251:]
        elif structure.pdb_code.index=='7UL5':
            ref_seq = ref_seq[:244]+ref_seq[245:253]+ref_seq[258:]
            temp_seq = temp_seq[:237]+temp_seq[239:242]+'LSGSR'+temp_seq[251:]
        elif structure.pdb_code.index=='7PP1':
            temp_seq = temp_seq[:128]+'P'+temp_seq[128:134]+temp_seq[135:161]+'L'+temp_seq[161:177]+temp_seq[178:]
        elif structure.pdb_code.index=='7RAN':
            temp_seq = temp_seq[:283]+'C'+temp_seq[283:287]+temp_seq[288:]
        elif structure.pdb_code.index=='7S0F':
            temp_seq = temp_seq[:247]+temp_seq[248:283]+'F'+temp_seq[283:]
        elif structure.pdb_code.index=='7VIH':
            temp_seq = temp_seq[:33]+'KL'+temp_seq[33:45]+temp_seq[47:]
        elif structure.pdb_code.index=='7VQX':
            temp_seq = temp_seq[:288]+'S'+temp_seq[288:297]+temp_seq[298:]
        elif structure.pdb_code.index in ['7VVK']:
            temp_seq = temp_seq[:4]+temp_seq[65:84]+temp_seq[4:65]+temp_seq[84:]
        elif structure.pdb_code.index in ['7VVL']:
            temp_seq = temp_seq[:4]+temp_seq[62:81]+temp_seq[4:62]+temp_seq[81:]
        elif structure.pdb_code.index in ['7VVM']:
            temp_seq = temp_seq[:4]+temp_seq[61:80]+temp_seq[4:61]+temp_seq[80:]
        elif structure.pdb_code.index=='7VVN':
            temp_seq = temp_seq[:6]+temp_seq[78:95]+temp_seq[6:67]+temp_seq[95:102]+11*'-'+temp_seq[102:365]+'T'+temp_seq[365:372]+temp_seq[373:403]+'T'+temp_seq[403:408]+temp_seq[409:]
        elif structure.pdb_code.index=='7VVO':
            temp_seq = temp_seq[:6]+temp_seq[79:99]+temp_seq[6:79]+temp_seq[99:]
        elif structure.pdb_code.index=='7W57':
            temp_seq = temp_seq[:192]+'P'+temp_seq[192:198]+temp_seq[199:]
        elif structure.pdb_code.index in ['7W6P']:
            temp_seq = temp_seq[:182]+'P'+temp_seq[182:200]+temp_seq[201:]
        elif structure.pdb_code.index=='7W7E':
            temp_seq = temp_seq[:182]+'P'+temp_seq[182:200]+temp_seq[201:242]+'R'+temp_seq[242:377]+temp_seq[378:]
        elif structure.pdb_code.index=='7WBJ':
            temp_seq = temp_seq[:288]+'S'+temp_seq[288:297]+temp_seq[298:]
        elif structure.pdb_code.index in ['7X8R']:
            temp_seq = temp_seq[:7]+temp_seq[89:113]+temp_seq[7:89]+temp_seq[113:]
        elif structure.pdb_code.index in ['7X8S']:
            temp_seq = temp_seq[:7]+temp_seq[90:114]+temp_seq[7:90]+temp_seq[114:]
        elif structure.pdb_code.index=='8HA0':
            temp_seq = temp_seq[:4]+temp_seq[58:78]+temp_seq[4:58]+temp_seq[78:]
        elif structure.pdb_code.index=='8HAF':
            temp_seq = temp_seq[:4]+temp_seq[53:78]+temp_seq[4:53]+temp_seq[78:]
        elif structure.pdb_code.index=='8HAO':
            temp_seq = temp_seq[:4]+temp_seq[57:78]+temp_seq[4:57]+temp_seq[78:]
        elif structure.pdb_code.index=='7T8X':
            temp_seq = temp_seq[:214]+'K'+temp_seq[214:237]+temp_seq[238:]
        elif structure.pdb_code.index in ['7ZBE','8A6C']:
            temp_seq = temp_seq[:228]+'T'+temp_seq[228:242]+temp_seq[243:]
        elif structure.pdb_code.index=='8FMZ':
            temp_seq = temp_seq[:172]+'A-'+temp_seq[174:]
        elif structure.pdb_code.index=='8ID4':
            temp_seq = temp_seq[:72]+'A--'+temp_seq[75:]
        elif structure.pdb_code.index=='7XJJ':
            temp_seq = temp_seq[:140]+'R'+temp_seq[140:146]+temp_seq[147:]
        elif structure.pdb_code.index=='8DZS':
            temp_seq = temp_seq[:247]+'S----'+temp_seq[252:]
        elif structure.pdb_code.index=='8G94':
            temp_seq = temp_seq[:36]+'I'+temp_seq[36:44]+temp_seq[45:]
        elif structure.pdb_code.index=='8IW1':
            temp_seq = temp_seq[:170]+'G'+temp_seq[170:188]+temp_seq[189:]
        elif structure.pdb_code.index in ['8IW4','8IWE']:
            temp_seq = temp_seq[:180]+'V-'+temp_seq[182:]
        elif structure.pdb_code.index in ['8JWY','8JWZ']:
            ref_seq = ref_seq[:221]+ref_seq[222:]
            temp_seq = temp_seq[:217]+temp_seq[218:]
        elif structure.pdb_code.index in ['7Y35','7Y36']:
            temp_seq = temp_seq[:31]+'R'+temp_seq[31:77]+temp_seq[78:]
        elif structure.pdb_code.index=='7YMJ':
            ref_seq = ref_seq[:215]+ref_seq[216:]
            temp_seq = temp_seq[:212]+temp_seq[214:217]+40*'-'+temp_seq[218:222]+temp_seq[223:227]+'D'+temp_seq[232:236]+'RITRLVL'+temp_seq[276:]
        elif structure.pdb_code.index=='8H0P':
            temp_seq = temp_seq[:244]+'N'+temp_seq[244:250]+temp_seq[251:]
        elif structure.pdb_code.index in ['8JCV','8JCX']:
            temp_seq = temp_seq[:641]+'I'+temp_seq[641:655]+temp_seq[656:]
        elif structure.pdb_code.index in ['8JD1']:
            temp_seq = temp_seq[:642]+'F'+temp_seq[642:654]+temp_seq[655:]
        elif structure.pdb_code.index in ['8JRV']:
            temp_seq = temp_seq[79:100]+temp_seq[:79]+temp_seq[100:]
        elif structure.pdb_code.index in ['8GTG','8GTM']:
            temp_seq = temp_seq[:117]+'T---'+temp_seq[121:]
        elif structure.pdb_code.index in ['8HTI','8J46','8W77']:
            ref_seq = ref_seq[:114]+ref_seq[115:145]+'IL'+ref_seq[149:160]+ref_seq[161:170]+ref_seq[171:215]+ref_seq[216:]
            temp_seq = temp_seq[:118]+temp_seq[119:150]+temp_seq[152:156]+temp_seq[157:168]+temp_seq[169:209]+temp_seq[210:]
        elif structure.pdb_code.index=='8J24':
            ref_seq = ref_seq[:245]+ref_seq[246:]
            temp_seq = temp_seq[:247]+temp_seq[248:]
        elif structure.pdb_code.index in ['7YFC','7FYD']:
            temp_seq = temp_seq[:203]+'R'+temp_seq[203:296]+temp_seq[297:]
        elif structure.pdb_code.index in ['8GGA','8GGP']:
            temp_seq = temp_seq[:236]+'E'+temp_seq[236:267]+temp_seq[268:]
        elif structure.pdb_code.index=='8GGB':
            temp_seq = temp_seq[:95]+'T-'+temp_seq[97:236]+'E'+temp_seq[236:267]+temp_seq[268:]
        elif structure.pdb_code.index=='8GGE':
            temp_seq = temp_seq[:186]+'E-'+temp_seq[188:]
        elif structure.pdb_code.index=='8GTI':
            temp_seq = temp_seq[:219]+'T---'+temp_seq[223:]
        elif structure.pdb_code.index in ['8HN8','8HOC']:
            ref_seq = ref_seq[:180]+ref_seq[181:]
            temp_seq = temp_seq[:182]+temp_seq[183:202]+'G'+temp_seq[202:222]+'-----'+temp_seq[230:291]+'HL'+temp_seq[291:]
        elif structure.pdb_code.index=='8IKH':
            temp_seq = temp_seq[:78]+'F---'+temp_seq[82:]
        elif structure.pdb_code.index in ['8J22','8J23']:
            ref_seq = ref_seq[:245]+ref_seq[246:]
            temp_seq = temp_seq[:247]+temp_seq[248:]
        elif structure.pdb_code.index=='8JXW':
            temp_seq = temp_seq[:205]+'H'+temp_seq[205:291]+temp_seq[292:]
        elif structure.pdb_code.index=='8KH5':
            temp_seq = temp_seq[:154]+'RT'+temp_seq[154:163]+temp_seq[165:]
        elif structure.pdb_code.index=='8PJK':
            temp_seq = temp_seq[:175]+'R'+temp_seq[175:180]+temp_seq[181:]
        elif structure.pdb_code.index=='8QJ2':
            temp_seq = temp_seq[:70]+'KN--'+temp_seq[74:]
        elif structure.pdb_code.index=='8QW4':
            temp_seq = temp_seq[:460]+'Y'+temp_seq[460:466]+temp_seq[467:]
        elif structure.pdb_code.index=='8SZF':
            temp_seq = temp_seq[:684]+'L'+temp_seq[684:702]+temp_seq[703:]
        elif structure.pdb_code.index=='8T3Q':
            temp_seq = temp_seq[:123]+'RG---'+temp_seq[128:]
        elif structure.pdb_code.index=='8TR2':
            temp_seq = temp_seq[:472]+'QT-G'+temp_seq[476:511]+'P---'+temp_seq[515:629]+'F--'+temp_seq[632:]
        elif structure.pdb_code.index in ['8TRC','8TRD']:
            temp_seq = temp_seq[:472]+'Q----'+temp_seq[477:]
        elif structure.pdb_code.index=='8UWL':
            temp_seq = temp_seq[:200]+'T'+temp_seq[200:245]+temp_seq[246:]
        elif structure.pdb_code.index=='8V6U':
            temp_seq = temp_seq[:283]+'C----'+temp_seq[288:]
        elif structure.pdb_code.index=='8WCB':
            temp_seq = temp_seq[:178]+'G-'+temp_seq[180:]
        elif structure.pdb_code.index=='8WKY':
            temp_seq = temp_seq[:94]+'TD--'+temp_seq[98:]
        elif structure.pdb_code.index=='8WPG':
            temp_seq = temp_seq[:683]+'L'+temp_seq[683:701]+temp_seq[702:]
        elif structure.pdb_code.index in ['8X79','8X7A']:
            temp_seq = temp_seq[:221]+'G'+temp_seq[221:230]+temp_seq[231:]
        elif structure.pdb_code.index=='8XQO':
            temp_seq = temp_seq[:160]+'N'+temp_seq[160:169]+temp_seq[170:215]+'K'+temp_seq[215:223]+temp_seq[224:]
        elif structure.pdb_code.index=='8XQP':
            temp_seq = temp_seq[:160]+'N'+temp_seq[160:169]+temp_seq[170:]
        elif structure.pdb_code.index in ['8YW3']:
            temp_seq = temp_seq[:105]+'S'+temp_seq[105:111]+temp_seq[112:]
        elif structure.pdb_code.index=='8YW4':
            temp_seq = temp_seq[:9]+temp_seq[87:105]+temp_seq[9:87]+temp_seq[105:]
        elif structure.pdb_code.index=='8YW5':
            temp_seq = temp_seq[:20]+'L'+temp_seq[20:29]+temp_seq[30:]
        elif structure.pdb_code.index=='8ZFJ':
            temp_seq = temp_seq[:262]+'C--'+temp_seq[265:]
        elif structure.pdb_code.index=='8ZSJ':
            temp_seq = temp_seq[:228]+'K'+temp_seq[228:243]+temp_seq[244:]
        elif structure.pdb_code.index=='9AVL':
            temp_seq = temp_seq[:123]+'N'+temp_seq[123:129]+temp_seq[130:]
        elif structure.pdb_code.index=='8UXY':
            ref_seq = ref_seq[:90]+ref_seq[92:108]+ref_seq[109:151]+ref_seq[152:200]+ref_seq[201:213]+ref_seq[214:]
            temp_seq = temp_seq[:84]+'I'+temp_seq[87:105]+temp_seq[106:149]+temp_seq[150:202]+temp_seq[203:211]+temp_seq[212:271]+'-YS'+temp_seq[274:]
        elif structure.pdb_code.index=='8UXV':
            ref_seq = ref_seq[:144]+ref_seq[145:210]+ref_seq[211:266]+ref_seq[267:]
            temp_seq = temp_seq[:146]+temp_seq[147:208]+temp_seq[209:268]+temp_seq[269:]
        elif structure.pdb_code.index=='8WVV':
            temp_seq = temp_seq[:542]+'S---'+temp_seq[546:673]+'P-----'+temp_seq[679:]
        elif structure.pdb_code.index=='8XWP':
            temp_seq = temp_seq[:235]+'L'+temp_seq[235:245]+temp_seq[246:]
        elif structure.pdb_code.index=='8YN4':
            temp_seq = temp_seq[:215]+'I'+temp_seq[215:226]+temp_seq[227:]
        elif structure.pdb_code.index in ['9JR2']:
            temp_seq = temp_seq[:5]+temp_seq[58:81]+temp_seq[5:58]+temp_seq[81:]
        elif structure.pdb_code.index in ['9JR3']:
            temp_seq = temp_seq[:4]+temp_seq[57:81]+temp_seq[4:57]+temp_seq[81:]


        for i, r in enumerate(ref_seq, 1): #loop over alignment to create lookups (track pos)
            if self.debug:
                print(i,r,temp_seq[i-1]) #print alignment for sanity check
            if r == "-":
                gaps += 1
            if r != "-":
                ref_positions[i] = [i-gaps,r]
            elif r == "-":
                ref_positions[i] = [None,'-']

            if temp_seq[i-1]=='-':
                unmapped_ref[i-gaps] = '-'

        gaps = 0
        for i, r in enumerate(temp_seq, 1): #make second lookup
            # print(i,r,ref_seq[i-1]) #print alignment for sanity check
            if r == "-":
                gaps += 1
            if r != "-":
                mapped_seq[i-gaps] = [r,ref_positions[i]]
                # if r!=ref_seq[i-1]:
                #     print('aa mismatch')
        # print("seg res not mapped",gaps)

        pdb = structure.pdb_data.pdb
        protein_conformation=structure.protein_conformation
        temp = ''
        check = 0
        errors = 0
        mismatch_seq = 0
        match_seq = 0
        not_matched = 0
        matched_by_pos = 0
        aa_mismatch = 0
        generic_change = 0
        generic_change_from = []

        pdblines_temp = pdb.splitlines()
        pdblines = []
        for line in pdblines_temp: #Get rid of all odd records
            if line.startswith('ATOM') or (line[17:20] in ['YCM','CSD','TYS','SEP'] and line.startswith('HETATM')):
                pdblines.append(line)
            # Only build residues for the first model
            if line.startswith('MODEL        2'):
                break
        pdblines.append('') #add a line to not "run out"
        rotamer_bulk = []
        rotamer_data_bulk = []
        residues_bulk = []
        if structure.pdb_code.index=='5LWE':
            seg_ends['5b'] = 209
            seg_ends['5e'] = 244
        # import pprint
        # pprint.pprint(wt_lookup)
        # pprint.pprint(mapped_seq)
        # pprint.pprint(unmapped_ref)
        # print('deletions: ',deletions)
        # print('removed: ',removed)
        for i,line in enumerate(pdblines):
            # print(line)
            if line.startswith('ATOM') or (line[17:20] in ['YCM','CSD','TYS','SEP'] and line.startswith('HETATM')):
                # if line[17:20] in ['YCM','CSD','TYS','SEP']: # sanity check for non-standard helix residues
                #     print(line)
                chain = line[21]
                if preferred_chain and chain!=preferred_chain: #If perferred is defined and is not the same as the current line, then skip
                    pass
                else:
                    nextline = pdblines[i+1]
                    residue_number = line[22:26].strip()
                    if (check==0 or nextline[22:26].strip()==check) and nextline.startswith('TER')==False and (nextline.startswith('ATOM')==True or nextline.startswith('HETATM')==True): #If this is either the begining or the same as previous line add to current rotamer
                        temp += line + "\n"
                        #print('same res',pdb.splitlines()[i+1])
                    else: #if this is a new residue
                        #print(pdb.splitlines()[i+1][22:26].strip(),check)
                        temp += line + "\n"
                        if structure.pdb_code.index in ['7E9H','8JD6'] and line[17:20]=='SEP':
                            continue
                        #(int(check.strip())<2000 or structure.pdb_code.index=="4PHU") and
                        if int(check.strip()) not in removed:
                            # print(line)
                            residue = Residue()
                            residue.sequence_number = int(check.strip())
                            residue.amino_acid = AA[residue_name.upper()]
                            residue.protein_conformation = protein_conformation
                            try:
                                seq_num_pos = pdbseq[chain][residue.sequence_number][0]
                            except:
                                # print('failed residue',pdb_path,residue.sequence_number)
                                temp = "" #start new line for rotamer
                                check = pdblines[i+1][22:26].strip()
                                continue
                            # print('hi',seq_num_pos,residue.sequence_number)
                            if seq_num_pos in mapped_seq:
                                # print(int(check.strip()),seq_num_pos) #sanity check
                                if mapped_seq[seq_num_pos][1][0]==None:
                                    # print('no match found') #sanity check
                                    # print(residue.sequence_number,residue.amino_acid) #sanity check
                                    residue.display_generic_number = None
                                    residue.generic_number = None
                                    residue.protein_segment = None
                                    not_matched +=1
                                else:
                                    wt_r = wt_lookup[mapped_seq[seq_num_pos][1][0]]
                                    # print(seq_num_pos, mapped_seq[seq_num_pos], wt_r, residue.sequence_number, residue.amino_acid) #sanity check for mapped resis
                                    if residue.sequence_number!=wt_r.sequence_number and residue.amino_acid!=wt_r.amino_acid and residue.sequence_number in wt_lookup: #if pos numbers not work -- see if the pos number might be in WT and unmapped
                                        if wt_lookup[residue.sequence_number].amino_acid==residue.amino_acid:
                                            if residue.sequence_number in unmapped_ref: #WT was not mapped, so could be it
                                                # print(residue.sequence_number,residue.amino_acid) #sanity check
                                                # print('wrongly matched, better match on pos+aa',residue.sequence_number,residue.amino_acid,wt_r.sequence_number,wt_r.amino_acid)
                                                wt_r = wt_lookup[residue.sequence_number]
                                                matched_by_pos +=1
                                                match_seq += 1
                                            else:
                                                mismatch_seq += 1
                                                ### REPLACE seq number with WT to fix odd PDB annotation. FIXME kinda dangerous, but best way to ensure consistent GN numbering
                                                residue.sequence_number = wt_r.sequence_number
                                                #print('could have been matched, but already aligned to another position',residue.sequence_number,residue.amino_acid,wt_r.sequence_number,wt_r.amino_acid)
                                        else:
                                            # print('WT pos not same AA, mismatch',residue.sequence_number,residue.amino_acid,wt_r.sequence_number,wt_r.amino_acid)
                                            wt_pdb_lookup.append(OrderedDict([('WT_POS',wt_r.sequence_number), ('PDB_POS',residue.sequence_number), ('AA','.')]))
                                            mismatch_seq += 1
                                            aa_mismatch += 1
                                    elif residue.sequence_number!=wt_r.sequence_number:
                                        # print('WT pos not same pos, mismatch',residue.sequence_number,residue.amino_acid,wt_r.sequence_number,wt_r.amino_acid)
                                        wt_pdb_lookup.append(OrderedDict([('WT_POS',wt_r.sequence_number), ('PDB_POS',residue.sequence_number), ('AA',wt_r.amino_acid)]))
                                        if structure.pdb_code.index not in ['4GBR','6C1R','6C1Q','7XBX','7F1Q','7ZLY','8JWY','8JWZ','8JMT','8TB7','8ITM','9D3G','9D3E']:
                                            if residue.sequence_number in unmapped_ref:
                                                # print('residue.sequence_number',residue.sequence_number,'not mapped though')
                                                if residue.amino_acid == wt_lookup[residue.sequence_number].amino_acid:
                                                    #print('they are same amino acid!')
                                                    wt_r = wt_lookup[residue.sequence_number]
                                                    mismatch_seq -= 1
                                        mismatch_seq += 1
                                        ### REPLACE seq number with WT to fix odd PDB annotation. FIXME kinda dangerous, but best way to ensure consistent GN numbering
                                        ### 2019.01.18 DISABLED underneat, to be sure that sequence number can be found in DB correctly.
                                        # residue.sequence_number = wt_r.sequence_number
                                    if residue.amino_acid!=wt_r.amino_acid:
                                        if debug: print('aa mismatch',residue.sequence_number,residue.amino_acid,wt_r.sequence_number,wt_r.amino_acid)
                                        aa_mismatch += 1

                                    else:
                                        match_seq += 1
                                    if wt_r.generic_number is not None:
                                        residue.display_generic_number = wt_r.display_generic_number
                                        residue.generic_number = wt_r.generic_number
                                    else:
                                        residue.display_generic_number = None
                                        residue.generic_number = None
                                        # print('no GN')

                                    residue.protein_segment = wt_r.protein_segment
                                    residue.missing_gn = False
                                    # print(residue, residue.protein_segment, residue.display_generic_number, wt_r, wt_r.protein_segment, wt_r.display_generic_number)
                                    if len(seg_ends):
                                        if residue.protein_segment.slug=='TM1':
                                            if seg_ends['1b']!='-' and seg_ends['1e']!='-':
                                                if residue.sequence_number<seg_ends['1b']:
                                                    residue.protein_segment = self.segments['N-term']
                                                    residue.display_generic_number = None
                                                    residue.generic_number = None
                                                elif residue.sequence_number>seg_ends['1e']:
                                                    residue.protein_segment = self.segments['ICL1']
                                                    residue.display_generic_number = None
                                                    residue.generic_number = None
                                        elif residue.protein_segment.slug=='ICL1':
                                            if seg_ends['i1b']!='-' and seg_ends['i1e']!='-':
                                                if residue.sequence_number<seg_ends['i1b'] and residue.sequence_number<=seg_ends['1e']:
                                                    residue.protein_segment = self.segments['TM1']
                                                # elif residue.sequence_number>seg_ends['i1e']:
                                                #     residue.protein_segment = self.segments['TM2']
                                                elif (residue.sequence_number>=seg_ends['i1b'] and residue.sequence_number<=seg_ends['i1e']) and residue.generic_number is None:
                                                    if debug: print("Missing GN in loop!",residue.sequence_number)
                                                    residue.missing_gn = True
                                                elif (residue.sequence_number<seg_ends['i1b'] and residue.sequence_number>seg_ends['i1e']):
                                                    if debug: print("Remove former GN for ",residue.generic_number)
                                                    residue.display_generic_number = None
                                                    residue.generic_number = None
                                            else:
                                                if residue.generic_number is not None:
                                                    if debug: print("Remove former GN for ",residue.generic_number)
                                                    residue.display_generic_number = None
                                                    residue.generic_number = None
                                                if residue.sequence_number<=seg_ends['1e']:
                                                    residue.protein_segment = self.segments['TM1']
                                                elif residue.sequence_number>=seg_ends['2b']:
                                                    residue.protein_segment = self.segments['TM2']
                                        elif residue.protein_segment.slug=='TM2':
                                            if seg_ends['2b']!='-' and seg_ends['2e']!='-':
                                                if residue.sequence_number<seg_ends['2b']:
                                                    residue.protein_segment = self.segments['ICL1']
                                                    residue.display_generic_number = None
                                                    residue.generic_number = None
                                                elif residue.sequence_number>seg_ends['2e']:
                                                    residue.protein_segment = self.segments['ECL1']
                                                    residue.display_generic_number = None
                                                    residue.generic_number = None
                                        elif residue.protein_segment.slug=='ECL1':
                                            if seg_ends['e1b']!='-' and seg_ends['e1e']!='-' and seg_ends['2e']!='-':
                                                if residue.sequence_number<seg_ends['e1b'] and residue.sequence_number<=seg_ends['2e']:
                                                    residue.protein_segment = self.segments['TM2']
                                                # elif residue.sequence_number>seg_ends['e1e']:
                                                #     residue.protein_segment = self.segments['TM3']
                                                elif (residue.sequence_number>=seg_ends['e1b'] and residue.sequence_number<=seg_ends['e1e']) and residue.generic_number is None:
                                                    if debug: print("Missing GN in loop!",residue.sequence_number)
                                                    residue.missing_gn = True
                                                elif (residue.sequence_number<seg_ends['e1b'] and residue.sequence_number>seg_ends['e1e']):
                                                    if debug: print("Remove former GN for ",residue.generic_number)
                                                    residue.display_generic_number = None
                                                    residue.generic_number = None
                                            else:
                                                if residue.generic_number is not None:
                                                    if debug: print("Remove former GN for ",residue.generic_number)
                                                    residue.display_generic_number = None
                                                    residue.generic_number = None
                                                if residue.sequence_number<=seg_ends['2e']:
                                                    residue.protein_segment = self.segments['TM2']
                                                elif residue.sequence_number>=seg_ends['3b']:
                                                    residue.protein_segment = self.segments['TM3']
                                        elif residue.protein_segment.slug=='TM3':
                                            if seg_ends['3b']!='-' and seg_ends['3e']!='-':
                                                if residue.sequence_number<seg_ends['3b']:
                                                    residue.protein_segment = self.segments['ECL1']
                                                    residue.display_generic_number = None
                                                    residue.generic_number = None
                                                elif residue.sequence_number>seg_ends['3e']:
                                                    residue.protein_segment = self.segments['ICL2']
                                                    residue.display_generic_number = None
                                                    residue.generic_number = None
                                        elif residue.protein_segment.slug=='ICL2':
                                            if seg_ends['i2b']!='-' and seg_ends['i2e']!='-':
                                                # if residue.sequence_number<seg_ends['i2b']:
                                                #     residue.protein_segment = self.segments['TM3']
                                                if residue.sequence_number>seg_ends['i2e'] and residue.sequence_number>=seg_ends['4b']:
                                                    residue.protein_segment = self.segments['TM4']
                                                elif (residue.sequence_number>=seg_ends['i2b'] and residue.sequence_number<=seg_ends['i2e']) and residue.generic_number is None:
                                                    if debug: print("Missing GN in loop!",residue.sequence_number)
                                                    residue.missing_gn = True
                                                elif (residue.sequence_number<seg_ends['i2b'] and residue.sequence_number>seg_ends['i2e']):
                                                    if debug: print("Remove former GN for ",residue.generic_number)
                                                    residue.display_generic_number = None
                                                    residue.generic_number = None
                                            else:
                                                if residue.generic_number is not None:
                                                    if debug: print("Remove former GN for ",residue.generic_number)
                                                    residue.display_generic_number = None
                                                    residue.generic_number = None
                                                if residue.sequence_number<=seg_ends['3e']:
                                                    residue.protein_segment = self.segments['TM3']
                                                elif residue.sequence_number>=seg_ends['4b']:
                                                    residue.protein_segment = self.segments['TM4']
                                        elif residue.protein_segment.slug=='TM4':
                                            if seg_ends['4b']!='-' and seg_ends['4e']!='-':
                                                if residue.sequence_number<seg_ends['4b']:
                                                    residue.protein_segment = self.segments['ICL2']
                                                    residue.display_generic_number = None
                                                    residue.generic_number = None
                                                elif residue.sequence_number>seg_ends['4e']:
                                                    residue.protein_segment = self.segments['ECL2']
                                                    residue.display_generic_number = None
                                                    residue.generic_number = None
                                        elif residue.protein_segment.slug=='ECL2':
                                            if seg_ends['e2b']!='-' and seg_ends['e2e']!='-' and seg_ends['4e']!='-':
                                                if residue.sequence_number<seg_ends['e2b'] and residue.sequence_number<=seg_ends['4e']:
                                                    residue.protein_segment = self.segments['TM4']
                                                elif residue.sequence_number>seg_ends['e2e'] and residue.sequence_number>=seg_ends['5b']:
                                                    residue.protein_segment = self.segments['TM5']
                                            else:
                                                if residue.sequence_number<=seg_ends['4e']:
                                                    residue.protein_segment = self.segments['TM4']
                                                elif residue.sequence_number>=seg_ends['5b']:
                                                    residue.protein_segment = self.segments['TM5']
                                        elif residue.protein_segment.slug=='TM5':
                                            if seg_ends['5b']!='-' and seg_ends['5e']!='-':
                                                if residue.sequence_number<seg_ends['5b']:
                                                    residue.protein_segment = self.segments['ECL2']
                                                    residue.display_generic_number = None
                                                    residue.generic_number = None
                                                elif residue.sequence_number>seg_ends['5e']:
                                                    residue.protein_segment = self.segments['ICL3']
                                                    residue.display_generic_number = None
                                                    residue.generic_number = None
                                        elif residue.protein_segment.slug=='ICL3':
                                            if residue.sequence_number<=seg_ends['5e']:
                                                residue.protein_segment = self.segments['TM5']
                                            elif residue.sequence_number>=seg_ends['6b']:
                                                residue.protein_segment = self.segments['TM6']
                                        elif residue.protein_segment.slug=='TM6':
                                            if seg_ends['6b']!='-' and seg_ends['6e']!='-':
                                                if residue.sequence_number<seg_ends['6b']:
                                                    residue.protein_segment = self.segments['ICL3']
                                                    residue.display_generic_number = None
                                                    residue.generic_number = None
                                                elif residue.sequence_number>seg_ends['6e']:
                                                    residue.protein_segment = self.segments['ECL3']
                                                    residue.display_generic_number = None
                                                    residue.generic_number = None
                                        elif residue.protein_segment.slug=='TM7':
                                            if seg_ends['7b']!='-' and seg_ends['7e']!='-':
                                                if residue.sequence_number<seg_ends['7b']:
                                                    residue.protein_segment = self.segments['ECL3']
                                                    residue.display_generic_number = None
                                                    residue.generic_number = None
                                                elif residue.sequence_number>seg_ends['7e']:
                                                    if seg_ends['8b']!='-':
                                                        if residue.sequence_number<seg_ends['8b']:
                                                            residue.protein_segment = self.segments['ICL4']
                                                            residue.display_generic_number = None
                                                            residue.generic_number = None
                                                    else:
                                                        residue.protein_segment = self.segments['C-term']
                                                        residue.display_generic_number = None
                                                        residue.generic_number = None
                                        elif residue.protein_segment.slug=='H8':
                                            if seg_ends['8b']!='-' and seg_ends['8e']!='-':
                                                if residue.sequence_number<seg_ends['8b']:
                                                    residue.protein_segment = self.segments['ICL4']
                                                    residue.display_generic_number = None
                                                    residue.generic_number = None
                                                elif residue.sequence_number>seg_ends['8e']:
                                                    residue.protein_segment = self.segments['C-term']
                                                    residue.display_generic_number = None
                                                    residue.generic_number = None
                                            else:
                                                residue.protein_segment = self.segments['C-term']
                                                residue.display_generic_number = None
                                                residue.generic_number = None
                                        elif residue.protein_segment.slug=='C-term':
                                            if seg_ends['8e']!='-':
                                                if residue.sequence_number<=seg_ends['8e']:
                                                    residue.protein_segment = self.segments['H8']

                                        residue.sequence_number = int(check.strip())
                                        #HAVE TO RESET SO IT FITS FOR INTERACTION SCRIPT

                                        if residue.protein_segment != wt_r.protein_segment:
                                            residue.display_generic_number = None
                                            residue.generic_number = None
                                            generic_change += 1
                                            # print(residue.protein_segment.slug[0:2])
                                            if residue.protein_segment.slug[0:2]=="TM" or 1==1:
                                                if debug: print(structure.protein_conformation.protein.entry_name,residue.amino_acid,"XTAL POS",residue.sequence_number, "WT POS",wt_r.sequence_number,"XTAL:",residue.protein_segment,"WT:",wt_r.protein_segment)

                                                # FIX ME
                                                # self.logger.info('No SEG ENDS info for {}'.format(structure.protein_conformation.protein.entry_name,residue.amino_acid,residue.sequence_number, wt_r.sequence_number,residue.protein_segment,wt_r.protein_segment))
                                                pass
                                    if structure.pdb_code.index=='5LWE' and 209<=residue.sequence_number<=244:
                                        residue.display_generic_number = None
                                        residue.generic_number = None
                                        residue.protein_segment = None
                                    # print(residue.sequence_number, "(",int(check.strip()),")",residue.amino_acid,residue.protein_segment,wt_r.amino_acid,wt_r.sequence_number) #sanity check

                            else:
                                #print('wierd error') #sanity check
                                residue.display_generic_number = None
                                residue.generic_number = None
                                residue.protein_segment = None



                            # print('aa ',residue.sequence_number,residue.amino_acid,residue.display_generic_number, residue.protein_segment)
                            if residue.protein_segment==None:
                                pass
                            else:
                                #print('inserted',residue.sequence_number) #sanity check
                                # residue.save()
                                residues_bulk.append(residue)
                                rotamer_data, created = PdbData.objects.get_or_create(pdb=temp)
                                #rotamer_data_bulk.append(PdbData(pdb=temp))
                                missing_atoms = False
                                if rotamer_data.pdb.startswith('COMPND'):
                                    lines = len(rotamer_data.pdb.split('\n'))-2
                                else:
                                    lines = len(rotamer_data.pdb.split('\n'))
                                if lines<atom_num_dict[residue.amino_acid]:
                                    missing_atoms = True
                                rotamer_data_bulk.append([rotamer_data, missing_atoms])
                                # rotamer, created = Rotamer.objects.get_or_create(residue=residue, structure=structure, pdbdata=rotamer_data)
                                #rotamer_bulk.append(Rotamer(residue=residue, structure=structure, pdbdata=rotamer_data))

                        temp = "" #start new line for rotamer
                        check = pdblines[i+1][22:26].strip()

                    check = pdblines[i+1][22:26].strip()
                chain = line[21]
                residue_name = line[17:20].title() #use title to get GLY to Gly so it matches

        ns = settings.DEFAULT_NUMBERING_SCHEME
        ns_obj = ResidueNumberingScheme.objects.get(slug=ns)
        scheme = structure.protein_conformation.protein.residue_numbering_scheme
        segments_present = []
        for res in residues_bulk:
            if res.protein_segment:
                if res.protein_segment.slug not in segments_present:
                    segments_present.append(res.protein_segment.slug)
                if res.generic_number==None and (res.protein_segment.category == "helix" or res.missing_gn): # residue.missing_gn
                    if (res.protein_segment==prev_segment):
                        # print(res, prev_gn, res.protein_segment)
                        gn_split = prev_gn.split("x")
                        new_gn = gn_split[0]+"x"+str(int(gn_split[1])+1)

                        display_split=prev_display.split("x")
                        seq_split = display_split[0].split(".")

                        new_display = seq_split[0]+"."+str(int(seq_split[1])+1)+"x"+str(int(display_split[1])+1)
                        new_equivalent = seq_split[0]+"x"+str(int(display_split[1])+1)

                        if debug: print("Added Generic Number for",res.sequence_number,": GN",new_gn," Display",new_display)

                        gn, created = ResidueGenericNumber.objects.get_or_create(
                                scheme=ns_obj, label=new_gn, protein_segment=res.protein_segment)
                        display_gn, created = ResidueGenericNumber.objects.get_or_create(
                                scheme=scheme, label=new_display, protein_segment=res.protein_segment)

                        try:
                            gn_equivalent, created = ResidueGenericNumberEquivalent.objects.get_or_create(
                                default_generic_number=gn,
                                scheme=scheme,
                                defaults={'label': new_equivalent})
                        except IntegrityError:
                            gn_equivalent = ResidueGenericNumberEquivalent.objects.get(
                                default_generic_number=gn,
                                scheme=scheme)

                        res.generic_number = gn
                        res.display_generic_number = display_gn

                        prev_gn = new_gn
                        prev_display = new_display
                        prev_segment = res.protein_segment
                else:
                    if res.generic_number:
                        prev_gn = res.generic_number.label
                        prev_display = res.display_generic_number.label
                    else:
                        prev_gn = None
                        prev_display = None
                    prev_segment = res.protein_segment

        for res in reversed(residues_bulk):
            if res.protein_segment:
                if res.generic_number==None and (res.protein_segment.category == "helix" or res.missing_gn): # residue.missing_gn
                    if (res.protein_segment==prev_segment):
                        gn_split = prev_gn.split("x")
                        new_gn = gn_split[0]+"x"+str(int(gn_split[1])-1)

                        display_split=prev_display.split("x")
                        seq_split = display_split[0].split(".")

                        new_display = seq_split[0]+"."+str(int(seq_split[1])-1)+"x"+str(int(display_split[1])-1)
                        new_equivalent = seq_split[0]+"x"+str(int(display_split[1])-1)

                        if debug: print("Added Generic Number for",res.sequence_number,": GN",new_gn," Display",new_display)

                        gn, created = ResidueGenericNumber.objects.get_or_create(
                                scheme=ns_obj, label=new_gn, protein_segment=res.protein_segment)
                        display_gn, created = ResidueGenericNumber.objects.get_or_create(
                                scheme=scheme, label=new_display, protein_segment=res.protein_segment)

                        try:
                            gn_equivalent, created = ResidueGenericNumberEquivalent.objects.get_or_create(
                                default_generic_number=gn,
                                scheme=scheme,
                                defaults={'label': new_equivalent})
                        except IntegrityError:
                            gn_equivalent = ResidueGenericNumberEquivalent.objects.get(
                                default_generic_number=gn,
                                scheme=scheme)

                        res.generic_number = gn
                        res.display_generic_number = display_gn

                        prev_gn = new_gn
                        prev_display = new_display
                        prev_segment = res.protein_segment
                else:
                    if res.generic_number:
                        prev_gn = res.generic_number.label
                        prev_display = res.display_generic_number.label
                    else:
                        prev_gn = None
                        prev_display = None
                    prev_segment = res.protein_segment
        bulked_res = Residue.objects.bulk_create(residues_bulk)
        #bulked_rot = PdbData.objects.bulk_create(rotamer_data_bulk)
        bulked_rot = rotamer_data_bulk

        rotamer_bulk = []
        for i,res in enumerate(bulked_res):
            rotamer_bulk.append(Rotamer(residue=res, structure=structure, pdbdata=bulked_rot[i][0],
                                        missing_atoms=bulked_rot[i][1]))

        Rotamer.objects.bulk_create(rotamer_bulk)
        #
        # for i in bulked:
        #     print(i.pk)
        if debug: print("WT",structure.protein_conformation.protein.parent.entry_name,"length",len(parent_seq),structure.pdb_code.index,'length',len(seq),len(mapped_seq),'mapped res',str(mismatch_seq+match_seq+aa_mismatch),'pos mismatch',mismatch_seq,'aa mismatch',aa_mismatch,'not mapped',not_matched,' mapping off, matched on pos,aa',matched_by_pos,"generic_segment_changes",generic_change)
        if (len(segments_present)<8 and 'H8' in segments_present) or len(segments_present)<7:
            print("Present helices:",segments_present)
            print("MISSING HELICES?!")
        if debug: print("===============**================")

        if not os.path.exists(os.sep.join([settings.DATA_DIR, 'structure_data', 'wt_pdb_lookup'])):
            os.makedirs(os.sep.join([settings.DATA_DIR, 'structure_data', 'wt_pdb_lookup']))
        wt_pdb_lookup_folder = os.sep.join([settings.DATA_DIR, 'structure_data', 'wt_pdb_lookup', str(structure) + '.json'])
        if len(wt_pdb_lookup)>0:
            with open(wt_pdb_lookup_folder, 'w') as f2:
                json.dump(wt_pdb_lookup, f2)
        elif os.path.exists(os.sep.join([settings.DATA_DIR, 'structure_data', 'wt_pdb_lookup', str(structure) + '.json'])):
            os.remove(os.sep.join([settings.DATA_DIR, 'structure_data', 'wt_pdb_lookup', str(structure) + '.json']))
        return None


    def build_contact_network(self, pdb_code):
        try:
            # interacting_pairs, distances  = compute_interactions(pdb_code, save_to_db=True)
            compute_interactions(pdb_code, protein=None, lig=None, do_interactions=True, do_complexes=False, do_peptide_ligand=True, save_to_db=True, file_input=False)
            # compute_interactions(pdb_code, do_interactions=True, do_peptide_ligand=True, save_to_db=True)
        except:
            self.logger.error('Error with computing interactions (%s)' % (pdb_code))
            return

    @staticmethod
    def parsecalculation(pdb_id, data, ligand_name, debug=True, ignore_ligand_preset=False):
        module_dir = '/tmp/interactions'
        web_resource = WebResource.objects.get(slug='pdb')
        web_link, _ = WebLink.objects.get_or_create(web_resource=web_resource, index=pdb_id)
        structure = Structure.objects.filter(pdb_code=web_link)
        if structure.exists():
            structure = Structure.objects.get(pdb_code=web_link)

            if structure.pdb_data is None:
                f = module_dir + "/pdbs/" + pdb_id + ".pdb"
                if os.path.isfile(f):
                    pdbdata, created = PdbData.objects.get_or_create(pdb=open(f, 'r').read())  # does this close the file?
                else:
                    print('quitting due to no pdb in filesystem')
                    quit()
                structure.pdb_data = pdbdata
                structure.save()

            protein = structure.protein_conformation
            lig_key = list(data.keys())[0]

            f = module_dir + "/results/" + pdb_id + "/interaction" + "/" + pdb_id + "_" + lig_key + ".pdb"
            if os.path.isfile(f):
                pdbdata, created = PdbData.objects.get_or_create(pdb=open(f, 'r').read())  # does this close the file?
                print("Found file" + f)
            else:
                print('quitting due to no pdb for fragment in filesystem', f)
                quit()

            lig_db_key = lig_key
            if lig_key!=ligand_name and len(lig_key)==3 and len(ligand_name)==5:
                lig_db_key = ligand_name
            struct_lig_interactions = StructureLigandInteraction.objects.filter(pdb_reference=lig_db_key, structure=structure, annotated=True) #, pdb_file=None
            if struct_lig_interactions.exists():  # if the annotated exists
                try:
                    struct_lig_interactions = struct_lig_interactions.get()
                    struct_lig_interactions.pdb_file = pdbdata
                    ligand = struct_lig_interactions.ligand
                except Exception as msg:
                    print('error with duplication structureligand',lig_db_key,msg)
            elif StructureLigandInteraction.objects.filter(pdb_reference=lig_db_key, structure=structure).exists():
                try:
                    struct_lig_interactions = StructureLigandInteraction.objects.filter(pdb_reference=lig_db_key, structure=structure).get()
                    struct_lig_interactions.pdb_file = pdbdata
                except StructureLigandInteraction.DoesNotExist: #already there
                    struct_lig_interactions = StructureLigandInteraction.objects.filter(pdb_reference=lig_db_key, structure=structure, pdb_file=pdbdata).get()
                ligand = struct_lig_interactions.ligand
            else:  # create ligand and pair
                print(pdb_id, "Skipping interactions with ", pdb_id)
                pass

            struct_lig_interactions.save()

            ResidueFragmentInteraction.objects.filter(structure_ligand_pair=struct_lig_interactions).delete()

            for interaction in data[lig_key]['interactions']:
                aa = interaction[0]
                if aa[-1] != structure.preferred_chain:
                    continue
                aa_single, pos, _ = regexaa(aa)
                residue = check_residue(protein, pos, aa_single)
                f = interaction[1]
                fragment, rotamer = extract_fragment_rotamer(f, residue, structure, ligand)
                if fragment is not None:
                    interaction_type, created = ResidueFragmentInteractionType.objects.get_or_create(
                                                slug=interaction[2],
                                                name=interaction[3],
                                                type=interaction[4], direction=interaction[5])
                    fragment_interaction, created = ResidueFragmentInteraction.objects.get_or_create(
                                                    structure_ligand_pair=struct_lig_interactions,
                                                    interaction_type=interaction_type,
                                                    fragment=fragment, rotamer=rotamer)
        else:
            pass

        # results = sorted(results, key=itemgetter(3), reverse=True)

        return data

    def main_func(self, positions, iterations, count, lock):
        # setting up processes
        # if not positions[1]:
        #     pdbs = self.parsed_structures[positions[0]:]
        # else:
        #     pdbs = self.parsed_structures[positions[0]:positions[1]]
        pdbs = self.parsed_structures.pdb_ids
        while count.value<len(pdbs):
            with lock:
                pdb_id = pdbs[count.value]
                count.value +=1

            sd = self.parsed_structures.structures[pdb_id]
            # is this a representative structure (will be used to guide structure-based alignments)?
            representative = False
            if sd['name'].upper() in self.xtal_representatives:#sd['representative']:
                representative = True

            # only process representative structures on first iteration
            # if not representative and iteration == 1:
            #     continue

            # skip representative structures on second iteration
            # if representative and iteration == 2:
            #     continue

            # is there a construct?
            # if 'construct' not in sd:
            #     self.logger.error('No construct specified, skipping!')
            #     continue

            self.logger.info('Building structure {}'.format(sd['name']))

            # does the construct exist?
            try:
                con = Protein.objects.get(entry_name=sd['name'].lower())
            except Protein.DoesNotExist:
                print('BIG ERROR Construct {} does not exists, skipping!'.format(sd['name'].lower()))
                self.logger.error('Construct {} does not exists, skipping!'.format(sd['name'].lower()))
                continue

            # create a structure record
            try:
                s = Structure.objects.get(protein_conformation__protein=con)

                # If update_flag is true then update existing structures
                # Otherwise only make new structures
                if not self.incremental_mode:
                    s = s.delete()
                    s = Structure()
                else:
                    if not s.build_check:
                        s = s.delete()
                        s = Structure()
                    else:
                        continue

            except Structure.DoesNotExist:
                s = Structure()

            s.representative = representative

            # protein state
            if 'state' not in sd:
                self.logger.warning('State not defined, using default state {}'.format(
                    settings.DEFAULT_PROTEIN_STATE))
                state = settings.DEFAULT_STATE.title()
            else:
                state = sd['state']
            state_slug = slugify(state)
            try:
                ps, created = ProteinState.objects.get_or_create(slug=state_slug, defaults={'name': state})
                if created:
                    self.logger.info('Created protein state {}'.format(ps.name))
            except IntegrityError:
                ps = ProteinState.objects.get(slug=state_slug)
            s.state = ps
            s.author_state = ps

            # # xtal activation value aka Delta Distance (Å)
            # if 'distance' not in sd:
            #     self.logger.warning('Delta distance not defined, using default value {}'.format(None))
            #     distance = None
            # else:
            #     distance = sd['distance']
            # s.distance = distance

            # protein conformation
            try:
                s.protein_conformation = ProteinConformation.objects.get(protein=con)
            except ProteinConformation.DoesNotExist:
                self.logger.error('Protein conformation for construct {} does not exists'.format(con))
                continue
            if s.protein_conformation.state is not state:
                ProteinConformation.objects.filter(protein=con).update(state=ps)

            # get the PDB file and save to DB
            sd['pdb'] = sd['name'].upper()
            if not os.path.exists(self.pdb_data_dir):
                os.makedirs(self.pdb_data_dir)

            pdb_path = os.sep.join([self.pdb_data_dir, sd['pdb'] + '.pdb'])
            if not os.path.isfile(pdb_path):
                self.logger.info('Fetching PDB file {}'.format(sd['pdb']))
                url = 'http://www.rcsb.org/pdb/files/%s.pdb' % sd['pdb']
                pdbdata_raw = urlopen(url).read().decode('utf-8')
                with open(pdb_path, 'w') as f:
                    f.write(pdbdata_raw)
            else:
                with open(pdb_path, 'r') as pdb_file:
                    pdbdata_raw = pdb_file.read()

            pdbdata, created = PdbData.objects.get_or_create(pdb=pdbdata_raw)
            s.pdb_data = pdbdata

            self.parsed_pdb = PDBParser(PERMISSIVE=True, QUIET=True).get_structure('ref', pdb_path)[0]

            # UPDATE HETSYN with its PDB reference instead + GRAB PUB DATE, PMID, DOI AND RESOLUTION
            hetsyn = {}
            hetsyn_reverse = {}
            for line in pdbdata_raw.splitlines():
                if line.startswith('HETSYN'):
                    m = re.match("HETSYN[\s]+([\w]{3})[\s]+(.+)",line) ### need to fix bad PDB formatting where col4 and col5 are put together for some reason -- usually seen when the id is +1000
                    if (m):
                        hetsyn[m.group(2).strip()] = m.group(1).upper()
                        hetsyn_reverse[m.group(1)] = m.group(2).strip().upper()
                if line.startswith('HETNAM'):
                    m = re.match("HETNAM[\s]+([\w]{3})[\s]+(.+)",line) ### need to fix bad PDB formatting where col4 and col5 are put together for some reason -- usually seen when the id is +1000
                    if (m):
                        hetsyn[m.group(2).strip()] = m.group(1).upper()
                        hetsyn_reverse[m.group(1)] = m.group(2).strip().upper()
                if line.startswith('REVDAT   1'):
                    sd['publication_date'] = line[13:22]
                if line.startswith('JRNL        PMID'):
                    sd['pubmed_id'] = line[19:].strip()
                if line.startswith('JRNL        DOI'):
                    sd['doi_id'] = line[19:].strip()

            if pdb_id=='8FEG':
                sd['doi_id'] = '10.1038/s41467-023-43718-w'

            if len(hetsyn) == 0:
                self.logger.info("PDB file contained NO hetsyn")

            with open(pdb_path,'r') as header:
                header_dict = parse_pdb_header(header)
            sd['publication_date'] = header_dict['release_date']
            if str(header_dict['resolution']).strip()!='None':
                sd['resolution'] = str(header_dict['resolution']).strip()
            sd['structure_method'] = header_dict['structure_method']

            # structure type
            if 'structure_method' in sd and sd['structure_method']:
                if sd['structure_method']=='unknown':
                    sd['structure_method'] = self.exp_method_dict[sd['method_from_file']]

                structure_type = sd['structure_method'].capitalize()
                structure_type_slug = slugify(sd['structure_method'])
                if sd['pdb']=='6ORV':
                    structure_type_slug = 'electron-microscopy'
                elif sd['pdb'] in ['6YVR','6Z4Q','6Z4S','6Z4V','6Z66','6Z8N','6ZA8','6ZIN','7B6W']:
                    structure_type_slug = 'x-ray-diffraction'

                try:
                    st, created = StructureType.objects.get_or_create(slug=structure_type_slug,
                        defaults={'name': structure_type})
                    if created:
                        self.logger.info('Created structure type {}'.format(st))
                except IntegrityError:
                    st = StructureType.objects.get(slug=structure_type_slug)
                s.structure_type = st
            else:
                self.logger.warning('No structure type specified in PDB file {}'.format(sd['pdb']))

            matched = 0
            if 'ligand' in sd and sd['ligand']:
                if isinstance(sd['ligand'], list):
                    ligands = sd['ligand']
                else:
                    ligands = [sd['ligand']]
                for ligand in ligands:
                    if 'name' in ligand:
                        if ligand['name'].upper() in hetsyn:
                            self.logger.info('Ligand {} matched to PDB records'.format(ligand['name']))
                            matched = 1
                            ligand['name'] = hetsyn[ligand['name'].upper()]
                        elif ligand['name'].upper() in hetsyn_reverse:
                            matched = 1

            if matched==0 and len(hetsyn)>0:
                self.logger.info('No ligand names found in HET in structure {}'.format(sd['pdb']))

            # REMOVE? can be used to dump structure files with updated ligands
            # yaml.dump(sd, open(source_file_path, 'w'), indent=4)

            # pdb code
            if 'pdb' in sd:
                web_resource = WebResource.objects.get(slug='pdb')
                s.pdb_code, created = WebLink.objects.get_or_create(index=sd['pdb'], web_resource=web_resource)
            else:
                self.logger.error('PDB code not specified for structure {}, skipping!'.format(sd['pdb']))
                continue

            # insert into plain text fields
            if 'preferred_chain' in sd:
                if '.' in sd['preferred_chain']:
                    pref_c = sd['preferred_chain'].split('.')[0]
                else:
                    pref_c = sd['preferred_chain']
                s.preferred_chain = pref_c
            else:
                self.logger.warning('Preferred chain not specified for structure {}'.format(sd['pdb']))
            if 'resolution' in sd:
                s.resolution = float(sd['resolution'])
            else:
                self.logger.warning('Resolution not specified for structure {}'.format(sd['pdb']))

            ### Publication date - if pdb file is incorrect, fetch from structures.csv
            if 'publication_date' in sd:
                s.publication_date = sd['publication_date']
                if int(s.publication_date[:4])<1990:
                    s.publication_date = sd['date_from_file']
                    print('WARNING: publication date for {} is incorrect ({}), switched to ({}) from structures.csv'.format(s, sd['publication_date'], sd['date_from_file']))
            else:
                self.logger.warning('Publication date not specified for structure {}'.format(sd['pdb']))

            # publication
            try:
                if 'doi_id' in sd:
                    s.publication = Publication.get_or_create_from_doi(sd['doi_id'])
                elif 'pubmed_id' in sd:
                    s.publication = Publication.get_or_create_from_pubmed(sd['pubmed_id'])
            except:
                self.logger.error('Error saving publication'.format(sd['pdb']))

            # if sd['pdb'] in self.xtal_seg_ends and not self.incremental_mode:
            s.annotated = True
            # else:
            #     s.annotated = False

            s.refined = False
            s.stats_text = None

            # save structure before adding M2M relations
            s.save()
            # StructureLigandInteraction.objects.filter(structure=s).delete()

            # ligands
            peptide_chain = ""
            if self.debug:
                print(sd)
            if 'ligand' in sd and sd['ligand'] and sd['ligand']!='None':
                if isinstance(sd['ligand'], list):
                    ligands = sd['ligand']
                else:
                    ligands = [sd['ligand']]
                for ligand in ligands:
                    l = False
                    peptide_chain = ""
                    if ligand['chain']!='':
                        peptide_chain = ligand['chain']
                        # ligand['name'] = 'pep'
                    if ligand['name'] and ligand['name'] != 'None': # some inserted as none.
                        ligand['type'] = ligand['type'].lower().strip()
                        # use annoted ligand type or default type
                        if ligand['type']:
                            lt, created = LigandType.objects.get_or_create(slug=slugify(ligand['type']),
                                defaults={'name': ligand['type']})
                        else:
                            lt, created = LigandType.objects.get_or_create(
                                slug=slugify(default_ligand_type), defaults={'name': default_ligand_type})
                        # set pdb reference for structure-ligand interaction
                        # if ligand['type'] in ['peptide','protein']:
                        #     pdb_reference = 'pep'
                        #     db_lig = Ligand.objects.filter(name=ligand['title'])
                        # else:
                        #     pdb_reference = ligand['name']
                        #     if ligand['name']=='RET':
                        #         db_lig = Ligand.objects.filter(name=ligand['title'])
                        #     else:
                        #         db_lig = Ligand.objects.filter(pdbe=ligand['name'])
                        #
                        # # check if ligand exists already
                        # if len(db_lig)>0:
                        #     l = db_lig[0]
                        # else:
                        # use ligand title, if specified
                        if 'title' in ligand and ligand['title']:
                            ligand_title = ligand['title']
                        else:
                            ligand_title = ligand['name']

                        # Adding the PDB three-letter code
                        ids = {}
                        pdb_reference = ligand['name']

                        # reformat pdb_reference to string if the 3 letter code is a number and gets formatted as float
                        if len(pdb_reference)>3 and '.' in pdb_reference:
                            pdb_reference = pdb_reference.split('.')[0]

                        if ligand['name'] != "pep" and ligand['name'] != "apo":
                            ids["pdb"] = ligand['name']

                        # use pubchem_id
                        if 'pubchemId' in ligand and ligand['pubchemId'] and ligand['pubchemId'] != 'None':
                            # TODO expand list and expand structure annotation
                            if isinstance(ligand['pubchemId'], str) and ligand['pubchemId'].startswith("gtoplig:"):
                                ids["gtoplig"] = ligand['pubchemId'][8:]
                            else:
                                ids["pubchem"] = ligand['pubchemId']
                        elif "pdb" in ids:
                            # match PDB via UniChem
                            uc_entries = match_id_via_unichem("pdb", ids["pdb"])
                            #print("UNICHEM matching", ids["pdb"])
                            for entry in uc_entries:
                                if entry["type"] not in ids:
                                    ids[entry["type"]] = entry["id"]
                        # sequence
                        if peptide_chain in self.parsed_pdb:
                            seq = ''
                            for res in self.parsed_pdb[peptide_chain]:
                                one_letter = Polypeptide.protein_letters_3to1.get(res.get_resname())
                                if not one_letter:
                                    if res.get_resname() in self.unnatural_amino_acids:
                                        one_letter = self.unnatural_amino_acids[res.get_resname()]
                                    else:
                                        print('WARNING: {} residue in structure {} is missing from unnatural amino acid definitions (data/protwis/gpcr/residue_data/unnatural_amino_acids.yaml)'.format(res, s))
                                        continue
                                seq+=one_letter
                            ids['sequence'] = seq

                        with lock:
                            l = get_or_create_ligand(ligand_title, ids, ligand['type'])
                        # Create LigandPeptideStructure object to store chain ID for peptide ligands - supposed to b TEMP
                        if ligand['type'] in ['peptide','protein']:
                            lps, created = LigandPeptideStructure.objects.get_or_create(structure=s, ligand=l, chain=peptide_chain)
                    else:
                        continue

                    # structure-ligand interaction
                    if l and ligand['role']:
                        role_slug = slugify(ligand['role'])
                        try:
                            lr, created = LigandRole.objects.get_or_create(slug=role_slug,
                            defaults={'name': ligand['role']})
                            if created:
                                self.logger.info('Created ligand role {}'.format(ligand['role']))
                        except IntegrityError:
                            lr = LigandRole.objects.get(slug=role_slug)

                        i, created = StructureLigandInteraction.objects.get_or_create(structure=s,
                            ligand=l, ligand_role=lr, annotated=True,
                            defaults={'pdb_reference': pdb_reference})
                        if i.pdb_reference != pdb_reference:
                            i.pdb_reference = pdb_reference
                            i.save()

            # protein anomalies
            anomaly_entry = self.xtal_anomalies[s.protein_conformation.protein.parent.entry_name]
            segment_codes = {'1':'TM1','12':'ICL1','2':'TM2','23':'ECL1','3':'TM3','34':'ICL2','4':'TM4','5':'TM5','6':'TM6','7':'TM7'}
            all_bulges, all_constrictions = OrderedDict(), OrderedDict()
            for key, val in anomaly_entry.items():
                if 'x' not in val:
                    continue
                if key[0] not in segment_codes:
                    continue
                segment = segment_codes[key.split('x')[0]]
                if len(key.split('x')[1])==3:
                    try:
                        all_bulges[segment] = all_bulges[segment]+[val]
                    except:
                        all_bulges[segment] = [val]
                else:
                    try:
                        all_constrictions[segment] = all_constrictions[segment]+[val]
                    except:
                        all_constrictions[segment] = [val]

            scheme = s.protein_conformation.protein.residue_numbering_scheme
            if len(all_bulges)>0:
                pa_slug = 'bulge'
                try:
                    pab, created = ProteinAnomalyType.objects.get_or_create(slug=pa_slug, defaults={
                        'name': 'Bulge'})
                    if created:
                        self.logger.info('Created protein anomaly type {}'.format(pab))
                except IntegrityError:
                    pab = ProteinAnomalyType.objects.get(slug=pa_slug)

                for segment, bulges in all_bulges.items():
                    for bulge in bulges:
                        try:
                            gn, created = ResidueGenericNumber.objects.get_or_create(label=bulge,
                                scheme=scheme, defaults={'protein_segment': ProteinSegment.objects.get(
                                slug=segment)})
                            if created:
                                self.logger.info('Created generic number {}'.format(gn))
                        except IntegrityError:
                            gn =  ResidueGenericNumber.objects.get(label=bulge, scheme=scheme)

                        try:
                            pa, created = ProteinAnomaly.objects.get_or_create(anomaly_type=pab,
                                generic_number=gn)
                            if created:
                                self.logger.info('Created protein anomaly {}'.format(pa))
                        except IntegrityError:
                            pa, created = ProteinAnomaly.objects.get(anomaly_type=pab, generic_number=gn)

                        s.protein_anomalies.add(pa)
            if len(all_constrictions)>0:
                pa_slug = 'constriction'
                try:
                    pac, created = ProteinAnomalyType.objects.get_or_create(slug=pa_slug, defaults={
                        'name': 'Constriction'})
                    if created:
                        self.logger.info('Created protein anomaly type {}'.format(pac))
                except IntegrityError:
                    pac = ProteinAnomalyType.objects.get(slug=pa_slug)

                for segment, constrictions in all_constrictions.items():
                    for constriction in constrictions:
                        try:
                            gn, created = ResidueGenericNumber.objects.get_or_create(label=constriction,
                                scheme=scheme, defaults={'protein_segment': ProteinSegment.objects.get(
                                slug=segment)})
                            if created:
                                self.logger.info('Created generic number {}'.format(gn))
                        except IntegrityError:
                            gn =  ResidueGenericNumber.objects.get(label=constriction, scheme=scheme)

                        try:
                            pa, created = ProteinAnomaly.objects.get_or_create(anomaly_type=pac,
                                generic_number=gn)
                            if created:
                                self.logger.info('Created protein anomaly {}'.format(pa))
                        except IntegrityError:
                            pa, created = ProteinAnomaly.objects.get(anomaly_type=pac, generic_number=gn)

                        s.protein_anomalies.add(pa)

            # stabilizing agents, FIXME - redesign this!
            # fusion proteins moved to constructs, use this for G-proteins and other agents?
            aux_proteins = []
            # if 'signaling_protein' in sd and sd['signaling_protein'] and sd['signaling_protein'] != 'None':
            #     aux_proteins.append('signaling_protein')
            if 'auxiliary_protein' in sd and sd['auxiliary_protein'] and sd['auxiliary_protein'] != 'None':
                aux_proteins.append('auxiliary_protein')
            for index in aux_proteins:
                if isinstance(sd[index], list):
                    aps = sd[index]
                else:
                    aps = [sd[index]]
                for aux_protein in aps:
                    aux_protein_slug = slugify(aux_protein)[:50]
                    try:
                        sa, created = StructureStabilizingAgent.objects.get_or_create(
                            slug=aux_protein_slug, defaults={'name': aux_protein})
                    except IntegrityError:
                        sa = StructureStabilizingAgent.objects.get(slug=aux_protein_slug)
                    s.stabilizing_agents.add(sa)

            # save structure
            s.save()
            #Delete previous interaction data to prevent errors.
            ResidueFragmentInteraction.objects.filter(structure_ligand_pair__structure=s).delete()
            #Remove previous Rotamers/Residues to prepare repopulate
            Fragment.objects.filter(structure=s).delete()
            Rotamer.objects.filter(structure=s).delete()
            Residue.objects.filter(protein_conformation=s.protein_conformation).delete()

            d = {}

            try:
                current = time.time()
                #protein = Protein.objects.filter(entry_name=s.protein_conformation).get()
                d = fetch_pdb_info(sd['pdb'],con, preferred_chain=s.preferred_chain)
                #delete before adding new
                #Construct.objects.filter(name=d['construct_crystal']['pdb_name']).delete()
                # add_construct(d)
                end = time.time()
                diff = round(end - current,1)
                self.logger.info('construction calculations done for {}. {} seconds.'.format(
                            s.protein_conformation.protein.entry_name, diff))
            except Exception as msg:
                print(msg)
                print('ERROR WITH CONSTRUCT FETCH {}'.format(sd['pdb']))
                self.logger.error('ERROR WITH CONSTRUCT FETCH for {}'.format(sd['pdb']))

                self.construct_errors.append(s)

            try:
                current = time.time()
                self.create_rotamers(s,pdb_path,d)
                # residue_errors = sbc.check_rotamers(s.pdb_code.index)
                # if len(residue_errors)>0:
                #     raise Exception('Error with rotamer check: {}'.format(residue_errors))
                end = time.time()
                diff = round(end - current,1)
                self.logger.info('Create resides/rotamers done for {}. {} seconds.'.format(
                            s.protein_conformation.protein.entry_name, diff))
            except Exception as msg:
                print(msg)
                print('ERROR WITH ROTAMERS {}'.format(sd['pdb']))
                self.logger.error('Error with rotamers for {}'.format(sd['pdb']))

                self.rotamer_errors.append(s)

            try:
                s.protein_conformation.generate_sites()
            except:
                pass

            if self.run_contactnetwork:
                try:
                    current = time.time()
                    self.build_contact_network(sd['pdb'])
                    end = time.time()
                    diff = round(end - current,1)
                    self.logger.info('Create contactnetwork done for {}. {} seconds.'.format(s.protein_conformation.protein.entry_name, diff))
                except Exception as msg:
                    print(msg)
                    print('ERROR WITH CONTACTNETWORK {}'.format(sd['pdb']))
                    self.logger.error('Error with contactnetwork for {}'.format(sd['pdb']))

                    self.contactnetwork_errors.append(s)

            for ligand in ligands:
                if ligand['type'].strip() in ['small molecule', 'protein', 'peptide'] and ligand['in_structure']:
                    try:
                        current = time.time()
                        peptide_chain = ""
                        if ligand['chain']!='':
                            peptide_chain = ligand['chain']
                        # mypath = '/tmp/interactions/results/' + sd['pdb'] + '/output'
                        # if not os.path.isdir(mypath):
                        #     #Only run calcs, if not already in temp
                        # runcalculation(sd['pdb'],peptide_chain)
                        data_results = runcalculation_2022(sd['pdb'], peptide_chain)
                        if 'NAG' in data_results:
                            del data_results['NAG']
                        self.parsecalculation(sd['pdb'], data_results, ligand['name'], False)
                        end = time.time()
                        diff = round(end - current,1)
                        print('Interaction calculations done for {}. {} seconds.'.format(
                                    s.protein_conformation.protein.entry_name, diff))
                        self.logger.info('Interaction calculations done for {}. {} seconds.'.format(
                                    s.protein_conformation.protein.entry_name, diff))
                    except Exception as msg:
                        print(msg)
                        # print(traceback.format_exc())
                        print('ERROR WITH INTERACTIONS {}'.format(sd['pdb']))
                        self.logger.error('Error parsing interactions output for {}'.format(sd['pdb']))

                        self.interaction_errors.append(s)

                        # print('{} done'.format(sd['pdb']))

            s.build_check = True
            s.save()
