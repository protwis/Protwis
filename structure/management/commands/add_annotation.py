from django.core.management.base import BaseCommand
from django.conf import settings

from protein.models import Protein
from residue.models import Residue
from structure.models import Structure
from structure.functions import ParseStructureCSV, X50Finder
from tools.management.commands.build_structure_angles import NonHetSelect
from contactnetwork.interaction import InteractingPair
from construct.functions import fetch_pdb_info, construct_structure_annotation_override
from build.management.commands.build_human_proteins import Command as Parse

from Bio.PDB import PDBParser, Polypeptide
import logging
import os
import yaml
import Bio
from collections import OrderedDict
import pprint
from datetime import datetime
from urllib.request import urlopen
import json
from copy import deepcopy



starttime = datetime.now()


class Command(BaseCommand):
    help = '''Adds segment end annotations for new structures from structure related csv files into the yaml files. 
              Run it on a complete db after parse_excel_annotations.py was run and the csv files got updated with the 
              updated GPCRdb_structure_info.xlsx'''

    logger = logging.getLogger(__name__)

    xtal_seg_end_file = os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'mod_xtal_segends.yaml'])
    nonxtal_seg_end_file = os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'non_xtal_segends.yaml'])
    anomalies_file = os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'all_anomalies.yaml'])
    with open(anomalies_file, 'r') as f2:
        anomalies = yaml.load(f2, Loader=yaml.FullLoader)
    pdb_data_dir = os.sep.join([settings.DATA_DIR, 'structure_data', 'pdbs'])

    sequence_file = os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'sequences.yaml'])
    with open(sequence_file, 'r') as f:
        gpcr_sequences = yaml.load(f, Loader=yaml.FullLoader)

    default_segends = {'1b':'-','1e':'-','i1b':'-','i1e':'-','2b':'-','2e':'-','e1b':'-','e1e':'-',
                       '3b':'-','3e':'-','i2b':'-','i2e':'-','4b':'-','4e':'-','e2b':'-','e2e':'-',
                       '5b':'-','5e':'-','6b':'-','6e':'-','7b':'-','7e':'-','8b':'-','8e':'-'}

    default_x50s = {'1x':'-','i1x':'-','2x':'-','e1x':'-','3x':'-','i2x':'-','4x':'-','e2x':'-','5x':'-','6x':'-','7x':'-','8x':'-'}

    x50_dict = {'1.50x50':'1x','12.50x50':'i1x','2.50x50':'2x','23.50x50':'e1x','3.50x50':'3x','34.50x50':'i2x','4.50x50':'4x',
                '45.50x50':'e2x','5.50x50':'5x','6.50x50':'6x','7.50x50':'7x','8.50x50':'8x'}

    def add_arguments(self, parser):
        parser.add_argument('--debug',
            action='store_true',
            dest='debug',
            default=False,
            help='Print info for debugging')
        parser.add_argument('--no_save',
            action='store_false',
            dest='no_save',
            default=True,
            help='Does not save annotation')
        parser.add_argument('-s', '--structure',
            dest='structure',
            help='Structure to annotate',
            nargs='+')
        parser.add_argument('--alphafold',
            action='store_true',
            dest='alphafold',
            default=False,
            help='Run annotation on Alphafold models')

    def handle(self, *args, **options):
        self.debug = options['debug']
        self.save_annotation = options['no_save']

        with open(self.nonxtal_seg_end_file, 'r') as f:
            self.nonxtal_seg_ends = yaml.load(f, Loader=yaml.FullLoader)

        out_segends, out_anomalies = {}, {}

        if options['alphafold']:
            if options['structure']:
                structure_input = options['structure']
                if len(structure_input)==1 and os.path.exists(structure_input[0]):
                    path_to_models = structure_input[0]
                    models = os.listdir(path_to_models)
                    p = Parse()
                    missing_segments = {}
                    outdated_models = []
                    c = 0

                    for m in models:
                        c+=1
                        print(m,c)
                        dssp = None
                        res_dict = OrderedDict()
                        accession = m.split('.')[0]
                        up = deepcopy(p.parse_uniprot_file(accession))

                        ### Running x50 finder
                        if up['entry_name'] not in self.nonxtal_seg_ends:
                            x = X50Finder(os.sep.join([path_to_models, m]), self.debug)
                            x50s_gn = x.run()

                            if self.debug:
                                print('Top blast hit:', x.top_hit)
                            temp_anomalies = {}
                            temp_anomalies['UniProt'] = up['entry_name']
                            temp_anomalies['Xtal Templ'] = up['entry_name']
                            temp_anomalies['Xtalised'] = 'NoXtal'
                            temp_anomalies['7x44'] = '7x44'
                            out_anomalies[up['entry_name']] = temp_anomalies
                        else:
                            continue

                        x50s = deepcopy(self.default_x50s)
                        if self.debug:
                            print(x50s_gn)
                        for chain, xs in x50s_gn.items():
                            for gn, x50 in xs.items():
                                x50s[self.x50_dict[gn]] = x50

                        segends = deepcopy(self.default_segends)

                        if up['entry_name'] not in self.gpcr_sequences:
                            self.add_sequence(up['entry_name'], False)

                        dssp = self.dssp(accession, x.biopdb[0][chain], x.biopdb)
    
                        struct_seq = ''
                        for d in dssp:
                            struct_seq+=d[1]

                        if up['sequence']!=struct_seq:
                            print(up['sequence'])
                            print(struct_seq)
                            outdated_models.append(up['entry_name'])
                            continue

                        ### Exceptions
                        if m=='Q8NGY7.pdb':
                            x50s['6x'] = 253
                            x50s['7x'] = '-' 
                            x50s['8x'] = '-'
                        elif m=='Q8NH09.pdb':
                            x50s['8x'] = '-'
                        elif m=='A6NMU1.pdb':
                            x50s['6x'] = 262
                            x50s['7x'] = 286
                            x50s['8x'] = 295
                        elif m=='O14718.pdb':
                            x50s['1x'] = 43
                        ###

                        for lab, x50 in x50s.items():
                            ### Exceptions
                            if m=='A6NET4.pdb' and lab=='5x':
                                segends[lab[0]+'b'] = 192
                                segends[lab[0]+'e'] = 226
                                continue
                            elif m=='A6NMU1.pdb' and lab=='7x':
                                segends[lab[0]+'b'] = 269
                                segends[lab[0]+'e'] = 291
                                continue

                            ### Helixes
                            if lab[0] in ['i','e']:
                                continue
                            b, e = '-', '-'
                            if x50=='-':
                                continue
                            if dssp[x50-1][2] not in ['H','G','I','T']:
                                x50s[lab] = '-'
                                continue
                            if int(lab[0])>1 and segends[str(int(lab[0])-1)+'e']=='-':
                                continue
                            ### End
                            end_list = dssp[x50-1:]
                            for i, d in enumerate(end_list):
                                if d[2] not in ['H','G','I']:
                                    if d[2]!='T' or (d[2]=='T' and end_list[i+1][2] not in ['H','G','I']):
                                        e = d[0]
                                        break
                            ### Start
                            if lab[0]=='8':
                                b = segends['7e']+1
                            else:
                                rev_list = sorted(dssp[:x50], key=lambda x: (-x[0]))
                                for i, d in enumerate(rev_list):
                                    if d[2] not in ['H','G','I'] and (int(lab[0])==1 or d[0]>segends[str(int(lab[0])-1)+'e']):
                                        if d[2]!='T' or (d[2]=='T' and rev_list[i+1][2] not in ['H','G','I']):
                                            b = d[0]
                                            break
                                if e!='-' and b=='-':
                                    b = segends[str(int(lab[0])-1)+'e']+1

                            segends[lab[0]+'b'] = b
                            segends[lab[0]+'e'] = e

                            
                        ### Exceptions
                        if m=='Q8NH73.pdb':
                            segends['1e'] = 50
                            segends['i1b'] = 51
                            segends['i1e'] = 54
                            segends['2b'] = 55
                            x50s['i1x'] = 53

                        ### ICL1
                        if x50s['i1x']!='-' and segends['1e']!='-' and segends['2b']!='-':
                            has_i1 = False
                            i1_range = range(segends['1e']+1, segends['2b'])
                            struct_i1 = [i1r for i1r in i1_range if i1r in res_dict]
                            non_helical, remove_list = self.get_non_helicals(i1_range, res_dict)
                            if len(non_helical)==len(i1_range):
                                for nh in non_helical:
                                    if self.check_H_bond(res_dict, nh, 'start'):
                                        has_i1 = True
                                        break
                            elif len(struct_i1)!=len(i1_range):
                                pass
                            else:
                                for nh in non_helical:
                                    if self.check_H_bond(res_dict, nh, 'start'):
                                        has_i1 = True
                                        break
                                if len(i1_range)-len(non_helical)>2:
                                    has_i1 = True

                            if has_i1:
                                segends['i1b'] = segends['1e']+1
                                segends['i1e'] = segends['2b']-1
                            else:
                                x50s['i1x'] = '-'

                        ### ECL1
                        if x50s['e1x']!='-' and segends['2e']!='-' and segends['3b']!='-':
                            segends['e1b'] = segends['2e']+1
                            segends['e1e'] = segends['3b']-1

                        ### ICL2
                        if x50s['i2x']!='-' and segends['3e']!='-' and segends['4b']!='-':
                            i2_range = range(segends['3e']+1, segends['4b'])
                            non_helical, remove_list = self.get_non_helicals(i2_range, res_dict)
                            if len(list(i2_range))-len(non_helical)>=4:
                                segends['i2b'] = segends['3e']+1
                                segends['i2e'] = segends['4b']-1
                            else:
                                x50s['i2x'] = '-'
                        else:
                            x50s['i2x'] = '-'
                        
                        ### ECL2
                        if x50s['e2x']!='-' and up['sequence'][x50s['e2x']-1]=='C' and x50s['3x']!='-' and 'C' in up['sequence'][x50s['3x']-26:x50s['3x']-23]:
                            segends['e2b'] = x50s['e2x']
                            segends['e2e'] = x50s['e2x']+2
                        else:
                            x50s['e2x'] = '-'

                        for l, s in segends.items():
                            x50s[l] = s

                        ### Check for no overlapping segends
                        ordered_segments = ['1','i1','2','e1','3','i2','4','e2','5','6','7','8']
                        out = OrderedDict()
                        if self.debug:
                            print(x50s)
                        for i, l in enumerate(ordered_segments):
                            if self.debug:
                                print(l)
                            if x50s[l+'x']!='-' and not x50s[l+'b']<=x50s[l+'x']<=x50s[l+'e']:
                                print('ERROR: number ascension error within segment {}'.format(l))
                                raise AssertionError
                            if i>0 and x50s[l+'b']!='-' and x50s[ordered_segments[i-1]+'e']!='-' and not x50s[ordered_segments[i-1]+'e']<x50s[l+'b']:
                                print('ERROR: segment start {} is not higher than segment end {}'.format(l, ordered_segments[i-1]))
                                raise AssertionError
                            ### Reorder
                            out[l+'b'] = x50s[l+'b']
                            out[l+'x'] = x50s[l+'x']
                            out[l+'e'] = x50s[l+'e']
                            ### Check missing
                            if l[0] not in ['i','e'] and x50s[l+'x']=='-':
                                if m not in missing_segments:
                                    missing_segments[m] = [l]
                                else:
                                    missing_segments[m].append(l)
                        
                        formatted_out = {}
                        for i,j in out.items():
                            formatted_out[i] = j

                        out_segends[up['entry_name']] = formatted_out

                    print('Missing segments')
                    pprint.pprint(missing_segments)
                    print('Outdated models')
                    print(outdated_models)

                    if self.save_annotation:
                        with open(self.nonxtal_seg_end_file, 'a') as f1:
                            yaml.dump(out_segends, f1, default_flow_style=False)

                        with open(self.anomalies_file, 'a') as f2:
                            yaml.dump(out_anomalies, f2, default_flow_style=False)
                else:
                    # either path is incorrect or can be individual accession input FIXME
                    pass

        else:
            self.parsed_structures = ParseStructureCSV()
            self.parsed_structures.parse_ligands()
            self.parsed_structures.parse_nanobodies()
            self.parsed_structures.parse_fusion_proteins()
            self.parsed_structures.parse_ramp()
            self.parsed_structures.parse_grk()
            self.parsed_structures.parse_parent_segends()

            if options['structure']:
                self.structures_to_annotate = options['structure']
            else:
                self.structures_to_annotate = []

            with open(self.xtal_seg_end_file, 'r') as f:
                self.xtal_seg_ends = yaml.safe_load(f)

            segends = OrderedDict()
            # mismatches = {}
            new_unique_receptor_structures = {}
            missing_ligand_info = []

            errors = {}

            for s, data in self.parsed_structures.structures.items():
                try:
                    ### New structures
                    if len(self.structures_to_annotate)>0 and s not in self.structures_to_annotate:
                        continue
                    ### Warnings for missing data
                    if 'ligand' not in self.parsed_structures.structures[s]:
                        print('WARNING: {} missing ligand annotation'.format(s))
                        missing_ligand_info.append(s)
                    else:
                        for l in data['ligand']:
                            if l['name']=='pep' and l['in_structure'] and l['chain']=='':
                                print('WARNING: {} {} peptide ligand missing chain ID'.format(s, l['title']))
                            if l['role']=='':
                                print('WARNING: {} {} ligand missing modality'.format(s, l['title']))
                    if s not in self.xtal_seg_ends or s in self.structures_to_annotate:
                        print(s)
                        segends[s] = deepcopy(self.default_segends)

                        self.download_pdb(s)

                        if data['protein'] not in self.gpcr_sequences:
                            self.add_sequence(data['protein'], True)
                            
                        structure = Bio.PDB.PDBParser(QUIET=True).get_structure(s, self.pdb_path)
                        parent_protein = Protein.objects.get(entry_name=data['protein'])
                        parent_residues = Residue.objects.filter(protein_conformation__protein=parent_protein)
                        parent_seq = parent_protein.sequence

                        # Check for fusion protein
                        fusion_present = False
                        deletions, removed = [], []
                        if 'auxiliary_protein' in data:
                            if len([a for a in data['auxiliary_protein'] if a in self.parsed_structures.fusion_proteins])>0:
                                fusion_present = True
                                try:
                                    d = fetch_pdb_info(s, parent_protein, True, preferred_chain=data['preferred_chain'])
                                    if self.debug:
                                        print('pdb info fetched', datetime.now()- starttime)
                                    if 'deletions' in d:
                                        for del_range in d['deletions']:
                                            if del_range['start']==146 and s=='4K5Y':
                                                #Manual fix for faulty 4K5Y annotation
                                                continue
                                            for i in range(del_range['start'],del_range['end']+1):
                                                deletions.append(i)
                                        #print("Annotation missing WT residues",d['deletions'])
                                    ## Remove segments that arent receptor (tags, fusion etc)
                                    if 'xml_segments' in d:
                                        for seg in d['xml_segments']:
                                            if seg[1]:
                                                # Odd rules to fit everything..
                                                # print(seg[1][0], entry_name)
                                                if seg[1][0]!=data['protein'] and seg[-1]!=True and seg[1][0]!='Uncharacterized protein' and 'receptor' not in seg[1][0]:
                                                    if seg[0].split("_")[1]==data['preferred_chain']:
                                                        for i in seg[6]:
                                                            removed.append(i)

                                    removed, deletions = construct_structure_annotation_override(s, removed, deletions)

                                    if self.debug:
                                        print('Deletions:')
                                        print(deletions)
                                        print('Removed:')
                                        print(removed)
                                except TypeError:
                                    print('ERROR: Fetch pdb data failed {} needs debug check'.format(s))

                        seq = ''
                        res_list = []
                        for res in structure[0][data['preferred_chain']]:
                            ### Skip residues with missing backbone atoms
                            if not res.has_id('N') or not res.has_id('CA') or not res.has_id('C') or not res.has_id('O'):
                                continue
                            ### Skip heteroatom residues - non-receptor residues
                            if res.get_id()[0].startswith('H'):
                                continue
                            try:
                                seq+=Polypeptide.protein_letters_3to1.get(res.get_resname())
                                res_list.append(res)
                            except KeyError:
                                pass

                        # crh = ClosestReceptorHomolog(data['protein'], force=True)
                        # crh.find_closest_receptor_homolog()
                        # max5sim = heapq.nlargest(5, [p.similarity for p in crh.all_proteins[1:]])
                        # max5prot = []
                        # for p in crh.all_proteins:
                        #     if p.protein not in max5prot and p.similarity in max5sim:
                        #         max5prot.append(p.protein)
                        # structures = Structure.objects.filter(protein_conformation__protein__parent__in=max5prot)

                        # print(len(structures))
                        # print(max5prot)

                        # Custom fixes
                        if s=='2I35':
                            seq = seq[:-2]
                        elif s=='3ZEV':
                            seq = seq[:-5]
                        elif s in ['8JD1','8JD2','8JD3','8JD4','8JD5','8IZB','8WPG','8WPU','8WRB']:
                            seq = seq[:-1]

                        if self.debug:
                            print(seq)

                        dssp = self.dssp(s, structure[0][data['preferred_chain']], structure)

                        if fusion_present or s in ['7V68','7V69','7V6A','7W6P','7W7E','8E9W','8E9X','8E9Y','8E9Z','8EA0','7T8X','7T90','7T94','7T96',
                                                   '7TRK','7TRP','7TRQ','7TRS','8IRU','8FX5','8W8R','8W8S','7V9L','8TZQ','8U02','8YN2']:
                            pw2 = Bio.pairwise2.align.localms(parent_seq, seq, 3, -3, -3.5, -1)
                        else:
                            pw2 = Bio.pairwise2.align.localms(parent_seq, seq, 3, -4, -5, -2)

                        ref_seq, temp_seq = str(pw2[0][0]), str(pw2[0][1])
                        ref_i, temp_i = 0, 0
                        res_dict = OrderedDict()
                        wt_pdb_lookup = {}

                        for i, r in enumerate(ref_seq, 1):
                            fusion_res = False
                            if self.debug:
                                print(i,r,temp_seq[i-1])
                            # if r=='-':
                            #     temp_i+=1
                            if temp_seq[i-1]!='-':
                                if s=='7ZI0':
                                    if fusion_present and res_list[temp_i].get_id()[1] in removed:
                                        fusion_res = True
                                else:
                                    if fusion_present and (res_list[temp_i].get_id()[1] in deletions or res_list[temp_i].get_id()[1] in removed):
                                        fusion_res = True
                                if not fusion_res:
                                    res_dict[parent_residues[ref_i].sequence_number] = [res_list[temp_i], dssp[temp_i]]
                                    wt_pdb_lookup[parent_residues[ref_i].sequence_number] = res_list[temp_i].get_id()[1]
                                temp_i+=1
                            if r!='-':
                                ref_i+=1

                        needs_lookup = False
                        if len(wt_pdb_lookup)>0:
                            needs_lookup = True

                        if self.debug:
                            pprint.pprint(wt_pdb_lookup)
                        try:
                            parent_segends = self.parsed_structures.parent_segends[parent_protein.entry_name]
                        except KeyError:
                            # Ortholog not in non_xtal_segends
                            parent_segends = OrderedDict()
                            for seg in parent_residues.order_by('protein_segment').values_list('protein_segment__slug', flat=True).distinct():
                                x50 = parent_residues.filter(protein_segment__slug=seg, generic_number__label__endswith='x50')
                                if len(x50)>0:
                                    x50 = x50[0].sequence_number
                                    b = parent_residues.filter(protein_segment__slug=seg, generic_number__isnull=False)[0].sequence_number
                                    e = parent_residues.filter(protein_segment__slug=seg, generic_number__isnull=False).reverse()[0].sequence_number
                                else:
                                    x50, b, e, = '-','-','-'
                                if seg.startswith('TM') or seg=='H8':
                                    prefix = seg[-1]
                                elif seg.startswith('ICL'):
                                    prefix = 'i'+seg[-1]
                                elif seg.startswith('ECL'):
                                    prefix = 'e'+seg[-1]
                                else:
                                    continue
                                parent_segends[prefix+'b'] = b
                                parent_segends[prefix+'x'] = x50
                                parent_segends[prefix+'e'] = e

                        if self.debug:
                            for i,j in res_dict.items():
                                print(i,j)
                            print(parent_segends)

                        ### Helices
                        for i in range(1,9):
                            if i==8 and parent_segends[str(i)+'x']=='-':
                                print('WARNING: no H8 annotated for wt {}'.format(parent_protein))
                                continue
                            parent_x50 = int(parent_segends[str(i)+'x'])
                            parent_start = int(parent_segends[str(i)+'b'])
                            parent_end = int(parent_segends[str(i)+'e'])

                            ### Start
                            start_range = range(parent_x50, parent_start-1, -1)
                            non_helical, remove_list = self.get_non_helicals(start_range, res_dict)

                            if len(non_helical)==0:
                                while parent_start in res_dict and res_dict[parent_start][1][2]=='H' and parent_start-1 in res_dict and res_dict[parent_start][0].get_id()[1]-res_dict[parent_start-1][0].get_id()[1]==1:
                                    parent_start-=1
                                start_range = range(parent_x50, parent_start-1, -1)
                                non_helical, remove_list = self.get_non_helicals(start_range, res_dict)
                            
                            ### Checking for start non-helical but backbone H-bond residues at start
                            for nh in non_helical:
                                if not self.check_H_bond(res_dict, nh, 'start'):
                                    remove_list.append(nh)
                            try:
                                start = min([k for k in start_range if k not in remove_list])
                            except ValueError:
                                closest_to_x50 = parent_x50
                                break_point, break_i = 50, 0
                                while closest_to_x50 not in res_dict and break_i<=break_point:
                                    closest_to_x50+=1
                                    break_i+=1
                                if break_i==break_point:
                                    start = '-'
                                else:
                                    start = closest_to_x50

                            ### End
                            end_range = range(parent_x50, parent_end+1)
                            non_helical, remove_list = self.get_non_helicals(end_range, res_dict)
                            if len(non_helical)==0:
                                while parent_end in res_dict and res_dict[parent_end][1][2]=='H' and parent_end+1 in res_dict and res_dict[parent_end+1][0].get_id()[1]-res_dict[parent_end][0].get_id()[1]==1:
                                    parent_end+=1
                                end_range = range(parent_x50, parent_end+1)
                                non_helical, remove_list = self.get_non_helicals(end_range, res_dict)

                            ### Checking for end non-helical but backbone H-bond residues at start
                            for nh in non_helical:
                                if not self.check_H_bond(res_dict, nh, 'end'):
                                    remove_list.append(nh)
                            try:
                                end = max([k for k in end_range if k not in remove_list])
                            except ValueError:
                                start = '-'
                                end = '-'

                            if i==8 and start!='-':
                                if start==segends[s]['7e'] and parent_x50-start==3:
                                    segends[s]['7e']-=1
                                elif start-segends[s]['7e']==2:
                                    segends[s]['7e'] = start-1
                            if needs_lookup:
                                if start in wt_pdb_lookup:
                                    start = wt_pdb_lookup[start]
                                if end in wt_pdb_lookup:
                                    end = wt_pdb_lookup[end]

                            ### Custom segends
                            if s=='7VVJ' and i==6:
                                start = 402
                                end = 425
                            elif s=='8JD4' and i==3:
                                start = 628
                                end = 655
                            elif s=='8K4N' and i==7:
                                start = 298
                                end = 326
                            elif s=='8WPG' and i==3:
                                start = 673
                                end = 699
                            elif s in ['8GGC','8GGF'] and i==7:
                                start = 304
                                end = 322
                            elif s=='8TQB' and i==3:
                                start = 640
                                end = 663
                            elif s=='8TR0' and i==3:
                                start = 640
                                end = 664
                            elif s=='8TR2' and i==3:
                                start = 641
                                end = 665
                            elif s=='8YW3' and i==6:
                                start = 346
                                end = 367
                            elif s in ['9C1P','9C2F'] and i==3:
                                start = 673
                                end = 698

                            if i<8 and (start=='-' or end=='-'):
                                print('WARNING: helix {} for {} {} has missing annotation'.format(i, s, parent_protein))

                            segends[s][str(i)+'b'] = start
                            segends[s][str(i)+'e'] = end

                            ### ICL1
                            if i==1 and parent_segends['i1b']!='-':
                                has_i1 = False
                                i1_range = range(int(parent_segends['i1b'])-1, int(parent_segends['i1e'])+1)
                                struct_i1 = [i1r for i1r in i1_range if i1r in res_dict]
                                non_helical, remove_list = self.get_non_helicals(i1_range, res_dict)
                                if len(non_helical)==len(i1_range):
                                    for nh in non_helical:
                                        if self.check_H_bond(res_dict, nh, 'start'):
                                            has_i1 = True
                                            break
                                elif len(struct_i1)!=len(i1_range):
                                    pass
                                else:
                                    for nh in non_helical:
                                        if self.check_H_bond(res_dict, nh, 'start'):
                                            has_i1 = True
                                            break
                                if has_i1:
                                    segends[s]['i1b'] = wt_pdb_lookup[int(parent_segends['i1b'])]
                                    segends[s]['i1e'] = wt_pdb_lookup[int(parent_segends['i1e'])]
                                    if segends[s]['1e']==segends[s]['i1b']:
                                        segends[s]['1e']-=1
                                    if segends[s]['2b']==segends[s]['i1e']:
                                        segends[s]['2b']+=1

                            ### ECL1
                            if i==2 and parent_segends['e1b']!='-':
                                e1_range = range(int(parent_segends['e1b']), int(parent_segends['e1e'])+1)
                                struct_e1 = [e for e in e1_range if e in res_dict]
                                if len(struct_e1)==len(e1_range):
                                    segends[s]['e1b'] = wt_pdb_lookup[int(parent_segends['e1b'])]
                                    segends[s]['e1e'] = wt_pdb_lookup[int(parent_segends['e1e'])]
                                    if segends[s]['2e']==segends[s]['e1b']:
                                        segends[s]['2e']-=1
                                    if segends[s]['3b']==segends[s]['e1e']:
                                        segends[s]['3b']+=1

                            ### ICL2
                            if i==3 and parent_segends['i2b']!='-':
                                i2_range = range(int(parent_segends['i2b']), int(parent_segends['i2e'])+1)
                                non_helical, remove_list = self.get_non_helicals(i2_range, res_dict)
                                if len(non_helical)<4 and len([i for i in i2_range if i in res_dict])==len(i2_range):
                                    segends[s]['i2b'] = wt_pdb_lookup[int(parent_segends['i2b'])]
                                    segends[s]['i2e'] = wt_pdb_lookup[int(parent_segends['i2e'])]
                                    if segends[s]['3e']==segends[s]['i2b']:
                                        segends[s]['3e']-=1
                                    if segends[s]['4b']==segends[s]['i2e']:
                                        segends[s]['4b']+=1

                            ### ECL2
                            if i==4 and parent_segends['e2b']!='-':
                                # e2_range = range(int(parent_segends['e2b']), int(parent_segends['e2e'])+1)
                                tm3_cys = parent_residues.get(generic_number__label='3x25').sequence_number
                                if parent_segends['e2b'] in res_dict and res_dict[parent_segends['e2b']][1][1]=='C' and tm3_cys in res_dict and res_dict[tm3_cys][1][1]=='C':
                                    segends[s]['e2b'] = wt_pdb_lookup[int(parent_segends['e2b'])]
                                    e2e = int(parent_segends['e2e'])
                                    e2e_counter = 0
                                    while e2e not in wt_pdb_lookup and e2e_counter<2:
                                        e2e-=1
                                        e2e_counter+=1
                                    segends[s]['e2e'] = wt_pdb_lookup[e2e]

                            parent_start, parent_x50, parent_end = None, None, None

                        if self.debug:
                            pprint.pprint(segends[s])
                        proteins_in_db = Structure.objects.filter(protein_conformation__protein__parent=parent_protein).exclude(structure_type__slug__startswith='af-')
                        if len(proteins_in_db)==0:
                            if parent_protein not in new_unique_receptor_structures:
                                new_unique_receptor_structures[parent_protein] = [s]
                            else:
                                new_unique_receptor_structures[parent_protein].append(s)
                        # print(new_unique_receptor_structures)

                        ### Check with done structures
                        # for seg, val in segends[s].items():
                        #     if val!=int(parent_segends[seg]) and s not in mismatches:
                        #         mismatches[s] = [seg, int(parent_segends[seg]), val]
                        #     elif val!=int(parent_segends[seg]):
                        #         mismatches[s].append([seg, int(parent_segends[seg]), val])

                        for s, data in segends.items():
                            formatted_out = {}
                            for seg, end in data.items():
                                formatted_out[seg] = end
                            self.xtal_seg_ends[s] = formatted_out

                        ### Sort
                        sorted_keys = list(self.xtal_seg_ends)
                        sorted_keys.sort()
                        sorted_dict = {i: self.xtal_seg_ends[i] for i in sorted_keys}
                        self.xtal_seg_ends = sorted_dict

                        ### Save to file
                        if self.save_annotation:
                            with open(self.xtal_seg_end_file, 'w') as f1:
                                yaml.dump(self.xtal_seg_ends, f1, default_flow_style=False)
                except Exception as msg:
                    errors[s] = msg
            
            print('Missing ligand info:')
            print(missing_ligand_info)
            # pprint.pprint(segends)
            print('New unique receptor structures')
            pprint.pprint(new_unique_receptor_structures)

            print('ERRORS:')
            pprint.pprint(errors)



    @staticmethod
    def get_non_helicals(res_range, residues):
        non_helical, remove_list = [], []
        for j in res_range:
            # print(j)
            if j not in residues:
                remove_list.append(j)
                continue
            # print(j, residues[j])
            if residues[j][1][2] not in ['H','G','I']:
                non_helical.append(j)
        return non_helical, remove_list

    @staticmethod
    def check_H_bond(residues, seqnum, start_or_end, cutoff=3.5):
        if start_or_end=='start':
            bond_partner_i = [3, 4, 5]
            this_atom = 'O'
            other_atom = 'N'
        elif start_or_end=='end':
            bond_partner_i = [-3, -4, -5]
            this_atom = 'N'
            other_atom = 'O'
        this_vector = residues[seqnum][0][this_atom].get_vector()
        for i in bond_partner_i:
            try:
                other_vector = residues[seqnum+i][0][other_atom].get_vector()
            except KeyError:
                continue
            d = other_vector-this_vector
            if start_or_end=='start':
                angle = InteractingPair.verify_hbond_angle(residues[seqnum+i][0], 'N', residues[seqnum][0], 'O')
            elif start_or_end=='end':
                angle = InteractingPair.verify_hbond_angle(residues[seqnum][0], 'N', residues[seqnum+i][0], 'O')
            # print(seqnum, seqnum+i, d.norm(), angle)
            if d.norm()<cutoff and angle:
                return True
        return False

    def add_sequence(self, protein, xtalised):
        url = 'https://rest.uniprot.org/uniprotkb/{}.fasta'.format(protein)
        uniprot_data = urlopen(url).read().decode('utf-8').split('\n')
        accession = uniprot_data[0].split('|')[1]
        sequence = ''.join(uniprot_data[1:])
        print(accession)
        print(sequence)
        if xtalised:
            x = 'Xtal'
        else:
            x = 'NoXtal'
        dic = {'Sequence':sequence, 'UniProt':protein, 'Xtalised':x}
        self.gpcr_sequences[protein.lower()] = dic
        ### Adding canonical sequence from UniProt to sequences.yaml
        if self.save_annotation:
            with open(os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'sequences.yaml']), 'w') as f:
                yaml.dump(self.gpcr_sequences, f)

    def download_pdb(self, pdb_code):
        self.pdb_path = os.sep.join([self.pdb_data_dir, pdb_code + '.pdb'])
        if not os.path.isfile(self.pdb_path):
            self.logger.info('Fetching PDB file {}'.format(pdb_code))
            url = 'http://www.rcsb.org/pdb/files/%s.pdb' % pdb_code
            pdbdata_raw = urlopen(url).read().decode('utf-8')
            with open(self.pdb_path, 'w') as f:
                f.write(pdbdata_raw)
        else:
            with open(self.pdb_path, 'r') as pdb_file:
                pdbdata_raw = pdb_file.read()

    def dssp(self, pdb_code, pchain, structure):
        # DSSP
        filename = "{}_temp.pdb".format(pdb_code)
        pdbio = Bio.PDB.PDBIO()
        pdbio.set_structure(pchain)
        pdbio.save(filename, NonHetSelect())
        if os.path.exists("/env/bin/dssp"):
            dssp = Bio.PDB.DSSP(structure[0], filename, dssp='/env/bin/dssp')
        elif os.path.exists("/env/bin/mkdssp"):
            dssp = Bio.PDB.DSSP(structure[0], filename, dssp='/env/bin/mkdssp')
        elif os.path.exists("/usr/local/bin/mkdssp"):
            dssp = Bio.PDB.DSSP(structure[0], filename, dssp='/usr/local/bin/mkdssp')
        # if self.debug:
        #     for d in dssp:
        #         print(d)
        os.remove(filename)
        return [i for i in dssp]
