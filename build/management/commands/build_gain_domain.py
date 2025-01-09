from django.core.management.base import BaseCommand, CommandError
from build.management.commands.base_build import Command as BaseBuild
from django.conf import settings
from django.db import connection
from django.db import IntegrityError

from protein.models import (Protein, ProteinConformation, ProteinState, ProteinFamily, ProteinAlias,
        ProteinSequenceType, Species, Gene, ProteinSource, ProteinSegment)
from residue.models import (ResidueNumberingScheme, ResidueGenericNumber, Residue, ResidueGenericNumberEquivalent)
from signprot.models import SignprotComplex, SignprotStructure, SignprotStructureExtraProteins
from common.models import WebResource, WebLink, Publication
from structure.models import StructureType, StructureStabilizingAgent, PdbData, Rotamer
from structure.functions import get_pdb_ids, create_structure_rotamer, fetch_signprot_data, build_signprot_struct
from common.tools import test_model_updates

import re
from Bio import pairwise2
from collections import OrderedDict
import logging
import shlex, subprocess
from io import StringIO
from Bio.PDB import PDBParser, PPBuilder, PDBIO, Polypeptide
from Bio import pairwise2
import pprint
import json
import yaml
import urllib
import django.apps
import traceback
import sys, os
import datetime
import csv


class Command(BaseBuild):

    #Setting the variables for the test tracking of the model upadates
    tracker = {}
    all_models = django.apps.apps.get_models()[6:]
    test_model_updates(all_models, tracker, initialize=True)
    gain_csv = os.sep.join([settings.DATA_DIR, "structure_data", "annotation", "GAIN-3.csv"])
    annotation_dict = {}
    # original
    # segment_labels = ['A.H1','A.h1h2','A.H2','A.h2h3','A.H3','A.h3h4','A.H4','A.h4h5','A.H5','A.h5h6','A.H6','A.h6h7','A.H7','A.h7h8','A.H8','A.h8s1','B.S1','B.s1s2','B.S2','B.s2s3','B.S3','B.s3s4','B.S4','B.s4s5','B.S5','B.s5s6','B.S6','B.s6s7','B.S7','B.s7s8','B.S8','B.s8s9','B.S9','B.s9s10','B.S10','B.s10s11','B.S11','B.s11s12','B.S12','B.s12gps','B.GPS','B.gpss13','B.S13','B.s13tm1']
    # updated
    segment_labels = ['A.H1','A.h1h2','A.H2','A.h2h3','A.H3','A.h3h4','A.H4','A.h4h5','A.H5','A.h5h6','A.H6','A.h6s1','B.S1','B.s1s2','B.S2','B.s2s3','B.S3','B.s3s4','B.S4','B.s4s5','B.S5','B.s5s6','B.S6','B.s6s7','B.S7','B.s7s8','B.S8','B.s8s9','B.S9','B.s9s10','B.S10','B.s10s11','B.S11','B.s11s12','B.S12','B.s12s13','B.S13','B.s13gps','B.GPS','B.gpss14','B.S14','B.s14tm1']
    segments = []
    for s in segment_labels:
        segments.append(ProteinSegment.objects.get(slug=s, proteinfamily='GPCR'))

    def add_arguments(self, parser):
        parser.add_argument("--debug", default=False, action="store_true", help="Debug mode")

    def handle(self, *args, **options):
        self.parse_csv()
        gn_scheme_parent, _ = ResidueNumberingScheme.objects.get_or_create(slug='gain', short_name='GAIN', name='GAIN', parent=None)
        gn_scheme, _ = ResidueNumberingScheme.objects.get_or_create(slug='gpcrdbgain', short_name='GPCRdb(GAIN)', name='GPCRdb(GAIN)', parent=gn_scheme_parent)
        for protein, ends in self.annotation_dict.items():
            residues = Residue.objects.filter(protein_conformation__protein=protein, protein_segment__slug='N-term')
            tm1_start_pos = Residue.objects.filter(protein_conformation__protein=protein, protein_segment__slug='TM1')[0].sequence_number
            print(protein)
            print(ends)
            print(self.segments)
            first_seg = None
            annotated_segments = []
            for s, pos in ends.items():
                if pos!='-':
                    if s.endswith('b'):
                        annotated_segments.append(s)
                    elif '+' in s:
                        annotated_segments.append(s[:-2])
            first_seg = annotated_segments[0]
            for seg in self.segments:
                ### Helices and sheets
                if seg.category!='loop':
                    print(seg)
                    ### Special treatment for GPS segment
                    if seg.slug=='B.GPS':
                        for gpsgn in ['B.GPS-2','B.GPS-1','B.GPS+1']:
                            gn, _ = ResidueGenericNumber.objects.get_or_create(label=gpsgn, protein_segment=seg, scheme=gn_scheme)
                            seqnum = int(ends[gpsgn])
                            res = residues.get(sequence_number=seqnum)
                            res.display_generic_number = gn
                            res.generic_number = gn
                            res.protein_segment = seg
                            res.save()
                            print(res.sequence_number, res.generic_number, res, res.protein_segment)
                        ### residues between GPS-1 and GPS+1
                        if int(ends['B.GPS+1'])-int(ends['B.GPS-1'])>1:
                            for i in range(int(ends['B.GPS-1'])+1, int(ends['B.GPS+1'])):
                                res = residues.get(sequence_number=i)
                                res.protein_segment = seg
                                res.save()
                                print(res.sequence_number, res.generic_number, res, res.protein_segment)
                    ### Other segments
                    elif seg.slug+'b' in ends:
                        if ends[seg.slug+'b']=='-' and self.segment_labels.index(seg.slug)<self.segment_labels.index(first_seg[:-1]):
                            continue
                        ### Missing segment
                        elif ends[seg.slug+'b']=='-':
                            print('Missing segment, assigning', self.segments[self.segments.index(seg)-1])
                            continue
                        start = int(ends[seg.slug+'b'])
                        x50 = int(ends[seg.slug+'x'])
                        end = int(ends[seg.slug+'e'])
                    
                        # if ':
                        #     start = int(ends['B.GPS-2'])
                        #     x50 = int(ends['B.GPS-1'])
                        #     end = int(ends['B.GPS+1'])-1
                        # elif seg.slug=='B.GPS+':
                        #     start = int(ends['B.GPS+1'])
                        #     x50 = start
                        #     end = start
                        start_lab = x50-start
                        ### Start
                        for num in range(start, x50+1):
                            gn, _ = ResidueGenericNumber.objects.get_or_create(label=seg.slug+'.'+str(50-start_lab), protein_segment=seg, scheme=gn_scheme)
                            print(num, seg.slug+'.'+str(50-start_lab), gn, _)
                            res = residues.get(sequence_number=num)
                            res.display_generic_number = gn
                            res.generic_number = gn
                            res.protein_segment = seg
                            print(res.sequence_number, res.generic_number, res, res.protein_segment)
                            res.save()
                            start_lab-=1
                        ### End
                        end_lab = 51
                        for num in range(x50+1,end+1):
                            gn, _ = ResidueGenericNumber.objects.get_or_create(label=seg.slug+'.'+str(end_lab), protein_segment=seg, scheme=gn_scheme)
                            # print(num, seg.slug+'.'+str(end_lab), gn, _)
                            try:
                                res = residues.get(sequence_number=num)
                            except Residue.DoesNotExist:
                                raise AssertionError('ERROR: {} at {} {} {} is annotated incorrectly'.format(protein, seg, gn, num))
                            res.display_generic_number = gn
                            res.generic_number = gn
                            res.protein_segment = seg
                            print(res.sequence_number, res.generic_number, res, res.protein_segment)
                            res.save()
                            end_lab+=1
                    # break
                ### Loops
                else:
                    print(seg)
                    if seg.slug=='B.gps':
                        end_tag = ''
                        b_tag = ''
                    elif seg.slug=='B.s13gps':
                        end_tag = 'e'
                        b_tag = ''
                    elif seg.slug=='B.gpss14':
                        end_tag = ''
                        b_tag = 'b'
                    else:
                        end_tag = 'e'
                        b_tag = 'b'
                    if self.segments[self.segments.index(seg)-1].slug+end_tag=='B.GPS':
                        prev_end = 'B.GPS+1'
                    # elif self.segments[self.segments.index(seg)-1].slug+end_tag=='B.GPS+':
                    #     prev_end = 'B.GPS+1'
                    else:
                        prev_end = self.segments[self.segments.index(seg)-1].slug+end_tag
                    if ends[prev_end]=='-':
                        continue
                    loop_index = self.segments.index(seg)
                    if self.segments[loop_index-1].slug+end_tag=='B.GPS':
                        l_prev_end = 'B.GPS+1'
                    # elif self.segments[loop_index-1].slug+end_tag=='B.GPS+':
                    #     l_prev_end = 'B.GPS+1'
                    else:
                        l_prev_end = self.segments[loop_index-1].slug+end_tag
                    loop_start = int(ends[l_prev_end])+1



                    ### s14tm1 domain end
                    if loop_index==len(self.segments)-1:
                        # print(seg, loop_start, tm1_start_pos)
                        for num in range(loop_start, tm1_start_pos):
                            res = residues.get(sequence_number=num)
                            res.protein_segment = seg
                            print(res.sequence_number, res.generic_number, res, res.protein_segment)
                            res.save()
                        continue
                    next_seg_offset = 1
                    ### In between loops
                    while self.segments[loop_index+next_seg_offset].slug+b_tag not in annotated_segments:
                        self.segments[loop_index+next_seg_offset].slug+b_tag
                        next_seg_offset+=1
                    if self.segments[loop_index+next_seg_offset].slug+b_tag=='B.GPS':
                        next_start = 'B.GPS-2'
                    # elif self.segments[loop_index+next_seg_offset].slug+b_tag=='B.GPS+':
                    #     next_start = 'B.GPS+1'
                    else:
                        next_start = self.segments[loop_index+next_seg_offset].slug+b_tag
                    # print(seg, loop_start, loop_index, int(ends[next_start])-1, next_seg_offset)
                    for num in range(loop_start, int(ends[next_start])):
                        res = residues.get(sequence_number=num)
                        res.protein_segment = seg
                        print(res.sequence_number, res.generic_number, res, res.protein_segment)
                        res.save()
            # break

        

    def parse_csv(self):
        proteins = Protein.objects.filter(family__parent__slug__startswith='003_', accession__isnull=False)
        with open(self.gain_csv, newline='') as csvfile:
            gaincsvreader = csv.reader(csvfile, delimiter=';', quotechar='|')
            for i, row in enumerate(gaincsvreader):
                if i==0:
                    header = row
                    continue
                ends = {}
                for j in range(2, len(header)):
                    if row[j]=='':
                        pos = '-'
                    else:
                        pos = row[j]
                    this_header = header[j].replace('.start','b').replace('.end','e').replace('.anchor','x')
                    ends[this_header] = pos
                self.annotation_dict[proteins.get(accession=row[0])] = ends
                # break