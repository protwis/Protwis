from django.shortcuts import render
from django.conf import settings
from django.views.generic import TemplateView, View
from django.http import HttpResponse, JsonResponse, HttpResponseRedirect
from django.db.models import Count, Q, Prefetch
from django import forms
from django.core.cache import cache
from django.views.decorators.cache import cache_page

from common.phylogenetic_tree import PhylogeneticTreeGenerator
from protein.models import Gene, ProteinSegment
from structure.models import Structure, StructureModel, StructureModelStatsRotamer, StructureModelSeqSim, StructureRefinedStatsRotamer, StructureRefinedSeqSim
from structure.functions import CASelector, SelectionParser, GenericNumbersSelector, SubstructureSelector, check_gn
from structure.assign_generic_numbers_gpcr import GenericNumbering
from structure.structural_superposition import ProteinSuperpose,FragmentSuperpose
from structure.forms import *
from interaction.models import ResidueFragmentInteraction,StructureLigandInteraction
from protein.models import Protein, ProteinFamily
from construct.models import Construct
from construct.functions import convert_ordered_to_disordered_annotation,add_construct
from common.views import AbsSegmentSelection,AbsReferenceSelection
from common.selection import Selection, SelectionItem
from common.extensions import MultiFileField
from common.models import ReleaseNotes

Alignment = getattr(__import__('common.alignment_' + settings.SITE_NAME, fromlist=['Alignment']), 'Alignment')

import inspect
import os
import time
import zipfile
import math
import json
import ast
from copy import deepcopy
from io import StringIO, BytesIO
from collections import OrderedDict
from Bio.PDB import PDBIO, PDBParser
from operator import itemgetter
import xlsxwriter

from email.mime.text import MIMEText
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart
from smtplib import SMTP
import smtplib
import sys

class_tree = {'001':'A','002':'B1','003':'B2','004':'C','005':'F','006':'T'}

class StructureBrowser(TemplateView):
    """
    Fetching Structure data for browser
    """

    template_name = "structure_browser.html"

    def get_context_data (self, **kwargs):

        context = super(StructureBrowser, self).get_context_data(**kwargs)
        try:
            context['structures'] = Structure.objects.filter(refined=False).select_related(
                "pdb_code__web_resource",
                "protein_conformation__protein__species",
                "protein_conformation__protein__source",
                "protein_conformation__protein__family__parent__parent__parent",
                "publication__web_link__web_resource").prefetch_related(
                "stabilizing_agents",
                "protein_conformation__protein__parent__endogenous_ligands__properities__ligand_type",
                Prefetch("ligands", queryset=StructureLigandInteraction.objects.filter(
                annotated=True).prefetch_related('ligand__properities__ligand_type', 'ligand_role')))
            context['refined'] = [i.pdb_code.index[:4] for i in Structure.objects.filter(refined=True)]
        except Structure.DoesNotExist as e:
            pass

        return context


class ServeHomologyModels(TemplateView):

    template_name = "homology_models.html"
    def get_context_data(self, **kwargs):
        context = super(ServeHomologyModels, self).get_context_data(**kwargs)
        try:
            context['structure_model'] = StructureModel.objects.all().select_related(
                "protein__family",
                "state",
                "protein__family__parent__parent__parent",
                "protein__species",
                "main_template__protein_conformation__protein__parent__family",
                "main_template__pdb_code")
        except StructureModel.DoesNotExist as e:
            pass

        return context


def HomologyModelDetails(request, modelname, state):
    """
    Show homology models details
    """
    modelname = modelname
    color_palette = ["orange","cyan","yellow","lime","fuchsia","green","teal","olive","thistle","grey","chocolate","blue","red","pink","maroon",]
    
    if state=='refined':
        model = Structure.objects.get(pdb_code__index=modelname+'_refined')
        model_main_template = Structure.objects.get(pdb_code__index=modelname)
        rotamers = StructureRefinedStatsRotamer.objects.filter(structure=model).prefetch_related("structure", "residue", "backbone_template", "rotamer_template").order_by('residue__sequence_number')
        main_template_seqsim = StructureRefinedSeqSim.objects.get(structure=model, template=model_main_template).similarity
    else:
        model = StructureModel.objects.get(protein__entry_name=modelname, state__slug=state)
        model_main_template = model.main_template
        rotamers = StructureModelStatsRotamer.objects.filter(homology_model=model).prefetch_related("homology_model", "residue", "backbone_template", "rotamer_template").order_by('residue__sequence_number')
        main_template_seqsim = StructureModelSeqSim.objects.get(homology_model=model, template=model_main_template).similarity

    backbone_templates, rotamer_templates = [],[]
    segments, segments_formatted, segments_out = {},{},{}
    bb_temps, r_temps = OrderedDict(), OrderedDict()
    bb_main, bb_alt, bb_none = 0,0,0
    sc_main, sc_alt, sc_none = 0,0,0

    for r in rotamers:
        if r.backbone_template not in backbone_templates and r.backbone_template!=None:
            backbone_templates.append(r.backbone_template)
            if r.backbone_template.protein_conformation.protein.parent not in bb_temps:
                bb_temps[r.backbone_template.protein_conformation.protein.parent] = [r.backbone_template]
            else:
                bb_temps[r.backbone_template.protein_conformation.protein.parent].append(r.backbone_template)
        if r.rotamer_template not in rotamer_templates and r.rotamer_template!=None:
            rotamer_templates.append(r.rotamer_template)
            if r.rotamer_template.protein_conformation.protein.parent not in r_temps:
                r_temps[r.rotamer_template.protein_conformation.protein.parent] = [r.rotamer_template]
            else:
                r_temps[r.rotamer_template.protein_conformation.protein.parent].append(r.rotamer_template)
        if r.backbone_template not in segments:
            segments[r.backbone_template] = [r.residue.sequence_number]
        else:
            segments[r.backbone_template].append(r.residue.sequence_number)
        if r.backbone_template==model_main_template:
            bb_main+=1
        elif r.backbone_template!=None:
            bb_alt+=1
        elif r.backbone_template==None:
            bb_none+=1
        if r.rotamer_template==model_main_template:
            sc_main+=1
        elif r.rotamer_template!=None:
            sc_alt+=1
        elif r.rotamer_template==None:
            sc_none+=1
    for s, nums in segments.items():
        for i, num in enumerate(nums):
            if i==0:
                segments_formatted[s] = [[num]]
            elif nums[i-1]!=num-1:
                if segments_formatted[s][-1][0]==nums[i-1]:
                    segments_formatted[s][-1] = '{}-{}'.format(segments_formatted[s][-1][0], nums[i-1])
                else:
                    segments_formatted[s][-1] = '{}-{}'.format(segments_formatted[s][-1][0], nums[i-1])
                segments_formatted[s].append([num])
                if i+1==len(segments[s]):
                    segments_formatted[s][-1] = '{}-{}'.format(segments_formatted[s][-1][0], segments_formatted[s][-1][0])
            elif i+1==len(segments[s]):
                segments_formatted[s][-1] = '{}-{}'.format(segments_formatted[s][-1][0], nums[i-1]+1)
        if len(nums)==1:
            segments_formatted[s] = ['{}-{}'.format(segments_formatted[s][0][0], segments_formatted[s][0][0])]

    colors = OrderedDict([(model_main_template,"darkorchid"), (None,"white")])
    i = 0
    for s, nums in segments_formatted.items():
        if len(nums)>1:
            text = ''
            for n in nums:
                text+='{} or '.format(n)
            segments_formatted[s] = text[:-4]
        else:
            segments_formatted[s] = segments_formatted[s][0]
        if s==model_main_template:
            pass
        elif s==None:
            segments_out["white"] = segments_formatted[s]
        else:
            segments_out[color_palette[i]] = segments_formatted[s]
            colors[s] = color_palette[i]
        i+=1
    template_list = []
    for b, temps in bb_temps.items():
        for i, t in enumerate(temps):
            t.color = colors[t]
            bb_temps[b][i] = t
            template_list.append(t.pdb_code.index)

    return render(request,'homology_models_details.html',{'model': model, 'modelname': modelname, 'rotamers': rotamers, 'backbone_templates': bb_temps, 'backbone_templates_number': len(backbone_templates),
                                                          'rotamer_templates': r_temps, 'rotamer_templates_number': len(rotamer_templates), 'color_residues': segments_out, 'bb_main': round(bb_main/len(rotamers)*100, 1),
                                                          'bb_alt': round(bb_alt/len(rotamers)*100, 1), 'bb_none': round(bb_none/len(rotamers)*100, 1), 'sc_main': round(sc_main/len(rotamers)*100, 1), 'sc_alt': round(sc_alt/len(rotamers)*100, 1),
                                                          'sc_none': round(sc_none/len(rotamers)*100, 1), 'main_template_seqsim': main_template_seqsim, 'template_list': template_list, 'model_main_template': model_main_template,
                                                          'state': state})

def ServeHomModDiagram(request, modelname, state):
    if state=='refined':
        model=Structure.objects.filter(pdb_code__index=modelname+'_refined')
    else:
        model=StructureModel.objects.filter(protein__entry_name=modelname, state__slug=state)
    if model.exists():
        model=model.get()
    else:
         quit() #quit!

    if model.pdb is None:
        quit()

    response = HttpResponse(model.pdb, content_type='text/plain')
    return response

def StructureDetails(request, pdbname):
    """
    Show structure details
    """
    pdbname = pdbname
    structures = ResidueFragmentInteraction.objects.values('structure_ligand_pair__ligand__name','structure_ligand_pair__pdb_reference','structure_ligand_pair__annotated').filter(structure_ligand_pair__structure__pdb_code__index=pdbname, structure_ligand_pair__annotated=True).annotate(numRes = Count('pk', distinct = True)).order_by('-numRes')
    resn_list = ''

    main_ligand = 'None'
    for structure in structures:
        if structure['structure_ligand_pair__annotated']:
            resn_list += ",\""+structure['structure_ligand_pair__pdb_reference']+"\""
            main_ligand = structure['structure_ligand_pair__pdb_reference']

    crystal = Structure.objects.get(pdb_code__index=pdbname)
    p = Protein.objects.get(protein=crystal.protein_conformation.protein)
    residues = ResidueFragmentInteraction.objects.filter(structure_ligand_pair__structure__pdb_code__index=pdbname, structure_ligand_pair__annotated=True).order_by('rotamer__residue__sequence_number')
    try:
        refined = Structure.objects.get(pdb_code__index=pdbname+'_refined')
    except:
        refined = False
    return render(request,'structure_details.html',{'pdbname': pdbname, 'structures': structures, 'crystal': crystal, 'protein':p, 'residues':residues, 'annotated_resn': resn_list, 'main_ligand': main_ligand, 'refined': refined})

def ServePdbDiagram(request, pdbname):
    structure=Structure.objects.filter(pdb_code__index=pdbname)
    if structure.exists():
        structure=structure.get()
    else:
         quit() #quit!

    if structure.pdb_data is None:
        quit()

    response = HttpResponse(structure.pdb_data.pdb, content_type='text/plain')
    return response


def ServePdbLigandDiagram(request,pdbname,ligand):
    pair = StructureLigandInteraction.objects.filter(structure__pdb_code__index=pdbname).filter(Q(ligand__properities__inchikey=ligand) | Q(ligand__name=ligand)).exclude(pdb_file__isnull=True).get()
    response = HttpResponse(pair.pdb_file.pdb, content_type='text/plain')
    return response

class StructureStatistics(TemplateView):
    """
    So not ready that EA wanted to publish it.
    """
    template_name = 'structure_statistics.html'

    def get_context_data (self, **kwargs):
        context = super().get_context_data(**kwargs)

        families = ProteinFamily.objects.all()
        lookup = {}
        for f in families:
            lookup[f.slug] = f.name

        all_structs = Structure.objects.all().prefetch_related('protein_conformation__protein__family').exclude(refined=True)
        all_complexes = all_structs.exclude(ligands=None)
        #FIXME G protein list is hard-coded for now. Table structure needs to be expanded for fully automatic approach.
        all_gprots = all_structs.filter(stabilizing_agents__slug='gs')
        all_active = all_structs.filter(protein_conformation__state__slug = 'active')        

        years = self.get_years_range(list(set([x.publication_date.year for x in all_structs])))
        
        unique_structs = Structure.objects.order_by('protein_conformation__protein__family__name', 'state',
            'publication_date', 'resolution').distinct('protein_conformation__protein__family__name').prefetch_related('protein_conformation__protein__family')
        unique_complexes = all_complexes.distinct('ligands', 'protein_conformation__protein__family__name')
        #FIXME G protein list is hard-coded for now. Table structure needs to be expanded for fully automatic approach.
        unique_gprots = unique_structs.filter(stabilizing_agents__slug='gs')
        unique_active = unique_structs.filter(protein_conformation__state__slug = 'active')
        
        #Stats
        struct_count = Structure.objects.all().annotate(Count('id'))
        struct_lig_count = Structure.objects.exclude(ligands=None)
                
        context['all_structures'] = len(all_structs)
        context['all_structures_by_class'] = self.count_by_class(all_structs, lookup)
        context['all_complexes'] = len(all_complexes)
        context['all_complexes_by_class'] = self.count_by_class(all_complexes, lookup)
        context['all_gprots'] = len(all_gprots)
        context['all_gprots_by_class'] = self.count_by_class(all_gprots, lookup)
        context['all_active'] = len(all_active)
        context['all_active_by_class'] = self.count_by_class(all_active, lookup)

        context['unique_structures'] = len(unique_structs)
        context['unique_structures_by_class'] = self.count_by_class(unique_structs, lookup)
        context['unique_complexes'] = len(unique_complexes)
        context['unique_complexes_by_class'] = self.count_by_class(unique_complexes, lookup)
        context['unique_gprots'] = len(unique_gprots)
        context['unique_gprots_by_class'] = self.count_by_class(unique_gprots, lookup)
        context['unique_active'] = len(unique_active)
        context['unique_active_by_class'] = self.count_by_class(unique_active, lookup)
        context['release_notes'] = ReleaseNotes.objects.all()[0]

        context['chartdata'] = self.get_per_family_cumulative_data_series(years, unique_structs, lookup)
        context['chartdata_y'] = self.get_per_family_data_series(years, unique_structs, lookup)
        context['chartdata_all'] = self.get_per_family_cumulative_data_series(years, all_structs, lookup)
        context['chartdata_reso'] = self.get_resolution_coverage_data_series(all_structs)
        #context['coverage'] = self.get_diagram_coverage()
        #{
        #    'depth': 3,
        #    'anchor': '#crystals'}
        tree = PhylogeneticTreeGenerator()
        class_a_data = tree.get_tree_data(ProteinFamily.objects.get(name='Class A (Rhodopsin)'))
        context['class_a_options'] = deepcopy(tree.d3_options)
        context['class_a_options']['anchor'] = 'class_a'
        context['class_a_options']['leaf_offset'] = 50
        context['class_a_options']['label_free'] = []
        context['class_a'] = json.dumps(class_a_data.get_nodes_dict('crystals'))
        class_b1_data = tree.get_tree_data(ProteinFamily.objects.get(name__startswith='Class B1 (Secretin)'))
        context['class_b1_options'] = deepcopy(tree.d3_options)
        context['class_b1_options']['anchor'] = 'class_b1'
        context['class_b1_options']['branch_trunc'] = 60
        context['class_b1_options']['label_free'] = [1,]
        context['class_b1'] = json.dumps(class_b1_data.get_nodes_dict('crystals'))
        class_b2_data = tree.get_tree_data(ProteinFamily.objects.get(name__startswith='Class B2 (Adhesion)'))
        context['class_b2_options'] = deepcopy(tree.d3_options)
        context['class_b2_options']['anchor'] = 'class_b2'
        context['class_b2_options']['label_free'] = [1,]
        context['class_b2'] = json.dumps(class_b2_data.get_nodes_dict('crystals'))
        class_c_data = tree.get_tree_data(ProteinFamily.objects.get(name__startswith='Class C (Glutamate)'))
        context['class_c_options'] = deepcopy(tree.d3_options)
        context['class_c_options']['anchor'] = 'class_c'
        context['class_c_options']['branch_trunc'] = 50
        context['class_c_options']['label_free'] = [1,]
        context['class_c'] = json.dumps(class_c_data.get_nodes_dict('crystals'))
        class_f_data = tree.get_tree_data(ProteinFamily.objects.get(name__startswith='Class F (Frizzled)'))
        context['class_f_options'] = deepcopy(tree.d3_options)
        context['class_f_options']['anchor'] = 'class_f'
        context['class_f_options']['label_free'] = [1,]
        #json.dump(class_f_data.get_nodes_dict('crystalized'), open('tree_test.json', 'w'), indent=4)
        context['class_f'] = json.dumps(class_f_data.get_nodes_dict('crystals'))
        class_t2_data = tree.get_tree_data(ProteinFamily.objects.get(name='Taste 2'))
        context['class_t2_options'] = deepcopy(tree.d3_options)
        context['class_t2_options']['anchor'] = 'class_t2'
        context['class_t2_options']['label_free'] = [1,]
        context['class_t2'] = json.dumps(class_t2_data.get_nodes_dict('crystals'))

        return context

    def get_families_dict(self, queryset, lookup):

        families = []
        for s in queryset:
            fid = s.protein_conformation.protein.family.slug.split("_")
            fname = lookup[fid[0]+"_"+fid[1]]
            cname = lookup[fid[0]]
            if fname not in families:
                families.append(fname)
        return families

    def count_by_class(self, queryset, lookup):

        #Ugly walkaround
        classes = [lookup[x] for x in lookup.keys() if x in ['001', '002', '003', '004', '005', '006']]
        records = []
        for s in queryset:
            fid = s.protein_conformation.protein.family.slug.split("_")
            cname = lookup[fid[0]]
            records.append(cname)

        tmp = OrderedDict()
        for x in sorted(classes):
            tmp[x] = records.count(x)

        return tmp

    def get_years_range(self, years_list):

        min_y = min(years_list)
        max_y = max(years_list)
        return range(min_y, max_y+1)


    def get_per_family_data_series(self, years, structures, lookup):
        """
        Prepare data for multiBarGraph of unique crystallized receptors. Returns data series for django-nvd3 wrapper.
        """
        families = self.get_families_dict(structures, lookup)
        series = []
        data = {}
        for year in years:
            for family in families:
                if family not in data.keys():
                    data[family] = []
                count = 0
                for structure in structures:
                    fid = structure.protein_conformation.protein.family.slug.split("_")
                    # if structure.protein_conformation.protein.get_protein_family() == family and structure.publication_date.year == year:
                    if lookup[fid[0]+"_"+fid[1]] == family and structure.publication_date.year == year:
                        count += 1
                data[family].append(count)
        for family in families:
            series.append({"values":
                [{
                    'x': years[i],
                    'y': j
                    } for i, j in enumerate(data[family])],
                "key": family,
                "yAxis": "1"})
        return json.dumps(series)


    def get_per_family_cumulative_data_series(self, years, structures, lookup):
        """
        Prepare data for multiBarGraph of unique crystallized receptors. Returns data series for django-nvd3 wrapper.
        """
        families = self.get_families_dict(structures, lookup)
        series = []
        data = {}
        for year in years:
            for family in families:
                if family not in data.keys():
                    data[family] = []
                count = 0
                for structure in structures:
                    fid = structure.protein_conformation.protein.family.slug.split("_")
                    # if structure.protein_conformation.protein.get_protein_family() == family and structure.publication_date.year == year:
                    if lookup[fid[0]+"_"+fid[1]] == family and structure.publication_date.year == year:
                        count += 1
                if len(data[family]) > 0:
                    data[family].append(count + data[family][-1])
                else:
                    data[family].append(count)
        for family in families:
            series.append({"values":
                [{
                    'x': years[i],
                    'y': j
                    } for i, j in enumerate(data[family])],
                "key": family,
                "yAxis": "1"})
        return json.dumps(series)


    def get_resolution_coverage_data_series(self, structures):
        """
        Prepare data for multiBarGraph of resolution coverage of available crystal structures.
        """
        #Resolutions boundaries
        reso_min = float(min([round(x.resolution, 1) for x in structures]))
        reso_max = float(max([round(x.resolution, 1) for x in structures]))
        step = (reso_max - reso_min)/10

        brackets = [reso_min + step*x for x in range(10)] + [reso_max]

        reso_count = []
        bracket_labels = []
        for idx, bracket in enumerate(brackets):
            if idx == 0:
                reso_count.append(len([x for x in structures if x.resolution <= bracket]))
                bracket_labels.append('< {:.1f}'.format(bracket))
            else:
                reso_count.append(len([x for x in structures if bracket-step < x.resolution <= bracket]))
                bracket_labels.append('{:.1f}-{:.1f}'.format(brackets[idx-1],bracket))

        return json.dumps([{"values": [{
                    'x': bracket_labels[i],
                    'y': j
                    } for i, j in enumerate(reso_count)],
                "key": 'Resolution coverage',
                "yAxis": "1"}])

    def get_diagram_coverage(self):
        """
        Prepare data for coverage diagram.
        """

        families = ProteinFamily.objects.all()
        lookup = {}
        for f in families:
            lookup[f.slug] = f.name.replace("receptors","")

        class_proteins = Protein.objects.filter(family__slug__startswith="00", source__name='SWISSPROT').prefetch_related('family').order_by('family__slug')

        coverage = OrderedDict()

        temp = OrderedDict([
                            ('name',''),
                            ('interactions', 0),
                            ('receptor_i', 0) ,
                            ('mutations' , 0),
                            ('receptor_m', 0),
                            ('mutations_an' , 0),
                            ('receptor_m_an', 0),
                            ('receptor_t',0),
                            ('children', OrderedDict()) ,
                            ('fraction_i',0),
                            ('fraction_m',0),
                            ('fraction_m_an',0)
                            ])

        for p in class_proteins:
            fid = p.family.slug.split("_")
            if fid[0] not in coverage:
                coverage[fid[0]] = deepcopy(temp)
                coverage[fid[0]]['name'] = lookup[fid[0]]
            if fid[1] not in coverage[fid[0]]['children']:
                coverage[fid[0]]['children'][fid[1]] = deepcopy(temp)
                coverage[fid[0]]['children'][fid[1]]['name'] = lookup[fid[0]+"_"+fid[1]]
            if fid[2] not in coverage[fid[0]]['children'][fid[1]]['children']:
                coverage[fid[0]]['children'][fid[1]]['children'][fid[2]] = deepcopy(temp)
                coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['name'] = lookup[fid[0]+"_"+fid[1]+"_"+fid[2]][:28]
            if fid[3] not in coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children']:
                coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]] = deepcopy(temp)
                coverage[fid[0]]['receptor_t'] += 1
                coverage[fid[0]]['children'][fid[1]]['receptor_t'] += 1
                coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['receptor_t'] += 1
                coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['name'] = p.entry_name.split("_")[0] #[:10]
                coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['receptor_t'] = 1


        class_interactions = ResidueFragmentInteraction.objects.filter(structure_ligand_pair__annotated=True).prefetch_related(
            'rotamer__residue__display_generic_number','interaction_type',
            'structure_ligand_pair__structure__protein_conformation__protein__parent__family',
            'structure_ligand_pair__ligand__properities',
            )


        score_copy = {'score': {'a':0,'i':0,'i_weight':0,'m':0,'m_weight':0,'s':0,'s_weight':0} , 'interaction' : {},'mutation': {}}

        # Replace above as fractions etc is not required and it was missing xtals that didnt have interactions.
        unique_structs = list(Structure.objects.order_by('protein_conformation__protein__parent', 'state',
            'publication_date', 'resolution').distinct('protein_conformation__protein__parent').prefetch_related('protein_conformation__protein__family'))

        for s in unique_structs:
            fid = s.protein_conformation.protein.family.slug.split("_")
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['receptor_i'] = 1
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['interactions'] += 1

        CSS_COLOR_NAMES = ["SteelBlue","SlateBlue","LightCoral","Orange","LightGreen","LightGray","PeachPuff","PaleGoldenRod"]

        tree = OrderedDict({'name':'GPCRs','children':[]})
        i = 0
        n = 0
        for c,c_v in coverage.items():
            c_v['name'] = c_v['name'].split("(")[0]
            if c_v['name'].strip() == 'Other GPCRs':
                continue
            children = []
            for lt,lt_v in c_v['children'].items():
                if lt_v['name'].strip() == 'Orphan' and c_v['name'].strip()=="Class A":
                    continue
                children_rf = []
                for rf,rf_v in lt_v['children'].items():
                    rf_v['name'] = rf_v['name'].split("<")[0]
                    if rf_v['name'].strip() == 'Taste 2':
                        continue
                    children_r = []
                    for r,r_v in rf_v['children'].items():
                        r_v['color'] = CSS_COLOR_NAMES[i]
                        r_v['sort'] = n
                        children_r.append(r_v)
                        n += 1
                    rf_v['children'] = children_r
                    rf_v['sort'] = n
                    rf_v['color'] = CSS_COLOR_NAMES[i]
                    children_rf.append(rf_v)
                lt_v['children'] = children_rf
                lt_v['sort'] = n
                lt_v['color'] = CSS_COLOR_NAMES[i]
                children.append(lt_v)
            c_v['children'] = children
            c_v['sort'] = n
            c_v['color'] = CSS_COLOR_NAMES[i]
            tree['children'].append(c_v)
            i += 1

        return json.dumps(tree)


    def get_diagram_crystals(self):
        """
        Prepare data for coverage diagram.
        """

        crystal_proteins = [x.protein_conformation.protein.parent for x in Structure.objects.order_by('protein_conformation__protein__parent', 'state',
            'publication_date', 'resolution').distinct('protein_conformation__protein__parent').prefetch_related('protein_conformation__protein__parent__family')]

        families = []
        for cryst_prot in crystal_proteins:
            families.append(cryst_prot.family)
            tmp = cryst_prot.family
            while tmp.parent is not None:
                tmp = tmp.parent
                families.append(tmp)
        lookup = {}
        for f in families:
            lookup[f.slug] = f.name.replace("receptors","")

        coverage = OrderedDict()
        temp = OrderedDict([
                            ('name',''),
                            ('interactions', 0),
                            ('receptor_i', 0) ,
                            ('mutations' , 0),
                            ('receptor_m', 0),
                            ('mutations_an' , 0),
                            ('receptor_m_an', 0),
                            ('receptor_t',0),
                            ('children', OrderedDict()) ,
                            ('fraction_i',0),
                            ('fraction_m',0),
                            ('fraction_m_an',0)
                            ])

        for p in crystal_proteins:
            fid = p.family.slug.split("_")
            if fid[0] not in coverage:
                coverage[fid[0]] = deepcopy(temp)
                coverage[fid[0]]['name'] = lookup[fid[0]]
            if fid[1] not in coverage[fid[0]]['children']:
                coverage[fid[0]]['children'][fid[1]] = deepcopy(temp)
                coverage[fid[0]]['children'][fid[1]]['name'] = lookup[fid[0]+"_"+fid[1]]
            if fid[2] not in coverage[fid[0]]['children'][fid[1]]['children']:
                coverage[fid[0]]['children'][fid[1]]['children'][fid[2]] = deepcopy(temp)
                coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['name'] = lookup[fid[0]+"_"+fid[1]+"_"+fid[2]][:28]
            if fid[3] not in coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children']:
                coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]] = deepcopy(temp)
                coverage[fid[0]]['receptor_t'] += 1
                coverage[fid[0]]['children'][fid[1]]['receptor_t'] += 1
                coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['receptor_t'] += 1
                coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['name'] = p.entry_name.split("_")[0] #[:10]
                coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['receptor_t'] = 1


        class_interactions = ResidueFragmentInteraction.objects.filter(structure_ligand_pair__annotated=True).prefetch_related(
            'rotamer__residue__display_generic_number','interaction_type',
            'structure_ligand_pair__structure__protein_conformation__protein__parent__family',
            'structure_ligand_pair__ligand__properities',
            )

        score_copy = {'score': {'a':0,'i':0,'i_weight':0,'m':0,'m_weight':0,'s':0,'s_weight':0} , 'interaction' : {},'mutation': {}}

        # Replace above as fractions etc is not required and it was missing xtals that didnt have interactions.
        unique_structs = list(Structure.objects.order_by('protein_conformation__protein__family__name', 'state',
            'publication_date', 'resolution').distinct('protein_conformation__protein__family__name').prefetch_related('protein_conformation__protein__family'))

        for s in unique_structs:
            fid = s.protein_conformation.protein.family.slug.split("_")
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['receptor_i'] = 1
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['interactions'] += 1

        CSS_COLOR_NAMES = ["SteelBlue","SlateBlue","LightCoral","Orange","LightGreen","LightGray","PeachPuff","PaleGoldenRod"]

        tree = OrderedDict({'name':'GPCRs','children':[]})
        i = 0
        n = 0
        for c,c_v in coverage.items():
            c_v['name'] = c_v['name'].split("(")[0]
            if c_v['name'].strip() == 'Other GPCRs':
                continue
            children = []
            for lt,lt_v in c_v['children'].items():
                if lt_v['name'].strip() == 'Orphan' and c_v['name'].strip()=="Class A":
                    continue
                children_rf = []
                for rf,rf_v in lt_v['children'].items():
                    rf_v['name'] = rf_v['name'].split("<")[0]
                    if rf_v['name'].strip() == 'Taste 2':
                        continue
                    children_r = []
                    for r,r_v in rf_v['children'].items():
                        r_v['color'] = CSS_COLOR_NAMES[i]
                        r_v['sort'] = n
                        children_r.append(r_v)
                        n += 1
                    rf_v['children'] = children_r
                    rf_v['sort'] = n
                    rf_v['color'] = CSS_COLOR_NAMES[i]
                    children_rf.append(rf_v)
                lt_v['children'] = children_rf
                lt_v['sort'] = n
                lt_v['color'] = CSS_COLOR_NAMES[i]
                children.append(lt_v)
            c_v['children'] = children
            c_v['sort'] = n
            c_v['color'] = CSS_COLOR_NAMES[i]
            tree['children'].append(c_v)
            i += 1

        return json.dumps(tree)


class GenericNumberingIndex(TemplateView):
    """
    Starting page of generic numbering assignment workflow.
    """
    template_name = 'common_structural_tools.html'

    #Left panel
    step = 1
    number_of_steps = 2
    documentation_url = settings.DOCUMENTATION_URL
    docs = 'structures.html#pdb-file-residue-numbering'
    title = "UPLOAD A PDB FILE"
    description = """
    Upload a pdb file to be annotated with generic numbers from GPCRdb.

    The numbers can be visualized in molecular viewers such as PyMOL, with scripts available with the output files.

    Once you have selected all your targets, click the green button.
    """

    #Input file form data
    header = "Select a file to upload:"
    upload_form_data = {
        "pdb_file": forms.FileField(),
        }
    form_code = forms.Form()
    form_code.fields = upload_form_data
    form_id = 'gn_pdb_file'
    url = '/structure/generic_numbering_results'
    mid_section = "upload_file_form.html"
    form_height = 200
    #Buttons
    buttons = {
        'continue' : {
            'label' : 'Assign generic numbers',
            'color' : 'success',
            },
        }


    def get_context_data (self, **kwargs):

        context = super(GenericNumberingIndex, self).get_context_data(**kwargs)
        # get attributes of this class and add them to the context
        context['form_code'] = str(self.form_code)
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]

        return context


#Class rendering results from generic numbers assignment
class GenericNumberingResults(TemplateView):

    template_name = 'common_structural_tools.html'

    #Left panel
    step = 1
    number_of_steps = 2
    title = "SELECT SUBSTRUCTURE"
    description = 'Download the desired substructures.'
    #Mid section
    mid_section = 'gn_results.html'
    #Buttons - none


    def post (self, request, *args, **kwargs):

        generic_numbering = GenericNumbering(StringIO(request.FILES['pdb_file'].file.read().decode('UTF-8',"ignore")))
        out_struct = generic_numbering.assign_generic_numbers()
        out_stream = StringIO()
        io = PDBIO()
        io.set_structure(out_struct)
        io.save(out_stream)
        if len(out_stream.getvalue()) > 0:
            request.session['gn_outfile'] = out_stream
            request.session['gn_outfname'] = request.FILES['pdb_file'].name
            self.success = True
        else:
            self.input_file = request.FILES['pdb_file'].name
            self.success = False

        context = super(GenericNumberingResults, self).get_context_data(**kwargs)
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]

        return render(request, self.template_name, context)



class GenericNumberingSelection(AbsSegmentSelection):
    """
    Segment selection for download of annotated substructure.
    """

    step = 2
    number_of_steps = 2

    docs = 'structures.html#pdb-file-residue-numbering'

    #Mid section
    #mid_section = 'segment_selection.html'

    #Right panel
    segment_list = True
    buttons = {
        'continue': {
            'label': 'Download substructure',
            'url': '/structure/generic_numbering_results/substr',
            'color': 'success',
        },
    }
    # OrderedDict to preserve the order of the boxes
    selection_boxes = OrderedDict([('reference', False),
        ('targets', False),
        ('segments', True),])


    def get_context_data(self, **kwargs):

        context = super(GenericNumberingSelection, self).get_context_data(**kwargs)

        simple_selection = self.request.session.get('selection', False)
        selection = Selection()
        if simple_selection:
            selection.importer(simple_selection)

        context['selection'] = {}
        context['selection']['site_residue_groups'] = selection.site_residue_groups
        context['selection']['active_site_residue_group'] = selection.active_site_residue_group
        for selection_box, include in self.selection_boxes.items():
            if include:
                context['selection'][selection_box] = selection.dict(selection_box)['selection'][selection_box]

        # get attributes of this class and add them to the context
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]
        return context


class GenericNumberingDownload(View):
    """
    Serve the (sub)structure depending on user's choice.
    """
    def get(self, request, *args, **kwargs):

        if self.kwargs['substructure'] == 'custom':
            return HttpResponseRedirect('/structure/generic_numbering_selection')

        simple_selection = self.request.session.get('selection', False)
        selection = Selection()
        if simple_selection:
            selection.importer(simple_selection)
        out_stream = StringIO()
        io = PDBIO()
        request.session['gn_outfile'].seek(0)
        gn_struct = PDBParser(PERMISSIVE=True, QUIET=True).get_structure(request.session['gn_outfname'], request.session['gn_outfile'])[0]

        if self.kwargs['substructure'] == 'full':
            io.set_structure(gn_struct)
            io.save(out_stream)

        if self.kwargs['substructure'] == 'substr':
            io.set_structure(gn_struct)
            io.save(out_stream, GenericNumbersSelector(parsed_selection=SelectionParser(selection)))

        root, ext = os.path.splitext(request.session['gn_outfname'])
        response = HttpResponse(content_type="chemical/x-pdb")
        response['Content-Disposition'] = 'attachment; filename="{}_GPCRDB.pdb"'.format(root)
        response.write(out_stream.getvalue())

        return response


#==============================================================================

#========================Superposition of structures===========================
#Class for starting page of superposition workflow
class SuperpositionWorkflowIndex(TemplateView):

    template_name = "common_structural_tools.html"

    #Left panel
    step = 1
    number_of_steps = 3
    documentation_url = settings.DOCUMENTATION_URL
    docs = 'structures.html#structure-superposition'
    title = "UPLOAD YOUR FILES"
    description = """
    Upload a pdb file for reference structure, and one or more files that will be superposed. You can also select the structures from crystal structure browser.

    Once you have uploaded/selected all your targets, click the green button.
    """

    header = "Upload or select your structures:"
    #
    upload_form_data = OrderedDict([
        ('ref_file', forms.FileField(label="Reference structure")),
        ('alt_files', MultiFileField(label="Structure(s) to superpose", max_num=10, min_num=1)),
        #('exclusive', forms.BooleanField(label='Download only superposed subset of atoms', widget=forms.CheckboxInput())),
        ])
    form_code = forms.Form()
    form_code.fields = upload_form_data
    form_id = 'superpose_files'
    url = '/structure/superposition_workflow_selection'
    mid_section = 'superposition_workflow_upload_file_form.html'

    #Buttons
    buttons = {
        'continue' : {
            'label' : 'Select segments',
            'color' : 'success',
            }
        }

    # OrderedDict to preserve the order of the boxes
    selection_boxes = OrderedDict([('reference', True),
        ('targets', True),
        ('segments', False)])

    def get_context_data (self, **kwargs):

        context = super(SuperpositionWorkflowIndex, self).get_context_data(**kwargs)

        # get selection from session and add to context
        # get simple selection from session
        simple_selection = self.request.session.get('selection', False)
        # print(simple_selection)
        # create full selection and import simple selection (if it exists)
        selection = Selection()
        if simple_selection:
            selection.importer(simple_selection)
        # print(self.kwargs.keys())
        #Clearing selections for fresh run
        if 'clear' in self.kwargs.keys():
            selection.clear('reference')
            selection.clear('targets')
            selection.clear('segments')
            if 'alt_files' in self.request.session.keys():
                del self.request.session['alt_files']
            if 'ref_file' in self.request.session.keys():
                del self.request.session['ref_file']
        context['selection'] = {}
        for selection_box, include in self.selection_boxes.items():
            if include:
                context['selection'][selection_box] = selection.dict(selection_box)['selection'][selection_box]

        # get attributes of this class and add them to the context
        context['form_code'] = str(self.form_code)
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]

        return context



#Class rendering selection box for sequence segments
class SuperpositionWorkflowSelection(AbsSegmentSelection):

    #Left panel
    step = 2
    number_of_steps = 3

    docs = 'structures.html#structure-superposition'

    #Mid section
    #mid_section = 'segment_selection.html'

    #Right panel
    segment_list = True
    buttons = {
        'continue': {
            'label': 'Superpose proteins',
            'url': '/structure/superposition_workflow_results',
            'color': 'success',
        },
    }
    # OrderedDict to preserve the order of the boxes
    selection_boxes = OrderedDict([('reference', False),
        ('targets', False),
        ('segments', True),])


    def post (self, request, *args, **kwargs):

        # create full selection and import simple selection (if it exists)
        simple_selection = request.session.get('selection', False)
        selection = Selection()
        if simple_selection:
            selection.importer(simple_selection)
        
        if 'ref_file' in request.FILES:
            request.session['ref_file'] = request.FILES['ref_file']
        if 'alt_files' in request.FILES:
            request.session['alt_files'] = request.FILES.getlist('alt_files')

        context = super(SuperpositionWorkflowSelection, self).get_context_data(**kwargs)
        context['selection'] = {}
        for selection_box, include in self.selection_boxes.items():
            if include:
                context['selection'][selection_box] = selection.dict(selection_box)['selection'][selection_box]


        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]

        return render(request, self.template_name, context)

    def get_context_data(self, **kwargs):

        self.buttons = {
            'continue': {
                'label': 'Download substructures',
                'url': '/structure/superposition_workflow_results/custom',
                'color': 'success',
            },
        }
        context = super(SuperpositionWorkflowSelection, self).get_context_data(**kwargs)

        simple_selection = self.request.session.get('selection', False)
        selection = Selection()
        if simple_selection:
            selection.importer(simple_selection)

        context['selection'] = {}
        context['selection']['site_residue_groups'] = selection.site_residue_groups
        context['selection']['active_site_residue_group'] = selection.active_site_residue_group
        for selection_box, include in self.selection_boxes.items():
            if include:
                context['selection'][selection_box] = selection.dict(selection_box)['selection'][selection_box]

        # get attributes of this class and add them to the context
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]
        return context


#Class rendering results from superposition workflow
class SuperpositionWorkflowResults(TemplateView):
    """
    Select download mode for the superposed structures. Full structures, superposed fragments only, select substructure.
    """

    template_name = 'common_structural_tools.html'

    #Left panel
    step = 3
    number_of_steps = 3
    title = "SELECT SUBSTRUCTURE"
    description = 'Download the desired substructures.'

    #Mid section
    mid_section = 'superposition_results.html'
    #Buttons - none


    def get_context_data (self, **kwargs):

        context = super(SuperpositionWorkflowResults, self).get_context_data(**kwargs)

        simple_selection = self.request.session.get('selection', False)
        selection = Selection()
        if simple_selection:
            selection.importer(simple_selection)

        if 'ref_file' in self.request.session.keys():
            ref_file = StringIO(self.request.session['ref_file'].file.read().decode('UTF-8'))
        elif selection.reference != []:
            ref_file = StringIO(selection.reference[0].item.get_cleaned_pdb())
        if 'alt_files' in self.request.session.keys():
            alt_files = [StringIO(alt_file.file.read().decode('UTF-8')) for alt_file in self.request.session['alt_files']]
        elif selection.targets != []:
            alt_files = [StringIO(x.item.get_cleaned_pdb()) for x in selection.targets if x.type in ['structure', 'structure_model', 'structure_model_Inactive', 'structure_model_Intermediate', 'structure_model_Active']]

        superposition = ProteinSuperpose(deepcopy(ref_file),alt_files, selection)
        out_structs = superposition.run()
        if 'alt_files' in self.request.session.keys():
            alt_file_names = [x.name for x in self.request.session['alt_files']]
        else:
            alt_file_names = []
            for x in selection.targets:
                if x.type=='structure':
                    alt_file_names.append('{}_{}.pdb'.format(x.item.protein_conformation.protein.entry_name, x.item.pdb_code.index))
                elif x.type=='structure_model' or x.type=='structure_model_Inactive' or x.type=='structure_model_Intermediate' or x.type=='structure_model_Active':
                    alt_file_names.append('Class{}_{}_{}_{}_GPCRdb.pdb'.format(class_tree[x.item.protein.family.slug[:3]], x.item.protein.entry_name, x.item.state.name, x.item.main_template.pdb_code.index))
        if len(out_structs) == 0:
            self.success = False
        elif len(out_structs) >= 1:
            io = PDBIO()
            self.request.session['alt_structs'] = {}
            for alt_struct, alt_file_name in zip(out_structs, alt_file_names):
                tmp = StringIO()
                io.set_structure(alt_struct)
                io.save(tmp)
                self.request.session['alt_structs'][alt_file_name] = tmp

            self.success = True

        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]
        return context

class SuperpositionWorkflowDownload(View):
    """
    Serve the (sub)structures depending on user's choice.
    """

    def get(self, request, *args, **kwargs):

        if self.kwargs['substructure'] == 'select':
            return HttpResponseRedirect('/structure/superposition_workflow_selection')

        io = PDBIO()
        out_stream = BytesIO()
        zipf = zipfile.ZipFile(out_stream, 'w')
        simple_selection = self.request.session.get('selection', False)
        selection = Selection()
        if simple_selection:
            selection.importer(simple_selection)
        self.alt_substructure_mapping = {}
        #reference
        if 'ref_file' in request.session.keys():
            self.request.session['ref_file'].file.seek(0)
            ref_struct = PDBParser(PERMISSIVE=True, QUIET=True).get_structure('ref', StringIO(self.request.session['ref_file'].file.read().decode('UTF-8')))[0]
            gn_assigner = GenericNumbering(structure=ref_struct)
            gn_assigner.assign_generic_numbers()
            self.ref_substructure_mapping = gn_assigner.get_substructure_mapping_dict()
            ref_name = self.request.session['ref_file'].name
        elif selection.reference != []:
            ref_struct = PDBParser(PERMISSIVE=True, QUIET=True).get_structure('ref', StringIO(selection.reference[0].item.get_cleaned_pdb()))[0]
            gn_assigner = GenericNumbering(structure=ref_struct)
            gn_assigner.assign_generic_numbers()
            self.ref_substructure_mapping = gn_assigner.get_substructure_mapping_dict()
            if selection.reference[0].type=='structure':
                ref_name = '{}_{}_ref.pdb'.format(selection.reference[0].item.protein_conformation.protein.entry_name, selection.reference[0].item.pdb_code.index)
            elif selection.reference[0].type=='structure_model' or selection.reference[0].type=='structure_model_Inactive' or selection.reference[0].type=='structure_model_Intermediate' or selection.reference[0].type=='structure_model_Active':
                ref_name = 'Class{}_{}_{}_{}_GPCRdb_ref.pdb'.format(class_tree[selection.reference[0].item.protein.family.slug[:3]], selection.reference[0].item.protein.entry_name, 
                                                                    selection.reference[0].item.state.name, selection.reference[0].item.main_template.pdb_code.index)

        alt_structs = {}
        for alt_id, st in self.request.session['alt_structs'].items():
            st.seek(0)
            alt_structs[alt_id] = PDBParser(PERMISSIVE=True, QUIET=True).get_structure(alt_id, st)[0]
            gn_assigner = GenericNumbering(structure=alt_structs[alt_id])
            gn_assigner.assign_generic_numbers()
            self.alt_substructure_mapping[alt_id] = gn_assigner.get_substructure_mapping_dict()

        if self.kwargs['substructure'] == 'full':

            io.set_structure(ref_struct)
            tmp = StringIO()
            io.save(tmp)
            zipf.writestr(ref_name, tmp.getvalue())

            for alt_name in self.request.session['alt_structs']:
                tmp = StringIO()
                io.set_structure(alt_structs[alt_name])
                io.save(tmp)
                zipf.writestr(alt_name, tmp.getvalue())

        elif self.kwargs['substructure'] == 'substr':

            consensus_gn_set = CASelector(SelectionParser(selection), ref_struct, alt_structs.values()).get_consensus_gn_set()
            io.set_structure(ref_struct)
            tmp = StringIO()
            io.save(tmp, GenericNumbersSelector(consensus_gn_set))
            zipf.writestr(ref_name, tmp.getvalue())
            for alt_name in self.request.session['alt_structs']:
                tmp = StringIO()
                io.set_structure(alt_structs[alt_name])
                io.save(tmp, GenericNumbersSelector(consensus_gn_set))
                zipf.writestr(alt_name, tmp.getvalue())

        elif self.kwargs['substructure'] == 'custom':

            io.set_structure(ref_struct)
            tmp = StringIO()
            io.save(tmp, SubstructureSelector(self.ref_substructure_mapping, parsed_selection=SelectionParser(selection)))
            
            zipf.writestr(ref_name, tmp.getvalue())
            for alt_name in self.request.session['alt_structs']:
                tmp = StringIO()
                io.set_structure(alt_structs[alt_name])
                io.save(tmp, SubstructureSelector(self.alt_substructure_mapping[alt_name], parsed_selection=SelectionParser(selection)))
                zipf.writestr(alt_name, tmp.getvalue())

        zipf.close()
        if len(out_stream.getvalue()) > 0:
            response = HttpResponse(content_type="application/zip")
            response['Content-Disposition'] = 'attachment; filename="Superposed_structures.zip"'
            response.write(out_stream.getvalue())

        if 'ref_file' in request.FILES:
            request.session['ref_file'] = request.FILES['ref_file']
        if 'alt_files' in request.FILES:
            request.session['alt_files'] = request.FILES.getlist('alt_files')


        return response


class FragmentSuperpositionIndex(TemplateView):

    template_name = 'common_structural_tools.html'

    #Left panel
    step = 1
    number_of_steps = 1

    documentation_url = settings.DOCUMENTATION_URL
    docs = 'sites.html#pharmacophore-generation'

    title = "SUPERPOSE FRAGMENTS OF CRYSTAL STRUCTURES"
    description = """
    The tool implements a fragment-based pharmacophore method, as published in <a href='http://www.ncbi.nlm.nih.gov/pubmed/25286328'>Fidom K, et al (2015)</a>. Interacting ligand moiety - residue pairs extracted from selected crystal structures of GPCRs are superposed onto the input pdb file based on gpcrdb generic residue numbers. Resulting aligned ligand fragments can be used for placement of pharmacophore features.

    Upload a pdb file you want to superpose the interacting moiety - residue pairs.

    Once you have selected all your targets, click the green button.
        """

    #Input file form data
    header = "Select a file to upload:"
    #Can't control the class properly - staying with the dirty explicit html code
    form_id='fragments'
    form_code = """
    Pdb file:<input id="id_pdb_file" name="pdb_file" type="file" /></br>
    Similarity:</br>
    <input id="similarity" name="similarity" type="radio" value="identical" /> Use fragments with identical residues</br>
    <input checked="checked" id="similarity" name="similarity" type="radio" value="similar" /> Use fragments with residues of similar properties</br>

    Fragments:</br>
    <input checked="checked" id="representative" name="representative" type="radio" value="closest" /> Use fragments from the evolutionary closest crystal structure</br>
    <input id="representative" name="representative" type="radio" value="any" /> Use all available fragments</br></br>
    State:<select id="id_state" name="state">
    <option value="active">Antagonist-bound structures</option>
    <option value="inactive">Agonist-bound structures</option>
    </select>
    """
    url = '/structure/fragment_superposition_results'
    mid_section = "upload_file_form.html"
    form_height = 350

    #Buttons
    buttons = {
        'continue' : {
            'label' : 'Retrieve fragments',
            'color' : 'success',
            },
        }


    def get_context_data (self, **kwargs):

        context = super(FragmentSuperpositionIndex, self).get_context_data(**kwargs)
        # get attributes of this class and add them to the context
        context['form_code'] = str(self.form_code)
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]

        return context



class FragmentSuperpositionResults(TemplateView):

    template_name = "common_structural_tools.html"

    #Left panel - blank
    #Mid section
    mid_section = 'fragment_superposition_results.html'
    #Buttons - none

    def post (self, request, *args, **kwargs):

        frag_sp = FragmentSuperpose(StringIO(request.FILES['pdb_file'].file.read().decode('UTF-8', 'ignore')),request.FILES['pdb_file'].name)
        superposed_fragments = []
        superposed_fragments_repr = []
        if request.POST['similarity'] == 'identical':
            if request.POST['representative'] == 'any':
                superposed_fragments = frag_sp.superpose_fragments()
            else:
                superposed_fragments_repr = frag_sp.superpose_fragments(representative=True, state=request.POST['state'])
                superposed_fragments = frag_sp.superpose_fragments()
        else:
            if request.POST['representative'] == 'any':
                superposed_fragments = frag_sp.superpose_fragments(use_similar=True)
            else:
                superposed_fragments_repr = frag_sp.superpose_fragments(representative=True, use_similar=True, state=request.POST['state'])
                superposed_fragments = frag_sp.superpose_fragments(use_similar=True)
        if superposed_fragments == []  and superposed_fragments_repr == []:
            self.message = "No fragments were aligned."
        else:
            io = PDBIO()
            out_stream = BytesIO()
            zipf = zipfile.ZipFile(out_stream, 'a', zipfile.ZIP_DEFLATED)
            for fragment, pdb_data in superposed_fragments:
                io.set_structure(pdb_data)
                tmp = StringIO()
                io.save(tmp)
                if request.POST['representative'] == 'any':
                    zipf.writestr(fragment.generate_filename(), tmp.getvalue())
                else:
                    zipf.writestr("all_fragments//{!s}".format(fragment.generate_filename()), tmp.getvalue())
            if superposed_fragments_repr != []:
                for fragment, pdb_data in superposed_fragments_repr:
                    io.set_structure(pdb_data)
                    tmp = StringIO()
                    io.save(tmp)
                    zipf.writestr("representative_fragments//{!s}".format(fragment.generate_filename()), tmp.getvalue())
            zipf.close()
            if len(out_stream.getvalue()) > 0:
                request.session['outfile'] = { 'interacting_moiety_residue_fragments.zip' : out_stream, }
                self.outfile = 'interacting_moiety_residue_fragments.zip'
                self.success = True
                self.zip = 'zip'
                self.message = '{:n} fragments were superposed.'.format(len(superposed_fragments))

        context = super(FragmentSuperpositionResults, self).get_context_data(**kwargs)
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]

        return render(request, self.template_name, context)


#==============================================================================
class TemplateTargetSelection(AbsReferenceSelection):
    """
    Starting point for template selection workflow. Target selection.
    """

    type_of_selection = 'reference'
    # Left panel
    description = 'Select a reference target by searching or browsing in the middle column.' \
        + '\n\nThe selected reference target will appear in the right column.' \
        + '\n\nOnce you have selected your reference target, either proceed with all TMs alignment ("Find template"' \
        + 'button) or specify the sequence segments manualy ("Advanced segment selection" button).'
    step = 1
    number_of_steps = 2
    redirect_on_select = False

    docs = 'structures.html#template-selection'

    # Mid section

    # Right panel
    buttons = OrderedDict()
    buttons['continue'] = {
        'label': 'Find template',
        'url': '/structure/template_browser',
        'color': 'success',
    }
    buttons['segments'] = {
        'label' : 'Advanced segment selection',
        'url' : '/structure/template_segment_selection',
        'color' : 'info',
    }

    selection_boxes = OrderedDict([('reference', True),
        ('targets', False),
        ('segments', False),])



#==============================================================================
class TemplateSegmentSelection(AbsSegmentSelection):
    """
    Advanced selection of sequence segments for template search.
    """
   #Left panel
    step = 2
    number_of_steps = 2

    docs = 'structures.html#template-selection'

    #Mid section
    #mid_section = 'segment_selection.html'

    #Right panel
    segment_list = True
    buttons = {
        'continue': {
            'label': 'Find template',
            'url': '/structure/template_browser',
            'color': 'success',
        },
    }
    # OrderedDict to preserve the order of the boxes
    selection_boxes = OrderedDict([('reference', True),
        ('targets', False),
        ('segments', True),])



#==============================================================================
class TemplateBrowser(TemplateView):
    """
    Fetching Structure data and ordering by similarity
    """

    template_name = "template_browser.html"

    def get_context_data (self, **kwargs):

        context = super(TemplateBrowser, self).get_context_data(**kwargs)

        # get simple selection from session
        simple_selection = self.request.session.get('selection', False)

        # make an alignment
        a = Alignment()
        a.ignore_alternative_residue_numbering_schemes = True

        # load the selected reference into the alignment
        a.load_reference_protein_from_selection(simple_selection)

        # fetch
        qs = Structure.objects.filter(refined=False).select_related(
            "pdb_code__web_resource",
            "protein_conformation__protein__species",
            "protein_conformation__protein__source",
            "protein_conformation__protein__family__parent__parent__parent",
            "publication__web_link__web_resource").prefetch_related(
            "stabilizing_agents",
            "protein_conformation__protein__parent__endogenous_ligands__properities__ligand_type",
            Prefetch("ligands", queryset=StructureLigandInteraction.objects.filter(
            annotated=True).prefetch_related('ligand__properities__ligand_type', 'ligand_role')))

        # Dirty but fast
        qsd = {}
        for st in qs:
            qsd[st.protein_conformation.protein.id] = st

        # add proteins to the alignment
        a.load_proteins([x.protein_conformation.protein for x in qs])

        if simple_selection.segments != []:
            a.load_segments_from_selection(simple_selection)
        else:
            a.load_segments(ProteinSegment.objects.filter(slug__in=['TM1', 'TM2', 'TM3', 'TM4','TM5','TM6', 'TM7']))

        a.build_alignment()
        a.calculate_similarity()

        context['structures'] = []
        for prot in a.proteins[1:]:
            try:
                context['structures'].append([prot.similarity, prot.identity, qsd[prot.protein.id]])
                del qsd[prot.protein.id]
            except KeyError:
                pass
        return context

#==============================================================================
class PDBClean(TemplateView):
    """
    Extraction, packing and serving out the pdb records selected via structure/template browser.
    """
    template_name = "pdb_download.html"

    def post(self, request, *args, **kwargs):
        context = super(PDBClean, self).get_context_data(**kwargs)

        class_dict = {'001':'A','002':'B1','003':'B2','004':'C','005':'F','006':'T','007':'O'}

        self.posted = True
        pref = True
        water = False
        hets = False

        if 'pref_chain' not in request.POST.keys():
            pref = False
        if 'water' in request.POST.keys():
            water = True
        if 'hets' in request.POST.keys():
            hets = True

        # get simple selection from session
        simple_selection = request.session.get('selection', False)
        selection = Selection()
        if simple_selection:
            selection.importer(simple_selection)
        out_stream = BytesIO()
        io = PDBIO()
        zipf = zipfile.ZipFile(out_stream, 'w', zipfile.ZIP_DEFLATED)
        if selection.targets != []:
            if selection.targets != [] and selection.targets[0].type == 'structure':
                for selected_struct in [x for x in selection.targets if x.type == 'structure']:
                    struct_name = '{}_{}.pdb'.format(selected_struct.item.protein_conformation.protein.parent.entry_name, selected_struct.item.pdb_code.index)
                    if hets:
                        lig_names = [x.pdb_reference for x in StructureLigandInteraction.objects.filter(structure=selected_struct.item, annotated=True)]
                    else:
                        lig_names = None
                    gn_assigner = GenericNumbering(structure=PDBParser(QUIET=True).get_structure(struct_name, StringIO(selected_struct.item.get_cleaned_pdb(pref, water, lig_names)))[0])
                    tmp = StringIO()
                    io.set_structure(gn_assigner.assign_generic_numbers())
                    request.session['substructure_mapping'] = gn_assigner.get_substructure_mapping_dict()
                    io.save(tmp)
                    zipf.writestr(struct_name, tmp.getvalue())
                    del gn_assigner, tmp
                for struct in selection.targets:
                    selection.remove('targets', 'structure', struct.item.id)
            elif selection.targets != [] and selection.targets[0].type == 'structure_model':
                for hommod in [x for x in selection.targets if x.type == 'structure_model']:
                    mod_name = 'Class{}_{}_{}_{}_{}_GPCRDB.pdb'.format(class_dict[hommod.item.protein.family.slug[:3]], hommod.item.protein.entry_name, 
                                                                                  hommod.item.state.name, hommod.item.main_template.pdb_code.index, hommod.item.version)
                    tmp = StringIO(hommod.item.pdb)
                    request.session['substructure_mapping'] = 'full'
                    zipf.writestr(mod_name, tmp.getvalue())
                    del tmp
                    rotamers = StructureModelStatsRotamer.objects.filter(homology_model=hommod.item).prefetch_related('residue','backbone_template','rotamer_template').order_by('residue__sequence_number')
                    stats_data = 'Segment,Sequence_number,Generic_number,Backbone_template,Rotamer_template\n'
                    for r in rotamers:
                        try:
                            gn = r.residue.generic_number.label
                        except:
                            gn = '-'
                        if r.backbone_template:
                            bt = r.backbone_template.pdb_code.index
                        else:
                            bt = '-'
                        if r.rotamer_template:
                            rt = r.rotamer_template.pdb_code.index
                        else:
                            rt = '-'
                        stats_data+='{},{},{},{},{}\n'.format(r.residue.protein_segment.slug, r.residue.sequence_number, gn, bt, rt)
                    stats_name = mod_name[:-3]+'templates.csv'
                    zipf.writestr(stats_name, stats_data)
                    del stats_data
                for mod in selection.targets:
                    selection.remove('targets', 'structure_model', mod.item.id)

            # export simple selection that can be serialized
            simple_selection = selection.exporter()

            request.session['selection'] = simple_selection
            request.session['cleaned_structures'] = out_stream

        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]

        if selection.targets != [] and selection.targets[0].type == 'structure_model':
            return zipf
        else:
            return render(request, self.template_name, context)


    def get_context_data (self, **kwargs):

        context = super(PDBClean, self).get_context_data(**kwargs)
        self.success = False
        self.posted = False

        # get simple selection from session
        simple_selection = self.request.session.get('selection', False)
        selection = Selection()
        if simple_selection:
            selection.importer(simple_selection)
        if selection.targets != []:
            self.success = True

        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]

        return context


class PDBSegmentSelection(AbsSegmentSelection):

    #Left panel
    step = 2
    number_of_steps = 2

    #Mid section
    #mid_section = 'segment_selection.html'

    #Right panel
    segment_list = True

    # OrderedDict to preserve the order of the boxes
    selection_boxes = OrderedDict([('reference', False),
        ('targets', False),
        ('segments', True),])

    def get_context_data(self, **kwargs):

        self.buttons = {
            'continue': {
                'label': 'Download substructures',
                'url': '/structure/pdb_download/custom',
                'color': 'success',
            },
        }
        context = super(PDBSegmentSelection, self).get_context_data(**kwargs)

        simple_selection = self.request.session.get('selection', False)
        selection = Selection()
        if simple_selection:
            selection.importer(simple_selection)

        context['selection'] = {}
        context['selection']['site_residue_groups'] = selection.site_residue_groups
        context['selection']['active_site_residue_group'] = selection.active_site_residue_group
        for selection_box, include in self.selection_boxes.items():
            if include:
                context['selection'][selection_box] = selection.dict(selection_box)['selection'][selection_box]

        # get attributes of this class and add them to the context
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]
        return context

#==============================================================================
class PDBDownload(View):
    """
    Serve the PDB (sub)structures depending on user's choice.
    """

    def get(self, request, hommods=False, *args, **kwargs):
        if self.kwargs['substructure'] == 'select':
            return HttpResponseRedirect('/structure/pdb_segment_selection')

        if self.kwargs['substructure'] == 'full':
            out_stream = request.session['cleaned_structures']

        elif self.kwargs['substructure'] == 'custom':
            simple_selection = request.session.get('selection', False)
            selection = Selection()
            if simple_selection:
                selection.importer(simple_selection)
            io = PDBIO()
            zipf_in = zipfile.ZipFile(request.session['cleaned_structures'], 'r')
            out_stream = BytesIO()
            zipf_out = zipfile.ZipFile(out_stream, 'w', zipfile.ZIP_DEFLATED)
            for name in zipf_in.namelist():
                tmp = StringIO()
                io.set_structure(PDBParser(QUIET=True).get_structure(name, StringIO(zipf_in.read(name).decode('utf-8')))[0])
                io.save(tmp, SubstructureSelector(request.session['substructure_mapping'], parsed_selection=SelectionParser(selection)))
                zipf_out.writestr(name, tmp.getvalue())

            zipf_in.close()
            zipf_out.close()
            del request.session['substructure_mapping']
        if len(out_stream.getvalue()) > 0:
            response = HttpResponse(content_type="application/zip")
            if hommods == False:
                response['Content-Disposition'] = 'attachment; filename="pdb_structures.zip"'
            else:
                response['Content-Disposition'] = 'attachment; filename="GPCRDB_homology_models.zip"'
            response.write(out_stream.getvalue())

        return response

#==============================================================================
def ConvertStructuresToProteins(request):
    "For alignment from structure browser"

    simple_selection = request.session.get('selection', False)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)
    if selection.targets != []:
        for struct in selection.targets:
            if 'refined' in struct.item.pdb_code.index:
                prot = struct.item.protein_conformation.protein
            else:
                prot = struct.item.protein_conformation.protein.parent
            selection.remove('targets', 'structure', struct.item.id)
            selection.add('targets', 'protein', SelectionItem('protein', prot))
        if selection.reference != []:
            selection.add('targets', 'protein', selection.reference[0])
    # export simple selection that can be serialized
    simple_selection = selection.exporter()

    # add simple selection to session
    request.session['selection'] = simple_selection

    return HttpResponseRedirect('/alignment/segmentselection')


def ConvertStructureModelsToProteins(request):
    "For alignment from homology model browser"
    simple_selection = request.session.get('selection', False)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)
    if selection.targets != []:
        for struct_mod in selection.targets:
            try:
                prot = struct_mod.item.protein
                selection.remove('targets', 'structure_model', struct_mod.item.id)
                selection.add('targets', 'protein', SelectionItem('protein', prot))
            except:
                prot = struct_mod.item.protein_conformation.protein.parent
                selection.remove('targets', 'structure', struct_mod.item.id)
                selection.add('targets', 'protein', SelectionItem('protein', prot))
        if selection.reference != []:
            selection.add('targets', 'protein', selection.reference[0])
    # export simple selection that can be serialized
    simple_selection = selection.exporter()

    # add simple selection to session
    request.session['selection'] = simple_selection

    return HttpResponseRedirect('/alignment/segmentselection')


def HommodDownload(request):
    "Download selected homology models in zip file"
    p = PDBClean()
    p.post(request)
    p1 = PDBDownload()
    p1.kwargs = {}
    p1.kwargs['substructure'] = 'full'
    p2 = p1.get(request, hommods=True)
    return p2

def SingleModelDownload(request, modelname, state, csv=False):
    "Download single homology model"
    class_dict = {'001':'A','002':'B1','003':'B2','004':'C','005':'F','006':'T','007':'O'}
    if state=='refined':
        hommod = Structure.objects.get(pdb_code__index=modelname+'_refined')
    else:
        hommod = StructureModel.objects.get(protein__entry_name=modelname, state__slug=state)
    if csv:
        if state=='refined':
            rotamers = StructureRefinedStatsRotamer.objects.filter(structure=hommod).prefetch_related("structure", "residue", "backbone_template", "rotamer_template").order_by('residue__sequence_number')
        else:
            rotamers = StructureModelStatsRotamer.objects.filter(homology_model=hommod).prefetch_related("homology_model", "residue", "backbone_template", "rotamer_template").order_by('residue__sequence_number')
        text_out = "Segment,Sequence_number,Generic_number,Backbone_template,Rotamer_template\n"
        for r in rotamers:
            if r.backbone_template:
                bt = r.backbone_template.pdb_code.index
            else:
                bt = '-'
            if r.rotamer_template:
                rt = r.rotamer_template.pdb_code.index
            else:
                rt = '-'
            if r.residue.generic_number:
                gn = r.residue.generic_number.label
            else:
                gn = '-'
            text_out+='{},{},{},{},{}\n'.format(r.residue.protein_segment.slug, r.residue.sequence_number, gn, bt, rt)
        response = HttpResponse(text_out, content_type="homology_models/csv")
        if state=='refined':
            file_name = 'Class{}_{}_{}_GPCRDB.templates.csv'.format(class_dict[hommod.protein_conformation.protein.family.slug[:3]], hommod.protein_conformation.protein.entry_name,
                                                                               hommod.pdb_code.index)
        else:
            file_name = 'Class{}_{}_{}_{}_{}_GPCRDB.templates.csv'.format(class_dict[hommod.protein.family.slug[:3]], hommod.protein.entry_name, 
                                                                                     hommod.state.name, hommod.main_template.pdb_code.index, hommod.version)
    else:
        if state=='refined':
            response = HttpResponse(hommod.pdb_data.pdb, content_type="homology_models/model")
        else:
            response = HttpResponse(hommod.pdb, content_type="homology_models/model")
        if state=='refined':
            file_name = 'Class{}_{}_{}_GPCRDB.pdb'.format(class_dict[hommod.protein_conformation.protein.family.slug[:3]], hommod.protein_conformation.protein.entry_name,
                                                                     hommod.pdb_code.index)
        else:
            file_name = 'Class{}_{}_{}_{}_{}_GPCRDB.pdb'.format(class_dict[hommod.protein.family.slug[:3]], hommod.protein.entry_name, 
                                                                           hommod.state.name, hommod.main_template.pdb_code.index, hommod.version)
    response['Content-Disposition'] = 'attachment; filename="{}"'.format(file_name)

    return response

def ServePdbOutfile (request, outfile, replacement_tag):

    root, ext = os.path.splitext(outfile)
    out_stream = request.session['outfile'][outfile]
    response = HttpResponse(content_type="chemical/x-pdb")
    response['Content-Disposition'] = 'attachment; filename="{}_{}.pdb"'.format(root, replacement_tag)
    response.write(out_stream.getvalue())

    return response


def ServeZipOutfile (request, outfile):

    out_stream = request.session['outfile'][outfile]
    response = HttpResponse(content_type="application/zip")
    response['Content-Disposition'] = 'attachment; filename="{}"'.format(outfile)
    response.write(out_stream.getvalue())

    return response

def RenderTrees(request):
    number = request.GET['number']
    tree = open(settings.STATICFILES_DIRS[0] +'/home/images/00'+number+'_tree.xml').read()
    legend = open(settings.STATICFILES_DIRS[0] +'/home/images/00'+number+'_legend.svg').read()
    context = {'tree':tree, 'leg':legend, 'num':number}
    return render(request, 'phylogenetic_trees.html', context)

def webform(request):
    form = construct_form()
    context = {'form':form}
    return render(request, 'web_form.html',context)

def webform_two(request, slug=None):
    context = {}
    if slug:
        c = Construct.objects.filter(name=slug).get()
        # print(c.json)
        # test = ast.literal_eval(c.json)
        # print(test)
        json_data = json.loads(c.json)
        if 'raw_data' not in json_data:
            json_data = convert_ordered_to_disordered_annotation(json_data)
        else:
            if 'csrfmiddlewaretoken' in json_data['raw_data']:
                del json_data['raw_data']['csrfmiddlewaretoken'] #remove to prevent errors

        context = {'edit':json.dumps(json_data)}
    return render(request, 'web_form_2.html',context)

def webformdata(request) :

    data = request.POST
    raw_data = deepcopy(data)
    purge_keys = ('Please Select','aamod_position','wt_aa','mut_aa','insert_pos_type','protein_type','deletion','csrfmiddlewaretoken')
    data = dict((k, v) for k, v in data.items() if v!='' and v!='Please Select') #remove empty
    deletions = []
    mutations = []
    contact_info= OrderedDict()
    construct_crystal=OrderedDict()
    auxiliary=OrderedDict()
    expression=OrderedDict()
    solubilization = OrderedDict()
    crystallization = OrderedDict()
    modifications = []
    aamod, aamod_start, aamod_end=[], [], []

    i=1
    error = 0
    error_msg = []
    for key,value in sorted(data.items()):
        try:
            if key.startswith('delet_start'):
                deletions.append({'start':value, 'end':data[key.replace('start','end')], 'origin':'user', 'type':'range'})
                data.pop(key, None)
                data.pop(key.replace('start','end'), None)
            elif key.startswith('ins_start'):
                deletions.append({'start':value, 'end':data[key.replace('start','end')], 'origin':'insertion'+key.replace('ins_start',''), 'type':'range'})
                data.pop(key, None)
                data.pop(key.replace('start','end'), None)
                data.pop(key.replace('ins_start',''), None)
            elif key.startswith(('deletion_single', 'insert_pos_single')):
                if key.startswith('insert_pos_single'):
                    deletions.append({'pos':value, 'origin':'insertion'+key.replace('insert_pos_single',''), 'type':'single'})
                    data.pop(key.replace('insert_pos_single',''), None)
                else:
                    deletions.append({'pos':value, 'origin':'user', 'type':'single'})
                data.pop(key, None)

            if key.startswith('aa_no'):
                pos_id = key.replace('aa_no','')
                if pos_id=='':
                    mut_id='1'
                else:
                    mut_id=pos_id.replace('_','')

                if 'mut_remark'+pos_id in data:
                    remark = data['mut_remark'+pos_id]
                else:
                    remark = ''

                mutations.append({'pos':value,'wt':data['wt_aa'+pos_id],'mut':data['mut_aa'+pos_id], 'type':data['mut_type'+pos_id], 'remark':remark})
                data.pop(key, None)
                data.pop('wt_aa'+pos_id, None)
                data.pop('mut_aa'+pos_id, None)
                data.pop('mut_type'+pos_id, None)

            if key.startswith(('date','name_cont', 'pi_name',
                'pi_address','address','url','pi_email' )):
                contact_info[key]=value
                data.pop(key, None)

            if key.startswith(('pdb', 'pdb_name',
                'uniprot','ligand_name', 'ligand_activity', 'ligand_conc', 'ligand_conc_unit','ligand_id','ligand_id_type')):
                construct_crystal[key]=value
                data.pop(key, None)

            if key.startswith('position'):
                pos_id = key.replace('position','')
                if pos_id=='':
                    aux_id='1'
                else:
                    aux_id=pos_id.replace('_','')

                if 'aux'+aux_id not in auxiliary:
                    auxiliary['aux'+aux_id] = {'position':value,'type':data['protein_type'+pos_id],'presence':data['presence'+pos_id]}

                    data.pop(key, None)
                    data.pop('protein_type'+pos_id, None)
                    data.pop('presence'+pos_id, None)

            if key.startswith(('tag', 'fusion_prot', 'signal', 'linker_seq','prot_cleavage', 'other_prot_cleavage' )):
                temp = key.split('_')
                if len(temp)==4:
                    pos_id = "_"+temp[3]
                    aux_id=pos_id.replace('_','')
                elif len(temp)==3:
                    pos_id = "_"+temp[2]
                    aux_id=pos_id.replace('_','')
                elif len(temp)==2 and temp[1].isdigit():
                    pos_id = "_"+temp[1]
                    aux_id=pos_id.replace('_','')
                else:
                    pos_id = ''
                    aux_id = '1'
                print(key,aux_id,pos_id)

                if 'aux'+aux_id not in auxiliary:
                    auxiliary['aux'+aux_id] = {'position':data['position'+pos_id],'type':data['protein_type'+pos_id],'presence':data['presence'+pos_id]}

                    data.pop('position'+pos_id, None)
                    data.pop('protein_type'+pos_id, None)
                    data.pop('presence'+pos_id, None)

                # if value=='Other':
                #     auxiliary['aux'+aux_id]['other'] = data['other_'+auxiliary['aux'+aux_id]['type']+pos_id]
                #     data.pop('other_'+auxiliary['aux'+aux_id]['type']+pos_id,None)

                auxiliary['aux'+aux_id]['subtype'] = value
                data.pop(key, None)

            if key.startswith(('expr_method', 'host_cell_type',
                    'host_cell', 'expr_remark','expr_other','other_host','other_host_cell' )):
                expression[key]=value
                data.pop(key, None)

            if key.startswith(('deterg_type','deterg_concentr','deterg_concentr_unit','solub_additive','additive_concentr','addit_concentr_unit','chem_enz_treatment','sol_remark')):
                solubilization[key]=value
                data.pop(key, None)

            elif key.startswith(('crystal_type','crystal_method','other_method','other_crystal_type',
                               'protein_concentr','protein_conc_unit','temperature','ph_single','ph',
                               'ph_range_one','ph_range_two','crystal_remark','lcp_lipid','lcp_add',
                               'lcp_conc','lcp_conc_unit','detergent','deterg_conc','deterg_conc_unit','lipid','lipid_concentr','lipid_concentr_unit',
                               'other_deterg','other_deterg_type', 'other_lcp_lipid','other_lipid')):
                crystallization[key]=value
                data.pop(key, None)

            if key.startswith('chemical_comp') and not key.startswith('chemical_comp_type'):

                if 'chemical_components' not in crystallization:
                    crystallization['chemical_components'] = []

                # print(key)    
                if key!='chemical_comp': #not first
                    comp_id = key.replace('chemical_comp','')
                else:
                    comp_id = ''

                crystallization['chemical_components'].append({'component':value,'type':data['chemical_comp_type'+comp_id],'value':data['concentr'+comp_id],'unit':data['concentr_unit'+comp_id]})
                data.pop(key, None)
                data.pop('concentr'+comp_id, None)
                data.pop('concentr_unit'+comp_id, None)
                data.pop('chemical_comp_type'+comp_id, None)


            if key.startswith('aamod') and not key.startswith('aamod_position') and not key.startswith('aamod_pair') and not key=='aamod_position' and not key=='aamod_single':
                if key!='aamod': #not first
                    mod_id = key.replace('aamod','')
                else:
                    mod_id = ''

                if data['aamod_position'+mod_id]=='single':
                    pos = ['single',data['aamod_single'+mod_id]]
                    data.pop('aamod_single'+mod_id, None)
                elif data['aamod_position'+mod_id]=='range':
                    pos = ['range',[data['aamod_start'+mod_id],data['aamod_end'+mod_id]]]
                    data.pop('aamod_start'+mod_id, None)
                    data.pop('aamod_end'+mod_id, None)
                elif data['aamod_position'+mod_id]=='pair':
                    pos = ['pair',[data['aamod_pair_one'+mod_id],data['aamod_pair_two'+mod_id]]]
                    data.pop('aamod_pair_one'+mod_id, None)
                    data.pop('aamod_pair_two'+mod_id, None)

                remark = ''
                if 'mod_remark'+mod_id in data:
                    remark = data['mod_remark'+mod_id]
                modifications.append({'type':value,'remark':remark,'position':pos })
                data.pop(key, None)
                data.pop('mod_remark'+mod_id, None)
                data.pop('aamod_position'+mod_id, None)

            if key.startswith(purge_keys):
                data.pop(key, None)
        except BaseException as e:
            error_msg.append(str(e))
            error = 1

    auxiliary = OrderedDict(sorted(auxiliary.items()))

    context = OrderedDict( [('contact_info',contact_info), ('construct_crystal',construct_crystal),
                           ('auxiliary' , auxiliary),  ('deletions',deletions), ('mutations',mutations),
                           ('modifications', modifications), ('expression', expression), ('solubilization',solubilization),
                           ('crystallization',crystallization),  ('unparsed',data),  ('raw_data',raw_data), ('error', error), ('error_msg',error_msg)] )

    add_construct(context)

    if error==0:
        dump_dir = '/protwis/construct_dump'
        # dump_dir = '/web/sites/files/construct_data' #for sites
        if not os.path.exists(dump_dir):
            os.makedirs(dump_dir)
        ts = int(time.time())
        json_data = context
        json.dump(json_data, open(dump_dir+"/"+str(ts)+"_"+construct_crystal['pdb']+".json", 'w'), indent=4, separators=(',', ': '))

        context['data'] = sorted(data.items())
        #context['data'] = sorted(raw_data.items())

        recipients = ['christian@munk.be']
        emaillist = [elem.strip().split(',') for elem in recipients]
        msg = MIMEMultipart()
        msg['Subject'] = 'GPCRdb: New webform data'
        msg['From'] = 'gpcrdb@gmail.com'
        msg['Reply-to'] = 'gpcrdb@gmail.com'

        msg.preamble = 'Multipart massage.\n'

        part = MIMEText("Hi, please find the attached file")
        msg.attach(part)

        part = MIMEApplication(open(str(dump_dir+"/"+str(ts)+"_"+construct_crystal['pdb']+".json"),"rb").read())
        part.add_header('Content-Disposition', 'attachment', filename=str(dump_dir+"/"+str(ts)+"_"+construct_crystal['pdb']+".json"))
        msg.attach(part)


        server = smtplib.SMTP("smtp.gmail.com:587")
        server.ehlo()
        server.starttls()
        server.login("gpcrdb@gmail.com", "gpcrdb2016")

        server.sendmail(msg['From'], emaillist , msg.as_string())

        context['filename'] = str(ts)+"_"+construct_crystal['pdb']

    return render(request, 'web_form_results.html', context)

def webform_download(request,slug):
    dump_dir = '/protwis/construct_dump'
    # dump_dir = '/web/sites/files/construct_data' #for sites
    file = dump_dir+"/"+str(slug)+".json"
    out_stream = open(file,"rb").read()
    response = HttpResponse(content_type="application/json")
    response['Content-Disposition'] = 'attachment; filename="{}"'.format(file)
    response.write(out_stream)
    return response
