from django.http import HttpResponse, JsonResponse
from django.shortcuts import render
from django.views.generic import TemplateView
from django.views.decorators.cache import cache_page
from django.views.decorators.csrf import csrf_exempt
from django.conf import settings
from django.db.models import Count, Case, When, Min, Q
from django.core.cache import cache
from django.contrib.postgres.aggregates import ArrayAgg

from common import definitions
Alignment = getattr(__import__('common.alignment_' + settings.SITE_NAME, fromlist=['Alignment']), 'Alignment')

from common.selection import SimpleSelection, Selection, SelectionItem
from ligand.models import AssayExperiment, BiasedData, BalancedLigands
from structure.models import Structure, StructureModel, StructureComplexModel
from protein.models import Protein, ProteinFamily, ProteinSegment, Species, ProteinSource, ProteinSet, ProteinCouplings
from residue.models import ResidueGenericNumber, ResidueNumberingScheme, ResidueGenericNumberEquivalent, ResiduePositionSet, Residue
from interaction.forms import PDBform
from signprot.models import SignprotStructure
from construct.tool import FileUploadForm

from svglib.svglib import SvgRenderer
from reportlab.graphics import renderPDF
from lxml import etree
import inspect
import html
import re
from collections import OrderedDict
from io import BytesIO
import xlsxwriter, xlrd
import time
import json
import urllib

default_schemes_excluded = ["cgn", "ecd", "can"]

def getLigandTable(receptor_id, browser_type):
    cache_key = "reference_table_" + str(receptor_id) + browser_type
    data_table = cache.get(cache_key)
    # data_table = None
    if data_table == None:
        ligands = list(BiasedData.objects.filter(receptor_id=receptor_id).values_list(
            "ligand__name",
            "endogenous_status",
            "ligand__ligand_type__name",
            "ligand__id",
            "ligand__inchikey",
            "ligand__smiles",).distinct())

        data_table = "<table id='uniprot_selection' class='uniprot_selection stripe compact'> \
            <thead>\
              <tr> \
                <th colspan=1>&nbsp;</th> \
                <th colspan=4>Ligands</th> \
                <th colspan=2>Data</th> \
              </tr> \
              <tr> \
                <th><br><br><input autocomplete='off' class='form-check-input' type='checkbox' onclick='return check_all_targets();'></th> \
                <th style=\"width; 100px;\">Ligand name<br>&nbsp;</th> \
                <th>2D structure<br>&nbsp;</th> \
                <th>Ligand type<br>&nbsp;</th> \
                <th>Physiological<br>&nbsp;</th> \
                <th>Publication count</th> \
                <th>Comparable ligands</th> \
              </tr> \
            </thead>\
            \n \
            <tbody>\n"

        #link_setup = "<a target=\"_blank\" href=\"{}\"><span class=\"glyphicon glyphicon-new-window btn-xs\"></span></a>"
        link_setup = "<a target=\"_blank\" href=\"{}\">{}</a>"
        # img_setup_inchi = "<img style=\"max-height: 300px; max-width: 300px;\" src=\"http://www.ebi.ac.uk/chembl/api/data/image/{}.svg\">"
        img_setup_smiles = "<img style=\"max-height: 300px; max-width: 300px;\" src=\"https://cactus.nci.nih.gov/chemical/structure/{}/image\">"

        # Prefetch and aggregate related data (ligand_ids and publication_ids)
        lig_pubs = list(BiasedData.objects.filter(receptor_id=receptor_id).values("ligand_id").
                        annotate(pubs_ids=ArrayAgg("publication_id", distinct=True)))
        lig_pubs = { item["ligand_id"]: item["pubs_ids"] for item in lig_pubs }
        pub_lig_ids = list(BiasedData.objects.filter(receptor_id=receptor_id).values("publication_id")
                         .annotate(lig_ids=ArrayAgg("ligand_id", distinct=True)))
        pub_lig_ids = { item["publication_id"]:item["lig_ids"] for item in pub_lig_ids }
        for p in ligands:
            pubs = lig_pubs[p[3]]
            compared = set()
            for pub_id in pubs:
                compared.update(pub_lig_ids[pub_id])

            t = {}
            t['ligandname'] = link_setup.format(str(p[3])+"/info", p[0])
            t['ligandtype'] = p[2]
            t['endogenous'] = p[1]
            t['publications'] = len(pubs)
            t['compared'] = len(compared) - 1
            t['2d_structure'] = img_setup_smiles.format(urllib.parse.quote(p[5])) if p[5] != None else "Image not available"
            data_table += "<tr> \
            <td data-sort=\"0\"><input autocomplete='off' class=\"form-check-input\" type=\"checkbox\" name=\"reference\" id=\"{}\" data-entry=\"{}\" entry-value=\"{}\"></td> \
            <td data-html=\"true\">{}</td> \
            <td>{}</td> \
            <td>{}</td> \
            <td>{}</td> \
            <td>{}</td> \
            <td>{}</td> \
            </tr> \n".format(
                p[0],
                p[3],
                p[3],
                t['ligandname'],
                t['2d_structure'],
                t['ligandtype'],
                ("No" if t['endogenous'] == None else t['endogenous']),
                t['publications'],
                t['compared'],
            )

        data_table += "</tbody></table>"
        cache.set(cache_key, data_table, 60*60*24*7)

    return data_table

def getLigandCountTable():
    data_table = cache.get("ligand_count_table")
    # data_table = None
    if data_table == None:
        proteins = Protein.objects.filter(sequence_type__slug="wt",
                                          family__slug__startswith="0").prefetch_related(
                                          # species__common_name="Human").prefetch_related(
            "family",
            "family__parent__parent__parent",
            "species"
        )
        # Acquired slugs
        # entry_names = [ p.entry_name for p in proteins ]
        drugtargets_approved = list(Protein.objects.filter(drugs__status="approved").values("entry_name").annotate(num_ligands=Count("drugs__name", distinct=True)))
        # drugtargets_approved = list(Protein.objects.filter(drugs__status="approved").values_list("entry_name", flat=True))
        approved = {}
        for entry in drugtargets_approved:
            approved[entry['entry_name']] = entry['num_ligands']
        drugtargets_trials = list(Protein.objects.filter(drugs__status__in=["in trial"],
                                                         drugs__clinicalstatus__in=["completed", "not open yet",
                                                                                    "ongoing", "recruiting",
                                                                                    "suspended"]).values(
            "entry_name").annotate(num_ligands=Count("drugs__name", distinct=True)))

        trials = {}
        for entry in drugtargets_trials:
            trials[entry['entry_name']] = entry['num_ligands']
        # ligand_set = list(AssayExperiment.objects.values_list("protein__family__slug", "protein__species_id__latin_name")\
        #     .annotate(num_ligands=Count("ligand", distinct=True)))

        ligand_set = list(AssayExperiment.objects.values("protein__entry_name")\
            .annotate(num_ligands=Count("ligand__pk", distinct=True)))

        ligand_count = {}
        for entry in ligand_set:
            ligand_count[entry["protein__entry_name"]] = entry["num_ligands"]

        data_table = "<table id='uniprot_selection' class='uniprot_selection stripe compact'> \
            <thead>\
              <tr> \
                <th colspan=1>&nbsp;</th> \
                <th colspan=6>Receptor classification</th> \
                <th colspan=1>Ligands</th> \
                <th colspan=2>Drugs</th> \
              </tr> \
              <tr> \
                <th><br><br><input autocomplete='off' class='form-check-input' type='checkbox' onclick='return check_all_targets();'></th> \
                <th>Class<br>&nbsp;</th> \
                <th>Ligand type<br>&nbsp;</th> \
                <th style=\"width; 100px;\">Family<br>&nbsp;</th> \
                <th>Species<br>&nbsp;</th> \
                <th style=\"color:red\">Receptor<br>(UniProt)</th> \
                <th style=\"color:red\">Receptor<br>(GtP)</th> \
                <th>Count</th> \
                <th>Approved</th> \
                <th>In clinical<br>trials</th> \
              </tr> \
            </thead>\
            \n \
            <tbody>\n"

        slug_list = []
        #link_setup = "<a target=\"_blank\" href=\"{}\"><span class=\"glyphicon glyphicon-new-window btn-xs\"></span></a>"
        link_setup = "<a target=\"_blank\" href=\"{}\">{}</a>"
        for p in proteins:
            # Do not repeat slugs (only unhuman proteins)
            # if p.family.slug in slug_list:
            if p.entry_name in slug_list:
                continue
            slug_list.append(p.entry_name)
            t = {}
            t['slug'] = p.family.slug
            t['entry_name'] = p.entry_name
            # Ligand count
            t['ligand_count'] = 0
            if t['entry_name'] in ligand_count and ligand_count[t['entry_name']] > 0:
                t['species'] = p.species
                t['id'] = p.id
                #t['ligand_count'] = link_setup.format("/ligand/target/all/" + t['slug'], ligand_count[t['entry_name']])
                t['ligand_count'] = ligand_count[t['entry_name']]
                t['accession'] = p.accession
                t['class'] = p.family.parent.parent.parent.short().split(' ')[0]
                t['ligandtype'] = p.family.parent.parent.short()
                t['family'] = p.family.parent.short()
                t['uniprot'] = p.entry_short()
                t['iuphar'] = p.family.name.replace("receptor", '').strip()

                # Web resource links
                #t['uniprot_link'] = ""
                #t['gtp_link'] = ""
                uniprot_links = p.web_links.filter(web_resource__slug='uniprot')
                if uniprot_links.count() > 0:
                    #t['uniprot_link'] = link_setup.format(p.web_links.filter(web_resource__slug='uniprot')[0])
                    t['uniprot'] = link_setup.format(uniprot_links[0], t['uniprot'])

                gtop_links = p.web_links.filter(web_resource__slug='gtop')
                if gtop_links.count() > 0:
                    #t['gtp_link'] = link_setup.format(p.web_links.filter(web_resource__slug='gtop')[0])
                    t['iuphar'] = link_setup.format(gtop_links[0], t['iuphar'])

                t['approved_target'] = approved[t['entry_name']] if t['entry_name'] in approved.keys() else 0
                t['clinical_target'] = trials[t['entry_name']] if t['entry_name'] in trials.keys() else 0

                data_table += "<tr> \
                <td data-sort=\"0\"><input autocomplete='off' class=\"form-check-input\" type=\"checkbox\" name=\"reference\" data-entry=\"{}\" entry-value=\"{}\"></td> \
                <td>{}</td> \
                <td>{}</td> \
                <td>{}</td> \
                <td>{}</td> \
                <td><span class=\"expand\">{}</span></td> \
                <td><span class=\"expand\">{}</span></td> \
                <td>{}</td> \
                <td>{}</td> \
                <td>{}</td> \
                </tr> \n".format(
                    t['slug'],
                    t['id'],
                    t['class'],
                    t['ligandtype'],
                    t['family'],
                    t['species'],
                    t['uniprot'],
                    t['iuphar'],
                    t['ligand_count'],
                    t['approved_target'],
                    t['clinical_target'],
                )

        data_table += "</tbody></table>"
        cache.set("ligand_count_table", data_table, 60*60*24*7)

    return data_table

def getTargetTable():
    data_table = cache.get("target_table")
    if data_table == None:
        proteins = Protein.objects.filter(sequence_type__slug="wt",
                                          family__slug__startswith="0",
                                          species__common_name="Human").prefetch_related(
            "family",
            "family__parent__parent__parent"
        )
        # Acquired slugs
        slug_list = [ p.family.slug for p in proteins ]

        # Acquire all targets that do not have a human ortholog
        missing_slugs = list(Protein.objects.filter(sequence_type__slug="wt", family__slug__startswith="0")\
                                         .exclude(family__slug__in=slug_list)\
                                         .distinct("family__slug")\
                                         .values_list("family__slug", flat=True))

        for i in missing_slugs:
            missing = Protein.objects.filter(family__slug=i)\
                                        .order_by("id")\
                                        .prefetch_related(
                "family",
                "family__parent__parent__parent"
            )
            proteins = proteins | missing[:1]

        pdbids = list(Structure.objects.all().exclude(structure_type__slug__startswith='af-').values_list("pdb_code__index", "protein_conformation__protein__family_id"))

        allpdbs = {}
        for pdb in pdbids:
            if pdb[1] not in allpdbs:
                allpdbs[pdb[1]] = [pdb[0]]
            else:
                allpdbs[pdb[1]].append(pdb[0])

        drugtargets_approved = list(Protein.objects.filter(drugs__status="approved").values_list("entry_name", flat=True))
        drugtargets_trials = list(Protein.objects.filter(drugs__status__in=["in trial"],
                                                         drugs__clinicalstatus__in=["completed", "not open yet",
                                                                                    "ongoing", "recruiting",
                                                                                    "suspended"]).values_list(
            "entry_name", flat=True))

        ligand_set = list(AssayExperiment.objects.values("protein__family__slug")\
            .annotate(num_ligands=Count("ligand", distinct=True)))

        ligand_count = {}
        for entry in ligand_set:
            ligand_count[entry["protein__family__slug"]] = entry["num_ligands"]

        # Filter data source to Guide to Pharmacology until other coupling transduction sources are "consolidated".
        couplings = ProteinCouplings.objects.filter(source="GuideToPharma").values_list("protein__entry_name",
                                                                                           "g_protein__name",
                                                                                           "transduction")

        signaling_data = {}
        for pairing in couplings:
            if pairing[0] not in signaling_data:
                signaling_data[pairing[0]] = {}
            signaling_data[pairing[0]][pairing[1]] = pairing[2]

        data_table = "<table id='uniprot_selection' class='uniprot_selection stripe compact'> \
            <thead>\
              <tr> \
                <th colspan=1>&nbsp;</th> \
                <th colspan=5>Receptor classification</th> \
                <th colspan=1>Ligands</th> \
                <th colspan=2>Structures</th> \
<!--                <th colspan=2>Drugs</th> -->\
                <th colspan=4>G protein coupling</th> \
              </tr> \
              <tr> \
                <th><br><br><input autocomplete='off' class='form-check-input' type='checkbox' onclick='return check_all_targets();'></th> \
                <th>Class<br>&nbsp;</th> \
                <th>Ligand type<br>&nbsp;</th> \
                <th style=\"width; 100px;\">Family<br>&nbsp;</th> \
                <th style=\"color:red\">Receptor<br>(UniProt)</th> \
                <th style=\"color:red\">Receptor<br>(GtP)</th> \
                <th>Count</th> \
                <th>Count</th> \
                <th>PDB(s)<br>&nbsp;</th> \
<!--                <th>Target of an approved drug</th> \
                <th>Target in clinical trials</th> --> \
                <th>Gs<br>&nbsp;</th> \
                <th>Gi/o<br>&nbsp;</th> \
                <th>Gq/11<br>&nbsp;</th> \
                <th>G12/13<br>&nbsp;</th> \
              </tr> \
            </thead>\
            \n \
            <tbody>\n"

        slug_list = []
        #link_setup = "<a target=\"_blank\" href=\"{}\"><span class=\"glyphicon glyphicon-new-window btn-xs\"></span></a>"
        link_setup = "<a target=\"_blank\" href=\"{}\">{}</a>"
        for p in proteins:
            # Do not repeat slugs (only unhuman proteins)
            if p.family.slug in slug_list:
                continue
            slug_list.append(p.family.slug)
            t = {}
            t['accession'] = p.accession
            t['name'] = p.entry_name.split("_")[0]
            t['slug'] = p.family.slug
            t['class'] = p.family.parent.parent.parent.short().split(' ')[0]
            t['ligandtype'] = p.family.parent.parent.short()
            t['family'] = p.family.parent.short()
            t['uniprot'] = p.entry_short()
            t['iuphar'] = p.family.name.replace("receptor", '').strip()

            # Web resource links
            #t['uniprot_link'] = ""
            #t['gtp_link'] = ""
            uniprot_links = p.web_links.filter(web_resource__slug='uniprot')
            if uniprot_links.count() > 0:
                t['uniprot'] = link_setup.format(uniprot_links[0], t['uniprot'])

            gtop_links = p.web_links.filter(web_resource__slug='gtop')
            if gtop_links.count() > 0:
                t['iuphar'] = link_setup.format(gtop_links[0], t['iuphar'])

            # Ligand count
            t['ligand_count'] = 0
            if t['slug'] in ligand_count:
                t['ligand_count'] = link_setup.format("/ligand/target/all/" + t['slug'], ligand_count[t['slug']])

            t['pdbid'] = t['pdbid_two'] = t['pdbid_tooltip'] = "-"
            t['pdb_count'] = 0
            if p.family_id in allpdbs:
                t['pdb_count'] = len(allpdbs[p.family_id])

                pdb_entries = allpdbs[p.family_id]
                pdb_entries.sort()
                t['pdbid'] = ",".join(pdb_entries)
                t['pdbid_two'] = ",".join(pdb_entries[:2])
                if len(allpdbs[p.family_id]) > 2:
                    t['pdbid_two'] += ",..."
                    n = 4 # Number of PDBs per line
                    pdb_sets = ["&nbsp;&nbsp;".join(pdb_entries[i:i + n]) for i in range(0, len(pdb_entries), n)]
                    t['pdbid_tooltip'] = "<br>".join(pdb_sets)

            t['approved_target'] = "Yes" if p.entry_name in drugtargets_approved else "No"
            t['clinical_target'] = "Yes" if p.entry_name in drugtargets_trials else "No"

            gprotein_families = ["Gs", "Gi/o", "Gq/11", "G12/13"]
            for gprotein in gprotein_families:
                if p.entry_name in signaling_data and gprotein in signaling_data[p.entry_name]:
                    t[gprotein] = signaling_data[p.entry_name][gprotein]
                else:
                    t[gprotein] = "-"

            data_table += "<tr> \
            <td data-sort=\"0\"><input autocomplete='off' class=\"form-check-input\" type=\"checkbox\" name=\"targets\" id=\"{}\" data-entry=\"{}\" data-human=\"{}\"></td> \
            <td>{}</td> \
            <td>{}</td> \
            <td>{}</td> \
            <td><span class=\"expand\">{}</span></td> \
            <td><span class=\"expand\">{}</span></td> \
            <td>{}</td> \
            <td>{}</td> \
            <td><span {} data-html=\"true\" data-placement=\"bottom\" title=\"{}\" data-search=\"{}\" >{}</span></td> \
            <!--<td>{}</td> \
            <td>{}</td>--> \
            <td>{}</td> \
            <td>{}</td> \
            <td>{}</td> \
            <td>{}</td> \
            </tr> \n".format(
                t['slug'],
                t['name'],
                ("No" if t['slug'] in missing_slugs else "Yes"),
                t['class'],
                t['ligandtype'],
                t['family'],
                t['uniprot'],
                t['iuphar'],
                t['ligand_count'],
                t['pdb_count'],
                ("data-toggle=\"tooltip\"" if t['pdbid_tooltip']!="-" else ""),
                t['pdbid_tooltip'],
                t['pdbid'],      # This one hidden used for search box.
                t['pdbid_two'],  # This one shown. Show only first two pdb's.
                t['approved_target'],
                t['clinical_target'],
                t[gprotein_families[0]].capitalize(),
                t[gprotein_families[1]].capitalize(),
                t[gprotein_families[2]].capitalize(),
                t[gprotein_families[3]].capitalize(),
            )

        data_table += "</tbody></table>"
        cache.set("target_table", data_table, 60*60*24*7)

    return data_table

def getReferenceTable(pathway, subtype):
    cache_key = "reference_table_" + pathway + "_" + subtype
    data_table = cache.get(cache_key)
    # data_table = None
    if data_table == None:
        #get all the proteins that are in biaseddata
        biased_proteins = list(BiasedData.objects.values_list("receptor_id__entry_name").distinct())

        biased_entry_names = [b[0] for b in biased_proteins]

        proteins = Protein.objects.filter(
                                          entry_name__in=biased_entry_names,
                                          sequence_type__slug="wt",
                                          family__slug__startswith="0").prefetch_related(
            "family",
            "family__parent__parent__parent",
            "species"
        )

        #Complete data
        totals = list(BiasedData.objects.values("receptor_id").annotate(total=Count("ligand_id", distinct=True)))

        if subtype == 'yes':
            physio_bias = list(BiasedData.objects.filter(subtype_biased__isnull=False).values("receptor_id").annotate(physio=Count("ligand_id", distinct=True)))
            balanced_refs = list(BalancedLigands.objects.filter(subtype_balanced=True).values("receptor_id").annotate(balanced=Count("ligand_id", distinct=True)))
            path_bias = list(BiasedData.objects.filter(pathway_subtype_biased__isnull=False).values("receptor_id").annotate(path=Count("ligand_id", distinct=True)))
        else:
            physio_bias = list(BiasedData.objects.filter(physiology_biased__isnull=False).values("receptor_id").annotate(physio=Count("ligand_id", distinct=True)))
            balanced_refs = list(BalancedLigands.objects.filter(subtype_balanced=False).values("receptor_id").annotate(balanced=Count("ligand_id", distinct=True)))
            path_bias = list(BiasedData.objects.filter(pathway_biased__isnull=False).values("receptor_id").annotate(path=Count("ligand_id", distinct=True)))

        ligand_tot = {}
        for entry in totals:
            if entry['receptor_id'] not in ligand_tot.keys():
                ligand_tot[entry['receptor_id']] = [entry['total'], '-', '-', '-']
        for entry in physio_bias:
            if entry['receptor_id'] in ligand_tot.keys():
                ligand_tot[entry['receptor_id']][1] = entry['physio']
        for entry in balanced_refs:
            if entry['receptor_id'] in ligand_tot.keys():
                ligand_tot[entry['receptor_id']][2] = entry['balanced']
        for entry in path_bias:
            if entry['receptor_id'] in ligand_tot.keys():
                ligand_tot[entry['receptor_id']][3] = entry['path']

        if pathway == "yes":
            data_table = "<table id='uniprot_selection' class='uniprot_selection stripe compact'> \
                <thead>\
                  <tr> \
                    <th colspan=1>&nbsp;</th> \
                    <th colspan=6>Receptor classification</th> \
                    <th colspan=1 style=\"border-left: 1px solid black; text-align:left\">Number of ligands</th> \
                  </tr> \
                  <tr> \
                    <th><br><br></th> \
                    <th>Class<br>&nbsp;</th> \
                    <th>Ligand type<br>&nbsp;</th> \
                    <th style=\"width; 100px;\">Family<br>&nbsp;</th> \
                    <th>Species<br>&nbsp;</th> \
                    <th style=\"color:red\">Receptor<br>(UniProt)</th> \
                    <th style=\"color:red\">Receptor<br>(GtP)</th> \
                    <th>Tested<br>(total)</th> \
                  </tr> \
                </thead>\
                \n \
                <tbody>\n"
        else:
            data_table = "<table id='uniprot_selection' class='uniprot_selection stripe compact'> \
                <thead>\
                  <tr> \
                    <th colspan=1>&nbsp;</th> \
                    <th colspan=6>Receptor classification</th> \
                    <th colspan=4 style=\"border-left: 1px solid black; text-align:left\">Number of ligands</th> \
                  </tr> \
                  <tr> \
                    <th><br><br></th> \
                    <th>Class<br>&nbsp;</th> \
                    <th>Ligand type<br>&nbsp;</th> \
                    <th style=\"width; 100px;\">Family<br>&nbsp;</th> \
                    <th>Species<br>&nbsp;</th> \
                    <th style=\"color:red\">Receptor<br>(UniProt)</th> \
                    <th style=\"color:red\">Receptor<br>(GtP)</th> \
                    <th>Tested<br>(total)</th> \
                    <th>Balanced<br>references</th> \
                    <th>Pathway<br>biased *</th> \
                    <th>Physiology<br>biased *</th> \
                  </tr> \
                </thead>\
                \n \
                <tbody>\n"

        # slug_list = []
        id_list = []
        link_setup = "<a target=\"_blank\" href=\"{}\">{}</a>"

        # NEW CODE
        for p in proteins:
            if p.id in id_list:
                continue
            id_list.append(p.id)
            t = {}
            t['id'] = p.id
            t['accession'] = p.accession
            t['name'] = p.entry_name.split("_")[0]
            t['slug'] = p.family.slug
            t['class'] = p.family.parent.parent.parent.short().split(' ')[0]
            t['ligandtype'] = p.family.parent.parent.short()
            t['species'] = p.species.latin_name
            t['family'] = p.family.parent.short()
            t['uniprot'] = p.entry_short()
            t['iuphar'] = p.family.name.replace("receptor", '').strip()

            uniprot_links = p.web_links.filter(web_resource__slug='uniprot')
            if uniprot_links.count() > 0:
                t['uniprot'] = link_setup.format(uniprot_links[0], t['uniprot'])

            gtop_links = p.web_links.filter(web_resource__slug='gtop')
            if gtop_links.count() > 0:
                t['iuphar'] = link_setup.format(gtop_links[0], t['iuphar'])

            # Ligand count
            t['ligand_count'] = 0

            if t['id'] in ligand_tot:
                t['ligand_count'] = link_setup.format("/ligand/target/all/" + t['slug'], ligand_tot[t['id']][0])
                if ligand_tot[t['id']][1] == '-':
                    t['biased_count'] = '-'
                    t['biased_span'] = ''
                else:
                    t['biased_count'] = link_setup.format("/ligand/target/all/" + t['slug'], ligand_tot[t['id']][1])
                    t['biased_span'] = ligand_tot[t['id']][1]
                if ligand_tot[t['id']][2] == '-':
                    t['balanced_refs'] = '-'
                    t['balanced_span'] = ''
                else:
                    t['balanced_refs'] = link_setup.format("/ligand/target/all/" + t['slug'], ligand_tot[t['id']][2])
                    t['balanced_span'] = ligand_tot[t['id']][2]
                if ligand_tot[t['id']][3] == '-':
                    t['pathway_count'] = '-'
                    t['pathway_span'] = ''
                else:
                    t['pathway_count'] = link_setup.format("/ligand/target/all/" + t['slug'], ligand_tot[t['id']][3])
                    t['pathway_span'] = ligand_tot[t['id']][3]

            if pathway == "yes":
                data_table += "<tr> \
                <td data-sort=\"0\"><input autocomplete='off' class=\"form-check-input\" type=\"checkbox\" name=\"reference\" id=\"{}\" data-entry=\"{}\" entry-value=\"{}\"></td> \
                <td>{}</td> \
                <td>{}</td> \
                <td>{}</td> \
                <td>{}</td> \
                <td><span class=\"expand\">{}</span></td> \
                <td><span class=\"expand\">{}</span></td> \
                <td style=\"border-left: 1px solid black; text-align:left\">{}</td> \
                </tr> \n".format(
                    t['slug'],
                    t['name'],
                    t['id'],
                    t['class'],
                    t['ligandtype'],
                    t['family'],
                    t['species'],
                    t['uniprot'],
                    t['iuphar'],
                    t['ligand_count'],
                )
            else:
                data_table += "<tr> \
                <td data-sort=\"0\"><input autocomplete='off' class=\"form-check-input\" type=\"checkbox\" name=\"reference\" id=\"{}\" data-entry=\"{}\" entry-value=\"{}\"></td> \
                <td>{}</td> \
                <td>{}</td> \
                <td>{}</td> \
                <td>{}</td> \
                <td><span class=\"expand\">{}</span></td> \
                <td><span class=\"expand\">{}</span></td> \
                <td style=\"border-left: 1px solid black; text-align:left\">{}</td> \
                <td data-search=\"{}\">{}</td> \
                <td data-search=\"{}\">{}</td> \
                <td data-search=\"{}\">{}</td> \
                </tr> \n".format(
                    t['slug'],
                    t['name'],
                    t['id'],
                    t['class'],
                    t['ligandtype'],
                    t['family'],
                    t['species'],
                    t['uniprot'],
                    t['iuphar'],
                    t['ligand_count'],
                    t['balanced_span'],
                    t['balanced_refs'],
                    t['pathway_span'],
                    t['pathway_count'],
                    t['biased_span'],
                    t['biased_count'],
                )

        data_table += "</tbody></table>"
        cache.set(cache_key, data_table, 60*60*24*7)

    return data_table

class AbsReferenceSelectionTable(TemplateView):

    """An abstract class for the table reference selection page used in many apps.
    To use it in another app, create a class-based view that extends this class
    """

    template_name = 'common/referenceselectiontable.html'

    type_of_selection = 'reference_table'
    selection_only_receptors = True
    step = 1
    number_of_steps = 2
    title = 'SELECT TARGET'
    description = 'Select target in the table (below) or browse the classification tree (right). You can select entire target' \
        + ' families or individual targets.\n\nOnce you have selected all your targets, click the green button.'
    documentation_url = settings.DOCUMENTATION_URL

    docs = False
    filters = True

    target_input = False
    import_export_box = True
    default_species = 'Human'
    default_slug = '000'
    default_subslug = '0'

    numbering_schemes = False
    search = False
    family_tree = True
    redirect_on_select = True
    filter_gprotein = False
    selection_heading = False
    buttons = {
        'continue': {
            'label': 'Continue to next step',
            'url': '#',
            'color': 'success',
        },
    }
    # OrderedDict to preserve the order of the boxes
    selection_boxes = OrderedDict([
        ('reference', True),
        ('targets', False),
        ('segments', False),
    ])

    # proteins and families
    #try - except block prevents manage.py from crashing - circular dependencies between protein - common
    # try:
    #     if ProteinFamily.objects.filter(slug=default_slug).exists():
    #         ppf = ProteinFamily.objects.get(slug=default_slug)
    #         pfs = ProteinFamily.objects.filter(parent=ppf.id).filter(slug__startswith=default_subslug)
    #         ps = Protein.objects.filter(family=ppf)
    #         psets = ProteinSet.objects.all().prefetch_related('proteins')
    #         tree_indent_level = []
    #         action = 'expand'
    #         # remove the parent family (for all other families than the root of the tree, the parent should be shown)
    #         del ppf

            # Load the target table data
    # table_data = getReferenceTable('no', 'no')
    # except Exception as e:
    #     pass

    # species
    sps = Species.objects.all()

    # numbering schemes
    gns = ResidueNumberingScheme.objects.exclude(slug=settings.DEFAULT_NUMBERING_SCHEME).exclude(slug__in=default_schemes_excluded)

    def get_context_data(self, **kwargs):
        """Get context from parent class

        (really only relevant for children of this class, as TemplateView does
        not have any context variables)
        """

        context = super().get_context_data(**kwargs)

        context["table_data"] = getReferenceTable('no', 'no')
        # get selection from session and add to context
        # get simple selection from session
        simple_selection = self.request.session.get('selection', False)

        # create full selection and import simple selection (if it exists)
        selection = Selection()

        # on the first page of a workflow, clear the selection (or dont' import from the session)
        if self.step is not 1:
            if simple_selection:
                selection.importer(simple_selection)

        # update session
        # receptor = Protein.objects.get(entry_name = simple_selection.reference)
        # context['selection']['receptor_id'] = selection.receptor.id
        simple_selection = selection.exporter()
        self.request.session['selection'] = simple_selection

        context['selection'] = {}
        for selection_box, include in self.selection_boxes.items():
            if include:
                context['selection'][selection_box] = selection.dict(selection_box)['selection'][selection_box]
        # if self.filters:
        #     context['selection']['species'] = selection.species
        #     context['selection']['annotation'] = selection.annotation

        # get attributes of this class and add them to the context
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]
        return context

class AbsTargetSelectionTable(TemplateView):
    """An abstract class for the tablew target selection page used in many apps.

    To use it in another app, create a class-based view that extends this class
    """

    template_name = 'common/targetselectiontable.html'

    type_of_selection = 'targets_table'
    selection_only_receptors = False
    import_export_box = True
    step = 1
    number_of_steps = 2
    title = 'SELECT TARGETS'
    description = 'Select targets in the table (below) or browse the classification tree (right). You can select entire target' \
        + ' families or individual targets.\n\nOnce you have selected all your targets, click the green button.'
    documentation_url = settings.DOCUMENTATION_URL

    docs = False
    filters = True

    target_input = False

    default_species = 'Human'
    default_slug = '000'
    default_subslug = '0'

    numbering_schemes = False
    search = False
    family_tree = True
    redirect_on_select = False
    filter_gprotein = False
    selection_heading = False
    buttons = {
        'continue': {
            'label': 'Continue to next step',
            'url': '#',
            'color': 'success',
        },
    }
    # OrderedDict to preserve the order of the boxes
    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])

    # proteins and families
    #try - except block prevents manage.py from crashing - circular dependencies between protein - common
    try:
        if ProteinFamily.objects.filter(slug=default_slug).exists():
            ppf = ProteinFamily.objects.get(slug=default_slug)
            pfs = ProteinFamily.objects.filter(parent=ppf.id).filter(slug__startswith=default_subslug)
            ps = Protein.objects.filter(family=ppf)
            psets = ProteinSet.objects.all().prefetch_related('proteins')
            tree_indent_level = []
            action = 'expand'
            # remove the parent family (for all other families than the root of the tree, the parent should be shown)
            del ppf

            # Load the target table data
            table_data = getTargetTable()
    except Exception as e:
        pass

    # species
    sps = Species.objects.all()

    # g proteins
    g_prots_slugs = ['100_001_002', '100_001_003', '100_001_001', '100_001_004', '100_001_005']
    gprots = ProteinFamily.objects.filter(slug__in=g_prots_slugs)
    # gprots = ProteinGProtein.objects.all()

    # numbering schemes
    gns = ResidueNumberingScheme.objects.exclude(slug=settings.DEFAULT_NUMBERING_SCHEME).exclude(slug__in=default_schemes_excluded)

    def get_context_data(self, **kwargs):
        """Get context from parent class

        (really only relevant for children of this class, as TemplateView does
        not have any context variables)
        """

        context = super().get_context_data(**kwargs)

        # get selection from session and add to context
        # get simple selection from session
        simple_selection = self.request.session.get('selection', False)

        # create full selection and import simple selection (if it exists)
        selection = Selection()

        # on the first page of a workflow, clear the selection (or dont' import from the session)
        if self.step is not 1:
            if simple_selection:
                selection.importer(simple_selection)

        # default species selection
        if self.default_species:
            sp = Species.objects.get(common_name=self.default_species)
            o = SelectionItem('species', sp)
            selection.species = [o]

        # update session
        simple_selection = selection.exporter()
        self.request.session['selection'] = simple_selection

        context['selection'] = {}
        for selection_box, include in self.selection_boxes.items():
            if include:
                context['selection'][selection_box] = selection.dict(selection_box)['selection'][selection_box]
        if self.filters:
            context['selection']['species'] = selection.species
            context['selection']['annotation'] = selection.annotation
            context['selection']['g_proteins'] = selection.g_proteins
            context['selection']['pref_g_proteins'] = selection.pref_g_proteins

        # get attributes of this class and add them to the context
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]
        return context

class AbsTargetSelection(TemplateView):
    """An abstract class for the target selection page used in many apps.

    To use it in another app, create a class-based view that extends this class.
    """

    template_name = 'common/targetselection.html'

    type_of_selection = 'targets'
    selection_only_receptors = False
    import_export_box = True
    step = 1
    number_of_steps = 2
    title = 'SELECT TARGETS'
    description = 'Select targets by searching or browsing in the middle column. You can select entire target' \
        + ' families or individual targets.\n\nYou can also enter the list of UNIPROT names of the targets (one per line) and click "Add targets" button to add those.\n\nSelected targets will appear in the right column, where you can edit' \
        + ' the list.\n\nOnce you have selected all your targets, click the green button.'
    documentation_url = settings.DOCUMENTATION_URL
    docs = False
    filters = True
    target_input = True
    default_species = 'Human'
    default_slug = '000'
    default_subslug = '0'
    numbering_schemes = False
    search = True
    family_tree = True
    filter_tableselect = False
    redirect_on_select = False
    filter_gprotein = False
    selection_heading = False
    buttons = {
        'continue': {
            'label': 'Continue to next step',
            'url': '#',
            'color': 'success',
        },
    }
    # OrderedDict to preserve the order of the boxes
    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])

    # proteins and families
    #try - except block prevents manage.py from crashing - circular dependencies between protein - common
    try:
        if ProteinFamily.objects.filter(slug=default_slug).exists():
            ppf = ProteinFamily.objects.get(slug=default_slug)
            pfs = ProteinFamily.objects.filter(parent=ppf.id).filter(slug__startswith=default_subslug)
            ps = Protein.objects.filter(family=ppf)
            psets = ProteinSet.objects.all().prefetch_related('proteins')
            tree_indent_level = []
            action = 'expand'
            # remove the parent family (for all other families than the root of the tree, the parent should be shown)
            del ppf
    except Exception as e:
        pass

    # species
    sps = Species.objects.all()

    # g proteins
    g_prots_slugs = ['100_001_002', '100_001_003', '100_001_001', '100_001_004', '100_001_005']
    gprots = ProteinFamily.objects.filter(slug__in=g_prots_slugs)
    # gprots = ProteinGProtein.objects.all()

    # numbering schemes
    gns = ResidueNumberingScheme.objects.exclude(slug=settings.DEFAULT_NUMBERING_SCHEME).exclude(slug__in=default_schemes_excluded)

    def get_context_data(self, **kwargs):
        """get context from parent class

        (really only relevant for children of this class, as TemplateView does
        not have any context variables)
        """

        context = super().get_context_data(**kwargs)

        # get selection from session and add to context
        # get simple selection from session
        simple_selection = self.request.session.get('selection', False)

        # create full selection and import simple selection (if it exists)
        selection = Selection()

        # on the first page of a workflow, clear the selection (or dont' import from the session)
        if self.step is not 1:
            if simple_selection:
                selection.importer(simple_selection)

        # default species selection
        if self.default_species:
            sp = Species.objects.get(common_name=self.default_species)
            o = SelectionItem('species', sp)
            selection.species = [o]

        # update session
        simple_selection = selection.exporter()
        self.request.session['selection'] = simple_selection

        context['selection'] = {}
        for selection_box, include in self.selection_boxes.items():
            if include:
                context['selection'][selection_box] = selection.dict(selection_box)['selection'][selection_box]
        if self.filters:
            context['selection']['species'] = selection.species
            context['selection']['annotation'] = selection.annotation
            context['selection']['g_proteins'] = selection.g_proteins
            context['selection']['pref_g_proteins'] = selection.pref_g_proteins

        # get attributes of this class and add them to the context
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]
        return context


class AbsReferenceSelection(AbsTargetSelection):
    type_of_selection = 'reference'
    step = 1
    number_of_steps = 3
    import_export_box = True
    title = 'SELECT A REFERENCE TARGET'
    description = 'Select a reference target by searching or browsing in the right column.\n\nThe reference will be compared to the targets you select later in the workflow.\n\nOnce you have selected your reference target, you will be redirected to the next step.'
    redirect_on_select = True
    selection_boxes = OrderedDict([
        ('reference', True),
        ('targets', False),
        ('segments', False),
    ])
    psets = [] # protein sets not applicable for this selection


class AbsBrowseSelection(AbsTargetSelection):
    type_of_selection = 'browse'
    step = 1
    number_of_steps = 1
    title = 'SELECT A TARGET OR FAMILY'
    description = 'Select a target or family by searching or browsing in the right column.'
    psets = [] # protein sets not applicable for this selection


class AbsSegmentSelection(TemplateView):
    """An abstract class for the segment selection page used in many apps.

    To use it in another app, create a class-based view for that app that extends this class
    """

    template_name = 'common/segmentselection.html'

    step = 2
    number_of_steps = 2
    title = 'SELECT SEQUENCE SEGMENTS'
    description = 'Select sequence segments in the middle column. You can expand helices and select individual' \
        + ' residues by clicking on the down arrows next to each helix.\n\nSelected segments will appear in the' \
        + ' right column, where you can edit the list.\n\nOnce you have selected all your segments, click the green' \
        + ' button. \n\n'
    documentation_url = settings.DOCUMENTATION_URL
    docs = False
    segment_list = True
    structure_upload = False
    upload_form = PDBform()
    position_type = 'residue'
    buttons = {
        'continue': {
            'label': 'Show alignment',
            'url': '/alignment/render',
            'color': 'success',
        },
    }
    # OrderedDict to preserve the order of the boxes
    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', True),
    ])

    try:
        rsets = ResiduePositionSet.objects.filter(protein_group='gpcr').prefetch_related('residue_position').order_by('set_type','name')
    except Exception as e:
        pass

    ss = ProteinSegment.objects.filter(partial=False, proteinfamily='GPCR').exclude(name__startswith='Fungal').exclude(name__startswith='ECD').prefetch_related('generic_numbers')
    ss_cats = ss.values_list('category').order_by('category').distinct('category')
    action = 'expand'

    amino_acid_groups = definitions.AMINO_ACID_GROUPS
    amino_acid_group_names = definitions.AMINO_ACID_GROUP_NAMES

    # Necessary for the site search functionality
    amino_acid_groups_old = definitions.AMINO_ACID_GROUPS_OLD
    amino_acid_group_names_old = definitions.AMINO_ACID_GROUP_NAMES_OLD

    def get_context_data(self, **kwargs):
        """get context from parent class

        (really only relevant for child classes of this class, as TemplateView
        does not have any context variables).
        """

        context = super().get_context_data(**kwargs)

        # get selection from session and add to context
        # get simple selection from session
        simple_selection = self.request.session.get('selection', False)

        # create full selection and import simple selection (if it exists)
        selection = Selection()
        if simple_selection:
            selection.importer(simple_selection)

        # context['selection'] = selection
        context['selection'] = {}
        context['selection']['site_residue_groups'] = selection.site_residue_groups
        context['selection']['active_site_residue_group'] = selection.active_site_residue_group
        for selection_box, include in self.selection_boxes.items():
            if include:
                context['selection'][selection_box] = selection.dict(selection_box)['selection'][selection_box]

        for f in selection.targets:
            if f.type=='family':
                family = get_gpcr_class(f.item)
                if family.name.startswith('Class D1'):
                    self.ss = ProteinSegment.objects.filter(partial=False, proteinfamily='GPCR').exclude(name__startswith='ECD').prefetch_related('generic_numbers')
                    self.ss_cats = self.ss.values_list('category').order_by('category').distinct('category')
                elif family.name.startswith('Class B'):
                    self.ss = ProteinSegment.objects.filter(partial=False, proteinfamily='GPCR').exclude(name__startswith='Class D1').prefetch_related('generic_numbers')
                    self.ss_cats = self.ss.values_list('category').order_by('category').distinct('category')
            elif f.type=='protein':
                if f.item.family.parent.parent.parent.name.startswith('Class D1'):
                    self.ss = ProteinSegment.objects.filter(partial=False, proteinfamily='GPCR').exclude(name__startswith='ECD').prefetch_related('generic_numbers')
                    self.ss_cats = self.ss.values_list('category').order_by('category').distinct('category')
                elif f.item.family.parent.parent.parent.name.startswith('Class B'):
                    self.ss = ProteinSegment.objects.filter(partial=False, proteinfamily='GPCR').exclude(name__startswith='Class D1').prefetch_related('generic_numbers')
                    self.ss_cats = self.ss.values_list('category').order_by('category').distinct('category')

        # get attributes of this class and add them to the context
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]

        # Set segment selection check if not set
        for button in context["buttons"]:
            if "onclick" not in context["buttons"][button]:
                context["buttons"][button]["onclick"] = 'return VerifyMinSegmentSelection()'

        return context


class AbsMiscSelection(TemplateView):
    """An abstract class for selection pages of other types than target- and segmentselection"""
    template_name = 'common/miscselection.html'
    step = 3
    number_of_steps = 3
    title = ''
    description = ''
    documentation_url = settings.DOCUMENTATION_URL
    docs = False
    buttons = {}
    tree_settings = False
    blast_input = False

    # OrderedDict to preserve the order of the boxes
    selection_boxes = OrderedDict([
        ('targets', True),
        ('segments', True),
    ])
    def get_context_data(self, **kwargs):
        """get context from parent class (really only relevant for child classes of this class, as TemplateView does
        not have any context variables)"""
        context = super().get_context_data(**kwargs)

        # get selection from session and add to context
        # get simple selection from session
        simple_selection = self.request.session.get('selection', False)

        # create full selection and import simple selection (if it exists)
        selection = Selection()

        # on the first page of a workflow, clear the selection (or dont' import from the session)
        if self.step is not 1:
            if simple_selection:
                selection.importer(simple_selection)

        context['selection'] = {}
        context['selection']['tree_settings'] = selection.tree_settings

        for selection_box, include in self.selection_boxes.items():
            if include:
                context['selection'][selection_box] = selection.dict(selection_box)['selection'][selection_box]

        # get attributes of this class and add them to the context
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]
        return context


def AddToSelection(request):
    """Receives a selection request, adds the selected item to session, and returns the updated selection"""
    selection_type = request.GET['selection_type']
    selection_subtype = request.GET['selection_subtype']
    selection_id = request.GET['selection_id']
    # get simple selection from session
    simple_selection = request.session.get('selection', False)

    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)

    # process selected object
    o = []
    if selection_type == 'reference' or selection_type == 'targets':
        # print(selection_type, selection_subtype, selection_id)
        if selection_subtype == 'protein':
            o.append(Protein.objects.get(pk=selection_id))
        if selection_subtype == 'protein_entry':
            o.append(Protein.objects.get(entry_name=selection_id))
            # print("Added {}".format(Protein.objects.get(entry_name=selection_id).name))

        elif selection_subtype == 'protein_set':
            selection_subtype = 'protein'
            pset = ProteinSet.objects.get(pk=selection_id)
            for protein in pset.proteins.all():
                o.append(protein)

        elif selection_subtype == 'family':
            o.append(ProteinFamily.objects.get(pk=selection_id))

        elif selection_subtype == 'set':
            o.append(ProteinSet.objects.get(pk=selection_id))
        # AddToSelection('reference', 'structure',  $(this).children().eq(11).text()+"_refined");
        # 4IAQ
        elif selection_subtype == 'structure':
            o.append(Structure.objects.get(pdb_code__index=selection_id.upper()))

        elif selection_subtype == 'structure_many':
            selection_subtype = 'structure'
            for pdb_code in selection_id.split(","):
                o.append(Structure.objects.get(pdb_code__index=pdb_code.upper()))

        elif selection_subtype == 'signprot_many':
            selection_subtype = 'structure'
            for pdb_code in selection_id.split(","):
                try:
                    o.append(SignprotStructure.objects.get(pdb_code__index=pdb_code.upper()))
                except SignprotStructure.DoesNotExist:
                    o.append(Structure.objects.get(pdb_code__index=pdb_code.upper()))

        elif selection_subtype == 'signprot':
            selection_subtype = 'structure'
            try:
                o.append(SignprotStructure.objects.get(pdb_code__index=selection_id.upper()))
            except SignprotStructure.DoesNotExist:
                o.append(Structure.objects.get(pdb_code__index=selection_id.upper()))

        elif selection_subtype == 'structure_complex_receptor':
            receptor, signprot = selection_id.split('-')
            try:
                scm = StructureComplexModel.objects.get(receptor_protein__entry_name=receptor, sign_protein__entry_name=signprot)
            except StructureComplexModel.DoesNotExist:
                scm = StructureComplexModel.objects.get(receptor_protein__parent__entry_name=selection_id, sign_protein__entry_name=signprot)
            o.append(scm)

        elif selection_subtype == 'structure_complex_signprot':
            o.append(StructureComplexModel.objects.filter(sign_protein__entry_name=selection_id)[0])

        elif selection_subtype == 'structure_model':
            state = selection_id.split('_')[-1]
            entry_name = '_'.join(selection_id.split('_')[:-1])
            try:
                sm = StructureModel.objects.get(protein__entry_name=entry_name.lower(), state__name=state)
            except StructureModel.DoesNotExist:
                sm = StructureModel.objects.get(protein__parent__entry_name=entry_name.lower(), state__name=state)
            o.append(sm)

        elif selection_subtype == 'structure_model_many':
            selection_subtype = 'structure_model'
            for model in selection_id.split(","):
                state = model.split('_')[-1]
                entry_name = '_'.join(model.split('_')[:-1])
                if state == 'refined':
                    o.append(StructureModel.objects.get(protein__entry_name=entry_name.lower()))
                else:
                    o.append(StructureModel.objects.get(protein__entry_name=entry_name.lower(), state__name=state))


    elif selection_type == 'segments':
        if selection_subtype == 'residue':
            o.append(ResidueGenericNumberEquivalent.objects.get(pk=selection_id))
        elif selection_subtype == 'residue_position_set':
            selection_subtype = 'residue'
            rset = ResiduePositionSet.objects.get(pk=selection_id)
            for residue in rset.residue_position.all():
                o.append(residue)
        elif selection_subtype == 'site_residue': # used in site search
            o.append(ResidueGenericNumberEquivalent.objects.get(pk=selection_id))

        else:
            o.append(ProteinSegment.objects.get(pk=selection_id))

    for obj in o:
        # add the selected item to the selection
        selection_object = SelectionItem(selection_subtype, obj)
        selection.add(selection_type, selection_subtype, selection_object)

    # export simple selection that can be serialized
    simple_selection = selection.exporter()
    # add simple selection to session
    request.session['selection'] = simple_selection

    # context
    context = selection.dict(selection_type)

    # template to load
    if selection_subtype == 'site_residue':
        template = 'common/selection_lists_sitesearch.html'
        amino_acid_groups = {
            'amino_acid_groups_old': definitions.AMINO_ACID_GROUPS_OLD,
            'amino_acid_group_names_old': definitions.AMINO_ACID_GROUP_NAMES_OLD}
        context.update(amino_acid_groups)
    else:
        template = 'common/selection_lists.html'

    # amino acid groups
    # print('+++++++++++++++++++++++++')
    # print(context['selection_type'])
    # print(context)
    # for c in context['selection']['targets']:
    #     print(c, c.type, c.item)
    return render(request, template, context)

def RemoveFromSelection(request):
    """Removes one selected item from the session"""
    selection_type = request.GET['selection_type']
    selection_subtype = request.GET['selection_subtype']
    selection_id = request.GET['selection_id']

    # get simple selection from session
    simple_selection = request.session.get('selection', False)

    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)

    # remove the selected item to the selection
    selection.remove(selection_type, selection_subtype, selection_id)
    # export simple selection that can be serialized
    simple_selection = selection.exporter()

    # add simple selection to session
    request.session['selection'] = simple_selection

    # context
    context = selection.dict(selection_type)

    # template to load
    if selection_subtype == 'site_residue':
        template = 'common/selection_lists_sitesearch.html'
        amino_acid_groups = {
            'amino_acid_groups_old': definitions.AMINO_ACID_GROUPS_OLD,
            'amino_acid_group_names_old': definitions.AMINO_ACID_GROUP_NAMES_OLD}
        context.update(amino_acid_groups)
    else:
        template = 'common/selection_lists.html'

    return render(request, template, context)

def ClearSelection(request):
    """Clears all selected items of the selected type from the session"""
    selection_type = request.GET['selection_type']

    # get simple selection from session
    simple_selection = request.session.get('selection', False)

    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)

    # remove the selected item to the selection
    selection.clear(selection_type)

    # export simple selection that can be serialized
    simple_selection = selection.exporter()

    # add simple selection to session
    request.session['selection'] = simple_selection

    return render(request, 'common/selection_lists.html', selection.dict(selection_type))

def CheckSelection(request):
    """Check all selected items of the selected type from the session"""
    selection_type = request.GET['selection_type']

    # get simple selection from session
    simple_selection = request.session.get('selection', False)

    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)

    if selection_type == 'reference':
        tot = len(selection.reference)
    else:
        tot = len(selection.targets)

    return JsonResponse({'total': tot})


def SelectRange(request):
    """Adds generic numbers within the given range"""

    selection_type = request.GET['selection_type']
    selection_subtype = request.GET['selection_subtype']
    range_start = request.GET['range_start']
    range_end = request.GET['range_end']

    # get simple selection from session
    simple_selection = request.session.get('selection', False)

    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)

    # process selected object
    o = []
    if selection_type == 'segments' and selection_subtype == 'residue':
        residue_nums = ResidueGenericNumberEquivalent.objects.all()
        for resn in residue_nums:
            if range_start < float(resn.label.replace('x','.')) < range_end:
                o.append(resn)

def ImportTargetSelection(request):
    """Adds source and proteins into session, returns info to update filters"""
    entry_names = request.GET['entry_names']
    proteins = Protein.objects.filter(entry_name__in=entry_names.split(',')).select_related('species', 'source')
    sources = set([p.source for p in proteins])
    species = set([p.species for p in proteins])
    found_entries = [p.entry_name for p in proteins]

    if len(sources)>1 or list(sources)[0].name!='SWISSPROT':
        source = 'All'
    else:
        source = 'Swissprot'
    species_ids = list(proteins.values_list('species__id', flat=True).distinct())
    slugs = list(proteins.values_list('family__slug', flat=True))
    out = {'slugs':slugs, 'species':species_ids, 'source':source, 'found_entries':found_entries}

    simple_selection = request.session.get('selection', False)

    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)

    for prot in proteins:
        selection_object = SelectionItem('protein', prot)
        selection.add('targets', 'protein', selection_object)

    for sp in species:
        selection_object = SelectionItem('species', sp)
        selection.add('species', 'species', selection_object)

    for ps in sources:
        selection_object = SelectionItem('annotation', ps)
        selection.add('annotation', 'annotation', selection_object)

    # export simple selection that can be serialized
    simple_selection = selection.exporter()

    # add simple selection to session
    request.session['selection'] = simple_selection

    return JsonResponse(json.dumps(out), safe=False)

def SelectFullSequence(request):
    """Adds all segments to the selection"""
    selection_type = request.GET['selection_type']

    # get simple selection from session
    simple_selection = request.session.get('selection', False)

    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)

    # get all segments
    if "protein_type" in request.GET:

        if request.GET['protein_type'] == 'gprotein':
            segmentlist = definitions.G_PROTEIN_SEGMENTS
            pfam = 'Alpha'
        else:
            segmentlist = definitions.ARRESTIN_SEGMENTS
            pfam = 'Arrestin'

        preserved = Case(*[When(slug=pk, then=pos) for pos, pk in enumerate(segmentlist['Full'])])
        segments = ProteinSegment.objects.filter(slug__in=segmentlist['Full'], partial=False, proteinfamily=pfam).order_by(preserved)

    else:
        add_class_b, add_class_d = False, False
        for f in selection.targets:
            if f.type=='family':
                family = get_gpcr_class(f.item)
                if family.name.startswith('Class D1'):
                    add_class_d = True
                elif family.name.startswith('Class B1'):
                    add_class_b = True
            elif f.type=='protein':
                if f.item.family.parent.parent.parent.name.startswith('Class D1'):
                    add_class_d = True
                elif f.item.family.parent.parent.parent.name.startswith('Class B1'):
                    add_class_b = True
        if add_class_d:
            segments = ProteinSegment.objects.filter(partial=False, proteinfamily='GPCR').exclude(name__startswith='ECD')
        elif add_class_b:
            segments = ProteinSegment.objects.filter(partial=False, proteinfamily='GPCR').exclude(name__startswith='Class D1')
        else:
            segments = ProteinSegment.objects.filter(partial=False, proteinfamily='GPCR').exclude(name__startswith='Class D1').exclude(name__startswith='ECD')

    for segment in segments:
        selection_object = SelectionItem(segment.category, segment)
        # add the selected item to the selection
        selection.add(selection_type, segment.category, selection_object)

    # export simple selection that can be serialized
    simple_selection = selection.exporter()

    # add simple selection to session
    request.session['selection'] = simple_selection

    return render(request, 'common/selection_lists.html', selection.dict(selection_type))

def SetTreeSelection(request):
    """Adds all alignable segments to the selection"""
    option_no = request.GET['option_no']
    option_id = request.GET['option_id']
    # get simple selection from session
    simple_selection = request.session.get('selection', False)
    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)
    selection.tree_settings[int(option_no)]=option_id
    simple_selection = selection.exporter()
    # add simple selection to session
    request.session['selection'] = simple_selection
    return render(request, 'common/tree_options.html', selection.dict('tree_settings'))

def SelectAlignableSegments(request):
    """Adds all alignable segments to the selection"""
    selection_type = request.GET['selection_type']

    # get simple selection from session
    simple_selection = request.session.get('selection', False)

    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)

    # get specific segments
    if "protein_type" in request.GET:
        if request.GET['protein_type'] == 'gprotein':
            segmentlist = definitions.G_PROTEIN_SEGMENTS
            pfam = 'Alpha'
        else:
            segmentlist = definitions.ARRESTIN_SEGMENTS
            pfam = 'Arrestin'

        preserved = Case(*[When(slug=pk, then=pos) for pos, pk in enumerate(segmentlist['Structured'])])
        segments = ProteinSegment.objects.filter(slug__in=segmentlist['Structured'], partial=False, proteinfamily=pfam).order_by(preserved)
    else:
        segments = ProteinSegment.objects.filter(partial=False, slug__startswith='TM')

    for segment in segments:
        selection_object = SelectionItem(segment.category, segment)
        # add the selected item to the selection
        selection.add(selection_type, segment.category, selection_object)

    # export simple selection that can be serialized
    simple_selection = selection.exporter()

    # add simple selection to session
    request.session['selection'] = simple_selection

    return render(request, 'common/selection_lists.html', selection.dict(selection_type))

def SelectAlignableResidues(request):
    """Adds all alignable residues to the selection"""
    selection_type = request.GET['selection_type']

    # get simple selection from session
    simple_selection = request.session.get('selection', False)

    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)
    if "protein_type" in request.GET:
        if request.GET['protein_type'] == 'gprotein':
            segmentlist = definitions.G_PROTEIN_SEGMENTS
            pfam = 'Gprotein'
        else:
            segmentlist = definitions.ARRESTIN_SEGMENTS
            pfam = 'Arrestin'

        preserved = Case(*[When(slug=pk, then=pos) for pos, pk in enumerate(segmentlist['Structured'])])
        segments = ProteinSegment.objects.filter(slug__in=segmentlist['Structured'], partial=False, proteinfamily=pfam).order_by(preserved)
    else:
        segments = ProteinSegment.objects.filter(proteinfamily='GPCR').order_by('pk')

    numbering_scheme_slug = 'false'

    # find the relevant numbering scheme (based on target selection)

    cgn = False
    seg_ids_all = []
    numbering_schemes = []
    if numbering_scheme_slug == 'cgn':
        cgn = True
    elif numbering_scheme_slug == 'false':
        if simple_selection and simple_selection.reference:
            if simple_selection.reference[0].type == 'family':
                proteins = Protein.objects.filter(family__slug__startswith=simple_selection.reference[0].item.slug)
                r_prot = proteins[0]
            elif simple_selection.reference[0].type == 'protein':
                r_prot = simple_selection.reference[0].item
            elif simple_selection.reference[0].type == 'structure':
                r_prot = simple_selection.reference[0].item.protein_conformation.protein
            elif simple_selection.reference[0].type == 'signprot':
                r_prot = simple_selection.reference[0].item.protein

            seg_ids_all = get_protein_segment_ids(r_prot, seg_ids_all)
            if r_prot.residue_numbering_scheme not in numbering_schemes:
                numbering_schemes.append(r_prot.residue_numbering_scheme)

        if simple_selection and simple_selection.targets:
            for t in simple_selection.targets:
                if t.type == 'family':
                    proteins = Protein.objects.filter(family__slug__startswith=t.item.slug)
                    t_prot = proteins[0]
                elif t.type == 'protein':
                    t_prot = t.item
                elif t.type == 'structure':
                    t_prot = t.item.protein_conformation.protein
                elif t.type == 'structure_model':
                    t_prot = t.item.protein
                elif t.type == 'signprot':
                    try:
                        t_prot = t.item.protein_conformation.protein
                    except AttributeError:
                        t_prot = t.item.protein

                seg_ids_all = get_protein_segment_ids(t_prot, seg_ids_all)
                if t_prot.residue_numbering_scheme not in numbering_schemes:
                    numbering_schemes.append(t_prot.residue_numbering_scheme)

        # Filter based on reference and target proteins
        filtered_segments = []
        for segment in segments:
            if segment.id in seg_ids_all:
                filtered_segments.append(segment)

        if len(numbering_schemes) == 0 and len(filtered_segments) == 0:
            numbering_schemes.append(ResidueNumberingScheme.objects.get(slug="gpcrdba"))
            filtered_segments = segments

        segments = filtered_segments
    else:
        numbering_schemes = [ResidueNumberingScheme.objects.get(slug=numbering_scheme_slug)]

    for segment in segments:
        if segment.fully_aligned:
            selection_object = SelectionItem(segment.category, segment)
            # add the selected item to the selection
            selection.add(selection_type, segment.category, selection_object)
        else:
            #if not fully aligned

            if ResidueGenericNumberEquivalent.objects.filter(
            default_generic_number__protein_segment=segment,
            scheme__in=numbering_schemes).exists():
                segment.only_aligned_residues = True
                selection_object = SelectionItem(segment.category, segment, properties={'only_aligned_residues':True})
                selection.add(selection_type, segment.category, selection_object)


    simple_selection = selection.exporter()
    # add simple selection to session
    request.session['selection'] = simple_selection

    return render(request, 'common/selection_lists.html', selection.dict('segments'))

def get_protein_segment_ids(protein, seg_ids_all):
    seg_ids = Residue.objects.filter(protein_conformation__protein=protein).order_by('protein_segment__id').distinct('protein_segment__id').values_list('protein_segment', flat=True)
    for s in seg_ids:
        if s not in seg_ids_all:
            seg_ids_all.append(s)
    return seg_ids_all

def ToggleFamilyTreeNode(request):
    """Opens/closes a node in the family selection tree"""
    action = request.GET['action']
    type_of_selection = request.GET['type_of_selection']

    node_id = request.GET['node_id']
    parent_tree_indent_level = int(request.GET['tree_indent_level'])
    tree_indent_level = []
    for i in range(parent_tree_indent_level+1):
        tree_indent_level.append(0)
    parent_tree_indent_level = tree_indent_level[:]
    del parent_tree_indent_level[-1]

    # session
    simple_selection = request.session.get('selection')
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)

    ppf = ProteinFamily.objects.get(pk=node_id)
    if action == 'expand':
        pfs = ProteinFamily.objects.filter(parent=node_id)

        # species filter
        species_list = []
        for species in selection.species:
            species_list.append(species.item)

        # annotation filter
        protein_source_list = []
        for protein_source in selection.annotation:
            protein_source_list.append(protein_source.item)

        # preferred g proteins filter
        pref_g_proteins_list = []
        for g_protein in selection.pref_g_proteins:
            pref_g_proteins_list.append(g_protein.item)


        # g proteins filter
        g_proteins_list = []
        for g_protein in selection.g_proteins:
            g_proteins_list.append(g_protein.item)

        if species_list:
            ps = Protein.objects.order_by('id').filter(family=ppf,
                species__in=(species_list),
                source__in=(protein_source_list)).order_by('source_id', 'id')
        else:
            ps = Protein.objects.order_by('id').filter(family=ppf,
                source__in=(protein_source_list)).order_by('source_id', 'id')

        if pref_g_proteins_list:
            proteins = [x.protein_id for x in ProteinCouplings.objects.filter(g_protein__in=g_proteins_list, transduction='primary')]
            ps = Protein.objects.order_by('id').filter(pk__in=proteins).filter(pk__in=ps)

        if g_proteins_list:
            proteins = [x.protein_id for x in ProteinCouplings.objects.filter(g_protein__in=g_proteins_list)]
            ps = Protein.objects.order_by('id').filter(pk__in=proteins).filter(pk__in=ps)

        # Excluding G protein Alpha subunit protein structure objects, e.g. 3sn6_a
        ps = ps.exclude(accession__isnull=True, family__parent__parent__name='Alpha')

        action = 'collapse'
    else:
        pfs = ps = {}
        action = 'expand'

    return render(request, 'common/selection_tree.html', {
        'action': action,
        'type_of_selection': type_of_selection,
        'ppf': ppf,
        'pfs': pfs,
        'ps': ps,
        'parent_tree_indent_level': parent_tree_indent_level,
        'tree_indent_level': tree_indent_level,
    })

def SelectionAnnotation(request):
    """Updates the selected level of annotation"""
    protein_source = request.GET['protein_source']

    if protein_source == 'All':
        pss = ProteinSource.objects.all()
    else:
        pss = ProteinSource.objects.filter(name=protein_source)

    # get simple selection from session
    simple_selection = request.session.get('selection', False)

    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)

    # reset the annotation selection
    selection.clear('annotation')

    # add the selected items to the selection
    for ps in pss:
        selection_object = SelectionItem('annotation', ps)
        selection.add('annotation', 'annotation', selection_object)

    # export simple selection that can be serialized
    simple_selection = selection.exporter()
    print('annotation',simple_selection)
    # add simple selection to session
    request.session['selection'] = simple_selection

    return render(request, 'common/selection_filters_annotation.html', selection.dict('annotation'))

def SelectionSpeciesPredefined(request):
    """Updates the selected species to predefined sets (Human and all)"""
    species = request.GET['species']

    # get simple selection from session
    simple_selection = request.session.get('selection', False)

    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)
    all_sps = Species.objects.all()
    sps = False
    if species == 'All':
        sps = []
    if species != 'All' and species:
        sps = Species.objects.filter(common_name=species)
    if sps != False:
        # reset the species selection
        selection.clear('species')

        # add the selected items to the selection
        for sp in sps:
            selection_object = SelectionItem('species', sp)
            selection.add('species', 'species', selection_object)

    # export simple selection that can be serialized
    simple_selection = selection.exporter()

    # add simple selection to session
    request.session['selection'] = simple_selection

    # add all species objects to context (for comparison to selected species)
    context = selection.dict('species')
    context['sps'] = all_sps

    return render(request, 'common/selection_filters_species.html', context)

def SelectionSpeciesToggle(request):
    """Updates the selected species arbitrary selections"""
    species_id = request.GET['species_id']

    all_sps = Species.objects.all()
    sps = Species.objects.filter(pk=species_id)
    # get simple selection from session
    simple_selection = request.session.get('selection', False)

    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)

    # add the selected items to the selection
    for sp in sps:
        exists = selection.remove('species', 'species', species_id)
        if not exists:
            selection_object = SelectionItem('species', sp)
            selection.add('species', 'species', selection_object)

    # export simple selection that can be serialized
    simple_selection = selection.exporter()

    # add simple selection to session
    request.session['selection'] = simple_selection

    # add all species objects to context (for comparison to selected species)
    context = selection.dict('species')
    context['sps'] = Species.objects.all()

    return render(request, 'common/selection_filters_species_selector.html', context)

def SelectionGproteinPredefined(request):
    """Updates the selected g proteins to predefined sets (Gi/Go and all)"""
    g_protein = request.GET['g_protein']
    preferred = request.GET['preferred']

    conversion = {
        'Gs': 'Gs',
        'Gi/o': 'Gi/o',
        'Gq/11': 'Gq/11',
        'G12/13': 'G12/13',
        'GPa1 family': 'GPa1 family'
    }
    # get simple selection from session
    simple_selection = request.session.get('selection', False)

    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)

    g_prots_slugs = ['100_001_002', '100_001_003', '100_001_001', '100_001_004', '100_001_005']
    all_gprots = ProteinFamily.objects.filter(slug__in=g_prots_slugs)
    # all_gprots = ProteinGProtein.objects.all()
    gprots = False
    if g_protein == 'All':
        gprots = []
    if g_protein != 'All' and g_protein:
        gprots = ProteinFamily.objects.filter(name=conversion[g_protein])
        # gprots = ProteinGProtein.objects.filter(name=g_protein)

    if gprots != False:
        # reset the g proteins selection
        if preferred == 'true':
            selection.clear('pref_g_proteins')
        else:
            selection.clear('g_proteins')

        # add the selected items to the selection
        for gprot in gprots:
            selection_object = SelectionItem('g_protein', gprot)
            if preferred == 'true':
                selection.add('pref_g_proteins', 'g_protein', selection_object)
            else:
                selection.add('g_proteins', 'g_protein', selection_object)

    # export simple selection that can be serialized
    simple_selection = selection.exporter()

    # add simple selection to session
    request.session['selection'] = simple_selection

    # add all species objects to context (for comparison to selected species)
    if preferred == 'true':
        context = selection.dict('pref_g_proteins')
    else:
        context = selection.dict('g_proteins')
    context['gprots'] = all_gprots
    # context['gprots'] = ProteinGProtein.objects.all()

    if preferred == 'true':
        return render(request, 'common/selection_filters_pref_gproteins.html', context)
    else:
        return render(request, 'common/selection_filters_gproteins.html', context)

def SelectionGproteinToggle(request):
    """Updates the selected g proteins arbitrary selections"""
    g_protein_id = request.GET['g_protein_id']
    preferred = request.GET['preferred']

    conversion = {
        1: 536,
        2: 545,
        3: 532,
        4: 550,
        5: 553
    }
    g_prots_slugs = ['100_001_002', '100_001_003', '100_001_001', '100_001_004', '100_001_005']
    all_gprots = ProteinFamily.objects.filter(slug__in=g_prots_slugs)
    # all_gprots = ProteinGProtein.objects.all()
    gprots = ProteinFamily.objects.filter(pk=conversion[g_protein_id])
    # gprots = ProteinGProtein.objects.filter(pk=g_protein_id)
    # print("'{}'".format(ProteinGProtein.objects.get(pk=g_protein_id).name))

    # get simple selection from session
    simple_selection = request.session.get('selection', False)

    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)

    # add the selected items to the selection
    for gprot in gprots:
        if preferred == 'true':
            # exists = selection.remove('pref_g_proteins', 'g_protein', conversion[g_protein_id])
            exists = selection.remove('pref_g_proteins', 'g_protein', g_protein_id)
        else:
            # exists = selection.remove('g_proteins', 'g_protein', conversion[g_protein_id])
            exists = selection.remove('g_proteins', 'g_protein', g_protein_id)
        if not exists:
            selection_object = SelectionItem('g_protein', gprot)
            if preferred == 'true':
                selection.add('pref_g_proteins', 'g_protein', selection_object)
            else:
                selection.add('g_proteins', 'g_protein', selection_object)

    # export simple selection that can be serialized
    simple_selection = selection.exporter()

    # add simple selection to session
    request.session['selection'] = simple_selection

    # add all species objects to context (for comparison to selected species)
    if preferred == 'true':
        context = selection.dict('pref_g_proteins')
    else:
        context = selection.dict('g_proteins')

    context['gprots'] = all_gprots
    # context['gprots'] = ProteinGProtein.objects.all()

    if preferred == 'true':
        # print(request.session['selection'])
        return render(request, 'common/selection_filters_pref_gproteins_selector.html', context)
    else:
        return render(request, 'common/selection_filters_gproteins_selector.html', context)

def ExpandSegment(request):
    """Expands a segment to show it's generic numbers"""
    segment_id = request.GET['segment_id']
    position_type = request.GET['position_type']
    numbering_scheme_slug = str(request.GET['numbering_scheme'])

    # get simple selection from session
    simple_selection = request.session.get('selection', False)

    # find the relevant numbering scheme (based on target selection)
    cgn = False
    if numbering_scheme_slug == 'cgn':
        cgn = True
    elif numbering_scheme_slug == 'false' and simple_selection:
        first_item = False
        if simple_selection.reference:
            first_item = simple_selection.reference[0]
        elif simple_selection.targets:
            first_item = simple_selection.targets[0]

        if first_item:
            if first_item.type == 'family':
                proteins = Protein.objects.filter(family__slug__startswith=first_item.item.slug)
                numbering_scheme = proteins[0].residue_numbering_scheme
            elif first_item.type == 'protein':
                numbering_scheme = first_item.item.residue_numbering_scheme
            elif first_item.type == 'structure':
                numbering_scheme = first_item.item.protein_conformation.protein.residue_numbering_scheme
        else:
            numbering_scheme = ResidueNumberingScheme.objects.get(slug="gpcrdba")
    elif numbering_scheme_slug:
            numbering_scheme = ResidueNumberingScheme.objects.get(slug=numbering_scheme_slug)
    else:
        numbering_scheme = ResidueNumberingScheme.objects.get(slug="gpcrdba")

    if cgn ==True:
        # fetch the generic numbers for CGN differently
        context = {}
        context['generic_numbers'] = ResidueGenericNumberEquivalent.objects.filter(
            default_generic_number__protein_segment__id=segment_id,
            scheme__slug='cgn').order_by('label')
        context['position_type'] = position_type
        context['scheme'] = ResidueNumberingScheme.objects.filter(slug='cgn')
        context['schemes'] = ResidueNumberingScheme.objects.filter(slug='cgn')
        context['segment_id'] = segment_id
    else:
        # fetch the generic numbers
        context = {}
        context['generic_numbers'] = ResidueGenericNumberEquivalent.objects.filter(
            default_generic_number__protein_segment__id=segment_id,
            scheme=numbering_scheme).order_by('label')
        context['position_type'] = position_type
        context['scheme'] = numbering_scheme
        context['schemes'] = ResidueNumberingScheme.objects.filter(parent__isnull=False)
        context['segment_id'] = segment_id

    return render(request, 'common/segment_generic_numbers.html', context)

def SelectionSchemesPredefined(request):
    """Updates the selected numbering_schemes to predefined sets (GPCRdb and All)"""
    numbering_schemes = request.GET['numbering_schemes']

    # get simple selection from session
    simple_selection = request.session.get('selection', False)

    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)

    all_gns = ResidueNumberingScheme.objects.exclude(slug=settings.DEFAULT_NUMBERING_SCHEME).exclude(slug__in=default_schemes_excluded)
    gns = False
    if numbering_schemes == 'All':
        if len(selection.numbering_schemes) == all_gns.count():
            gns = []
        else:
            gns = all_gns

    # reset the species selection
    selection.clear('numbering_schemes')

    # add the selected items to the selection
    for gn in gns:
        selection_object = SelectionItem('numbering_schemes', gn)
        selection.add('numbering_schemes', 'numbering_schemes', selection_object)

    # export simple selection that can be serialized
    simple_selection = selection.exporter()

    # add simple selection to session
    request.session['selection'] = simple_selection

    # add all species objects to context (for comparison to selected species)
    context = selection.dict('numbering_schemes')
    context['gns'] = all_gns

    return render(request, 'common/selection_filters_numbering_schemes.html', context)

def SelectionSchemesToggle(request):
    """Updates the selected numbering schemes arbitrary selections"""
    numbering_scheme_id = request.GET['numbering_scheme_id']
    gns = ResidueNumberingScheme.objects.filter(pk=numbering_scheme_id).exclude(slug__in=default_schemes_excluded)

    # get simple selection from session
    simple_selection = request.session.get('selection', False)

    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)

    # add the selected items to the selection
    for gn in gns:
        exists = selection.remove('numbering_schemes', 'numbering_schemes', numbering_scheme_id)
        if not exists:
            selection_object = SelectionItem('numbering_schemes', gn)
            selection.add('numbering_schemes', 'numbering_schemes', selection_object)

    # export simple selection that can be serialized
    simple_selection = selection.exporter()

    # add simple selection to session
    request.session['selection'] = simple_selection

    # add all species objects to context (for comparison to selected species)
    context = selection.dict('numbering_schemes')
    context['gns'] = ResidueNumberingScheme.objects.exclude(slug=settings.DEFAULT_NUMBERING_SCHEME).exclude(slug__in=default_schemes_excluded)

    return render(request, 'common/selection_filters_numbering_schemes.html', context)

def UpdateSiteResidueFeatures(request):
    """Updates the selected features of a site residue"""
    selection_type = request.GET['selection_type']
    selection_subtype = request.GET['selection_subtype']
    selection_id = request.GET['selection_id']

    o = []
    if selection_type == 'reference' or selection_type == 'targets':
        if selection_subtype == 'protein':
            o.append(Protein.objects.get(pk=selection_id))
        elif selection_subtype == 'protein_set':
            selection_subtype = 'protein'
            pset = ProteinSet.objects.get(pk=selection_id)
            for protein in pset.proteins.all():
                o.append(protein)
        elif selection_subtype == 'family':
            o.append(ProteinFamily.objects.get(pk=selection_id))
        elif selection_subtype == 'set':
            o.append(ProteinSet.objects.get(pk=selection_id))
        elif selection_subtype == 'structure':
            o.append(Protein.objects.get(entry_name=selection_id.lower()))
    elif selection_type == 'segments':
        if selection_subtype == 'residue':
            o.append(ResidueGenericNumberEquivalent.objects.get(pk=selection_id))
        elif selection_subtype == 'residue_with_properties': # used in site search
            o.append(ResidueGenericNumberEquivalent.objects.get(pk=selection_id))
        else:
            o.append(ProteinSegment.objects.get(pk=selection_id))

    # get simple selection from session
    simple_selection = request.session.get('selection', False)

    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)

    for obj in o:
        # add the selected item to the selection
        selection_object = SelectionItem(selection_subtype, obj)
        selection.add(selection_type, selection_subtype, selection_object)

    # export simple selection that can be serialized
    simple_selection = selection.exporter()

    # add simple selection to session
    request.session['selection'] = simple_selection

    return render(request, 'common/selection_lists.html', selection.dict(selection_type))

def SelectResidueFeature(request):
    """Receives a selection request, add a feature selection to an item, and returns the updated selection"""
    selection_type = request.GET['selection_type']
    selection_subtype = request.GET['selection_subtype']
    selection_id = int(request.GET['selection_id'])
    feature = request.GET['feature']

    # get simple selection from session
    simple_selection = request.session.get('selection', False)

    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)

    # process selected object
    if selection_type == 'segments' and selection_subtype == 'site_residue':
        for obj in selection.segments:
            if int(obj.item.id) == selection_id:
                obj.properties['feature'] = feature
                obj.properties['amino_acids'] = ','.join(definitions.AMINO_ACID_GROUPS_OLD[feature])
                break

    # export simple selection that can be serialized
    simple_selection = selection.exporter()

    # add simple selection to session
    request.session['selection'] = simple_selection

    # context
    context = selection.dict(selection_type)

    # amino acid groups
    amino_acid_groups = {
        'amino_acid_groups_old': definitions.AMINO_ACID_GROUPS_OLD,
        'amino_acid_group_names_old': definitions.AMINO_ACID_GROUP_NAMES_OLD }
    context.update(amino_acid_groups)

    # template to load
    template = 'common/selection_lists_sitesearch.html'

    return render(request, template, context)

def AddResidueGroup(request):
    """Receives a selection request, creates a new residue group, and returns the updated selection"""
    selection_type = request.GET['selection_type']

    # get simple selection from session
    simple_selection = request.session.get('selection', False)

    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)

    # add a group

    selection.site_residue_groups.append([])

    # make the new group active
    selection.active_site_residue_group = len(selection.site_residue_groups)

    # export simple selection that can be serialized
    simple_selection = selection.exporter()

    # add simple selection to session
    request.session['selection'] = simple_selection

    # context
    context = selection.dict(selection_type)

    # amino acid groups
    amino_acid_groups = {
        'amino_acid_groups_old': definitions.AMINO_ACID_GROUPS_OLD,
        'amino_acid_group_names_old': definitions.AMINO_ACID_GROUP_NAMES_OLD,
        'amino_acid_groups': definitions.AMINO_ACID_GROUPS,
        'amino_acid_group_names': definitions.AMINO_ACID_GROUP_NAMES }
    context.update(amino_acid_groups)

    # template to load
    template = 'common/selection_lists_sitesearch.html'

    return render(request, template, context)

def SelectResidueGroup(request):
    """Receives a selection request, updates the active residue group, and returns the updated selection"""
    selection_type = request.GET['selection_type']
    group_id = int(request.GET['group_id'])

    # get simple selection from session
    simple_selection = request.session.get('selection', False)

    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)

    # update the selected group
    selection.active_site_residue_group = group_id

    # export simple selection that can be serialized
    simple_selection = selection.exporter()

    # add simple selection to session
    request.session['selection'] = simple_selection

    # context
    context = selection.dict(selection_type)

    # amino acid groups
    amino_acid_groups = {
        'amino_acid_groups_old': definitions.AMINO_ACID_GROUPS_OLD,
        'amino_acid_group_names_old': definitions.AMINO_ACID_GROUP_NAMES_OLD,
        'amino_acid_groups': definitions.AMINO_ACID_GROUPS,
        'amino_acid_group_names': definitions.AMINO_ACID_GROUP_NAMES }
    context.update(amino_acid_groups)

    # template to load
    template = 'common/selection_lists_sitesearch.html'

    return render(request, template, context)

def RemoveResidueGroup(request):
    """Receives a selection request, removes a residue group, and returns the updated selection"""
    selection_type = request.GET['selection_type']
    group_id = int(request.GET['group_id'])

    # get simple selection from session
    simple_selection = request.session.get('selection', False)

    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)

    # find all positions in the group and delete them
    for position in selection.segments:
        if position.type == 'site_residue':
            if position.properties['site_residue_group'] == group_id:
                selection.remove(selection_type, position.type, position.item.id)
            else:
                if position.properties['site_residue_group'] > group_id:
                    position.properties['site_residue_group'] -= 1

    # export simple selection that can be serialized
    simple_selection = selection.exporter()

    # add simple selection to session
    request.session['selection'] = simple_selection

    # context
    context = selection.dict(selection_type)

    # amino acid groups
    amino_acid_groups = {
        'amino_acid_groups_old': definitions.AMINO_ACID_GROUPS_OLD,
        'amino_acid_group_names_old': definitions.AMINO_ACID_GROUP_NAMES_OLD,
        'amino_acid_groups': definitions.AMINO_ACID_GROUPS,
        'amino_acid_group_names': definitions.AMINO_ACID_GROUP_NAMES }
    context.update(amino_acid_groups)

    # template to load
    template = 'common/selection_lists_sitesearch.html'

    return render(request, template, context)

def SetGroupMinMatch(request):
    """Receives a selection request, sets a minimum match for a group, and returns the updated selection"""
    selection_type = request.GET['selection_type']
    group_id = int(request.GET['group_id'])
    min_match = int(request.GET['min_match'])

    # get simple selection from session
    simple_selection = request.session.get('selection', False)

    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)

    # update the group
    selection.site_residue_groups[group_id - 1][0] = min_match

    # export simple selection that can be serialized
    simple_selection = selection.exporter()

    # add simple selection to session
    request.session['selection'] = simple_selection

    # context
    context = selection.dict(selection_type)

    # amino acid groups
    amino_acid_groups = {
        'amino_acid_groups_old': definitions.AMINO_ACID_GROUPS_OLD,
        'amino_acid_group_names_old': definitions.AMINO_ACID_GROUP_NAMES_OLD,
        'amino_acid_groups': definitions.AMINO_ACID_GROUPS,
        'amino_acid_group_names': definitions.AMINO_ACID_GROUP_NAMES }
    context.update(amino_acid_groups)

    # template to load
    template = 'common/selection_lists_sitesearch.html'

    return render(request, template, context)

def VerifyMinimumSelection(request):
    """Receives a selection request, checks if the minimum # has been selected"""
    selection_type = request.GET['selection_type']
    minimum = int(request.GET['minimum'])

    # get simple selection from session
    simple_selection = request.session.get('selection', False)

    if selection_type == "receptors":
        # Borrow function from alignment to check receptor count
        a = Alignment()
        a.load_proteins_from_selection(simple_selection)

        if len(a.proteins) >= minimum:
            return HttpResponse(True)
        else:
            return HttpResponse(False)
    else:
        return HttpResponse("Not supported", 404)

def ResiduesDownload(request):

    simple_selection = request.session.get('selection', False)

    outstream = BytesIO()
    wb = xlsxwriter.Workbook(outstream, {'in_memory': True})
    worksheet = wb.add_worksheet()
    row_count = 0

    for position in simple_selection.segments:
        if position.type == 'residue':
            worksheet.write_row(row_count, 0, ['residue', position.item.scheme.slug, position.item.label])
            row_count += 1
        elif position.type == 'helix':
            worksheet.write_row(row_count, 0, ['helix', '', position.item.slug])
            row_count += 1
    wb.close()
    outstream.seek(0)
    response = HttpResponse(outstream.read(), content_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
    response['Content-Disposition'] = "attachment; filename=segment_selection.xlsx"

    return response

def ResiduesUpload(request):
    """Receives a file containing generic residue positions along with numbering scheme and adds those to the selection."""

    # get simple selection from session
    simple_selection = request.session.get('selection', False)

    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)

    selection_type = 'segments'
    if request.FILES == {}:
        return render(request, 'common/selection_lists.html', '')

    #Overwriting the existing selection
    selection.clear(selection_type)

    #The excel file
    o = []
    workbook = xlrd.open_workbook(file_contents=request.FILES['xml_file'].read())
    worksheets = workbook.sheet_names()
    for worksheet_name in worksheets:
        worksheet = workbook.sheet_by_name(worksheet_name)
        for row in worksheet.get_rows():
            if len(row) < 2:
                continue
            try:
                if row[0].value == 'residue':
                    position = ResidueGenericNumberEquivalent.objects.get(label=row[2].value, scheme__slug=row[1].value)
                    o.append(position)
                elif row[0].value == 'helix':
                    o.append(ProteinSegment.objects.get(slug=row[2].value))
            except Exception as msg:
                print(msg)
                continue
    for obj in o:
        # add the selected item to the selection
        if obj.__class__.__name__ == 'ResidueGenericNumberEquivalent':
            selection_subtype = 'residue'
        else:
            selection_subtype = 'helix'
        selection_object = SelectionItem(selection_subtype, obj)
        selection.add(selection_type, selection_subtype, selection_object)

    # export simple selection that can be serialized
    simple_selection = selection.exporter()

    # add simple selection to session
    request.session['selection'] = simple_selection

    return render(request, 'common/selection_lists.html', selection.dict(selection_type))

@csrf_exempt
def ReadTargetInput(request):
    """Receives the data from the input form and adds the listed targets to the selection"""

    # get simple selection from session
    simple_selection = request.session.get('selection', False)

    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)

    selection_type = 'targets'
    selection_subtype = 'protein'

    if request.POST == {}:
        return render(request, 'common/selection_lists.html', '')

    # Process input names
    up_names = request.POST['input-targets'].split('\r')
    for up_name in up_names:
        up_name = up_name.strip()
        obj = None
        if "_" in up_name: # Maybe entry name
            selection_subtype = 'protein'
            try:
                obj = Protein.objects.get(entry_name=up_name.lower())
            except:
                obj = None
        else: # Maybe accession code
            selection_subtype = 'protein'
            try:
                obj = Protein.objects.get(accession=up_name.upper())
            except:
                obj = None

        # Try id
        if obj == None and (up_name.isnumeric()):
            selection_subtype = 'protein'
            try:
                obj = up_name
            except:
                obj = None

        # Try slugs
        if obj == None and (up_name.isnumeric() or "_" in up_name):
            selection_subtype = 'family'
            try:
                obj = ProteinFamily.objects.get(slug=up_name)
            except:
                obj = None

        if obj != None:
            selection_object = SelectionItem(selection_subtype, obj)
            selection.add(selection_type, selection_subtype, selection_object)

    # export simple selection that can be serialized
    simple_selection = selection.exporter()
    # add simple selection to session
    request.session['selection'] = simple_selection

    # context
    context = selection.dict(selection_type)

    return render(request, 'common/selection_lists.html', context)

@csrf_exempt
def ReadReferenceInput(request):
    """Receives the data from the input form and adds the selected reference to the selection"""

    # get simple selection from session
    simple_selection = request.session.get('selection', False)
    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)
    selection_type = 'reference'
    selection_subtype = 'protein'

    if request.POST == {}:
        return render(request, 'common/selection_lists.html', '')

    # Process input names
    up_names = request.POST['input-targets'].split('\r')
    for up_name in up_names:
        up_name = up_name.strip()
        obj = None
        if "_" in up_name: # Maybe entry name
            selection_subtype = 'protein'
            try:
                obj = Protein.objects.get(entry_name=up_name.lower())
            except:
                obj = None
        else: # Maybe accession code
            selection_subtype = 'protein'
            try:
                obj = Protein.objects.get(accession=up_name.upper())
            except:
                obj = None

        # # Try id
        if obj == None and (up_name.isnumeric()):
            selection_subtype = 'protein'
            try:
                obj = up_name
            except:
                obj = None

        # Try slugs
        if obj == None and (up_name.isnumeric() or "_" in up_name):
            selection_subtype = 'family'
            try:
                obj = ProteinFamily.objects.get(slug=up_name)
            except:
                obj = None

        if obj != None:
            selection_object = SelectionItem(selection_subtype, obj)
            selection.add(selection_type, selection_subtype, selection_object)

    # export simple selection that can be serialized
    simple_selection = selection.exporter()


    # add simple selection to session
    request.session['selection'] = simple_selection
    # context
    context = selection.dict(selection_type)

    return render(request, 'common/selection_lists.html', context)

@csrf_exempt
def ExportExcelSuggestions(request):
    """Convert json file to excel file"""
    headers = ['reference','review', 'protein', 'mutation_pos', 'generic', 'mutation_from', 'mutation_to',
        'ligand_name', 'ligand_idtype', 'ligand_id', 'ligand_class',
        'exp_type', 'exp_func',  'exp_wt_value',  'exp_wt_unit','exp_mu_effect_sign', 'exp_mu_effect_type', 'exp_mu_effect_value',
        'exp_mu_effect_qual', 'exp_mu_effect_ligand_prop',  'exp_mu_ligand_ref', 'opt_type', 'opt_wt',
        'opt_mu', 'opt_sign', 'opt_percentage', 'opt_qual','opt_agonist', 'added_date'
         ] #'added_by',



    # icl2_end
    # thermo
    # icl2_start
    # nterm
    # cterm
    # icl3_start
    # removals
    # icl3_end

    sheets = ['nterm','icl2_start','icl2_end','icl3_start','icl3_end','cterm','thermo','removals']

    headers = {}
    headers['nterm'] = ['From TM1','hits','Homology levels','Fusions']
    headers['icl2_start'] = ['Start position','hits','Homology levels','PDB codes']
    headers['icl2_end'] = ['End position','hits','Homology levels','PDB codes']
    headers['icl3_start'] = ['Start position','hits','Homology levels','PDB codes']
    headers['icl3_end'] = ['End position','hits','Homology levels','PDB codes']
    headers['cterm'] = ['From TM8','hits','Homology levels']
    headers['thermo'] = ['Generic position','Mutation','hits','methods']
    headers['removals'] = ['segment','mutation','type','subtype']

    data = request.POST['d']
    data = json.loads(data)

    #EXCEL SOLUTION
    output = BytesIO()
    workbook = xlsxwriter.Workbook(output)

    for name in sheets:
        values = data[name]
        if len(values)==0:
            continue
        worksheet = workbook.add_worksheet(name)
        col = 0
        for h in headers[name]:
            worksheet.write(0, col, h)
            col += 1
        row = 1
        for d in values:
            col = 0
            #print(d)
            # for h in headers:
            #     worksheet.write(row, col, d[h])
            #     col += 1
            for c in d:
                #print(c)
                if isinstance(c, list):
                    c = ",".join(c)
                worksheet.write(row, col, c)
                col += 1
            row += 1
    workbook.close()
    output.seek(0)
    xlsx_data = output.read()

    ts = time.time()

    cache.set(ts,xlsx_data,30)
    response = HttpResponse(ts)
    return response

@csrf_exempt
def ExportExcelModifications(request):
    """Convert json file to excel file"""
    headers = ['#','type', 'method', 'range', 'info','insert_location','order','from','to','sequence','fixed','extra']

    #EXCEL SOLUTION
    output = BytesIO()
    workbook = xlsxwriter.Workbook(output)

    if 'd' in request.POST:
        data = request.POST['d']
        data = json.loads(data)

        worksheet = workbook.add_worksheet("modifications")
        row = 1
        index = {}
        col = 0
        for h in headers:
            worksheet.write(0, col, h)
            index[h] = col
            col += 1
        number = 0
        for mod in data:
            worksheet.write(row, 0, str(number))
            number += 1
            for m,v in mod.items():
                if isinstance(v, list) and m=='range' and len(v)>1:
                    #v = ",".join(str(x) for x in v)
                    v = str(v[0])+"-"+str(v[-1]) #first and last
                elif isinstance(v, list):
                    v = ",".join(str(x) for x in v)
                if m in index:
                    worksheet.write(row, index[m], str(v))
                else:
                    print('No column for '+m)
            row += 1
    elif 'm' in request.POST:
        datas = request.POST['m']
        datas = json.loads(datas)
        worksheet = workbook.add_worksheet("modifications")
        row = 0
        for i, data in enumerate(datas):

            worksheet.write(row, 0, "Construct Number #"+str(i+1))
            row += 1
            index = {}
            col = 0
            for h in headers:
                worksheet.write(row, col, h)
                index[h] = col
                col += 1
            number = 0
            row += 1
            for mod in data:
                if len(mod)>0:
                    worksheet.write(row, 0, str(number))
                    number += 1
                    for m,v in mod.items():
                        if isinstance(v, list) and m=='range' and len(v)>1:
                            #v = ",".join(str(x) for x in v)
                            v = str(v[0])+"-"+str(v[-1]) #first and last
                        elif isinstance(v, list):
                            v = ",".join(str(x) for x in v)
                        if m in index:
                            worksheet.write(row, index[m], str(v))
                        else:
                            print('No column for '+m)
                row += 1
            row += 1


    worksheet2 = workbook.add_worksheet("FASTA SEQUENCES")
    # worksheet2.write(0, 0, "sequences")
    row = 0

    sequences = request.POST['s']
    sequences = json.loads(sequences)
    for s in sequences:
        # worksheet2.write(row, 0, "Modifications # used")
        # worksheet2.write(row, 1, ' '.join(s[0]))
        # row += 1
        worksheet2.write(row, 0, ">" + s[1] + "[Modifications:" + ' '.join(s[0]) +"]")
        row += 1
        #worksheet2.write(row, 0, "Sequence")
        worksheet2.write(row, 0, s[2])
        row += 1

    worksheet3 = workbook.add_worksheet("FASTA SEQUENCES BLOCK")
    # worksheet2.write(0, 0, "sequences")
    row = 0

    sequences = request.POST['s']
    sequences = json.loads(sequences)
    for s in sequences:
        # worksheet2.write(row, 0, "Modifications # used")
        # worksheet2.write(row, 1, ' '.join(s[0]))
        # row += 1
        worksheet3.write(row, 0, ">" + s[1] + "[Modifications:" + ' '.join(s[0]) +"]")
        row += 1
        #worksheet2.write(row, 0, "Sequence")
        worksheet3.write(row, 0, s[3])
        row += 1

    workbook.close()
    output.seek(0)
    xlsx_data = output.read()

    ts = time.time()

    cache.set(ts,xlsx_data,30)
    response = HttpResponse(ts)
    return response

def ExportExcelDownload(request, ts, entry_name):
    """Convert json file to excel file"""

    response = HttpResponse(cache.get(ts),content_type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')
    response['Content-Disposition'] = 'attachment; filename='+entry_name+'.xlsx' #% 'mutations'

    return response

@csrf_exempt
def ImportExcel(request, **response_kwargs):
    """Receives excel, outputs json"""
    o = []
    if request.method == 'POST':
        form = FileUploadForm(data=request.POST, files=request.FILES)
        if form.is_valid():
            workbook = xlrd.open_workbook(file_contents=request.FILES['file_source'].read())
            worksheets = workbook.sheet_names()
            for worksheet_name in worksheets:
                worksheet = workbook.sheet_by_name(worksheet_name)
                # Only load modifications
                if worksheet_name!='modifications':
                    continue
                num_rows = worksheet.nrows - 1
                num_cells = worksheet.ncols - 1
                curr_row = 0 #skip first, otherwise -1
                while curr_row < num_rows:
                    curr_row += 1
                    #row = worksheet.row(curr_row)
                    curr_cell = -1
                    temprow = []
                    if worksheet.cell_value(curr_row, 0) == '': #if empty
                        continue
                    while curr_cell < num_cells:
                        curr_cell += 1
                        #cell_type = worksheet.cell_type(curr_row, curr_cell)
                        cell_value = worksheet.cell_value(curr_row, curr_cell)
                        temprow.append(cell_value)
                    o.append(temprow)
            jsondata = json.dumps(o)
            response_kwargs['content_type'] = 'application/json'
            return HttpResponse(jsondata, **response_kwargs)

        else:
            pass
    jsondata = json.dumps(o)
    response_kwargs['content_type'] = 'application/json'
    return HttpResponse(jsondata, **response_kwargs)

@csrf_exempt
def ConvertSVG(request, **response_kwargs):
    """Receives SVG, outputs PDF"""

    pdf_content = ""
    filename = "converted.pdf"
    if request.method == 'POST':
        # Obtain the POST data
        filename = request.POST.get('filename')
        svg_content = request.POST.get('dataUrl')

        # unescape the SVG string
        svg_content = html.unescape(svg_content)

        ## Removing inline tspans as super/subscript breaks the layout
        clean = re.compile('</?tspan.*?>')
        svg_content = re.sub(clean, '', svg_content)

        # Build XML tree
        svg = etree.fromstring(svg_content)

        # Render in PDF
        renderer = SvgRenderer("")
        drawing = renderer.render(svg)
        pdf_content = renderPDF.drawToString(drawing)

    response = HttpResponse(content_type="application/pdf")
    response['Content-Disposition'] = 'attachment; filename="' + filename + '"'
    response.write(pdf_content)
    return response

def get_gpcr_class(item):
    while item.parent.parent!=None:
        item = item.parent
    return item

@csrf_exempt
@cache_page(60*60*24*7)
def TargetTableData(request):
    """
    Creates a table for selection of targets.

    The goal is to to offer an alternative to the togglefamilytreenode
    alternative already present in the selection logic.
    Here we do server-side rendering of the table. This is a common trick to optimize response to client.
    See the following relevant video from google chrome.
    <https://youtu.be/ff4fgQxPaO0>
    """

    return HttpResponse(getTargetTable())
