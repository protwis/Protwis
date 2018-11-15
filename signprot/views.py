from django.shortcuts import get_object_or_404, render
from django.http import HttpResponse
from django.http import JsonResponse
from django.core.cache import cache
from django.views.decorators.cache import cache_page

from protein.models import Protein, ProteinConformation, ProteinAlias, ProteinFamily, Gene, ProteinGProtein, ProteinGProteinPair
from residue.models import Residue, ResiduePositionSet

from structure.models import Structure
from mutation.models import MutationExperiment
from common.selection import Selection
from common.diagrams_gpcr import DrawSnakePlot
from common.diagrams_gprotein import DrawGproteinPlot
from common.diagrams_arrestin import DrawArrestinPlot

from signprot.models import SignprotStructure, SignprotBarcode, SignprotInteractions

from common import definitions
from collections import OrderedDict
from common.views import AbsTargetSelection

import json
# Create your views here.
class BrowseSelection(AbsTargetSelection):
    step = 1
    number_of_steps = 1
    psets = False
    filters = True
    filter_gprotein = True

    type_of_selection = 'browse_gprot'

    description = 'Select a G protein or family by searching or browsing in the right column.'
    description = 'Select a G protein (family) by searching or browsing in the middle. The selection is viewed to' \
        + ' the right.'
    docs = 'receptors.html'
    target_input=False

    selection_boxes = OrderedDict([
        ('reference', False), ('targets', True),
        ('segments', False),
    ])
    try:
        ppf_g = ProteinFamily.objects.get(slug="100_000")
        # ppf_a = ProteinFamily.objects.get(slug="200_000")
        # pfs = ProteinFamily.objects.filter(parent__in=[ppf_g.id,ppf_a.id])
        pfs = ProteinFamily.objects.filter(parent__in=[ppf_g.id])
        ps = Protein.objects.filter(family__in=[ppf_g]) # ,ppf_a
        tree_indent_level = []
        # action = 'expand'
        # remove the parent family (for all other families than the root of the tree, the parent should be shown)
        # del ppf_g
        # del ppf_a
    except Exception as e:
        pass

@cache_page(60*60*24*2) # 2 days caching
def GProtein(request):

    name_of_cache = 'gprotein_statistics'

    context = cache.get(name_of_cache)

    if context==None:

        context = OrderedDict()
        i=0

        gproteins = ProteinGProtein.objects.all().prefetch_related('proteingproteinpair_set')
        slugs = ['001','002','004','005']
        slug_translate = {'001':"ClassA", '002':"ClassB1",'004':"ClassC", '005':"ClassF"}
        selectivitydata = {}
        for slug in slugs:
            jsondata = {}
            for gp in gproteins:
                # ps = gp.proteingproteinpair_set.all()
                ps = gp.proteingproteinpair_set.filter(protein__family__slug__startswith=slug)

                if ps:
                    jsondata[str(gp)] = []
                    for p in ps:
                        if str(p.protein.entry_name).split('_')[0].upper() not in selectivitydata:
                            selectivitydata[str(p.protein.entry_name).split('_')[0].upper()] = []
                        selectivitydata[str(p.protein.entry_name).split('_')[0].upper()].append(str(gp))
                        # print(p.protein.family.parent.parent.parent)
                        jsondata[str(gp)].append(str(p.protein.entry_name)+'\n')

                    jsondata[str(gp)] = ''.join(jsondata[str(gp)])

            context[slug_translate[slug]] = jsondata

        context["selectivitydata"] = selectivitydata



    return render(request, 'signprot/gprotein.html', context)

@cache_page(60*60*24*2)
def familyDetail(request, slug):
    # get family
    pf = ProteinFamily.objects.get(slug=slug)

    # get family list
    ppf = pf
    families = [ppf.name]
    while ppf.parent.parent:
        families.append(ppf.parent.name)
        ppf = ppf.parent
    families.reverse()

    # number of proteins
    proteins = Protein.objects.filter(family__slug__startswith=pf.slug, sequence_type__slug='wt')
    no_of_proteins = proteins.count()
    no_of_human_proteins = Protein.objects.filter(family__slug__startswith=pf.slug, species__id=1,
        sequence_type__slug='wt').count()
    list_proteins = list(proteins.values_list('pk',flat=True))

    # get structures of this family
    structures = SignprotStructure.objects.filter(origin__family__slug__startswith=slug
        )

    mutations = MutationExperiment.objects.filter(protein__in=proteins).prefetch_related('residue__generic_number',
                                'exp_qual', 'ligand')

    mutations_list = {}
    for mutation in mutations:
        if not mutation.residue.generic_number: continue #cant map those without display numbers
        if mutation.residue.generic_number.label not in mutations_list: mutations_list[mutation.residue.generic_number.label] = []
        if mutation.ligand:
            ligand = mutation.ligand.name
        else:
            ligand = ''
        if mutation.exp_qual:
            qual = mutation.exp_qual.qual
        else:
            qual = ''
        mutations_list[mutation.residue.generic_number.label].append([mutation.foldchange,ligand.replace("'", "\\'"),qual])

    # Update to consensus sequence in protein confirmation!
    try:
        pc = ProteinConformation.objects.filter(protein__family__slug=slug, protein__sequence_type__slug='consensus')
    except ProteinConformation.DoesNotExist:
        pc = ProteinConformation.objects.get(protein__family__slug=slug, protein__species_id=1,
            protein__sequence_type__slug='wt')

    residues = Residue.objects.filter(protein_conformation=pc).order_by('sequence_number').prefetch_related(
        'protein_segment', 'generic_number', 'display_generic_number')

    jsondata = {}
    jsondata_interaction = {}
    for r in residues:
        if r.generic_number:
            if r.generic_number.label in mutations_list:
                jsondata[r.sequence_number] = [mutations_list[r.generic_number.label]]
            if r.generic_number.label in interaction_list:
                jsondata_interaction[r.sequence_number] = interaction_list[r.generic_number.label]

    # process residues and return them in chunks of 10
    # this is done for easier scaling on smaller screens
    chunk_size = 10
    r_chunks = []
    r_buffer = []
    last_segment = False
    border = False
    title_cell_skip = 0
    for i, r in enumerate(residues):
        # title of segment to be written out for the first residue in each segment
        segment_title = False

        # keep track of last residues segment (for marking borders)
        if r.protein_segment.slug != last_segment:
            last_segment = r.protein_segment.slug
            border = True

        # if on a border, is there room to write out the title? If not, write title in next chunk
        if i == 0 or (border and len(last_segment) <= (chunk_size - i % chunk_size)):
            segment_title = True
            border = False
            title_cell_skip += len(last_segment) # skip cells following title (which has colspan > 1)

        if i and i % chunk_size == 0:
            r_chunks.append(r_buffer)
            r_buffer = []

        r_buffer.append((r, segment_title, title_cell_skip))

        # update cell skip counter
        if title_cell_skip > 0:
            title_cell_skip -= 1
    if r_buffer:
        r_chunks.append(r_buffer)

    context = {'pf': pf, 'families': families, 'structures': structures, 'no_of_proteins': no_of_proteins,
        'no_of_human_proteins': no_of_human_proteins, 'mutations':mutations, 'r_chunks': r_chunks, 'chunk_size': chunk_size}

    return render(request, 'signprot/family_details.html', context)

class TargetSelection(AbsTargetSelection):
    step = 1
    number_of_steps = 1
    filters = False
    psets = False
    target_input = False
    redirect_on_select = True
    type_of_selection = 'ginterface'
    title = 'SELECT TARGET for Gs INTERFACE'
    description = 'Select a reference target by searching or browsing.' \
        + '\n\nThe Gs interface from adrb2 (PDB: 3SN6) will be superposed onto the selected target.' \
        + '\n\nAn interaction browser for the adrb2 Gs interface will be given for comparison"'

    # template_name = 'common/targetselection.html'

    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])

    buttons = {
        'continue': {
            'label': 'Continue to next step',
            'url': '#',
            'color': 'success',
        },
    }

@cache_page(60*60*24*2)
def Ginterface(request, protein = None):

    residuelist = Residue.objects.filter(protein_conformation__protein__entry_name=protein).prefetch_related('protein_segment','display_generic_number','generic_number')
    SnakePlot = DrawSnakePlot(
                residuelist, "Class A (Rhodopsin)", protein, nobuttons=1)

    # TEST
    gprotein_residues = Residue.objects.filter(protein_conformation__protein__entry_name='gnaz_human').prefetch_related('protein_segment','display_generic_number','generic_number')
    gproteinplot = DrawGproteinPlot(
                gprotein_residues, "Gprotein", protein)

    crystal = Structure.objects.get(pdb_code__index="3SN6")
    aa_names = definitions.AMINO_ACID_GROUP_NAMES_OLD
    names_aa = dict(zip(aa_names.values(),aa_names.keys()))
    names_aa['Polar (S/T)'] = 'pol_short'
    names_aa['Polar (N/Q/H)'] = 'pol_long'

    residues_browser = [{'pos': 135, 'aa': 'I', 'gprotseg': "H5",'segment': 'TM3', 'ligand': 'Gs', 'type': aa_names['hp'], 'gpcrdb': '3.54x54', 'gpnum': 'G.H5.16', 'gpaa': 'Q384', 'availability': 'interacting'},{'pos': 136, 'aa': 'T', 'gprotseg': "H5",'segment': 'TM3', 'ligand': 'Gs', 'type': 'Polar (S/T)', 'gpcrdb': '3.55x55', 'gpnum': 'G.H5.12', 'gpaa': 'R380', 'availability': 'interacting'},{'pos': 139, 'aa': 'F', 'gprotseg': "H5",'segment': 'ICL2', 'ligand': 'Gs', 'type': 'Aromatic', 'gpcrdb': '34.51x51', 'gpnum': 'G.H5.8', 'gpaa': 'F376', 'availability': 'interacting'},{'pos': 139, 'aa': 'F', 'gprotseg': "S1",'segment': 'ICL2', 'ligand': 'Gs', 'type': 'Aromatic', 'gpcrdb': '34.51x51', 'gpnum': 'G.S1.2', 'gpaa': 'H41', 'availability': 'interacting'},{'pos': 141, 'aa': 'Y', 'gprotseg': "H5",'segment': 'ICL2', 'ligand': 'Gs', 'type': 'Aromatic', 'gpcrdb': '34.53x53', 'gpnum': 'G.H5.19', 'gpaa': 'H387', 'availability': 'interacting'},{'pos': 225, 'aa': 'E', 'gprotseg': "H5",'segment': 'TM5', 'ligand': 'Gs', 'type': 'Negative charge', 'gpcrdb': '5.64x64', 'gpnum': 'G.H5.12', 'gpaa': 'R380', 'availability': 'interacting'},{'pos': 225, 'aa': 'E', 'gprotseg': "H5",'segment': 'TM5', 'ligand': 'Gs', 'type': 'Negative charge', 'gpcrdb': '5.64x64', 'gpnum': 'G.H5.16', 'gpaa': 'Q384', 'availability': 'interacting'},{'pos': 229, 'aa': 'Q', 'gprotseg': "H5",'segment': 'TM5', 'ligand': 'Gs', 'type': 'Polar (N/Q/H)', 'gpcrdb': '5.68x68', 'gpnum': 'G.H5.13', 'gpaa': 'D381', 'availability': 'interacting'},{'pos': 229, 'aa': 'Q', 'gprotseg': "H5",'segment': 'TM5', 'ligand': 'Gs', 'type': 'Polar (N/Q/H)', 'gpcrdb': '5.68x68', 'gpnum': 'G.H5.16', 'gpaa': 'Q384', 'availability': 'interacting'},{'pos': 229, 'aa': 'Q', 'gprotseg': "H5",'segment': 'TM5', 'ligand': 'Gs', 'type': 'Polar (N/Q/H)', 'gpcrdb': '5.68x68', 'gpnum': 'G.H5.17', 'gpaa': 'R385', 'availability': 'interacting'},{'pos': 274, 'aa': 'T', 'gprotseg': "H5",'segment': 'TM6', 'ligand': 'Gs', 'type': 'Polar (S/T)', 'gpcrdb': '6.36x36', 'gpnum': 'G.H5.24', 'gpaa': 'E392', 'availability': 'interacting'},{'pos': 328, 'aa': 'R', 'gprotseg': "H5",'segment': 'TM7', 'ligand': 'Gs', 'type': 'Positive charge', 'gpcrdb': '7.55x55', 'gpnum': 'G.H5.24', 'gpaa': 'E392', 'availability': 'interacting'}, {'pos': 232, 'aa': 'K', 'segment': 'TM5', 'ligand': 'Gs', 'type': 'Positive charge', 'gpcrdb': '5.71x71', 'gprotseg': "H5", 'gpnum': 'G.H5.13', 'gpaa': 'D381', 'availability': 'interacting'}]

    # accessible_gn = ['3.50x50', '3.53x53', '3.54x54', '3.55x55', '34.50x50', '34.51x51', '34.53x53', '34.54x54', '5.61x61', '5.64x64', '5.65x65', '5.67x67', '5.68x68', '5.71x71', '5.72x72', '5.74x74', '5.75x75', '6.29x29', '6.32x32', '6.33x33', '6.36x36', '6.37x37', '7.55x55', '8.48x48', '8.49x49']

    accessible_gn = ['3.50x50', '3.53x53', '3.54x54', '3.55x55', '3.56x56', '34.50x50', '34.51x51', '34.52x52', '34.53x53', '34.54x54', '34.55x55', '34.56x56', '34.57x57', '5.61x61', '5.64x64', '5.65x65', '5.66x66', '5.67x67', '5.68x68', '5.69x69', '5.71x71', '5.72x72', '5.74x74', '5.75x75', '6.25x25', '6.26x26', '6.28x28', '6.29x29', '6.32x32', '6.33x33', '6.36x36', '6.37x37', '6.40x40', '7.55x55', '7.56x56', '8.47x47', '8.48x48', '8.49x49', '8.51x51']

    exchange_table = OrderedDict([('hp', ('V','I', 'L', 'M')),
                                 ('ar', ('F', 'H', 'W', 'Y')),
                                 ('pol_short', ('S', 'T')), # Short/hydroxy
                                 ('pol_long', ('N', 'Q', 'H')), # Amino-like (both donor and acceptor
                                 ('neg', ('D', 'E')),
                                 ('pos', ('K', 'R'))])

    interacting_gn = []

    accessible_pos = list(residuelist.filter(display_generic_number__label__in=accessible_gn).values_list('sequence_number', flat=True))

    # Which of the Gs interacting_pos are conserved?
    GS_none_equivalent_interacting_pos = []
    GS_none_equivalent_interacting_gn = []

    for interaction in residues_browser:
        interacting_gn.append(interaction['gpcrdb'])
        gs_b2_interaction_type_long = (next((item['type'] for item in residues_browser if item['gpcrdb'] == interaction['gpcrdb']), None))

        interacting_aa = residuelist.filter(display_generic_number__label__in=[interaction['gpcrdb']]).values_list('amino_acid', flat=True)

        if interacting_aa:
            interaction['aa'] = interacting_aa[0]
            pos = residuelist.filter(display_generic_number__label__in=[interaction['gpcrdb']]).values_list('sequence_number', flat=True)[0]
            interaction['pos'] = pos

            feature = names_aa[gs_b2_interaction_type_long]

            if interacting_aa[0] not in exchange_table[feature]:
                GS_none_equivalent_interacting_pos.append(pos)
                GS_none_equivalent_interacting_gn.append(interaction['gpcrdb'])

    GS_equivalent_interacting_pos = list(residuelist.filter(display_generic_number__label__in=interacting_gn).values_list('sequence_number', flat=True))

    gProteinData = ProteinGProteinPair.objects.filter(protein__entry_name=protein)

    primary = []
    secondary = []

    for entry in gProteinData:
        if entry.transduction == 'primary':
            primary.append((entry.g_protein.name.replace("Gs","G<sub>s</sub>").replace("Gi","G<sub>i</sub>").replace("Go","G<sub>o</sub>").replace("G11","G<sub>11</sub>").replace("G12","G<sub>12</sub>").replace("G13","G<sub>13</sub>").replace("Gq","G<sub>q</sub>").replace("G","G&alpha;"),entry.g_protein.slug))
        elif entry.transduction == 'secondary':
            secondary.append((entry.g_protein.name.replace("Gs","G<sub>s</sub>").replace("Gi","G<sub>i</sub>").replace("Go","G<sub>o</sub>").replace("G11","G<sub>11</sub>").replace("G12","G<sub>12</sub>").replace("G13","G<sub>13</sub>").replace("Gq","G<sub>q</sub>").replace("G","G&alpha;"),entry.g_protein.slug))


    return render(request, 'signprot/ginterface.html', {'pdbname': '3SN6', 'snakeplot': SnakePlot, 'gproteinplot': gproteinplot, 'crystal': crystal, 'interacting_equivalent': GS_equivalent_interacting_pos, 'interacting_none_equivalent': GS_none_equivalent_interacting_pos, 'accessible': accessible_pos, 'residues': residues_browser, 'mapped_protein': protein, 'interacting_gn': GS_none_equivalent_interacting_gn, 'primary_Gprotein': set(primary), 'secondary_Gprotein': set(secondary)} )

def ajaxInterface(request, slug, **response_kwargs):

    name_of_cache = 'ajaxInterface_'+slug

    jsondata = cache.get(name_of_cache)

    if jsondata == None:

        if slug == "arrs_human":
            rsets = ResiduePositionSet.objects.get(name="Arrestin interface")
        else:
            rsets = ResiduePositionSet.objects.get(name="Gprotein Barcode")
        # residues = Residue.objects.filter(protein_conformation__protein__entry_name=slug, display_generic_number__label=residue.label)

        jsondata = {}
        positions = []
        for x, residue in enumerate(rsets.residue_position.all()):
            try:
                pos = str(list(Residue.objects.filter(protein_conformation__protein__entry_name=slug, display_generic_number__label=residue.label))[0])
            except:
                print("Protein has no residue position at", residue.label)
            a = pos[1:]

            jsondata[a] = [5, 'Receptor interface position', residue.label]

        jsondata = json.dumps(jsondata)

    cache.set(name_of_cache, jsondata, 60*60*24*2) #two days timeout on cache

    response_kwargs['content_type'] = 'application/json'

    return HttpResponse(jsondata, **response_kwargs)

def ajaxBarcode(request, slug, cutoff, **response_kwargs):

    name_of_cache = 'ajaxBarcode_'+slug+cutoff

    jsondata = cache.get(name_of_cache)

    if jsondata == None:
        jsondata = {}

        selectivity_pos = list(SignprotBarcode.objects.filter(protein__entry_name=slug, seq_identity__gte=cutoff).values_list('residue__display_generic_number__label', flat=True))

        conserved = list(SignprotBarcode.objects.filter(protein__entry_name=slug, paralog_score__gte=cutoff, seq_identity__gte=cutoff).prefetch_related('residue__display_generic_number').values_list('residue__display_generic_number__label', flat=True))

        na_data = list(SignprotBarcode.objects.filter(protein__entry_name=slug, seq_identity=0, paralog_score=0).values_list('residue__display_generic_number__label', flat=True))

        all_positions = Residue.objects.filter(protein_conformation__protein__entry_name=slug).prefetch_related('display_generic_number')

        for res in all_positions:
            cgn = str(res.generic_number)
            res = str(res.sequence_number)
            if cgn in conserved:
                jsondata[res] = [0, 'Conserved', cgn]
            elif cgn in selectivity_pos and cgn not in conserved:
                jsondata[res] = [1, 'Selectivity determining', cgn]
            elif cgn in na_data:
                jsondata[res] = [3, 'NA', cgn]
            else:
                jsondata[res] = [2, 'Evolutionary neutral', cgn]

        jsondata = json.dumps(jsondata)
        response_kwargs['content_type'] = 'application/json'

        cache.set(name_of_cache, jsondata, 60*60*24*2) #two days timeout on cache

    return HttpResponse(jsondata, **response_kwargs)

@cache_page(60*60*24*2)
def StructureInfo(request, pdbname):
    """
    Show structure details
    """
    protein = Protein.objects.get(signprotstructure__PDB_code=pdbname)

    crystal = SignprotStructure.objects.get(PDB_code=pdbname)

    return render(request,'signprot/structure_info.html',{'pdbname': pdbname, 'protein': protein, 'crystal': crystal})

# @cache_page(60*60*24*2)
def signprotdetail(request, slug):
    # get protein

    slug = slug.lower()
    p = Protein.objects.prefetch_related('web_links__web_resource').get(entry_name=slug, sequence_type__slug='wt')

    # get family list
    pf = p.family
    families = [pf.name]
    while pf.parent.parent:
        families.append(pf.parent.name)
        pf = pf.parent
    families.reverse()

    # get protein aliases
    aliases = ProteinAlias.objects.filter(protein=p).values_list('name', flat=True)

    # get genes
    genes = Gene.objects.filter(proteins=p).values_list('name', flat=True)
    gene = genes[0]
    alt_genes = genes[1:]

    # get structures of this signal protein
    structures = SignprotStructure.objects.filter(origin=p)

    # mutations
    mutations = MutationExperiment.objects.filter(protein=p)


    # get residues
    pc = ProteinConformation.objects.get(protein=p)

    residues = Residue.objects.filter(protein_conformation=pc).order_by('sequence_number').prefetch_related(
        'protein_segment', 'generic_number', 'display_generic_number')

    # process residues and return them in chunks of 10
    # this is done for easier scaling on smaller screens
    chunk_size = 10
    r_chunks = []
    r_buffer = []
    last_segment = False
    border = False
    title_cell_skip = 0
    for i, r in enumerate(residues):
        # title of segment to be written out for the first residue in each segment
        segment_title = False

        # keep track of last residues segment (for marking borders)
        if r.protein_segment.slug != last_segment:
            last_segment = r.protein_segment.slug
            border = True

        # if on a border, is there room to write out the title? If not, write title in next chunk
        if i == 0 or (border and len(last_segment) <= (chunk_size - i % chunk_size)):
            segment_title = True
            border = False
            title_cell_skip += len(last_segment) # skip cells following title (which has colspan > 1)

        if i and i % chunk_size == 0:
            r_chunks.append(r_buffer)
            r_buffer = []

        r_buffer.append((r, segment_title, title_cell_skip))

        # update cell skip counter
        if title_cell_skip > 0:
            title_cell_skip -= 1
    if r_buffer:
        r_chunks.append(r_buffer)
    context = {'p': p, 'families': families, 'r_chunks': r_chunks, 'chunk_size': chunk_size, 'aliases': aliases,
        'gene': gene, 'alt_genes': alt_genes, 'structures': structures, 'mutations': mutations}

    return render(request, 'signprot/signprot_details.html', context)

def InteractionMatrix(request):
    from django.db.models import F

    dataset = {
        '3sn6' : [
        ['R','F',139,'34.51x51','A','V',217, ["hydrophobic", "van-der-waals"]],
        ['R','Q',229,'5.68x68','A','R',385, ["hydrophobic", "polar-backbone-sidechain", "polar-sidechain-sidechain", "h-bond"]],
        ['R','R',63,'12.49x49','B','D',312, ["polar-backbone-sidechain"]],
        ['R','E',225,'5.64x64','A','Q',384, ["polar-sidechain-sidechain", "van-der-waals"]],
        ['R','E',62,'12.48x48','B','D',312, ["polar-sidechain-sidechain"]],
        ['R','R',239,'-','A','D',343, ["polar-sidechain-backbone"]],
        ['R','F',139,'34.51x51','A','I',383, ["hydrophobic"]],
        ['R','R',239,'-','A','L',346, ["hydrophobic"]],
        ['R','S',143,'34.55x55','A','A',39, ["hydrophobic", "van-der-waals"]],
        ['R','F',139,'34.51x51','A','C',379, ["hydrophobic"]],
        ['R','I',233,'5.72x72','A','L',394, ["hydrophobic"]],
        ['R','L',275,'6.37x37','A','L',393, ["hydrophobic", "van-der-waals"]],
        ['R','Q',229,'5.68x68','A','L',388, ["hydrophobic"]],
        ['R','I',135,'3.54x54','A','L',388, ["hydrophobic", "van-der-waals"]],
        ['R','T',274,'6.36x36','A','L',393, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
        ['R','I',135,'3.54x54','A','Q',384, ["polar-backbone-sidechain", "van-der-waals"]],
        ['R','E',225,'5.64x64','A','R',380, ["polar-sidechain-sidechain"]],
        ['R','K',232,'5.71x71','A','D',381, ["hydrophobic", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
        ['R','I',233,'5.72x72','A','R',385, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
        ['R','Q',229,'5.68x68','A','D',381, ["polar-sidechain-backbone", "polar-sidechain-sidechain"]],
        ['R','R',239,'-','A','T',350, ["hydrophobic", "polar-sidechain-sidechain"]],
        ['R','I',135,'3.54x54','A','Y',391, ["hydrophobic"]],
        ['R','K',235,'5.74x74','A','D',323, ["polar-backbone-sidechain"]],
        ['R','A',271,'6.33x33','A','L',393, ["hydrophobic", "van-der-waals"]],
        ['R','A',226,'5.65x65','A','L',388, ["hydrophobic", "van-der-waals"]],
        ['R','I',135,'3.54x54','A','H',387, ["hydrophobic", "polar-backbone-sidechain"]],
        ['R','F',139,'34.51x51','A','R',380, ["hydrophobic", "van-der-waals"]],
        ['R','T',274,'6.36x36','A','Y',391, ["polar-sidechain-backbone"]],
        ['R','Q',142,'34.54x54','A','I',383, ["hydrophobic"]],
        ['R','R',228,'5.67x67','A','D',381, ["polar-sidechain-sidechain"]],
        ['R','P',138,'34.50x50','A','R',380, ["hydrophobic"]],
        ['R','T',136,'3.55x55','A','R',380, ["polar-backbone-sidechain", "van-der-waals"]],
        ['R','F',139,'34.51x51','A','F',376, ["edge-to-face", "face-to-edge", "hydrophobic", "van-der-waals"]],
        ['R','L',230,'5.69x69','A','L',394, ["hydrophobic", "van-der-waals"]],
        ['R','P',138,'34.50x50','A','I',383, ["hydrophobic", "van-der-waals"]],
        ['R','R',228,'5.67x67','A','Q',384, ["polar-sidechain-sidechain"]],
        ['R','I',233,'5.72x72','A','Y',358, ["hydrophobic", "van-der-waals"]],
        ['R','I',135,'3.54x54','A','L',393, ["hydrophobic", "van-der-waals"]],
        ['R','Q',142,'34.54x54','A','H',387, ["polar-sidechain-sidechain"]],
        ['R','P',138,'34.50x50','A','H',387, ["hydrophobic"]],
        ['R','Y',141,'34.53x53','A','H',387, ["edge-to-face", "face-to-edge", "pi-cation", "hydrophobic"]],
        ['R','D',130,'3.49x49','A','Y',391, ["polar-sidechain-sidechain"]],
        ['R','Q',229,'5.68x68','A','Q',384, ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
        ['R','R',131,'3.50x50','A','Y',391, ["cation-pi", "hydrophobic", "van-der-waals"]],
        ['R','A',134,'3.53x53','A','H',387, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
        ['R','T',274,'6.36x36','A','E',392, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
        ['R','R',239,'-','A','R',347, ["hydrophobic", "polar-sidechain-backbone"]],
        ['R','K',140,'34.52x52','A','R',380, ["polar-sidechain-sidechain"]],
        ['R','V',222,'5.61x61','A','L',393, ["hydrophobic", "van-der-waals"]],
        ['R','F',139,'34.51x51','A','H',41, ["edge-to-face", "face-to-edge", "pi-cation", "hydrophobic"]],
        ['R','P',138,'34.50x50','A','Q',384, ["hydrophobic", "van-der-waals"]],
        ],
        '4x1h' : [
        ['A','A',246,'6.29x29','C','F',350, ["hydrophobic", "van-der-waals"]],
        ['A','V',139,'3.54x54','C','V',340, ["hydrophobic", "water-mediated"]],
        ['A','V',138,'3.53x53','C','D',343, ["hydrophobic", "polar-backbone-sidechain", "water-mediated"]],
        ['A','T',242,'6.25x25','C','F',350, ["hydrophobic", "van-der-waals"]],
        ['A','M',309,'7.56x56','C','G',348, ["water-mediated", "water-mediated"]],
        ['A','I',305,'7.52x52','C','G',348, ["water-mediated"]],
        ['A','Q',312,'8.49x49','C','S',346, ["water-mediated"]],
        ['A','K',141,'3.56x56','C','V',340, ["water-mediated"]],
        ['A','N',73,'2.40x40','C','G',348, ["water-mediated"]],
        ['A','T',229,'5.64x64','C','V',340, ["hydrophobic"]],
        ['A','E',249,'6.32x32','C','F',350, ["hydrophobic"]],
        ['A','V',250,'6.33x33','C','L',349, ["hydrophobic"]],
        ['A','A',233,'5.68x68','C','V',340, ["hydrophobic"]],
        ['A','N',73,'2.40x40','C','C',347, ["water-mediated", "water-mediated"]],
        ['A','V',139,'3.54x54','C','L',344, ["hydrophobic"]],
        ['A','L',72,'2.39x39','C','C',347, ["water-mediated", "water-mediated"]],
        ['A','T',70,'2.37x37','C','S',346, ["water-mediated"]],
        ['A','A',246,'6.29x29','C','L',344, ["hydrophobic"]],
        ['A','N',73,'2.40x40','C','S',346, ["water-mediated", "water-mediated"]],
        ['A','L',72,'2.39x39','C','S',346, ["hydrophobic", "van-der-waals", "water-mediated", "water-mediated"]],
        ['A','K',245,'6.28x28','C','F',350, ["hydrophobic"]],
        ['A','V',230,'5.65x65','C','L',344, ["hydrophobic"]],
        ['A','N',310,'8.47x47','C','S',346, ["water-mediated"]],
        ['A','A',246,'6.29x29','C','L',341, ["hydrophobic", "van-der-waals"]],
        ['A','N',310,'8.47x47','C','C',347, ["polar-sidechain-backbone", "water-mediated"]],
        ['A','K',311,'8.48x48','C','F',350, ["polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
        ['A','T',242,'6.25x25','C','L',341, ["hydrophobic"]],
        ['A','M',253,'6.36x36','C','G',348, ["water-mediated"]],
        ['A','R',135,'3.50x50','C','C',347, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals", "water-mediated", "water-mediated"]],
        ['A','R',135,'3.50x50','C','S',346, ["water-mediated"]],
        ['A','R',135,'3.50x50','C','L',349, ["hydrophobic", "van-der-waals"]],
        ['A','K',141,'3.56x56','C','D',343, ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals", "water-mediated"]],
        ['A','N',310,'8.47x47','C','G',348, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals", "water-mediated", "water-mediated", "water-mediated", "water-mediated"]],
        ['A','L',226,'5.61x61','C','L',344, ["hydrophobic"]],
        ['A','V',139,'3.54x54','C','D',343, ["hydrophobic", "water-mediated"]],
        ['A','V',250,'6.33x33','C','L',344, ["hydrophobic"]],
        ['A','L',226,'5.61x61','C','L',349, ["hydrophobic"]],
        ['A','V',138,'3.53x53','C','V',340, ["water-mediated"]],
        ['A','T',243,'6.26x26','C','L',341, ["hydrophobic"]],
        ['A','P',71,'2.38x38','C','S',346, ["water-mediated"]],
        ['A','Y',306,'7.53x53','C','G',348, ["water-mediated", "water-mediated"]],
        ['A','Q',312,'8.49x49','C','C',347, ["water-mediated"]],
        ['A','Q',312,'8.49x49','C','G',348, ["water-mediated"]],
        ['A','R',135,'3.50x50','C','G',348, ["water-mediated", "water-mediated"]],
        ],
        '5g53' : [
        ['A','K',227,'6.29x29','C','L',394, ["hydrophobic", "polar-sidechain-sidechain"]],
        ['A','L',235,'6.37x37','C','L',393, ["hydrophobic", "van-der-waals"]],
        ['A','Q',207,'5.68x68','C','L',388, ["hydrophobic"]],
        ['A','M',211,'5.72x72','C','Y',358, ["hydrophobic", "polar-sidechain-sidechain"]],
        ['A','Q',207,'5.68x68','C','Q',384, ["polar-sidechain-sidechain", "h-bond"]],
        ['A','A',203,'5.64x64','C','Q',384, ["hydrophobic", "polar-backbone-sidechain"]],
        ['A','L',110,'34.51x51','C','V',217, ["hydrophobic", "van-der-waals"]],
        ['A','I',106,'3.54x54','C','Q',384, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
        ['A','I',200,'5.61x61','C','L',388, ["hydrophobic"]],
        ['A','P',109,'34.50x50','C','R',380, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
        ['A','R',293,'8.48x48','C','E',392, ["hydrophobic", "polar-backbone-sidechain"]],
        ['A','Q',207,'5.68x68','C','Y',360, ["polar-sidechain-sidechain"]],
        ['A','I',200,'5.61x61','C','L',393, ["hydrophobic"]],
        ['A','L',208,'5.69x69','C','L',394, ["hydrophobic", "van-der-waals"]],
        ['A','A',204,'5.65x65','C','L',388, ["hydrophobic"]],
        ['A','L',110,'34.51x51','C','R',380, ["hydrophobic"]],
        ['A','P',109,'34.50x50','C','Q',384, ["hydrophobic"]],
        ['A','Q',207,'5.68x68','C','R',385, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
        ['A','Q',207,'5.68x68','C','D',381, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
        ['A','R',111,'34.52x52','C','R',380, ["polar-sidechain-sidechain"]],
        ['A','I',106,'3.54x54','C','L',388, ["hydrophobic", "van-der-waals"]],
        ['A','Y',112,'34.53x53','C','H',387, ["edge-to-face", "face-to-edge", "pi-cation", "hydrophobic"]],
        ['A','L',110,'34.51x51','C','I',383, ["hydrophobic"]],
        ['A','I',106,'3.54x54','C','Y',391, ["hydrophobic"]],
        ['A','R',107,'3.55x55','C','R',380, ["polar-backbone-sidechain", "van-der-waals"]],
        ['A','Q',210,'5.71x71','C','D',381, ["polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
        ['A','L',110,'34.51x51','C','C',379, ["hydrophobic"]],
        ['A','L',110,'34.51x51','C','H',41, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
        ['A','A',203,'5.64x64','C','L',388, ["hydrophobic"]],
        ['A','L',110,'34.51x51','C','F',219, ["hydrophobic"]],
        ['A','R',291,'7.56x56','C','E',392, ["polar-backbone-sidechain"]],
        ['A','R',102,'3.50x50','C','Y',391, ["cation-pi", "hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
        ['A','I',108,'3.56x56','C','R',380, ["polar-backbone-sidechain"]],
        ['A','A',231,'6.33x33','C','L',393, ["hydrophobic"]],
        ['A','L',110,'34.51x51','C','F',376, ["hydrophobic", "van-der-waals"]],
        ['A','P',109,'34.50x50','C','I',383, ["hydrophobic", "van-der-waals"]],
        ['A','R',296,'8.51x51','C','E',392, ["polar-sidechain-sidechain", "van-der-waals"]],
        ['A','R',111,'34.52x52','C','V',217, ["hydrophobic", "van-der-waals"]],
        ['A','R',291,'7.56x56','C','Y',391, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
        ['A','I',106,'3.54x54','C','H',387, ["hydrophobic"]],
        ['A','A',105,'3.53x53','C','H',387, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
        ['A','R',111,'34.52x52','C','D',215, ["polar-sidechain-backbone", "van-der-waals"]],
        ['B','R',111,'34.52x52','D','D',215, ["polar-sidechain-backbone", "van-der-waals"]],
        ['B','I',106,'3.54x54','D','L',388, ["hydrophobic", "van-der-waals"]],
        ['B','I',106,'3.54x54','D','Y',391, ["hydrophobic"]],
        ['B','L',110,'34.51x51','D','I',383, ["hydrophobic"]],
        ['B','L',110,'34.51x51','D','H',41, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
        ['B','L',110,'34.51x51','D','F',376, ["hydrophobic", "van-der-waals"]],
        ['B','P',109,'34.50x50','D','Q',384, ["hydrophobic"]],
        ['B','G',114,'34.55x55','D','A',39, ["hydrophobic", "van-der-waals"]],
        ['B','R',296,'8.51x51','D','E',392, ["polar-sidechain-sidechain", "van-der-waals"]],
        ['B','Y',112,'34.53x53','D','H',387, ["edge-to-face", "face-to-edge", "pi-cation", "hydrophobic"]],
        ['B','L',110,'34.51x51','D','R',380, ["hydrophobic"]],
        ['B','A',231,'6.33x33','D','L',393, ["hydrophobic", "van-der-waals"]],
        ['B','I',106,'3.54x54','D','Q',384, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
        ['B','Q',207,'5.68x68','D','L',388, ["hydrophobic"]],
        ['B','A',204,'5.65x65','D','L',393, ["hydrophobic", "van-der-waals"]],
        ['B','A',203,'5.64x64','D','Q',384, ["hydrophobic", "polar-backbone-sidechain"]],
        ['B','R',102,'3.50x50','D','Y',391, ["hydrophobic", "polar-backbone-sidechain"]],
        ['B','A',105,'3.53x53','D','H',387, ["hydrophobic", "polar-backbone-sidechain"]],
        ['B','I',200,'5.61x61','D','L',393, ["hydrophobic", "van-der-waals"]],
        ['B','A',203,'5.64x64','D','L',388, ["hydrophobic"]],
        ['B','Q',207,'5.68x68','D','R',385, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
        ['B','L',110,'34.51x51','D','V',217, ["hydrophobic", "van-der-waals"]],
        ['B','L',235,'6.37x37','D','L',393, ["hydrophobic", "van-der-waals"]],
        ['B','A',204,'5.65x65','D','L',388, ["hydrophobic"]],
        ['B','P',109,'34.50x50','D','I',383, ["hydrophobic", "van-der-waals"]],
        ['B','A',105,'3.53x53','D','Y',391, ["van-der-waals"]],
        ['B','I',106,'3.54x54','D','H',387, ["hydrophobic"]],
        ['B','R',111,'34.52x52','D','V',217, ["hydrophobic", "van-der-waals"]],
        ['B','I',108,'3.56x56','D','R',380, ["polar-backbone-sidechain"]],
        ['B','R',291,'7.56x56','D','E',392, ["polar-backbone-sidechain"]],
        ['B','L',110,'34.51x51','D','F',219, ["hydrophobic"]],
        ['B','R',111,'34.52x52','D','R',380, ["polar-sidechain-sidechain"]],
        ['B','P',109,'34.50x50','D','R',380, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
        ['B','N',113,'34.54x54','D','A',39, ["hydrophobic", "van-der-waals"]],
        ['B','L',110,'34.51x51','D','C',379, ["hydrophobic"]],
        ['B','R',293,'8.48x48','D','E',392, ["polar-backbone-sidechain", "van-der-waals"]],
        ['B','R',107,'3.55x55','D','R',380, ["polar-backbone-sidechain", "van-der-waals"]],
        ['B','I',200,'5.61x61','D','L',388, ["hydrophobic"]],
        ['B','Q',207,'5.68x68','D','Q',384, ["polar-sidechain-sidechain"]],
        ['B','Q',207,'5.68x68','D','D',381, ["polar-sidechain-backbone"]],
        ],
        '5uz7' : [
        ['R','Y',243,'3.53x53','A','Y',391, ["van-der-waals"]],
        ['R','N',396,'8.48x48','A','E',392, ["polar-backbone-sidechain"]],
        ['R','T',330,'-','A','Y',358, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
        ['R','R',180,'2.46x46','A','Q',390, ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
        ['R','L',348,'6.45x45','A','L',393, ["hydrophobic"]],
        ['R','Q',415,'8.67x67','B','V',307, ["polar-sidechain-backbone", "van-der-waals"]],
        ['R','F',253,'-','A','V',217, ["hydrophobic"]],
        ['R','T',254,'-','A','H',387, ["polar-sidechain-sidechain"]],
        ['R','N',396,'8.48x48','A','R',356, ["polar-sidechain-sidechain"]],
        ['R','F',253,'-','A','H',41, ["polar-backbone-sidechain", "van-der-waals"]],
        ['R','V',249,'3.59x59','A','Q',384, ["polar-backbone-sidechain"]],
        ['R','R',180,'2.46x46','A','Y',391, ["cation-pi", "hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
        ['R','I',248,'3.58x58','A','Q',384, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
        ['R','H',184,'2.50x50','A','Y',391, ["hydrophobic"]],
        ['R','Q',415,'8.67x67','B','Q',44, ["hydrophobic", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
        ['R','V',252,'-','A','R',380, ["hydrophobic", "van-der-waals"]],
        ['R','K',326,'5.64x64','A','R',385, ["polar-backbone-sidechain"]],
        ['R','Q',408,'8.60x60','B','H',311, ["polar-sidechain-backbone", "van-der-waals"]],
        ['R','I',411,'8.63x63','B','V',307, ["hydrophobic"]],
        ['R','I',248,'3.58x58','A','H',387, ["hydrophobic"]],
        ['R','V',252,'-','A','I',383, ["hydrophobic", "van-der-waals"]],
        ['R','H',331,'-','A','Y',358, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
        ['R','K',326,'5.64x64','A','R',380, ["polar-sidechain-sidechain"]],
        ['R','C',394,'7.60x60','A','E',392, ["polar-backbone-sidechain"]],
        ['R','K',326,'5.64x64','A','Q',384, ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
        ['R','L',244,'3.54x54','A','Y',391, ["hydrophobic"]],
        ['R','L',247,'3.57x57','A','H',387, ["hydrophobic", "van-der-waals"]],
        ['R','R',404,'8.56x56','B','D',312, ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
        ['R','V',252,'-','A','Q',384, ["hydrophobic", "van-der-waals"]],
        ['R','T',345,'6.42x42','A','L',393, ["hydrophobic", "van-der-waals"]],
        ['R','E',329,'-','A','R',385, ["polar-backbone-sidechain"]],
        ['R','I',248,'3.58x58','A','L',388, ["hydrophobic"]],
        ['R','L',323,'5.61x61','A','L',388, ["hydrophobic", "van-der-waals"]],
        ['R','M',327,'5.65x65','A','L',394, ["hydrophobic", "van-der-waals"]],
        ['R','Q',408,'8.60x60','B','A',309, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
        ['R','L',247,'3.57x57','A','Y',391, ["hydrophobic"]],
        ['R','L',348,'6.45x45','A','E',392, ["hydrophobic"]],
        ['R','Q',408,'8.60x60','B','G',310, ["polar-sidechain-backbone"]],
        ],
        '5vai' : [
        ['R','Y',402,'7.57x57','A','E',392, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
        ['R','K',415,'8.56x56','B','D',291, ["polar-sidechain-backbone"]],
        ['R','N',338,'-','A','Y',360, ["hydrophobic", "van-der-waals"]],
        ['R','L',255,'3.58x58','A','R',380, ["polar-backbone-sidechain"]],
        ['R','R',419,'8.60x60','B','N',293, ["polar-sidechain-sidechain"]],
        ['R','R',419,'8.60x60','B','H',311, ["polar-sidechain-backbone"]],
        ['R','L',251,'3.54x54','A','Y',391, ["hydrophobic", "van-der-waals"]],
        ['R','H',180,'2.50x50','A','Y',391, ["cation-pi", "hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
        ['R','R',170,'12.48x48','B','R',52, ["polar-sidechain-sidechain"]],
        ['R','L',255,'3.58x58','A','Q',384, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
        ['R','H',171,'12.49x49','B','D',312, ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
        ['R','K',334,'5.64x64','A','R',385, ["polar-backbone-sidechain", "van-der-waals"]],
        ['R','R',419,'8.60x60','B','G',310, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
        ['R','L',339,'-','A','L',394, ["hydrophobic"]],
        ['R','R',264,'4.40x41','A','Q',35, ["polar-backbone-sidechain"]],
        ['R','N',406,'8.47x47','A','E',392, ["hydrophobic", "polar-backbone-sidechain", "polar-sidechain-sidechain", "van-der-waals"]],
        ['R','A',256,'3.59x59','A','Q',384, ["polar-backbone-sidechain"]],
        ['R','S',352,'6.41x41','A','E',392, ["polar-sidechain-backbone"]],
        ['R','R',419,'8.60x60','B','A',309, ["hydrophobic", "van-der-waals"]],
        ['R','L',359,'6.48x48','A','Y',391, ["hydrophobic", "van-der-waals"]],
        ['R','E',408,'8.49x49','A','Q',390, ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
        ['R','K',342,'-','A','T',350, ["polar-sidechain-backbone", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
        ['R','L',254,'3.57x57','A','H',387, ["hydrophobic", "van-der-waals"]],
        ['R','N',407,'8.48x48','A','E',392, ["polar-backbone-sidechain"]],
        ['R','E',262,'4.38x39','A','Q',35, ["polar-backbone-sidechain", "van-der-waals"]],
        ['R','L',356,'6.45x45','A','L',393, ["hydrophobic"]],
        ['R','R',419,'8.60x60','B','F',292, ["polar-sidechain-backbone"]],
        ['R','V',405,'7.60x60','A','E',392, ["polar-backbone-sidechain", "van-der-waals"]],
        ['R','S',352,'6.41x41','A','L',394, ["polar-sidechain-backbone"]],
        ['R','N',338,'-','A','C',359, ["polar-sidechain-backbone"]],
        ['R','R',176,'2.46x46','A','Q',390, ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
        ['R','R',176,'2.46x46','A','Y',391, ["hydrophobic"]],
        ['R','K',415,'8.56x56','B','F',292, ["cation-pi", "hydrophobic", "van-der-waals"]],
        ['R','L',356,'6.45x45','A','Y',391, ["hydrophobic"]],
        ['R','L',339,'-','A','Y',358, ["hydrophobic", "van-der-waals"]],
        ['R','K',334,'5.64x64','A','D',381, ["hydrophobic", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
        ['R','R',419,'8.60x60','B','D',312, ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
        ['R','E',412,'8.53x53','B','D',312, ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
        ['R','N',406,'8.47x47','A','Q',390, ["polar-sidechain-backbone"]],
        ['R','T',353,'6.42x42','A','L',393, ["hydrophobic"]],
        ['R','A',256,'3.59x59','A','R',380, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
        ['R','S',261,'4.37x38','A','Q',35, ["hydrophobic", "polar-backbone-sidechain", "polar-sidechain-sidechain", "van-der-waals"]],
        ['R','Q',263,'4.39x40','A','Q',35, ["polar-backbone-sidechain", "van-der-waals"]],
        ['R','V',331,'5.61x61','A','L',388, ["hydrophobic", "van-der-waals"]],
        ['R','Q',263,'4.39x40','A','Q',31, ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
        ['R','L',339,'-','A','R',385, ["polar-backbone-sidechain"]],
        ['R','S',352,'6.41x41','A','L',393, ["polar-sidechain-backbone", "van-der-waals"]],
        ['R','L',401,'7.56x56','A','E',392, ["polar-backbone-sidechain", "van-der-waals"]],
        ],
        '6b3j' : [
        ['R','N',407,'8.48x48','A','E',392, ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
        ['R','K',334,'5.64x64','A','L',388, ["hydrophobic"]],
        ['R','K',415,'8.56x56','B','D',312, ["polar-sidechain-sidechain", "h-bond"]],
        ['R','R',176,'2.46x46','A','Y',391, ["cation-pi", "hydrophobic"]],
        ['R','L',255,'3.58x58','A','Q',384, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
        ['R','V',331,'5.61x61','A','L',393, ["hydrophobic"]],
        ['R','E',423,'8.64x64','B','Q',44, ["hydrophobic"]],
        ['R','S',261,'-','A','Q',35, ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
        ['R','L',251,'3.54x54','A','Y',391, ["hydrophobic", "van-der-waals"]],
        ['R','V',331,'5.61x61','A','L',394, ["hydrophobic"]],
        ['R','L',255,'3.58x58','A','L',388, ["hydrophobic"]],
        ['R','S',258,'-','A','I',383, ["hydrophobic"]],
        ['R','R',176,'2.46x46','A','Q',390, ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain"]],
        ['R','V',259,'-','A','V',217, ["hydrophobic"]],
        ['R','K',334,'5.64x64','A','R',385, ["hydrophobic", "polar-backbone-sidechain"]],
        ['R','R',348,'6.37x37','A','L',394, ["hydrophobic"]],
        ['R','Y',250,'3.53x53','A','Y',391, ["van-der-waals"]],
        ['R','S',352,'6.41x41','A','L',394, ["polar-sidechain-backbone"]],
        ['R','K',334,'5.64x64','A','Q',384, ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
        ['R','V',331,'5.61x61','A','L',388, ["hydrophobic"]],
        ['R','H',171,'12.49x49','B','D',312, ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
        ['R','R',419,'8.60x60','B','A',309, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
        ['R','L',255,'3.58x58','A','H',387, ["hydrophobic"]],
        ['R','E',262,'4.38x39','A','R',38, ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
        ['R','V',327,'5.57x57','A','L',393, ["hydrophobic"]],
        ['R','L',254,'3.57x57','A','H',387, ["hydrophobic", "van-der-waals"]],
        ['R','H',180,'2.50x50','A','Y',391, ["hydrophobic"]],
        ['R','K',334,'5.64x64','A','D',381, ["hydrophobic", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
        ['R','L',356,'6.45x45','A','L',393, ["hydrophobic"]],
        ['R','S',352,'6.41x41','A','L',393, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
        ['R','E',262,'4.38x39','A','K',34, ["polar-sidechain-sidechain"]],
        ],
        '6cmo' : [
        ['R','M',309,'7.56x56','A','G',352, ["van-der-waals"]],
        ['R','V',139,'3.54x54','A','L',348, ["hydrophobic"]],
        ['R','K',141,'3.56x56','A','D',193, ["polar-sidechain-sidechain"]],
        ['R','S',240,'-','A','E',318, ["hydrophobic", "polar-backbone-sidechain", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
        ['R','K',66,'12.48x48','B','D',312, ["polar-sidechain-sidechain"]],
        ['R','V',139,'3.54x54','A','N',347, ["hydrophobic"]],
        ['R','K',311,'8.48x48','A','K',349, ["polar-sidechain-backbone"]],
        ['R','T',242,'6.25x25','A','D',315, ["polar-sidechain-backbone"]],
        ['R','E',239,'-','A','E',318, ["polar-backbone-sidechain"]],
        ['R','T',242,'6.25x25','A','F',354, ["hydrophobic", "van-der-waals"]],
        ['R','Q',237,'5.72x72','A','D',341, ["hydrophobic"]],
        ['R','R',135,'3.50x50','A','L',353, ["hydrophobic"]],
        ['R','R',147,'34.55x55','A','R',32, ["hydrophobic", "polar-sidechain-backbone"]],
        ['R','A',246,'6.29x29','A','F',354, ["hydrophobic", "van-der-waals"]],
        ['R','E',249,'6.32x32','A','L',353, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
        ['R','A',246,'6.29x29','A','L',348, ["hydrophobic"]],
        ['R','K',311,'8.48x48','A','L',353, ["polar-sidechain-backbone"]],
        ['R','S',240,'-','A','K',345, ["polar-sidechain-sidechain", "h-bond"]],
        ['R','Q',237,'5.72x72','A','Y',320, ["polar-sidechain-sidechain"]],
        ['R','E',239,'-','A','Y',320, ["hydrophobic", "van-der-waals"]],
        ['R','E',249,'6.32x32','A','F',354, ["hydrophobic"]],
        ['R','K',311,'8.48x48','A','F',354, ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
        ['R','A',241,'6.24x24','A','E',318, ["polar-backbone-sidechain"]],
        ['R','K',245,'6.28x28','A','F',354, ["hydrophobic"]],
        ['R','N',310,'8.47x47','A','G',352, ["hydrophobic", "polar-sidechain-backbone"]],
        ['R','M',253,'6.36x36','A','L',353, ["hydrophobic", "van-der-waals"]],
        ['R','K',311,'8.48x48','A','G',352, ["polar-sidechain-backbone"]],
        ['R','R',147,'34.55x55','A','A',31, ["polar-sidechain-backbone"]],
        ['R','V',250,'6.33x33','A','L',353, ["hydrophobic"]],
        ['R','T',243,'6.26x26','A','D',341, ["polar-sidechain-sidechain"]],
        ],
        '6d9h' : [
        ['R','Q',210,'5.68x68','A','I',345, ["hydrophobic", "van-der-waals"]],
        ['R','K',224,'6.25x25','A','E',319, ["polar-sidechain-sidechain"]],
        ['R','R',291,'7.56x56','A','G',353, ["hydrophobic"]],
        ['R','I',207,'5.65x65','A','L',349, ["hydrophobic"]],
        ['R','Q',38,'12.48x48','B','F',335, ["hydrophobic", "van-der-waals"]],
        ['R','K',231,'6.32x32','A','F',355, ["polar-sidechain-backbone"]],
        ['R','L',113,'34.51x51','A','L',195, ["hydrophobic"]],
        ['R','R',108,'3.53x53','A','N',348, ["hydrophobic", "polar-backbone-sidechain", "polar-sidechain-sidechain"]],
        ['R','Q',210,'5.68x68','A','D',342, ["polar-sidechain-sidechain", "van-der-waals"]],
        ['R','P',112,'34.50x50','A','N',348, ["polar-backbone-sidechain"]],
        ['R','L',211,'5.69x69','A','K',346, ["hydrophobic"]],
        ['R','K',213,'5.71x71','A','D',342, ["polar-sidechain-sidechain"]],
        ['R','K',228,'6.29x29','A','F',355, ["cation-pi", "hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
        ['R','L',236,'6.37x37','A','L',354, ["hydrophobic", "van-der-waals"]],
        ['R','V',203,'5.61x61','A','L',354, ["hydrophobic"]],
        ['R','K',224,'6.25x25','A','D',316, ["hydrophobic", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
        ['R','I',207,'5.65x65','A','L',354, ["hydrophobic"]],
        ['R','V',109,'3.54x54','A','L',349, ["hydrophobic", "van-der-waals"]],
        ['R','R',108,'3.53x53','A','D',351, ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
        ['R','L',113,'34.51x51','A','T',341, ["hydrophobic"]],
        ['R','K',294,'8.49x49','A','D',351, ["hydrophobic", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
        ['R','P',112,'34.50x50','A','I',344, ["hydrophobic", "van-der-waals"]],
        ['R','P',112,'34.50x50','A','I',345, ["hydrophobic", "van-der-waals"]],
        ['R','L',113,'34.51x51','A','I',344, ["hydrophobic"]],
        ['R','R',105,'3.50x50','A','L',354, ["hydrophobic"]],
        ['R','I',232,'6.33x33','A','F',355, ["hydrophobic"]],
        ['R','I',292,'8.47x47','A','G',353, ["hydrophobic"]],
        ['R','R',108,'3.53x53','A','C',352, ["polar-sidechain-sidechain", "van-der-waals"]],
        ['R','K',228,'6.29x29','A','D',316, ["polar-sidechain-sidechain"]],
        ['R','N',37,'1.60x60','B','D',312, ["polar-sidechain-sidechain"]],
        ['R','F',45,'2.40x40','A','D',351, ["van-der-waals"]],
        ['R','L',113,'34.51x51','A','F',337, ["hydrophobic"]],
        ['R','D',42,'2.37x37','A','D',351, ["polar-sidechain-sidechain"]],
        ['R','R',105,'3.50x50','A','C',352, ["hydrophobic", "polar-sidechain-backbone"]],
        ['R','K',301,'8.56x56','B','D',312, ["polar-sidechain-backbone"]],
        ['R','I',292,'8.47x47','A','C',352, ["hydrophobic", "van-der-waals"]],
        ['R','I',232,'6.33x33','A','L',354, ["hydrophobic", "van-der-waals"]],
        ['R','K',301,'8.56x56','B','F',292, ["cation-pi", "hydrophobic", "van-der-waals"]],
        ],
        '6dde' : [
        ['R','R',263,'-','A','I',319, ["polar-sidechain-backbone"]],
        ['R','R',263,'-','A','Y',320, ["cation-pi", "hydrophobic", "van-der-waals"]],
        ['R','V',173,'34.51x51','A','F',336, ["hydrophobic", "van-der-waals"]],
        ['R','D',177,'34.55x55','A','R',32, ["polar-backbone-sidechain", "polar-sidechain-sidechain"]],
        ['R','L',259,'5.65x65','A','I',344, ["hydrophobic"]],
        ['R','V',169,'3.54x54','A','L',348, ["hydrophobic", "van-der-waals"]],
        ['R','L',176,'34.54x54','A','I',343, ["hydrophobic", "van-der-waals"]],
        ['R','L',176,'34.54x54','A','L',194, ["hydrophobic"]],
        ['R','V',173,'34.51x51','A','L',194, ["hydrophobic"]],
        ['R','M',264,'-','A','D',341, ["hydrophobic", "polar-backbone-sidechain"]],
        ['R','P',172,'34.50x50','A','I',343, ["hydrophobic", "van-der-waals"]],
        ['R','L',176,'34.54x54','A','R',32, ["hydrophobic", "polar-backbone-sidechain"]],
        ['R','R',277,'6.32x32','A','L',353, ["polar-sidechain-backbone"]],
        ['R','T',103,'2.39x39','A','C',351, ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
        ['R','V',173,'34.51x51','A','D',193, ["hydrophobic", "van-der-waals"]],
        ['R','L',259,'5.65x65','A','L',348, ["hydrophobic"]],
        ['R','D',340,'8.47x47','A','G',352, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
        ['R','M',281,'6.36x36','A','L',353, ["hydrophobic", "van-der-waals"]],
        ['R','I',278,'6.33x33','A','F',354, ["hydrophobic", "van-der-waals"]],
        ['R','E',270,'6.25x25','A','D',315, ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
        ['R','M',255,'5.61x61','A','L',353, ["hydrophobic"]],
        ['R','T',103,'2.39x39','A','D',350, ["polar-sidechain-backbone"]],
        ['R','V',262,'5.68x68','A','D',341, ["hydrophobic", "van-der-waals"]],
        ['R','K',271,'6.26x26','A','K',314, ["polar-sidechain-backbone", "van-der-waals"]],
        ['R','R',165,'3.50x50','A','L',353, ["hydrophobic"]],
        ['R','R',165,'3.50x50','A','C',351, ["polar-sidechain-sidechain", "van-der-waals"]],
        ['R','R',179,'34.57x57','A','N',347, ["polar-sidechain-sidechain", "van-der-waals"]],
        ['R','V',262,'5.68x68','A','I',344, ["hydrophobic"]],
        ['R','P',172,'34.50x50','A','T',340, ["hydrophobic"]],
        ['R','A',168,'3.53x53','A','N',347, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
        ['R','I',278,'6.33x33','A','L',348, ["hydrophobic"]],
        ['R','M',264,'-','A','T',316, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
        ['R','K',271,'6.26x26','A','D',315, ["hydrophobic", "van-der-waals"]],
        ['R','I',278,'6.33x33','A','L',353, ["hydrophobic"]],
        ['R','P',172,'34.50x50','A','I',344, ["hydrophobic", "van-der-waals"]],
        ['R','M',264,'-','A','K',345, ["hydrophobic"]],
        ['R','S',268,'6.23x23','A','D',315, ["polar-backbone-sidechain"]],
        ['R','R',258,'5.64x64','A','I',344, ["hydrophobic", "van-der-waals"]],
        ['R','R',263,'-','A','D',341, ["polar-backbone-sidechain"]],
        ['R','R',182,'4.40x40','A','R',24, ["hydrophobic", "polar-sidechain-sidechain"]],
        ],
        '6ddf' : [
        ['R','E',341,'8.48x48','A','F',354, ["hydrophobic"]],
        ['R','R',263,'-','A','I',319, ["polar-sidechain-backbone", "van-der-waals"]],
        ['R','V',173,'34.51x51','A','F',336, ["hydrophobic", "van-der-waals"]],
        ['R','L',176,'34.54x54','A','I',343, ["hydrophobic"]],
        ['R','L',176,'34.54x54','A','L',194, ["hydrophobic"]],
        ['R','M',264,'-','A','D',341, ["polar-backbone-sidechain"]],
        ['R','R',277,'6.32x32','A','F',354, ["polar-sidechain-backbone", "van-der-waals"]],
        ['R','P',172,'34.50x50','A','I',343, ["hydrophobic", "van-der-waals"]],
        ['R','R',277,'6.32x32','A','L',353, ["polar-sidechain-backbone", "van-der-waals"]],
        ['R','V',173,'34.51x51','A','L',194, ["hydrophobic"]],
        ['R','V',173,'34.51x51','A','D',193, ["hydrophobic"]],
        ['R','L',259,'5.65x65','A','I',344, ["hydrophobic"]],
        ['R','D',340,'8.47x47','A','G',352, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
        ['R','K',174,'34.52x52','A','D',193, ["polar-sidechain-sidechain"]],
        ['R','E',341,'8.48x48','A','L',353, ["polar-sidechain-backbone"]],
        ['R','R',165,'3.50x50','A','C',351, ["hydrophobic", "polar-sidechain-backbone"]],
        ['R','V',169,'3.54x54','A','L',348, ["hydrophobic", "van-der-waals"]],
        ['R','P',172,'34.50x50','A','N',347, ["polar-backbone-sidechain"]],
        ['R','P',172,'34.50x50','A','T',340, ["hydrophobic", "van-der-waals"]],
        ['R','A',168,'3.53x53','A','N',347, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
        ['R','I',278,'6.33x33','A','L',348, ["hydrophobic"]],
        ['R','K',271,'6.26x26','A','D',315, ["hydrophobic", "van-der-waals"]],
        ['R','M',255,'5.61x61','A','L',353, ["hydrophobic"]],
        ['R','K',271,'6.26x26','A','K',317, ["polar-sidechain-backbone"]],
        ['R','R',263,'-','A','D',341, ["hydrophobic", "polar-backbone-sidechain"]],
        ['R','I',278,'6.33x33','A','L',353, ["hydrophobic"]],
        ['R','V',262,'5.68x68','A','I',344, ["hydrophobic"]],
        ['R','R',263,'-','A','Y',320, ["cation-pi", "hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
        ['R','D',177,'34.55x55','A','R',32, ["polar-sidechain-sidechain"]],
        ['R','E',341,'8.48x48','A','G',352, ["polar-sidechain-backbone", "van-der-waals"]],
        ['R','P',172,'34.50x50','A','I',344, ["hydrophobic", "van-der-waals"]],
        ['R','A',168,'3.53x53','A','C',351, ["van-der-waals"]],
        ['R','D',164,'3.49x49','A','C',351, ["polar-sidechain-sidechain"]],
        ['R','L',176,'34.54x54','A','R',32, ["hydrophobic"]],
        ['R','T',103,'2.39x39','A','C',351, ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain"]],
        ['R','L',259,'5.65x65','A','L',348, ["hydrophobic", "van-der-waals"]],
        ['R','R',263,'-','A','E',318, ["van-der-waals"]],
        ['R','M',281,'6.36x36','A','L',353, ["hydrophobic", "van-der-waals"]],
        ['R','I',278,'6.33x33','A','F',354, ["hydrophobic", "van-der-waals"]],
        ['R','E',270,'6.25x25','A','D',315, ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
        ['R','T',103,'2.39x39','A','D',350, ["polar-sidechain-backbone", "van-der-waals"]],
        ['R','V',262,'5.68x68','A','D',341, ["hydrophobic", "van-der-waals"]],
        ['R','K',271,'6.26x26','A','K',314, ["polar-sidechain-backbone", "van-der-waals"]],
        ['R','R',165,'3.50x50','A','L',353, ["hydrophobic"]],
        ['R','R',179,'34.57x57','A','N',347, ["polar-sidechain-sidechain"]],
        ['R','D',340,'8.47x47','A','L',353, ["polar-sidechain-backbone"]],
        ['R','R',179,'34.57x57','A','C',351, ["polar-sidechain-sidechain", "van-der-waals"]],
        ['R','M',264,'-','A','T',316, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
        ['R','M',264,'-','A','K',345, ["hydrophobic"]],
        ['R','S',268,'6.23x23','A','D',315, ["polar-backbone-sidechain", "polar-sidechain-sidechain", "van-der-waals"]],
        ['R','R',258,'5.64x64','A','I',344, ["hydrophobic", "van-der-waals"]],
        ],
        '6g79' : [
        ['S','L',316,'6.37x37','A','L',353, ["hydrophobic"]],
        ['S','A',235,'5.65x65','A','L',348, ["hydrophobic"]],
        ['S','R',238,'5.68x68','A','D',341, ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
        ['S','R',308,'6.29x29','A','N',316, ["polar-sidechain-backbone"]],
        ['S','T',315,'6.36x36','A','L',353, ["hydrophobic", "polar-sidechain-backbone"]],
        ['S','V',155,'34.51x51','A','T',340, ["hydrophobic"]],
        ['S','R',147,'3.50x50','A','C',351, ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
        ['S','A',150,'3.53x53','A','N',347, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
        ['S','S',372,'7.56x56','A','C',351, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
        ['S','I',151,'3.54x54','A','N',347, ["hydrophobic"]],
        ['S','K',311,'6.32x32','A','G',352, ["polar-sidechain-backbone"]],
        ['S','A',154,'34.50x50','A','I',344, ["hydrophobic", "van-der-waals"]],
        ['S','I',151,'3.54x54','A','I',344, ["hydrophobic"]],
        ['S','A',312,'6.33x33','A','L',353, ["hydrophobic", "van-der-waals"]],
        ['S','R',238,'5.68x68','A','L',348, ["hydrophobic"]],
        ['S','T',315,'6.36x36','A','G',352, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
        ['S','I',239,'5.69x69','A','Y',354, ["hydrophobic"]],
        ['S','I',231,'5.61x61','A','L',353, ["hydrophobic", "van-der-waals"]],
        ['S','R',238,'5.68x68','A','I',344, ["hydrophobic", "van-der-waals"]],
        ['S','S',372,'7.56x56','A','G',352, ["hydrophobic"]],
        ['S','V',155,'34.51x51','A','F',336, ["hydrophobic"]],
        ['S','R',238,'5.68x68','A','Y',354, ["hydrophobic"]],
        ['S','A',154,'34.50x50','A','I',343, ["hydrophobic"]],
        ['S','R',308,'6.29x29','A','Y',354, ["cation-pi", "hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
        ['S','K',311,'6.32x32','A','Y',354, ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
        ['S','A',235,'5.65x65','A','L',353, ["hydrophobic"]],
        ['S','R',161,'34.57x57','A','N',347, ["polar-sidechain-sidechain"]],
        ['S','N',373,'8.47x47','A','G',350, ["polar-sidechain-backbone"]],
        ['S','I',151,'3.54x54','A','L',348, ["hydrophobic", "van-der-waals"]],
        ],
        '6gdg' : [
        ['A','Q',207,'5.68x68','D','L',384, ["hydrophobic"]],
        ['A','L',110,'34.51x51','D','I',373, ["hydrophobic"]],
        ['A','Y',112,'34.53x53','D','H',377, ["edge-to-face", "face-to-edge", "pi-cation", "hydrophobic"]],
        ['A','Q',207,'5.68x68','D','Q',374, ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
        ['A','R',102,'3.50x50','D','Y',381, ["cation-pi", "hydrophobic", "van-der-waals"]],
        ['A','N',36,'12.49x49','B','D',312, ["hydrophobic", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
        ['A','R',111,'34.52x52','D','K',216, ["hydrophobic", "polar-sidechain-backbone"]],
        ['A','M',211,'5.72x72','D','Y',348, ["hydrophobic"]],
        ['A','L',110,'34.51x51','D','H',41, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
        ['A','P',109,'34.50x50','D','R',370, ["hydrophobic"]],
        ['A','I',106,'3.54x54','D','H',377, ["hydrophobic", "van-der-waals"]],
        ['A','S',35,'12.48x48','B','D',333, ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
        ['A','L',110,'34.51x51','D','V',217, ["hydrophobic"]],
        ['A','S',234,'6.36x36','D','L',383, ["hydrophobic", "van-der-waals"]],
        ['A','E',294,'8.49x49','D','Q',380, ["polar-backbone-sidechain"]],
        ['A','L',110,'34.51x51','D','C',369, ["hydrophobic"]],
        ['A','I',106,'3.54x54','D','Q',374, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
        ['A','S',35,'12.48x48','B','D',312, ["hydrophobic"]],
        ['A','S',35,'12.48x48','B','F',335, ["hydrophobic", "van-der-waals"]],
        ['A','I',106,'3.54x54','D','Y',381, ["hydrophobic"]],
        ['A','I',200,'5.61x61','D','L',383, ["hydrophobic", "van-der-waals"]],
        ['A','N',34,'1.60x60','B','D',312, ["polar-sidechain-sidechain"]],
        ['A','L',208,'5.69x69','D','L',384, ["hydrophobic"]],
        ['A','Q',210,'5.71x71','D','D',371, ["polar-sidechain-sidechain"]],
        ['A','S',234,'6.36x36','D','E',382, ["polar-sidechain-backbone"]],
        ['A','L',110,'34.51x51','D','R',370, ["hydrophobic"]],
        ['A','R',111,'34.52x52','D','D',215, ["polar-sidechain-backbone"]],
        ['A','R',291,'7.56x56','D','Y',381, ["polar-sidechain-backbone"]],
        ['A','R',293,'8.48x48','D','Q',380, ["polar-backbone-sidechain"]],
        ['A','L',110,'34.51x51','D','F',366, ["hydrophobic", "van-der-waals"]],
        ['A','Q',207,'5.68x68','D','D',371, ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
        ['A','R',111,'34.52x52','D','V',217, ["hydrophobic", "van-der-waals"]],
        ['A','Q',207,'5.68x68','D','L',378, ["hydrophobic"]],
        ['A','L',235,'6.37x37','D','L',383, ["hydrophobic", "van-der-waals"]],
        ['A','A',204,'5.65x65','D','L',378, ["hydrophobic"]],
        ['A','I',106,'3.54x54','D','L',378, ["hydrophobic", "van-der-waals"]],
        ['A','A',231,'6.33x33','D','L',383, ["hydrophobic", "van-der-waals"]],
        ['A','H',230,'6.32x32','D','E',382, ["polar-sidechain-backbone"]],
        ['A','N',113,'34.54x54','D','A',39, ["hydrophobic"]],
        ['A','P',109,'34.50x50','D','Q',374, ["hydrophobic", "van-der-waals"]],
        ['A','R',293,'8.48x48','D','E',382, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
        ['A','A',105,'3.53x53','D','H',377, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
        ['A','G',114,'34.55x55','D','A',39, ["hydrophobic"]],
        ['A','A',203,'5.64x64','D','L',378, ["hydrophobic"]],
        ['A','Q',207,'5.68x68','D','R',375, ["polar-sidechain-backbone"]],
        ['A','R',107,'3.55x55','D','R',370, ["polar-backbone-sidechain", "polar-sidechain-sidechain", "van-der-waals"]],
        ['A','R',291,'7.56x56','D','E',382, ["polar-backbone-sidechain", "polar-sidechain-backbone", "van-der-waals"]],
        ['A','N',113,'34.54x54','D','R',38, ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
        ['A','P',109,'34.50x50','D','I',373, ["hydrophobic", "van-der-waals"]],
        ['A','Q',38,'12.51x51','B','R',52, ["polar-sidechain-sidechain", "van-der-waals"]],
        ],
        '6e3y' : [
        ['R','V',243,'3.60x60','A','R',380, ["polar-backbone-sidechain"]],
        ['R','L',240,'3.57x57','A','Y',391, ["hydrophobic"]],
        ['R','V',242,'3.59x59','A','R',380, ["polar-backbone-sidechain", "van-der-waals"]],
        ['R','I',241,'3.58x58','A','H',387, ["hydrophobic"]],
        ['R','G',389,'8.48x48','A','E',392, ["polar-backbone-sidechain", "van-der-waals"]],
        ['R','E',248,'-','A','A',39, ["hydrophobic"]],
        ['R','K',319,'5.64x64','A','Q',384, ["hydrophobic", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
        ['R','T',323,'-','A','Y',358, ["hydrophobic"]],
        ['R','L',240,'3.57x57','A','H',387, ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
        ['R','F',246,'-','A','C',379, ["hydrophobic"]],
        ['R','L',316,'5.61x61','A','L',388, ["hydrophobic", "van-der-waals"]],
        ['R','K',333,'6.37x37','A','L',394, ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
        ['R','I',241,'3.58x58','A','L',388, ["hydrophobic"]],
        ['R','R',336,'6.40x40','A','L',393, ["polar-sidechain-backbone"]],
        ['R','I',241,'3.58x58','A','Y',391, ["hydrophobic"]],
        ['R','R',336,'6.40x40','A','L',394, ["polar-sidechain-sidechain"]],
        ['R','K',319,'5.64x64','A','D',381, ["hydrophobic", "van-der-waals"]],
        ['R','F',246,'-','A','R',380, ["pi-cation", "hydrophobic", "van-der-waals"]],
        ['R','N',388,'8.47x47','A','E',392, ["polar-sidechain-sidechain", "van-der-waals"]],
        ['R','K',319,'5.64x64','A','L',388, ["hydrophobic"]],
        ['R','S',168,'12.49x49','B','D',312, ["polar-sidechain-sidechain"]],
        ['R','R',173,'2.46x46','A','Y',391, ["cation-pi", "hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
        ['R','R',173,'2.46x46','A','Q',390, ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
        ['R','E',248,'-','A','R',38, ["hydrophobic", "van-der-waals"]],
        ['R','F',246,'-','A','I',383, ["hydrophobic", "van-der-waals"]],
        ['R','R',397,'8.56x56','B','H',311, ["polar-sidechain-backbone"]],
        ['R','R',397,'8.56x56','B','D',312, ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
        ['R','H',177,'2.50x50','A','Y',391, ["hydrophobic"]],
        ['R','S',168,'12.49x49','B','F',335, ["hydrophobic"]],
        ['R','L',320,'5.65x65','A','L',394, ["hydrophobic", "van-der-waals"]],
        ['R','V',242,'3.59x59','A','Q',384, ["polar-backbone-sidechain"]],
        ['R','V',245,'-','A','I',383, ["hydrophobic", "van-der-waals"]],
        ['R','F',387,'7.60x60','A','E',392, ["polar-backbone-sidechain"]],
        ['R','I',241,'3.58x58','A','L',393, ["hydrophobic"]],
        ['R','R',336,'6.40x40','A','E',392, ["polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
        ['R','L',237,'3.54x54','A','Y',391, ["hydrophobic", "van-der-waals"]],
        ['R','F',246,'-','A','F',219, ["hydrophobic"]],
        ['R','I',241,'3.58x58','A','Q',384, ["polar-backbone-sidechain", "van-der-waals"]],
        ['R','K',319,'5.64x64','A','R',385, ["polar-backbone-sidechain", "van-der-waals"]],
        ['R','F',246,'-','A','H',41, ["pi-cation", "hydrophobic"]],
        ['R','L',341,'6.45x45','A','L',393, ["hydrophobic", "van-der-waals"]],
        ['R','I',312,'5.57x57','A','L',393, ["hydrophobic"]],
        ['R','V',245,'-','A','H',387, ["polar-backbone-sidechain"]],
        ['R','L',316,'5.61x61','A','L',393, ["hydrophobic"]],
        ['R','F',246,'-','A','V',217, ["hydrophobic", "van-der-waals"]],
        ['R','K',167,'12.48x48','B','R',52, ["hydrophobic", "polar-backbone-sidechain"]],
        ['R','V',245,'-','A','R',380, ["hydrophobic", "van-der-waals"]],
        ['R','F',246,'-','A','F',376, ["hydrophobic", "van-der-waals"]],
        ['R','R',397,'8.56x56','B','F',292, ["cation-pi"]],
        ['R','V',245,'-','A','Q',384, ["hydrophobic", "van-der-waals"]],
        ['R','L',316,'5.61x61','A','L',394, ["hydrophobic", "van-der-waals"]],
        ]
}

    complex_info = [
        {
            'pdb_id': '3sn6',
            'receptor': 'beta2',
            'gprotein': 'GNAS2_BOVIN',
            'alternative_gprotein': 'GNAS2_HUMAN'
            },
        {
            'pdb_id': '4x1h',
            'receptor': '',
            'gprotein': '',
            'alternative_protein': ''
            },
        {
            'pdb_id': '5g53',
            'receptor': '',
            'gprotein': '',
            'alternative_protein': ''
            },
        {
            'pdb_id': '5g53',
            'receptor': '',
            'gprotein': '',
            'alternative_protein': ''
            },
        {
            'pdb_id': '5uz7',
            'receptor': '',
            'gprotein': '',
            'alternative_protein': ''
            },
        {
            'pdb_id': '5vai',
            'receptor': '',
            'gprotein': '',
            'alternative_protein': ''
            },
        {
            'pdb_id': '6b3j',
            'receptor': '',
            'gprotein': '',
            'alternative_protein': ''
            },
        {
            'pdb_id': '6cmo',
            'receptor': '',
            'gprotein': '',
            'alternative_protein': ''
            },
        {
            'pdb_id': '6d9h',
            'receptor': '',
            'gprotein': '',
            'alternative_protein': ''
            },
        {
            'pdb_id': '6dde',
            'receptor': '',
            'gprotein': '',
            'alternative_protein': ''
            },
        {
            'pdb_id': '6ddf',
            'receptor': '',
            'gprotein': '',
            'alternative_protein': ''
            },
        {
            'pdb_id': '6g79',
            'receptor': '',
            'gprotein': '',
            'alternative_protein': ''
            },
        {
            'pdb_id': '6gdg',
            'receptor': '',
            'gprotein': '',
            'alternative_protein': ''
            },
    ]

    ps = ProteinConformation.objects.filter(
        protein__sequence_type__slug='wt',
        protein__species__common_name="Human",
        protein__family__slug__startswith='00',  # receptors, no gproteins
        # structure__refined=False
        ).values(
            name = F('protein__name'),
            entry_name = F('protein__entry_name'),
            pdb_id = F('structure__pdb_code__index'),
            rec_id = F('protein__id'),
            protein_family = F('protein__family__parent__name'),
            protein_class = F('protein__family__parent__parent__parent__name'),
            ligand = F('protein__endogenous_ligands__properities__ligand_type__name')
        )

    # data = Structure.objects.filter(
    #         # protein_conformation__protein__sequence_type__slug='wt',
    #         protein_conformation__protein__species__common_name="Human",
    #         protein_conformation__protein__family__slug__startswith='00',
    #         refined=False
    #     ).select_related(
    #         "state",
    #         "pdb_code__web_resource",
    #         "protein_conformation__protein__species",
    #         "protein_conformation__protein__source",
    #         "protein_conformation__protein__family__parent__parent__parent",
    #         "publication__web_link__web_resource"
    #     ).prefetch_related(
    #         "stabilizing_agents",
    #         "construct__crystallization__crystal_method",
    #         "protein_conformation__protein__parent__endogenous_ligands__properities__ligand_type",
    #         "protein_conformation__site_protein_conformation__site"
    #     )

    # ps = []
    # for s in data:
    #     r = {}
    #     r['seq_slug'] = s.protein_conformation.protein.sequence_type.slug
    #     r['pdb_id'] = s.pdb_code.index
    #     r['rec_id'] = s.protein_conformation.protein_id
    #     r['name'] = s.protein_conformation.protein.parent.name
    #     r['entry_name'] = s.protein_conformation.protein.parent.entry_name
    #     r['protein'] = s.protein_conformation.protein.parent.entry_short()
    #     r['protein_long'] = s.protein_conformation.protein.parent.short()
    #     r['protein_family'] = s.protein_conformation.protein.parent.family.parent.short()
    #     r['class'] = s.protein_conformation.protein.parent.family.parent.parent.parent.short()
    #     r['species'] = s.protein_conformation.protein.species.common_name
    #     r['date'] = str(s.publication_date)
    #     r['state'] = s.state.name
    #     r['representative'] = 'Yes' if s.representative else 'No'
    #     ps.append(r)

    names = set(pi['entry_name'] for pi in ps)
    residuelist = Residue.objects.filter(
            protein_conformation__protein__entry_name__in=names,
            # protein_conformation__structure__refined=False
            ).prefetch_related(
                'protein_conformation'
            ).values(
                pdb_id = F('protein_conformation__structure__pdb_code__index'),
                rec_id = F('protein_conformation__protein__id'),
                name = F('protein_conformation__protein__name'),
                rec_aa = F('amino_acid'),
                rec_gn = F('display_generic_number__label')
            )

    interactions_metadata = complex_info
    context = {
        'interactions': dataset,
        'interactions_metadata': interactions_metadata,
        'ps': json.dumps(list(ps)),
        'rs': json.dumps(list(residuelist)),
        }

    return render(request, 'signprot/matrix.html', context)

def IMSequenceSignature(request):

    import time
    t1 = time.time()

    import re
    from itertools import chain

    from protein.models import Protein
    from protein.models import ProteinSegment
    from residue.models import ResidueGenericNumberEquivalent
    from seqsign.sequence_signature import SignatureMatch
    from seqsign.sequence_signature import SequenceSignature

    # example data
    segments = list(ProteinSegment.objects.filter(proteinfamily='GPCR'))
    pos_set = ["5ht2c_human", "acm4_human", "drd1_human"]
    neg_set = ["agtr1_human", "ednrb_human", "gnrhr_human"]

    # receive data

    # get pos/neg set objects
    pos_set = Protein.objects.filter(entry_name__in=pos_set).select_related('residue_numbering_scheme', 'species')
    neg_set = Protein.objects.filter(entry_name__in=neg_set).select_related('residue_numbering_scheme', 'species')

    # res numbers
    # segments = []
    # # label looks like "2x51"
    # gen_object = ResidueGenericNumberEquivalent.objects.get(label=s, scheme__slug='gpcrdb')
    # segments.append(gen_object)
    # # a.load_segments(gen_list)

    # Calculate Sequence Signature
    signature = SequenceSignature()
    signature.setup_alignments(segments, pos_set, neg_set)
    signature.calculate_signature()

    # process data for return
    signature_data = signature.prepare_display_data()

    # FEATURES
    feats = [feature for feature in signature_data['a_pos'].features]
    len_feats = len(feats)

    trans = {
        'N-term': 'N',
        'TM1': 1,
        'ICL1': 12,
        'TM2': 2,
        'ECL1': 23,
        'TM3': 3,
        'ICL2': 34,
        'TM4': 4,
        'ECL2': 45,
        'TM5': 5,
        'ICL3': 56,
        'TM6': 6,
        'ECL3': 67,
        'TM7': 7,
        'ICL4': 78,
        'H8': 8,
        'C-term': 'C',
    }

    # GET GENERIC NUMBERS
    generic_numbers = []
    for _, segments in signature_data['common_generic_numbers'].items():
        for elem, num in segments.items():
            gnl = []
            for x, dn in num.items():
                if dn != '':
                    rexp = r'(?<=<b>)\d{1,}|\.?\d{2,}[\-?\d{2,}]*|x\d{2,}'
                    gn = re.findall(rexp, dn)
                else:
                    gn = ''.join([str(trans[elem]), '.', str(x)])
                gnl.append(''.join(gn))
            generic_numbers.append(gnl)


    # FEATURE FREQUENCIES
    signature_features = {}
    x = 0
    for i, feature in enumerate(signature_data['feats_signature']):
        for j, segment in enumerate(feature):
            for k, freq in enumerate(segment):
                # freq0: score
                # freq1: level of conservation
                # freq2: a - b explanation
                try:
                    signature_features[x] = {
                        'feature': str(feats[i]),
                        'gn': str(generic_numbers[j][k]),
                        'freq': int(freq[0]),
                        'cons': int(freq[1]),
                        'expl': str(freq[2]),
                    }
                    x += 1
                except IndexError as e:
                    print(e)


    # SIGNATURE CONSENSUS
    generic_numbers_flat = list(chain.from_iterable(generic_numbers))
    sigcons = {}
    x = 0
    for segment, cons in signature_data['signature_consensus'].items():
        for i, pos in enumerate(cons):
            # pos0: Code
            # pos1: Name
            # pos2: Score
            # pos3: level of conservation
            sigcons[x] = {
                'gn': str(generic_numbers_flat[x]),
                'code': str(pos[0]),
                'name': str(pos[1]),
                'score': int(pos[2]),
                'cons': int(pos[3]),
            }
            x += 1

    # define list of features to keep
    # subset results
    # pass back to front
    res = {
        'cons': sigcons,
        'feat': signature_features,
    }

    t2 = time.time()
    print('Runtime: {}'.format((t2-t1)*1000.0))

    return JsonResponse(res, safe=False)

