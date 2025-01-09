from django.shortcuts import render, redirect
from django.conf import settings
from django.http import HttpResponse
from django.views.generic import TemplateView
from django.db.models import Case, When
from django.core.cache import cache
from django.core.cache import caches

try:
    cache_alignment = caches['alignments']
except:
    cache_alignment = cache

from alignment.functions import get_proteins_from_selection
from common import definitions
from common.selection import Selection, SelectionItem
from common.views import AbsTargetSelection, AbsTargetSelectionTable
from common.views import AbsSegmentSelection
from common.views import AbsMiscSelection
from structure.functions import BlastSearch
from protwis.context_processors import site_title

# from common.alignment_SITE_NAME import Alignment
Alignment = getattr(__import__('common.alignment_' + settings.SITE_NAME, fromlist=['Alignment']), 'Alignment')
from protein.models import Protein, ProteinSegment, ProteinFamily, ProteinSet
from residue.models import ResidueNumberingScheme, ResiduePositionSet
from seqsign.sequence_signature import SequenceSignature, signature_score_excel

from collections import OrderedDict
from copy import deepcopy
import hashlib
import inspect
from io import BytesIO
import itertools
import json
import numpy as np
import os
import xlsxwriter
import xlrd


# class TargetSelection(AbsTargetSelection):
#     step = 1
#     number_of_steps = 2
#     filter_tableselect = False
#     docs = 'sequences.html#structure-based-alignments'
#     selection_boxes = OrderedDict([
#         ('reference', False),
#         ('targets', True),
#         ('segments', False),
#     ])
#     buttons = {
#         'continue': {
#             'label': 'Continue to next step',
#             'url': '/alignment/segmentselection',
#             'color': 'success',
#         },
#     }

class TargetSelection(AbsTargetSelectionTable):
    step = 1
    number_of_steps = 2
    filter_tableselect = False
    docs = 'sequences.html#structure-based-alignments'
    title = "SELECT RECEPTORS"
    description = 'Select receptors in the table (below) or browse the classification tree (right). You can select entire' \
        + ' families or individual receptors.\n\nOnce you have selected all your receptors, click the green button.'
    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])
    buttons = {
        'continue': {
            'label': 'Next',
            'onclick': "submitSelection('/alignment/segmentselection');",
            'color': 'success',
        },
    }

class TargetSelectionGprotein(AbsTargetSelection):
    step = 1
    number_of_steps = 2
    psets = False
    filters = True
    filter_gprotein = True

    docs = 'sequences.html#structure-based-alignments'

    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])
    buttons = {
        'continue': {
            'label': 'Continue to next step',
            'url': '/alignment/segmentselectiongprot',
            'color': 'success',
        },
    }
    try:
        if ProteinFamily.objects.filter(slug="100_001").exists():
            ppf = ProteinFamily.objects.get(slug="100_001")
            pfs = ProteinFamily.objects.filter(parent=ppf.id)
            ps = Protein.objects.filter(family=ppf)

            tree_indent_level = []
            action = 'expand'
            # remove the parent family (for all other families than the root of the tree, the parent should be shown)
            del ppf
    except Exception as e:
        pass

class TargetSelectionArrestin(AbsTargetSelection):
    step = 1
    number_of_steps = 2
    psets = False
    filters = True
    filter_gprotein = True

    docs = 'sequences.html#structure-based-alignments'

    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])
    buttons = {
        'continue': {
            'label': 'Continue to next step',
            'url': '/alignment/segmentselectionarrestin',
            'color': 'success',
        },
    }
    try:
        if ProteinFamily.objects.filter(slug="200_000").exists():
            ppf = ProteinFamily.objects.get(slug="200_000")
            pfs = ProteinFamily.objects.filter(parent=ppf.id)
            ps = Protein.objects.filter(family=ppf)

            tree_indent_level = []
            action = 'expand'
            # remove the parent family (for all other families than the root of the tree, the parent should be shown)
            del ppf
    except Exception as e:
        pass

class SegmentSelection(AbsSegmentSelection):
    step = 2
    number_of_steps = 2
    docs = 'sequences.html#structure-based-alignments'
    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', True),
    ])
    buttons = {
        'continue': {
            'label': 'Show alignment',
            'url': '/alignment/render',
            'color': 'success',
        },
    }

class SegmentSelectionGprotein(AbsSegmentSelection):
    step = 2
    number_of_steps = 2
    docs = 'sequences.html#structure-based-alignments'
    description = 'Select sequence segments in the middle column for G proteins. You can expand every structural element and select individual' \
        + ' residues by clicking on the down arrows next to each helix, sheet or loop.\n\n You can select the full sequence or show all structured regions at the same time.\n\nSelected segments will appear in the' \
        + ' right column, where you can edit the list.\n\nOnce you have selected all your segments, click the green' \
        + ' button.'

    template_name = 'common/segmentselection.html'

    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', True),
    ])
    buttons = {
        'continue': {
            'label': 'Show alignment',
            'url': '/alignment/render',
            'color': 'success',
        },
    }

    position_type = 'gprotein'
    rsets = ResiduePositionSet.objects.filter(name__in=['Gprotein Barcode', 'YM binding site']).prefetch_related('residue_position')

    ss = ProteinSegment.objects.filter(partial=False, proteinfamily='Alpha').prefetch_related('generic_numbers')
    ss_cats = ss.values_list('category').order_by('category').distinct('category')


class SegmentSelectionArrestin(AbsSegmentSelection):
    step = 2
    number_of_steps = 2
    docs = 'sequences.html#structure-based-alignments'
    description = 'Select sequence segments in the middle column for beta and visual arrestins. You can expand every structural element and select individual' \
        + ' residues by clicking on the down arrows next to each helix, sheet or loop.\n\n You can select the full sequence or show all structured regions at the same time.\n\nSelected segments will appear in the' \
        + ' right column, where you can edit the list.\n\nOnce you have selected all your segments, click the green' \
        + ' button.'

    template_name = 'common/segmentselection.html'

    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', True),
    ])
    buttons = {
        'continue': {
            'label': 'Show alignment',
            'url': '/alignment/render',
            'color': 'success',
        },
    }

    position_type = 'arrestin'

    ## Add some Arrestin specific positions
    rsets = ResiduePositionSet.objects.filter(name__in=['Arrestin interface']).prefetch_related('residue_position')

    ## ProteinSegment for different proteins
    ss = ProteinSegment.objects.filter(partial=False, proteinfamily='Arrestin').prefetch_related('generic_numbers')
    ss_cats = ss.values_list('category').order_by('category').distinct('category')


class BlastSearchInput(AbsMiscSelection):
    step = 1
    number_of_steps = 1
    docs = 'sequences.html#similarity-search-blast'
    title = 'BLAST search'
    description = 'Enter a sequence into the text box and press the green button.'
    buttons = {
        'continue': {
            'label': 'BLAST',
            'onclick': 'document.getElementById(\'form\').submit()',
            'color': 'success',
        },
    }
    selection_boxes = {}
    blast_input = True


class BlastSearchResults(TemplateView):
    """
    An interface for blast similarity search of the input sequence.
    """
    template_name="blast/blast_search_results.html"

    def post(self, request, *args, **kwargs):

        if 'human' in request.POST.keys():
            blast = BlastSearch(blastdb=os.sep.join([settings.STATICFILES_DIRS[0], 'blast', 'protwis_human_blastdb']), top_results=50)
            blast_out = blast.run(request.POST['input_seq'])
        else:
            blast = BlastSearch(top_results=50)
            blast_out = blast.run(request.POST['input_seq'])

        context = {}
        context['results'] = [(Protein.objects.get(pk=x[0]), x[1]) for x in blast_out]
        context["input"] = request.POST['input_seq']

        return render(request, self.template_name, context)


def render_alignment(request):
    # get the user selection from session
    simple_selection = request.session.get('selection', False)
    if simple_selection == False or not simple_selection.targets:
        return redirect("/alignment/targetselection")

    # create an alignment object
    a = Alignment()

    # load data from selection into the alignment
    # only show wildtype protein entries if selection type is family
    if len([t for t in simple_selection.targets if t.type=='family'])>0:
        a.load_proteins_from_selection(simple_selection, only_wildtype=True)
    else:
        a.load_proteins_from_selection(simple_selection)
    a.load_segments_from_selection(simple_selection)

    key = "ALIGNMENT_" + a.get_hash()
    return_html = cache_alignment.get(key)

    if return_html==None:
        # build the alignment data matrix
        check = a.build_alignment()
        if check == 'Too large':
            return render(request, 'alignment/error.html', {'proteins': len(a.proteins), 'residues':a.number_of_residues_total})

        # calculate consensus sequence + amino acid and feature frequency
        a.calculate_statistics()
        num_of_sequences = len(a.proteins)
        num_residue_columns = len(a.positions) + len(a.segments)
        segment_headers_to_hide = [seg for seg, posis in a.segments.items() if len(posis)==0]

        return_html = render(request, 'alignment/alignment.html', {'a': a, 'num_of_sequences': num_of_sequences,
            'num_residue_columns': num_residue_columns, 'segment_headers_to_hide': segment_headers_to_hide})

    cache_alignment.set(key, return_html, 60*60*24*7) #set alignment cache one week

    return return_html

def render_family_alignment(request, slug):
    # create an alignment object
    a = Alignment()

    # fetch proteins and segments
    proteins = Protein.objects.filter(family__slug__startswith=slug, sequence_type__slug='wt')

    if len(proteins)>50 and len(slug.split("_"))<4:
        # If alignment is going to be too big, only pick human.
        proteins = Protein.objects.filter(family__slug__startswith=slug, sequence_type__slug='wt', species__latin_name='Homo sapiens')

    if slug.startswith('100'):

        gsegments = definitions.G_PROTEIN_SEGMENTS

        preserved = Case(*[When(slug=pk, then=pos) for pos, pk in enumerate(gsegments['Full'])])
        segments = ProteinSegment.objects.filter(slug__in=gsegments['Full'], partial=False).order_by(preserved)

    elif slug.startswith('200'):
        arrsegments = definitions.ARRESTIN_SEGMENTS

        preserved = Case(*[When(slug=pk, then=pos) for pos, pk in enumerate(arrsegments['Full'])])
        segments = ProteinSegment.objects.filter(slug__in=arrsegments['Full'], partial=False).order_by(preserved)

    else:
        segments = ProteinSegment.objects.filter(partial=False, proteinfamily='GPCR')
        if len(proteins)>50:
            # if a lot of proteins, exclude some segments
            segments = ProteinSegment.objects.filter(partial=False, proteinfamily='GPCR').exclude(slug__in=['N-term','C-term'])
        if len(proteins)>200:
            # if many more proteins exluclude more segments
            segments = ProteinSegment.objects.filter(partial=False, proteinfamily='GPCR').exclude(slug__in=['N-term','C-term']).exclude(category='loop')

    protein_ids = []
    for p in proteins:
        protein_ids.append(p.pk)
    protein_list = ','.join(str(x) for x in sorted(protein_ids))

    #create unique proteins_id
    segments_ids = []
    for s in segments:
        segments_ids.append(s.slug)
    segments_list = ','.join(str(x) for x in sorted(segments_ids))

    # Store proteins and segments as selection to enable Fasta/Excel/CSV downloads
    selection = Selection()
    for prot in proteins:
        selection.add('targets', 'protein', SelectionItem('protein', prot))
    for segment in segments:
        selection.add('segments', 'protein_segment', SelectionItem('protein_segment', segment))
    request.session['selection'] = selection.exporter()

    s = str(protein_list+"_"+segments_list)
    key = "ALIGNMENT_"+hashlib.md5(s.encode('utf-8')).hexdigest()
    return_html = cache_alignment.get(key)
    if return_html==None:
        # load data into the alignment
        a.load_proteins(proteins)
        a.load_segments(segments)

        # build the alignment data matrix
        a.build_alignment()

        # calculate consensus sequence + amino acid and feature frequency
        a.calculate_statistics()

        num_of_sequences = len(a.proteins)
        num_residue_columns = len(a.positions) + len(a.segments)

        return_html = render(request, 'alignment/alignment.html', {'a': a, 'num_of_sequences': num_of_sequences,
        'num_residue_columns': num_residue_columns})

    #update it if used
    cache_alignment.set(key,return_html, 60*60*24*7) #set alignment cache one week

    return return_html

def render_fasta_alignment(request):
    # get the user selection from session
    simple_selection = request.session.get('selection', False)

    # create an alignment object
    a = Alignment()
    a.show_padding = False

    # load data from selection into the alignment
    a.load_proteins_from_selection(simple_selection)
    a.load_segments_from_selection(simple_selection)

    # build the alignment data matrix
    a.build_alignment()

    response = render(request, 'alignment/alignment_fasta.html', context={'a': a}, content_type='text/fasta')
    response['Content-Disposition'] = "attachment; filename=" + site_title(request)["site_title"] + "_alignment.fasta"
    return response

def render_fasta_family_alignment(request, slug):
    # create an alignment object
    a = Alignment()
    a.show_padding = False

    # fetch proteins and segments
    proteins = Protein.objects.filter(family__slug__startswith=slug, sequence_type__slug='wt')
    segments = ProteinSegment.objects.filter(partial=False)

    # load data into the alignment
    a.load_proteins(proteins)
    a.load_segments(segments)

    # build the alignment data matrix
    a.build_alignment()

    response = render(request, 'alignment/alignment_fasta.html', context={'a': a}, content_type='text/fasta')
    response['Content-Disposition'] = "attachment; filename=" + site_title(request)["site_title"] + "_alignment.fasta"
    return response

def render_csv_alignment(request):
    # get the user selection from session
    simple_selection = request.session.get('selection', False)

    # create an alignment object
    a = Alignment()
    a.show_padding = False

    # load data from selection into the alignment
    a.load_proteins_from_selection(simple_selection)
    a.load_segments_from_selection(simple_selection)

    # build the alignment data matrix
    a.build_alignment()

    # calculate consensus sequence + amino acid and feature frequency
    a.calculate_statistics()

    response = render(request, 'alignment/alignment_csv.html', context={'a': a}, content_type='text/csv')
    response['Content-Disposition'] = "attachment; filename=" + site_title(request)["site_title"] + "_alignment.csv"
    return response

# Excel download based on seq. signature tool
def render_alignment_excel(request):

    # Grab all targets
    targets = request.session.get('selection', False)

    # create placeholder seq signature
    signature = SequenceSignature()
    signature.setup_alignments_from_selection(targets, targets)

    # calculate the signture
    signature.calculate_signature()
    signature.calculate_zscales_signature()

    outstream = BytesIO()
    wb = xlsxwriter.Workbook(outstream, {'in_memory': True})

    # Sequence alignment of targets
    signature.prepare_excel_worksheet(
        wb,
        'Alignment',
        'positive',
        'alignment'
    )

    # Residue properties stats
    signature.prepare_excel_worksheet(
        wb,
        'Property_conservation',
        'positive',
        'features'
    )
    # Z-scales
    signature.zscales_excel(
        wb,
        "Z-scales",
        'positive'
    )
    wb.close()
    outstream.seek(0)
    response = HttpResponse(
        outstream.read(),
        content_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )
    response['Content-Disposition'] = "attachment; filename=" + site_title(request)["site_title"] + "_alignment.xlsx"

    return response
