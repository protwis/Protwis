from django.shortcuts import render, redirect
from django.conf import settings
from django.http import HttpResponse


from alignment.functions import get_proteins_from_selection
from common.selection import Selection
#from common.views import AbsTargetSelection
from common.views import AbsTargetSelectionTable
from common.views import AbsSegmentSelection
from seqsign.sequence_signature import SequenceSignature, SignatureMatch, signature_score_excel

Alignment = getattr(__import__('common.alignment_' + settings.SITE_NAME, fromlist=['Alignment']), 'Alignment')

from collections import OrderedDict
from copy import deepcopy
from io import BytesIO
import xlsxwriter


class PosTargetSelection(AbsTargetSelectionTable):
    step = 1
    number_of_steps = 4
    title = "SELECT RECEPTORS"
    description = "Select receptors in the table (below) or browse the classification tree (right). You can select entire" \
        + " families or individual receptors.\n\nOnce you have selected all your receptors, click the green button."
    docs = "sequences.html#structure-based-alignments"
    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])
    buttons = {
        "continue": {
            "label": "Next",
            "onclick": "submitSelection('/seqsign/savepos');",
            "color": "success",
        },
    }

# class PosTargetSelection(AbsTargetSelection):
#     step = 1
#     number_of_steps = 4
#     docs = 'sequences.html#structure-based-alignments'
#     selection_boxes = OrderedDict([
#         ('reference', False),
#         ('targets', True),
#         ('segments', False),
#     ])
#     buttons = {
#         'continue': {
#             'label': 'Continue to next step',
#             'url': '/seqsign/savepos',
#             'color': 'success',
#         },
#     }


def preserve_targets(request):

    request.session['targets_pos'] = deepcopy(request.session.get('selection', False))
    # get simple selection from session
    simple_selection = request.session.get('selection', False)

    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)
        selection.clear('targets')
        # export simple selection that can be serialized
        simple_selection = selection.exporter()
        # add simple selection to session
        request.session['selection'] = simple_selection

    return redirect('/seqsign/negativegroupselection',)


class NegTargetSelection(AbsTargetSelectionTable):

    default_species = 'Human'
    step = 2
    number_of_steps = 4
    title = "SELECT RECEPTORS"
    description = "Select receptors in the table (below) or browse the classification tree (right). You can select entire" \
        + " families or individual receptors.\n\nOnce you have selected all your receptors, click the green button."
    docs = "sequences.html#structure-based-alignments"
    selection_boxes = OrderedDict([
        ("reference", False),
        ("targets", True),
        ("segments", False),
    ])
    buttons = {
        "continue": {
            "label": "Next",
            "onclick": "submitSelection('/seqsign/segmentselectionsignature');",
            "color": "success",
        },
    }

# class NegTargetSelection(AbsTargetSelection):
#
#     default_species = 'Human'
#     step = 2
#     number_of_steps = 4
#     docs = 'sequences.html#structure-based-alignments'
#     selection_boxes = OrderedDict([
#         ('reference', False),
#         ('targets', True),
#         ('segments', False),
#     ])
#     buttons = {
#         'continue': {
#             'label': 'Continue to next step',
#             'url': '/seqsign/segmentselectionsignature',
#             'color': 'success',
#         },
#     }


class SegmentSelectionSignature(AbsSegmentSelection):
    step = 3
    number_of_steps = 4

    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', False),
        ('segments', True),
    ])
    buttons = {
        'continue': {
            'label': 'Calculate sequence signature',
            'url': '/seqsign/render_signature',
            'color': 'success',
        },
    }

def render_reordered(request, group):

    #grab the selections from session data
    #targets set #1
    ss_pos = request.session.get('targets_pos', False)
    #targets set #2
    ss_neg = request.session.get('selection', False)

    aln = Alignment()

    if group == 'positive':
        aln.load_proteins_from_selection(ss_pos)
    elif group == 'negative':
        aln.load_proteins_from_selection(ss_neg)

    aln.load_segments_from_selection(ss_neg)
    aln.build_alignment()
    aln.calculate_statistics()
    return render(request, 'seqsign/alignment_reordered.html', context={
        'aln': aln,
        'num_residue_columns': len(aln.positions) + len(aln.segments)
        })


def render_signature(request):

    # grab the selections from session data

    # targets set #1
    ss_pos = request.session.get('targets_pos', False)
    # targets set #2
    ss_neg = request.session.get('selection', False)

    # setup signature
    signature = SequenceSignature()
    signature.setup_alignments_from_selection(ss_pos, ss_neg)
    # calculate the signature
    signature.calculate_signature()

    # calculate the Z-scores signatures
    signature.calculate_zscales_signature()

    # save for later
    # signature_map = feats_delta.argmax(axis=0)
    request.session['signature'] = signature.prepare_session_data()
    request.session.modified = True

    signature.segment_headers_to_hide = [seg for seg, posis in signature.common_segments.items() if len(posis)==0]

    return_html = render(
        request,
        'sequence_signature.html',
        signature.prepare_display_data()
        )

    return return_html

def render_signature_excel(request):

    # version #2 - 5 sheets with separate pieces of signature outline

    # step 1 - repeat the data preparation for a sequence signature

    # targets set #1
    ss_pos = request.session.get('targets_pos', False)
    # targets set #2
    ss_neg = request.session.get('selection', False)

    signature = SequenceSignature()
    signature.setup_alignments_from_selection(ss_pos, ss_neg)

    # calculate the signture
    signature.calculate_signature()
    signature.calculate_zscales_signature()

    outstream = BytesIO()
    # wb = xlsxwriter.Workbook('excel_test.xlsx', {'in_memory': False})
    wb = xlsxwriter.Workbook(outstream, {'in_memory': True})
    # Feature stats for signature
    signature.prepare_excel_worksheet(
        wb,
        'signature_properties',
        'signature',
        'features'
    )
    # Signature Z-scales
    signature.zscales_excel(
        wb,
        "signature_zscales",
        'signature'
    )
    # Feature stats for positive group alignment
    signature.prepare_excel_worksheet(
        wb,
        'protein_set1_properties',
        'positive',
        'features'
    )
    # Positive group alignment
    signature.prepare_excel_worksheet(
        wb,
        'protein_set1_aln',
        'positive',
        'alignment'
    )
    # Positive group Z-scales
    signature.zscales_excel(
        wb,
        "protein_set1_zscales",
        'positive'
    )
    # Feature stats for negative group alignment
    signature.prepare_excel_worksheet(
        wb,
        'protein_set2_properties',
        'negative',
        'features'
    )
    # Negative group alignment
    signature.prepare_excel_worksheet(
        wb,
        'protein_set2_aln',
        'negative',
        'alignment'
    )
    # Negative group Z-scales
    signature.zscales_excel(
        wb,
        "protein_set2_zscales",
        'negative'
    )
    signature.per_gn_signature_excel(wb)

    wb.close()
    outstream.seek(0)
    response = HttpResponse(
        outstream.read(),
        content_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )
    response['Content-Disposition'] = "attachment; filename=sequence_signature.xlsx"

    return response

def render_signature_match_scores(request, cutoff):

    signature_data = request.session.get('signature')

    # targets set #1
    ss_pos = request.session.get('targets_pos', False)
    # targets set #2
    ss_neg = request.session.get('selection', False)

    signature_match = SignatureMatch(
        signature_data['common_positions'],
        signature_data['numbering_schemes'],
        signature_data['common_segments'],
        signature_data['diff_matrix'],
        get_proteins_from_selection(ss_pos),
        get_proteins_from_selection(ss_neg),
        cutoff = int(cutoff)
    )
    signature_match.score_protein_class(get_proteins_from_selection(ss_pos)[0].family.slug[:3])

    # Temp disabling Excel due to issue with serializing SeqSign object
    # request.session['signature_match'] = {
    #     'scores': dict(signature_match.protein_report),
    #     'scores_pos': dict(signature_match.scores_pos),
    #     'scores_neg': dict(signature_match.scores_neg),
    #     'protein_signatures': dict(signature_match.protein_signatures),
    #     'signatures_pos': dict(signature_match.signatures_pos),
    #     'signatures_neg': dict(signature_match.signatures_neg),
    #     'signature_filtered': dict(signature_match.signature_consensus),
    #     'relevant_gn': dict(signature_match.relevant_gn),
    #     'relevant_segments': dict(signature_match.relevant_segments),
    #     'numbering_schemes': signature_match.schemes,
    # }

    response = render(
        request,
        'signature_match.html',
        {'scores': signature_match}
        )
    return response

def render_signature_match_excel(request):

    scores_data = request.session.get('signature_match', False)

    outstream = BytesIO()
    wb = xlsxwriter.Workbook(outstream, {'in_memory': True})

    signature_score_excel(
        wb,
        scores_data['scores'],
        scores_data['protein_signatures'],
        scores_data['signature_filtered'],
        scores_data['relevant_gn'],
        scores_data['relevant_segments'],
        scores_data['numbering_schemes'],
        scores_data['scores_pos'],
        scores_data['scores_neg'],
        scores_data['signatures_pos'],
        scores_data['signatures_neg'],

    )
    wb.close()
    outstream.seek(0)
    response = HttpResponse(
        outstream.read(),
        content_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )
    response['Content-Disposition'] = "attachment; filename=sequence_signature_protein_scores.xlsx"

    return response
