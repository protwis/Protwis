"""
A module for generating sequence signatures for the given two sets of proteins.
"""
from django.conf import settings
from django.core import exceptions

from alignment.functions import strip_html_tags, get_format_props
Alignment = getattr(__import__(
    'common.alignment_' + settings.SITE_NAME,
    fromlist=['Alignment']
    ), 'Alignment')

from common.definitions import AMINO_ACIDS, AMINO_ACID_GROUPS, AMINO_ACID_GROUP_NAMES
from protein.models import Protein, ProteinConformation
from residue.models import Residue


from collections import OrderedDict
from copy import deepcopy
import numpy as np
import re
import time

class SequenceSignature:
    """
    A class handling the sequence signature.
    """

    def __init__(self):

        self.aln_pos = Alignment()
        self.aln_neg = Alignment()

        self.features_normalized_pos = OrderedDict()
        self.features_normalized_neg = OrderedDict()
        self.features_frequency_difference = OrderedDict()
        self.features_frequency_diff_display = []

        self.freq_cutoff = 30
        self.common_gn = OrderedDict()

    def setup_alignments(self, segments, protein_set_positive = None, protein_set_negative = None):

        if protein_set_positive:
            self.aln_pos.load_proteins(protein_set_positive)
        if protein_set_negative:
            self.aln_neg.load_proteins(protein_set_negative)

        # In case positive and negative sets come from different classes
        # unify the numbering schemes
        self.common_schemes = self.merge_numbering_schemes()
        self.aln_pos.numbering_schemes = self.common_schemes
        self.aln_neg.numbering_schemes = self.common_schemes
        # now load the segments and generic numbers
        self.aln_pos.load_segments(segments)
        self.aln_neg.load_segments(segments)

        self.aln_pos.build_alignment()
        self.aln_neg.build_alignment()

        self.common_gn = deepcopy(self.aln_pos.generic_numbers)
        for scheme in self.aln_neg.numbering_schemes:
            for segment in self.aln_neg.segments:
                for pos in self.aln_neg.generic_numbers[scheme[0]][segment].items():
                    if pos[0] not in self.common_gn[scheme[0]][segment].keys():
                        self.common_gn[scheme[0]][segment][pos[0]] = pos[1]
                self.common_gn[scheme[0]][segment] = OrderedDict(sorted(
                    self.common_gn[scheme[0]][segment].items(),
                    key=lambda x: x[0].split('x')
                    ))
        self.common_segments = dict([
            (x, sorted(list(set(self.aln_pos.segments[x]) | set(self.aln_neg.segments[x])), key=lambda x: x.split('x'))) for x in self.aln_neg.segments
        ])
        # tweaking alignment
        self._update_alignment(self.aln_pos)
        self.aln_pos.calculate_statistics()
        # tweaking consensus seq
        self._update_consensus_sequence(self.aln_pos)

        # tweaking negative alignment
        self._update_alignment(self.aln_neg)
        self.aln_neg.calculate_statistics()
        # tweaking consensus seq
        self._update_consensus_sequence(self.aln_neg)

    def _update_alignment(self, alignment):

        for prot in alignment.proteins:
            for seg, resi in prot.alignment.items():
                consensus = []
                aln_list = [x[0] for x in resi]
                aln_dict = dict([
                    (x[0], x) for x in resi
                ])
                for pos in self.common_segments[seg]:
                    if pos not in aln_list:
                        consensus.append([pos, False, '-', 0])
                    else:
                        consensus.append(aln_dict[pos])
                prot.alignment[seg] = consensus

    def _update_consensus_sequence(self, alignment):

        for seg, resi in alignment.consensus.items():
            consensus = OrderedDict()
            aln_list = [x for x in resi.keys()]
            aln_dict = dict([
                (x, resi[x]) for x in resi.keys()
            ])
            for pos in self.common_segments[seg]:
                if pos not in aln_list:
                    consensus[pos] = ['_', 0, 100]
                else:
                    consensus[pos] = aln_dict[pos]
            alignment.consensus[seg] = consensus

    def _convert_feature_stats(self, fstats, aln):

        tmp_fstats = []
        for row in range(len(AMINO_ACID_GROUPS.keys())):
            tmp_row = []
            print(aln.feature_stats[row])
            for segment in self.common_segments:
                print(fstats[segment][row])
                tmp_row.append([[
                    str(x),
                    str(int(x/10)),
                ] for x in fstats[segment][row]])
            print(tmp_row)
            tmp_fstats.append(tmp_row)
        aln.feature_stats = tmp_fstats

    def setup_alignments_from_selection(self, positive_selection, negative_selection):
        """
        The function gathers necessary information from provided selections
        and runs the calculations of the sequence alignments independently for
        both protein sets. It also finds the common set of residue positions.

        Arguments:
            positive_selection {Selection} -- selection containing first group of proteins
            negative_selection {[type]} -- selction containing second group of proteins along with the user-selcted sequence segments for the alignment
        """

        self.aln_pos.load_proteins_from_selection(positive_selection)
        self.aln_neg.load_proteins_from_selection(negative_selection)

        # local segment list
        segments = []

        # read selection
        for segment in negative_selection.segments:
            segments.append(segment)

        self.setup_alignments(segments)

    def calculate_signature(self):
        """
        Calculates the feature frequency difference between two protein sets.
        Generates the full differential matrix as well as maximum difference for a position (for scatter plot).
        """
        for sid, segment in enumerate(self.aln_neg.segments):
            self.features_normalized_pos[segment] = np.array(
                [[x[0] for x in feat[sid]] for feat in self.aln_pos.feature_stats],
                dtype='int'
                )
            self.features_normalized_neg[segment] = np.array(
                [[x[0] for x in feat[sid]] for feat in self.aln_neg.feature_stats],
                dtype='int'
                )

        for segment in self.aln_neg.segments:
            #TODO: get the correct default numering scheme from settings
            for idx, res in enumerate(self.common_gn[self.common_schemes[0][0]][segment].keys()):
                if res not in self.aln_pos.generic_numbers[self.common_schemes[0][0]][segment].keys():
                    self.features_normalized_pos[segment] = np.insert(self.features_normalized_pos[segment], idx, 0, axis=1)
                    # Set 100% occurence for a gap feature
                    self.features_normalized_pos[segment][-1, idx] = 100
                elif res not in self.aln_neg.generic_numbers[self.common_schemes[0][0]][segment].keys():
                    self.features_normalized_neg[segment] = np.insert(self.features_normalized_neg[segment], idx, 0, axis=1)
                    # Set 100% occurence for a gap feature
                    self.features_normalized_neg[segment][-1, idx] = 100

            # now the difference
            self.features_frequency_difference[segment] = np.subtract(
                self.features_normalized_pos[segment],
                self.features_normalized_neg[segment]
                )

        self._convert_feature_stats(self.features_normalized_pos, self.aln_pos)
        self._convert_feature_stats(self.features_normalized_neg, self.aln_neg)

        # Version with display data
        for row in range(len(AMINO_ACID_GROUPS.keys())):
            tmp_row = []
            for segment in self.aln_neg.segments:
                #first item is the real value,
                # second is the assignmnent of color (via css)
                # 0 - red, 5 - yellow, 10 - green
                #third item is a tooltip
                tmp_row.append([[
                    x,
                    int(x/20)+5,
                    "{} - {}".format(
                        self.features_normalized_pos[segment][row][y],
                        self.features_normalized_neg[segment][row][y]
                        )
                    ] for y, x in enumerate(self.features_frequency_difference[segment][row])])
            self.features_frequency_diff_display.append(tmp_row)

        self.signature = OrderedDict([(x, []) for x in self.aln_neg.segments])
        for segment in self.aln_neg.segments:
            tmp = np.array(self.features_frequency_difference[segment])
            signature_map = np.absolute(tmp).argmax(axis=0)
            self.signature[segment] = []
            for col, pos in enumerate(list(signature_map)):
                self.signature[segment].append([
                    list(AMINO_ACID_GROUPS.keys())[pos],
                    list(AMINO_ACID_GROUP_NAMES.values())[pos],
                    self.features_frequency_difference[segment][pos][col],
                    int(self.features_frequency_difference[segment][pos][col]/20)+5
                ])

        features_pos = OrderedDict()
        features_neg = OrderedDict()
        self.features_consensus_pos = OrderedDict([(x, []) for x in self.aln_neg.segments])
        self.features_consensus_neg = OrderedDict([(x, []) for x in self.aln_neg.segments])
        for sid, segment in enumerate(self.aln_neg.segments):
            features_pos[segment] = np.array(
                [[x[0] for x in feat[sid]] for feat in self.aln_pos.feature_stats],
                dtype='int'
                )
            features_neg[segment] = np.array(
                [[x[0] for x in feat[sid]] for feat in self.aln_neg.feature_stats],
                dtype='int'
                )
            features_cons_pos = np.absolute(features_pos[segment]).argmax(axis=0)
            features_cons_neg = np.absolute(features_neg[segment]).argmax(axis=0)

            for col, pos in enumerate(list(features_cons_pos)):
                self.features_consensus_pos[segment].append([
                    list(AMINO_ACID_GROUPS.keys())[pos],
                    list(AMINO_ACID_GROUP_NAMES.values())[pos],
                    features_pos[segment][pos][col],
                    int(features_pos[segment][pos][col]/20)+5
                ])
            for col, pos in enumerate(list(features_cons_neg)):
                self.features_consensus_neg[segment].append([
                    list(AMINO_ACID_GROUPS.keys())[pos],
                    list(AMINO_ACID_GROUP_NAMES.values())[pos],
                    features_neg[segment][pos][col],
                    int(features_neg[segment][pos][col]/20)+5
                ])
        self._convert_feature_stats(self.features_normalized_pos, self.aln_pos)
        self._convert_feature_stats(self.features_normalized_neg, self.aln_neg)

    def prepare_display_data(self):

        options = {
            'num_residue_columns': len(sum([[x for x in self.common_gn[self.common_schemes[0][0]][segment]] for segment in self.aln_neg.segments], [])),
            'num_of_sequences_pos': len(self.aln_pos.proteins),
            'num_residue_columns_pos': len(self.aln_pos.positions),
            'num_of_sequences_neg': len(self.aln_neg.proteins),
            'num_residue_columns_neg': len(self.aln_neg.positions),
            'common_segments': self.common_segments,
            'common_generic_numbers': self.common_gn,
            'feats_signature': self.features_frequency_diff_display,
            'signature_consensus': self.signature,
            'feats_cons_pos': self.features_consensus_pos,
            'feats_cons_neg': self.features_consensus_neg,
            'a_pos': self.aln_pos,
            'a_neg': self.aln_neg,
        }
        return options

    def prepare_session_data(self):

        session_signature = {
            'common_positions': self.common_gn,
            'diff_matrix': self.features_frequency_difference,
            'numbering_schemes': self.common_schemes,
            'segments': self.aln_neg.segments
        }
        return session_signature

    def merge_numbering_schemes(self):
        """
        Extract all of the numbering schemes used for a set of proteins.

        Arguments:
            proteins {selection} -- A set of proteins to analyze
        """

        numbering_schemes = {}
        for prot in self.aln_pos.proteins + self.aln_neg.proteins:
            if prot.protein.residue_numbering_scheme.slug not in numbering_schemes:
                rnsn = prot.protein.residue_numbering_scheme.name
                numbering_schemes[prot.protein.residue_numbering_scheme.slug] = rnsn
        # order and convert numbering scheme dict to tuple
        return sorted(numbering_schemes.items(), key=lambda x: x[0])

    def prepare_excel_worksheet(self, workbook, worksheet_name, aln='positive', data='alignment'):
        """
        A function saving alignment data subset into the excel spreadsheet.
        It adds a worksheet to an existing workbook and saves only a selected subset of alignment data.
        For a complete save of the alignment it needs to be wrapped with additional code.

        The outline of the excel worksheet is similar to the one of html page.
        First column shows nunbering schemes, protein list, etc
        The frequency data start from column B

        Arguments:
            workbook {xlrsxwriter.Workbook} -- object to add workseet to
            worksheet_name {string} -- name for the new workseet

        Keyword Arguments:
            alignment {string} -- alignment to extract data from.
                                    Possible choices: positive, negative, signature
            data {string} -- data type to save to workshet: 'alignment' or 'features' frequencies
        """

        props = AMINO_ACID_GROUP_NAMES.values()
        worksheet = workbook.add_worksheet(worksheet_name)

        if aln == 'positive':
            numbering_schemes = self.aln_pos.numbering_schemes
            generic_numbers_set = self.aln_pos.generic_numbers
            alignment = self.aln_pos
            if data == 'features':
                data_block = self.aln_pos.feature_stats
        elif aln == 'negative':
            numbering_schemes = self.aln_neg.numbering_schemes
            generic_numbers_set = self.aln_neg.generic_numbers
            alignment = self.aln_neg
            if data == 'features':
                data_block = self.aln_neg.feature_stats
        else:
            numbering_schemes = self.common_schemes
            generic_numbers_set = self.common_gn
            if data == 'features':
                data_block = self.features_frequency_diff_display

        # First column, numbering schemes
        for row, scheme in enumerate(numbering_schemes):
            worksheet.write(1 + 3*row, 0, scheme[1])

        # First column, stats
        if data == 'features':
            for offset, prop in enumerate(props):
                worksheet.write(1 + 3 * len(numbering_schemes) + offset, 0, prop)

        # First column, protein list (for alignment) and line for consensus sequence
        else:
            for offset, prot in enumerate(alignment.proteins):
                worksheet.write(
                    1 + 3 * len(numbering_schemes) + offset,
                    0,
                    prot.protein.entry_name
                )
            worksheet.write(
                1 + len(numbering_schemes) + len(alignment.proteins),
                0,
                'CONSENSUS'
                )

        # Second column and on
        # Segments
        offset = 0
        for segment in generic_numbers_set[numbering_schemes[0][0]].keys():
            worksheet.merge_range(
                0,
                1 + offset,
                0,
                len(generic_numbers_set[numbering_schemes[0][0]][segment]) + offset - 1,
                segment
            )
            offset += len(generic_numbers_set[numbering_schemes[0][0]][segment])

        # Generic numbers
        # for row, item in enumerate(generic_numbers_set.items()):
        for row, item in enumerate(numbering_schemes):
            scheme = item[0]
            offset = 1
            for sn, gn_list in generic_numbers_set[scheme].items():
                for col, gn_pair in enumerate(gn_list.items()):
                    try:
                        tm, bw, gpcrdb = re.split('\.|x', strip_html_tags(gn_pair[1]))
                    except:
                        tm, bw, gpcrdb = ('', '', '')
                    worksheet.write(
                        1 + 3 * row,
                        col + offset,
                        tm
                    )
                    worksheet.write(
                        2 + 3 * row,
                        col + offset,
                        bw
                    )
                    worksheet.write(
                        3 + 3*row,
                        col + offset,
                        gpcrdb
                    )
                offset += len(gn_list.items())

        # Stats
        if data == 'features':
            offset = 1 + 3 * len(numbering_schemes)

            for row, prop in enumerate(data_block):
                col_offset = 0
                for segment in prop:
                    for col, freq in enumerate(segment):
                        cell_format = workbook.add_format(get_format_props(freq[1]))
                        worksheet.write(
                            offset + row,
                            1 + col + col_offset,
                            freq[0] if isinstance(freq[0], int) else int(freq[0]),
                            cell_format
                        )
                    col_offset += len(segment)
            col_offset = 0
            for segment, cons_feat in self.signature.items():
                for col, chunk in enumerate(cons_feat):
                    worksheet.write(
                        offset + len(AMINO_ACID_GROUPS),
                        1 + col + col_offset,
                        chunk[0]
                    )
                    cell_format = workbook.add_format(get_format_props(int(chunk[2]/20)+5))
                    worksheet.write(
                        1 + offset + len(AMINO_ACID_GROUPS),
                        1 + col + col_offset,
                        chunk[2],
                        cell_format
                    )
                col_offset += len(cons_feat)
        # Alignment
        else:
            offset = 1 + 3 * len(alignment.numbering_schemes)

            for row, data in enumerate(alignment.proteins):
                col_offset = 0
                for segment, sequence in data.alignment.items():
                    for col, res in enumerate(sequence):
                        cell_format = workbook.add_format(get_format_props(res=res[2]))
                        worksheet.write(
                            offset + row,
                            1 + col + col_offset,
                            res[2],
                            cell_format
                        )
                    col_offset += len(sequence)
            # Consensus sequence
            row = 1 + 3 * len(alignment.numbering_schemes) + len(alignment.proteins)
            col_offset = 0
            for segment, sequence in alignment.consensus.items():
                for col, data in enumerate(sequence.items()):
                    res = data[1]
                    cell_format = workbook.add_format(get_format_props(res=res[0]))
                    worksheet.write(
                        row,
                        1 + col + col_offset,
                        res[0],
                        cell_format
                    )
                    cell_format = workbook.add_format(get_format_props(res[1]))
                    worksheet.write(
                        row + 1,
                        1 + col + col_offset,
                        res[2],
                        cell_format
                    )

                col_offset += len(sequence.items())


class SignatureMatch():

    def __init__(self, common_positions, numbering_schemes, segments, difference_matrix, protein_set, cutoff=40):

        self.cutoff = cutoff
        self.common_gn = common_positions
        self.schemes = numbering_schemes
        self.segments = segments
        self.diff_matrix = difference_matrix
        self.signature_filtered = OrderedDict()
        self.protein_set = protein_set
        self.relevant_gn = OrderedDict([(x[0], OrderedDict()) for x in self.schemes])
        self.relevant_segments = []

        self.find_relevant_gns()

        self.residue_to_feat = dict(
            [(x, set()) for x in AMINO_ACIDS.keys()]
            )
        for fidx, feat in enumerate(AMINO_ACID_GROUPS.items()):
            for res in feat[1].split(','):
                self.residue_to_feat[res].add(fidx)


    def find_relevant_gns(self):

        signature_consensus = OrderedDict()
        for segment in self.segments:
            segment_consensus = []
            signature_map = np.absolute(self.diff_matrix[segment]).argmax(axis=0)
            for col, pos in enumerate(list(signature_map)):
                if self.diff_matrix[segment][pos][col] > self.cutoff:
                    segment_consensus.append(self.diff_matrix[segment][ : , col])
                    for scheme in self.schemes:
                        gnum = list(self.common_gn[scheme[0]][segment].items())[col]
                        try:
                            self.relevant_gn[scheme[0]][segment][gnum[0]] = gnum[1]
                        except:
                            self.relevant_gn[scheme[0]][segment] = OrderedDict()
                            self.relevant_gn[scheme[0]][segment][gnum[0]] = gnum[1]

            segment_consensus = np.array(segment_consensus).T
            if segment_consensus != []:
                signature_consensus[segment] = segment_consensus
        self.signature_filtered = signature_consensus
        self.relevant_segments = signature_consensus.keys()


    def score_protein_class(self, pclass_slug='001'):

        a = time.time()
        protein_scores = {}
        class_proteins = Protein.objects.filter(
            species__common_name='Human',
            family__slug__startswith=pclass_slug
            ).exclude(
                id__in=[x.id for x in self.protein_set]
                )
        class_a_pcf = ProteinConformation.objects.order_by('protein__family__slug',
            'protein__entry_name').filter(protein__in=class_proteins, protein__sequence_type__slug='wt').exclude(protein__entry_name__endswith='-consensus')
        for pcf in class_a_pcf:
            start = time.time()
            protein_scores[pcf] = self.score_protein(pcf)
            end = time.time()
            print("Time elapsed for {}: ".format(pcf.protein.entry_name), end - start)
        b = time.time()
        print("Total time: ", b - a)
        return sorted(protein_scores.items(), key=lambda x: x[1], reverse=True)


    def score_protein(self, pcf):

        prot_score = 0.0
        for segment in self.relevant_segments:
            resi = Residue.objects.filter(
                protein_segment__slug=segment,
                protein_conformation=pcf,
                generic_number__label__in=self.relevant_gn[self.schemes[0][0]][segment].keys(),
                )
            for idx, pos in enumerate(self.relevant_gn[self.schemes[0][0]][segment].keys()):
                dslice = self.diff_matrix[segment][ : , idx]
                try:
                    res = resi.get(generic_number__label=pos).amino_acid
                    n = set(dslice.nonzero)
                    p = self.residue_to_feat[res]
                    prot_score += np.sum(dslice[p & n])
                    prot_score -= np.sum(dslice[n - p])
                except Exception as e:
                    prot_score -= np.sum(dslice)
        return prot_score/100

    def prepare_session_data(self):
        session_signature = {
            'common_positions': self.common_gn,
            'diff_matrix': self.diff_matrix,
            'numbering_schemes': self.schemes,
            'segments': self.segments,
            'relevant_gn': self.relevant_gn,
            'relevant_segments': self.relevant_segments,
            'signature_filtered': self.signature_filtered,
        }
        return session_signature

class ScoreBreakdown():

    def __init__(self, pcf, cutoff, **kwargs):

        self.pcf = pcf
        self.cutoff = cutoff
        self.common_gn = kwargs['common_positions']
        self.schemes = kwargs['numbering_schemes']
        self.segments = kwargs['segments']
        self.diff_matrix = kwargs['diff_matrix']
        self.relevant_gn = kwargs['relevant_gn']
        self.relevant_segments = kwargs['relevant_segments']
        self.score_breakdown = kwargs['signature_filtered']
        # self.score_breakdown = OrderedDict()
        # self.relevant_gn = OrderedDict([(x[0], OrderedDict()) for x in self.schemes])

    def find_relevant_gns(self):

        signature_consensus = OrderedDict()
        for segment in self.segments:
            segment_consensus = []
            signature_map = np.absolute(self.diff_matrix[segment]).argmax(axis=0)
            for col, pos in enumerate(list(signature_map)):
                if self.diff_matrix[segment][pos][col] > self.cutoff:
                    segment_consensus.append(self.diff_matrix[segment][ : , col])
                    for scheme in self.schemes:
                        gnum = list(self.common_gn[scheme[0]][segment].items())[col]
                        try:
                            self.relevant_gn[scheme[0]][segment][gnum[0]] = gnum[1]
                        except:
                            self.relevant_gn[scheme[0]][segment] = OrderedDict()
                            self.relevant_gn[scheme[0]][segment][gnum[0]] = gnum[1]

            segment_consensus = np.array(segment_consensus).T
            if segment_consensus != []:
                signature_consensus[segment] = segment_consensus
        self.score_breakdown = signature_consensus

    def prepare_display_data(self):

        display_score_breakdown = []
        for fidx, feat in enumerate(AMINO_ACID_GROUPS.items()):
            display_score_breakdown.append([])
            for segment in self.score_breakdown.keys():
                tmp = []
                rs = Residue.objects.filter(
                    protein_segment__slug=segment,
                    protein_conformation__id = self.pcf
                    )
                for idx, pos in enumerate(self.relevant_gn[self.schemes[0][0]][segment].keys()):
                    try:
                        res = rs.get(generic_number__label=pos)
                    except (exceptions.ObjectDoesNotExist, exceptions.MultipleObjectsReturned):
                        res = None
                    val = self.score_breakdown[segment][fidx][idx]
                    if not res:
                        tmp.append([-val, "red"] if val != 0 else [val, "white"])
                        continue
                    if res.amino_acid in feat[1]:
                        if val != 0:
                            tmp.append([val, "green"])
                        else:
                            tmp.append([val, "white"])
                    else:
                        if val != 0:
                            tmp.append([-val, "red"])
                        else:

                            tmp.append([val, "white"])
                display_score_breakdown[fidx].append(tmp)
        signature_consensus = OrderedDict([
            (x,
            np.sum(
                self.score_breakdown[x],
                axis=0)/100
            ) for x in self.score_breakdown.keys()
            ])
        display_data = {
            'common_numbering_schemes': self.schemes,
            'common_generic_numbers': self.relevant_gn,
            'common_segments': self.score_breakdown.keys(),
            'feats_signature': display_score_breakdown,
            'signature_consensus': signature_consensus,

            'num_residue_columns': OrderedDict(
                [(x, len(signature_consensus[x])) for x in self.score_breakdown.keys()]
                )
        }
        return display_data
