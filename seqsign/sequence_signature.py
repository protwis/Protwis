"""
A module for generating sequence signatures for the given two sets of proteins.
"""
from django.conf import settings
#from django.core import exceptions

from alignment.functions import strip_html_tags, get_format_props, prepare_aa_group_preference, zscales_color_scale
Alignment = getattr(__import__(
    'common.alignment_' + settings.SITE_NAME,
    fromlist=['Alignment']
    ), 'Alignment')

from common.definitions import AA_ZSCALES, AMINO_ACIDS, AMINO_ACID_GROUPS, AMINO_ACID_GROUP_NAMES, AMINO_ACID_GROUP_PROPERTIES, ZSCALES
from protein.models import Protein, ProteinConformation
from residue.models import Residue
from signprot.models import SignprotComplex


from collections import OrderedDict
from copy import deepcopy
import numpy as np
from operator import itemgetter
import re
from scipy.stats import t
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

        self.features_consensus_pos = OrderedDict()
        self.features_consensus_neg = OrderedDict()

        self.freq_cutoff = 30
        self.common_gn = OrderedDict()
        self.common_segments = OrderedDict()
        self.common_schemes = {}
        self.selection_segments = OrderedDict()

        self.signature = OrderedDict()

        self.zscales_signature = OrderedDict()

        self.feature_preference = prepare_aa_group_preference()
        self.group_lengths = dict([
            (x, len(y)) for x,y in enumerate(AMINO_ACID_GROUPS.values())
        ])
        self.default_column = np.array([((y == '-') and 100) or 0 for y in AMINO_ACID_GROUPS.keys()])

        self.segment_headers_to_hide = []


    def _assign_preferred_features(self, signature, segment, ref_matrix):

        new_signature = []
        for pos, argmax in enumerate(signature):
            new_signature.append(self._calculate_best_feature(pos, segment, argmax, ref_matrix))
        return new_signature


    def _calculate_best_feature(self, pos, segment, argmax, ref_matrix):

        tmp = self.feature_preference[argmax]
        equiv_feat = np.where(np.isin(ref_matrix[segment][:, pos], ref_matrix[segment][argmax][pos]))[0]
        pref_feat = argmax
        min_len = self.group_lengths[argmax]

        for efeat in equiv_feat:
            if efeat in tmp and self.group_lengths[efeat] < min_len:
                pref_feat = efeat
                min_len = self.group_lengths[efeat]
            # when two features have the same aa count, take the one from positive set
            elif efeat in tmp and self.group_lengths[efeat] == min_len:
                if ref_matrix[segment][pref_feat][pos] < 0 and ref_matrix[segment][efeat][pos] > 0:
                    pref_feat = efeat
        return pref_feat

    def setup_alignments(self, segments, protein_set_positive=None, protein_set_negative=None):
        """Setup (fetch and normalize) the data necessary for calculation of the signature.

        Arguments:
            segments {list} -- List of segments to calculate the signature from

        Keyword Arguments:
            protein_set_positive {list} -- list of Protein objects - a positive (reference) set (default: {None})
            protein_set_negative {list} -- list of Protein objects - a negative set (default: {None})
        """

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

        # SM: "Borrowing" the selected segments from Alignment object to avoid extracting Selection object
        # SM: segments are for now the same for both positive and negative set
        self.selection_segments = self.aln_pos.segments

        self.aln_pos.build_alignment()
        self.aln_neg.build_alignment()

        for x in self.aln_neg.segments.keys():
            if x not in self.aln_pos.segments.keys():
                self.aln_pos.segments[x] = {}
                self.aln_pos.generic_numbers[self.common_schemes[0][0]][x] = {}

        for x in self.aln_pos.segments.keys():
            if x not in self.aln_neg.segments.keys():
                self.aln_neg.segments[x] = {}
                self.aln_pos.generic_numbers[self.common_schemes[0][0]][x] = {}

        self.common_gn = deepcopy(self.aln_pos.generic_numbers)
        for scheme in self.aln_neg.numbering_schemes:
            for segment in self.aln_neg.segments:
                for pos in self.aln_neg.generic_numbers[scheme[0]][segment].items():
                    if segment not in self.common_gn[scheme[0]].keys():
                        self.common_gn[scheme[0]][segment] = {}
                    if pos[0] not in self.common_gn[scheme[0]][segment].keys():
                        self.common_gn[scheme[0]][segment][pos[0]] = pos[1]
                self.common_gn[scheme[0]][segment] = OrderedDict(sorted(
                    self.common_gn[scheme[0]][segment].items(),
                    key=lambda x: x[0].split('x')
                    ))

        self.common_segments = OrderedDict([
            (x, sorted(list(set(self.aln_pos.segments[x]) | set(self.aln_neg.segments[x])), key=lambda x: x.split('x'))) for x in self.aln_neg.segments
        ])

        # tweaking alignment
        self.aln_pos.calculate_statistics()
        self._update_alignment(self.aln_pos)
        # tweaking consensus seq
        self._update_consensus_sequence(self.aln_pos)

        # tweaking negative alignment
        self.aln_neg.calculate_statistics()
        self._update_alignment(self.aln_neg)
        # tweaking consensus seq
        self._update_consensus_sequence(self.aln_neg)

    def setup_alignments_signprot(self, segments, protein_set_positive=None,
            protein_set_negative=None, ignore_in_alignment=None):
        """Setup (fetch and normalize) the data necessary for calculation of the signature.

        Arguments:
            segments {list} -- List of segments to calculate the signature from

        Keyword Arguments:
            protein_set_positive {list} -- list of Protein objects - a positive (reference) set (default: {None})
            protein_set_negative {list} -- list of Protein objects - a negative set (default: {None})
        """


        if protein_set_positive:
            self.aln_pos.load_proteins(protein_set_positive)
        if protein_set_negative:
            self.aln_neg.load_proteins(protein_set_negative)

        # In case positive and negative sets come from different classes
        # unify the numbering schemes
        if protein_set_positive and protein_set_negative:
            self.common_schemes = self.merge_numbering_schemes()
            self.aln_pos.numbering_schemes = self.common_schemes
            self.aln_neg.numbering_schemes = self.common_schemes

        # now load the segments and generic numbers
        if protein_set_positive:
            self.common_schemes = self.merge_numbering_schemes()
            self.aln_pos.load_segments(segments)
            self.aln_pos.build_alignment()
        if protein_set_negative:
            self.common_schemes = self.merge_numbering_schemes()
            self.aln_neg.load_segments(segments)
            self.aln_neg.build_alignment()

        self.common_gn = deepcopy(self.aln_pos.generic_numbers)
        if protein_set_negative:
            for scheme in self.aln_neg.numbering_schemes:
                for segment in self.aln_neg.segments:
                    for pos in self.aln_neg.generic_numbers[scheme[0]][segment].items():
                        if pos[0] not in self.common_gn[scheme[0]][segment].keys():
                            self.common_gn[scheme[0]][segment][pos[0]] = pos[1]
                    self.common_gn[scheme[0]][segment] = OrderedDict(sorted(
                        self.common_gn[scheme[0]][segment].items(),
                        key=lambda x: x[0].split('x')
                        ))
            self.common_segments = OrderedDict([
                (x, sorted(list(set(self.aln_pos.segments[x]) | set(self.aln_neg.segments[x])), key=lambda x: x.split('x'))) for x in self.aln_neg.segments
            ])
        else:
            for scheme in self.aln_pos.numbering_schemes:
                for segment in self.aln_pos.segments:
                    for pos in self.aln_pos.generic_numbers[scheme[0]][segment].items():
                        if pos[0] not in self.common_gn[scheme[0]][segment].keys():
                            self.common_gn[scheme[0]][segment][pos[0]] = pos[1]
                    self.common_gn[scheme[0]][segment] = OrderedDict(sorted(
                        self.common_gn[scheme[0]][segment].items(),
                        key=lambda x: x[0].split('x')
                        ))
            self.common_segments = OrderedDict([
                (x, sorted(list(set(self.aln_pos.segments[x])), key=lambda x:
                    x.split('x'))) for x in self.aln_pos.segments
            ])

        if protein_set_positive:
            # tweaking alignment
            if ignore_in_alignment:
                self.aln_pos.calculate_statistics(ignore_in_alignment)
            else:
                self.aln_pos.calculate_statistics()
            self._update_alignment(self.aln_pos)
            # tweaking consensus seq
            self._update_consensus_sequence(self.aln_pos)

        if protein_set_negative:
            # tweaking negative alignment
            if ignore_in_alignment:
                self.aln_neg.calculate_statistics(ignore_in_alignment)
            else:
                self.aln_neg.calculate_statistics()
            self._update_alignment(self.aln_neg)
            # tweaking consensus seq
            self._update_consensus_sequence(self.aln_neg)

    def _update_alignment(self, alignment):

        # print(alignment.segments.keys())
        for prot in alignment.proteins:
            for seg, resi in prot.alignment.items():
                if seg not in self.common_segments.keys():
                    continue
                consensus = []
                aln_list = [x[0] for x in resi]
                aln_dict = dict([
                    (x[0], x) for x in resi
                ])
                for pos in self.common_segments[seg]:
                    if pos not in aln_list:
                        consensus.append([pos, False, '_', 0])
                    else:
                        consensus.append(aln_dict[pos])
                prot.alignment[seg] = consensus

    def _update_consensus_sequence(self, alignment):

        for seg, resi in alignment.consensus.items():
            if seg not in self.common_segments.keys():
                continue
            consensus = OrderedDict()
            aln_list = [x for x in resi.keys()]
            aln_dict = dict([
                (x, resi[x]) for x in resi.keys()
            ])
            for pos in self.common_segments[seg]:
                if pos not in aln_list:
                    consensus[pos] = ['-', 0, 100]
                else:
                    consensus[pos] = aln_dict[pos]
            alignment.consensus[seg] = consensus

    def _convert_feature_stats(self, fstats, aln):

        tmp_fstats = []
        for row in range(len(AMINO_ACID_GROUPS.keys())):
            tmp_row = []
            for segment in self.common_segments:
                tmp_row.append([[
                    str(x),
                    str(int(x/10)), #if x != 0 else -1,
                ] for x in fstats[segment][row]])
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

    def setup_alignments_from_selection_onesided(self, positive_selection):
        """
        The function gathers necessary information from provided selections
        and runs the calculations of the sequence alignments.

        Arguments:
            positive_selection {Selection} -- selection containing first group of proteins
        """

        self.aln_pos.load_proteins_from_selection(positive_selection)

        # local segment list
        segments = []

        # read selection
        for segment in positive_selection.segments:
            segments.append(segment)

        self.setup_alignments(segments)

    def calculate_signature(self):
        """
        Calculates the feature frequency difference between two protein sets.
        Generates the full differential matrix as well as maximum difference for a position (for scatter plot).
        """

        for sid, segment in enumerate(self.selection_segments):
            if segment not in self.common_segments:
                continue
            self.features_normalized_pos[segment] = np.array(
                [[x[0] for x in feat[sid]] for feat in self.aln_pos.feature_stats],
                dtype='int'
                )
            self.features_normalized_neg[segment] = np.array(
                [[x[0] for x in feat[sid]] for feat in self.aln_neg.feature_stats],
                dtype='int'
                )

        for segment in self.common_segments:
            #TODO: get the correct default numering scheme from settings
            for idx, res in enumerate(self.common_gn[self.common_schemes[0][0]][segment].keys()):
                if res not in self.aln_pos.generic_numbers[self.common_schemes[0][0]][segment].keys():
                    self.features_normalized_pos[segment] = np.insert(self.features_normalized_pos[segment], idx, self.default_column, axis=1)
                elif res not in self.aln_neg.generic_numbers[self.common_schemes[0][0]][segment].keys():
                    self.features_normalized_neg[segment] = np.insert(self.features_normalized_neg[segment], idx, self.default_column, axis=1)

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
                    int(x/20)+5 if x!= 0 else -1,
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
            #signature_map = tmp.argmax(axis=0)
            # Update mapping to prefer features with fewer amino acids
            signature_map = self._assign_preferred_features(signature_map, segment, self.features_frequency_difference)

            self.signature[segment] = []
            for col, pos in enumerate(list(signature_map)):
                self.signature[segment].append([
                    list(AMINO_ACID_GROUP_PROPERTIES.values())[pos]['display_name_short'] if self.features_frequency_difference[segment][pos][col] > 0 else '-' + list(AMINO_ACID_GROUP_PROPERTIES.values())[pos]['display_name_short'], # latest implementation of NOT... properties
                    list(AMINO_ACID_GROUP_NAMES.values())[pos] if self.features_frequency_difference[segment][pos][col] > 0 else "Not " + list(AMINO_ACID_GROUP_NAMES.values())[pos], # latest implementation of NOT... properties
                    self.features_frequency_difference[segment][pos][col],
                    int(self.features_frequency_difference[segment][pos][col]/20)+5,
                    list(AMINO_ACID_GROUP_PROPERTIES.values())[pos]['length'],
                    list(AMINO_ACID_GROUPS.keys())[pos],
                    list(AMINO_ACID_GROUP_PROPERTIES.values())[pos]['display_name_short']
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
            features_cons_pos = features_pos[segment].argmax(axis=0)
            features_cons_pos = self._assign_preferred_features(features_cons_pos, segment, features_pos)
            features_cons_neg = features_neg[segment].argmax(axis=0)
            features_cons_neg = self._assign_preferred_features(features_cons_neg, segment, features_neg)

            for col, pos in enumerate(list(features_cons_pos)):
                self.features_consensus_pos[segment].append([
                    # list(AMINO_ACID_GROUPS.keys())[pos],
                    list(AMINO_ACID_GROUP_PROPERTIES.values())[pos]['display_name_short'],
                    list(AMINO_ACID_GROUP_NAMES.values())[pos],
                    features_pos[segment][pos][col],
                    int(features_pos[segment][pos][col]/20)+5,
                    list(AMINO_ACID_GROUP_PROPERTIES.values())[pos]['length'],
                    list(AMINO_ACID_GROUPS.keys())[pos]
                ])
            for col, pos in enumerate(list(features_cons_neg)):
                self.features_consensus_neg[segment].append([
                    # list(AMINO_ACID_GROUPS.keys())[pos],
                    list(AMINO_ACID_GROUP_PROPERTIES.values())[pos]['display_name_short'],
                    list(AMINO_ACID_GROUP_NAMES.values())[pos],
                    features_neg[segment][pos][col],
                    int(features_neg[segment][pos][col]/20)+5,
                    list(AMINO_ACID_GROUP_PROPERTIES.values())[pos]['length'],
                    list(AMINO_ACID_GROUPS.keys())[pos]
                ])
        self._convert_feature_stats(self.features_normalized_pos, self.aln_pos)
        self._convert_feature_stats(self.features_normalized_neg, self.aln_neg)

    def calculate_signature_onesided(self):
        """
        Calculates the feature frequency for one protein set.
        Generates the full differential matrix as well as maximum difference for a position (for scatter plot).
        """
        for sid, segment in enumerate(self.aln_pos.segments):
            self.features_normalized_pos[segment] = np.array(
                [[x[0] for x in feat[sid]] for feat in self.aln_pos.feature_stats],
                dtype='int'
                )

        for segment in self.aln_pos.segments:
            #TODO: get the correct default numering scheme from settings
            for idx, res in enumerate(self.common_gn[self.common_schemes[0][0]][segment].keys()):
                if res not in self.aln_pos.generic_numbers[self.common_schemes[0][0]][segment].keys():
                    self.features_normalized_pos[segment] = np.insert(self.features_normalized_pos[segment], idx, self.default_column, axis=1)

            # now the difference
            self.features_frequency_difference[segment] = np.subtract(
                self.features_normalized_pos[segment],
                0
                )

        self._convert_feature_stats(self.features_normalized_pos, self.aln_pos)

        # Version with display data
        for row in range(len(AMINO_ACID_GROUPS.keys())):
            tmp_row = []
            for segment in self.aln_pos.segments:
                #first item is the real value,
                # second is the assignmnent of color (via css)
                # 0 - red, 5 - yellow, 10 - green
                #third item is a tooltip
                tmp_row.append([[
                    x,
                    int(x/20)+5 if x!= 0 else -1,
                    "{} - {}".format(
                        self.features_normalized_pos[segment][row][y],
                        0
                        )
                    ] for y, x in enumerate(self.features_frequency_difference[segment][row])])
            self.features_frequency_diff_display.append(tmp_row)

        self.signature = OrderedDict([(x, []) for x in self.aln_pos.segments])
        for segment in self.aln_pos.segments:
            tmp = np.array(self.features_frequency_difference[segment])
            signature_map = np.absolute(tmp).argmax(axis=0)
            #signature_map = tmp.argmax(axis=0)
            # Update mapping to prefer features with fewer amino acids
            signature_map = self._assign_preferred_features(signature_map, segment, self.features_frequency_difference)

            self.signature[segment] = []
            for col, pos in enumerate(list(signature_map)):
                self.signature[segment].append([
                    list(AMINO_ACID_GROUP_PROPERTIES.values())[pos]['display_name_short'] if self.features_frequency_difference[segment][pos][col] > 0 else '-' + list(AMINO_ACID_GROUP_PROPERTIES.values())[pos]['display_name_short'], # latest implementation of NOT... properties
                    list(AMINO_ACID_GROUP_NAMES.values())[pos] if self.features_frequency_difference[segment][pos][col] > 0 else "Not " + list(AMINO_ACID_GROUP_NAMES.values())[pos], # latest implementation of NOT... properties
                    self.features_frequency_difference[segment][pos][col],
                    int(self.features_frequency_difference[segment][pos][col]/20)+5,
                    list(AMINO_ACID_GROUP_PROPERTIES.values())[pos]['length'],
                    list(AMINO_ACID_GROUPS.keys())[pos],
                    list(AMINO_ACID_GROUP_PROPERTIES.values())[pos]['display_name_short']
                ])

        features_pos = OrderedDict()
        self.features_consensus_pos = OrderedDict([(x, []) for x in self.aln_pos.segments])
        for sid, segment in enumerate(self.aln_pos.segments):
            features_pos[segment] = np.array(
                [[x[0] for x in feat[sid]] for feat in self.aln_pos.feature_stats],
                dtype='int'
                )
            features_cons_pos = features_pos[segment].argmax(axis=0)
            features_cons_pos = self._assign_preferred_features(features_cons_pos, segment, features_pos)

            for col, pos in enumerate(list(features_cons_pos)):
                self.features_consensus_pos[segment].append([
                    # list(AMINO_ACID_GROUPS.keys())[pos],
                    list(AMINO_ACID_GROUP_PROPERTIES.values())[pos]['display_name_short'],
                    list(AMINO_ACID_GROUP_NAMES.values())[pos],
                    features_pos[segment][pos][col],
                    int(features_pos[segment][pos][col]/20)+5,
                    list(AMINO_ACID_GROUP_PROPERTIES.values())[pos]['length'],
                    list(AMINO_ACID_GROUPS.keys())[pos]
                ])
        self._convert_feature_stats(self.features_normalized_pos, self.aln_pos)

    def calculate_zscales_signature(self):
        """
        Calculates the Z-scales (Z1-Z5) difference between two protein sets for each GN residue position
        Generates the full difference matrix and calculates the relevance (P-value) for each z-scale & position combination.
        """
        # Prepare zscales for both sets
        self.aln_pos.calculate_zscales()
        self.aln_neg.calculate_zscales()

        # Difference + p-value calculation for shared residues
        ZSCALES.sort()
        for zscale in ZSCALES:
            self.zscales_signature[zscale] = OrderedDict()
            for segment in self.aln_pos.zscales[ZSCALES[0]].keys():
                self.zscales_signature[zscale][segment] = OrderedDict()

                all_keys = set(self.aln_pos.zscales[zscale][segment].keys()).union(self.aln_neg.zscales[zscale][segment].keys())
                shared_keys = set(self.aln_pos.zscales[zscale][segment].keys()).intersection(self.aln_neg.zscales[zscale][segment].keys())
                for entry in sorted(all_keys):
                    if entry in shared_keys:
                        var1 = self.aln_pos.zscales[zscale][segment][entry]
                        var2 = self.aln_neg.zscales[zscale][segment][entry]

                        # Welch's t-test
                        # One-liner alternative : sed = np.sqrt(var1[1]**2.0/var1[2] + var2[1]**2.0/var2[2])
                        # se1 = var1[1]/np.sqrt(var1[2])
                        # se2 = var2[1]/np.sqrt(var2[2])
                        # sed = np.sqrt(se1**2.0 + se2**2.0)

                        # Student t-test assuming similar variance different sample sizes
                        sed = 0
                        if var1[2] > 0 and var2[2] > 0 and (var1[2]+var2[2] - 2) > 0:
                            sed = np.sqrt(((var1[2] - 1) * var1[1]**2.0 + (var2[2]-1)*var2[1]**2.0)/(var1[2]+var2[2]-2)) * np.sqrt(1/var1[2] + 1/var2[2])

                        t_value = 1
                        p = 100
                        color = -1
                        if sed != 0:
                            mean_diff = var1[0] - var2[0]
                            t_value = mean_diff / sed

                            # Grab P-value
                            df = var1[2] + var2[2] - 2
                            p = (1.0 - t.cdf(abs(t_value), df)) * 2.0

                            # Coloring based on statistical significance
                            #if p <= 0.05:
                            #    color = int(round(10 - 9 * p/0.05, 0))
                            #else:
                            #    color = 0

                            # Coloring difference Z-scale means when statistically significant
                            if p <= 0.05 and abs(mean_diff) > 0.6:
                                color = round(mean_diff / 4 * 5, 0)
                                if abs(color) > 5:
                                    color = color/abs(color) * 5
                                color = int(color + 5)

                        tooltip = entry + " ("+ zscale + ")<br>" + \
                                  "Set 1: " + str(round(var1[0], 2)) + " ± " + str(round(var1[1], 2)) + " (" + str(var1[2]) + ")</br>" + \
                                  "Set 2: " + str(round(var2[0], 2)) + " ± " + str(round(var2[1], 2)) + " (" + str(var2[2]) + ")</br>"
                        if p > 0.001:
                            tooltip += "P-value:  {0:.3f}".format(p)
                        else:
                            tooltip += "P-value:  {0:.2E}".format(p)

                        self.zscales_signature[zscale][segment][entry] = [round(var1[0]-var2[0],1), color, tooltip] # diff, P-value, tooltip
                    else:
                        tooltip = entry + "<br/>Set 1: GAP<br/>"
                        if entry in self.aln_pos.zscales[zscale][segment]:
                            var1 = self.aln_pos.zscales[zscale][segment][entry]
                            tooltip = entry + "<br/>Set 1: " + str(round(var1[0], 2)) + " ± " + str(round(var1[1], 2)) + " (" + str(var1[2]) + ")</br>"

                        if entry in self.aln_neg.zscales[zscale][segment]:
                            var2 = self.aln_neg.zscales[zscale][segment][entry]
                            tooltip += "Set 2: " + str(round(var2[0], 2)) + " ± " + str(round(var2[1], 2)) + " (" + str(var2[2]) + ")</br>"
                        else:
                            tooltip += "Set 2: GAP<br/>"

                        self.zscales_signature[zscale][segment][entry] = ["-", -1, tooltip] # diff, P-value, tooltip

    def suggest_mutations(self, ref_prot):
        """
        A function that returns a list of suggestions for point mutations based on sequence signature.
        No cut-off on signature strength is applied here, the function iterates over all residues in the signature.

        Args:
            ref_prot ([Protein]): A protein to suggest the mutations for.
        """
        residues =  { x.generic_number.label: x for x in Residue.objects.filter(protein_conformation=ref_prot) if x.generic_number is not None }
        output = OrderedDict() # A list of lists in the values
        for s, gns in self.common_segments.items():
            for pos, resi in enumerate(gns):
                if resi not in residues:
                    continue
                segment = residues[resi].protein_segment.slug
                signature_prop_group = self.signature[s][pos][5]
                if residues[resi].amino_acid in AMINO_ACID_GROUPS[signature_prop_group]:
                    continue
                signature_strength = self.signature[s][pos][2]
                mutation = ", ".join(AMINO_ACID_GROUPS[signature_prop_group])
                #Adds a list with mutation suggestion. Fields:
                # GN, WT aa, Suggested aa, Signature strength, Set 1/2
                try:
                    output[segment].append(
                        [resi,
                        residues[resi].amino_acid,
                        mutation,
                        abs(signature_strength), "Set 1 " if signature_strength >= 0 else "Set 2"
                        ])
                except KeyError:
                    output[segment] = [ [resi, residues[resi].amino_acid, mutation, abs(signature_strength), "Set 1 " if signature_strength >= 0 else "Set 2"] ]

        return output

    def prepare_display_data_onesided(self):

        options = {
            'common_segments': self.common_segments,
            'common_generic_numbers': self.common_gn,
            'signature_consensus': self.signature,
            'feats_cons_pos': self.features_consensus_pos,
            'a_pos': self.aln_pos,
        }

        return options

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
            'zscales_signature': self.zscales_signature,
            'feats_cons_pos': self.features_consensus_pos,
            'feats_cons_neg': self.features_consensus_neg,
            'a_pos': self.aln_pos,
            'a_neg': self.aln_neg,
            'segment_headers_to_hide': self.segment_headers_to_hide,
        }

        return options

    def prepare_session_data(self):

        session_signature = {
            'common_positions': self.common_gn,
            'diff_matrix': self.features_frequency_difference,
            'numbering_schemes': self.common_schemes,
            'common_segments': self.common_segments,
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
                try:
                    #New way of breaking down the numbering scheme
                    rnsn_parent = prot.protein.residue_numbering_scheme.parent.short_name
                except Exception as msg:
                    rnsn_parent = ''
                numbering_schemes[prot.protein.residue_numbering_scheme.slug] = (rnsn, rnsn_parent)
        # order and convert numbering scheme dict to tuple
        return sorted([(x[0], x[1][0], x[1][1]) for x in numbering_schemes.items()], key=itemgetter(0))

    # def apply_cutoff(self, cutoff=0):

    #     matrix_consensus = OrderedDict()
    #     for segment in self.segments:
    #         # print(segment)
    #         segment_consensus = []
    #         signature_map = np.absolute(self.features_frequency_difference[segment]).argmax(axis=0)
    #         for col, pos in enumerate(list(signature_map)):
    #             if abs(self.features_frequency_difference[segment][pos][col]) > self.cutoff:
    #                 segment_consensus.append(self.features_frequency_difference[segment][ : , col])
    #                 for scheme in self.schemes:
    #                     gnum = list(self.common_gn[scheme[0]][segment].items())[col]
    #                     try:
    #                         self.relevant_gn[scheme[0]][segment][gnum[0]] = gnum[1]
    #                     except:
    #                         self.relevant_gn[scheme[0]][segment] = OrderedDict()
    #                         self.relevant_gn[scheme[0]][segment][gnum[0]] = gnum[1]

    #         segment_consensus = np.array(segment_consensus).T
    #         if segment_consensus != []:
    #             matrix_consensus[segment] = segment_consensus
    #     self.signature_matrix_filtered = matrix_consensus
    #     self.relevant_segments = OrderedDict([
    #         (
    #             x[0],
    #             self.relevant_gn[self.schemes[0][0]][x[0]].keys()
    #         ) for x in self.signature_matrix_filtered.items()
    #     ])
    #     signature = OrderedDict([(x[0], []) for x in matrix_consensus.items()])
    #     for segment in self.relevant_segments:
    #         signature_map = np.absolute(self.signature_matrix_filtered[segment]).argmax(axis=0)
    #         tmp = np.array(self.signature_matrix_filtered[segment])
    #         for col, pos in enumerate(list(signature_map)):
    #             signature[segment].append([
    #                 list(AMINO_ACID_GROUPS.keys())[pos],
    #                 list(AMINO_ACID_GROUP_NAMES.values())[pos],
    #                 tmp[pos][col],
    #                 int(tmp[pos][col]/20)+5
    #             ])
    #     self.signature_consensus = signature

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
            alignment = self.aln_pos
            if data == 'features':
                data_block = self.aln_pos.feature_stats
                feat_consensus = self.features_consensus_pos
        elif aln == 'negative':
            alignment = self.aln_neg
            if data == 'features':
                data_block = self.aln_neg.feature_stats
                feat_consensus = self.features_consensus_neg
        else:
            if data == 'features':
                data_block = self.features_frequency_diff_display
                feat_consensus = self.signature
        numbering_schemes = self.common_schemes
        generic_numbers_set = self.common_gn

        # First column, numbering schemes
        for row, scheme in enumerate(numbering_schemes):
            # worksheet.write(1 + 3*row, 0, scheme[1])
            worksheet.write(1 + 3 * row, 0, 'Residue number')
            worksheet.write(2 + 3 * row, 0, 'Sequence-based ({})'.format(scheme[2]))
            worksheet.write(3 + 3 * row, 0, 'Structure-based (GPCRdb)')

        # First column, stats
        if data == 'features':
            worksheet.write(1 + 3 * len(numbering_schemes) - 1, 1, "Length")
            #for offset, prop in enumerate(props):
            offset = 0
            for key, prop in AMINO_ACID_GROUP_NAMES.items():
                worksheet.write(1 + 3 * len(numbering_schemes) + offset, 0, prop)
                worksheet.write(1 + 3 * len(numbering_schemes) + offset, 1, AMINO_ACID_GROUP_PROPERTIES[key]["length"])
                offset += 1

            if aln == 'signature':
                worksheet.write(
                    1 + 3 * len(numbering_schemes) + len(props),
                    0,
                    'Signature consensus'
                    )
            else:
                worksheet.write(
                    1 + 3 * len(numbering_schemes) + len(props),
                    0,
                    'Prop/AA consensus'
                    )
            worksheet.write(
                2 + 3 * len(numbering_schemes) + len(props),
                0,
                'Length'
            )
        # First column, protein list (for alignment) and line for consensus sequence
        else:
            for offset, prot in enumerate(alignment.proteins):
                worksheet.write(
                    1 + 3 * len(numbering_schemes) + offset,
                    0,
                    prot.protein.entry_name
                )
            worksheet.write(
                1 + 3 * len(numbering_schemes) + len(alignment.proteins),
                0,
                'Seq consensus'
                )

        # Second column and on
        # Segments
        offset = 0
        col_offset = 1
        if data == 'features':
            col_offset = 2
        for segment in generic_numbers_set[numbering_schemes[0][0]].keys():
            worksheet.merge_range(
                0,
                col_offset + offset,
                0,
                len(generic_numbers_set[numbering_schemes[0][0]][segment]) + offset - 1,
                segment
            )
            offset += len(generic_numbers_set[numbering_schemes[0][0]][segment])

        # Generic numbers
        # for row, item in enumerate(generic_numbers_set.items()):
        for row, item in enumerate(numbering_schemes):
            scheme = item[0]
            offset = col_offset
            for _, gn_list in generic_numbers_set[scheme].items():
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
                col_offset = 1
                for segment in prop:
                    for col, freq in enumerate(segment):
                        # if aln == 'signature':
                        #     cell_format = workbook.add_format(get_format_props(freq[1] if freq[0] != 0 else 5))
                        # else:
                        if aln == 'signature':
                            cell_format = workbook.add_format(get_format_props(freq[1]))
                        else:
                            cell_format = workbook.add_format(get_format_props(freq_gs=freq[1]))

                        worksheet.write(
                            offset + row,
                            1 + col + col_offset,
                            freq[0] if isinstance(freq[0], int) else int(freq[0]),
                            cell_format
                        )
                    col_offset += len(segment)
            col_offset = 1
            for segment, cons_feat in feat_consensus.items():
                for col, chunk in enumerate(cons_feat):
                    if aln == 'signature':
                        cell_format = workbook.add_format(get_format_props(feat=chunk[-1].replace('\u03b1', 'a')))
                    else:
                        cell_format = workbook.add_format(get_format_props(feat=chunk[0].replace('\u03b1', 'a')))
                    #Property group
                    worksheet.write(
                        offset + len(AMINO_ACID_GROUPS),
                        1 + col + col_offset,
                        chunk[0],
                        cell_format
                    )
                    #Length of prop
                    worksheet.write(
                        1 + offset + len(AMINO_ACID_GROUPS),
                        1 + col + col_offset,
                        chunk[4],
                    )
                    if aln == 'signature':
                        cell_format = workbook.add_format(get_format_props(int(chunk[2]/20)+5))
                    else:
                        cell_format = workbook.add_format(get_format_props(int(chunk[2]/10))) #if chunk[2] != 0 else get_format_props(-1))
                    #Percentages
                    worksheet.write(
                        2 + offset + len(AMINO_ACID_GROUPS),
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

    def per_gn_signature_excel(self, workbook, worksheet_name='SignByCol'):

        per_gn_signature = []
        for segment in self.common_segments:
            for pos, item in enumerate(self.signature[segment]):
                gn = list(self.common_gn[self.common_schemes[0][0]][segment].keys())[pos]
                if 'x' not in gn:
                    continue # skip positions without a generic number

                prop = AMINO_ACID_GROUP_PROPERTIES[item[5]]['display_name_short']
                length = AMINO_ACID_GROUP_PROPERTIES[item[5]]['length']
                prop_name = item[1]
                score = abs(item[2])

                per_gn_signature.append([gn, score, prop, length, prop_name])

        worksheet = workbook.add_worksheet(worksheet_name)
        worksheet.write_row(0, 0, ['Pos', 'Score', 'Prop/AA', 'Length', 'Property'])
        for row, pos in enumerate(sorted(per_gn_signature, key=lambda x: x[1], reverse=True)):
            worksheet.write_row(
                row + 1,
                0,
                pos,
            )


    def zscales_excel(self, workbook, worksheet_name, aln):

        worksheet = workbook.add_worksheet(worksheet_name)
        if aln == 'signature':
            generic_numbers_set = self.common_gn
            segments = self.common_segments
            zscales = self.zscales_signature
        if aln == 'positive':
            generic_numbers_set = self.aln_pos.generic_numbers
            segments = self.aln_pos.segments
            zscales = self.aln_pos.zscales
        if aln == 'negative':
            generic_numbers_set = self.aln_neg.generic_numbers
            segments = self.aln_neg.segments
            zscales = self.aln_neg.zscales

        numbering_schemes = self.common_schemes

        # First column, numbering schemes
        for row, scheme in enumerate(numbering_schemes):
            worksheet.write(1 + 3 * row, 0, 'Residue number')
            worksheet.write(2 + 3 * row, 0, 'Sequence-based ({})'.format(scheme[2]))
            worksheet.write(3 + 3 * row, 0, 'Structure-based (GPCRdb)')

        # Segments
        offset = 0
        col_offset = 1

        # First ccolumn - labels
        for segment in generic_numbers_set[numbering_schemes[0][0]].keys():
            worksheet.merge_range(
                0,
                col_offset + offset,
                0,
                len(generic_numbers_set[numbering_schemes[0][0]][segment]) + offset - 1,
                segment
            )
            offset += len(generic_numbers_set[numbering_schemes[0][0]][segment])

        # Generic numbers
        for row, item in enumerate(numbering_schemes):
            scheme = item[0]
            offset = col_offset
            for _, gn_list in generic_numbers_set[scheme].items():
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


        # Setting starting point for writing z-scales into excel sheet
        offset = 1 + 3 * len(numbering_schemes)


        # This part is split into alignments and signature since zscale items have different format
        # (4 positions for alignment vs 3 positions for signature zscales)
        # Signature
        if aln == 'signature':

            for z_id, zscale in enumerate(ZSCALES):
                # Z-scale label
                worksheet.write(
                    offset + z_id,
                    0,
                    zscale
                )
                for segment, gns in segments.items():
                    for g_id, gn in enumerate(gns):
                        # Unpack the values form z-scales
                        zs_diff, p_val, tooltip = zscales[zscale][segment][gn]
                        # Colour cell based on p-value
                        cell_format = workbook.add_format(get_format_props(p_val))
                        worksheet.write(
                            offset + z_id, # row
                            col_offset + g_id, # column
                            zs_diff, # Z-scale value
                            cell_format # cell formatting
                        )
                        # Stuff from tooltip goes into comments
                        worksheet.write_comment(
                            offset + z_id, # row
                            col_offset + g_id, # column
                            tooltip.replace('</br>', '\n')
                        )
                    col_offset += len(generic_numbers_set[numbering_schemes[0][0]][segment])
                col_offset = 1
        else:
            for z_id, zscale in enumerate(ZSCALES):
                # Z-scales label
                worksheet.write(
                    offset + z_id,
                    0,
                    zscale
                )
                for segment, gns in segments.items():
                    for g_id, gn in enumerate(gns):
                        res = zscales[zscale][segment][gn]
                        # Unpack the values form z-scales
                        zs_mean, _, _, tooltip = res # mean, std_dev, obs_count, tooltip
                        # Colour cell based on p-value
                        cell_format = workbook.add_format(get_format_props(freq_gs=zscales_color_scale(res)))
                        worksheet.write(
                            offset + z_id, # row
                            col_offset + g_id, # column
                            "%.2f" % zs_mean, # Z-scale value
                            cell_format # cell formatting
                        )
                        # Stuff from tooltip goes into comments
                        worksheet.write_comment(
                            offset + z_id, # row
                            col_offset + g_id, # column
                            tooltip
                        )
                    col_offset += len(generic_numbers_set[numbering_schemes[0][0]][segment])
                col_offset = 1



class SignatureMatch():

    def __init__(self, common_positions, numbering_schemes, segments, difference_matrix,
            protein_set_pos=None, protein_set_neg=None, cutoff=40, signprot=False):

        self.cutoff = cutoff
        self.norm = 0.0
        self.common_gn = common_positions
        self.schemes = numbering_schemes
        self.segments = segments
        self.diff_matrix = difference_matrix
        self.signature_matrix_filtered = OrderedDict()
        self.signature_consensus = OrderedDict()
        if protein_set_pos and protein_set_neg:
            self.protein_set = protein_set_pos + protein_set_neg
        elif protein_set_pos:
            self.protein_set = protein_set_pos
        elif protein_set_neg:
            self.protein_set = protein_set_neg
        self.protein_set_pos = protein_set_pos
        self.protein_set_neg = protein_set_neg
        self.relevant_gn = OrderedDict([(x[0], OrderedDict()) for x in self.schemes])
        self.relevant_segments = OrderedDict()
        self.scored_proteins = []
        self.protein_report = OrderedDict()
        self.protein_signatures = OrderedDict()
        self.feature_preference = prepare_aa_group_preference()

        self.group_lengths = dict([
            (x, len(y)) for x,y in enumerate(AMINO_ACID_GROUPS.values())
        ])

        self.find_relevant_gns()

        self.residue_to_feat = dict(
            [(x, set()) for x in AMINO_ACIDS.keys()]
            )
        for fidx, feat in enumerate(AMINO_ACID_GROUPS.items()):
            for res in feat[1]:
                try:
                    self.residue_to_feat[res].add(fidx)
                except KeyError:
                    self.residue_to_feat['-'].add(fidx)

        self._find_norm()
        if protein_set_pos:
            self.scores_pos, self.signatures_pos, self.scored_proteins_pos = self.score_protein_set(self.protein_set_pos, signprot)
        if protein_set_neg:
            self.scores_neg, self.signatures_neg, self.scored_proteins_neg = self.score_protein_set(self.protein_set_neg, signprot)


    def _assign_preferred_features(self, signature, segment, ref_matrix):

        new_signature = []
        for pos, argmax in enumerate(signature):
            new_signature.append(self._calculate_best_feature(pos, segment, argmax, ref_matrix))
        return new_signature

    def _calculate_best_feature(self, pos, segment, argmax, ref_matrix):

        tmp = self.feature_preference[argmax]
        equiv_feat = np.where(np.isin(ref_matrix[segment][:, pos], ref_matrix[segment][argmax][pos]))[0]
        pref_feat = argmax
        min_len = self.group_lengths[argmax]

        for efeat in equiv_feat:
            if efeat in tmp and self.group_lengths[efeat] < min_len:
                pref_feat = efeat
                min_len = self.group_lengths[efeat]
            # when two features have the same aa count, take the one from positive set
            elif efeat in tmp and self.group_lengths[efeat] == min_len:
                if ref_matrix[segment][pref_feat][pos] < 0 and ref_matrix[segment][efeat][pos] > 0:
                    pref_feat = efeat
        return pref_feat

    def _find_norm(self):

        norm = 0.0
        for segment in self.relevant_segments:
            norm += np.sum(np.amax(np.absolute(self.signature_matrix_filtered[segment]), axis=0))
        self.norm = norm

    def find_relevant_gns(self):
        """
        Find the set of generic residue positions meeting the cutoff.
        """

        matrix_consensus = OrderedDict()
        for segment in self.segments:
            segment_consensus = []
            #signature_map = self.diff_matrix[segment].argmax(axis=0)
            signature_map = np.absolute(self.diff_matrix[segment]).argmax(axis=0)
            # Update mapping to prefer features with fewer amino acids
            signature_map = self._assign_preferred_features(signature_map, segment, self.diff_matrix)
            for col, pos in enumerate(list(signature_map)):
                if abs(self.diff_matrix[segment][pos][col]) >= self.cutoff:
                    segment_consensus.append(self.diff_matrix[segment][ : , col])
                    for scheme in self.schemes:
                        gnum = list(self.common_gn[scheme[0]][segment].items())[col]
                        try:
                            self.relevant_gn[scheme[0]][segment][gnum[0]] = gnum[1]
                        except KeyError:
                            self.relevant_gn[scheme[0]][segment] = OrderedDict()
                            self.relevant_gn[scheme[0]][segment][gnum[0]] = gnum[1]
            segment_consensus = np.array(segment_consensus).T

            if segment_consensus.shape != (0,):
                matrix_consensus[segment] = segment_consensus
        self.signature_matrix_filtered = matrix_consensus
        self.relevant_segments = OrderedDict([
            (
                x[0],
                self.relevant_gn[self.schemes[0][0]][x[0]].keys()
            ) for x in self.signature_matrix_filtered.items()
        ])

        signature = OrderedDict([(x[0], []) for x in matrix_consensus.items()])
        for segment in self.relevant_segments:
            # signature_map = self.signature_matrix_filtered[segment].argmax(axis=0)
            signature_map = np.absolute(self.signature_matrix_filtered[segment]).argmax(axis=0)
            signature_map = self._assign_preferred_features(signature_map, segment, self.signature_matrix_filtered)
            tmp = np.array(self.signature_matrix_filtered[segment])
            for col, pos in enumerate(list(signature_map)):
                signature[segment].append([
                    # list(AMINO_ACID_GROUPS.keys())[pos],
                    list(AMINO_ACID_GROUP_PROPERTIES.values())[pos]['display_name_short'],
                    list(AMINO_ACID_GROUP_NAMES.values())[pos],
                    tmp[pos][col],
                    int(tmp[pos][col]/20)+5,
                    list(AMINO_ACID_GROUP_PROPERTIES.values())[pos]['length'],
                ])
        self.signature_consensus = signature


    def score_protein_class(self, pclass_slug='001', signprot=False):

        start = time.time()
        protein_scores = {}
        protein_signature_match = {}
        class_proteins = Protein.objects.filter(
            species__common_name='Human',
            family__slug__startswith=pclass_slug
            ).exclude(
                id__in=[x.id for x in self.protein_set]
                )

        if signprot:
            complex_objs = SignprotComplex.objects.prefetch_related('structure__protein_conformation__protein').values_list('structure__protein_conformation__protein__parent_id', flat=True)
            class_proteins = class_proteins.exclude(id__in=complex_objs)

        class_a_pcf = ProteinConformation.objects.order_by(
            'protein__family__slug',
            'protein__entry_name'
            ).filter(
                protein__in=class_proteins,
                protein__sequence_type__slug='wt'
            ).exclude(protein__entry_name__endswith='-consensus').prefetch_related('protein','protein__family__parent','protein__species')

        relevant_gns_total = []
        for segment in  self.relevant_segments:
            for idx, pos in enumerate(self.relevant_gn[self.schemes[0][0]][segment].keys()):
                relevant_gns_total.append(pos)

        resi = Residue.objects.filter(
            protein_conformation__in=class_a_pcf,
            generic_number__label__in=relevant_gns_total
            ).prefetch_related('generic_number','protein_conformation')
        resi_dict_all = {}
        for r in resi:
            if r.generic_number:
                pcf = r.protein_conformation.pk
                if pcf not in resi_dict_all:
                    resi_dict_all[pcf] = {}
                resi_dict_all[pcf][r.generic_number.label] = r

        for pcf in class_a_pcf:
            # p_start = time.time()
            score, nscore, signature_match = self.score_protein(pcf, resi_dict_all)
            protein_scores[pcf] = (score, nscore)
            protein_signature_match[pcf] = signature_match
            # p_end = time.time()
            # print("Time elapsed for {}: ".format(pcf.protein.entry_name), p_end - p_start)
        end = time.time()
        self.protein_report = OrderedDict(sorted(protein_scores.items(), key=lambda x: x[1][0], reverse=True))
        for prot in self.protein_report.items():
            self.protein_signatures[prot[0]] = protein_signature_match[prot[0]]
        self.scored_proteins = list(self.protein_report.keys())
        print("Total time: ", end - start)


    def score_protein_set(self, protein_set, signprot=False):

        start = time.time()
        protein_scores = {}
        protein_signature_match = {}

        seq_type_slug=['wt']
        if signprot:
            seq_type_slug.append('mod')

        pcfs = ProteinConformation.objects.order_by(
            'protein__family__slug',
            'protein__entry_name'
            ).filter(
                protein__in=protein_set,
                protein__sequence_type__slug__in=seq_type_slug
                ).exclude(protein__entry_name__endswith='-consensus').prefetch_related('protein')

        relevant_gns_total = []
        for segment in  self.relevant_segments:
            for idx, pos in enumerate(self.relevant_gn[self.schemes[0][0]][segment].keys()):
                relevant_gns_total.append(pos)

        resi = Residue.objects.filter(
            protein_conformation__in=pcfs,
            generic_number__label__in=relevant_gns_total
            ).prefetch_related('generic_number','protein_conformation')
        resi_dict_all = {}
        for r in resi:
            if r.generic_number:
                pcf = r.protein_conformation.pk
                if pcf not in resi_dict_all:
                    resi_dict_all[pcf] = {}
                resi_dict_all[pcf][r.generic_number.label] = r

        for pcf in pcfs:
            # p_start = time.time()
            score, nscore, signature_match = self.score_protein(pcf,resi_dict_all)
            protein_scores[pcf] = (score, nscore)
            protein_signature_match[pcf] = signature_match
            # p_end = time.time()
            # print("Time elapsed for {}: ".format(pcf.protein.entry_name), p_end - p_start)
        end = time.time()
        protein_report = OrderedDict(sorted(protein_scores.items(), key=lambda x: x[1][0], reverse=True))
        protein_signatures = OrderedDict()
        for prot in protein_report.items():
            protein_signatures[prot[0]] = protein_signature_match[prot[0]]
        scored_proteins = list(protein_report.keys())
        print("Total time: ", end - start)

        return (protein_report, protein_signatures, scored_proteins)

    def score_protein(self, pcf,resi_dict_all):
        prot_score = 0.0
        #norm = 0.0
        consensus_match = OrderedDict([(x, []) for x in self.relevant_segments])

        if resi_dict_all == None or pcf.pk not in resi_dict_all:
            relevant_gns_total = []
            for segment in  self.relevant_segments:
                for idx, pos in enumerate(self.relevant_gn[self.schemes[0][0]][segment].keys()):
                    relevant_gns_total.append(pos)
            resi = Residue.objects.filter(
                protein_conformation=pcf,
                generic_number__label__in=relevant_gns_total
                ).prefetch_related('generic_number')
            resi_dict = {}
            for r in resi:
                if r.generic_number:
                    resi_dict[r.generic_number.label] = r
        else:
            resi_dict = resi_dict_all[pcf.pk]
        for segment in self.relevant_segments:
            tmp = []
            # signature_map = self.signature_matrix_filtered[segment].argmax(axis=0)
            signature_map = np.absolute(self.signature_matrix_filtered[segment]).argmax(axis=0)
            signature_map = self._assign_preferred_features(signature_map, segment, self.signature_matrix_filtered)

            #norm += np.sum(np.amax(self.signature_matrix_filtered[segment], axis=0))

            for idx, pos in enumerate(self.relevant_gn[self.schemes[0][0]][segment].keys()):
                feat = signature_map[idx]
                feat_abr = list(AMINO_ACID_GROUPS.keys())[feat]
                feat_name = list(AMINO_ACID_GROUP_NAMES.values())[feat]
                val = self.signature_matrix_filtered[segment][feat][idx]
                if pos in resi_dict:
                    res = resi_dict[pos]
                    if feat in self.residue_to_feat[res.amino_acid]:
                        if val > 0:
                            prot_score += val
                        tmp.append([
                            feat_abr,
                            feat_name,
                            val,
                            #"green",
                            "#808080",
                            res.amino_acid, pos
                            ]) if val > 0 else tmp.append([
                                feat_abr,
                                feat_name,
                                val,
                                "white",
                                res.amino_acid,
                                pos
                                ])
                    else:
                        #David doesn't want the negative values in the score
                        # prot_score -= val
                        if val < 0:
                            prot_score -= val #if a receptor does NOT have the negative property, add the score
                        tmp.append([
                            feat_abr,
                            feat_name,
                            val,
                            # "red",
                            "white",
                            res.amino_acid,
                            pos
                            ]) if val > 0 else tmp.append([
                                feat_abr,
                                feat_name,
                                val,
                                # "green",
                                "#808080",
                                res.amino_acid,
                                pos
                                ])
                else:
                    if feat_name == 'Gap':
                        tmp.append([
                            feat_abr,
                            feat_name,
                            val,
                            # "green",
                            "#808080",
                            '-',
                            pos
                            ]) if val > 0 else tmp.append([
                                feat_abr,
                                feat_name,
                                val,
                                "white",
                                '-',
                                pos
                                ])
                        prot_score += val
                    else:
                        #David doesn't want the negative values in the score
                        #prot_score -= val
                        tmp.append([
                            feat_abr,
                            feat_name,
                            val,
                            # "red",
                            "white",
                            '-',
                            pos
                            ]) if val > 0 else tmp.append([
                                feat_abr,
                                feat_name,
                                val,
                                "white",
                                '-',
                                pos
                                ])
            consensus_match[segment] = tmp
        return (prot_score/100, prot_score/self.norm*100, consensus_match)

def signature_score_excel(workbook, scores, protein_signatures, signature_filtered, relevant_gn, relevant_segments, numbering_schemes, scores_positive=None, scores_negative=None, signatures_positive=None, signatures_negative=None):

    worksheet = workbook.add_worksheet('scored_proteins')
    #wrap = workbook.add_format({'text_wrap': True})


    # First column, numbering schemes
    for row, scheme in enumerate(numbering_schemes):
        worksheet.write(1 + 3 * row, 4, 'Residue number')
        worksheet.write(2 + 3 * row, 4, 'Sequence-based ({})'.format(scheme[2]))
        worksheet.write(3 + 3*row, 4, 'Structure-based (GPCRdb)')
    worksheet.write(2 + 3 * len(numbering_schemes), 0, 'UniProt')
    worksheet.write(2 + 3 * len(numbering_schemes), 1, 'Receptor name (IUPHAR)')
    worksheet.write(2 + 3 * len(numbering_schemes) ,2, 'Receptor family')
    worksheet.write(2 + 3 * len(numbering_schemes) ,3, 'Ligand class')
    # Score header
    worksheet.write(2 + 3 * len(numbering_schemes), 4, 'Score')
    #worksheet.write(1, 3, 'Normalized score')

    offset = 0
    # Segments
    for segment, resi in relevant_segments.items():
        worksheet.merge_range(
            0,
            5 + offset,
            0,
            5 + len(resi) + offset - 1,
            segment
        )
        offset += len(resi)

    # Generic numbers
    # for row, item in enumerate(generic_numbers_set.items()):
    for row, item in enumerate(numbering_schemes):
        scheme = item[0]
        offset = 1
        for _, gn_list in relevant_gn[scheme].items():
            for col, gn_pair in enumerate(gn_list.items()):
                try:
                    tm, bw, gpcrdb = re.split('\.|x', strip_html_tags(gn_pair[1]))
                except:
                    tm, bw, gpcrdb = ('', '', '')
                worksheet.write(
                    1 + 3 * row,
                    4 + col + offset,
                    tm
                )
                worksheet.write(
                    2 + 3 * row,
                    4 + col + offset,
                    bw
                )
                worksheet.write(
                    3 + 3*row,
                    4 + col + offset,
                    gpcrdb
                )
            offset += len(gn_list.items())

    # Line for sequence signature
    worksheet.write(
        1 + 3 * len(numbering_schemes),
        4,
        'CONSENSUS'
        )
    col_offset = 0
    for segment, cons_feat in signature_filtered.items():
        for col, chunk in enumerate(cons_feat):
            worksheet.write(
                1 + 3 * len(numbering_schemes),
                5 + col + col_offset,
                chunk[0]
            )
            cell_format = workbook.add_format(get_format_props(int(chunk[2]/20)+5))
            worksheet.write(
                2 + 3 * len(numbering_schemes),
                5 + col + col_offset,
                chunk[2],
                cell_format
            )
        col_offset += len(cons_feat)



    # Score lines
    row_offset = 0
    for protein, score in scores.items():
        worksheet.write(
            3 + 3 * len(numbering_schemes) + row_offset,
            0,
            protein.protein.entry_name,
        )
        worksheet.write(
            3 + 3 * len(numbering_schemes) + row_offset,
            1,
            protein.protein.name
        )
        worksheet.write(
            3 + 3 * len(numbering_schemes) + row_offset,
            2,
            protein.protein.family.parent.name,
        )
        worksheet.write(
            3 + 3 * len(numbering_schemes) + row_offset,
            3,
            protein.protein.family.parent.parent.name,
        )
        # worksheet.write(
        #     3 + 3 * len(numbering_schemes) + row_offset,
        #     2,
        #     score[0],
        # )
        worksheet.write(
            3 + 3 * len(numbering_schemes) + row_offset,
            4,
            score[1],
        )
        col_offset = 0
        for segment, data in protein_signatures[protein].items():
            for col, res in enumerate(data):
                cell_format = workbook.add_format({'bg_color': res[3],})
                worksheet.write(
                    3 + 3 * len(numbering_schemes) + row_offset,
                    5 + col + col_offset,
                    res[4],
                    cell_format
                )
            col_offset += len(data)
        row_offset += 1


    static_offset = 3 + 3 * len(numbering_schemes) + len(protein_signatures.items())
    #Scores for positive set (if specified)
    if scores_positive:
        worksheet.write(
            static_offset,
            0,
            'Protein set 1'
        )
        static_offset += 1

        row_offset = 0
        for protein, score in scores_positive.items():
            worksheet.write(
                static_offset + row_offset,
                0,
                protein.protein.entry_name,
            )
            worksheet.write(
                static_offset + row_offset,
                1,
                protein.protein.name
            )
            worksheet.write(
                static_offset + row_offset,
                2,
                protein.protein.family.parent.name,
            )
            worksheet.write(
                static_offset + row_offset,
                3,
                protein.protein.family.parent.parent.name,
            )
            # worksheet.write(
            #     3 + 3 * len(numbering_schemes) + row_offset,
            #     2,
            #     score[0],
            # )
            worksheet.write(
                static_offset + row_offset,
                4,
                score[1],
            )
            col_offset = 0
            for segment, data in signatures_positive[protein].items():
                for col, res in enumerate(data):
                    cell_format = workbook.add_format({'bg_color': res[3],})
                    worksheet.write(
                        static_offset + row_offset,
                        5 + col + col_offset,
                        res[4],
                        cell_format
                    )
                col_offset += len(data)
            row_offset += 1
        static_offset += len(scores_positive.items())



    #Scores for negative set (if specified)
    if scores_negative:
        worksheet.write(
            static_offset,
            0,
            'Protein set 2'
        )
        static_offset += 1

        row_offset = 0
        for protein, score in scores_negative.items():
            worksheet.write(
                static_offset + row_offset,
                0,
                protein.protein.entry_name,
            )
            worksheet.write(
                static_offset + row_offset,
                1,
                protein.protein.name
            )
            worksheet.write(
                static_offset + row_offset,
                2,
                protein.protein.family.parent.name,
            )
            worksheet.write(
                static_offset + row_offset,
                3,
                protein.protein.family.parent.parent.name,
            )
            # worksheet.write(
            #     3 + 3 * len(numbering_schemes) + row_offset,
            #     2,
            #     score[0],
            # )
            worksheet.write(
                static_offset + row_offset,
                4,
                score[1],
            )
            col_offset = 0
            for segment, data in signatures_negative[protein].items():
                for col, res in enumerate(data):
                    cell_format = workbook.add_format({'bg_color': res[3],})
                    worksheet.write(
                        static_offset + row_offset,
                        5 + col + col_offset,
                        res[4],
                        cell_format
                    )
                col_offset += len(data)
            row_offset += 1
