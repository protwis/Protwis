import hashlib
import logging
import os
from collections import OrderedDict
from copy import deepcopy
from operator import itemgetter

import numpy as np

from alignment.functions import prepare_aa_group_preference
from Bio.Align import substitution_matrices
from common.definitions import *
from django.conf import settings
from django.core.cache import cache, caches
from django.db.models import Q
from protein.models import (Protein, ProteinConformation, ProteinFamily,
                            ProteinSegment, ProteinState)
from residue.functions import dgn, ggn
from residue.models import (Residue, ResidueGenericNumber, ResidueNumberingScheme)
# from structure.functions import StructureSeqNumOverwrite
from signprot.models import SignprotComplex
from structure.models import Rotamer, Structure

try:
    cache_alignments = caches['alignment_core']
except:
    cache_alignments = cache


class Alignment:
    """A class representing a protein sequence alignment, with or without a reference sequence."""
    def __init__(self):
        self.reference = False
        self.proteins = []
        self.unique_proteins = []
        self.non_matching_proteins = [] # proteins that do not match user specified site definitions
        self.segments = OrderedDict()
        self.segments_only_alignable = []
        self.numbering_schemes = {}
        self.generic_numbers = OrderedDict()
        self.generic_number_objs = {}
        self.positions = set()
        self.consensus = OrderedDict()
        self.forced_consensus = OrderedDict() # consensus sequence where all conflicts are solved by rules
        self.full_consensus = [] # consensus sequence with full residue objects
        self.similarity_matrix = OrderedDict()
        self.amino_acids = []
        self.amino_acid_stats = []
        self.aa_count = OrderedDict()
        self.aa_count_with_protein = OrderedDict()
        self.features = []
        self.features_combo = []
        self.feature_stats = []
        self.feat_consensus = OrderedDict()
        self.default_numbering_scheme = ResidueNumberingScheme.objects.get(slug=settings.DEFAULT_NUMBERING_SCHEME)
        self.states = [settings.DEFAULT_PROTEIN_STATE] # inactive, active etc
        self.use_residue_groups = False
        self.ignore_alternative_residue_numbering_schemes = False # set to true if no numbering is to be displayed
        self.residues_to_delete = []
        self.normalized_scores = OrderedDict()
        self.stats_done = False
        self.zscales = OrderedDict()

        # refers to which ProteinConformation attribute to order by (identity, similarity or similarity score)
        self.order_by = 'similarity'

        # name of custom segment, where individually selected postions are collected
        self.custom_segment_label = 'Custom'
        self.interaction_segment_label = 'I'

        # gap symbols
        self.gaps = ['-', '_']

        # when true, gaps at the beginning or end of a segment have a different symbol than other gaps
        self.show_padding = True

        # prepare the selection order of the property groups for calculation of the feature consensus
        self.feature_preference = prepare_aa_group_preference()
        self.group_lengths = dict([
            (x, len(y)) for x,y in enumerate(AMINO_ACID_GROUPS.values())
        ])

    def __str__(self):
        return str(self.__dict__)

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

    def load_reference_protein(self, protein):
        """Loads a protein into the alignment as a reference."""
        self.reference = True

        # fetch the selected conformations of the protein
        # FIXME take many conformational states into account
        try:
            pconf = ProteinConformation.objects.get(protein=protein)
        except ProteinConformation.DoesNotExist:
            raise Exception ('Protein conformation {} not found for protein {}'.format(self.states[0],
                                                                                       protein.entry_name))

        self.proteins.insert(0, pconf)
        self.update_numbering_schemes()
        self.stats_done = False

    def load_reference_protein_from_selection(self, simple_selection):
        """Read user selection and add selected reference protein."""
        if simple_selection and simple_selection.reference:
            self.load_reference_protein(simple_selection.reference[0].item)
        self.stats_done = False

    def load_proteins(self, proteins):
        """Load a list of protein objects into the alignment."""
        # fetch all conformations of selected proteins
        # FIXME only show inactive?
        protein_conformations = ProteinConformation.objects.order_by('protein__family__slug',
                                                                     'protein__entry_name').filter(protein__in=proteins).select_related('protein__residue_numbering_scheme',
                                                                                                                                        'protein__species', 'state')
        pconfs = OrderedDict()
        for pconf in protein_conformations:
            pconf_label = pconf.__str__()
            if pconf_label not in pconfs:
                pconfs[pconf_label] = {}
            pconfs[pconf_label] = pconf

        for pconf_label, pconf in pconfs.items():
            self.proteins.append(pconf)
        self.update_numbering_schemes()
        self.stats_done = False

    def load_proteins_from_selection(self, simple_selection, only_wildtype = False):
        """Read user selection and add selected proteins."""
        # local protein list
        proteins = []

        # flatten the selection into individual proteins
        for target in simple_selection.targets:
            if target.type == 'protein':
                proteins.append(target.item)
            elif target.type == 'family':
                # species filter
                species_list = []
                for species in simple_selection.species:
                    species_list.append(species.item)

                # annotation filter
                protein_source_list = []
                for protein_source in simple_selection.annotation:
                    protein_source_list.append(protein_source.item)

                family_proteins = Protein.objects.filter(family__slug__startswith=target.item.slug, source__in=(protein_source_list))
                if species_list:
                    family_proteins = family_proteins.filter(species__in=(species_list))

                if only_wildtype:
                    family_proteins = family_proteins.filter(sequence_type__slug="wt")

                family_proteins.select_related('residue_numbering_scheme', 'species')
                for fp in family_proteins:
                    proteins.append(fp)

        # load protein list
        self.load_proteins(proteins)
        self.stats_done = False

    def load_segments(self, selected_segments):
        selected_residue_positions = []
        segment_lookup = {}
        segment_lookup_positions = {}
        unsorted_segments = OrderedDict()

        ## Check segments and find all ResidueGenericNumber to lessen queries afterwards
        normal_segments = []
        segment_positions_lookup = OrderedDict()
        for s in selected_segments:
            if hasattr(s, 'item'):
                selected_segment = s.item
            else:
                selected_segment = s
            # handle residue positions separately
            if selected_segment.__class__.__name__ == 'ResidueGenericNumberEquivalent':
                continue
            normal_segments.append(selected_segment)
            segment_positions_lookup[selected_segment.slug] = []
        if normal_segments:
            # If normal segments, prepare all ResidueGenericNumber for less overhead
            segment_positions = ResidueGenericNumber.objects.filter(protein_segment__in=normal_segments,
                                                                    scheme=self.default_numbering_scheme).select_related('protein_segment').order_by('label')
            for segment_residue in segment_positions:
                segment = segment_residue.protein_segment.slug
                segment_positions_lookup[segment].append(segment_residue)

        for s in selected_segments:
            if hasattr(s, 'item'):
                selected_segment = s.item
                if hasattr(selected_segment, 'only_aligned_residues'):
                    self.segments_only_alignable.append(selected_segment.slug)
            else:
                selected_segment = s

            # handle residue positions separately
            if selected_segment.__class__.__name__ == 'ResidueGenericNumberEquivalent':
                segment_residue = [selected_segment]

                # if in site search, add residue group to selected item
                if hasattr(s, 'properties') and 'site_residue_group' in s.properties:
                    segment_residue.append(s.properties['site_residue_group'])
                selected_residue_positions.append(segment_residue)
                continue

            unsorted_segments[selected_segment.pk] = []
            segment_lookup[selected_segment.pk] = selected_segment.slug
            segment_lookup_positions[selected_segment.pk] = segment_positions_lookup[selected_segment.slug]
            for segment_residue in segment_positions_lookup[selected_segment.slug]:
                unsorted_segments[selected_segment.pk].append(segment_residue.label)

        # Use PK values of segments to sort them before making alignment to ensure logical order
        sorted_segments = sorted(unsorted_segments)
        for s in sorted_segments:
            self.segments[segment_lookup[s]] = unsorted_segments[s]
            # generic numbers in the schemes of all selected proteins
            self.load_generic_numbers(segment_lookup[s], segment_lookup_positions[s])

        # combine individual residue positions into a custom segment
        if selected_residue_positions:
            if self.use_residue_groups:
                interaction_groups = {}
                for residue_position in selected_residue_positions:
                    segment_slug = self.interaction_segment_label + str(residue_position[1])
                    if segment_slug not in self.segments:
                        self.segments[segment_slug] = []
                    if segment_slug not in interaction_groups:
                        interaction_groups[segment_slug] = []
                    self.segments[segment_slug].append(residue_position[0].default_generic_number.label)
                    interaction_groups[segment_slug].append(residue_position)
                for segment_slug, group_residues in interaction_groups.items():
                    self.load_generic_numbers(segment_slug, group_residues)
            else:
                segment_slug = self.custom_segment_label
                if segment_slug not in self.segments:
                    self.segments[segment_slug] = []
                for residue_position in selected_residue_positions:
                    self.segments[segment_slug].append(residue_position[0].default_generic_number.label)
                self.load_generic_numbers(segment_slug, selected_residue_positions)
        self.stats_done = False

    def load_segments_from_selection(self, simple_selection):
        """Read user selection and add selected protein segments/residue positions."""
        # local segment list
        segments = []

        # read selection
        for segment in simple_selection.segments:
            segments.append(segment)

        # load segment positions
        self.load_segments(segments)
        self.stats_done = False

    def load_generic_numbers(self, segment_slug, residues):
        """Loads generic numbers in the schemes of all selected proteins."""
        for ns in self.numbering_schemes:
            if ns[0] not in self.generic_numbers:
                self.generic_numbers[ns[0]] = OrderedDict()
            self.generic_numbers[ns[0]][segment_slug] = OrderedDict()
            for segment_residue in residues:
                if segment_slug == self.custom_segment_label or self.use_residue_groups:
                    residue_position = segment_residue[0].default_generic_number.label
                else:
                    residue_position = segment_residue.label
                self.generic_numbers[ns[0]][segment_slug][residue_position] = []
        self.stats_done = False

    def update_numbering_schemes(self):
        """Update numbering scheme list."""
        self.numbering_schemes = {}
        for pc in self.proteins:
            if pc.protein.residue_numbering_scheme.slug not in self.numbering_schemes:
                rnsn = pc.protein.residue_numbering_scheme.name
                try:
                    #New way of breaking down the numbering scheme
                    rnsn_parent = pc.protein.residue_numbering_scheme.parent.short_name
                except:
                    if pc.protein.residue_numbering_scheme.slug in ['cgn','can']:
                        rnsn_parent = pc.protein.residue_numbering_scheme.short_name
                    else:
                        rnsn_parent = ''
                self.numbering_schemes[pc.protein.residue_numbering_scheme.slug] = (rnsn, rnsn_parent)

        # order and convert numbering scheme dict to tuple
        self.numbering_schemes = sorted([(x[0], x[1][0], x[1][1]) for x in self.numbering_schemes.items()], key=itemgetter(0))
        self.stats_done = False

    def get_hash(self):
        # create unique hash key for alignment combo
        protein_ids = sorted(set([ str(protein.id) for protein in self.proteins ]))
        segment_ids = sorted(set( self.segments ))
        hash_key = "|" + "-".join(protein_ids)
        hash_key += "|" + "-".join(segment_ids)
        if 'Custom' in self.segments:
            hash_key += "|" + "-".join(sorted(set( self.segments['Custom'] )))
        hash_key += "|" + "-".join(sorted(set([ scheme[0] for scheme in self.numbering_schemes])))
        hash_key += "|" + "-".join(sorted(set(self.segments_only_alignable)))
        hash_key += "|" + str(self.ignore_alternative_residue_numbering_schemes)
        hash_key += "|" + str(self.custom_segment_label)
        hash_key += "|" + str(self.use_residue_groups)

        return hashlib.md5(hash_key.encode('utf-8')).hexdigest()

    # AJK: point for optimization - primary bottleneck (#1 cleaning, #2 last for-loop in this function)
    def build_alignment(self):
        """Fetch selected residues from DB and build an alignment."""
        # AJK: prevent prefetching all data for large alignments before checking #residues (DB + memory killer)
        rs = Residue.objects.filter(protein_segment__slug__in=self.segments, protein_conformation__in=self.proteins)


        self.number_of_residues_total = len(rs)
        if self.number_of_residues_total>125000: #300 receptors, 400 residues limit
            return "Too large"

        # AJK: performance boost -> Internal caching (not for very small alignments)
        # AJK: note -> ideally we would have fully cached alignments and only select the relevant segments afterwards.
        cache_key = "ALIGNMENTS_"+self.get_hash()

        # cache_alignments.set(cache_key, 0, 0)
        if self.number_of_residues_total < 2500 or not cache_alignments.has_key(cache_key):
            # fetch segment residues
            if not self.ignore_alternative_residue_numbering_schemes and len(self.numbering_schemes) > 1:
                rs = Residue.objects.filter(
                    protein_segment__slug__in=self.segments, protein_conformation__in=self.proteins).prefetch_related(
                    'protein_conformation__protein', 'protein_conformation__state', 'protein_segment',
                    'generic_number__scheme', 'display_generic_number__scheme', 'alternative_generic_numbers__scheme')
            else:
                rs = Residue.objects.filter(
                    protein_segment__slug__in=self.segments, protein_conformation__in=self.proteins).prefetch_related(
                    'protein_conformation__protein', 'protein_conformation__state', 'protein_segment',
                    'generic_number__scheme', 'display_generic_number__scheme')

            # If segment flagged to only include the alignable residues, exclude the ones with no GN
            for s in self.segments_only_alignable:
                rs = rs.exclude(protein_segment__slug=s, generic_number=None)

            # fetch individually selected residues (Custom segment)
            crs = {}
            for segment in self.segments:
                if segment == self.custom_segment_label or self.use_residue_groups:
                    if not self.ignore_alternative_residue_numbering_schemes and len(self.numbering_schemes) > 1:
                        crs[segment] = Residue.objects.filter(
                            generic_number__label__in=self.segments[segment],
                            protein_conformation__in=self.proteins).prefetch_related(
                            'protein_conformation__protein', 'protein_conformation__state', 'protein_segment',
                            'generic_number__scheme', 'display_generic_number__scheme', 'alternative_generic_numbers__scheme')
                    else:
                        crs[segment] = Residue.objects.filter(
                            generic_number__label__in=self.segments[segment],
                            protein_conformation__in=self.proteins).prefetch_related(
                            'protein_conformation__protein', 'protein_conformation__state', 'protein_segment',
                            'generic_number__scheme', 'display_generic_number__scheme')

            # create a dict of proteins, segments and residues
            proteins = {}
            segment_counters = {}
            aligned_residue_encountered = {}
            fusion_protein_inserted = {}
            for r in rs:
                ps = r.protein_segment.slug

                # identifiers for protein/state
                pcid = r.protein_conformation.protein.entry_name + "-" + r.protein_conformation.state.slug

                # update protein dict
                if pcid not in proteins:
                    proteins[pcid] = {}
                if ps not in proteins[pcid]:
                    proteins[pcid][ps] = {}

                # update aligned residue tracker
                if pcid not in aligned_residue_encountered:
                    aligned_residue_encountered[pcid] = {}
                if ps not in aligned_residue_encountered[pcid]:
                    aligned_residue_encountered[pcid][ps] = False

                # what part of the segment is this? There are 4 possibilities:
                # 1. The aligned part (for both fully and partially aligned segments)
                # 2. The part before the aligned part in a partially aligned segment
                # 3. The part after the aligned part in a partially aligned segment
                # 4. An unaligned segment (then there is only one part)
                if r.generic_number:
                    segment_part = 1
                elif ps in settings.REFERENCE_POSITIONS and not aligned_residue_encountered[pcid][ps]:
                    segment_part = 2
                elif ps in settings.REFERENCE_POSITIONS and aligned_residue_encountered[pcid][ps]:
                    segment_part = 3
                else:
                    segment_part = 4

                # update segment counters
                if pcid not in segment_counters:
                    segment_counters[pcid] = {}
                if segment_part == 3:
                    part_ps = ps + '_after'
                else:
                    part_ps = ps
                if part_ps not in segment_counters[pcid]:
                    segment_counters[pcid][part_ps] = 1
                else:
                    segment_counters[pcid][part_ps] += 1


                # update fusion protein tracker
                if pcid not in fusion_protein_inserted:
                    fusion_protein_inserted[pcid] = {}
                if ps not in fusion_protein_inserted[pcid]:
                    fusion_protein_inserted[pcid][ps] = False

                # user generic numbers as keys for aligned segments
                if r.generic_number:
                    proteins[pcid][ps][r.generic_number.label] = r

                    # register the presence of an aligned residue
                    aligned_residue_encountered[pcid][ps] = True
                # use custom keys for non-aligned segments
                else:
                    # label prefix + index
                    # Unaligned segments should be split in the middle, with the first part "left aligned", and the second
                    # "right aligned". If there is an aligned part of the segment, it goes in the middle.
                    if segment_part == 2:
                        prefix = '00-'
                    elif segment_part == 3:
                        prefix = 'zz-'
                    else:
                        prefix = '01-'

                    # Note that there is not enough information to assign correct indicies to "right aligned" residues, but
                    # those are corrected below
                    index = str("%04d" % (segment_counters[pcid][part_ps],))

                    # position label
                    pos_label =  prefix + ps + "-" + index

                    # residue
                    proteins[pcid][ps][pos_label] = r

            # correct alignment of split segments
            for pcid, segments in proteins.items():
                for ps, positions in segments.items():
                    pos_num = 1
                    pos_num_after = 1
                    for pos_label in sorted(positions):
                        res_obj = proteins[pcid][ps][pos_label]
                        right_align = False
                        # In a "normal", non split, unaligned segment, is this past the middle?
                        if (pos_label.startswith('01-')
                                and res_obj.protein_segment.category != 'terminus'
                                and pos_num > (segment_counters[pcid][ps] / 2 + 0.5)):
                            right_align = True
                        # In an partially aligned segment (prefixed with 00), where conserved residues are lacking, treat
                        # as an unaligned segment
                        elif (pos_label.startswith('00-')
                              and not aligned_residue_encountered[pcid][ps]
                              and pos_num > (segment_counters[pcid][ps] / 2 + 0.5)
                              or res_obj.protein_segment.slug == 'N-term'):
                            right_align = True
                        # In an N-terminus, always right align everything
                        elif pos_label.startswith('01-') and res_obj.protein_segment.slug == 'N-term':
                            right_align = True

                        if right_align:
                            # if so, "right align" from here using a zz prefixed label
                            updated_index = 'zz' + pos_label[2:]
                            proteins[pcid][ps][updated_index] = proteins[pcid][ps].pop(pos_label)
                            pos_label = updated_index

                        if pos_label.startswith('zz-'):
                            segment_label_after = ps + '_after' # parts after a partly aligned segment start with zz
                            if segment_label_after in segment_counters[pcid]:
                                segment_length = segment_counters[pcid][segment_label_after]
                                counter = pos_num_after

                            # this might be the "second part" of an unaligned segment, e.g.
                            # AAAA----AAAAA
                            # AAAAAAAAAAAAA
                            else:
                                segment_length = segment_counters[pcid][ps]
                                counter = pos_num

                            updated_index = pos_label[:-4] + str(9999 - (segment_length - counter))
                            proteins[pcid][ps][updated_index] = proteins[pcid][ps].pop(pos_label)
                            pos_label = updated_index
                            pos_num_after += 1
                        if pos_label not in self.segments[ps]:
                            self.segments[ps].append(pos_label)
                        pos_num += 1

            # individually selected residues (Custom segment)
            for segment in self.segments:
                if segment == self.custom_segment_label or self.use_residue_groups:
                    for r in crs[segment]:
                        ps = segment
                        pcid = r.protein_conformation.protein.entry_name + "-" + r.protein_conformation.state.slug
                        if pcid not in proteins:
                            proteins[pcid] = {}
                        if ps not in proteins[pcid]:
                            proteins[pcid][ps] = {}
                        proteins[pcid][ps][r.generic_number.label] = r

            # remove split segments from segment list and order segment positions
            for segment, positions in self.segments.items():
                s = segment.split("_")
                if len(s) > 1:
                    del self.segments[segment]
                else:
                    sorted_segment = []
                    if segment=='B.GPS':
                        if 'B.GPS-2' in self.segments[segment]:
                            sorted_segment.append('B.GPS-2')
                        if 'B.GPS-1' in self.segments[segment]:
                            sorted_segment.append('B.GPS-1')
                        for gn in sorted(self.segments[segment]):
                            if not gn.endswith('-2') and not gn.endswith('-1') and not gn.endswith('+1'):
                                sorted_segment.append(gn)
                        if 'B.GPS+1' in self.segments[segment]:
                            sorted_segment.append('B.GPS+1')
                    else:
                        for gn in sorted(self.segments[segment], key=lambda x: x.split('x')):
                            sorted_segment.append(gn)
                    self.segments[segment] = sorted_segment

            self.unique_proteins = list(set(self.proteins))
            for pc in self.unique_proteins:
                row = OrderedDict()
                pc.alignment_list = []
                for segment, positions in self.segments.items():
                    s = []
                    first_residue_found = False

                    # counters to keep track of gaps at the end of a segment
                    gap_counter = 0
                    position_counter = 1

                    # numbering scheme
                    ns_slug = pc.protein.residue_numbering_scheme.slug

                    # loop all positions in this segment
                    for pos in positions:
                        try:
                            # find the residue record from the dict defined above
                            pcid = pc.protein.entry_name + "-" + pc.state.slug
                            r = proteins[pcid][segment][pos]

                            # add position to the list of positions that are not empty
                            if pos not in self.positions:
                                self.positions.add(pos)

                            # add display number to list of display numbers for this position
                            if r.display_generic_number:
                                if pos not in self.generic_numbers[ns_slug][segment]:
                                    self.generic_numbers[ns_slug][segment][pos] = []
                                if r.display_generic_number.label not in self.generic_numbers[ns_slug][segment][pos]:
                                    self.generic_numbers[ns_slug][segment][pos].append(r.display_generic_number.label)
                            else:
                                if pos not in self.generic_numbers[ns_slug][segment]:
                                    self.generic_numbers[ns_slug][segment][pos] = []

                            # add display numbers for other numbering schemes of selected proteins
                            if (not self.ignore_alternative_residue_numbering_schemes and len(self.numbering_schemes) > 1):
                                if r.generic_number:
                                    for arn in r.alternative_generic_numbers.all():
                                        for ns in self.numbering_schemes:
                                            if (arn.scheme.slug == ns[0] and arn.scheme.slug != ns_slug):
                                                self.generic_numbers[arn.scheme.slug][segment][pos].append(arn.label)
                                                break
                                else:
                                    for ns in self.numbering_schemes:
                                        if pos not in self.generic_numbers[ns[0]][segment] and ns[0] != ns_slug:
                                            self.generic_numbers[ns[0]][segment][pos] = []

                            # append the residue to the matrix
                            if r.generic_number:
                                # s.append([pos, r.display_generic_number.label, r.amino_acid,
                                #   r.display_generic_number.scheme.short_name, r.sequence_number])

                                s.append([pos, r.display_generic_number.label, r.amino_acid,
                                          r.display_generic_number.scheme.short_name, r.sequence_number, r.generic_number.label])


                                # update generic residue object dict
                                if pos not in self.generic_number_objs:
                                    self.generic_number_objs[pos] = r.display_generic_number
                            else:
                                s.append([pos, "", r.amino_acid, "", r.sequence_number])

                            first_residue_found = True

                            # reset gap counter
                            gap_counter = 0
                        except:
                            if self.show_padding:
                                padding_symbol = '_'
                            else:
                                padding_symbol = '-'
                            if first_residue_found:
                                s.append([pos, False, '-', 0])

                                # update gap counter
                                gap_counter += 1

                                # if this is the last residue and there are gaps and the end of the segment, update them to
                                # end gaps
                                if self.show_padding:
                                    if (position_counter) == len(positions):
                                        for i in range(gap_counter):
                                            s[len(positions)-(i+1)][2] = padding_symbol
                            else:
                                s.append([pos, False, padding_symbol, 0])

                        # update position counter
                        position_counter += 1
                    row[segment] = s
                    pc.alignment_list.append(s)
                pc.alignment = row
            #                pc.alignment_list = row_list # FIXME redundant, remove when dependecies are removed

            self.sort_generic_numbers()
            self.merge_generic_numbers()
            self.clear_empty_positions()
            # self.clear_empty_segments()

            if self.number_of_residues_total >= 2500:
                self.calculate_statistics()
                cache_data = {'unique_proteins': self.unique_proteins,
                              'amino_acids': self.amino_acids,
                              'amino_acid_stats': self.amino_acid_stats,
                              'aa_count': self.aa_count,
                              'aa_count_with_protein': self.aa_count_with_protein,
                              'gaps': self.gaps,
                              'feat_consensus': self.feat_consensus,
                              'features': self.features,
                              'features_combo': self.features_combo,
                              'feature_stats': self.feature_stats,
                              'consensus': self.consensus,
                              'forced_consensus': self.forced_consensus,
                              'full_consensus': self.full_consensus,
                              'generic_number_objs': self.generic_number_objs,
                              'generic_numbers': self.generic_numbers,
                              #                            'numbering_schemes': self.numbering_schemes,
                              'positions': self.positions,
                              'segments': self.segments,
                              'zscales': self.zscales}
                cache_alignments.set(cache_key, cache_data, 60*60*24*14)
        else:
            cache_data = cache_alignments.get(cache_key)
            self.unique_proteins = cache_data['unique_proteins']
            self.amino_acids = cache_data['amino_acids']
            self.amino_acid_stats = cache_data['amino_acid_stats']
            self.aa_count = cache_data['aa_count']
            self.aa_count_with_protein = cache_data['aa_count_with_protein']
            self.gaps = cache_data['gaps']
            self.feat_consensus = cache_data['feat_consensus']
            self.features = cache_data['features']
            self.features_combo = cache_data['features_combo']
            self.feature_stats = cache_data['feature_stats']
            self.consensus = cache_data['consensus']
            self.forced_consensus = cache_data['forced_consensus']
            self.full_consensus = cache_data['full_consensus']
            self.generic_number_objs = cache_data['generic_number_objs']
            self.generic_numbers = cache_data['generic_numbers']
            #            self.numbering_schemes = cache_data['numbering_schemes']
            self.positions = cache_data['positions']
            self.segments = cache_data['segments']
            self.zscales = cache_data['zscales']
            self.stats_done = True

        # Adapt alignment to order in current self.proteins
        self.proteins = [prot2 for prot1 in self.proteins for prot2 in self.unique_proteins if prot1.id==prot2.id]

    def remove_non_generic_numbers_from_alignment(self):
        """Remove all positions without a generic number from the protein alignment property."""
        to_delete = {}
        for prot in self.proteins:
            for seg in prot.alignment:
                if seg not in to_delete:
                    to_delete[seg] = []
                for i, res in enumerate(self.proteins[0].alignment[seg]):
                    if 'x' not in res[0] and i not in to_delete[seg]:
                        to_delete[seg].append(i)
        for prot in self.proteins:
            for s, ids in to_delete.items():
                offset = 0
                for i in ids:
                    prot.alignment[s].pop(i-offset)
                    offset+=1

    def clear_empty_positions(self):
        """Remove empty columns from the segments and matrix."""
        # segments

        # AJK optimized cleaning alignment - deepcopy not required and faster removal
        #generic_numbers = deepcopy(self.generic_numbers) # deepcopy is required because the dictionary changes during the loop
        for ns, segments in self.generic_numbers.items():
            for segment, positions in segments.items():
                self.generic_numbers[ns][segment] = OrderedDict([(pos, self.generic_numbers[ns][segment][pos]) for pos in self.generic_numbers[ns][segment].keys() if pos in self.positions ])
                self.segments[segment] = [ pos for pos in self.segments[segment] if pos in self.positions ]

                # Deprecated
        #                for pos, dn in positions.items():
        #                    if pos not in self.positions:
        # remove position from generic numbers dict
        #                        del self.generic_numbers[ns][segment][pos]

        # remove position from segment dict
        #                        if pos in self.segments[segment]:
        #                            self.segments[segment].remove(pos)

        # proteins
        # AJK optimized cleaning alignment - deepcopy not required and faster removal
        #proteins = deepcopy(self.proteins) # deepcopy is required because the list changes during the loop
        for i, protein in enumerate(self.unique_proteins):
            for j, s in protein.alignment.items():
                self.unique_proteins[i].alignment[j] = [ p for p in s if p[0] in self.positions ]

                # Deprecated
    #                for p in s:
    #                    if p[0] not in self.positions:
    #                        self.proteins[i].alignment[j].remove(p)

    def clear_empty_segments(self):
        # clear empty segments
        clear_list = []
        for segment, resis in self.segments.items():
            if len(resis)==0:
                clear_list.append(segment)
        for segment in clear_list:
            del self.segments[segment]
            for ns, segments in self.generic_numbers.items():
                del self.generic_numbers[ns][segment]
            # have to remove segments from protein.alignment for statistics
            for p in self.unique_proteins:
                del p.alignment[segment]

    def merge_generic_numbers(self):
        """Check whether there are many display numbers for each position, and merge them if there are."""
        # deepcopy is required because the dictionary changes during the loop
        # generic_numbers = deepcopy(self.generic_numbers)
        for ns, segments in self.generic_numbers.items():
            for segment, positions in segments.items():
                for pos, dns in positions.items():
                    if not dns: # don't format if there are no numbers
                        self.generic_numbers[ns][segment][pos] = ""
                    elif len(dns) == 1:
                        self.generic_numbers[ns][segment][pos] = self.format_generic_number(dns[0])
                    else:
                        self.generic_numbers[ns][segment][pos] = '-'.join(dns)

    def sort_generic_numbers(self):
        # remove split segments from generic numbers list
        for ns, gs in self.generic_numbers.items():
            for segment, positions in gs.items():
                s = segment.split("_")
                if len(s) > 1:
                    del self.generic_numbers[ns][segment]
                else:
                    ordered_generic_numbers = OrderedDict()

                    if segment=='B.GPS':
                        if 'B.GPS-2' in self.generic_numbers[ns][segment]:
                            ordered_generic_numbers['B.GPS-2'] = self.generic_numbers[ns][segment]['B.GPS-2']
                        if 'B.GPS-1' in self.generic_numbers[ns][segment]:
                            ordered_generic_numbers['B.GPS-1'] = self.generic_numbers[ns][segment]['B.GPS-1']
                        for gn in sorted(self.generic_numbers[ns][segment]):
                            if not gn.endswith('-2') and not gn.endswith('-1') and not gn.endswith('+1'):
                                ordered_generic_numbers[gn] = self.generic_numbers[ns][segment][gn]

                        if 'B.GPS+1' in self.generic_numbers[ns][segment]:
                            ordered_generic_numbers['B.GPS+1'] = self.generic_numbers[ns][segment]['B.GPS+1']
                    else:
                        for gn in sorted(self.generic_numbers[ns][segment], key=lambda x: x.split('x')):
                            ordered_generic_numbers[gn] = self.generic_numbers[ns][segment][gn]
                    self.generic_numbers[ns][segment] = ordered_generic_numbers

    def format_generic_number(self, generic_number):
        """A placeholder for an instance specific function."""
        return generic_number

    def calculate_statistics(self, ignore={}):
        """Calculate consensus sequence and amino acid and feature frequency."""
        if not self.stats_done:
            feature_count = OrderedDict()
            most_freq_aa = OrderedDict()
            amino_acids = OrderedDict([(a, 0) for a in AMINO_ACIDS]) # from common.definitions
            self.amino_acids = list(AMINO_ACIDS.keys())

            # AJK: optimized for speed - split into multiple steps
            for i, p in enumerate(self.unique_proteins):
                entry_name = p.protein.entry_name
                # print(entry_name)
                for j, s in p.alignment.items():
                    if i == 0:
                        self.aa_count[j] = OrderedDict()
                    for p in s:
                        generic_number = p[0]
                        display_generic_number = p[1]
                        amino_acid = p[2]

                        ignore_list = ignore.get(generic_number, [])

                        # Indicate gap and collect statistics
                        if amino_acid in self.gaps:
                            amino_acid = '-'

                        # Skip when unknown amino acid type
                        elif amino_acid == 'X':
                            continue

                        # skip when position is on the ignore list
                        if entry_name in ignore_list:
                            continue

                        if ignore and amino_acid in self.gaps:
                            continue

                        # Init counters
                        if generic_number not in self.aa_count[j]:
                            self.aa_count[j][generic_number] = amino_acids.copy()
                            #                        if generic_number in self.generic_number_objs:
                            self.aa_count_with_protein[generic_number] = {}

                        # update amino acid counter for this generic number
                        self.aa_count[j][generic_number][amino_acid] += 1

                        #if generic_number in self.generic_number_objs:
                        if amino_acid not in self.aa_count_with_protein[generic_number]:
                            self.aa_count_with_protein[generic_number][amino_acid] = {entry_name}
                        elif entry_name not in self.aa_count_with_protein[generic_number][amino_acid]:
                            self.aa_count_with_protein[generic_number][amino_acid].add(entry_name)


            # AJK: Only need for update at the end of loop
            # 1. Update most frequent amino_acids for this generic number
            # 2. Update feature counter for this generic number
            features = OrderedDict([(a, 0) for a in AMINO_ACID_GROUPS])
            self.features_combo = [(x, y['display_name_short'], y['length']) for x,y in zip(list(AMINO_ACID_GROUP_NAMES.values()), list(AMINO_ACID_GROUP_PROPERTIES.values()))]
            self.features = list(AMINO_ACID_GROUP_NAMES.values())

            for j in self.aa_count:
                most_freq_aa[j] = OrderedDict()
                feature_count[j] = OrderedDict()
                for generic_number in self.aa_count[j]:
                    feature_count[j][generic_number] = features.copy()
                    most_freq_aa[j][generic_number] = [[], 0]
                    for amino_acid in self.aa_count[j][generic_number]:
                        # Most frequent AA
                        if self.aa_count[j][generic_number][amino_acid] > most_freq_aa[j][generic_number][1]:
                            most_freq_aa[j][generic_number] = [[amino_acid], self.aa_count[j][generic_number][amino_acid]]
                        elif self.aa_count[j][generic_number][amino_acid] == most_freq_aa[j][generic_number][1]:
                            if amino_acid not in most_freq_aa[j][generic_number][0]:
                                most_freq_aa[j][generic_number][0].append(amino_acid)

                        # create property frequency
                        if self.aa_count[j][generic_number][amino_acid] > 0:
                            for feature in AMINO_ACID_GROUPS_AA[amino_acid]:
                                feature_count[j][generic_number][feature] += self.aa_count[j][generic_number][amino_acid]

            # merge the amino acid counts into a consensus sequence
            num_proteins = len(self.unique_proteins)
            sequence_counter = 1
            for i, s in most_freq_aa.items():
                self.consensus[i] = OrderedDict()
                self.forced_consensus[i] = OrderedDict()
                if i == 'Custom' and s == 'x':
                    sorted_res = sorted(s, key=lambda x: (x.split("x")[0], x.split("x")[1]))
                elif i == 'Custom' and s == '.':
                    sorted_res = sorted(s, key=lambda x: (x.split(".")[0], x.split(".")[1]))
                else:
                    sorted_res = sorted(s)
                for p in sorted_res:
                    r = s[p]
                    conservation = str(round(r[1]/num_proteins*100))
                    if len(conservation) == 1:
                        cons_interval = '0'
                    else:
                        # the intervals are defined as 0-10, where 0 is 0-9, 1 is 10-19 etc. Used for colors.
                        cons_interval = conservation[:-1]

                    # forced consensus sequence uses the first residue to break ties
                    self.forced_consensus[i][p] = r[0][0]

                    # consensus sequence displays + in tie situations
                    num_freq_aa = len(r[0])
                    if num_freq_aa == 1:
                        # Use raw data
                        self.consensus[i][p] = [
                            r[0][0],
                            cons_interval,
                            round(r[1]/num_proteins*100),
                            ""
                        ]
                    elif num_freq_aa > 1 and ignore:
                        self.consensus[i][p] = [
                            r[0][0],
                            cons_interval,
                            round(r[1]/num_proteins*100),
                            ", ".join(r[0])
                        ]
                    elif num_freq_aa > 1:
                        self.consensus[i][p] = [
                            '+',
                            cons_interval,
                            round(r[1]/num_proteins*100),
                            ", ".join(r[0])
                        ]

                    # create a residue object full consensus
                    res = Residue()
                    res.sequence_number = sequence_counter
                    if p in self.generic_number_objs:
                        res.display_generic_number = self.generic_number_objs[p]
                    res.family_generic_number = p
                    res.segment_slug = i
                    res.amino_acid = r[0][0]
                    res.frequency = self.consensus[i][p][2]
                    self.full_consensus.append(res)

                    # update sequence counter
                    sequence_counter += 1

            # process amino acid frequency
            for i, amino_acid in enumerate(AMINO_ACIDS):
                self.amino_acid_stats.append([])
                j = 0
                for segment, segment_num in self.aa_count.items():
                    self.amino_acid_stats[i].append([])
                    k = 0
                    # if segment == 'Custom' and segment_num == 'x':
                    #     sorted_res = sorted(segment_num, key=lambda x: (x.split("x")[0], x.split("x")[1]))
                    # elif segment == 'Custom' and segment_num == '.':
                    #     sorted_res = sorted(segment_num, key=lambda x: (x.split(".")[0], x.split(".")[1]))
                    # else:
                    #sorted_res = segment_num
                    for gn in self.generic_numbers[self.numbering_schemes[0][0]][segment]:
                        aas = segment_num[gn]
                        self.amino_acid_stats[i][j].append([])
                        for aa, freq in aas.items():
                            if aa == amino_acid:
                                frequency = str(round(freq/num_proteins*100))
                                if len(frequency) == 1:
                                    freq_interval = '0'
                                else:
                                    # intervals defined in the same way as for the consensus sequence
                                    freq_interval = frequency[:-1]
                                self.amino_acid_stats[i][j][k] = [frequency, freq_interval]
                        k += 1
                    j += 1

            # create index and prepare stats array
            index_AAG = {}
            for i, feature in enumerate(AMINO_ACID_GROUPS):
                index_AAG[feature] = i
                self.feature_stats.append([])
                for segment in feature_count:
                    self.feature_stats[i].append([])

            j = 0
            for segment, segment_num in feature_count.items():
                k = 0

                # first_gn = list(segment_num.items())[0]
                # if segment == 'Custom' and 'x' in first_gn:
                #     sorted_res = sorted(segment_num, key=lambda x: (x.split("x")[0], x.split("x")[1]))
                # elif segment == 'Custom' and '.' in first_gn:
                #     sorted_res = sorted(segment_num, key=lambda x: (x.split(".")[0], x.split(".")[1]))
                # else:

                for gn in self.generic_numbers[self.numbering_schemes[0][0]][segment]:
                    fs = segment_num[gn]

                    for f, freq in fs.items():
                        # find feature key
                        i = index_AAG[f]

                        frequency = str(round(freq/num_proteins*100))
                        if len(frequency) == 1:
                            freq_interval = '0'
                        else:
                            # intervals defined in the same way as for the consensus sequence
                            freq_interval = frequency[:-1]

                        self.feature_stats[i][j].append([])
                        self.feature_stats[i][j][k] = [frequency, freq_interval]
                    k += 1
                j += 1

            # process feature frequency
            feats = OrderedDict()
            self.feat_consensus = OrderedDict([(x, []) for x in self.segments])
            for sid, segment in enumerate(self.segments):

                # Feature_stats is accessed via feature, segment, residue
                feats[segment] = np.array(
                    [[x[0] for x in feat[sid]] for feat in self.feature_stats],
                    dtype='int'
                )

                feat_cons_tmp = feats[segment].argmax(axis=0)
                feat_cons_tmp = self._assign_preferred_features(feat_cons_tmp, segment, feats)
                for col, pos in enumerate(list(feat_cons_tmp)):
                    self.feat_consensus[segment].append([
                        list(AMINO_ACID_GROUP_PROPERTIES.values())[pos]['display_name_short'],
                        list(AMINO_ACID_GROUP_NAMES.values())[pos],
                        feats[segment][pos][col],
                        int(feats[segment][pos][col]/20)+5,
                        list(AMINO_ACID_GROUP_PROPERTIES.values())[pos]['length'],
                        list(AMINO_ACID_GROUPS.keys())[pos]
                    ])
                #print(self.feat_consensus[segment])

            del feats
            if 'feat_cons_tmp' in locals():
                del feat_cons_tmp

            self.calculate_zscales(True)
            self.stats_done = True

    def calculate_aa_count_per_generic_number(self):
        """Small function to return a dictionary of display_generic_number and the frequency of each AA."""
        generic_lookup_aa_freq = {}
        num_proteins = len(self.unique_proteins)
        for j, a in self.aa_count.items():
            for g, p in a.items():
                for aa, c in p.items():
                    if g in self.generic_number_objs:
                        if self.generic_number_objs[g].label in generic_lookup_aa_freq:
                            generic_lookup_aa_freq[self.generic_number_objs[g].label][aa] = round(c/num_proteins*100)
                        else:
                            generic_lookup_aa_freq[self.generic_number_objs[g].label] = {aa: round(c/num_proteins*100) }
        return generic_lookup_aa_freq

    def calculate_similarity(self, normalized=False):
        """Calculate the sequence identity/similarity of every selected protein compared to a selected reference."""
        for i, protein in enumerate(self.proteins):
            # skip the first row, as it is the reference
            if i == 0:
                continue

            # calculate identity, similarity and similarity score to the reference
            if not normalized:
                calc_values = self.pairwise_similarity(self.proteins[0], self.proteins[i])
                # update the protein
                if calc_values:
                    self.proteins[i].identity = calc_values[0]
                    self.proteins[i].similarity = calc_values[1]
                    self.proteins[i].similarity_score = calc_values[2]
            else:
                self.pairwise_similarity_normalized(self.proteins[0], self.proteins[i])

        # calculate normalized identity, similarity and similarity score, removes columns where reference is gapped, gaps in templates are removed from the specific pairwise alignment
        if normalized:
            i = 1
            for protein_2, values in self.normalized_scores.items():
                for j in range(0, 3):
                    for gn in list(values[j].keys()):
                        if gn in self.residues_to_delete:
                            del self.normalized_scores[protein_2][j][gn]
                identities = list(self.normalized_scores[protein_2][0].values())
                similarities = list(self.normalized_scores[protein_2][1].values())
                similarity_scores = list(self.normalized_scores[protein_2][2].values())
                identity = "{:10.0f}".format(sum(identities) / len(identities) * 100)
                similarity = "{:10.0f}".format(sum(similarities) / len(similarities) * 100)
                similarity_score = sum(similarity_scores)
                self.proteins[i].identity = identity
                self.proteins[i].similarity = similarity
                self.proteins[i].similarity_score = similarity_score
                i+=1

        # order protein list by similarity score
        ref = self.proteins.pop(0)
        order_by_value = int(getattr(self.proteins[0], self.order_by))
        if order_by_value:
            # Added secondary sorting by score
            self.proteins.sort(key=lambda x: (getattr(x, self.order_by), getattr(x, "similarity_score")), reverse=True)
        self.proteins.insert(0, ref)

    def calculate_similarity_matrix(self):
        """Calculate a matrix of sequence identity/similarity for every selected protein."""
        # Init results matrix
        self.similarity_matrix = OrderedDict()
        for i, protein in enumerate(self.proteins):
            protein_key = protein.protein.entry_name
            protein_name = "[" + protein.protein.species.common_name + "] " + protein.protein.name
            self.similarity_matrix[protein_key] = {'name': protein_name, 'values': [None] * len(self.proteins)}

        # similarity comparisons
        for i, protein in enumerate(self.proteins):
            protein_key = protein.protein.entry_name
            self.similarity_matrix[protein_key]['values'][i] = ['-', '-']

            for k in range(i+1, len(self.proteins)):
                # calculate identity, similarity and similarity score to the reference
                calc_values = self.pairwise_similarity(self.proteins[i], self.proteins[k])

                # Similarity
                value = calc_values[1].strip()
                if int(value) < 10:
                    color_class = 0
                else:
                    color_class = str(value)[:-1]
                self.similarity_matrix[self.proteins[k].protein.entry_name]['values'][i] = [value, color_class]

                # Identity
                value = calc_values[0].strip()
                if int(value) < 10:
                    color_class = 0
                else:
                    color_class = str(value)[:-1]
                self.similarity_matrix[protein_key]['values'][k] = [value, color_class]

    def calculate_zscales(self, from_stats = False):
        """Calculate Z-scales distribution for current alignment set."""
        if not self.stats_done:
            # Check if alignment statistics need to be calculated
            if len(self.aa_count) == 0 and not from_stats:
                self.calculate_statistics()
            else:
                # Prepare Z-scales per segment/GN position
                self.zscales = OrderedDict([ (zscale, OrderedDict()) for zscale in ZSCALES ])

                # Calculates distribution per GN position
                for segment in self.aa_count:
                    for zscale in ZSCALES:
                        self.zscales[zscale][segment] = OrderedDict()
                    for generic_number in self.aa_count[segment]:
                        zscale_position = { zscale: [] for zscale in ZSCALES }

                        for amino_acid in self.aa_count[segment][generic_number]:
                            if amino_acid != "-" and self.aa_count[segment][generic_number][amino_acid] > 0:
                                for key in range(len(ZSCALES)):
                                    # Frequency AA at this position * value
                                    if amino_acid in AA_ZSCALES:
                                        zscale_position[ZSCALES[key]].extend([AA_ZSCALES[amino_acid][key]] * self.aa_count[segment][generic_number][amino_acid])

                        # store average + stddev + count + display
                        for zscale in ZSCALES:
                            if len(zscale_position[zscale]) == 1:
                                display = tooltip = str(round(zscale_position[zscale][0], 2)) + " ± " + str(0) + " (1)"
                                self.zscales[zscale][segment][generic_number] = [zscale_position[zscale][0], 0, 1, display]
                            else:
                                z_mean = np.mean(zscale_position[zscale])
                                z_std = np.std(zscale_position[zscale], ddof=1)
                                z_count = len(zscale_position[zscale])
                                display = tooltip = str(round(z_mean,2)) + " ± " + str(round(z_std, 2)) + " (" + str(z_count) + ")"
                                self.zscales[zscale][segment][generic_number] = [z_mean, z_std, z_count, display]

    def evaluate_sites(self, request):
        """Evaluate which user selected site definitions match each protein sequence."""
        # get simple selection from session
        simple_selection = request.session.get('selection', False)

        # format site definititions
        site_defs = {}
        for position in simple_selection.segments:
            if position.type == 'site_residue' and position.properties['site_residue_group']:
                group_id = position.properties['site_residue_group']
                if group_id not in site_defs:
                    # min match is the value of the first item in each groups list (hence the [0])
                    # site_defs example:
                    # {
                    #     1: {
                    #         'min_match': 2,
                    #         'positions': {
                    #             {'3x51': 'hbd', '6x50': 'pos'}
                    #         },
                    #         'amino_acids': {
                    #             {'3x51': ['H', 'K', 'N', 'Q', 'R', 'S', 'T', 'W', 'Y'], '6x50': ['H', 'K', 'R']}
                    #         },
                    #     }
                    # }
                    site_defs[group_id] = {'min_match': simple_selection.site_residue_groups[group_id -1][0],
                                           'positions': {},
                                           'amino_acids': {},
                                           }

                site_defs[group_id]['positions'][position.item.label] = position.properties['feature']
                if 'amino_acids' in position.properties:
                    site_defs[group_id]['amino_acids'][position.item.label] = position.properties['amino_acids']

        # go through all proteins and match against site definitions
        for protein in self.proteins:
            for k, segment in enumerate(protein.alignment.values(), start = 1):
                num_matched = 0
                min_match = site_defs[k]['min_match']
                for position in segment:
                    # position example: ['6x49', '6.49x49', 'L', 'GPCRdb(A)', 282, 282]
                    try:
                        if position[2] in site_defs[k]['amino_acids'][position[0]]:
                            num_matched += 1
                            if num_matched >= min_match:
                                break
                    except:
                        if position[2] in AMINO_ACID_GROUPS[site_defs[k]['positions'][position[0]]]:
                            num_matched += 1
                            if num_matched >= min_match:
                                break
                else:
                    # if the protein sequence does not match the definitions, store it in non_matching_proteins
                    self.non_matching_proteins.append(protein)
                    break

        # remove non-matching proteins from protein list
        self.proteins = [p for p in self.proteins if p not in self.non_matching_proteins]

    def pairwise_similarity(self, protein_1, protein_2):
        """Calculate the identity, similarity and similarity score between a pair of proteins."""
        identityscore = 0
        similarityscore = 0
        totalcount = 0
        totalsimilarity = 0
        scoring_matrix = substitution_matrices.load("BLOSUM62")
        for j, s in protein_2.alignment.items():
            for k, p in enumerate(s):
                reference_residue = protein_1.alignment[j][k][2]
                protein_residue = protein_2.alignment[j][k][2]
                if not (reference_residue in self.gaps and protein_residue in self.gaps):
                    totalcount += 1
                    # identity
                    if protein_residue == reference_residue:
                        identityscore += 1

                    # similarity
                    if not (reference_residue in self.gaps or protein_residue in self.gaps):
                        pair = (protein_residue, reference_residue)
                        similarity = self.score_match(pair, scoring_matrix)
                        totalsimilarity += similarity
                        if similarity > 0:
                            similarityscore += 1


        # format the calculated values
        #if identityscore and similarityscore:
        if totalcount:
            identity = "{:10.0f}".format(identityscore / totalcount * 100)
            similarity = "{:10.0f}".format(similarityscore / totalcount * 100)
            similarity_score = totalsimilarity

            return identity, similarity, similarity_score
        else:
            # return False
            # NOTE returning F results in fatal errors. No aligned residues: return -1
            return "{:10.0f}".format(-1), "{:10.0f}".format(-1), 0

    def pairwise_similarity_normalized(self, protein_1, protein_2):
        """Calculate the identity, similarity and similarity score between a pair of proteins but delete gaps and normalize.

        Used for finding closest receptor homologue with crystal structure.
        """
        identities = OrderedDict()
        similarities = OrderedDict()
        similarity_scores = OrderedDict()
        ref_seq, temp_seq = '',''
        scoring_matrix = substitution_matrices.load("BLOSUM62")
        for j, s in protein_2.alignment.items():
            for k, p in enumerate(s):
                reference_residue = protein_1.alignment[j][k][2]
                protein_residue = protein_2.alignment[j][k][2]
                ref_seq+=reference_residue
                temp_seq+=protein_residue
                if not (reference_residue in self.gaps and protein_residue in self.gaps):
                    # identity
                    if protein_residue == reference_residue:
                        identities[p[0]] = 1
                    else:
                        identities[p[0]] = 0

                    # similarity
                    if reference_residue in self.gaps:
                        if p[0] not in self.residues_to_delete:
                            self.residues_to_delete.append(p[0])
                    elif protein_residue in self.gaps:
                        del identities[p[0]]
                        pass
                    else:
                        pair = (protein_residue, reference_residue)
                        similarity = self.score_match(pair, scoring_matrix)
                        if similarity > 0:
                            similarities[p[0]] = 1
                        else:
                            similarities[p[0]] = 0
                        similarity_scores[p[0]] = similarity
                else:
                    if p[0] not in self.residues_to_delete:
                        self.residues_to_delete.append(p[0])
        self.normalized_scores[protein_2] = [identities, similarities, similarity_scores]

    def score_match(self, pair, matrix):
        if pair not in matrix:
            return matrix[(tuple(reversed(pair)))]
        else:
            return matrix[pair]


class AlignedReferenceTemplate(Alignment):
    """ Creates a structure based alignment between reference protein and target proteins that are made up from the
        best available unique structures. It marks the best match as the main template structure.

        @param reference_protein: Protein object of reference protein. \n
        @param segments: list of segment ids to be considered in the alignment, e.g. ['TM1','TM2']. \n
        @param query_states: list of activation sites considered. \n
        @param order_by: str of ordering the aligned proteins. Identity, similarity or simscore.
        @param provide_main_temlpate_structure: Structure object, use only when aligning loops and when the main
        template is already known.
    """

    def __init__(self):
        super(AlignedReferenceTemplate, self).__init__()
        self.reference_dict = OrderedDict()
        self.template_dict = OrderedDict()
        self.alignment_dict = OrderedDict()
        self.code_dict = {'ICL1':'12x50','ECL1':'23x50','ICL2':'34x50'}
        self.loop_partial_except_list = {'ICL1':[],'ECL1':[],'ICL2':[],'ECL2':[],'ECL2_1':['3UZA','3UZC','3RFM','6G79'],
                                         'ECL2_mid':[],'ECL2_2':[],'ICL3':['3VW7'],'ECL3':[],'ICL4':[]}
        self.seq_nums_overwrite_cutoff_dict = {'4PHU':2000, '4LDL':1000, '4LDO':1000, '4QKX':1000, '5JQH':1000, '5TZY':2000, '5KW2':2000, '6D26':2000, '6D27':2000, '6CSY':1000}
        self.seq_num_overwrite_files = [i.split('.')[0] for i in os.listdir(os.sep.join([settings.DATA_DIR, 'structure_data','wt_pdb_lookup']))]

    def run_hommod_alignment(self, reference_protein, segments, query_states, order_by, provide_main_template_structure=None,
                             provide_similarity_table=None, main_pdb_array=None, provide_alignment=None, only_output_alignment=None,
                             complex_model=False, signprot=None, force_main_temp=False, core_alignment=None):
        self.logger = logging.getLogger('homology_modeling')
        self.segment_labels = segments
        self.complex = complex_model
        self.signprot = signprot
        self.force_main_temp = force_main_temp
        self.core_alignment = core_alignment
        self.main_temp_ban_list = ['opsd_todpa', 'adrb1_melga']
        self.gpcr_class = ProteinFamily.objects.get(name=reference_protein.get_protein_class())
        if len(str(reference_protein))==4:
            self.reference_protein = Protein.objects.get(entry_name=reference_protein.parent)
            self.revise_xtal = str(reference_protein)
        else:
            self.reference_protein = Protein.objects.get(entry_name=reference_protein)
            self.revise_xtal = None
        self.provide_alignment = provide_alignment
        if provide_main_template_structure==None and provide_similarity_table==None:
            self.query_states = query_states
            self.order_by = order_by
            self.load_reference_protein(self.reference_protein)
            if only_output_alignment!=None:
                self.load_proteins([only_output_alignment])
            else:
                self.load_proteins_by_structure()
            self.load_segments(ProteinSegment.objects.filter(slug__in=segments))
            self.build_alignment()
            self.calculate_similarity()
            self.reference_protein = self.proteins[0]
            self.main_template_protein = None
            self.ordered_proteins = []
            if only_output_alignment!=None:
                return self.proteins
        if provide_main_template_structure==None:
            self.main_template_structure = None
            self.provide_main_template_structure = False
        else:
            self.main_template_structure = provide_main_template_structure
            self.provide_main_template_structure = True
        segment_type = [str(x)[:2] for x in segments]
        if provide_similarity_table==None:
            self.provide_similarity_table = None
        else:
            self.provide_similarity_table = provide_similarity_table
        if main_pdb_array!=None:
            self.main_pdb_array = main_pdb_array
        if 'TM' in segment_type:
            self.similarity_table = self.create_helix_similarity_table()
        elif 'IC' in segment_type or 'EC' in segment_type and 'TM' not in segment_type:
            self.loop_table = OrderedDict()
            self.similarity_table = self.create_loop_similarity_table()
        if self.main_template_structure==None:
            self.changes_on_db = False
            self.main_template_structure = self.get_main_template()

    def local_pairwise_alignment(self, reference, template, segment):
        """Local pairwise alignment."""
        self.load_reference_protein(reference)
        self.load_proteins(template)
        self.load_segments(ProteinSegment.objects.get(slug__in=segment))
        self.build_alignment()
        return self.enhance_alignment(self.proteins[0], self.proteins[1])

    def __repr__(self):
        return '<AlignedReferenceTemplate: Ref: {} ; Temp: {}>'.format(self.reference_protein.protein.entry_name,
                                                                       self.main_template_structure)

    def load_proteins_by_structure(self):
        """Loads proteins into alignment based on available structures in the database."""
        if self.reference_protein.family.parent.parent.parent.slug=='003':
            template_family = ProteinFamily.objects.get(slug='002')
        elif self.reference_protein.family.parent.parent.parent.slug=='007':
            template_family = ProteinFamily.objects.get(slug='001')
        else:
            template_family = self.reference_protein.family.parent.parent.parent
        if self.reference_protein.family.parent.parent.parent.name == 'Class B2 (Adhesion)':
            self.structures_data = Structure.objects.filter(
                state__name__in=self.query_states).filter(Q(protein_conformation__protein__parent__family__parent__parent__parent=template_family) |
                                                          Q(protein_conformation__protein__parent__family__parent__parent__parent=self.reference_protein.family.parent.parent.parent)
                                                          ).order_by('protein_conformation__protein__parent','resolution').filter(annotated=True
                                                          ).exclude(structure_type__slug__startswith='af-').distinct()
        else:
            self.structures_data = Structure.objects.filter(
                state__name__in=self.query_states, protein_conformation__protein__parent__family__parent__parent__parent=
                template_family).order_by('protein_conformation__protein__parent',
                                          'resolution').filter(annotated=True
                                          ).exclude(structure_type__slug__startswith='af-').distinct()
        if self.revise_xtal==None:
            if self.force_main_temp:
                main_st = Structure.objects.get(pdb_code__index=self.force_main_temp.upper())
                if main_st.protein_conformation.protein.parent.entry_name in self.main_temp_ban_list:
                    self.main_temp_ban_list.remove(main_st.protein_conformation.protein.parent.entry_name)
            self.structures_data = self.structures_data.exclude(protein_conformation__protein__parent__entry_name__in=self.main_temp_ban_list)
        self.load_proteins([target.protein_conformation.protein.parent for target in self.structures_data])

    def get_main_template(self):
        """Returns main template structure after checking for matching helix start and end positions."""
        if self.force_main_temp:
            st = Structure.objects.get(pdb_code__index=self.force_main_temp.upper())
            # if self.core_alignment and st.pdb_code.index in self.seq_num_overwrite_files:
            #     self.overwrite_db_seq_nums(st, st.pdb_code.index)
            self.main_template_protein = [i for i in self.ordered_proteins if i.protein==st.protein_conformation.protein.parent][0]
            return st
        i = 1
        if self.complex:
            complex_templates = self.get_template_from_gprotein(self.signprot)
        for st in self.similarity_table:
            if st.pdb_code.index in ['5LWE','4Z9G'] and st.protein_conformation.protein.parent==self.ordered_proteins[i].protein:
                i+=1
                continue
            # only use complex main template in table signprot_complex
            if self.complex and st.pdb_code.index not in complex_templates:
                i+=1
                continue
            if st.protein_conformation.protein.parent==self.ordered_proteins[i].protein:
                self.main_template_protein = self.ordered_proteins[i]
                # if self.core_alignment and st.pdb_code.index in self.seq_num_overwrite_files:
                #     self.overwrite_db_seq_nums(st, st.pdb_code.index)
                return st

    def get_template_from_gprotein(self, signprot):
        gprotein = Protein.objects.get(entry_name=signprot)
        templates = SignprotComplex.objects.filter(protein=gprotein, structure__protein_conformation__protein__family__slug__startswith=self.gpcr_class.slug).exclude(beta_protein__isnull=True).values_list('structure__pdb_code__index', flat=True)
        if len(templates)==0:
            subfamily = Protein.objects.filter(family__parent=gprotein.family.parent).exclude(entry_name=gprotein.entry_name)
            templates = SignprotComplex.objects.filter(protein__in=subfamily, structure__protein_conformation__protein__family__slug__startswith=self.gpcr_class.slug).exclude(beta_protein__isnull=True).values_list('structure__pdb_code__index', flat=True)
        if len(templates)==0:
            templates = SignprotComplex.objects.all().exclude(beta_protein__isnull=True).values_list('structure__pdb_code__index', flat=True)
        return templates

    def overwrite_db_seq_nums(self, structure, cutoff):
        ssno = StructureSeqNumOverwrite(structure)
        ssno.seq_num_overwrite('pdb')
        self.changes_on_db = True
        self.logger.info('Structure {} residue table sequence number overwrite pdb to wt'.format(structure))

    def create_helix_similarity_table(self):
        """Creates an ordered dictionary of structure objects, where templates are sorted by similarity and resolution."""
        temp_list = []
        self.ordered_proteins = [self.proteins[0]]
        similarity_table = OrderedDict()
        for protein in self.proteins:
            try:
                matches = self.structures_data.filter(protein_conformation__protein__parent=protein.protein)
                for m in matches:
                    if m.protein_conformation.protein.parent==self.reference_protein.protein and int(protein.similarity)==0:
                        continue
                    if (m, int(protein.similarity), float(m.resolution), protein, m.representative) not in temp_list:
                        temp_list.append((m, int(protein.similarity), float(m.resolution), protein, m.representative))
            except:
                pass
        sorted_list = sorted(temp_list, key=lambda x: (-x[1],-x[4],x[2]))
        if self.revise_xtal!=None:
            temp_list = []
            for i in sorted_list:
                if self.revise_xtal==i[0].pdb_code.index.lower():
                    temp_list.append(i)
                    break
            for i in sorted_list:
                if self.revise_xtal!=i[0].pdb_code.index.lower():
                    temp_list.append(i)
            sorted_list = temp_list
        for i in sorted_list:
            similarity_table[i[0]] = i[1]
            self.ordered_proteins.append(i[3])
        return similarity_table

    def create_loop_similarity_table(self):
        """
        Creates an ordered dictionary of structure objects, where templates are sorted by similarity and resolution.

        Only templates that have the same loop length as the reference are considered.
        """
        temp_list, temp_list1, temp_list2, temp_list_mid = [],[],[],[]
        similarity_table = OrderedDict()
        self.main_template_protein = self.main_template_structure.protein_conformation.protein.parent
        ref_seq = Residue.objects.filter(protein_conformation__protein=self.reference_protein,
                                         protein_segment__slug=self.segment_labels[0])
        x50_ref = False
        for i in ref_seq:
            try:
                if i.generic_number.label[-3:]=='x50':
                    x50_ref = True
                    break
            except:
                pass
        prot_conf = ProteinConformation.objects.get(protein=self.reference_protein)
        segment_order = []
        for i in list(Residue.objects.filter(protein_conformation=prot_conf)):
            if i.protein_segment.slug not in segment_order:
                segment_order.append(i.protein_segment.slug)
        prev_seg = segment_order[segment_order.index(self.segment_labels[0])-1]
        next_seg = segment_order[segment_order.index(self.segment_labels[0])+1]
        if prev_seg=='C-term':
            orig_before_gns = []
        else:
            orig_before_gns = [i.replace('.','x') for i in list(self.main_pdb_array[prev_seg].keys())[-4:]]
        orig_after_gns = [j.replace('.','x') for j in list(self.main_pdb_array[next_seg].keys())[:4]]

        if len(orig_before_gns)==0:
            last_before_gn = None
        else:
            last_before_gn = orig_before_gns[-1]
        first_after_gn = orig_after_gns[0]

        if self.segment_labels[0]=='ECL2':
            try:
                ref_ECL2 = self.ECL2_slicer(ref_seq)
            except:
                ref_ECL2 = None

        for struct, similarity in self.provide_similarity_table.items():
            protein = struct.protein_conformation.protein.parent
            if protein==self.main_template_protein:
                main_temp_seq = Residue.objects.filter(protein_conformation=struct.protein_conformation,
                                                       protein_segment__slug=self.segment_labels[0])
                parent = ProteinConformation.objects.get(protein=struct.protein_conformation.protein.parent)
                main_temp_parent = Residue.objects.filter(protein_conformation=parent,
                                                          protein_segment__slug=self.segment_labels[0])
                try:
                    if self.segment_labels[0]=='ECL2' and ref_ECL2!=None:
                        main_temp_ECL2 = self.ECL2_slicer(main_temp_seq)
                        main_parent_ECL2 = self.ECL2_slicer(main_temp_parent)
                        rota = [x for x in Rotamer.objects.filter(structure=struct, residue__in=main_temp_ECL2[1]) if x.pdbdata.pdb.startswith('COMPND')==False]
                        rota = [x for x in rota if x.pdbdata.pdb[21:22] in x.structure.preferred_chain]
                        if len(rota)==3:
                            temp_list_mid.append((struct, 3, similarity, float(struct.resolution), protein, struct.representative))
                        if len(main_temp_ECL2[0])-1<=len(ref_ECL2[0])<=len(main_temp_ECL2[0])+1 and len(main_temp_ECL2[0])==len(main_parent_ECL2[0]):
                            temp_list1.append((struct, len(ref_ECL2[0]), similarity, float(struct.resolution),protein, struct.representative))
                        if len(ref_ECL2[2])==len(main_temp_ECL2[2]) and len(main_temp_ECL2[2])==len(main_parent_ECL2[2]):
                            temp_list2.append((struct, len(ref_ECL2[2]), similarity, float(struct.resolution), protein, struct.representative))

                        # Allow for partial main loop template
                        if len(main_parent_ECL2[0])-1<=len(ref_ECL2[0])<=len(main_parent_ECL2[0])+1 and [i.sequence_number for i in main_temp_ECL2[0]]!=[i.sequence_number for i in main_parent_ECL2[0]]:
                            if abs(len(main_parent_ECL2[0])-len(main_temp_ECL2[0]))<=len(main_parent_ECL2[0])/2:
                                # Ignore templates in loop partial except list
                                if struct.pdb_code.index not in self.loop_partial_except_list['ECL2_1']:
                                    temp_list1.append((struct, len(ref_ECL2[0]), 0, float(struct.resolution), protein, struct.representative))
                        if len(ref_ECL2[2])==len(main_parent_ECL2[2]) and [i.sequence_number for i in main_temp_ECL2[2]]!=[i.sequence_number for i in main_parent_ECL2[2]]:
                            if abs(len(main_parent_ECL2[2])-len(main_temp_ECL2[2]))<=len(main_parent_ECL2[2])/2:
                                temp_list2.append((struct, len(ref_ECL2[2]), 0, float(struct.resolution), protein, struct.representative))

                    else:
                        raise Exception()
                except:
                    if len(main_temp_seq)==0:
                        continue
                    if ((len(ref_seq)==len(main_temp_seq) and len(main_temp_seq)==len(main_temp_parent) and
                         [i.sequence_number for i in main_temp_seq]==[i.sequence_number for i in main_temp_parent]) or
                            self.segment_labels[0] in self.provide_alignment.reference_dict):
                        if len(main_temp_seq)!=len(main_temp_parent):
                            temp_list.append((struct, len(ref_seq), 0, float(struct.resolution), protein, struct.representative))
                        else:
                            similarity_table[self.main_template_structure] = self.provide_similarity_table[
                                self.main_template_structure]
                            temp_list.append((struct, len(main_temp_seq), similarity, float(struct.resolution), protein, struct.representative))
                    # Allow for partial main loop template
                    elif (len(ref_seq)>=len(main_temp_parent) and len(main_temp_parent)>len(main_temp_seq) and
                          [i.sequence_number for i in main_temp_seq]!=[i.sequence_number for i in main_temp_parent]):
                        if x50_ref==False and len(ref_seq)!=len(main_temp_parent):
                            continue
                        if self.segment_labels[0]=='ICL3':
                            if len(main_temp_parent)<=10:
                                temp_list.append((struct, len(ref_seq), 0, float(struct.resolution), protein, struct.representative))
                        else:
                            temp_list.append((struct, len(ref_seq), 0, float(struct.resolution), protein, struct.representative))
            else:
                temp_length, temp_length1, temp_length2 = [],[],[]
                try:
                    alt_last_gn = Residue.objects.get(protein_conformation=struct.protein_conformation,
                                                      display_generic_number__label=dgn(last_before_gn,
                                                                                        struct.protein_conformation))
                    alt_first_gn= Residue.objects.get(protein_conformation=struct.protein_conformation,
                                                      display_generic_number__label=dgn(first_after_gn,
                                                                                        struct.protein_conformation))
                    temp_length = alt_first_gn.sequence_number-alt_last_gn.sequence_number-1
                    alt_seq = Residue.objects.filter(protein_conformation=struct.protein_conformation, protein_segment__slug=self.segment_labels[0])
                    if self.segment_labels[0]=='ECL2' and ref_ECL2!=None:
                        alt_ECL2 = self.ECL2_slicer(alt_seq)
                        alt_rota = [x for x in Rotamer.objects.filter(structure=struct, residue__in=alt_ECL2[1]) if x.pdbdata.pdb.startswith('COMPND')==False]
                        if len(alt_rota)==3:
                            temp_list_mid.append((struct, 3, similarity, float(struct.resolution), protein, struct.representative))
                        if len(ref_ECL2[0])==len(alt_ECL2[0]) and len(ref_ECL2[2])==len(alt_ECL2[2]):
                            temp_length1 = len(alt_ECL2[0])
                            temp_length2 = len(alt_ECL2[2])
                        elif len(ref_ECL2[0])==len(alt_ECL2[0]) and len(ref_ECL2[2])!=len(alt_ECL2[2]):
                            temp_length1 = len(alt_ECL2[0])
                            temp_length2 = -1
                        elif len(ref_ECL2[0])!=len(alt_ECL2[0]) and len(ref_ECL2[2])==len(alt_ECL2[2]):
                            temp_length1 = -1
                            temp_length2 = len(alt_ECL2[2])
                        elif len(ref_ECL2[0])!=len(alt_ECL2[0]) and len(ref_ECL2[2])!=len(alt_ECL2[2]):
                            temp_length1 = -1
                            temp_length2 = -1
                    elif len(ref_seq)!=len(alt_seq):
                        continue
                    before_nums = list(range(alt_last_gn.sequence_number-3, alt_last_gn.sequence_number+1))
                    after_nums = list(range(alt_first_gn.sequence_number, alt_first_gn.sequence_number+4))
                    alt_before8 = Residue.objects.filter(protein_conformation__protein=protein,
                                                         sequence_number__in=before_nums)
                    alt_after8 = Residue.objects.filter(protein_conformation__protein=protein,
                                                        sequence_number__in=after_nums)
                    alt_before_gns = [ggn(r.display_generic_number.label) for r in alt_before8]
                    alt_after_gns = [ggn(r.display_generic_number.label) for r in alt_after8]
                    if orig_before_gns==alt_before_gns and orig_after_gns==alt_after_gns:
                        pass
                    else:
                        raise Exception()
                except:
                    temp_length, temp_length1, temp_length2 = -1,-1,-1
                temp_list.append((struct, temp_length, similarity, float(struct.resolution), protein, struct.representative))


                if self.segment_labels[0]=='ECL2' and ref_ECL2!=None:
                    temp_list1.append((struct, temp_length1, similarity, float(struct.resolution), protein, struct.representative))
                    temp_list2.append((struct, temp_length2, similarity, float(struct.resolution), protein, struct.representative))
        if self.segment_labels[0]=='ECL2' and ref_ECL2!=None:
            ECL2_1 = self.order_sim_table(temp_list1, ref_ECL2[0], OrderedDict(), ECL2_part='_1')
            ECL2_mid = self.order_sim_table(temp_list_mid, ref_ECL2[1], OrderedDict(), x50_ref, ECL2_part='_mid')
            ECL2_2 = self.order_sim_table(temp_list2, ref_ECL2[2], OrderedDict(), ECL2_part='_2')
            self.loop_table = OrderedDict([('ECL2_1',ECL2_1),('ECL2_mid',ECL2_mid),('ECL2_2',ECL2_2)])
            if len(ECL2_mid)==0:
                self.loop_table=None
            return self.loop_table
        else:
            # if self.segment_labels[0]!='ICL2':
            #     temp_list = temp_list[1:]
            return self.order_sim_table(temp_list, ref_seq, OrderedDict(), x50_ref)

    def order_sim_table(self, temp_list, ref_seq, similarity_table, x50_ref=None, ECL2_part=''):
        alt_temps_gn = []
        if self.segment_labels[0]!='ECL2' or self.segment_labels[0]=='ECL2' and x50_ref==True:
            for entry in temp_list:
                res_list = [i for i in list(Residue.objects.filter(protein_conformation=entry[0].protein_conformation,
                                                                   protein_segment__slug=self.segment_labels[0]))]
                for i in res_list:
                    try:
                        if i.generic_number.label[-3:]=='x50' and x50_ref==True:
                            alt_temps_gn.append(entry)
                    except:
                        pass

        alt_temps = [entry for entry in temp_list if entry[1]==len(ref_seq)]
        sorted_list_gn = sorted(alt_temps_gn, key=lambda x: (-x[2],-x[5],x[3]))
        sorted_list = sorted(alt_temps, key=lambda x: (-x[2],-x[5],x[3]))
        combined = sorted_list_gn+sorted_list
        main_in = [x for x in combined if x[0]==self.main_template_structure]
        if len(combined)>0 and combined[0][0]!=self.main_template_structure and len(main_in)>0:
            combined = main_in+combined
        if self.revise_xtal!=None:
            temp_list = []
            for i in combined:
                if self.revise_xtal==i[0].pdb_code.index.lower():
                    temp_list.append(i)
                    break
            main_t = None
            if self.revise_xtal.upper() in self.loop_partial_except_list[self.segment_labels[0]+ECL2_part]:
                main_t = combined[0]
                temp_list = []
            main_t_added = False
            for i in combined:
                if self.revise_xtal!=i[0].pdb_code.index.lower():
                    if main_t!=None and main_t_added==False and i[2]==0:
                        temp_list.append(main_t)
                        main_t_added = True
                    temp_list.append(i)
            combined = temp_list
        for i in combined:
            similarity_table[i[0]] = i[2]
        try:
            # self.main_template_protein = combined[0][4]
            # self.main_template_structure = combined[0][0]
            self.loop_table = similarity_table
        except:
            self.main_template_protein = None
            self.main_template_structure = None
            self.loop_table = None
            return None
        return similarity_table

    def ECL2_slicer(self, queryset):
        x50 = queryset.get(generic_number__label='45x50').sequence_number
        queryset_l = list(queryset)
        ECL2_1 = [i for i in queryset_l if i.sequence_number<x50]
        ECL2_mid = [i for i in queryset_l if x50<=i.sequence_number<x50+3]
        ECL2_2 = [i for i in queryset_l if i.sequence_number>=x50+3]
        if len(ECL2_mid)<3:
            raise AssertionError()
        return[ECL2_1,ECL2_mid,ECL2_2]

    def enhance_alignment(self, reference, template, keep_all=False):
        """Creates an alignment between reference and main_template where matching residues are depicted with the one-letter residue code, mismatches with '.', gaps with '-', gaps due to shorter sequences with 'x'."""
        for ref_seglab, temp_seglab in zip(reference.alignment, template.alignment):
            if 'TM' in ref_seglab or ref_seglab in ['ICL1','ECL1','ICL2','ECL2','H8']:
                ref_segment_dict,temp_segment_dict,align_segment_dict = OrderedDict(), OrderedDict(), OrderedDict()
                for ref_position, temp_position in zip(reference.alignment[ref_seglab],template.alignment[temp_seglab]):
                    ### Check res by res alignment
                    # print(ref_position, temp_position)
                    if ref_position[1]!=False and temp_position[1]!=False and ref_position[1]!='' and temp_position!='':
                        bw, gn = ref_position[1].split('x')
                        gen_num = '{}x{}'.format(bw.split('.')[0],gn)
                        ref_segment_dict[gen_num]=ref_position[2]
                        temp_segment_dict[gen_num]=temp_position[2]
                        if ref_position[2]==temp_position[2]:
                            align_segment_dict[gen_num]=ref_position[2]
                        else:
                            align_segment_dict[gen_num]='.'
                    elif ref_position[1]=='' and temp_position[1]=='':
                        ref_segment_dict[str(ref_position[4])]=ref_position[2]
                        temp_segment_dict[str(temp_position[4])]=temp_position[2]
                        if ref_position[2]==temp_position[2]:
                            align_segment_dict[str(ref_position[4])]=ref_position[2]
                        else:
                            align_segment_dict[ref_position[4]]='.'
                    elif ref_position[1]!=False and temp_position[1]==False and ref_position[1]!='':
                        bw, gn = ref_position[1].split('x')
                        gen_num = '{}x{}'.format(bw.split('.')[0],gn)
                        ref_segment_dict[gen_num]=ref_position[2]
                        if temp_position[2]=='-':
                            temp_segment_dict[gen_num]='-'
                            align_segment_dict[gen_num]='-'
                        elif temp_position[2]=='_':
                            temp_segment_dict[gen_num]='x'
                            align_segment_dict[gen_num]='x'
                    elif ref_position[1]=='' and temp_position[1]==False:
                        ref_segment_dict[str(ref_position[4])]=ref_position[2]
                        if temp_position[2]=='-':
                            temp_segment_dict[temp_position[0]]='-'
                            align_segment_dict[temp_position[0]]='-'
                        elif temp_position[2]=='_':
                            temp_segment_dict[temp_position[0]]='x'
                            align_segment_dict[temp_position[0]]='x'
                    elif ref_position[2]=='-' and temp_position[1]!=False and temp_position[1]!='':
                        bw, gn = temp_position[1].split('x')
                        gen_num = '{}x{}'.format(bw.split('.')[0],gn)
                        ref_segment_dict[gen_num]='-'
                        temp_segment_dict[gen_num]=temp_position[2]
                        align_segment_dict[gen_num]='-'
                    elif ref_position[2]=='-' and temp_position[1]=='':
                        ref_segment_dict[ref_position[0]]='-'
                        temp_segment_dict[str(temp_position[4])]=temp_position[2]
                        align_segment_dict[ref_position[0]]='-'
                    elif ref_position[2]=='_' and temp_position[1]=='':
                        ref_segment_dict[ref_position[0]]='x'
                        temp_segment_dict[str(temp_position[4])]=temp_position[2]
                        align_segment_dict[ref_position[0]]='x'
                    elif ref_position[2]=='_' and temp_position[1]!=False:
                        bw, gn = temp_position[1].split('x')
                        gen_num = '{}x{}'.format(bw.split('.')[0],gn)
                        ref_segment_dict[gen_num]='x'
                        temp_segment_dict[gen_num]=temp_position[2]
                        align_segment_dict[gen_num]='x'

                self.reference_dict[ref_seglab] = ref_segment_dict
                self.template_dict[ref_seglab] = temp_segment_dict
                self.alignment_dict[ref_seglab] = align_segment_dict

        if keep_all==False:
            delete_r = set()
            delete_t = set()
            delete_a = set()
            for r_seglab, t_seglab, a_seglab in zip(self.reference_dict,self.template_dict,self.alignment_dict):
                if r_seglab in ['ICL1','ECL1','ICL2']:
                    if len(list(self.reference_dict[r_seglab].keys()))==0:
                        well_aligned = False
                    else:
                        well_aligned = True
                        if (self.code_dict[r_seglab] not in self.reference_dict[r_seglab]
                                or self.code_dict[r_seglab] not in self.template_dict[t_seglab]):
                            well_aligned = False
                        x50_present = False
                        for i in self.reference_dict[r_seglab]:
                            if 'x50' in i and self.reference_dict[r_seglab][i]!='-':
                                x50_present = True
                                break
                        if x50_present==False:
                            well_aligned = False
                        if well_aligned==False:
                            delete_r.add(r_seglab)
                            delete_t.add(t_seglab)
                            delete_a.add(a_seglab)
                elif r_seglab=='ECL2':
                    delete_r.add(r_seglab)
                    delete_t.add(t_seglab)
                    delete_a.add(a_seglab)

            for i in delete_r:
                del self.reference_dict[i]
            for i in delete_t:
                del self.template_dict[i]
            for i in delete_a:
                del self.alignment_dict[i]

        return self


class GProteinAlignment(Alignment):
    def __init__(self):
        super(GProteinAlignment, self).__init__()
        self.reference_dict = OrderedDict()
        self.template_dict = OrderedDict()
        self.alignment_dict = OrderedDict()

    def run_alignment(self, reference_protein, template_protein, segments=None, calculate_similarity=False):
        self.load_reference_protein(reference_protein)
        self.load_proteins([template_protein])
        if segments:
            gprotein_segments = ProteinSegment.objects.filter(proteinfamily='Alpha', name__in=segments)
        else:
            gprotein_segments = ProteinSegment.objects.filter(proteinfamily='Alpha')
        self.load_segments(gprotein_segments)
        self.build_alignment()
        if calculate_similarity:
            self.calculate_similarity()
        self.enhance_alignment(self.proteins[0], self.proteins[1])

    def enhance_alignment(self, reference, template):
        for ref_seg, temp_seg in zip(reference.alignment, template.alignment):
            if ref_seg in ['S10','S11','S12','S13','S14','S15','S16','S17','S18','S19']:
                continue
            ref_segment_dict, temp_segment_dict, align_segment_dict = OrderedDict(), OrderedDict(), OrderedDict()
            for ref_pos, temp_pos in zip(reference.alignment[ref_seg], template.alignment[temp_seg]):
                if ref_pos[1] and temp_pos[1] and ref_pos[1]==temp_pos[1]:
                    ref_segment_dict[ref_pos[0]] = ref_pos[2]
                    temp_segment_dict[temp_pos[0]] = temp_pos[2]
                    if ref_pos[2]==temp_pos[2]:
                        align_segment_dict[ref_pos[0]] = ref_pos[2]
                    else:
                        align_segment_dict[ref_pos[0]] = '.'
                elif not ref_pos[1] and temp_pos[1]:
                    ref_segment_dict[ref_pos[0]] = '-'
                    temp_segment_dict[temp_pos[0]] = temp_pos[2]
                    align_segment_dict[ref_pos[0]] = '-'
                elif ref_pos[1] and not temp_pos[1]:
                    ref_segment_dict[ref_pos[0]] = ref_pos[2]
                    temp_segment_dict[temp_pos[0]] = '-'
                    align_segment_dict[ref_pos[0]] = '-'
                elif not ref_pos[1] and not temp_pos[1]:
                    ref_segment_dict[ref_pos[0]] = '-'
                    temp_segment_dict[temp_pos[0]] = '-'
                    align_segment_dict[ref_pos[0]] = '-'
            self.reference_dict[ref_seg] = ref_segment_dict
            self.template_dict[temp_seg] = temp_segment_dict
            self.alignment_dict[ref_seg] = align_segment_dict
        return self


class ClosestReceptorHomolog():
    """Finds the closest receptor homolog that has a structure. Uses the pairwise_similarity_normalized function that deletes gaps."""
    def __init__(self, protein, protein_segments=['TM1','TM2','TM3','TM4','TM5','TM6','TM7','H8'], normalized=True):
        self.protein = protein
        self.protein_segments = protein_segments
        self.normalized = normalized

        # TOFIX this is a bad workaround to select similar receptor family
        # slugs can change and the selection updates based on the available structures
        self.family_mapping = {'001':'001','002':'002','003':['002','003'],'004':'004','005':'005','006':'006','007':'001','008':['001','002','003','004','005','006']}
        self.all_proteins = []

    def find_closest_receptor_homolog(self):
        a = Alignment()
        p = Protein.objects.get(entry_name=self.protein)
        exclusion_list = ['opsd_todpa', 'g1sgd4_rabit', 'us28_hcmva', 'q08bg4_danre', 'q9wtk1_cavpo', 'q80km9_hcmv', 'q98sw5_xenla', 'b1b1u5_9arac']
        if self.protein in exclusion_list:
            exclusion_list.remove(self.protein)
        this_structs = Structure.objects.filter(protein_conformation__protein__parent__entry_name=self.protein)
        if len(this_structs)>0:
            return this_structs[0].protein_conformation.protein.parent
        else:
            if p.family.slug[:3]=='008':
                structures = Structure.objects.all().exclude(structure_type__slug__startswith='af-').exclude(annotated=False).exclude(protein_conformation__protein__parent__entry_name__in=exclusion_list)
            elif isinstance(self.family_mapping[p.family.slug[:3]], list):
                structures = []
                for slug in self.family_mapping[p.family.slug[:3]]:
                    structures+=list(Structure.objects.filter(protein_conformation__protein__parent__family__slug__istartswith=slug).exclude(
                    annotated=False).exclude(structure_type__slug__startswith='af-').exclude(protein_conformation__protein__parent__entry_name__in=exclusion_list))
            else:
                structures = Structure.objects.filter(protein_conformation__protein__parent__family__slug__istartswith=self.family_mapping[p.family.slug[:3]]).exclude(
                    annotated=False).exclude(structure_type__slug__startswith='af-').exclude(protein_conformation__protein__parent__entry_name__in=exclusion_list)
            a.load_reference_protein(p)
            structure_proteins = []
            for i in structures:
                if i.protein_conformation.protein.parent not in structure_proteins:
                    structure_proteins.append(i.protein_conformation.protein.parent)
            a.load_proteins(structure_proteins)
            a.load_segments(ProteinSegment.objects.filter(slug__in=self.protein_segments))
            a.build_alignment()
            a.calculate_similarity(normalized=self.normalized)
            self.all_proteins = a.proteins
            max_sim, max_id, max_i = 0, 0, 1
            for i, p in enumerate(self.all_proteins):
                if int(p.similarity)>max_sim:
                    max_sim = int(p.similarity)
                    max_id = int(p.identity)
                    max_i = i
                elif int(p.similarity)==max_sim and int(p.identity)>max_id:
                    max_sim = int(p.similarity)
                    max_id = int(p.identity)
                    max_i = i
            return a.proteins[max_i].protein
