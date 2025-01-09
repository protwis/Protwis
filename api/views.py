from rest_framework import views, generics, schemas
from rest_framework.decorators import api_view, renderer_classes
from rest_framework.response import Response
from rest_framework.parsers import FileUploadParser
from rest_framework.renderers import JSONRenderer

from rest_framework_swagger.views import get_swagger_view
from rest_framework_swagger.renderers import OpenAPIRenderer, SwaggerUIRenderer

from django.template.loader import render_to_string
from django.db.models import Prefetch, Q, Min, Count

from interaction.models import ResidueFragmentInteraction
from mutation.models import MutationRaw
from protein.models import Protein, ProteinFamily, Species, ProteinSegment
from ligand.models import LigandID, AssayExperiment, Ligand, LigandPeptideStructure, Endogenous_GTP
from residue.models import Residue, ResidueGenericNumberEquivalent
from structure.models import Structure, StructureExtraProteins
from structure.assign_generic_numbers_gpcr import GenericNumbering
from structure.sequence_parser import SequenceParser
from api.serializers import (ProteinSerializer, ProteinFamilySerializer, SpeciesSerializer, ResidueSerializer,
                             ResidueExtendedSerializer, StructureLigandInteractionSerializer,
                             ComplexInteractionSerializer,
                             StructurePeptideLigandInteractionSerializer,
                             MutationSerializer, ReceptorListSerializer, GuidetoPharmacologySerializer, EndogenousLigandSerializer)
from api.renderers import PDBRenderer
from common.alignment import Alignment
from common.definitions import AMINO_ACIDS, AMINO_ACID_GROUPS
from drugs.models import Drugs
from contactnetwork.models import InteractionPeptide, Interaction, InteractingPeptideResiduePair

from io import StringIO
import numpy as np
from Bio.PDB import PDBIO, PDBParser, parse_pdb_header
from collections import OrderedDict

# FIXME add
# getMutations
# numberPDBfile

class JSONOpenAPIRenderer(OpenAPIRenderer):
    media_type = 'application/json'

@api_view()
@renderer_classes([SwaggerUIRenderer, OpenAPIRenderer, JSONOpenAPIRenderer])
def schema_view(request):
    generator = schemas.SchemaGenerator(title='API calls')
    return Response(generator.get_schema(request=request))

class ProteinDetail(generics.RetrieveAPIView):

    """
    Get a single protein instance by entry name
    \n/protein/{entry_name}/
    \n{entry_name} is a protein identifier from Uniprot, e.g. adrb2_human
    """

    queryset = Protein.objects.filter(sequence_type__slug="wt").prefetch_related('family', 'species', 'source', 'residue_numbering_scheme', 'genes')
    serializer_class = ProteinSerializer
    lookup_field = 'entry_name'


class ProteinByAccessionDetail(ProteinDetail):

    """
    Get a single protein instance by accession
    \n/protein/accession/{accession}/
    \n{accession} is a protein identifier from Uniprot, e.g. P07550
    """

    lookup_field = 'accession'


class ReceptorList(generics.ListAPIView):

    """
    Get a list of all receptor proteins in GPCRdb (source: SWISSPROT)
    \n/receptorlist/
    """

    queryset = Protein.objects.filter(accession__isnull=False, parent__isnull=True, family__slug__startswith="0", source__name="SWISSPROT").prefetch_related('family','endogenous_gtp_set')
    serializer_class = ReceptorListSerializer

class ProteinFamilyList(generics.ListAPIView):

    """
    Get a list of protein families
    \n/proteinfamily/
    """

    queryset = ProteinFamily.objects.all().prefetch_related('parent')
    serializer_class = ProteinFamilySerializer


class ProteinFamilyDetail(generics.RetrieveAPIView):

    """
    Get a single protein family instance
    \n/proteinfamily/{slug}/
    \n{slug} is a protein family identifier, e.g. 001_001_001
    """

    queryset = ProteinFamily.objects.all().prefetch_related("parent")
    serializer_class = ProteinFamilySerializer
    lookup_field = 'slug'


class ProteinFamilyChildrenList(generics.ListAPIView):

    """
    Get a list of child families of a protein family
    \n/proteinfamily/children/{slug}/
    \n{slug} is a protein family identifier, e.g. 001_001_001
    """

    serializer_class = ProteinFamilySerializer

    def get_queryset(self):
        family = self.kwargs.get('slug')
        queryset = ProteinFamily.objects.all().prefetch_related("parent")
        return queryset.filter(parent__slug=family)


class ProteinFamilyDescendantList(generics.ListAPIView):

    """
    Get a list of descendant families of a protein family
    \n/proteinfamily/descendants/{slug}/
    \n{slug} is a protein family identifier, e.g. 001_001_001
    """

    serializer_class = ProteinFamilySerializer

    def get_queryset(self):
        family = self.kwargs.get('slug')
        queryset = ProteinFamily.objects.all().prefetch_related("parent")
        return queryset.filter(Q(slug__startswith=family) & ~Q(slug=family))


class ProteinsInFamilyList(generics.ListAPIView):

    """
    Get a list of proteins in a protein family
    \n/proteinfamily/proteins/{slug}/
    \n{slug} is a protein family identifier, e.g. 001_001_001
    """

    serializer_class = ProteinSerializer

    def get_queryset(self):
        queryset = Protein.objects.all()
        family = self.kwargs.get('slug')

        return queryset.filter(sequence_type__slug='wt', family__slug__startswith=family)\
                    .prefetch_related('family', 'species', 'source', 'residue_numbering_scheme', 'genes')


class ProteinsInFamilySpeciesList(generics.ListAPIView):

    """
    Get a list of proteins in a protein family
    \n/proteinfamily/proteins/{slug}/{species}
    \n{slug} is a protein family identifier, e.g. 001_001_001
    \n{latin_name} is a species identifier from Uniprot, e.g. Homo sapiens
    """

    serializer_class = ProteinSerializer

    def get_queryset(self):
        queryset = Protein.objects.all()
        family = self.kwargs.get('slug')
        species = self.kwargs.get('latin_name')
        return queryset.filter(sequence_type__slug='wt', family__slug__startswith=family,
                               species__latin_name__iexact=species).prefetch_related('family',
                               'species', 'source', 'residue_numbering_scheme', 'genes')


class ResiduesList(generics.ListAPIView):

    """
    Get a list of residues of a protein
    \n/residues/{entry_name}/
    \n{entry_name} is a protein identifier from Uniprot, e.g. adrb2_human
    """

    serializer_class = ResidueSerializer

    def get_queryset(self):
        queryset = Residue.objects.all()
        #protein_conformation__protein__sequence_type__slug='wt',
        return queryset.filter(
            protein_conformation__protein__entry_name=self.kwargs.get('entry_name')).prefetch_related('display_generic_number','protein_segment','alternative_generic_numbers')


class ResiduesExtendedList(ResiduesList):

    """
    Get a list of residues of a protein, including alternative generic numbers
    \n/residues/extended/{entry_name}/
    \n{entry_name} is a protein identifier from Uniprot, e.g. adrb2_human
    """

    serializer_class = ResidueExtendedSerializer


class SpeciesList(generics.ListAPIView):

    """
    Get a list of species
    \n/species/
    """

    queryset = Species.objects.all()
    serializer_class = SpeciesSerializer

class GtPIDList(generics.ListAPIView):

    """
    Get a list of Guide to Pharmacology ligand IDs
    \n/ligands/gtop_ids/
    """

    queryset = LigandID.objects.filter(web_resource_id__slug="gtoplig").values('index', 'ligand__name')
    serializer_class = GuidetoPharmacologySerializer

class EndogenousLigands(generics.ListAPIView):

    """
    Get a list of endogenous ligands - receptor pairs
    \n/ligands/endogenousligands/
    """

    queryset = Endogenous_GTP.objects.all().prefetch_related('ligand_id', 'ligand_id_ligand_type_id',
                                                             'receptor_id','receptor_id__family_id').values('endogenous_status', 'potency_ranking',
                                                                                                            'ligand_id__name', 'ligand_id__sequence',
                                                                                                            'ligand_id__ligand_type_id__slug', 'receptor_id__entry_name')
    serializer_class = EndogenousLigandSerializer


class PeptideList(views.APIView):

    """
    Get a list of peptide ligands with detailed info
    \n/ligands/peptides/
    """
    def get(self, request):
        peptides_ids = Ligand.objects.filter(ligand_type_id=2).values_list('id', flat=True)
        peptides_data = Ligand.objects.filter(ligand_type_id=2).values_list('id', flat=True).values_list('name', 'sequence')
        peptides_with_structures = LigandPeptideStructure.objects.filter(ligand_id__in=peptides_ids).prefetch_related('structure', 'structure_id__web_link',
                                                                                           'structure__protein_conformation', 'structure__protein_conformation__protein',
                                                                                           'structure__protein_conformation__protein__family').values_list(
                                                                                           'ligand_id__name', 'ligand_id__sequence', 'structure_id__pdb_code_id__index',
                                                                                           'structure_id__protein_conformation_id__protein_id__parent_id__entry_name',
                                                                                           'structure_id__protein_conformation_id__protein_id__family__parent__name',
                                                                                           'structure_id__protein_conformation_id__protein_id__family__parent__parent__parent__name',
                                                                                           'chain')
        s = []
        for peptide in peptides_with_structures:
            peptide_info = {
                'Peptide name': peptide[0],
                'Sequence': peptide[1],
                'Sequence length': len(peptide[1]) if peptide[1] is not None else None,
                'PDB': peptide[2],
                'Chain': peptide[6],
                'Receptor': peptide[3],
                'Family': peptide[4],
                'GPCR Class': peptide[5]
            }
            s.append(peptide_info)
        for peptide in peptides_data:
            peptide_info = {
                'Peptide name': peptide[0],
                'Sequence': peptide[1],
                'Sequence length': len(peptide[1]) if peptide[1] != None else None,
                'PDB': '-',
                'Chain': '-',
                'Receptor': '-',
                'Family': '-',
                'GPCR Class': '-'
            }
            s.append(peptide_info)


        return Response(s)

class SpeciesDetail(generics.RetrieveAPIView):

    """
    Get a single species instance
    \n/species/{latin_name}/
    \n{latin_name} is a species identifier from Uniprot, e.g. Homo sapiens
    """

    queryset = Species.objects.all()
    serializer_class = SpeciesSerializer
    lookup_field = 'latin_name'


class NumberPDBStructureView(views.APIView):

    """
    WRITEME
    """

    pass

class StructureModelsList(views.APIView):

    """
    Get a list of GPCR - peptide alphafold model structures
    \n/structure/alphafold_peptide_models
    """

    def get(self, request, pdb_code=None, entry_name=None, representative=None):

        structures = Structure.objects.filter(structure_type__slug__startswith='af-peptide')

        structures = structures.prefetch_related('protein_conformation__protein__parent__species', 'pdb_code',
            'protein_conformation__protein__parent__family', 'protein_conformation__protein__parent__species',
            'protein_conformation__protein__parent__family__parent__parent__parent',
            'publication__web_link', 'publication__web_link__web_resource', 'structure_type',
            'structureligandinteraction_set__ligand',
            'structureligandinteraction_set__ligand__ligand_type',
            'structureligandinteraction_set__ligand_role','state')

        # convert objects to a list of dictionaries
        # normal serializers can not be used because of abstraction of tables (e.g. protein_conformation)
        s = []
        for structure in structures:
            # essential fields
            structure_data = {
                'pdb_code': structure.pdb_code.index,
                'protein': structure.protein_conformation.protein.entry_name,
                'class': structure.protein_conformation.protein.family.parent.parent.parent.name,
                'family': structure.protein_conformation.protein.family.slug,
                'species': structure.protein_conformation.protein.species.latin_name,
                'preferred_chain': structure.preferred_chain,
                'resolution': structure.resolution,
                'publication_date': structure.publication_date,
                'type': structure.structure_type.name,
                'state': structure.state.name,
                'distance': structure.distance,
            }

            # publication
            if structure.publication:
                structure_data['publication'] = structure.publication.web_link.__str__()
            else:
                structure_data['publication'] = None

            # ligand
            ligands = []
            #for interaction in structure.structureligandinteraction_set.filter(annotated=True): # does this cancel prefetch?
            for interaction in structure.structureligandinteraction_set.all():
                if interaction.annotated:
                    ligand = {}
                    if interaction.ligand.name:
                        ligand['name'] = interaction.ligand.name
                    if interaction.ligand.ligand_type and interaction.ligand.ligand_type.name:
                        ligand['type'] = interaction.ligand.ligand_type.name
                    if interaction.ligand_role and interaction.ligand_role.name:
                        ligand['function'] = interaction.ligand_role.name
                    if interaction.ligand.pdbe:
                        ligand['PDB'] = interaction.ligand.pdbe
                    if interaction.ligand.smiles:
                        ligand['SMILES'] = interaction.ligand.smiles
                    if ligand:
                        ligands.append(ligand)
            structure_data['ligands'] = ligands

            s.append(structure_data)

        # if a structure is selected, return a single dict rather then a list of dicts
        if len(s) == 1:
            s = s[0]

        return Response(s)


class StructureList(views.APIView):

    """
    Get a list of structures
    \n/structure/
    """

    def get(self, request, pdb_code=None, entry_name=None, representative=None):
        if pdb_code:
            structures = Structure.objects.filter(pdb_code__index=pdb_code)
        elif entry_name and representative:
            structures = Structure.objects.filter(protein_conformation__protein__parent__entry_name=entry_name,
                representative=True).exclude(structure_type__slug__startswith='af-')
        elif entry_name:
            structures = Structure.objects.filter(protein_conformation__protein__parent__entry_name=entry_name).exclude(structure_type__slug__startswith='af-')
        elif representative:
            structures = Structure.objects.filter(representative=True).exclude(structure_type__slug__startswith='af-')
        else:
            structures = Structure.objects.all().exclude(structure_type__slug__startswith='af-')

        structures = structures.prefetch_related('protein_conformation__protein__parent__species', 'pdb_code',
            'protein_conformation__protein__parent__family', 'protein_conformation__protein__parent__species',
            'protein_conformation__protein__parent__family__parent__parent__parent',
            'publication__web_link', 'publication__web_link__web_resource', 'structure_type',
            'structureligandinteraction_set__ligand',
            'structureligandinteraction_set__ligand__ligand_type',
            'structureligandinteraction_set__ligand_role',
            'signprot_complex', 'signprot_complex__protein', 'signprot_complex__protein__parent',
            'signprot_complex__beta_protein', 'signprot_complex__gamma_protein',
            'signprot_complex__protein__parent__family', 'signprot_complex__protein__parent__species', 'state',
            Prefetch('extra_proteins', queryset=StructureExtraProteins.objects.filter(wt_protein__family__slug__startswith="200")),
            'extra_proteins__wt_protein')


        # convert objects to a list of dictionaries
        # normal serializers can not be used because of abstraction of tables (e.g. protein_conformation)
        s = []
        for structure in structures:
            # essential fields
            structure_data = {
                'pdb_code': structure.pdb_code.index,
                'protein': structure.protein_conformation.protein.parent.entry_name,
                'class': structure.protein_conformation.protein.parent.family.parent.parent.parent.name,
                'family': structure.protein_conformation.protein.parent.family.slug,
                'species': structure.protein_conformation.protein.parent.species.latin_name,
                'preferred_chain': structure.preferred_chain,
                'resolution': structure.resolution,
                'publication_date': structure.publication_date,
                'type': structure.structure_type.name,
                'state': structure.state.name,
                'distance': structure.distance,
            }

            # publication
            if structure.publication:
                structure_data['publication'] = structure.publication.web_link.__str__()
            else:
                structure_data['publication'] = None

            # ligand
            ligands = []
            #for interaction in structure.structureligandinteraction_set.filter(annotated=True): # does this cancel prefetch?
            for interaction in structure.structureligandinteraction_set.all():
                if interaction.annotated:
                    ligand = {}
                    if interaction.ligand.name:
                        ligand['name'] = interaction.ligand.name
                    if interaction.ligand.ligand_type and interaction.ligand.ligand_type.name:
                        ligand['type'] = interaction.ligand.ligand_type.name
                    if interaction.ligand_role and interaction.ligand_role.name:
                        ligand['function'] = interaction.ligand_role.name
                    if interaction.pdb_reference and interaction.pdb_reference != "pep":
                        ligand['PDB'] = interaction.pdb_reference
                    elif interaction.ligand.pdbe:
                        ligand['PDB'] = interaction.ligand.pdbe
                    if interaction.ligand.smiles:
                        ligand['SMILES'] = interaction.ligand.smiles
                    if ligand:
                        ligands.append(ligand)
            structure_data['ligands'] = ligands

            # signalling protein
            if structure.signprot_complex:
                sign_prot = {'type': 'G protein', 'data': {}}
                sign_prot['data']['entity1'] = {'entry_name':structure.signprot_complex.protein.entry_name, 'chain':structure.signprot_complex.alpha}
                if structure.signprot_complex.beta_protein:
                    sign_prot['data']['entity2'] = {'entry_name':structure.signprot_complex.beta_protein.entry_name, 'chain':structure.signprot_complex.beta_chain}
                if structure.signprot_complex.gamma_protein:
                    sign_prot['data']['entity3'] = {'entry_name':structure.signprot_complex.gamma_protein.entry_name, 'chain':structure.signprot_complex.gamma_chain}
                structure_data['signalling_protein'] = sign_prot
            if len(structure.extra_proteins.all())>0:
                for arrestin in structure.extra_proteins.all():
                    structure_data['signalling_protein'] = {'type': 'Arrestin', 'data': {'entity1':{'entry_name':arrestin.wt_protein.entry_name, 'chain':arrestin.chain}}}

            s.append(structure_data)

        # if a structure is selected, return a single dict rather then a list of dicts
        if len(s) == 1:
            s = s[0]

        return Response(s)

    def get_structures(self, pdb_code=None, representative=None):
        return Structure.objects.all().exclude(structure_type__slug__startswith='af-')

class RepresentativeStructureList(StructureList):

    """
    Get a list of representative structures (one for each protein and activation state)
    \n/structure/representative/
    """


class StructureListProtein(StructureList):

    """
    Get a list of structures of a protein
    \n/structure/protein/{entry_name}
    """


class RepresentativeStructureListProtein(StructureList):

    """
    Get a list of representative structures of a protein (one for each activation state)
    \n/structure/protein/{entry_name}/representative/
    """


class StructureDetail(StructureList):

    """
    Get a single structure instance
    \n/structure/{pdb_code}/
    \n{pdb_code} is a structure identifier from the Protein Data Bank, e.g. 2RH1
    """
    def get_structures(self, pdb_code=None, representative=None):
        return Structure.objects.filter(pdb_code__index=pdb_code)

class StructureAccessionHuman(views.APIView):

    """
    Get a list of all (human) UniProt accession codes for which receptor types an experimental structure is available
    \n/structure/accession_codes_human/
    """

    def get(self, request):
        unique_slugs = list(Structure.objects.filter(protein_conformation__protein__family__slug__startswith="0").exclude(structure_type__slug__startswith='af-')\
            .order_by('protein_conformation__protein__family__slug').values_list('protein_conformation__protein__family__slug', flat=True).distinct())
        accession_codes = list(Protein.objects.filter(family__slug__in=unique_slugs, sequence_type__slug='wt', species__latin_name='Homo sapiens')\
                            .values_list('entry_name', flat=True))
        accessions = [x.split("_")[0].upper() for x in accession_codes]
        return Response(accessions)

class FamilyAlignment(views.APIView):

    """
    Get a full sequence alignment of a protein family including a consensus sequence.
    Note that this method only includes Swiss-Prot sequences in the alignment.
    \n/alignment/family/{slug}/
    \n{slug} is a protein family identifier, e.g. 001_001_001
    """

    def get(self, request, slug=None, segments=None, latin_name=None, statistics=False, include_trembl=False):
        if slug is not None:
            # Check for specific species
            if latin_name is not None:
                ps = Protein.objects.filter(sequence_type__slug='wt', family__slug__startswith=slug,
                    species__latin_name__iexact=latin_name)
            else:
                if not include_trembl:
                    ps = Protein.objects.filter(sequence_type__slug='wt', family__slug__startswith=slug, source__id=1)
                else:
                    ps = Protein.objects.filter(sequence_type__slug='wt', family__slug__startswith=slug)

            # take the numbering scheme from the first protein
            #s_slug = Protein.objects.get(entry_name=ps[0]).residue_numbering_scheme_id
            s_slug = ps[0].residue_numbering_scheme_id

            protein_family = ps[0].family.slug[:3]

            gen_list = []
            segment_list = []
            if segments is not None:
                input_list = segments.split(",")
                # fetch a list of all segments

                protein_segments = ProteinSegment.objects.filter(partial=False).values_list('slug', flat=True)
                for s in input_list:
                    # add to segment list
                    if s in protein_segments:
                        segment_list.append(s)
                    # get generic numbering object for generic positions
                    else:
                        # make sure the query works for all positions
                        gen_object = ResidueGenericNumberEquivalent.objects.get(label=s, scheme__id=s_slug)
                        gen_object.properties = {}
                        gen_list.append(gen_object)

                # fetch all complete protein_segments
                ss = ProteinSegment.objects.filter(slug__in=segment_list, partial=False)
            else:
                ss = ProteinSegment.objects.filter(partial=False)

            if int(protein_family) < 100:
                ss = [ s for s in ss if s.proteinfamily == 'GPCR']
            elif protein_family == "100":
                ss = [ s for s in ss if s.proteinfamily == 'Gprotein']
            elif protein_family == "200":
                ss = [ s for s in ss if s.proteinfamily == 'Arrestin']

            # create an alignment object
            a = Alignment()
            a.show_padding = False

            # load data from selection into the alignment
            a.load_proteins(ps)

            # load generic numbers and TMs seperately
            if gen_list:
                a.load_segments(gen_list)
            a.load_segments(ss)

            # build the alignment data matrix
            a.build_alignment()

            a.calculate_statistics()

            residue_list = []
            for aa in a.full_consensus:
                residue_list.append(aa.amino_acid)

            # render the fasta template as string
            response = render_to_string('alignment/alignment_fasta.html', {'a': a}).split("\n")

            # convert the list to a dict
            ali_dict = OrderedDict({})
            protein = "reference"
            for row in response:
                if row.strip() != "":
                    if row.startswith(">"):
                        protein = row[1:]
                    else:
                        ali_dict[protein] = row
                        protein = False

            ali_dict["CONSENSUS"] = "".join(residue_list)

            # render statistics for output
            if statistics == True:
                feat = {}
                for i, feature in enumerate(AMINO_ACID_GROUPS):
                    feature_stats = a.feature_stats[i]
                    feature_stats_clean = []
                    for d in feature_stats:
                        sub_list = [x[0] for x in d]
                        feature_stats_clean.append(sub_list) # remove feature frequencies
                    # print(feature_stats_clean)
                    feat[feature] = [item for sublist in feature_stats_clean for item in sublist]

                for i, AA in enumerate(AMINO_ACIDS):
                    feature_stats = a.amino_acid_stats[i]
                    feature_stats_clean = []
                    for d in feature_stats:
                        sub_list = [x[0] for x in d]
                        feature_stats_clean.append(sub_list) # remove feature frequencies
                    # print(feature_stats_clean)
                    feat[AA] = [item for sublist in feature_stats_clean for item in sublist]

                ali_dict["statistics"] = feat

            return Response(ali_dict)

class FamilyAlignmentAll(FamilyAlignment):

    """
    Get a full sequence alignment of a protein family including both Swiss-Prot
    and TrEMBL sequences. Note that this method allows for the alignment of up
    to a few hundred sequences, larger alignments will result in an error.
    \n/alignment/family_all/{slug}
    \n{slug} is a protein family identifier, e.g. 001_001_001_001
    """

class FamilyAlignmentPartial(FamilyAlignment):

    """
    Get a partial sequence alignment of a protein family
    \n/alignment/family/{slug}/{segments}/
    \n{slug} is a protein family identifier, e.g. 001_001_001
    \n{segments} is a comma separated list of protein segment identifiers and/ or
    generic GPCRdb numbers, e.g. TM2,TM3,ECL2,4x50
    """

class FamilyAlignmentSpecies(FamilyAlignment):

    """
    Get a full sequence alignment of a protein family
    \n/alignment/family/{slug}//{species}
    \n{slug} is a protein family identifier, e.g. 001_001_001
    \n{species} is a species identifier from Uniprot, e.g. Homo sapiens
    """

class FamilyAlignmentPartialSpecies(FamilyAlignment):

    """
    Get a partial sequence alignment of a protein family
    \n/alignment/family/{slug}/{segments}/{species}
    \n{slug} is a protein family identifier, e.g. 001_001_001
    \n{segments} is a comma separated list of protein segment identifiers and/ or
    generic GPCRdb numbers, e.g. TM2,TM3,ECL2,4x50
    \n{species} is a species identifier from Uniprot, e.g. Homo sapiens
    """


class ProteinSimilaritySearchAlignment(views.APIView):

    """
    Get a segment sequence alignment of two or more proteins ranked by similarity
    \n/alignment/similarity/{proteins}/{segments}/
    \n{proteins} is a comma separated list of protein identifiers, e.g. adrb2_human,5ht2a_human,cxcr4_human,
    where the first protein is the query protein and the following are compared to it. Now PDB IDs can also be
    used to align structure sequences with wild type sequences, e.g. adrb2_human,3SN6
    \n{segments} is a comma separated list of protein segment identifiers and/ or
    generic GPCRdb numbers, e.g. TM2,TM3,ECL2,4x50
    """

    def get(self, request, proteins=None, segments=None):
        if proteins is not None:
            protein_list = proteins.split(",")
            # first in API should be reference
            reference = Protein.objects.filter(entry_name__iexact=protein_list[0])
            q_list = Q()
            for q in [Q(entry_name__iexact=p) for p in protein_list[1:]]:
                q_list |= q
            ps = Protein.objects.filter(q_list)

            # take the numbering scheme from the first protein
            s_slug = reference[0].residue_numbering_scheme_id

            protein_family = ps[0].family.slug[:3]

            gen_list = []
            segment_list = []
            if segments is not None:
                input_list = segments.split(",")
                # fetch a list of all segments
                protein_segments = ProteinSegment.objects.filter(partial=False).values_list('slug', flat=True)
                for s in input_list:
                    # add to segment list
                    if s in protein_segments:
                        segment_list.append(s)
                    # get generic numbering object for generic positions
                    else:
                        # make sure the query works for all positions
                        gen_object = ResidueGenericNumberEquivalent.objects.get(label=s, scheme__id=s_slug)
                        gen_object.properties = {}
                        gen_list.append(gen_object)

                # fetch all complete protein_segments
                ss = ProteinSegment.objects.filter(slug__in=segment_list, partial=False)
            else:
                ss = ProteinSegment.objects.filter(partial=False)

            if int(protein_family) < 100:
                ss = [ s for s in ss if s.proteinfamily == 'GPCR']
            elif protein_family == "100":
                ss = [ s for s in ss if s.proteinfamily == 'Alpha']
            elif protein_family == "200":
                ss = [ s for s in ss if s.proteinfamily == 'Arrestin']

            # create an alignment object
            a = Alignment()
            a.show_padding = False

            # load data from API into the alignment
            a.load_reference_protein(reference[0])
            a.load_proteins(ps)

            # load generic numbers and TMs seperately
            if gen_list:
                a.load_segments(gen_list)
            a.load_segments(ss)

            # build the alignment data matrix
            a.build_alignment()

            # calculate identity and similarity of each row compared to the reference
            a.calculate_similarity()

            # render the fasta template as string
            response = render_to_string('alignment/alignment_fasta.html', {'a': a}).split("\n")

            # convert the list to a dict
            ali_dict = {}
            k = False
            num = 0
            for i, row in enumerate(response):
                if row.startswith(">"):
                    k = row[1:]
                elif k:
                    # add the query as 100 identical/similar to the beginning (like on the website)
                    if num == 0:
                        a.proteins[num].identity = 100
                        a.proteins[num].similarity = 100
                    # order dict after custom list
                    keyorder = ["similarity","identity","AA"]
                    ali_dict[k] = {"AA": row, "identity": int(str(a.proteins[num].identity).replace(" ","")),
                    "similarity": int(str(a.proteins[num].similarity).replace(" ",""))}
                    ali_dict[k] = OrderedDict(sorted(ali_dict[k].items(), key=lambda t: keyorder.index(t[0])))
                    num+=1
                    k = False
            ali_dict_ordered = OrderedDict(sorted(ali_dict.items(), key=lambda x: x[1]['similarity'], reverse=True))
            return Response(ali_dict_ordered)

class ProteinAlignment(views.APIView):

    """
    Get a full sequence alignment of two or more proteins
    \n/alignment/protein/{proteins}/
    \n{proteins} is a comma separated list of protein identifiers, e.g. adrb2_human,5ht2a_human. PDB IDs can also be
    used to align structure sequences with wild type sequences, e.g. adrb2_human,3SN6
    """

    def get(self, request, proteins=None, segments=None, statistics=False):
        if proteins is not None:
            protein_list = proteins.split(",")
            q_list = Q()
            for q in [Q(entry_name__iexact=p) for p in protein_list]:
                q_list |= q
            ps = Protein.objects.filter(q_list)

            # take the numbering scheme from the first protein
            #s_slug = Protein.objects.get(entry_name=protein_list[0]).residue_numbering_scheme_id
            s_slug = ps[0].residue_numbering_scheme_id

            protein_family = ps[0].family.slug[:3]

            gen_list = []
            segment_list = []
            if segments is not None:
                input_list = segments.split(",")

                # fetch a list of all segments
                protein_segments = ProteinSegment.objects.filter(partial=False).values_list('slug', flat=True)
                for s in input_list:
                    # add to segment list
                    if s in protein_segments:
                        segment_list.append(s)
                    # get generic numbering object for generic positions
                    else:
                        gen_object = ResidueGenericNumberEquivalent.objects.get(label=s, scheme__id=s_slug)
                        gen_object.properties = {}
                        gen_list.append(gen_object)

                # fetch all complete protein_segments
                ss = ProteinSegment.objects.filter(slug__in=segment_list, partial=False)
            else:
                ss = ProteinSegment.objects.filter(partial=False)

            if int(protein_family) < 100:
                ss = [ s for s in ss if s.proteinfamily == 'GPCR']
            elif protein_family == "100":
                ss = [ s for s in ss if s.proteinfamily == 'Alpha']
            elif protein_family == "200":
                ss = [ s for s in ss if s.proteinfamily == 'Arrestin']

            # create an alignment object
            a = Alignment()
            a.show_padding = False

            # load data from selection into the alignment
            a.load_proteins(ps)

            # load generic numbers and TMs seperately
            if gen_list:
                a.load_segments(gen_list)
            a.load_segments(ss)

            # build the alignment data matrix
            a.build_alignment()

            # calculate statistics
            if statistics == True:
                a.calculate_statistics()

            # render the fasta template as string
            response = render_to_string('alignment/alignment_fasta.html', {'a': a}).split("\n")

            # convert the list to a dict
            ali_dict = {}
            k = False
            for row in response:
                if row.startswith(">"):
                    k = row[1:]
                elif k:
                    ali_dict[k] = row
                    k = False

            # render statistics for output
            if statistics == True:
                feat = {}
                for i, feature in enumerate(AMINO_ACID_GROUPS):
                    feature_stats = a.feature_stats[i]
                    feature_stats_clean = []
                    for d in feature_stats:
                        sub_list = [x[0] for x in d]
                        feature_stats_clean.append(sub_list) # remove feature frequencies
                    # print(feature_stats_clean)
                    feat[feature] = [item for sublist in feature_stats_clean for item in sublist]

                for i, AA in enumerate(AMINO_ACIDS):
                    feature_stats = a.amino_acid_stats[i]
                    feature_stats_clean = []
                    for d in feature_stats:
                        sub_list = [x[0] for x in d]
                        feature_stats_clean.append(sub_list) # remove feature frequencies
                    # print(feature_stats_clean)
                    feat[AA] = [item for sublist in feature_stats_clean for item in sublist]

                ali_dict["statistics"] = feat

            return Response(ali_dict)

class ProteinAlignmentStatistics(ProteinAlignment):

    """
    Add a /statics at the end of an alignment in order to
    receive an additional residue property statistics output e.g.:
    \n/alignment/protein/{proteins}/{segments}/statistics
    \n{proteins} is a comma separated list of protein identifiers, e.g. adrb2_human,5ht2a_human
    \n{segments} is a comma separated list of protein segment identifiers and/ or
    generic GPCRdb numbers, e.g. TM2,TM3,ECL2,4x50
    """

class ProteinAlignmentPartial(ProteinAlignment):

    """
    Get a partial sequence alignment of two or more proteins
    \n/alignment/protein/{proteins}/{segments}/
    \n{proteins} is a comma separated list of protein identifiers, e.g. adrb2_human,5ht2a_human
    \n{segments} is a comma separated list of protein segment identifiers and/ or
    generic GPCRdb numbers, e.g. TM2,TM3,ECL2,4x50
    """


class StructureTemplate(views.APIView):

    """
    Get the most similar structure template for a protein using a 7TM alignment
    \n/structure/template/{entry_name}/
    \n{entry_name} is a protein identifier from Uniprot, e.g. adrb2_human
    """

    def get(self, request, entry_name=None, segments=None):
        if entry_name is not None:
            ref = Protein.objects.get(sequence_type__slug='wt', entry_name=entry_name)

            structures =  Structure.objects.order_by('protein_conformation__protein__parent', 'state',
                'resolution').distinct('protein_conformation__protein__parent', 'state')

            ps = []
            for structure in structures:
                ps.append(structure.protein_conformation.protein.parent)

            if segments is not None:
                input_list = segments.split(",")
                ss = ProteinSegment.objects.filter(slug__in=input_list, partial=False)
            else:
                ss = ProteinSegment.objects.filter(partial=False, category='helix')

            # create an alignment object
            a = Alignment()
            a.show_padding = False

            # load data from selection into the alignment
            a.load_reference_protein(ref)
            a.load_proteins(ps)
            a.load_segments(ss)

            # build the alignment data matrix
            a.build_alignment()

            # calculate identity and similarity of each row compared to the reference
            a.calculate_similarity()

            # return the entry_name of the closest template
            return Response(a.proteins[1].protein.entry_name)


class StructureTemplatePartial(StructureTemplate):

    """
    Get the most similar structure template for a protein using a partial alignment
    \n/structure/template/{entry_name}/{segments}/
    \n{entry_name} is a protein identifier from Uniprot, e.g. adrb2_human
    \n{segments} is a comma separated list of protein segment identifiers, e.g. TM3,TM5,TM6
    """


class StructureAssignGenericNumbers(views.APIView):

    """
    Assign generic residue numbers (Ballesteros-Weinstein and GPCRdb schemes) to an uploaded pdb file.
    \n/structure/assign_generic_numbers\n
    e.g.
    curl -X POST -F "pdb_file=@myfile.pdb" https://gpcrdb.org/services/structure/assign_generic_numbers
    """
    parser_classes = (FileUploadParser,)
    renderer_classes = (PDBRenderer, )

    def post(self, request):

        # root, ext = os.path.splitext(request._request.FILES['pdb_file'].name)
        generic_numbering = GenericNumbering(StringIO(request._request.FILES['pdb_file'].file.read().decode('UTF-8', "ignore")))
        out_struct = generic_numbering.assign_generic_numbers()
        out_stream = StringIO()
        io = PDBIO()
        io.set_structure(out_struct)
        io.save(out_stream)
        print(len(out_stream.getvalue()))
        # filename="{}_GPCRdb.pdb".format(root)
        return Response(out_stream.getvalue())


class StructureSequenceParser(views.APIView):

    """
    Analyze the uploaded pdb structure listing auxiliary proteins, mutations, deletions and insertions.
    \n/structure/structure/parse_pdb\n
    e.g.
    curl -X POST -F "pdb_file=@myfile.pdb" https://gpcrdb.org/services/structure/parse_pdb
    """
    parser_classes = (FileUploadParser,)
    renderer_classes = (JSONRenderer, )

    def post(self, request):
        # root, ext = os.path.splitext(request._request.FILES['pdb_file'].name)
        pdb_file = StringIO(request._request.FILES['pdb_file'].file.read().decode('UTF-8', "ignore"))
        header = parse_pdb_header(pdb_file)
        parser = SequenceParser(pdb_file)

        json_data = OrderedDict()
        json_data["header"] = header
        json_data.update(parser.get_fusions())
        json_data.update(parser.get_mutations())
        json_data.update(parser.get_deletions())

        return Response(json_data)


class StructureLigandInteractions(generics.ListAPIView):

    """
    Get a list of interactions between structure and ligand
    \n/structure/{pdb_code}/interaction/
    \n{pdb_code} is a structure identifier from the Protein Data Bank, e.g. 2RH1
    """
    serializer_class = StructureLigandInteractionSerializer

    def get_queryset(self):
        pdb_code = self.kwargs.get('pdb_code')

        queryset = ResidueFragmentInteraction.objects.filter(structure_ligand_pair__structure__pdb_code__index__iexact=pdb_code,
                                                             structure_ligand_pair__annotated=True).prefetch_related('structure_ligand_pair__structure__pdb_code',
                                                             'interaction_type',
                                                             'fragment__residue__generic_number',
                                                             'fragment__residue__display_generic_number',
                                                             ).order_by('fragment__residue__sequence_number')
        return queryset


class ComplexInteractions(generics.ListAPIView):
    """
    Get a list of interactions between structure and G protein/GPCR complex
    \n/structure/{id}/interaction/
    \n{id} is a structure identifier from either the Protein Data Bank, e.g. 5UZ7,
    \na modified PDB identifier to get the refined structure interactions, e.g. 5UZ7_refined,
    \nor a GproteinDb AlphaFold2-Multimer complex model identifier, e.g. AFM_CALCR_HUMAN_GNAS2_HUMAN
    """
    serializer_class = ComplexInteractionSerializer

    def get_queryset(self):
        struct_id = self.kwargs.get('id')

        queryset = Interaction.objects.filter(
            interacting_pair__referenced_structure__pdb_code__index=struct_id,
            interacting_pair__res2_id__protein_conformation__protein__family__slug__startswith='100'
        ).values(
            'interacting_pair__res1_id__protein_conformation__protein__entry_name',
            'interacting_pair__res1_id__protein_conformation__protein',
            'interacting_pair__res1_id__protein_conformation__protein__parent',
            'interacting_pair__res1_id__sequence_number',
            'interacting_pair__res1_id__display_generic_number__label',
            'interacting_pair__res2_id__protein_conformation__protein',
            'interacting_pair__res2_id__protein_conformation__protein__parent',
            'interacting_pair__res2_id__sequence_number',
            'interacting_pair__res2_id__display_generic_number__label',
            'interaction_type',
            'interaction_level'
        ).distinct().order_by('interacting_pair__res2__display_generic_number')

        return queryset


class StructurePeptideLigandInteractions(generics.ListAPIView):

    """
    Get a list of interactions between structure and peptide ligand
    \n/structure/{value}/peptideinteraction/
    \n{value} can be a structure identifier from the Protein Data Bank, e.g. 5VBL
    \n{value} can also be a protein identifier from Uniprot, e.g. cxcr4_human
    \n{value} can also be a protein identifier from Uniprot, e.g. P61073
    \nThe inserted value will be queried in the following order: PDB code --> UniProt entry name --> UniProt accession
    \nBy default, UniProt values (entry name and accession) will be queried to AlphaFold Models interaction data
    """
    serializer_class = StructurePeptideLigandInteractionSerializer

    def get_queryset(self):
        value = self.kwargs.get('value')
        #trying different inputs: pdb_code, entry_name, accession_number
        queryset = InteractionPeptide.objects.filter(interacting_peptide_pair__peptide__structure__pdb_code__index__iexact=value)
        if len(queryset) == 0:
            queryset = InteractionPeptide.objects.filter(interacting_peptide_pair__peptide__structure__protein_conformation__protein__entry_name__iexact=value)
        if len(queryset) == 0: #check also the parent if it's a PDB structure
            queryset = InteractionPeptide.objects.filter(interacting_peptide_pair__peptide__structure__protein_conformation__protein__parent__entry_name__iexact=value)
        if len(queryset) == 0:
            queryset = InteractionPeptide.objects.filter(interacting_peptide_pair__peptide__structure__protein_conformation__protein__accession__iexact=value)
        if len(queryset) == 0: #check also the parent if it's a PDB structure
            queryset = InteractionPeptide.objects.filter(interacting_peptide_pair__peptide__structure__protein_conformation__protein__parent__accession__iexact=value)

        queryset = queryset.values('interacting_peptide_pair__peptide__structure__pdb_code__index',
                            'interacting_peptide_pair__peptide__ligand__name',
                            'interacting_peptide_pair__peptide__chain',
                            'interacting_peptide_pair__peptide_amino_acid',
                            'interacting_peptide_pair__peptide_amino_acid_three_letter',
                            'interacting_peptide_pair__peptide_sequence_number',
                            'interacting_peptide_pair__receptor_residue__amino_acid',
                            'interacting_peptide_pair__receptor_residue__sequence_number',
                            'interacting_peptide_pair__receptor_residue__display_generic_number__label',
                            'interacting_peptide_pair_id',
                            'interaction_type', 'interaction_level').order_by("interacting_peptide_pair__peptide_sequence_number").distinct(
                            ).annotate(
                                interaction_count=Count('interaction_type')
                            ).order_by('interacting_peptide_pair__peptide_sequence_number','interacting_peptide_pair__receptor_residue__sequence_number')

        return queryset

class MutantList(generics.ListAPIView):

    """
    Get a list of mutants of single protein instance by entry name
    \n/mutant/{entry_name}/
    \n{entry_name} is a protein identifier from Uniprot, e.g. adrb2_human
    """
    serializer_class = MutationSerializer

    def get_queryset(self):
        queryset = MutationRaw.objects.all()
        return queryset.filter(protein=self.kwargs.get('entry_name'))

class DrugList(views.APIView):

    """
    Get a list of drugs for a single protein instance by entry name
    \n/drugs/{proteins}/
    \n{entry_name} is a protein identifier from Uniprot, e.g. adrb2_human
    """
    def get(self, request, entry_name=None):

        drugs = Drugs.objects.filter(target__entry_name=entry_name).distinct()

        druglist = []
        for drug in drugs:
            drugname = drug.name
            drugtype = drug.drugtype
            clinical = drug.clinicalstatus
            phasedate = drug.phasedate
            if clinical != '-':
                status = drug.status + ' (' + drug.clinicalstatus + ', ' + phasedate + ')'
            else:
                status = drug.status
            approval = drug.approval
            indication = drug.indication
            moa = drug.moa
            novelty = drug.novelty
            druglist.append({'name':drugname, 'approval': approval, 'indication': indication, 'status':status, 'drugtype':drugtype, 'moa':moa, 'novelty': novelty})

        return Response(druglist)

class LigandList(views.APIView):

    """
    Get a list of ligands for a single protein instance by entry name
    \n/ligands/{prot_name}/
    \n{prot_name} is a protein identifier from Uniprot, e.g. adrb2_human
    """
    def get(self, request, prot_name=None):
        ligands = AssayExperiment.objects.filter(protein__entry_name=prot_name).prefetch_related('ligand','protein','publication').order_by('ligand_id__name')
        ligands = ligands.values('value_type',
                                 'standard_relation',
                                 'standard_activity_value',
                                 'assay_description',
                                 'assay_type',
                                 'p_activity_value',
                                 'p_activity_ranges',
                                 'source',
                                 'ligand__name',
                                 'ligand__ligand_type__name',
                                 'protein__species__common_name',
                                 'protein__entry_name',
                                 'protein__name',
                                 'ligand__rotatable_bonds',
                                 'ligand__smiles',
                                 'protein',
                                 'publication__web_link__index',
                                 'publication__web_link__web_resource__url',
                                 'affinity',
                                 'potency')

        ligandlist = []
        for compound in ligands:
            protein = compound['protein__name']
            lig_name = compound['ligand__name']
            lig_type = compound['ligand__ligand_type__name'].replace('-',' ')
            smiles = compound['ligand__smiles']
            p_activity_ranges = compound['p_activity_ranges']
            p_activity_value = compound['p_activity_value']
            value_type = compound['value_type']
            assay_type = compound['assay_type']
            assay_description = compound['assay_description']
            standard_value = compound['standard_activity_value']
            source = compound['source']
            if compound['publication__web_link__index'] is not None:
                DOI = compound['publication__web_link__web_resource__url'].strip('$index') + compound['publication__web_link__index']
            else:
                DOI = 'Not Available'

            ligandlist.append({'Protein name':protein,
                               'Ligand name': lig_name,
                               'Ligand type': lig_type,
                               'Smiles': smiles,
                               'Activity ranges (P)': p_activity_ranges,
                               'Activity value (P)': p_activity_value,
                               'Activity value (Standard)': standard_value,
                               'Value type': value_type,
                               'Assay type': assay_type,
                               'Assay description': assay_description,
                               'Source': source,
                               'DOI': DOI})

        return Response(ligandlist)

class HelixBoxView(views.APIView):

    """
    Get SVG source code for a protein's helix box plot
    \n/plot/helixbox/{entry_name}/
    \n{entry_name} is a protein identifier from Uniprot, e.g. adrb2_human
    """

    def get(self, request, entry_name=None):
        if entry_name is not None:
            p = Protein.objects.get(entry_name=entry_name)

            return Response(str(p.get_helical_box()).split("\n"))


class SnakePlotView(views.APIView):

    """
    Get SVG source code for a protein's snake plot
    \n/plot/snake/{entry_name}/
    \n{entry_name} is a protein identifier from Uniprot, e.g. adrb2_human
    """

    def get(self, request, entry_name=None):
        if entry_name is not None:
            p = Protein.objects.get(entry_name=entry_name)

            return Response(str(p.get_snake_plot()).split("\n"))

class PeptideInteractionCADistances(views.APIView):

    """
    Returns the distance between the C-alpha of receptor residue and C-alpha of peptide ligand
    \nfrom the given interaction id
    \n/structure/peptideinteraction/ca_distances/{interaction_id}/
    \n{interaction_id} is a peptide-ligand interaction ID derived from GPCRdb API
    """

    def get(self, request, interaction_id=None):
        one_to_three = {
            "A": "ALA",
            "R": "ARG",
            "N": "ASN",
            "D": "ASP",
            "C": "CYS",
            "Q": "GLN",
            "E": "GLU",
            "G": "GLY",
            "H": "HIS",
            "I": "ILE",
            "L": "LEU",
            "K": "LYS",
            "M": "MET",
            "F": "PHE",
            "P": "PRO",
            "S": "SER",
            "T": "THR",
            "W": "TRP",
            "Y": "TYR",
            "V": "VAL",
        }
        if interaction_id is not None:
            interaction_record = InteractingPeptideResiduePair.objects.filter(id=test_id)[0]
            conformation_id = interaction_record.receptor_residue.protein_conformation.id
            peptide_amino, peptide_pos, receptor_residue = interaction_record.peptide_amino_acid_three_letter, interaction_record.peptide_sequence_number, interaction_record.receptor_residue
            struct = Structure.objects.filter(protein_conformation=conformation_id)[0]

            # Create a StringIO object to simulate reading from a file-like object
            pdb_io = StringIO(struct.pdb_data.pdb)

            # Define the residue names and residue numbers for the two residues
            residue1_name = one_to_three[receptor_residue.amino_acid]
            residue1_number = receptor_residue.sequence_number
            residue2_name = peptide_amino
            residue2_number = peptide_pos

            # Create a PDB parser
            parser = PDBParser(QUIET=True)

            # Parse the PDB file
            pdb_struct = parser.get_structure("protein", pdb_io)

            # Iterate through the structure and find the Cα atoms of the specified residues
            ca1 = None
            ca2 = None
            cb = None
            for model in pdb_struct:
                for chain in model:
                    for res in chain:
                        if res.get_resname() == residue1_name and res.get_id()[1] == residue1_number:
                            ca1 = res["CA"].get_coord()
                            cb = res["CB"].get_coord()
                        elif res.get_resname() == residue2_name and res.get_id()[1] == residue2_number:
                            ca2 = res["CA"].get_coord()


            # Check if both Cα atoms were found
            if ca1 is not None and ca2 is not None:
                # Calculate the Euclidean distance between the Cα atoms
                distance = np.linalg.norm(ca1 - ca2)
                # Calculate vectors for ca1-ca2 (ligand Cα to receptor Cα) and ca2-CB (ligand Cα to receptor Carbon-Beta)
                vector_ca1_ca2 = ca2 - ca1
                vector_ca2_CB = cb - ca2
                # Calculate the dot product between the vectors
                dot_product = np.dot(vector_ca1_ca2, vector_ca2_CB)
                # Calculate the magnitudes (norms) of the vectors
                magnitude_ca1_ca2 = np.linalg.norm(vector_ca1_ca2)
                magnitude_ca2_CB = np.linalg.norm(vector_ca2_CB)
                # Calculate the cosine of the angle between the vectors using the dot product and magnitudes
                cosine_angle = dot_product / (magnitude_ca1_ca2 * magnitude_ca2_CB)
                # Calculate the angle in radians using the arccosine function (cosine rule)
                angle_radians = np.arccos(cosine_angle)
                # Convert the angle from radians to degrees
                angle_degrees = np.degrees(angle_radians)
                result = [{'Ca distance': distance,
                          'Ca-Ca-CB angle (degrees)': angle_degrees,
                          'Ca-Ca-CB angle (radians)': angle_radians}]
                return Response(result)
