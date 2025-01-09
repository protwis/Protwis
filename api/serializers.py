from rest_framework import serializers

from interaction.models import ResidueFragmentInteraction
from ligand.models import Endogenous_GTP, LigandID
from mutation.models import MutationRaw
from protein.models import Protein, ProteinConformation, ProteinFamily, Species, ProteinSource, ProteinSegment
from residue.models import Residue, ResidueNumberingScheme, ResidueGenericNumber
from structure.models import Structure, StructureComplexModel
from contactnetwork.models import InteractionPeptide, Interaction


class ProteinSerializer(serializers.ModelSerializer):
    family = serializers.SlugRelatedField(read_only=True, slug_field='slug')
    species = serializers.StringRelatedField(read_only=True)
    source = serializers.StringRelatedField(read_only=True)
    residue_numbering_scheme = serializers.SlugRelatedField(read_only=True, slug_field='short_name')
    genes = serializers.StringRelatedField(many=True)
    class Meta:
        model = Protein
        fields = ('entry_name', 'name', 'accession', 'family', 'species', 'source', 'residue_numbering_scheme',
            'sequence','genes')


class ProteinFromConformationSerializer(serializers.ModelSerializer):
    protein = serializers.SlugRelatedField(read_only=True, slug_field='entry_name')
    class Meta:
        model = ProteinConformation
        fields = ('protein', )


class ParentProteinFamilySerializer(serializers.ModelSerializer):
    class Meta:
        model = ProteinFamily
        fields = ('slug', 'name')


class ProteinFamilySerializer(serializers.ModelSerializer):
    parent = ParentProteinFamilySerializer(read_only=True)
    class Meta:
        model = ProteinFamily
        fields = ('slug', 'name', 'parent')

class EndogenousGTPSerializer(serializers.ModelSerializer):
    name = serializers.ReadOnlyField(source='ligand.name')

    class Meta:
        model = Endogenous_GTP
        fields = ('name', )

class ReceptorListSerializer(serializers.ModelSerializer):
    subfamily = serializers.ReadOnlyField(source='family.name')
    ligand_type = serializers.ReadOnlyField(source='family.parent.parent.name')
    receptor_family = serializers.ReadOnlyField(source='family.parent.name')
    receptor_class = serializers.ReadOnlyField(source='family.parent.parent.parent.name')
    endogenous_ligands = EndogenousGTPSerializer(source='endogenous_gtp_set', read_only=True, many=True)
    species = serializers.StringRelatedField(read_only=True)

    class Meta:
        model = Protein
        fields = ('entry_name', 'name', 'accession', 'receptor_class', 'receptor_family', 'ligand_type', 'subfamily', 'receptor_class', 'endogenous_ligands', 'species', 'sequence')


class SpeciesSerializer(serializers.ModelSerializer):
    class Meta:
        model = Species
        fields = ('latin_name', 'common_name')

class GuidetoPharmacologySerializer(serializers.ModelSerializer):
    gtp_ligand_id = serializers.ReadOnlyField(source='index')
    ligand_name = serializers.ReadOnlyField(source='ligand__name')
    class Meta:
        model = LigandID
        fields = ('gtp_ligand_id', 'ligand_name')

class EndogenousLigandSerializer(serializers.ModelSerializer):
    endogenous_status = serializers.ReadOnlyField()
    potency_ranking = serializers.ReadOnlyField()
    ligand_name = serializers.ReadOnlyField(source='ligand_id__name')
    sequence = serializers.ReadOnlyField(source='ligand_id__sequence')
    ligand_type = serializers.ReadOnlyField(source='ligand_id__ligand_type_id__slug')
    receptor = serializers.ReadOnlyField(source='receptor_id__entry_name')
    class Meta:
        model = Endogenous_GTP
        fields = ('receptor', 'ligand_name', 'sequence', 'sequence', 'ligand_type', 'endogenous_status', 'potency_ranking')

class ResidueNumberingSchemeSerializer(serializers.ModelSerializer):
    class Meta:
        model = ResidueNumberingScheme
        fields = ('slug', 'short_name')


class ResidueGenericNumberSerializer(serializers.ModelSerializer):
    scheme = serializers.SlugRelatedField(read_only=True, slug_field='short_name')
    class Meta:
        model = ResidueGenericNumber
        fields = ('scheme', 'label')


class ResidueSerializer(serializers.ModelSerializer):
    protein_segment = serializers.StringRelatedField(read_only=True)
    display_generic_number = serializers.StringRelatedField(read_only=True)
    class Meta:
        model = Residue
        fields = ('sequence_number', 'amino_acid', 'protein_segment', 'display_generic_number')


class ResidueExtendedSerializer(serializers.ModelSerializer):
    protein_segment = serializers.StringRelatedField(read_only=True)
    display_generic_number = serializers.StringRelatedField(read_only=True)
    alternative_generic_numbers = ResidueGenericNumberSerializer(read_only=True, many=True)
    class Meta:
        model = Residue
        fields = ('sequence_number', 'amino_acid', 'protein_segment', 'display_generic_number',
            'alternative_generic_numbers')


class StructureSerializer(serializers.ModelSerializer):
    pdb_code = serializers.SlugRelatedField(read_only=True, slug_field='index')
    protein_conformation = ProteinFromConformationSerializer()
    class Meta:
        model = Structure
        fields = ('pdb_code', 'resolution', 'protein_conformation')


class StructureLigandInteractionSerializer(serializers.ModelSerializer):
    pdb_code = serializers.ReadOnlyField(source='structure_ligand_pair.structure.pdb_code.index')
    ligand_name = serializers.ReadOnlyField(source='structure_ligand_pair.ligand.name')
    amino_acid = serializers.ReadOnlyField(source='fragment.residue.amino_acid')
    sequence_number = serializers.IntegerField(read_only=True, source='fragment.residue.sequence_number')
    display_generic_number = serializers.ReadOnlyField(source='fragment.residue.display_generic_number.label')
    interaction_type = serializers.ReadOnlyField(source='interaction_type.name')
    class Meta:
        model = ResidueFragmentInteraction
        fields = ('pdb_code', 'ligand_name', 'amino_acid', 'sequence_number',
                  'display_generic_number', 'interaction_type')


class StructurePeptideLigandInteractionSerializer(serializers.ModelSerializer):
    pdb_code = serializers.ReadOnlyField(source='interacting_peptide_pair__peptide__structure__pdb_code__index')
    ligand_name = serializers.ReadOnlyField(source='interacting_peptide_pair__peptide__ligand__name')
    ligand_chain = serializers.ReadOnlyField(source='interacting_peptide_pair__peptide__chain')
    peptide_amino_acid = serializers.ReadOnlyField(source='interacting_peptide_pair__peptide_amino_acid')
    peptide_amino_acid_three_letter = serializers.ReadOnlyField(source='interacting_peptide_pair__peptide_amino_acid_three_letter')
    peptide_residue_number = serializers.ReadOnlyField(source='interacting_peptide_pair__peptide_sequence_number')
    receptor_amino_acid = serializers.ReadOnlyField(source='interacting_peptide_pair__receptor_residue__amino_acid')
    receptor_residue_number = serializers.ReadOnlyField(source='interacting_peptide_pair__receptor_residue__sequence_number')
    receptor_residue_generic_number = serializers.ReadOnlyField(source='interacting_peptide_pair__receptor_residue__display_generic_number__label')
    interaction_type = serializers.ReadOnlyField()
    interaction_level = serializers.ReadOnlyField()
    interaction_count = serializers.ReadOnlyField()
    interacting_peptide_pair_id = serializers.ReadOnlyField()

    class Meta:
        model = InteractionPeptide
        fields = ('pdb_code', 'ligand_name', 'ligand_chain', 'peptide_amino_acid', 'peptide_amino_acid_three_letter', 'peptide_residue_number',
                  'receptor_amino_acid', 'receptor_residue_number', 'receptor_residue_generic_number', 'interacting_peptide_pair_id', 'interaction_type', 'interaction_level', 'interaction_count')


class ComplexInteractionSerializer(serializers.ModelSerializer):
    id = serializers.ReadOnlyField(
        source='interacting_pair__referenced_structure__pdb_code__index')

    receptor = serializers.ReadOnlyField(
        source='interacting_pair__res1_id__protein_conformation__protein__parent__entry_name')
    receptor_generic_residue_number = serializers.ReadOnlyField(source='interacting_pair__res1_id__display_generic_number__label')
    receptor_residue_number = serializers.ReadOnlyField(source='interacting_pair__res1_id__sequence_number')

    gprotein = serializers.ReadOnlyField(
        source='interacting_pair__res1_id__protein_conformation__protein__parent__entry_name')
    gprotein_generic_residue_number = serializers.ReadOnlyField(source='interacting_pair__res2_id__display_generic_number__label')
    gprotein_residue_number = serializers.ReadOnlyField(source='interacting_pair__res2_id__sequence_number')

    interaction_type = serializers.ReadOnlyField()
    interaction_level = serializers.ReadOnlyField()

    class Meta:
        model = Interaction
        fields = ('id',
                  'receptor', 'receptor_generic_residue_number', 'receptor_residue_number',
                  'gprotein', 'gprotein_generic_residue_number', 'gprotein_residue_number',
                  'interaction_type', 'interaction_level')


class MutationSerializer(serializers.ModelSerializer):
    class Meta:
        model = MutationRaw
        fields = ('reference', 'protein', 'mutation_pos', 'mutation_from', 'mutation_to',
            'ligand_name', 'ligand_idtype', 'ligand_id', 'ligand_class',
            'exp_type', 'exp_func',  'exp_wt_value',  'exp_wt_unit','exp_mu_effect_sign', 'exp_mu_effect_type', 'exp_mu_effect_value',
            'exp_fold_change',
            'exp_mu_effect_qual', 'exp_mu_effect_ligand_prop',  'exp_mu_ligand_ref', 'opt_receptor_expression', 'opt_basal_activity',
            'opt_gain_of_activity', 'opt_ligand_emax', 'opt_agonist')
