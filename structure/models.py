from django.db import models

from io import StringIO
from Bio.PDB import PDBIO

class Structure(models.Model):
    protein_conformation = models.ForeignKey('protein.ProteinConformation')
    structure_type = models.ForeignKey('StructureType')
    pdb_code = models.ForeignKey('common.WebLink')
    state = models.ForeignKey('protein.ProteinState')
    publication = models.ForeignKey('common.Publication', null=True)
    ligands = models.ManyToManyField('ligand.Ligand', through='interaction.StructureLigandInteraction')
    protein_anomalies = models.ManyToManyField('protein.ProteinAnomaly')
    stabilizing_agents = models.ManyToManyField('StructureStabilizingAgent')
    preferred_chain = models.CharField(max_length=20)
    resolution = models.DecimalField(max_digits=5, decimal_places=3)
    publication_date = models.DateField()
    pdb_data = models.ForeignKey('PdbData', null=True) #allow null for now, since dump file does not contain.
    representative = models.BooleanField(default=False)


    def __str__(self):
        return self.pdb_code.index

    def get_cleaned_pdb(self, pref_chain=True, remove_waters=True, ligands_to_keep=None, remove_aux=False, aux_range=5.0):
        
        tmp = []
        for line in self.pdb_data.pdb.split('\n'):
            save_line = False
            if pref_chain:
                if (line.startswith('ATOM') or line.startswith('HET')) and line[21] == self.preferred_chain[0]:
                    save_line = True
            else:
                save_line = True
            if remove_waters and line.startswith('HET') and line[17:20] == 'HOH':
                save_line = False
            if ligands_to_keep and line.startswith('HET'):
                if line[17:20] != 'HOH' and line[17:20] in ligands_to_keep:
                    save_line = True
                elif line[17:20] != 'HOH':
                    save_line=False
            if save_line:
                tmp.append(line)

        return '\n'.join(tmp)

                        
    def get_preferred_chain_pdb(self):

        tmp = []
        for line in self.pdb_data.pdb.split('\n'):
            # http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
            if (line.startswith('ATOM') or line.startswith('HET')) and line[21] == self.preferred_chain[0]:
                tmp.append(line)
        return '\n'.join(tmp)

    class Meta():
        db_table = 'structure'


class StructureModel(models.Model):
    protein = models.ForeignKey('protein.Protein')
    state = models.ForeignKey('protein.ProteinState')
    main_template = models.ForeignKey('structure.Structure')
    pdb = models.TextField()
    version = models.DecimalField(max_digits=2, decimal_places=1)
    
    def __repr__(self):
        return '<HomologyModel: '+str(self.protein.entry_name)+' '+str(self.state)+'>'

    class Meta():
        db_table = 'structure_model'      


class StructureModelStatsRotamer(models.Model):
    homology_model = models.ForeignKey('structure.StructureModel')
    segment = models.CharField(max_length=6)
    sequence_number = models.IntegerField()
    generic_number = models.CharField(max_length=5, null=True)
    rotamer_template = models.ForeignKey('structure.Structure', null=True)
    model_segment = models.ForeignKey('structure.StructureModelStatsSegment')

    def __repr__(self):
        return '<StructureModelStatsRotamer: seqnum '+str(self.sequence_number)+' hommod '+str(self.homology_model.protein)+'>'

    class Meta():
        db_table = 'structure_model_stats_rotamer'


class StructureModelStatsSegment(models.Model):
    homology_model = models.ForeignKey('structure.StructureModel')
    segment = models.CharField(max_length=6)
    start = models.IntegerField()
    end = models.IntegerField()
    start_gn = models.CharField(max_length=5, null=True)
    end_gn = models.CharField(max_length=5, null=True)
    backbone_template = models.ForeignKey('structure.Structure', null=True)

    def __repr__(self):
        return '<StructureModelStatsSegment: '+str(self.start)+'-'+str(self.end)+' hommod '+str(self.homology_model.protein)+'>'

    class Meta():
        db_table = 'structure_model_stats_segment'


class StructureModelRMSD(models.Model):
    homology_model = models.ForeignKey('structure.StructureModel')
    pdb = models.CharField(max_length=4)
    service = models.CharField(max_length=10)
    version = models.DecimalField(max_digits=2, decimal_places=1)
    date = models.DateField(null=True)
    overall_all = models.DecimalField(null=True, max_digits=3, decimal_places=1)
    overall_backbone = models.DecimalField(null=True, max_digits=3, decimal_places=1)
    TM_all = models.DecimalField(null=True, max_digits=3, decimal_places=1)
    TM_backbone = models.DecimalField(null=True, max_digits=3, decimal_places=1)

    def __repr__(self):
        return '<StructureModelRMSD: '+self.service+' hommod '+str(self.homology_model.protein)+'>'

    class Meta():
        db_table = 'structure_model_rmsd'


class StructureType(models.Model):
    slug = models.SlugField(max_length=20, unique=True)
    name = models.CharField(max_length=100)

    def __str__(self):
        return self.name

    class Meta():
        db_table = "structure_type"


class StructureStabilizingAgent(models.Model):
    slug = models.SlugField(max_length=50, unique=True)
    name = models.CharField(max_length=100)

    def __str__(self):
        return self.name

    class Meta():
        db_table = "structure_stabilizing_agent"


class PdbData(models.Model):
    pdb = models.TextField()

    def __str__(self):
        return self.pdb

    class Meta():
        db_table = "structure_pdb_data"


class Rotamer(models.Model):
    residue = models.ForeignKey('residue.Residue')
    structure = models.ForeignKey('structure.Structure')
    pdbdata = models.ForeignKey('PdbData')

    def __str__(self):
        return '{} {}{}'.format(self.structure.pdb_code.index, self.residue.amino_acid, self.residue.sequence_number)

    class Meta():
        db_table = "structure_rotamer"


class Fragment(models.Model):
    residue = models.ForeignKey('residue.Residue')
    ligand = models.ForeignKey('ligand.Ligand')
    structure = models.ForeignKey('structure.Structure')
    pdbdata = models.ForeignKey('PdbData')

    def __str__(self):
        return '{} {}{} {}'.format(self.structure.pdb_code.index, self.residue.amino_acid,
            self.residue.sequence_number, self.ligand.name)

    class Meta():
        db_table = "structure_fragment"


class StructureSegment(models.Model):
    structure = models.ForeignKey('Structure')
    protein_segment = models.ForeignKey('protein.ProteinSegment')
    start = models.IntegerField()
    end = models.IntegerField()

    def __str__(self):
        return self.structure.pdb_code.index + " " + self.protein_segment.slug

    class Meta():
        db_table = "structure_segment"


class StructureSegmentModeling(models.Model):
    """Annotations of segment borders that are observed in exp. structures, and can be used for modeling.
    This class is indentical to StructureSegment, but is kept separate to avoid confusion."""
    structure = models.ForeignKey('Structure')
    protein_segment = models.ForeignKey('protein.ProteinSegment')
    start = models.IntegerField()
    end = models.IntegerField()

    def __str__(self):
        return self.structure.pdb_code.index + " " + self.protein_segment.slug

    class Meta():
        db_table = "structure_segment_modeling"


class StructureCoordinates(models.Model):
    structure = models.ForeignKey('Structure')
    protein_segment = models.ForeignKey('protein.ProteinSegment')
    description = models.ForeignKey('StructureCoordinatesDescription')

    def __str__(self):
        return "{} {} {}".format(self.structure.pdb_code.index, self.protein_segment.slug, self.description.text)

    class Meta():
        db_table = "structure_coordinates"


class StructureCoordinatesDescription(models.Model):
    text = models.CharField(max_length=200, unique=True)

    def __str__(self):
        return self.text

    class Meta():
        db_table = "structure_coordinates_description"


class StructureEngineering(models.Model):
    structure = models.ForeignKey('Structure')
    protein_segment = models.ForeignKey('protein.ProteinSegment')
    description = models.ForeignKey('StructureEngineeringDescription')

    def __str__(self):
        return "{} {} {}".format(self.structure.pdb_code.index, self.protein_segment.slug, self.description.text)

    class Meta():
        db_table = "structure_engineering"


class StructureEngineeringDescription(models.Model):
    text = models.CharField(max_length=200, unique=True)

    def __str__(self):
        return self.text

    class Meta():
        db_table = "structure_engineering_description"
