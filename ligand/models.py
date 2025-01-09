from django.db import models
from django.utils.text import slugify
from django.db import IntegrityError
from django.db import connection

from django_rdkit import models as rdkit_models
from django.contrib.postgres.indexes import GistIndex

from common.models import WebResource
from common.models import WebLink, Publication
from common.tools import fetch_from_web_api
from string import Template
from structure.models import Structure
from urllib.request import urlopen, quote

import json
import yaml
import logging


class Ligand(models.Model):
    name = models.TextField()
    pdbe = models.CharField(max_length=5, null=True)
    ligand_type = models.ForeignKey('LigandType', null=True, on_delete=models.CASCADE)
    ambiguous_alias = models.BooleanField(default = False)
    uniprot = models.CharField(max_length=35, null=True)

    # structure definition
    smiles = models.TextField(null=True)
    inchikey = models.CharField(max_length=27, null=True, unique=True)
    clean_inchikey = models.CharField(max_length=27, null=True)
    sequence = models.CharField(max_length=1000, null=True)

    # Ligand properties
    mw = models.DecimalField(max_digits=15, decimal_places=3, null=True)
    rotatable_bonds = models.SmallIntegerField(null=True)
    hacc = models.SmallIntegerField(null=True)
    hdon = models.SmallIntegerField(null=True)
    logp = models.DecimalField(max_digits=10, decimal_places=3, null=True)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'ligand'

class CustomLigandMolManager(models.Manager):
    def truncate_table(self):
        cursor = connection.cursor()
        table_name = self.model._meta.db_table
        sql = 'TRUNCATE TABLE "{0}" CASCADE'.format(table_name)
        cursor.execute(sql)


class LigandMol(models.Model):
    ligand = models.OneToOneField('Ligand', null=False, on_delete=models.CASCADE)
    molecule = rdkit_models.MolField()

    objects = models.Manager()  # The default manager.
    custom_objects = CustomLigandMolManager()  # The custom manager.
    class Meta:
        indexes = [
            GistIndex(fields=['molecule']),
        ]

class CustomLigandFingerprintManager(models.Manager):
    def truncate_table(self):
        cursor = connection.cursor()
        table_name = self.model._meta.db_table
        sql = 'TRUNCATE TABLE "{0}" CASCADE'.format(table_name)
        cursor.execute(sql)

class LigandFingerprint(models.Model):
    ligand = models.OneToOneField('Ligand', null=False, on_delete=models.CASCADE)
    mfp2 = rdkit_models.BfpField(null=True) #Morgan fingerprints

    objects = models.Manager()  # The default manager.
    custom_objects = CustomLigandFingerprintManager()  # The custom manager.
    class Meta:
        indexes = [
            GistIndex(fields=['mfp2']),
        ]




    # def load_by_gtop_id(self, ligand_name, gtop_id, ligand_type):
    #     logger = logging.getLogger('build')
    #
    #     # get the data from cache or web services
    #     cache_dir = ['guidetopharmacology', 'ligands']
    #     url = 'http://www.guidetopharmacology.org/services/ligands/$index'
    #     gtop = fetch_from_web_api(url, gtop_id, cache_dir)
    #
    #     if gtop:
    #         # get name from response
    #         ligand_name = gtop['name']
    #         if ligand_name == '11-<i>cis</i>-retinal':
    #             ligand_name = 'retinal'
    #
    #     # does a ligand by this name already exists?
    #     try:
    #         existing_ligand = Ligand.objects.get(
    #             name=ligand_name, canonical=True)
    #         return existing_ligand
    #     except Ligand.DoesNotExist:
    #         web_resource = False
    #
    #         if gtop_id:
    #             # gtoplig webresource
    #             web_resource = WebResource.objects.get(slug='gtoplig')
    #
    #         return self.update_ligand(ligand_name, {}, ligand_type, web_resource, gtop_id)
    #
    # def load_from_pubchem(self, lookup_type, pubchem_id, ligand_type, ligand_title=False, pdbe=None):
    #     logger = logging.getLogger('build')
    #     # if ligand title is specified, use that as the name
    #     if ligand_title:
    #         ligand_name = ligand_title
    #
    #     # otherwise, fetch ligand name from pubchem
    #     else:
    #         # check cache
    #         cache_dir = ['pubchem', 'cid', 'synonyms']
    #         url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{}/$index/synonyms/json'.format(
    #             lookup_type)
    #         pubchem = fetch_from_web_api(url, pubchem_id, cache_dir)
    #         ##print (pubchem)
    #
    #         # get name from response
    #         try:
    #             ligand_name = pubchem['InformationList']['Information'][0]['Synonym'][0]
    #         except:
    #             # Some compounds do not have a name but are still a valid pubchem entry. (Peptides)
    #             logger.warning(
    #                 'Ligand {} does not have a name in PubChem'.format(pubchem_id))
    #             ligand_name = lookup_type + ' ' + pubchem_id
    #             # return None
    #
    #     # fetch ligand properties from pubchem
    #     properties = {}
    #
    #     # check cache
    #     cache_dir = ['pubchem', 'cid', 'property']
    #     url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{}/$index/property/CanonicalSMILES,InChIKey,MolecularWeight,HBondDonorCount,HBondAcceptorCount,XLogP,RotatableBondCount/json'.format(
    #         lookup_type)
    #     pubchem = fetch_from_web_api(url, pubchem_id, cache_dir)
    #     # get properties from response
    #     if pubchem == False:
    #         logger.warning('Ligand {} not found in PubChem'.format(pubchem_id))
    #         return None
    #
    #     if pubchem['PropertyTable']['Properties'][0]:
    #         if 'HBondAcceptorCount' in pubchem['PropertyTable']['Properties'][0]:
    #             properties['hacc'] = pubchem['PropertyTable']['Properties'][0]['HBondAcceptorCount']
    #         if 'HBondDonorCount' in pubchem['PropertyTable']['Properties'][0]:
    #             properties['hdon'] = pubchem['PropertyTable']['Properties'][0]['HBondDonorCount']
    #         if 'XLogP' in pubchem['PropertyTable']['Properties'][0]:
    #             properties['logp'] = pubchem['PropertyTable']['Properties'][0]['XLogP']
    #         if 'RotatableBondCount' in pubchem['PropertyTable']['Properties'][0]:
    #             properties['rotatable_bonds'] = pubchem['PropertyTable']['Properties'][0]['RotatableBondCount']
    #         if 'MolecularWeight' in pubchem['PropertyTable']['Properties'][0]:
    #             properties['mw'] = pubchem['PropertyTable']['Properties'][0]['MolecularWeight']
    #     try:
    #
    #         properties['smiles'] = pubchem['PropertyTable']['Properties'][0]['CanonicalSMILES']
    #         properties['inchikey'] = pubchem['PropertyTable']['Properties'][0]['InChIKey']
    #
    #     except:
    #         logger.warning('Ligand {} not found in PubChem'.format(pubchem_id))
    #         return None
    #
    #     # pubchem webresource
    #     web_resource = WebResource.objects.get(slug='pubchem')
    #     #print (web_resource)
    #
    #     # does a ligand with this canonical name already exist
    #     try:
    #         return Ligand.objects.get(name=ligand_name, canonical=True)
    #         # FIXME check inchikey
    #     except Ligand.DoesNotExist:
    #         pass  # continue
    #
    #     # does a (canonical) ligand with this inchikey already exist?
    #     try:
    #         existing_lp = LigandProperities.objects.get(
    #             inchikey=properties['inchikey'])
    #         self.properities = existing_lp
    #         self.name = ligand_name
    #         self.canonical = False
    #         self.ambigious_alias = False
    #         self.pdbe = pdbe
    #
    #         try:
    #             self.save()
    #             return self
    #         except IntegrityError:
    #             return Ligand.objects.get(name=ligand_name, canonical=False)
    #     except LigandProperities.DoesNotExist:
    #         return self.update_ligand(ligand_name, properties, ligand_type, web_resource, pubchem_id, pdbe)
    #
    # def update_ligand(self, ligand_name, properties, ligand_type, web_resource=False, web_resource_index=False, pdbe=None):
    #     lp = LigandProperities.objects.create(ligand_type=ligand_type)
    #
    #     # assign properties
    #     for prop in properties:
    #         setattr(lp, prop, properties[prop])
    #
    #     # assign web link
    #     if web_resource and web_resource_index:
    #         try:
    #             wl, created = WebLink.objects.get_or_create(
    #                 index=web_resource_index, web_resource=web_resource)
    #         except IntegrityError:
    #             wl = WebLink.objects.get(
    #                 index=web_resource_index, web_resource=web_resource)
    #         lp.web_links.add(wl)
    #
    #     # try saving the properties, catch IntegrityErrors due to concurrent processing
    #     try:
    #         lp.save()
    #     except IntegrityError:
    #         lp = LigandProperities.objects.get(inchikey=properties['inchikey'])
    #
    #     self.name = ligand_name
    #     self.canonical = True
    #     self.ambigious_alias = False
    #     self.properities = lp
    #     self.pdbe = pdbe
    #
    #     try:
    #         self.save()
    #         return self
    #     except IntegrityError:
    #         return Ligand.objects.get(name=ligand_name, canonical=True)
    #
    # def load_by_name(self, name):
    #     logger = logging.getLogger('build')
    #     # fetch ligand info from pubchem - start by getting name and 'canonical' name
    #     pubchem_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/' + \
    #         name + '/synonyms/TXT'
    #     # if inchikey has been added use this -- especially to avoid updating a wrong inchikey to a synonym.
    #     if self.properities.inchikey:
    #         pubchem_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/InchiKey/' + \
    #             self.properities.inchikey + '/synonyms/TXT'
    #     try:
    #         req = urlopen(pubchem_url)
    #         pubchem = req.read().decode('UTF-8').splitlines()
    #         pubchem_name = pubchem[0]
    #     except:  # name not matched in pubchem
    #         if self.properities.inchikey:  # if inchikey has been added for check this
    #             pubchem_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/InchiKey/' + \
    #                 self.properities.inchikey + '/synonyms/TXT'
    #             try:
    #                 req = urlopen(pubchem_url)
    #                 pubchem = req.read().decode('UTF-8').splitlines()
    #                 pubchem_name = pubchem[0]
    #             except:  # name not matched in pubchem - exit cos something is wrong
    #                 logger.info(
    #                     'Ligand not found by InchiKey in pubchem: ' + str(self.properities.inchikey))
    #                 return
    #         else:  # if name not found and no inchikey, then no point in looking further
    #             logger.info(
    #                 'Ligand not found in pubchem by name (Consider renaming): ' + str(name))
    #             return
    #
    #     pubchem_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/' + \
    #         quote(pubchem_name) + '/property/CanonicalSMILES,InChIKey/json'
    #
    #     if self.properities.inchikey:
    #         pubchem_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inChiKey/' + \
    #             self.properities.inchikey + '/property/CanonicalSMILES,InChIKey/json'
    #
    #     try:
    #         req = urlopen(pubchem_url)
    #         pubchem = json.loads(req.read().decode('UTF-8'))
    #     except:  # JSON failed
    #         return
    #
    #     # weblink
    #     pubchem_id = pubchem['PropertyTable']['Properties'][0]['CID']
    #     try:
    #         web_resource = WebResource.objects.get(slug='pubchem')
    #     except:
    #         # abort if pdb resource is not found
    #         raise Exception('PubChem resource not found, aborting!')
    #
    #     pubchem_inchikey = ''
    #     pubchem_smiles = ''
    #
    #     # SMILES
    #     pubchem_smiles = pubchem['PropertyTable']['Properties'][0]['CanonicalSMILES']
    #
    #     # InChIKey
    #     pubchem_inchikey = pubchem['PropertyTable']['Properties'][0]['InChIKey']
    #
    #     try:  # now that we have inchikey, try and see if it exists in DB
    #         existing_lp = LigandProperities.objects.get(
    #             inchikey=pubchem_inchikey)
    #         self.properities = existing_lp
    #     except:
    #         wl, created = WebLink.objects.get_or_create(
    #             index=pubchem_id, web_resource=web_resource)
    #         self.properities.web_links.add(wl)
    #         # ligand type
    #         self.properities.ligand_type, created = LigandType.objects.get_or_create(
    #             slug='small-molecule', defaults={'name': 'Small molecule'})
    #         self.properities.inchikey = pubchem_inchikey
    #         self.properities.smiles = pubchem_smiles
    #         self.properities.save()
    #
    #     if pubchem_name.lower() != name.lower():  # if not canonical name
    #         logger.info("Updating canonical flag to Pubchem. PubChem canonical: " +
    #                     pubchem_name + ". DB canonical: " + name)
    #         self.canonical = False
    #         try:
    #             self.save()
    #         except IntegrityError:
    #             logger.error("FAILED SAVING LIGAND, duplicate?")
    #         # NEED TO CHECK BY INCHI KEY - SOME CANONICAL NAMES HAVE MANY ICHIKEYS (DOXEPIN)
    #         canonical_entry = Ligand.objects.filter(
    #             name=pubchem_name, properities__inchikey=pubchem_inchikey)
    #         if canonical_entry.exists():
    #             return
    #         else:  # insert the 'canonical' entry
    #             try:
    #                 canonical_entry = Ligand()
    #                 canonical_entry.name = pubchem_name
    #                 canonical_entry.canonical = True
    #                 canonical_entry.properities = self.properities
    #                 canonical_entry.save()
    #             except IntegrityError:
    #                 logger.error(
    #                     "FAILED SAVING CANONICAL LIGAND, duplicate? " + pubchem_name + " " + name)
    #                 print("FAILED SAVING CANONICAL LIGAND, duplicate? " +
    #                       pubchem_name + " " + name)

# Dedicated WebLink-like model to relieve pressure of the WL model and be more creative
class LigandID(models.Model):
    ligand = models.ForeignKey(Ligand, related_name='ids', on_delete = models.CASCADE)
    index = models.TextField(null = False)
    web_resource = models.ForeignKey(WebResource, on_delete = models.CASCADE, blank=True, null=True)

    def __str__(self):
        return Template(str(self.web_resource)).substitute(index=self.index)

class LigandImage(models.Model):
    ligand = models.ForeignKey(Ligand, related_name='image', on_delete = models.CASCADE)
    image = models.TextField(null = False)

    def __str__(self):
        # TODO convert into base64 encoding by default for in-line integration
        return self.image

class LigandType(models.Model):
    slug = models.SlugField(max_length=20, unique=True)
    name = models.CharField(max_length=100)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'ligand_type'

class LigandPeptideStructure(models.Model):
    structure = models.ForeignKey(
        'structure.Structure', on_delete=models.CASCADE, null=True)
    ligand = models.ForeignKey('ligand.Ligand', on_delete=models.CASCADE)
    chain = models.CharField(max_length=20)
    model = models.ForeignKey(
        'structure.StructureModel', on_delete=models.CASCADE, null=True)

    def __str__(self):
        return '<PeptideLigand: {} {} {}>'.format(self.structure, self.ligand, self.chain)

    class Meta():
        db_table = "ligand_peptide_structure"

class LigandRole(models.Model):
    slug = models.SlugField(max_length=50, unique=True)
    name = models.CharField(max_length=100)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'ligand_role'

class AssayExperiment(models.Model):
    ligand = models.ForeignKey('Ligand', on_delete=models.CASCADE)
    protein = models.ForeignKey('protein.Protein', on_delete=models.CASCADE)
    assay_type = models.CharField(max_length=10)
    assay_description = models.TextField(max_length=1000, null=True)
    standard_activity_value = models.CharField(max_length=20, null=True)
    p_activity_value = models.CharField(max_length=100, null=True) #Only 1 value, median/max fot GTP (p activity)
    p_activity_ranges = models.CharField(max_length=40, null=True) #If we have ranges (GtP)
    standard_relation = models.CharField(max_length=10)
    value_type = models.CharField(max_length=50, null=True)
    source = models.CharField(max_length=50, null=True)
    publication = models.ManyToManyField(Publication)
    document_chembl_id = models.CharField(max_length=100, null=True)
    affinity = models.CharField(max_length=100, null=True)
    potency = models.CharField(max_length=100, null=True)
    count_affinity_test = models.CharField(max_length=10, null=True)
    count_potency_test = models.CharField(max_length=10, null=True)
    reference_ligand = models.CharField(max_length=300, null=True)


class LigandVendors(models.Model):
    slug = models.SlugField(max_length=100, unique=True)
    name = models.CharField(max_length=200, default='')
    url = models.TextField(null=True)


class LigandVendorLink(models.Model):
    vendor = models.ForeignKey('LigandVendors', on_delete=models.CASCADE)
    ligand = models.ForeignKey('Ligand', null=True, related_name='vendors', on_delete=models.CASCADE)
    url = models.CharField(max_length=400)  # SourceRecordURL
    external_id = models.CharField(max_length=500, null=True)  # RegistryID => often does not match vendor ID


class Endogenous_GTP(models.Model):
    ligand = models.ForeignKey(Ligand, on_delete=models.CASCADE)
    ligand_species = models.ForeignKey('protein.Species', on_delete=models.CASCADE, null=True)
    ligand_action = models.ForeignKey('ligand.LigandRole', on_delete=models.CASCADE, null=True)
    endogenous_status = models.CharField(max_length=200, null=True)
    potency_ranking = models.FloatField(max_length=60, null=True)
    receptor = models.ForeignKey('protein.Protein', on_delete=models.CASCADE)
    pec50 = models.CharField(max_length=200, null=True)
    pKi = models.CharField(max_length=200, null=True)
    pic50 = models.CharField(max_length=200, null=True)
    pKd = models.CharField(max_length=200, null=True)
    publication = models.ManyToManyField(Publication)

# Biased Signalling Data
class BiasedData(models.Model):
    ligand = models.ForeignKey(Ligand, on_delete=models.CASCADE) #LINK
    publication = models.ForeignKey(Publication, on_delete=models.CASCADE) #LINK
    experiment = models.CharField(max_length=60, null=True)
    endogenous_status = models.CharField(max_length=60, null=True) #need to fetch from endogenous ligand browser now fetching from the excel datasheet (Link?)
    receptor = models.ForeignKey('protein.Protein', on_delete=models.CASCADE) #LINK
    receptor_isoform = models.CharField(max_length=60, null=True)
    active_receptor_complex = models.CharField(max_length=60, null=True)
    cell_line = models.CharField(max_length=60, null=True)
    tissue = models.CharField(max_length=60, null=True)
    species = models.CharField(max_length=60, null=True)
    primary_effector_family = models.CharField(max_length=60, null=True)
    primary_effector_subtype = models.CharField(max_length=60, null=True)
    molecule_1 = models.CharField(max_length=60, null=True)
    molecule_2 = models.CharField(max_length=60, null=True)
    measured_process = models.CharField(max_length=60, null=True)
    pathway_level = models.CharField(max_length=60, null=True)
    assay_type = models.CharField(max_length=60, null=True)
    EC50 = models.FloatField(max_length=60, null=True)
    EC50_sign = models.CharField(max_length=60, null=True)
    qualitative_activity = models.CharField(max_length=60, null=True)
    Emax = models.FloatField(max_length=60, null=True)
    Emax_sign = models.CharField(max_length=60, null=True)
    Tau_KA = models.FloatField(max_length=60, null=True)
    delta_Tau_KA = models.FloatField(max_length=60, null=True)
    time_resolved = models.CharField(max_length=60, null=True)
    pathway_preferred = models.CharField(max_length=60, null=True)  #pathway preference
    subtype_biased = models.CharField(max_length=60, null=True)     #subtype biased ligands
    physiology_biased = models.CharField(max_length=60, null=True)  #biased ligands
    pathway_biased = models.CharField(max_length=60, null=True)     #balanced ligands
    pathway_subtype_biased = models.CharField(max_length=60, null=True)     #balanced subtype ligands

class BalancedLigands(models.Model):
    ligand = models.ForeignKey(Ligand, on_delete=models.CASCADE) #LINK
    receptor = models.ForeignKey('protein.Protein', on_delete=models.CASCADE) #LINK
    first_pathway = models.CharField(max_length=60, null=True)
    second_pathway = models.CharField(max_length=60, null=True)
    delta_logEmaxEC50 = models.FloatField(max_length=60, null=True)
    delta_logTauKA = models.FloatField(max_length=60, null=True)
    subtype_balanced = models.BooleanField(default=False)
    publication = models.ForeignKey(Publication, on_delete=models.CASCADE) #LINK


# Pathways part - start
class BiasedPathways(models.Model):
    submission_author = models.CharField(max_length=50)
    ligand = models.ForeignKey(Ligand, on_delete=models.CASCADE, null=True)
    publication = models.ForeignKey(Publication, on_delete=models.CASCADE)
    receptor = models.ForeignKey('protein.Protein', on_delete=models.CASCADE)
    chembl = models.CharField(max_length=40, null=True)
    effect_type = models.CharField(max_length=200, null=True)
    relevance = models.CharField(max_length=200, null=True)
    signalling_protein = models.CharField(max_length=50, null=True)


class BiasedPathwaysAssay(models.Model):
    biased_pathway = models.ForeignKey(
        BiasedPathways, related_name='biased_pathway',
        on_delete=models.CASCADE, null=True
    )
    pathway_outcome_high = models.CharField(max_length=200)
    pathway_outcome_summary = models.CharField(max_length=300, null=True)
    pathway_outcome_detail = models.CharField(max_length=300, null=True)
    experiment_pathway_distinction = models.CharField(
        max_length=200, null=True)
    experiment_system = models.CharField(max_length=40, null=True)
    experiment_outcome_method = models.CharField(max_length=200, null=True)
