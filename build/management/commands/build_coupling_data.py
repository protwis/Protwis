from django.conf import settings
from django.core.management.base import BaseCommand
from django.db.models import Q

from build.management.commands.build_ligand_functions import get_or_create_ligand
from build.management.commands.parse_excel_annotations import Command as ParseExcel
from common.models import Publication, WebLink, WebResource
from common.tools import test_model_updates
from protein.models import Protein, ProteinFamily, ProteinCouplings, Biosensor
from ligand.models import Ligand

import csv
import logging
import os
import sys
from itertools import islice
from collections import defaultdict
import django.apps
import xlrd


class Command(BaseCommand):
    help = 'Build G proteins'

    # source files
    iupharcoupling_file = os.sep.join([settings.DATA_DIR, 'g_protein_data', 'iuphar_coupling_data.csv'])
    master_file = os.sep.join([settings.DATA_DIR, 'g_protein_data', 'GPCR-G_protein_couplings.xlsx'])
    #Setting the variables for the test tracking of the model upadates
    tracker = {}
    all_models = django.apps.apps.get_models()[6:]
    test_model_updates(all_models, tracker, initialize=True)
    logger = logging.getLogger(__name__)

    def handle(self, *args, **options):
        self.purge_coupling_data()
        self.logger.info('PASS: purge_coupling_data')
        self.create_iuphar_couplings()
        test_model_updates(self.all_models, self.tracker, check=True)
        self.logger.info('PASS: create_iuphar_couplings')
        self.create_data_couplings()
        test_model_updates(self.all_models, self.tracker, check=True)
        self.logger.info('PASS: create_data_couplings')

    def purge_coupling_data(self):
        """DROP data from the protein_gprotein_pair table."""
        try:
            ProteinCouplings.objects.all().delete()
        except Exception as msg:
            self.logger.warning('Existing protein couplings cannot be deleted' + str(msg))

    def create_iuphar_couplings(self, filenames=False):
        """
        Function to add IUPHAR receptor-G protein couplings.

        The function reads a iupharcoupling_file, which comes from parsing the Guide_to_Pharmacology.
        """
        self.logger.info('CREATING: IUPHAR couplings')
        # print("PROCESSING: IUPHAR couplings")

        translation = {'Gs family': ['Gs', '100_001_001'],
                       'Gi/Go family': ['Gi/o', '100_001_002'],
                       'Gq/G11 family': ['Gq/11', '100_001_003'],
                       'G12/G13 family': ['G12/13', '100_001_004'],
                       'GPa1 family': ['GPa1 family', '100_001_005']
                       }

        filenames = [self.iupharcoupling_file]

        source = "GuideToPharma"
        for filename in filenames:
            filepath = self.iupharcoupling_file
            self.logger.info('Reading filename ' + filename)
            pub_years = defaultdict(int)
            pub_years_protein = defaultdict(set)

            with open(filepath, 'r') as f:
                reader = csv.reader(f)
                for row in islice(reader, 1, None):  # skip first line
                    entry_name = row[3]
                    primary = row[11]
                    secondary = row[12]
                    primary_pubmed = row[15]
                    secondary_pubmed = row[16]
                    # fetch protein
                    try:
                        p = Protein.objects.get(entry_name=entry_name)
                    except Protein.DoesNotExist:
                        self.logger.warning('Protein not found for entry_name {}'.format(entry_name))
                        print("protein not found for ", entry_name)
                        continue

                    primary = primary.replace("G protein (identity unknown)", "None")  # replace none
                    primary = primary.split(", ")

                    secondary = secondary.replace("G protein (identity unknown)", "None")  # replace none
                    secondary = secondary.split(", ")

                    if primary == 'None' and secondary == 'None':
                        print('no data for ', entry_name)
                        continue

                    try:
                        for gp in primary:
                            if gp in ['', 'None', '_-arrestin', 'Arrestin',
                                      'G protein independent mechanism']:  # skip bad ones
                                continue
                            g = ProteinFamily.objects.get_or_create(name=translation[gp][0], slug=translation[gp][1])[0]
                            gpair = ProteinCouplings(protein=p, g_protein=g, transduction='primary', source=source)
                            gpair.save()

                            for pmid in primary_pubmed.split("|"):
                                try:
                                    test = int(pmid)
                                except:
                                    continue
                                try:
                                    pub = Publication.objects.get(web_link__index=pmid,
                                                                  web_link__web_resource__slug='pubmed')
                                except Publication.DoesNotExist as e:
                                    pub = Publication()
                                    pub.update_from_pubmed_data(index=pmid)
                                    pub.save()
                                gpair.references.add(pub)
                                gpair.save()

                    except Exception as e:
                        print("error in primary assignment", p, gp, e)

                    try:
                        for gp in secondary:
                            if gp in ['None', '_-arrestin', 'Arrestin', 'G protein independent mechanism',
                                      '']:  # skip bad ones
                                continue
                            if gp in primary:  # sip those that were already primary
                                continue
                            g = ProteinFamily.objects.get_or_create(name=translation[gp][0], slug=translation[gp][1])[0]
                            gpair = ProteinCouplings(protein=p, g_protein=g, transduction='secondary', source=source)
                            gpair.save()

                            for pmid in secondary_pubmed.split("|"):
                                try:
                                    test = int(pmid)
                                except:
                                    continue

                                try:
                                    pub = Publication.objects.get(web_link__index=pmid,
                                                                  web_link__web_resource__slug='pubmed')
                                except Publication.DoesNotExist as e:
                                    pub = Publication()
                                    pub.update_from_pubmed_data(index=pmid)
                                    pub.save()
                                gpair.references.add(pub)
                                gpair.save()
                    except Exception as e:
                        print("error in secondary assignment", p, gp, e)


        self.logger.info('COMPLETED CREATING IUPHAR couplings')

    def assess_ligand_id(self, spreadsheet):
        """
        Function that fetches the ligands data from the master file for the different sources
        and compare them to the LigandIDs sheet. A dictionary is created that contains the information
        to be further parsed via the get_or_make_ligand function. Data from Bouvier and Inoue, processed
        by David Gloriam.
        """
        #rading data by the ligand sheets
        Bligands = spreadsheet.sheet_by_name("B-Ligands")
        Iligands = spreadsheet.sheet_by_name("I-Ligands")
        IDs = spreadsheet.sheet_by_name("LigandIDs")
        #counting sheets rows
        Blig_rows = Bligands.nrows
        Ilig_rows = Iligands.nrows
        IDs_rows = IDs.nrows

        #third value is the column number bearing the ligand name
        #bouvier --> ligand renamed - Bouvier (9)
        #inoue   --> ligand short (author)    (2)
        tupled_sources = [(Bligands, Blig_rows, 9),(Iligands, Ilig_rows, 2)]

        def parseValue(s, type='int'):
            """
            Function to return an int or NA from the ligand data values
            """
            if s == '':
                return 'NA'
            else:
                if type == 'int':
                    return int(s)
                else:
                    return str(s)

        IDsdata = {}
        # sources_data = {}
        #create a dictionary of the complete info in the LigandIDs
        #sheet, using GtP - Ligand ID column as key and setting the
        #GtP - PubChem CID, GtP - Ligand Name, GtP - SMILES, Accession Number
        #total 5 values
        for i in range(1, IDs_rows): #skip the header row
            if IDs.cell_value(i, 4) not in IDsdata.keys():
                IDsdata[IDs.cell_value(i, 4)] = [IDs.cell_value(i, 3), parseValue(IDs.cell_value(i, 7)), parseValue(IDs.cell_value(i, 11), 'str'), IDs.cell_value(i, 1)]

        #Now need to parse the Bouvier and Inoue ligand data
        #search for matches and provide expanded info for get_or_make_ligand function
        for triple in tupled_sources:
            for i in range(1, triple[1]):
                if triple[0].cell_value(i, 5) not in IDsdata.keys():
                    #if the ligand is not in the IDs sheet, save the UniProt and the ligand name
                    IDsdata[triple[0].cell_value(i, 5)] = [triple[0].cell_value(i, 0), triple[0].cell_value(i, triple[2])]

        #this creates a dictionary with all the ligands data from the sources
        #with PubChem CID, SMILES and Ligand Name where available.
        #Otherwise we have the ligand name from the sheet
        return IDsdata


    def assess_variants(self, dictionary, iterator, column_variants):
        """
        Function that assesses the variants of the G-Protein coupling or Arrestin coupling.
        For now variants assessed are: Isoform 2 of Go and Arrestin without GRK, but more can be further
        implemented, thus we need to keep this function highly flexible.
        """
        if iterator in column_variants:
            dictionary['variant'] = column_variants[iterator]
        else:
            dictionary['variant'] = 'Regular'

    def assess_type(self, accession_id):
        """
        Function that assess the type of the ligand tested making a call to the Protein model.
        It is used to provide precise info to the get_or_make_ligand function in terms of
        ligand_type for the ligand. Usually used for Peptide or Protein calls. Defaults as Peptide call.
        """
        call = list(Protein.objects.filter(accession=accession_id).values_list("family__parent__parent__name"))
        label = call[0][0].split(' ')[0].lower()
        if (label != 'peptide') and (label != 'protein'):
            label = 'peptide'
        return label

# this is the new function that will overwrite the function read_coupling
# reading all different sheets of the single file
    def read_all_coupling(self, filenames=False):
        """
        Yet another function to read G-protein coupling data coming in Excel files.
        The idea is that now the format will hopefully be fixed in the same way for data
        coming from different groups. For now the data comes from Bouvier, Inoue and Roth
        and has been processed by David Gloriam.
        """
        book = xlrd.open_workbook(filenames)

        bouvier = book.sheet_by_name("B-import")
        inoue = book.sheet_by_name("I-import")
        roth = book.sheet_by_name("R-import")
        B_rows = bouvier.nrows
        I_rows = inoue.nrows
        R_rows = roth.nrows

        combinations = [(bouvier, B_rows, 'Bouvier'), (inoue, I_rows, 'Inoue'), (roth, R_rows, 'Roth')]
        #variable that sets the number of subunits parsed
        #in the spreadsheet. If more subunits are introduced
        #please update this number
        subunits = 17
        #each key if the section of data in the spreadsheet
        #values are: start column, ArrB2_NO_GRK column, GoA column and GoB column
        variants = ["isoform 1 (GoA)", "isoform 2 (GoB)", "no GRK"]
        variant_indices = [31, 32, 40]
        start_index = 26
        columns = {
                    'logemaxec50':   {"start": start_index, "variants": dict(zip(variant_indices, variants))},
                    'pec50deg':     {"start": start_index+subunits*1, "variants": dict(zip([i+subunits*1 for i in variant_indices], variants))},
                    'emaxdeg':      {"start": start_index+subunits*2, "variants": dict(zip([i+subunits*2 for i in variant_indices], variants))},
                    'stddeg':       {"start": start_index+subunits*3, "variants": dict(zip([i+subunits*3 for i in variant_indices], variants))}
                    }

        data = {}
        """data is a dictionary and must have a format:
        {'<protein>':
            {'<gproteinsubunit>':
                {'logemaxec50': <logemaxec50>,
                 'pec50deg': <pec50deg>,
                 'emaxdeg': <emaxdeg>,
                 'stddeg': <stddeg>}
            }
        }"""

        ligands = self.assess_ligand_id(book)

        #for each sheet of data from souces
        for tuple in combinations:
            data[tuple[2]] = {}
            #for each row, set protein as dict key
            #and fetch fixed columns info
            for i in range(2, tuple[1]):
                protein = tuple[0].cell_value(i, 0)
                protein_dict = {}
                protein_dict['ligand_id'] = tuple[0].cell_value(i, 4)
                if tuple[0].cell_value(i, 5) == 'Physiological':
                    protein_dict['ligand_physiological'] = True
                else:
                    protein_dict['ligand_physiological'] = False
                #for each block of data parse the sheet
                #and retrieve the associated info
                for key in columns.keys():
                    for j in range(columns[key]["start"], (columns[key]["start"] + subunits)):
                        gproteinsubunit = tuple[0].cell_value(1, j).split("\n")[-1]
                        if gproteinsubunit not in protein_dict.keys():
                            protein_dict[gproteinsubunit] = {}
                        textType = tuple[0].cell(i, j).ctype
                        if textType == 5:
                            protein_dict[gproteinsubunit][key] = None
                        else:
                            protein_dict[gproteinsubunit][key] = tuple[0].cell_value(i, j)
                        self.assess_variants(protein_dict[gproteinsubunit], j, columns[key]["variants"])
                #apply temporary dict to master dict
                data[tuple[2]][protein] = protein_dict

        return data, ligands

    def read_new_coupling(self, filepath):
        pe = ParseExcel()
        data_raw = pe.parse_excel(filepath, True)

        ligands = self.assess_ligand_id(pe.workbook)

        signprots = Protein.objects.filter(Q(family__slug__istartswith='100_001') | Q(family__slug__istartswith='200'), species__common_name='Human', accession__isnull=False).prefetch_related('family', 'family__parent')

        ### Predefined subtype columns where variants (isoform, GRK) are considered separate entries
        ###### Needs updating if new subtypes, variants (column headers) are added to master excel
        signprot_dict = {### G proteins
                         'GsL': signprots.get(entry_name='gnas2_human'), 'Gi1': signprots.get(entry_name='gnai1_human'), 'Gi2': signprots.get(entry_name='gnai2_human'), 
                         'Gi3': signprots.get(entry_name='gnai3_human'),
                         'GoA': signprots.get(entry_name='gnao_human'), 'GoB': signprots.get(entry_name='gnao_human'), 'Gz': signprots.get(entry_name='gnaz_human'),
                         'Gq': signprots.get(entry_name='gnaq_human'), 'G11': signprots.get(entry_name='gna11_human'), 'G14': signprots.get(entry_name='gna14_human'), 
                         'G15': signprots.get(entry_name='gna15_human'), 'G12': signprots.get(entry_name='gna12_human'), 'G13': signprots.get(entry_name='gna13_human'),
                         'GsS': signprots.get(entry_name='gnas2_human'), 'Golf': signprots.get(entry_name='gnal_human'), 'Ggust': signprots.get(entry_name='gnat3_human'),
                         ### Arrestins
                         'β1/ GRK2': signprots.get(entry_name='arrb1_human'), 'β2': signprots.get(entry_name='arrb2_human'), 'β2/ GRK2': signprots.get(entry_name='arrb2_human')}

        data = {}
        for source, source_metadata in data_raw['Metadata'].items():
            data[source] = {'metadata':source_metadata}
            for receptor, row in data_raw[source].items():
                if 'Receptor (UniProt)' in row:
                    receptor = row['Receptor (UniProt)']
                elif 'UniProt' in row:
                    receptor = row['UniProt']
                # receptor_variant = None
                ### Variant info - check for receptor variants, e.g. with RAMP
                if 'RAMP' in receptor:
                    # rec, receptor_variant = receptor.split('+')
                    receptor = receptor.replace('+','/')
                data[source][receptor] = {}
                ### Loop through
                for sp_key, sp in signprot_dict.items():
                    ### Skip is G prot subtype not in this source
                    if sp_key not in row:
                        continue
                    data[source][receptor][sp_key] = {}
                    ### Family data
                    if 'Rank {} fam'.format(sp.family.parent.name) in row:
                        data[source][receptor][sp_key]['fam_rank'] = int(row['Rank {} fam'.format(sp.family.parent.name)])
                    if "% of 1' family {}Fam".format(sp.family.parent.name) in row:
                        data[source][receptor][sp_key]['fam_percent'] = round(row["% of 1' family {}Fam".format(sp.family.parent.name)])
                    ###### Family values
                    if '{}-fam'.format(sp.family.parent.name) in row and source_metadata['Parameter']=='log(Emax/EC50)':
                        data[source][receptor][sp_key]['fam_logemaxec50'] = round(row['{}-fam'.format(sp.family.parent.name)],1)

                        # data[source][receptor][sp_key]['fam_emax'] = round(row[f'Emax {sp.family.parent.name}'])
                    if '{}-fam'.format(sp.family.parent.name) in row and source_metadata['Parameter']=='Activation rate (s-1)':
                        data[source][receptor][sp_key]['fam_kon_mean'] = round(row['{}-fam'.format(sp.family.parent.name)],1)
                    if '{}-fam'.format(sp.family.parent.name) in row and source_metadata['Parameter'] in ['Efficacy','Econstitutive']:
                        data[source][receptor][sp_key]['fam_deltaGDP_conc'] = round(row['{}-fam'.format(sp.family.parent.name)],2)
                    ### Subtype data
                    if "% of 1' subtype "+sp_key in row:
                        if row["% of 1' subtype "+sp_key]=='':
                            data[source][receptor][sp_key]['subtype_percent'] = None
                        else:
                            data[source][receptor][sp_key]['subtype_percent'] = round(row["% of 1' subtype "+sp_key])
                    ###### Subtype values
                    if sp_key in row and source_metadata['Parameter']=='log(Emax/EC50)':
                        data[source][receptor][sp_key]['logemaxec50'] = round(row[sp_key],1)

                        try:
                            data[source][receptor][sp_key]['fam_emax'] = round(float(row[f'Emax {sp_key.replace(" ", "")}']), 1)
                        except:
                            new_a = row[f'Emax {sp_key.replace(" ", "")}']
                            print(f'ERROR new_a: {new_a}, {sp.family.parent.name}, sp_key: {sp_key}')
                            data[source][receptor][sp_key]['fam_emax'] = 0
                        try:
                            data[source][receptor][sp_key]['fam_pec50'] = round(float(row[f'pEC50 {sp_key.replace(" ","")}']), 1)
                        except:
                            new_b = row[f'pEC50 {sp_key.replace(" ", "")}']
                            print(f'ERROR new_b: {new_b}, {sp.family.parent.name}, spkey: {sp_key}')
                            data[source][receptor][sp_key]['fam_pec50'] = 0
                        a = row[f'Emax {sp_key.replace(" ", "")}']
                        print(f'Emax_score: {sp_key}')

                    if sp_key in row and source_metadata['Parameter']=='Activation rate (s-1)':
                        data[source][receptor][sp_key]['kon_mean'] = round(row[sp_key],1)

                    if sp_key in row and source_metadata['Parameter'] in ['Efficacy','Econstitutive']:
                        if row[sp_key]=='':
                            data[source][receptor][sp_key]['deltaGDP_conc'] = None
                        else:
                            data[source][receptor][sp_key]['deltaGDP_conc'] = round(row[sp_key], 2)
                                         
                    ### Additional data
                    ### Ligand info
                    if 'Ligand ID' in row and 'ligand_id' not in data[source][receptor]:
                        data[source][receptor]['ligand_id'] = row['Ligand ID']
                    if 'Phys/ Surr' in row and 'ligand_physiological' not in data[source][receptor]:
                        if row['Phys/ Surr'] in ['Phys','Physiological']:
                            data[source][receptor]['ligand_physiological'] = True
                        else:
                            data[source][receptor]['ligand_physiological'] = False
                    if 'Ligand ID type' in row and row['Ligand ID type']=='_FASTA seq (peptide)':
                        del data[source][receptor]['ligand_id']
                    ### DOI
                    if 'DOI' in row:
                        data[source][receptor]['doi'] = row['DOI']
                    ### Variant info - separated by | when multi-variant
                    variant = None
                    if sp_key=='GoA':
                        variant = self.add_divider(variant, 'isoform 1 (GoA)')
                    elif sp_key=='GoB':
                        variant = self.add_divider(variant, 'isoform 2 (GoB)')
                    elif 'GRK' in sp_key:
                        variant = self.add_divider(variant, '+ GRK')
                    elif sp_key=='GsL':
                        variant = self.add_divider(variant, 'isoform 1 (GsL)')
                    elif sp_key=='GsS':
                        variant = self.add_divider(variant, 'isoform 2 (GsS)')
                    if variant:
                        data[source][receptor][sp_key]['variant'] = variant
                    else:
                        data[source][receptor][sp_key]['variant'] = 'Regular'
                    data[source][receptor][sp_key]['prot_obj'] = sp

        return data, ligands

    @staticmethod
    def add_divider(prev_variant, variant):
        if prev_variant:
            return prev_variant+' | {}'.format(variant)
        else:
            return variant

    def create_data_couplings(self):
        """This function adds all coupling data coming from the master Excel file."""
        self.logger.info('CREATE data couplings')
        print('CREATE data couplings')
        # read source files
        filepath = self.master_file
        self.logger.info('Reading file ' + filepath)
        #new function for reading data from the master spreadsheet
        data, ligands = self.read_new_coupling(filepath)
        lookup = {}
        lig_lookup = {}
        bulk = []
        publications = {}
        ######### MODIFIED #########

        for source, data_vals in data.items():
            # print("PROCESSING: "+str(source)+" DATA")
            ### Split lab and biosensor
            source = source.split('-')[0]
            # if source=='Lambert':
            #     import pprint
            #     pprint.pprint(data_vals)
            ### Build biosensor
            biosensor, _ = Biosensor.objects.get_or_create(name=data_vals['metadata']['Biosensor'], downstream_steps=int(data_vals['metadata']['Down- stream steps']), parameter=data_vals['metadata']['Parameter'])
            for entry_name, couplings in data_vals.items():
                if entry_name=='metadata':
                    continue
                # if it has / then pick first, since it gets same protein
                entry_split = entry_name.split("/")
                entry_name = entry_split[0]
                if len(entry_split)>1:
                    other_protein = entry_split[1]
                else:
                    other_protein = None
                # print("PROCESSING: "+str(entry_name)+" ROW")
                # Fetch protein
                try:
                    p = Protein.objects.get(entry_name=entry_name.lower()+"_human")
                except Protein.DoesNotExist:
                    self.logger.warning('Protein not found for entry_name {}'.format(entry_name))
                    print("protein not found for ", entry_name)
                    continue

                gproteins = list(couplings.keys())

                for header in gproteins:
                    # print("PROCESSING: "+str(header)+" COLUMN")
                    if isinstance(couplings[header], dict):
                        gp = couplings[header]['prot_obj']
                    else:
                        continue

                    # Assume they are there.
                    if gp.family.slug not in lookup:
                        # print("SEARCHING FOR FAMILY OF SLUG: " + str("_".join(gp.family.slug.split("_")[:3])))
                        lookup[gp.family.slug] = ProteinFamily.objects.get(slug="_".join(gp.family.slug.split("_")[:3]))

                    g = lookup[gp.family.slug]
                    # print("PROCESSING: LIGANDS INFORMATION OF " + str(ligands[couplings['ligand_id']][0]) + " FOR "+ str(g))
                    #ligand here should be fetched via the function get_or_make_ligand
                    #(0) Ligand Name, (1) PubChem CID, (2) SMILES, (3) Accession Number
                    if 'ligand_physiological' not in couplings:
                        lig_phys_surr = None
                    else:
                        lig_phys_surr = couplings['ligand_physiological']
                    if 'ligand_id' not in couplings:
                        l = None
                    else:
                        if couplings['ligand_id'] in lig_lookup:
                            l = lig_lookup[couplings['ligand_id']]
                        elif len(ligands[couplings['ligand_id']]) == 4:
                            if ligands[couplings['ligand_id']][1] != 'NA':
                                ids = {"pubchem" : ligands[couplings['ligand_id']][1]}
                                l = get_or_create_ligand(ligands[couplings['ligand_id']][0], ids)
                            elif ligands[couplings['ligand_id']][2] != 'NA':
                                ids = {"smiles" : ligands[couplings['ligand_id']][2]}
                                l = get_or_create_ligand(ligands[couplings['ligand_id']][0], ids)
                            else:
                                #make the call to search for peptide/protein
                                # print("NO INFO IN THE DATABASE. ADDING NEW LIGAND.")
                                label = self.assess_type(ligands[couplings['ligand_id']][3])
                                ids = {}
                                if ligands[couplings['ligand_id']][3] != 'NA':
                                    ids = {"uniprot": ligands[couplings['ligand_id']][3]}

                                l = get_or_create_ligand(ligands[couplings['ligand_id']][0], ids, label)
                        else:
                            #make the call to search for peptide/protein
                            # print("NO INFO IN THE DATABASE. ADDING NEW LIGAND.")
                            label = self.assess_type(ligands[couplings['ligand_id']][1])
                            l = get_or_create_ligand(ligands[couplings['ligand_id']][0], {}, label)

                        lig_lookup[couplings['ligand_id']] = l
                    #here it fills the data in ProteinCouplings model
                    #new data in dictionary: {'ligand_name': 'Serotonin', 'ligand_id': 5.0, 'ligand_type': 'Physiological'}
                    if couplings[header]:
                        ### None is set to fields where there are no values in the data
                        ###### Needs updating if new fields are addedK
                        vals = {k: None for k in ['variant','fam_emax', 'fam_pec50','logemaxec50','subtype_percent','kon_mean','deltaGDP_conc','fam_rank','fam_percent','fam_logemaxec50','fam_kon_mean','fam_deltaGDP_conc']}
                        for key in vals:
                            if key in couplings[header]:
                                vals[key] = couplings[header][key]
                                print(f'HEADERS: {couplings[header]}')
                                print(f'Couplings: {key}: {couplings[header][key]}')
                        # print(p,g,gp,l,source,vals)
                        gpair = ProteinCouplings(protein=p, g_protein=g, g_protein_subunit=gp, ligand=l, source=source, physiological_ligand=lig_phys_surr,
                                                 variant=vals['variant'],
                                                 logemaxec50=vals['logemaxec50'],
                                                 percent_of_primary_subtype=vals['subtype_percent'],
                                                 kon_mean=vals['kon_mean'],
                                                 deltaGDP_conc=vals['deltaGDP_conc'],
                                                 family_rank=vals['fam_rank'],
                                                 percent_of_primary_family=vals['fam_percent'],
                                                 logemaxec50_family=vals['fam_logemaxec50'],
                                                 pec50=vals['fam_pec50'],
                                                 emax=vals['fam_emax'],
                                                 kon_mean_family=vals['fam_kon_mean'],
                                                 deltaGDP_conc_family=vals['fam_deltaGDP_conc'],
                                                 other_protein=other_protein,
                                                 biosensor=biosensor)

                        try:
                            pub = Publication.objects.get(web_link__index=couplings['doi'])
                        except Publication.DoesNotExist:
                            pub = Publication()
                            pub.get_or_create_from_doi(couplings['doi'])
                        if l:
                            lig = l.pk
                        else:
                            lig = None
                        publications['{}|{}|{}|{}|{}'.format(source, biosensor.name, p.entry_name, lig, other_protein)] = pub
                                                 
                        # print("APPENDING CALL TO BULK")
                        bulk.append(gpair)

        ProteinCouplings.objects.bulk_create(bulk)

        ### Add publication many to many after ProteinCouplings bulk create
        for key, pub in publications.items():
            s, b, p, l, o = key.split('|')
            if l=='None':
                l = None
            if o=='None':
                o = None
            
            couplings = ProteinCouplings.objects.filter(source=s, biosensor__name=b, protein__entry_name=p, ligand__pk=l, other_protein=o)
            for c in couplings:
                c.references.add(pub)
                c.save()

        self.logger.info('COMPLETED CREATE data couplings')

    def create_data_couplings_depracated(self):
        """This function adds all coupling data coming from the master Excel file."""
        self.logger.info('CREATE data couplings')

        # read source files
        filepath = self.master_file
        self.logger.info('Reading file ' + filepath)
        #new function for reading data from the master spreadsheet
        data, ligands = self.read_all_coupling(filepath)
        lookup = {}
        lig_lookup = {}
        bulk = []
        ######### MODIFIED #########

        for source in data.keys():
            print("PROCESSING: "+str(source)+" DATA")
            for entry_name, couplings in data[source].items():
                # if it has / then pick first, since it gets same protein
                entry_name = entry_name.split("/")[0]
                # print("PROCESSING: "+str(entry_name)+" ROW")
                # Fetch protein
                try:
                    p = Protein.objects.get(entry_name=entry_name.lower()+"_human")
                except Protein.DoesNotExist:
                    self.logger.warning('Protein not found for entry_name {}'.format(entry_name))
                    print("protein not found for ", entry_name)
                    continue

                gproteins = list(couplings.keys())[2:]

                for header in gproteins:
                    # print("PROCESSING: "+str(header)+" COLUMN")
                    gprotein = header.split('/')[0]
                    if gprotein not in lookup:
                        gp = Protein.objects.filter(family__name=gprotein, species__common_name="Human")[0]
                        lookup[gprotein] = gp
                    else:
                        gp = lookup[gprotein]

                    # Assume they are there.
                    if gp.family.slug not in lookup:
                        # print("SEARCHING FOR FAMILY OF SLUG: " + str("_".join(gp.family.slug.split("_")[:3])))
                        lookup[gp.family.slug] = ProteinFamily.objects.get(slug="_".join(gp.family.slug.split("_")[:3]))

                    g = lookup[gp.family.slug]
                    # print("PROCESSING: LIGANDS INFORMATION OF " + str(ligands[couplings['ligand_id']][0]) + " FOR "+ str(g))
                    #ligand here should be fetched via the function get_or_make_ligand
                    #(0) Ligand Name, (1) PubChem CID, (2) SMILES, (3) Accession Number
                    if couplings['ligand_id'] in lig_lookup:
                        l = lig_lookup[couplings['ligand_id']]
                    elif len(ligands[couplings['ligand_id']]) == 4:
                        if ligands[couplings['ligand_id']][1] != 'NA':
                            ids = {"pubchem" : ligands[couplings['ligand_id']][1]}
                            l = get_or_create_ligand(ligands[couplings['ligand_id']][0], ids)
                        elif ligands[couplings['ligand_id']][2] != 'NA':
                            ids = {"smiles" : ligands[couplings['ligand_id']][2]}
                            l = get_or_create_ligand(ligands[couplings['ligand_id']][0], ids)
                        else:
                            #make the call to search for peptide/protein
                            # print("NO INFO IN THE DATABASE. ADDING NEW LIGAND.")
                            label = self.assess_type(ligands[couplings['ligand_id']][3])
                            ids = {}
                            if ligands[couplings['ligand_id']][3] != 'NA':
                                ids = {"uniprot": ligands[couplings['ligand_id']][3]}

                            l = get_or_create_ligand(ligands[couplings['ligand_id']][0], ids, label)
                    else:
                        #make the call to search for peptide/protein
                        # print("NO INFO IN THE DATABASE. ADDING NEW LIGAND.")
                        label = self.assess_type(ligands[couplings['ligand_id']][1])
                        l = get_or_create_ligand(ligands[couplings['ligand_id']][0], {}, label)

                    lig_lookup[couplings['ligand_id']] = l
                    #here it fills the data in ProteinCouplings model
                    #new data in dictionary: {'ligand_name': 'Serotonin', 'ligand_id': 5.0, 'ligand_type': 'Physiological'}
                    # print("CREATING CALL TO MODEL")
                    if (couplings[header]['logemaxec50'] != None) or (couplings[header]['pec50deg'] != None) or (couplings[header]['emaxdeg'] != None) or (couplings[header]['stddeg'] != None):
                        gpair = ProteinCouplings(protein=p,
                                                   g_protein=g,
                                                   ligand=l,
                                                   variant=couplings[header]['variant'],
                                                   source=source,
                                                   logemaxec50=couplings[header]['logemaxec50'],
                                                   pec50=couplings[header]['pec50'],
                                                   emax=couplings[header]['emax'],
                                                   stand_dev=couplings[header]['stddeg'],
                                                   physiological_ligand=couplings['ligand_physiological'],
                                                   g_protein_subunit=gp)
                        # print("APPENDING CALL TO BULK")
                        bulk.append(gpair)

        ProteinCouplings.objects.bulk_create(bulk)

        self.logger.info('COMPLETED CREATE data couplings')
