from Bio.Blast import NCBIXML, NCBIWWW
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.PDBIO import Select
import Bio.PDB.Polypeptide as polypeptide
from Bio.PDB.AbstractPropertyMap import AbstractPropertyMap
from Bio.PDB.Polypeptide import CaPPBuilder, is_aa
try:
    from Bio.PDB.vectors import rotaxis
except:
    from Bio.PDB import rotaxis
from Bio.Align import PairwiseAligner

from django.conf import settings
from common.alignment import Alignment
from common.tools import urlopen_with_retry
from common.models import WebResource, WebLink, Publication
from protein.models import Protein, ProteinSegment, ProteinConformation, ProteinState
from residue.functions import dgn
from residue.models import Residue, ResidueGenericNumberEquivalent
from structure.models import Structure, Rotamer, PdbData, StructureStabilizingAgent, StructureType
from signprot.models import SignprotStructure
from ligand.models import Endogenous_GTP

from subprocess import Popen, PIPE
from io import StringIO
import pandas as pd
import os
import sys
import tempfile
import logging
import math
import urllib
from collections import OrderedDict
import Bio.PDB as PDB
from Bio import SeqIO
import csv
import numpy
import json
import yaml
from urllib.request import urlopen, Request
from datetime import datetime
import re


logger = logging.getLogger("protwis")

ATOM_FORMAT_STRING="%s%5i %-4s%c%3s %c%4i%c   %8.3f%8.3f%8.3f%s%6.2f      %4s%2s%2s\n"

#==============================================================================
# I have put it into separate class for the sake of future uses
class BlastSearch(object):


    def __init__ (self, blast_path='blastp',
        blastdb=os.sep.join([settings.STATICFILES_DIRS[0], 'blast', 'protwis_blastdb']), top_results=1):

        self.blast_path = blast_path
        self.blastdb = blastdb
        #typicaly top scored result is enough, but for sequences with missing
        #residues it is better to use more results to avoid getting sequence of
        #e.g.  different species
        self.top_results = top_results

    #takes Bio.Seq sequence as an input and returns a list of tuples with the
    #alignments
    def run (self, input_seq):

        output = []
        #Windows has problems with Popen and PIPE
        if sys.platform == 'win32':
            tmp = tempfile.NamedTemporaryFile()
            logger.debug("Running Blast with sequence: {}".format(input_seq))
            tmp.write(bytes(str(input_seq) + '\n', 'latin1'))
            tmp.seek(0)
            blast = Popen('%s -db %s -outfmt 5' % (self.blast_path, self.blastdb), universal_newlines=True, stdin=tmp,
                stdout=PIPE, stderr=PIPE)
            (blast_out, blast_err) = blast.communicate()
        else:
            #Rest of the world:
            blast = Popen('%s -db %s -outfmt 5' % (self.blast_path, self.blastdb), universal_newlines=True, shell=True,
                stdin=PIPE, stdout=PIPE, stderr=PIPE)
            (blast_out, blast_err) = blast.communicate(input=str(input_seq))

            # Log the BLAST output and errors for inspection
            logger.debug("BLAST Output: {}".format(blast_out))
            if blast_err:
                logger.error("BLAST Error: {}".format(blast_err))

            # Check if BLAST output is empty
            if not blast_out.strip():
                logger.error("No output returned from BLAST command.")
                return []

        if len(blast_err) != 0:
            logger.debug(blast_err)
        if blast_out!='\n':
            result = NCBIXML.read(StringIO(blast_out))
            for aln in result.alignments[:self.top_results]:
                logger.debug("Looping over alignments, current hit: {}".format(aln.hit_id))
                output.append((aln.hit_id, aln))
        return output
#==============================================================================

class BlastSearchOnline(object):

    def __init__(self, blast_program='blastp', db='swissprot', top_results=1):
        self.blast_program = blast_program
        self.db = db
        self.top_results = top_results
        pass

    def run(self, input_seq):
        output = []

        result = NCBIXML.read(NCBIWWW.qblast(self.blast_program, self.db, input_seq, auto_format='xml'))
        for aln in result.alignments[:self.top_results]:
            logger.debug("Looping over alignments, current hit: {}".format(aln.hit_id))
            output.append((aln.hit_id, aln))
        return output

    def get_uniprot_entry_name (self, up_id):

        #get the 'entry name' field for given uniprot id
        #'entry name' is a protein id in gpcrdb
        url = "http://www.uniprot.org/uniprot/?query=accession:%s&columns=entry name&format=tab" %up_id

        #used urllib, urllib2 throws an error here for some reason
        try:
            response = urlopen(url) #nosec
            page = response.readlines()
            return page[1].strip().lower()

        except urllib.HTTPError as error:
            print(error)
            return ''

#==============================================================================

#stores information about alignments and b-w numbers
class MappedResidue(object):

    def __init__ (self, res_num, res_name):

        self.number = res_num
        self.name = res_name
        self.pos_in_aln = 0
        self.mapping = {}
        self.bw = 0.
        self.gpcrdb = 0.
        self.gpcrdb_id = 0
        self.segment = ''
        self.display = ''
        self.residue_record = None

    def add_bw_number (self, bw_number=''):

        self.bw = bw_number

    def add_segment (self, segment=''):

        self.segment = segment

    def add_display_number (self, display = ''):

        self.display = display

    def add_gpcrdb_number_id (self, gpcrdb_number_id=''):

        self.gpcrdb_id = gpcrdb_number_id

    def add_residue_record (self, residue_record= None):

        self.residue_record = residue_record

    def add_gpcrdb_number (self, gpcrdb_number=''):

        #PDB format does not allow fractional part longer than 2 digits
        #so numbers x.xx1 are negative
        if len(gpcrdb_number.split('.')[1]) > 2:
          self.gpcrdb = '-' + gpcrdb_number[:4].replace('x', '.')
        else:
          self.gpcrdb = gpcrdb_number.replace('x', '.')

#==============================================================================

#turns selection into actual residues
class SelectionParser(object):

    def __init__ (self, selection, full_helix_name=False):

        self.generic_numbers = []
        self.helices = []
        self.substructures = []

        for segment in selection.segments:
            logger.debug('Segments in selection: {}'.format(segment))
            if segment.type == 'helix':
                if full_helix_name:
                    self.helices.append(segment.item.slug)
                else:
                    try:
                        self.helices.append(int(segment.item.slug[-1]))
                    except ValueError:
                        self.helices.append(segment.item.slug)
            elif segment.type == 'residue':
                self.generic_numbers.append(segment.item.label.replace('x','.'))
            else:
                self.substructures.append(segment.item.slug)
        logger.debug("Helices selected: {}; Residues: {}; Other substructures:{}".format(self.helices, self.generic_numbers, self.substructures))


#==============================================================================
class GenericNumbersSelector(Select):

    def __init__(self, generic_numbers=[], helices=[], parsed_selection=None):

        self.generic_numbers = generic_numbers
        self.helices = helices
        if parsed_selection:
            self.generic_numbers=parsed_selection.generic_numbers
            self.helices = parsed_selection.helices

    def accept_residue(self, residue):

        try:
            if "{:.2f}".format(residue['CA'].get_bfactor()) in self.generic_numbers:
                return 1
            if -8.1 < residue['CA'].get_bfactor() < 0 and "{:.3f}".format(-residue['CA'].get_bfactor() + 0.001) in self.generic_numbers:
                return 1
            if -8.1 < residue['CA'].get_bfactor() < 8.1 and int(math.floor(abs(residue['CA'].get_bfactor()))) in self.helices:
                return 1
        except:
            return 0
        return 0

#==============================================================================
class SubstructureSelector(Select):

    def __init__(self, segment_mapping, parsed_selection=None):

        self.residues = []

        for tm in parsed_selection.helices:
            if tm==8:
                seg_str = 'H'
            else:
                seg_str = 'TM'
            segname = seg_str+str(tm)
            if segname in segment_mapping:
                self.residues.extend(segment_mapping[segname])
        for substr in parsed_selection.substructures:
            if substr in segment_mapping:
                self.residues.extend(segment_mapping[substr])

    def accept_residue(self, residue):

        try:
            if int(residue.id[1]) in self.residues:
                return 1
        except:
            return 0
        return 0

#==============================================================================
class CASelector(object):

    def __init__ (self, parsed_selection, ref_pdbio_struct, alt_structs):

        self.ref_atoms = []
        self.alt_atoms = {}

        self.selection = parsed_selection
        try:
            self.ref_atoms.extend(self.select_generic_numbers(ref_pdbio_struct))
            self.ref_atoms.extend(self.select_helices(ref_pdbio_struct))
        except Exception as msg:
            logger.warning("Can't select atoms from the reference structure!\n{!s}".format(msg))

        for alt_struct in alt_structs:
            try:
                self.alt_atoms[alt_struct.id] = []
                self.alt_atoms[alt_struct.id].extend(self.select_generic_numbers(alt_struct))
                self.alt_atoms[alt_struct.id].extend(self.select_helices(alt_struct))

            except Exception as msg:
                logger.warning("Can't select atoms from structure {!s}\n{!s}".format(alt_struct.id, msg))


    def select_generic_numbers (self, structure):

        if self.selection.generic_numbers == []:
            return []

        atom_list = []

        for chain in structure:
            for res in chain:
                try:
                    if "{:.2f}".format(res['CA'].get_bfactor()) in self.selection.generic_numbers:
                        atom_list.append(res['CA'])
                    if -8.1 < res['CA'].get_bfactor() < 0 and "{:.3f}".format(-res['CA'].get_bfactor() + 0.001) in self.selection.generic_numbers:
                        atom_list.append(res['CA'])
                except:
                    continue

        if atom_list == []:
            logger.warning("No atoms with given generic numbers {} for {!s}".format(self.selection.generic_numbers, structure.id))
        return atom_list


    def select_helices (self, structure):

        if self.selection.helices == []:
            return []

        atom_list = []
        for chain in structure:
            for res in chain:
                try:
                    if -8.1 < res['CA'].get_bfactor() < 8.1 and int(math.floor(abs(res['CA'].get_bfactor()))) in self.selection.helices:
                        atom_list.append(res['CA'])
                except Exception:
                    logger.warning("No CA for residues {!s} in {!s}".format(res.id, structure.id))
                    continue

        if atom_list == []:
            logger.warning("No atoms with given generic numbers for  {!s}".format(structure.id))

        return atom_list


    def get_consensus_atom_sets (self, alt_id):

        tmp_ref = []
        tmp_alt = []

        for ref_at in self.ref_atoms:
            for alt_at in self.alt_atoms[alt_id]:
                if ref_at.get_bfactor() == alt_at.get_bfactor():
                    tmp_ref.append(ref_at)
                    tmp_alt.append(alt_at)

        if len(tmp_ref) != len(tmp_alt):
            return ([], [])

        return (tmp_ref, tmp_alt)


    def get_consensus_gn_set (self):

        gn_list = []
        for alt_id in self.alt_atoms.keys():
            tmp_ref, tmp_alt = self.get_consensus_atom_sets(alt_id)

            for ref_ca in tmp_ref:
                for alt_ca in tmp_alt:
                    if ref_ca.get_bfactor() == alt_ca.get_bfactor():
                        if 0 < ref_ca.get_bfactor() < 8.1 and "{:.2f}".format(ref_ca.get_bfactor()) in self.selection.generic_numbers:
                            gn_list.append("{:.2f}".format(ref_ca.get_bfactor()))
                        if -8.1 < ref_ca.get_bfactor() < 0 and "{:.3f}".format(-ref_ca.get_bfactor() + 0.001) in self.selection.generic_numbers:
                            gn_list.append("{:.3f}".format(-ref_ca.get_bfactor() + 0.001))
                        if 0 < ref_ca.get_bfactor() < 8.1 and int(math.floor(abs(ref_ca.get_bfactor()))) in self.selection.helices:
                            gn_list.append("{:.2f}".format(ref_ca.get_bfactor()))
                        if -8.1 < ref_ca.get_bfactor() < 0 and int(math.floor(abs(ref_ca.get_bfactor()))) in self.selection.helices:
                            gn_list.append("{:.3f}".format(-ref_ca.get_bfactor() + 0.001))
        return gn_list


    def get_ref_atoms (self):

        return self.ref_atoms


    def get_alt_atoms (self, alt_id):

        try:
            return self.alt_atoms[alt_id]
        except:
            return []


    def get_alt_atoms_all (self):

        return self.alt_atoms


#==============================================================================
class BackboneSelector():
    """
    Selects backbone atoms from reference structure and rotamer for Superposition
    """

    similarity_dict = {
        "fragment_residue" : 0,
        "interaction_type" : 1,
        "target_residue" : 2
        }

    #similarity_rules = [[['H', 'F', 'Y', 'W'], ['AEF', 'AFF'], ['H', 'F', 'Y', 'W']],
    #    [['Y'], ['AFE'], ['F']],
    #    [['S', 'T'], ['HBA', 'HBD'], ['S', 'T']],]

    #interaction_type slugs changed along the way
    similarity_rules = [[['H', 'F', 'Y', 'W'], ['aro_ef_protein', 'aro_ff'], ['H', 'F', 'Y', 'W']],
        [['Y'], ['aro_fe_protein'], ['F']],
        [['S', 'T'], ['polar_acceptor_protein', 'polar_donor_protein'], ['S', 'T']],]

    def __init__(self, ref_pdbio_struct, fragment, use_similar=False):

        self.ref_atoms = []
        self.alt_atoms = []

        self.ref_atoms = self.select_ref_atoms(fragment, ref_pdbio_struct, use_similar)
        self.alt_atoms = self.select_alt_atoms(PDBParser(PERMISSIVE=True, QUIET=True).get_structure('ref', StringIO(str(fragment.rotamer.pdbdata)))[0])


    def select_ref_atoms(self, fragment, ref_pdbio_struct, use_similar=False):

        for chain in ref_pdbio_struct:
            for res in chain:
                try:
                    gn = self.get_generic_number(res)
                    if gn == fragment.rotamer.residue.display_generic_number.label:
                        if polypeptide.three_to_one(res.resname) == fragment.rotamer.residue.amino_acid:
                            return [res['CA'], res['N'], res['O']]
                        else:
                            if use_similar:
                                for rule in self.similarity_rules:
                                    if polypeptide.three_to_one(res.resname) in rule[self.similarity_dict["target_residue"]] and \
                                        fragment.rotamer.residue.amino_acid in rule[self.similarity_dict["target_residue"]] and \
                                        fragment.interaction_type.slug in rule[self.similarity_dict["interaction_type"]]:
                                        return [res['CA'], res['N'], res['O']]
                        # else:
                        #     if fragment.interaction_type.slug not in ['acc', 'hyd']:
                        #         return [res['CA'], res['N'], res['O']]
                except Exception as msg:
                    continue
        return []


    def select_alt_atoms(self, rotamer_pdbio_struct):

        for chain in rotamer_pdbio_struct:
            for res in chain:
                try:
                    return [res['CA'], res['N'], res['O']]
                except:
                    continue
        return []


    def get_generic_number(self, res):

        if 'CA' not in res:
            return 0.0
        if 0 < res['CA'].get_bfactor() < 8.1:
            return "{:.2f}x{!s}".format(res['N'].get_bfactor(), self._get_fraction_string(res['CA'].get_bfactor()))
        if -8.1 < res['CA'].get_bfactor() < 0:
            return "{:.2f}x{!s}".format(res['N'].get_bfactor(),  self._get_fraction_string(res['CA'].get_bfactor() - 0.001))
        return 0.0

    #TODO: Is this function really neccessary?
    def _get_fraction_string(self, number):

        if number > 0:
            return "{:.2f}".format(number).split('.')[1]
        else:
            return "{:.3f}".format(number).split('.')[1]


    def get_ref_atoms(self):

        return self.ref_atoms


    def get_alt_atoms(self):

        return self.alt_atoms


#==============================================================================
def check_gn (pdb_struct):

    for chain in pdb_struct:
        for residue in chain:
            try:
                if -8.1 < residue['CA'].get_bfactor() < 8.1:
                    return True
            except:
                continue
    return False


#==============================================================================
def get_segment_template (protein, segments=['TM1', 'TM2', 'TM3', 'TM4','TM5','TM6', 'TM7'], state=None):

    a = Alignment()
    a.load_reference_protein(protein)
    #You are so gonna love it...
    if state:
        a.load_proteins([x.protein_conformation.protein.parent for x in list(Structure.objects.order_by('protein_conformation__protein__parent','resolution').exclude(protein_conformation__protein=protein.id, protein_conformation__state=state, structure_type__slug__startswith='af-'))])
    else:
        a.load_proteins([x.protein_conformation.protein.parent for x in list(Structure.objects.order_by('protein_conformation__protein__parent','resolution').exclude(protein_conformation__protein=protein.id, structure_type__slug__startswith='af-'))])
    a.load_segments(ProteinSegment.objects.filter(slug__in=segments))
    a.build_alignment()
    a.calculate_similarity()

    return a.proteins[1]


#==============================================================================
def fetch_template_structure (template_protein):

    return Structure.objects.get(protein_conformation__protein__parent=template_protein.entry_name).exclude(structure_type__slug__startswith='af-')


#==============================================================================
def extract_pdb_data(residue):
    """Returns PDB string of a given residue"""
    pdb_string = ''
    hetfield, resseq, icode=residue.get_id()
    resname=residue.get_resname()
    segid=residue.get_segid()
    atom_number = 1
    for atom in residue:
        pdb_string += get_atom_line(atom, hetfield, segid, atom_number, resname, resseq, icode, residue.get_parent().get_id())
        atom_number += 1
    return pdb_string


#==============================================================================
# def convert_csv_to_xlsx(self, csv, separator):
#     wb = Workbook()
#     sheet = wb.active

#     CSV_SEPARATOR = separator

#     with open(csv) as f:
#         reader = csv.reader(f)
#         for r, row in enumerate(reader):
#             for c, col in enumerate(row):
#                 for idx, val in enumerate(col.split(CSV_SEPARATOR)):
#                     cell = sheet.cell(row=r+1, column=idx+1)
#                     cell.value = val

#     wb.save("my_file.xlsx")



#==============================================================================
def get_atom_line(atom, hetfield, segid, atom_number, resname, resseq, icode, chain_id, charge="  "):
    """Returns an ATOM PDB string."""
    if hetfield!=" ":
        record_type="HETATM"
    else:
        record_type="ATOM  "
    if atom.element:
        element = atom.element.strip().upper()
        element = element.rjust(2)
    else:
        element = "  "
    name=atom.get_fullname()
    altloc=atom.get_altloc()
    x, y, z=atom.get_coord()
    bfactor=atom.get_bfactor()
    occupancy=atom.get_occupancy()
    try:
        occupancy_str = "%6.2f" % occupancy
    except TypeError:
        if occupancy is None:
            occupancy_str = " " * 6
            import warnings
            from Bio import BiopythonWarning
            warnings.warn("Missing occupancy in atom %s written as blank" % repr(atom.get_full_id()), BiopythonWarning)
        else:
            raise TypeError("Invalid occupancy %r in atom %r" % (occupancy, atom.get_full_id()))
        pass
    args=(record_type, atom_number, name, altloc, resname, chain_id, resseq, icode, x, y, z, occupancy_str, bfactor, segid, element, charge)
    return ATOM_FORMAT_STRING % args

#==============================================================================
class HSExposureCB(AbstractPropertyMap):
    """
    Abstract class to calculate Half-Sphere Exposure (HSE). #GP: Biopython class repurposed
    The HSE can be calculated based on the CA-CB vector, or the pseudo CB-CA
    vector based on three consecutive CA atoms. This is done by two separate
    subclasses.
    """
    def __init__(self, model, radius, offset=0, hse_up_key='HSE_U', hse_down_key='HSE_D', angle_key=None, check_chain_breaks=False,
                 check_knots=False, receptor=None, signprot=None,  restrict_to_chain=[], check_hetatoms=False):
        """
        @param model: model
        @type model: L{Model}

        @param radius: HSE radius
        @type radius: float

        @param offset: number of flanking residues that are ignored in the calculation
        of the number of neighbors
        @type offset: int

        @param hse_up_key: key used to store HSEup in the entity.xtra attribute
        @type hse_up_key: string

        @param hse_down_key: key used to store HSEdown in the entity.xtra attribute
        @type hse_down_key: string

        @param angle_key: key used to store the angle between CA-CB and CA-pCB in
        the entity.xtra attribute
        @type angle_key: string
        """
        assert(offset>=0)
        # For PyMOL visualization
        self.ca_cb_list=[]
        ppb=CaPPBuilder()
        ppl=ppb.build_peptides(model)
        hse_map={}
        hse_list=[]
        hse_keys=[]
        ### GP
        if model.get_id()!=0:
            model = model[0]
        residues_in_pdb,residues_with_proper_CA=[],[]
        if check_chain_breaks is True:
            for chain in model:
                for res in chain:
                    if is_aa(res):
                        residues_in_pdb.append(res.get_id()[1])
        het_resis, het_resis_close, het_resis_clash = [], [], []
        for chain in model:
            for res in chain:
                if res.get_id()[0]!=' ':
                    het_resis.append(res)
        self.clash_pairs = []
        self.chain_breaks = []

        if check_knots:
            possible_knots = PossibleKnots(receptor, signprot)
            knot_resis = possible_knots.get_resnums()
            self.remodel_resis = {}
        if len(restrict_to_chain)>0:
            restricted_ppl = []
            for p in ppl:
                if p[0].get_parent().get_id() in restrict_to_chain:
                    restricted_ppl.append(p)
            ppl = restricted_ppl
        ###########
        for pp1 in ppl:
            for i in range(0, len(pp1)):
                residues_with_proper_CA.append(pp1[i].get_id()[1])
                if i==0:
                    r1=None
                else:
                    r1=pp1[i-1]
                r2=pp1[i]
                if i==len(pp1)-1:
                    r3=None
                else:
                    r3=pp1[i+1]
                # This method is provided by the subclasses to calculate HSE
                result=self._get_cb(r1, r2, r3)
                if result is None:
                    # Missing atoms, or i==0, or i==len(pp1)-1
                    continue
                pcb, angle=result
                hse_u=0
                hse_d=0
                ca2=r2['CA'].get_vector()
                residue_up=[]   ### GP
                residue_down=[] ### GP
                for pp2 in ppl:
                    for j in range(0, len(pp2)):
                        try:
                            if r2.get_id()[1]-1!=r1.get_id()[1] or r2.get_id()[1]+1!=r3.get_id()[1]:
                                pass
                            else:
                                raise Exception
                        except:
                            if pp1 is pp2 and abs(i-j)<=offset:
                                # neighboring residues in the chain are ignored
                                continue
                        ro=pp2[j]
                        if not is_aa(ro) or not ro.has_id('CA'):
                            continue
                        cao=ro['CA'].get_vector()
                        d=(cao-ca2)
                        if d.norm()<radius:
                            if d.angle(pcb)<(math.pi/2):
                                hse_u+=1
                                ### GP
                                # Puts residues' names in a list that were found in the upper half sphere
                                residue_up.append(ro)

                                ### end of GP code
                            else:
                                hse_d+=1
                                ### GP
                                # Puts residues' names in a list that were found in the lower half sphere
                                residue_down.append(ro)
                                ### end of GP code
                res_id=r2.get_id()
                chain_id=r2.get_parent().get_id()
                # Fill the 3 data structures
                hse_map[(chain_id, res_id)]=(hse_u, hse_d, angle)
                hse_list.append((r2, (residue_up, residue_down, hse_u, hse_d, angle)))
                ### GP residue_up and residue_down added to hse_list
                hse_keys.append((chain_id, res_id))
                # Add to xtra
                r2.xtra[hse_up_key]=hse_u
                r2.xtra[hse_down_key]=hse_d
                if angle_key:
                    r2.xtra[angle_key]=angle

                ### GP checking for knots
                if check_knots:
                    for knot in knot_resis:
                        if knot[0][1]==pp1[i].get_id()[1] and knot[0][0]==pp1[i].get_parent().get_id():
                            # print(pp1[i].get_parent().get_id(),pp1[i]) #print reference
                            for r in residue_up:
                                if r.get_parent().get_id()==knot[1][0] and r.get_id()[1] in knot[1][1]:
                                    # print('close: ', r.get_parent().get_id(),r) #print res within radius
                                    resi_range = [knot[1][1][0], knot[1][1][-1]]
                                    if knot[1][0] not in self.remodel_resis:
                                        self.remodel_resis[knot[1][0]] = [resi_range]
                                    else:
                                        if resi_range not in self.remodel_resis[knot[1][0]]:
                                            self.remodel_resis[knot[1][0]].append(resi_range)

                ### GP checking for atom clashes
                include_prev, include_next = False, False
                try:
                    if pp1[i].get_id()[1]-1!=pp1[i-1].get_id()[1]:
                        include_prev = True
                except:
                    include_prev = False
                try:
                    if pp1[i].get_id()[1]+1!=pp1[i+1].get_id()[1]:
                        include_next = True
                except:
                    include_next = False
                for atom in pp1[i]:
                    ref_vector = atom.get_vector()
                    for other_res in residue_up:
                        try:
                            if other_res==pp1[i-1] and include_prev is False:
                                continue
                            elif len(pp1)>=i+1 and other_res==pp1[i+1] and include_next is False:
                                continue
                            else:
                                raise Exception
                        except:
                            for other_atom in other_res:
                                other_vector = other_atom.get_vector()
                                d = other_vector-ref_vector
                                if d.norm()<1.5:
                                    if len(str(pp1[i]['CA'].get_bfactor()).split('.')[1])==1:
                                        clash_res1 = float(str(pp1[i]['CA'].get_bfactor())+'0')
                                    else:
                                        clash_res1 = pp1[i]['CA'].get_bfactor()
                                    if len(str(other_res['CA'].get_bfactor()).split('.')[1])==1:
                                        clash_res2 = float(str(other_res['CA'].get_bfactor())+'0')
                                    else:
                                        clash_res2 = other_res['CA'].get_bfactor()
                                    self.clash_pairs.append([(clash_res1, pp1[i].get_id()[1]), (clash_res2, other_res.get_id()[1])])
                    if check_hetatoms:
                        for het_res in het_resis:
                            for het_atom in het_res:
                                het_atom_vector = het_atom.get_vector()
                                d = het_atom_vector-ref_vector
                                if d.norm()<6:
                                    if d.norm()<1:
                                        het_resis_clash.append(het_res)
                                    het_resis_close.append(het_res)
        ### GP checking HETRESIS to remove if not interacting with AAs
        self.hetresis_to_remove = []
        for i in het_resis:
            if i not in self.hetresis_to_remove and i not in het_resis_close or i in het_resis_clash:
                self.hetresis_to_remove.append(i)
        # self.hetresis_to_remove = [i for i in het_resis if i not in het_resis_close or i in het_resis_clash]
        if check_chain_breaks:
            for r in residues_in_pdb:
                if r not in residues_with_proper_CA:
                    self.chain_breaks.append(r)


    def _get_cb(self, r1, r2, r3):
        """
        Method to calculate CB-CA vector.

        @param r1, r2, r3: three consecutive residues (only r2 is used)
        @type r1, r2, r3: L{Residue}
        """
        if r2.get_resname()=='GLY':
            return self._get_gly_cb_vector(r2), 0.0
        else:
            if r2.has_id('CB') and r2.has_id('CA'):
                vcb=r2['CB'].get_vector()
                vca=r2['CA'].get_vector()
                return (vcb-vca), 0.0
        return None

    def _get_gly_cb_vector(self, residue):
        """
        Return a pseudo CB vector for a Gly residue.
        The pseudoCB vector is centered at the origin.

        CB coord=N coord rotated over -120 degrees
        along the CA-C axis.
        """
        try:
            n_v=residue["N"].get_vector()
            c_v=residue["C"].get_vector()
            ca_v=residue["CA"].get_vector()
        except:
            return None
        # center at origin
        n_v=n_v-ca_v
        c_v=c_v-ca_v
        # rotation around c-ca over -120 deg
        rot=rotaxis(-math.pi*120.0/180.0, c_v)
        cb_at_origin_v=n_v.left_multiply(rot)
        # move back to ca position
        cb_v=cb_at_origin_v+ca_v
        # This is for PyMol visualization
        self.ca_cb_list.append((ca_v, cb_v))
        return cb_at_origin_v


class PossibleKnots():
    def __init__(self, receptor, signprot):
        self.receptor = Protein.objects.get(entry_name=receptor)
        self.signprot = Protein.objects.get(entry_name=signprot)
        self.possible_knots = {'ICL3-H4':[['R','A'],['G.H4.11','G.H4.14','G.H4.15','G.h4s6.01']],
                               'G.h1ha-H.hehf':[['A','A'],['H.HE.08', 'H.hdhe.05']]}
        self.output = []

    def get_resnums(self):
        if self.receptor and self.signprot:
            for knot_label, values in self.possible_knots.items():
                chain1, chain2 = values[0]
                if not self.receptor.accession:
                    self.receptor = self.receptor.parent
                region1 = list(Residue.objects.filter(protein_conformation__protein=self.receptor, protein_segment__slug=knot_label.split('-')[0]).values_list('sequence_number', flat=True))
                if len(region1)==0:
                    region1 = list(Residue.objects.filter(protein_conformation__protein=self.signprot, protein_segment__slug=knot_label.split('-')[0]).values_list('sequence_number', flat=True))
                if len(region1)==0:
                    raise AssertionError('Protein segment slug error for loop knot: No residues found for {} in {} - {}'.format(knot_label.split('-')[0], self.receptor, self.signprot))
                if knot_label=='G.h1ha-H.hehf':
                    for i in range(0,3):
                        region1.append(region1[-1]+1)
                for r in values[1]:
                    region2 = Residue.objects.get(protein_conformation__protein=self.signprot, display_generic_number__label=r)
                    self.output.append([[chain2,region2.sequence_number],[chain1,region1]])
        return self.output


class PdbChainSelector():
    def __init__(self, pdb_code, protein):
        self.pdb_code = pdb_code
        try:
            protein.entry_name
        except:
            protein = Protein.objects.get(entry_name=protein)
        self.protein = protein
        self.chains = []
        self.dssp_dict = OrderedDict()
        self.dssp_info = OrderedDict()
        self.aux_residues = []

    def run_dssp(self):
        pdb = PDB.PDBList()
        pdb.retrieve_pdb_file(self.pdb_code, pdir='./', file_format="pdb")
        p = PDB.PDBParser()
        f = 'pdb{}.ent'.format(self.pdb_code.lower())
        wt_residues = [i for i in Residue.objects.filter(protein_conformation__protein=self.protein).exclude(protein_segment__slug__in=['N-term','C-term'])]
        gn_residues = [i.sequence_number for i in wt_residues if i.generic_number and i.protein_segment.slug not in ['ECL1','ECL2','ICL3','ECL3']]
        structure = p.get_structure(self.pdb_code, f)
        for chain in structure[0]:
            ch = chain.get_id()
            self.chains.append(ch)
            self.dssp_dict[ch] = OrderedDict()
            self.dssp_info[ch] = OrderedDict([('H',0),('B',0),('E',0),('G',0),('I',0),('T',0),('S',0),('-',0)])
        if len(self.dssp_dict)>1:
            dssp = PDB.DSSP(structure[0], f, dssp='/env/bin/dssp')
            for key in dssp.keys():
                if int(key[1][1]) in gn_residues:
                    self.dssp_dict[key[0]][key[1][1]] = dssp[key]
                    self.dssp_info[key[0]][dssp[key][2]] = self.dssp_info[key[0]][dssp[key][2]]+1
        os.remove(f)

    def get_seqnums_by_secondary_structure(self, chain, secondary_structure):
        ''' Returns list of sequence numbers that match the secondary structural property. \n
            H: alpha helix, B: isolated beta-bridge, E: strand, G: 3-10 helix, I: pi-helix, T: turn, S: bend, -: Other
        '''
        if secondary_structure not in ['H','B','E','G','I','T','S','-']:
            raise ValueError('Incorrect secondary_structure input. Input can be H, B, E, G, I, T, S, -')
        output = []
        for seqnum, val in self.dssp_dict[chain].items():
            if val[2]==secondary_structure:
                output.append(seqnum)
        return output

    def select_chain(self):
        num_helix_res = []
        seq_lengths = []
        for c, val in self.dssp_info.items():
            seq_length = 0
            num_helix_res.append(val['H']+val['G']+val['I'])
            for s, num in val.items():
                seq_length+=num
            seq_lengths.append(seq_length)

        max_res = num_helix_res[0]
        max_i = 0
        for i in range(1,len(num_helix_res)):
            if num_helix_res[i]>max_res:
                if num_helix_res[i]-max_res>=seq_lengths[max_i]-seq_lengths[i]:
                    max_res = num_helix_res[i]
                    max_i = i
            elif num_helix_res[i]<max_res:
                if max_res-num_helix_res[i]<seq_lengths[i]-seq_lengths[max_i]:
                    max_res = num_helix_res[i]
                    max_i = i
            elif num_helix_res[i]==max_res:
                if seq_lengths[max_i]<seq_lengths[i]:
                    max_res = num_helix_res[i]
                    max_i = i
        return self.chains[max_i]


class PdbStateIdentifier():
    def __init__(self, structure, tm2_gn='2x41', tm6_gn='6x38', tm3_gn='3x44', tm7_gn='7x52', inactive_cutoff=2, intermediate_cutoff=7.15):
        self.structure_type = None

        try:
            if structure.protein_conformation.protein.parent is None:
                raise Exception
            self.structure = structure
            self.structure_type = 'structure'
            family = structure.protein_conformation.protein.family
        except:
            try:
                structure.protein_conformation.protein
                self.structure = structure
                self.structure_type = 'refined'
                family = structure.protein_conformation.protein.family
            except:
                try:
                    structure.protein
                    self.structure = structure
                    self.structure_type = 'hommod'
                    family = structure.protein.family
                except:
                    structure.receptor_protein
                    self.structure = structure
                    self.structure_type = 'complex'
                    family = structure.receptor_protein.family
        if tm2_gn=='2x41' and tm6_gn=='6x38' and tm3_gn=='3x44' and tm7_gn=='7x52' and inactive_cutoff==2 and intermediate_cutoff==7.15:
            if family.slug.startswith('002') or family.slug.startswith('003'):
                tm6_gn, tm7_gn = '6x33', '7x51'
                inactive_cutoff, intermediate_cutoff = 2.5, 5.5
            elif family.slug.startswith('004'):
                inactive_cutoff, intermediate_cutoff = 5, 7.15
        self.tm2_gn, self.tm6_gn, self.tm3_gn, self.tm7_gn = tm2_gn, tm6_gn, tm3_gn, tm7_gn
        self.inactive_cutoff = inactive_cutoff
        self.intermediate_cutoff = intermediate_cutoff
        self.state = None
        self.activation_value = None
        self.line = False

    def run(self):
        if self.structure_type=='structure':
            self.parent_prot_conf = ProteinConformation.objects.get(protein=self.structure.protein_conformation.protein.parent)
            ssno = StructureSeqNumOverwrite(self.structure)
            ssno.seq_num_overwrite('pdb')
        elif self.structure_type=='refined':
            self.parent_prot_conf = ProteinConformation.objects.get(protein=self.structure.protein_conformation.protein)
        elif self.structure_type=='hommod':
            self.parent_prot_conf = ProteinConformation.objects.get(protein=self.structure.protein)
        elif self.structure_type=='complex':
            self.parent_prot_conf = ProteinConformation.objects.get(protein=self.structure.receptor_protein)
        # class A and T
        if self.parent_prot_conf.protein.family.slug.startswith('001') or self.parent_prot_conf.protein.family.slug.startswith('007'):
            tm6 = self.get_residue_distance(self.tm2_gn, self.tm6_gn)
            tm7 = self.get_residue_distance(self.tm3_gn, self.tm7_gn)
            print(tm6, tm7, tm6-tm7)
            if tm6 is not False and tm7 is not False:
                self.activation_value = tm6-tm7
                if self.activation_value<self.inactive_cutoff:
                    self.state = ProteinState.objects.get(slug='inactive')
                elif self.inactive_cutoff<=self.activation_value<=self.intermediate_cutoff:
                    self.state = ProteinState.objects.get(slug='intermediate')
                elif self.activation_value>self.intermediate_cutoff:
                    self.state = ProteinState.objects.get(slug='active')
        # class B
        elif self.parent_prot_conf.protein.family.slug.startswith('002') or self.parent_prot_conf.protein.family.slug.startswith('003'):
            tm2_gn_b = ResidueGenericNumberEquivalent.objects.get(default_generic_number__label=self.tm2_gn, scheme__short_name='GPCRdb(B)').label
            tm6_gn_b = ResidueGenericNumberEquivalent.objects.get(default_generic_number__label=self.tm6_gn, scheme__short_name='GPCRdb(B)').label
            tm3_gn_b = ResidueGenericNumberEquivalent.objects.get(default_generic_number__label=self.tm3_gn, scheme__short_name='GPCRdb(B)').label
            tm7_gn_b = ResidueGenericNumberEquivalent.objects.get(default_generic_number__label=self.tm7_gn, scheme__short_name='GPCRdb(B)').label

            tm6 = self.get_residue_distance(tm2_gn_b, tm6_gn_b)
            tm7 = self.get_residue_distance(tm3_gn_b, tm7_gn_b)
            if tm6 is not False and tm7 is not False:
                self.activation_value = tm6-tm7
                if self.activation_value<self.inactive_cutoff:
                    self.state = ProteinState.objects.get(slug='inactive')
                elif self.inactive_cutoff<=self.activation_value<=self.intermediate_cutoff:
                    self.state = ProteinState.objects.get(slug='intermediate')
                elif self.activation_value>self.intermediate_cutoff:
                    self.state = ProteinState.objects.get(slug='active')
        # class C
        elif self.parent_prot_conf.protein.family.slug.startswith('004'):
            tm2_gn_c = ResidueGenericNumberEquivalent.objects.get(default_generic_number__label=self.tm2_gn, scheme__short_name='GPCRdb(C)').label
            tm6_gn_c = ResidueGenericNumberEquivalent.objects.get(default_generic_number__label=self.tm6_gn, scheme__short_name='GPCRdb(C)').label
            tm3_gn_c = ResidueGenericNumberEquivalent.objects.get(default_generic_number__label=self.tm3_gn, scheme__short_name='GPCRdb(C)').label
            tm7_gn_c = ResidueGenericNumberEquivalent.objects.get(default_generic_number__label=self.tm7_gn, scheme__short_name='GPCRdb(C)').label

            tm6 = self.get_residue_distance(tm2_gn_c, tm6_gn_c)
            tm7 = self.get_residue_distance(tm3_gn_c, tm7_gn_c)
            if tm6 is not False and tm7 is not False:
                self.activation_value = tm6-tm7
                if self.activation_value<self.inactive_cutoff:
                    self.state = ProteinState.objects.get(slug='inactive')
                elif self.inactive_cutoff<=self.activation_value<=self.intermediate_cutoff:
                    self.state = ProteinState.objects.get(slug='intermediate')
                elif self.activation_value>self.intermediate_cutoff:
                    self.state = ProteinState.objects.get(slug='active')
        # class D
        #########
        # class F
        elif self.parent_prot_conf.protein.family.slug.startswith('006'):
            tm2_gn_f = ResidueGenericNumberEquivalent.objects.get(default_generic_number__label=self.tm2_gn, scheme__short_name='GPCRdb(F)').label
            tm6_gn_f = ResidueGenericNumberEquivalent.objects.get(default_generic_number__label=self.tm6_gn, scheme__short_name='GPCRdb(F)').label
            tm3_gn_f = ResidueGenericNumberEquivalent.objects.get(default_generic_number__label=self.tm3_gn, scheme__short_name='GPCRdb(F)').label
            tm7_gn_f = ResidueGenericNumberEquivalent.objects.get(default_generic_number__label=self.tm7_gn, scheme__short_name='GPCRdb(F)').label

            tm6 = self.get_residue_distance(tm2_gn_f, tm6_gn_f)
            tm7 = self.get_residue_distance(tm3_gn_f, tm7_gn_f)
            if tm6 is not False and tm7 is not False:
                self.activation_value = tm6-tm7
                if self.activation_value<0:
                    self.state = ProteinState.objects.get(slug='inactive')
                elif 0<=self.activation_value<=2:
                    self.state = ProteinState.objects.get(slug='intermediate')
                elif self.activation_value>2:
                    self.state = ProteinState.objects.get(slug='active')
        else:
            print('{} is not class A,B,C,F'.format(self.structure))
        if self.structure_type=='structure':
            ssno.seq_num_overwrite('pdb')

    def get_residue_distance(self, residue1, residue2):
        try:
            res1 = Residue.objects.get(protein_conformation__protein=self.structure.protein_conformation.protein.parent, display_generic_number__label=dgn(residue1, self.parent_prot_conf))
            res2 = Residue.objects.get(protein_conformation__protein=self.structure.protein_conformation.protein.parent, display_generic_number__label=dgn(residue2, self.parent_prot_conf))
            print(res1, res1.id, res2, res2.id)
            try:
                rota1 = Rotamer.objects.filter(structure=self.structure, residue__sequence_number=res1.sequence_number)
                if len(rota1)==0:
                    raise Exception
            except:
                rota1 = Rotamer.objects.filter(structure=self.structure, residue__display_generic_number__label=dgn(residue1, self.structure.protein_conformation))
            rota1 = right_rotamer_select(rota1, self.structure.preferred_chain[0])
            try:
                rota2 = Rotamer.objects.filter(structure=self.structure, residue__sequence_number=res2.sequence_number)
                if len(rota2)==0:
                    raise Exception
            except:
                rota2 = Rotamer.objects.filter(structure=self.structure, residue__display_generic_number__label=dgn(residue2, self.structure.protein_conformation))
            rota2 = right_rotamer_select(rota2, self.structure.preferred_chain[0])
            rotas = [rota1, rota2]
            io1 = StringIO(rotas[0].pdbdata.pdb)
            rota_struct1 = PDB.PDBParser(QUIET=True).get_structure('structure', io1)[0]
            io2 = StringIO(rotas[1].pdbdata.pdb)
            rota_struct2 = PDB.PDBParser(QUIET=True).get_structure('structure', io2)[0]

            for chain1, chain2 in zip(rota_struct1, rota_struct2):
                for r1, r2 in zip(chain1, chain2):
                    # print(self.structure, r1.get_id()[1], r2.get_id()[1], self.calculate_CA_distance(r1, r2), self.structure.state.name)
                    line = '{},{},{},{},{}\n'.format(self.structure, self.structure.state.name, round(self.calculate_CA_distance(r1, r2), 2), r1.get_id()[1], r2.get_id()[1])
                    self.line = line
                    return self.calculate_CA_distance(r1, r2)
        except:
            try:
                res1 = Residue.objects.get(protein_conformation=self.parent_prot_conf, display_generic_number__label=dgn(residue1, self.parent_prot_conf))
                res2 = Residue.objects.get(protein_conformation=self.parent_prot_conf, display_generic_number__label=dgn(residue2, self.parent_prot_conf))
                if self.structure_type=='refined':
                    pdb_data = self.structure.pdb_data.pdb
                elif self.structure_type=='hommod':
                    pdb_data = self.structure.pdb_data.pdb
                io = StringIO(pdb_data)
                struct = PDB.PDBParser(QUIET=True).get_structure('structure', io)[0]
                for chain in struct:
                    r1 = chain[res1.sequence_number]
                    r2 = chain[res2.sequence_number]
                    print(self.structure, r1.get_id()[1], r2.get_id()[1], self.calculate_CA_distance(r1, r2), self.structure.state.name)
                    line = '{},{},{},{},{}\n'.format(self.structure, self.structure.state.name, round(self.calculate_CA_distance(r1, r2), 2), r1.get_id()[1], r2.get_id()[1])
                    self.line = line
                    return self.calculate_CA_distance(r1, r2)

            except:
                print('Error: {} no matching rotamers ({}, {})'.format(self.structure.pdb_code.index, residue1, residue2))
                return False

    def calculate_CA_distance(self, residue1, residue2):
        diff_vector = residue1['CA'].get_coord()-residue2['CA'].get_coord()
        return numpy.sqrt(numpy.sum(diff_vector * diff_vector))


class StructureSeqNumOverwrite():
    def __init__(self, structure):
        self.structure = structure
        path = os.sep.join([settings.DATA_DIR, 'structure_data','wt_pdb_lookup', '{}.json'.format(self.structure.pdb_code.index)])
        if os.path.isfile(path):
            with open(path, 'r') as lookup_file:
                self.lookup = json.load(lookup_file)
            self.wt_pdb_table, self.pdb_wt_table = OrderedDict(), OrderedDict()
            for i in self.lookup:
                self.wt_pdb_table[i['WT_POS']] = i['PDB_POS']
                self.pdb_wt_table[i['PDB_POS']] = i['WT_POS']
        else:
            self.lookup = OrderedDict()
            self.wt_pdb_table = OrderedDict()
            self.pdb_wt_table = OrderedDict()

    def seq_num_overwrite(self, overwrite_target):
        ''' Overwrites Residue object sequence numbers in GPCRDB
            @param overwrite_target: 'pdb' if converting pdb to wt, 'wt' if the other way around
        '''
        resis = Residue.objects.filter(protein_conformation=self.structure.protein_conformation)
        if overwrite_target=='pdb':
            target_dict = self.pdb_wt_table
        elif overwrite_target=='wt':
            target_dict = self.wt_pdb_table
        for r in resis:
            if r.sequence_number in target_dict:
                r.sequence_number = int(target_dict[r.sequence_number])
                r.save()

class ParseAFModelsCSV():
    def __init__(self):
        self.complexes = []
        self.structures = {}
        with open(os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'AF_Models.csv']), newline='') as csvfile:
            structs = csv.reader(csvfile, delimiter=',')
            next(structs, None)
            for s in structs:
                #0:complex 1:finished 2:receptor 3:ligand 4:outside_found
                #5:chosen_rank 6:score 7:relaxed 8:closer_pep_term 9:has_x
                protein = s[0].split('_')[0]
                ligand = s[3]
                organism = s[0].split('-')[0].split('_')[1]
                # TODO: THIS IS HARDCODED BUT NEED TO BE CHANGED
                rankedpdbname = s[0] + '-rank' + s[5] + '.pdb'
                location = os.sep.join([settings.DATA_DIR, 'structure_data', 'af_peptide', rankedpdbname]) #'/protwis/data/protwis/gpcr/structure_data/af_peptide/'+s[0]+'-rank'+s[5]+'.pdb'
                cmpx = s[0] #protein+'_'+ligand
                self.complexes.append(cmpx)
                self.structures[cmpx]= {'protein':s[0].split('-')[0], 'name':protein.lower(), 'state':'Not defined', 'peptide_id': ligand, 'preferred_chain':'A', 'model':'af-peptide', 'location':location, 'organism':organism}

    def parse_ligands(self):
        endogenous_peptides = pd.read_csv(os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'Peptides_GPRCdb_dump.csv']))
        for complex in self.structures.keys():
            receptor, pep_id = complex.split('_')[0], complex.split('-')[1]
            self.structures[complex]['ligand'] = []
            try:
                role = endogenous_peptides.loc[(endogenous_peptides['Receptor'].str.match(receptor)) & (endogenous_peptides['GtP ID'] == int(pep_id)), 'Role'].values[0]
                name = endogenous_peptides.loc[(endogenous_peptides['Receptor'].str.match(receptor)) & (endogenous_peptides['GtP ID'] == int(pep_id)), 'Ligand Name'].values[0]
                type = endogenous_peptides.loc[(endogenous_peptides['Receptor'].str.match(receptor)) & (endogenous_peptides['GtP ID'] == int(pep_id)), 'Type'].values[0]
                gpcrdb_id = endogenous_peptides.loc[(endogenous_peptides['Receptor'].str.match(receptor)) & (endogenous_peptides['GtP ID'] == int(pep_id)), 'GPCRdb ID'].values[0]
                self.structures[complex]['ligand'].append(
                                {'chain':'B',
                                 'name':pep_id,
                                 'role':role,
                                 'title':name,
                                 'type':type,
                                 'gpcrdb id': gpcrdb_id,
                                 'in_structure': True})
            except IndexError:
                print('Cannot find information for complex {}'.format(complex))

class ParseAFComplexModels():

    residue_to_one_letter = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }

    def __init__(self):
        self.data_dir = os.sep.join([settings.DATA_DIR, 'structure_data', 'AlphaFold_multimer'])
        self.filedirs = os.listdir(self.data_dir)
        self.complexes = {}

        for f in self.filedirs:
            metrics_file = os.sep.join([self.data_dir, f, f+'_metrics.csv'])
            if not os.path.exists(metrics_file):
                metrics_file = os.sep.join([self.data_dir, f, f+'.csv'])
            metrics = [row for row in csv.DictReader(open(metrics_file, 'r'))][0]
            location = os.sep.join([self.data_dir, f, f+'.pdb'])

            parts = f.split('-')
            receptor = parts[0]
            if len(parts) == 3:  # Case with peptide
                peptide_id = parts[1]
                peptide = "-" + peptide_id
                signprot = parts[2]
                chain_e_sequence = self.get_ligand_sequence(location, 'E')
            else:  # Case without peptide
                peptide = None
                signprot = parts[1]
                chain_e_sequence = None

            # Grab model date/version from pdb file
            with open(location, 'r') as model_file:
                line = model_file.readlines()[0]
                print('LINE')
                print(line)
                date_re = re.search('HEADER[A-Z\S\D]+(\d{4}-\d{2}-\d{2})', line)
                try:
                    model_date = date_re.group(1)
                except AttributeError:
                    current_date = datetime.now().date()
                    model_date = current_date.strftime('%Y-%m-%d')

            # Check signprot type

            if signprot.startswith('gna'):

                # Check if model has full heterotrimer
                if 'gbb1_human' in f:
                    signprot = signprot.split('_')[0]+'_human'
                    beta_gamma = True
                else:
                    beta_gamma = False

            else:
                beta_gamma = False

            self.complexes[f'{receptor}{peptide}-{signprot}'] = {
                'receptor': receptor,
                'peptide': peptide,
                'signprot': signprot,
                'publication_date': model_date,
                'location': location,
                'model': 'af-signprot',
                'preferred_chain': 'A',
                'PTM': metrics['ptm'],
                'iPTM': metrics['iptm'],
                'PAE_mean': metrics['pae_mean'],
                'chain_e_sequence': chain_e_sequence
            }

            # Check for type of signprot
            if signprot.startswith('gna'):
                # Check if model has full heterotrimer
                if 'gbb1_human' in f:
                    signprot = signprot.split('_')[0] + '_human'
                    beta_gamma = True
                else:
                    beta_gamma = False
            else:
                self.complexes[f'{receptor}{peptide}-{signprot}']['model'] = 'af-arrestin'
                beta_gamma = False

            self.complexes[f'{receptor}{peptide}-{signprot}']['beta_gamma'] = beta_gamma

            # self.complexes[f'{receptor}{peptide}-{signprot}'] = complex_info

    def get_ligand_sequence(self, pdb_file, chain_id):
        sequence = ""
        with open(pdb_file, 'r') as file:
            for line in file:
                if line.startswith("ATOM") and line[21] == chain_id:
                    residue = line[17:20]
                    res_seq = int(line[22:26])
                    if len(sequence) < res_seq:
                        sequence += self.residue_to_one_letter.get(residue, 'X')
        return sequence

class ParseRFAAModels():
    def __init__(self):
        self.data_dir = os.sep.join([settings.DATA_DIR, 'structure_data', 'RFAA'])
        self.filedirs = os.listdir(self.data_dir)
        self.complexes = {}
        for f in self.filedirs:
            parts = f.split('-', 1)
            receptor = parts[0]
            ligand_secs = parts[1].split('[')
            ligand = ligand_secs[0]
            inchikey = ligand_secs[1].replace(']', '')
            model = 'af-rfaa-sm'

            metrics_file = os.sep.join([self.data_dir, f, receptor + '-' + ligand +'_metrics.csv'])
            metrics = [row for row in csv.DictReader(open(metrics_file, 'r'))][0]
            location = os.sep.join([self.data_dir, f, f+'.pdb'])
            ### Grab model date/version from pdb file
            with open(location, 'r') as model_file:
                line = model_file.readlines()[0]
                date_re = re.search('HEADER[A-Z\S\D]+(\d{4}-\d{2}-\d{2})', line)
                try:
                    model_date = date_re.group(1)
                except:
                    model_date = '2024-05-15'
            self.complexes[f'{receptor}-{inchikey}'] = {
                'receptor':receptor, 
                'ligand_id':inchikey,
                'publication_date':model_date, 
                'location':location, 
                'model': model, 
                'preferred_chain':'A', 
                'pae_7tm':metrics['pae_7tm'], 
                'plddt_mean':metrics['plddt_sm_mean'], 
                }


class ParseStructureCSV():
    def __init__(self):
        self.pdb_ids = []
        self.structures = {}
        self.parent_segends = {}
        with open(os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'structures.csv']), newline='') as csvfile:
            structures = csv.reader(csvfile, delimiter='\t')
            next(structures, None)
            for s in structures:
                self.pdb_ids.append(s[0])
                # Chain ID format check - remove .0 from number type chain ID e.g. 8JCU
                if '.' in s[5]:
                    s[5] = s[5][0]
                self.structures[s[0]]= {'protein':s[1], 'name':s[0].lower(), 'state':s[4], 'preferred_chain':s[5], 'resolution':s[3], 'date_from_file':s[7], 'method_from_file':s[2]}
        self.fusion_proteins = []

    def __str__(self):
        return '<ParsedStructures: {} entries>'.format(len(self.pdb_ids))

    def parse_ligands(self):
        with open(os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'ligands.csv']), newline='') as csvfile:
            ligands = csv.reader(csvfile, delimiter='\t')
            next(ligands, None)
            for ligand in ligands:
                if 'ligand' not in self.structures[ligand[0]]:
                    self.structures[ligand[0]]['ligand'] = []
                in_structure = True
                if ligand[8]!='':
                    in_structure = False
                self.structures[ligand[0]]['ligand'].append({'chain':ligand[1], 'name':ligand[2], 'pubchemId':ligand[3], 'role':ligand[4], 'title':ligand[5], 'type': ligand[6], 'in_structure': in_structure})

    def parse_nanobodies(self):
        self.parse_aux_file('nanobodies.csv')

    def parse_fusion_proteins(self):
        self.parse_aux_file('fusion_proteins.csv')

    def parse_ramp(self):
        self.parse_aux_file('ramp.csv')

    def parse_grk(self):
        self.parse_aux_file('grk.csv')

    def parse_g_proteins(self):
        with open(os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'g_proteins.csv']), newline='') as csvfile:
            g_proteins = csv.reader(csvfile, delimiter='\t')
            next(g_proteins, None)
            for g in g_proteins:
                self.structures[g[0]]['g_protein'] = {'alpha_uniprot': g[1], 'alpha_chain': g[2], 'beta_uniprot': g[3], 'beta_chain': g[4], 'gamma_uniprot': g[5], 'gamma_chain': g[6], 'note': g[7]}

    def parse_arrestins(self):
        with open(os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'arrestins.csv']), newline='') as csvfile:
            arrestins = csv.reader(csvfile, delimiter='\t')
            next(arrestins, None)
            for a in arrestins:
                self.structures[a[0]]['arrestin'] = {'protein': a[1], 'chain': a[2], 'note': a[3]}

    def parse_aux_file(self, aux_csv):
        with open(os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', aux_csv]), newline='') as csvfile:
            aux = csv.reader(csvfile, delimiter='\t')
            next(aux, None)
            for a in aux:
                if 'auxiliary_protein' not in self.structures[a[0]]:
                    self.structures[a[0]]['auxiliary_protein'] = []
                self.structures[a[0]]['auxiliary_protein'].append(a[1])
                if aux_csv=='fusion_proteins.csv':
                    if a[1] not in self.fusion_proteins:
                        self.fusion_proteins.append(a[1])

    @staticmethod
    def parse_yaml_file(yaml_file):
        with open(os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', yaml_file]), 'r') as yfile:
            y = yaml.safe_load(yfile)
        return y

    def parse_parent_segends(self):
        self.parent_segends = self.parse_yaml_file('non_xtal_segends.yaml')


class StructureBuildCheck():
    local_annotation_dir = os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation'])
    local_wt_pdb_lookup_dir = os.sep.join([settings.DATA_DIR, 'structure_data', 'wt_pdb_lookup'])
    local_g_protein_chimeras_gapped = os.sep.join([settings.DATA_DIR, 'g_protein_data', 'g_protein_chimeras_gapped.fasta'])

    def __init__(self):
        with open(self.local_annotation_dir+'/mod_xtal_segends.yaml', 'r') as f:
            self.segends_dict = yaml.safe_load(f)
        self.pdbs = ParseStructureCSV().pdb_ids
        self.wt_pdb_lookup_files = [i.split('.')[0] for i in os.listdir(self.local_wt_pdb_lookup_dir)]
        self.missing_seg = []
        self.start_error = []
        self.end_error = []
        self.helix_length_error = []
        self.duplicate_residue_error = {}
        self.g_protein_chimeras = SeqIO.to_dict(SeqIO.parse(open(self.local_g_protein_chimeras_gapped), "fasta"))
        self.g_prot_test_exceptions = {}
        for i, j in self.g_protein_chimeras.items():
            if '|' not in i and '-' in j.seq:
                self.g_prot_test_exceptions[i] = j.seq.count('-')

    def check_structures(self):
        for pdb in self.pdbs:
            try:
                Structure.objects.get(pdb_code__index=pdb)
            except Structure.DoesNotExist:
                print('Error: {} Structure object has not been built'.format(pdb))

    def check_duplicate_residues(self, structs, check_on='sequence_number'):
        for s in structs:
            if hasattr(s, 'protein_conformation'):
                resis = Residue.objects.filter(protein_conformation=s.protein_conformation)
            else:
                if s.protein.family.slug.startswith('100'):
                    tag = '_a'
                elif s.protein.family.slug.startswith('200'):
                    tag = '_arrestin'
                resis = Residue.objects.filter(protein_conformation__protein__entry_name=s.structure.pdb_code.index.lower()+tag)
            if len(resis)!=len(resis.values_list(check_on).distinct()):
                for r in resis:
                    try:
                        if hasattr(s, 'protein_conformation'):
                            if check_on=='sequence_number':
                                Residue.objects.get(protein_conformation=s.protein_conformation, sequence_number=r.sequence_number)
                            elif check_on=='display_generic_number':
                                Residue.objects.get(protein_conformation=s.protein_conformation, display_generic_number=r.display_generic_number)
                        else:
                            if check_on=='sequence_number':
                                re = Residue.objects.get(protein_conformation__protein__entry_name=s.structure.pdb_code.index.lower()+tag, sequence_number=r.sequence_number)
                            elif check_on=='display_generic_number':
                                re = Residue.objects.get(protein_conformation__protein__entry_name=s.structure.pdb_code.index.lower()+tag, display_generic_number=r.display_generic_number)
                    except Residue.MultipleObjectsReturned:
                        if s in self.duplicate_residue_error:
                            self.duplicate_residue_error[s].append(r)
                        else:
                            self.duplicate_residue_error[s] = [r]

    def check_segment_ends(self, structure):
        key = structure.pdb_code.index
        if key in self.segends_dict:
            parent_residues = Residue.objects.filter(protein_conformation__protein=structure.protein_conformation.protein.parent)
            structure_residues = Residue.objects.filter(protein_conformation=structure.protein_conformation)
            annotation = self.segends_dict[key]
            for i in range(1,9):
                if annotation[str(i)+'e']=='-':
                    continue
                else:
                    anno_b = int(annotation[str(i)+'b'])
                    anno_e = int(annotation[str(i)+'e'])
                    if i==8:
                        seg = 'H8'
                    else:
                        seg = 'TM'+str(i)
                    parent_seg_resis = parent_residues.filter(protein_segment__slug=seg)
                    seg_resis = structure_residues.filter(protein_segment__slug=seg)
                    if len(seg_resis)==0:
                        self.missing_seg.append([structure, seg, anno_b, anno_e])
                        continue
                    if i<8 and len(seg_resis)<5:
                        self.helix_length_error.append([structure, seg, len(seg_resis)])
                    if seg_resis[0].sequence_number!=anno_b:
                        if parent_seg_resis[0].sequence_number!=seg_resis[0].sequence_number:
                            if structure.pdb_code.index in self.wt_pdb_lookup_files:
                                with open(os.sep.join([self.local_wt_pdb_lookup_dir, structure.pdb_code.index+'.json']), 'r') as wt_pdb:
                                    lookup = json.load(wt_pdb)
                                    found_s = False
                                    for i in lookup:
                                        if i['PDB_POS']==seg_resis[0].sequence_number:
                                            found_s = True
                                            if i['WT_POS']!=parent_seg_resis[0].sequence_number:
                                                self.start_error.append([structure, seg, seg_resis[0].sequence_number, anno_b])
                                                break
                                            else:
                                                break
                                    if not found_s:
                                        self.start_error.append([structure, seg, seg_resis[0].sequence_number, anno_b])
                            else:
                                self.start_error.append([structure, seg, seg_resis[0].sequence_number, anno_b])
                    if seg_resis.reverse()[0].sequence_number!=anno_e:
                        if parent_seg_resis.reverse()[0].sequence_number!=seg_resis.reverse()[0].sequence_number:
                            if structure.pdb_code.index in self.wt_pdb_lookup_files:
                                with open(os.sep.join([self.local_wt_pdb_lookup_dir, structure.pdb_code.index+'.json']), 'r') as wt_pdb:
                                    lookup = json.load(wt_pdb)
                                    found_e = False
                                    for i in lookup:
                                        if i['PDB_POS']==seg_resis.reverse()[0].sequence_number:
                                            found_e = True
                                            if i['WT_POS']!=parent_seg_resis.reverse()[0].sequence_number:
                                                self.end_error.append([structure, seg, seg_resis.reverse()[0].sequence_number, anno_e])
                                                break
                                            else:
                                                break
                                    if not found_e:
                                        self.end_error.append([structure, seg, seg_resis.reverse()[0].sequence_number, anno_e])
                            else:
                                self.end_error.append([structure, seg, seg_resis.reverse()[0].sequence_number, anno_e])
        else:
            print('Warning: {} not annotated'.format(key))

    def check_signprot_struct_residues(self, signprot_complex):
        pdb = PDBParser(PERMISSIVE=True, QUIET=True).get_structure('struct', StringIO(str(signprot_complex.structure.pdb_data.pdb)))[0]
        if signprot_complex.protein.family.slug.startswith('100'):
            resis = Residue.objects.filter(protein_conformation=ProteinConformation.objects.get(protein__entry_name=signprot_complex.structure.pdb_code.index.lower()+'_a'))
            H5 = resis.filter(protein_segment__slug='G.H5')
            if len(H5)<10:
                self.missing_seg.append(['G.H5', signprot_complex.structure, len(H5)])
        elif signprot_complex.protein.family.slug.startswith('200'):
            resis = Residue.objects.filter(protein_conformation=ProteinConformation.objects.get(protein__entry_name=signprot_complex.structure.pdb_code.index.lower()+'_arrestin'))
        ppb = PDB.PPBuilder()
        seq = ""
        for pp in ppb.build_peptides(pdb[signprot_complex.alpha], aa_only=False):
            seq += str(pp.get_sequence())
        ### Skip if difference in length comes from unmappable residues from structure chimera
        if signprot_complex.protein.entry_name in self.g_prot_test_exceptions and len(seq)-len(resis)<=self.g_prot_test_exceptions[signprot_complex.protein.entry_name]:
            return 0
        ### Print structures where residues were not built
        if len(seq)!=len(resis):
            print(signprot_complex.structure, len(pdb[signprot_complex.alpha]), len(resis))


class ModelRotamer(object):
    def __init__(self):
        self.backbone_template = None
        self.rotamer_template = None


class X50Finder():
    '''Find corresponding x50 positions based on best BLAST hit in db
    '''
    def __init__(self, uniprot_file, debug=False):
        self.uniprot_file = uniprot_file
        p = PDBParser()
        self.biopdb = p.get_structure('structure', self.uniprot_file)
        self.top_hit = None
        self.debug = debug

    def get_sequence_from_structure(self):
        structure_seq = {}
        for chain in self.biopdb[0]:
            chain_id = chain.get_id()
            structure_seq[chain_id] = ''
            for res in chain:
                if "CA" in res and res.id[0]==" ":
                    structure_seq[chain_id]+=polypeptide.three_to_one(res.get_resname())
        return structure_seq

    def run(self):
        structure_seq = self.get_sequence_from_structure()
        self.blast_search = BlastSearch(top_results=25)
        out_x50 = {}
        for chain, seq in structure_seq.items():
            out_x50[chain] = {}
            blast_output = self.blast_search.run(seq)
            found_good_match = False
            for b in blast_output:
                ref = Protein.objects.get(id=int(b[0]))
                aligner = PairwiseAligner()
                aligner.open_gap_score = -2
                aligner.extend_gap_score = -1
                aligner.mode = 'global'
                pw2 = aligner.align(ref.sequence, seq)

                alignment_fragments = [i for i in str(pw2[0][1]).split('-') if i!='' and len(i)<10]
                if len(alignment_fragments)<4:
                    found_good_match = True
                    break
            self.top_hit = ref

            if not found_good_match:
                print('ERROR: no good pairwise alignment for {}'.format(self.uniprot_file.split('/')[-1]))

            ref_seq, temp_seq = str(pw2[0][0]), str(pw2[0][1])
            if self.debug:
                for r, t in zip(ref_seq, temp_seq):
                    print(r,t)

            x50s = Residue.objects.filter(protein_conformation__protein=ref, display_generic_number__label__endswith='x50')
            indeces = {}
            for x in x50s:
                c, i = 1, 0
                for r in ref_seq:
                    if r!='-':
                        c+=1
                    i+=1
                    if c==x.sequence_number:
                        break
                indeces[x.display_generic_number.label] = i

            for gn, i in indeces.items():
                before = len(temp_seq[:i+1].replace('-',''))
                after = len(temp_seq[i:].replace('-',''))
                if after==0 or (temp_seq[i]=='-' and temp_seq[i+1]=='-' and temp_seq[i+2]=='-'):
                    out_x50[chain][gn] = '-'
                else:
                    out_x50[chain][gn] = before

        return out_x50


def update_template_source(template_source, keys, struct, segment, just_rot=False):
    ''' Update the template_source dictionary with structure info for backbone and rotamers.
    '''
    for k in keys:
        if just_rot is True:
            try:
                template_source[segment][k][1] = struct
            except:
                pass
        else:
            try:
                template_source[segment][k][0] = struct
                template_source[segment][k][1] = struct
            except:
                pass
    return template_source

def compare_and_update_template_source(template_source, segment, signprot_pdb_array, i, cgn, template_source_key, segs_for_alt_complex_struct, alt_complex_struct, main_structure):
    if cgn in signprot_pdb_array[segment]:
        if segment in segs_for_alt_complex_struct and signprot_pdb_array[segment][cgn]!='x':
            update_template_source(template_source, [str(template_source_key)], alt_complex_struct, segment)
        elif signprot_pdb_array[segment][cgn]!='x':
            update_template_source(template_source, [str(template_source_key)], main_structure, segment)
        else:
            update_template_source(template_source, [str(template_source_key)], None, segment)
    else:
        update_template_source(template_source, [str(template_source_key)], None, segment)
    return template_source


def right_rotamer_select(rotamer, chain=None):
    ''' Filter out compound rotamers.
    '''
    if len(rotamer)>1:
        for i in rotamer:
            if i.pdbdata.pdb.startswith('COMPND') is False:
                if chain!=None and chain==i.pdbdata.pdb[21]:
                    rotamer = i
                    break
    else:
        rotamer=rotamer[0]
    return rotamer


def get_pdb_ids(uniprot_id):
    pdb_list = []
    data = { "query": {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute" : "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                    "operator":"in",
                    "value":[ uniprot_id ]
                }
            },
            "request_options":{
                "return_all_hits": True
            },
            "return_type": "entry" }
    url = 'https://search.rcsb.org/rcsbsearch/v2/query'
    req = Request(url)
    req.add_header('Content-Type', 'application/json; charset=utf-8')
    jsondata = json.dumps(data)
    jsondataasbytes = jsondata.encode('utf-8')
    req.add_header('Content-Length', len(jsondataasbytes))
    response = urlopen_with_retry(req, jsondataasbytes)
    if response is False or response is None:
        logger.warning(f"ERROR: no valid response from RCSB for accession code {uniprot_id}")
        return []

    rr = response.read()
    if len(rr)==0:
        return []
    out = json.loads(rr)
    for i in out['result_set']:
        pdb_list.append(i['identifier'])

    response.close()
    return pdb_list


def create_structure_rotamer(PDB_residue, residue_object, structure):
    atom_num_dict = {'E':9, 'S':6, 'Y':12, 'G':4, 'A':5, 'V':7, 'M':8, 'L':8, 'I':8, 'T':7, 'F':11, 'H':10, 'K':9,
                 'D':8, 'C':6, 'R':11, 'P':7, 'Q':9, 'N':8, 'W':14, '-':0}
    out_stream = StringIO()
    io = PDBIO()
    io.set_structure(PDB_residue)
    io.save(out_stream)
    pdbdata = PdbData.objects.get_or_create(pdb=out_stream.getvalue())[0]
    missing_atoms = atom_num_dict[polypeptide.three_to_one(PDB_residue.get_resname())] > len(PDB_residue.get_unpacked_list())
    rot = Rotamer(missing_atoms=missing_atoms, pdbdata=pdbdata, residue=residue_object, structure=structure)
    return rot

def fetch_signprot_data(pdb, protein, beta_uniprots=[], gamma_uniprots=[]):
    data = {}

    response = urllib.request.urlopen("https://data.rcsb.org/rest/v1/core/entry/{}".format(pdb))
    json_data = json.loads(response.read())
    response.close()

    data["method"] = json_data["exptl"][0]["method"]
    if data["method"].startswith("THEORETICAL") or data["method"] in ["SOLUTION NMR","SOLID-STATE NMR"]:
        return None
    if "citation" in json_data and "pdbx_database_id_doi" in json_data["citation"]:
        data["doi"] = json_data["citation"]["pdbx_database_id_doi"]
    else:
        data["doi"] = None
    if "pubmed_id" in json_data["rcsb_entry_container_identifiers"]:
        data["pubmedId"] = json_data["rcsb_entry_container_identifiers"]["pubmed_id"]
    else:
        data["pubmedId"] = None

    data["release_date"] = json_data["rcsb_accession_info"]["initial_release_date"][:10]
    data["resolution"] = json_data["rcsb_entry_info"]["resolution_combined"][0]
    entities_num = len(json_data["rcsb_entry_container_identifiers"]["polymer_entity_ids"])
    data["alpha"] = protein.accession
    data["alpha_chain"] = None
    data["alpha_coverage"] = None
    data["beta"] = None
    data["beta_chain"] = None
    data["gamma"] = None
    data["gamma_chain"] = None
    data["other"] = []
    for i in range(1,entities_num+1):
        response = urllib.request.urlopen("https://data.rcsb.org/rest/v1/core/polymer_entity/{}/{}".format(pdb, i))
        json_data = json.loads(response.read())
        response.close()
        if "uniprot_ids" in json_data["rcsb_polymer_entity_container_identifiers"]:
            for j, u_id in enumerate(json_data["rcsb_polymer_entity_container_identifiers"]["uniprot_ids"]):
                if u_id+".txt" in beta_uniprots:
                    data["beta"] = u_id
                    data["beta_chain"] = json_data["rcsb_polymer_entity_container_identifiers"]["auth_asym_ids"][j][0]
                elif u_id+".txt" in gamma_uniprots:
                    data["gamma"] = u_id
                    data["gamma_chain"] = json_data["rcsb_polymer_entity_container_identifiers"]["auth_asym_ids"][j][0]
                elif u_id==protein.accession:
                    data["alpha"] = u_id
                    data["alpha_coverage"] = json_data["entity_poly"]["rcsb_sample_sequence_length"]
                    try:
                        data["alpha_chain"] = json_data["rcsb_polymer_entity_container_identifiers"]["auth_asym_ids"][j][0]
                    except IndexError as e:
                        data["alpha_chain"] = json_data["rcsb_polymer_entity_container_identifiers"]["auth_asym_ids"][j-1][0]
                else:
                    if json_data["rcsb_polymer_entity"]["pdbx_description"] not in data["other"]:
                        data["other"].append(json_data["rcsb_polymer_entity"]["pdbx_description"])
        else:
            if json_data["rcsb_polymer_entity"]["pdbx_description"] not in data["other"]:
                data["other"].append(json_data["rcsb_polymer_entity"]["pdbx_description"])
    return data


def build_signprot_struct(protein, pdb, data):
    ss = SignprotStructure()
    pdb_code, p_c = WebLink.objects.get_or_create(index=pdb, web_resource=WebResource.objects.get(slug="pdb"))
    pub_date = data["release_date"]
    # Structure type
    if "x-ray" in data["method"].lower():
        structure_type_slug = "x-ray-diffraction"
    elif "electron" in data["method"].lower():
        structure_type_slug = "electron-microscopy"
    else:
        structure_type_slug = "-".join(data["method"].lower().split(" "))
    try:
        structure_type = StructureType.objects.get(slug=structure_type_slug)
    except StructureType.DoesNotExist as e:
        structure_type, c = StructureType.objects.get_or_create(slug=structure_type_slug, name=data["method"])
        # self.logger.info("Created StructureType:"+str(structure_type))

    # Publication
    pub = None
    if data["doi"]:
        try:
            pub = Publication.get_or_create_from_doi(data["doi"])
        except:
            # 2nd try (in case of paralellization clash)
            pub = Publication.get_or_create_from_doi(data["doi"])
    elif data["pubmedId"]:
        try:
            pub = Publication.get_or_create_from_pubmed(data["pubmedId"])
        except:
            # 2nd try (in case of paralellization clash)
            pub = Publication.get_or_create_from_pubmed(data["pubmedId"])

    # PDB data
    url = 'https://www.rcsb.org/pdb/files/{}.pdb'.format(pdb)
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req) as response:
        pdbdata_raw = response.read().decode('utf-8')
    pdbdata_object = PdbData.objects.get_or_create(pdb=pdbdata_raw)[0]
    ss.pdb_code = pdb_code
    ss.structure_type = structure_type
    ss.resolution = data["resolution"]
    ss.publication_date = pub_date
    ss.publication = pub
    ss.protein = protein
    ss.pdb_data = pdbdata_object
    if len(SignprotStructure.objects.filter(pdb_code=ss.pdb_code))>0:
        # self.logger.warning('SignprotStructure {} already created, skipping'.format(pdb_code))
        return 0
    ss.save()
    # Stabilizing agent
    for o in data["other"]:
        if len(o)>75:
            continue
        if o=="REGULATOR OF G-PROTEIN SIGNALING 14":
            o = "Regulator of G-protein signaling 14"
        elif o=="Nanobody 35":
            o = "Nanobody-35"
        elif o=="ADENYLATE CYCLASE, TYPE V":
            o = "Adenylate cyclase, type V"
        elif o=="1-phosphatidylinositol-4,5-bisphosphate phosphodiesterase beta-3":
            o = "1-phosphatidylinositol 4,5-bisphosphate phosphodiesterase beta-3"
        elif o=="FAB30 HEAVY CHAIN":
            o = "Fab30 Heavy Chain"
        elif o=="FAB30 LIGHT CHAIN":
            o = "Fab30 Light Chain"
        elif o=="VASOPRESSIN V2 RECEPTOR PHOSPHOPEPTIDE":
            o = "Vasopressin V2 receptor phosphopeptide"
        elif o=="AP-2 COMPLEX SUBUNIT BETA-1":
            o = "AP-2 complex subunit Beta-1"
        stabagent, sa_created = StructureStabilizingAgent.objects.get_or_create(slug=o.replace(" ","-").replace(" ","-"), name=o)
        ss.stabilizing_agents.add(stabagent)
    ss.save()
    return ss

def flip_residue(atoms, atom_type):
    for a in atoms:
        if a.get_id()==atom_type+'1':
            one_index = atoms.index(a)
            one_coords = a.get_coord()
        elif a.get_id()==atom_type+'2':
            atoms[one_index].coord = a.get_coord()
            a.coord = one_coords
    return atoms

def run_residue_flip(atoms, atom_types=None):
    if not atom_types:
        atom_types = ['CD','CE','CG','OE','OD','NH']
    for at in atom_types:
        atoms = flip_residue(atoms, at)
    return atoms

def atoms_to_dict(atom_list):
    prev_res = 0
    atom_resis = OrderedDict()
    for a in atom_list:
        if a.get_parent().get_id()[1]!=prev_res:
            atom_resis[a.get_parent().get_id()[1]] = [a]
        else:
            atom_resis[a.get_parent().get_id()[1]].append(a)
        prev_res = a.get_parent().get_id()[1]
    return atom_resis
