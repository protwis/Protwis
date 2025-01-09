from build.management.commands.base_build import Command as BaseBuild
from structure.models import Structure, StructureModel
from Bio.PDB import PDBIO, PDBParser, Select
from Bio.PDB.Model import Model
from structure.assign_generic_numbers_gpcr import GenericNumbering as as_gn
from io import StringIO
import os
import shutil
from protwis import settings

class ResidueBFactorSelect(Select):
    """
    A selection class for filtering residues based on the B-factor of their CA atoms.
    Only residues with CA atom B-factors within the specified range are selected.
    """
    def __init__(self, bfactor_range=(-8.0, 8.0)):
        """
        Initializes the selection object with the specified B-factor range.

        Parameters:
        -----------
        bfactor_range : tuple, optional
            The inclusive range of B-factors to select residues. Default is (-8.0, 8.0).
            Negative range is to catch bulges.
        """
        self.bfactor_range = bfactor_range

    def accept_residue(self, residue):
        """
        Checks if a residue should be accepted based on the B-factor of its CA atom.

        Parameters:
        -----------
        residue : Bio.PDB.Residue
            The residue to check.

        Returns:
        --------
        bool
            True if the residue should be accepted (CA atom's B-factor is within the range), False otherwise.
        """
        for atom in residue:
            if atom.get_id() == 'CA':
                bfactor = atom.get_bfactor()
                if self.bfactor_range[0] < bfactor < self.bfactor_range[1]:
                    return True
        return False

def extract_preferred_chain(pdb_data, preferred_chain, pdb_code):
    """
    Extracts the preferred chain from the PDB data and returns a model containing only that chain.

    Parameters:
    -----------
    pdb_data : str
        The PDB data as a string.
    preferred_chain : str
        The identifier of the preferred chain to extract.
    pdb_code : str
        The PDB code for identification purposes.

    Returns:
    --------
    Bio.PDB.Model.Model or None
        A model containing only the preferred chain, or None if the chain is not found.
    """
    parser = PDBParser(QUIET=True)
    pdb_io = StringIO(pdb_data)
    structure = parser.get_structure(pdb_code, pdb_io)

    # Get the first model
    model = structure[0]  # Assuming only one model

    # Create a new model with the same id
    new_model = Model(model.id)

    # Extract the preferred chain
    if preferred_chain in model:
        chain = model[preferred_chain]
        new_model.add(chain)
        return new_model
    else:
        print(f"Chain {preferred_chain} not found in structure {pdb_code}")
        return None

def create_db_dir(structure_type, output_dir):
    """
    Creates an output directory for the specified structure type.

    Parameters:
    -----------
    structure_type : str
        The structure type ('raw', 'ref', or 'af').
    output_dir : str
        The base output directory.

    Returns:
    --------
    str
        The path to the created output directory.
    """
    db_dir = f'{structure_type}_foldseek_db_trim'
    output_db_dir = os.path.join(output_dir, db_dir)
    if os.path.exists(output_db_dir):
        shutil.rmtree(output_db_dir)
    os.mkdir(output_db_dir)
    return output_db_dir

def process_structure(pdb_code, pdb_data, output_filename, preferred_chain=None):
    """
    Processes a single structure: parses PDB data, assigns generic numbers, saves annotated structure.

    Parameters:
    -----------
    pdb_code : str
        The PDB code or identifier of the structure.
    pdb_data : str
        The PDB data as a string.
    output_filename : str
        The filename to save the annotated structure.
    preferred_chain : str, optional
        The identifier of the preferred chain to extract. If None, the entire structure is used.
    """
    parser = PDBParser(QUIET=True)
    pdb_io = StringIO(pdb_data)
    structure = parser.get_structure(pdb_code, pdb_io)

    if preferred_chain:
        # Extract the preferred chain
        model = structure[0]  # Assuming only one model
        if preferred_chain in model:
            chain = model[preferred_chain]
            # Create a new model with the chain
            new_model = Model(model.id)
            new_model.add(chain)
            structure_to_use = new_model
        else:
            print(f"Chain {preferred_chain} not found in structure {pdb_code}")
            return
    else:
        # Use the entire structure
        structure_to_use = structure[0]  # Assuming first model

    # Assign generic numbers
    gn = as_gn(structure=structure_to_use, pdb_code=pdb_code)
    gn.assign_generic_numbers()
    annotated_structure = gn.get_annotated_structure()

    # Save annotated structure
    io = PDBIO()
    io.set_structure(annotated_structure)
    io.save(output_filename, ResidueBFactorSelect())
    print(f"Saved selected residues to {output_filename}")

def process_structures(structure_type, output_dir):
    """
    Processes structures of the given structure_type and saves the annotated structures.

    Parameters:
    -----------
    structure_type : str
        The structure type to process ('raw', 'ref', or 'af').
    output_dir : str
        The base output directory.
    """
    output_db_dir = create_db_dir(structure_type, output_dir)

    if structure_type == 'raw':
        # Query experimental structures from the database
        exp_structures = Structure.objects.filter(
            structure_type__slug__in=['x-ray-diffraction', 'electron-microscopy', 'electron-crystallography']
        ).values_list(
            'pdb_code__index',
            'pdb_data__pdb',
            'preferred_chain'
        )

        for pdb_code, pdb_data, preferred_chain in exp_structures:
            output_filename = f"{output_db_dir}/{pdb_code}_raw_info.pdb"
            print(f"Processing structure {pdb_code} with preferred chain {preferred_chain}")
            process_structure(pdb_code, pdb_data, output_filename, preferred_chain)

    elif structure_type == 'af':
        af_structures = StructureModel.objects.filter(main_template_id__isnull=True).values_list(
            'protein__entry_name',
            'pdb_data__pdb',
            'state__slug'
        )

        for pdb_code, pdb_data, state in af_structures:
            output_filename = f"{output_db_dir}/{pdb_code}_{state}_af_info.pdb"
            print(f"Processing structure {pdb_code}")
            process_structure(pdb_code, pdb_data, output_filename)

    elif structure_type == 'ref':
        inactive_structures = StructureModel.objects.filter(main_template_id__isnull=False).values_list(
            'protein__entry_name',
            'pdb_data__pdb'
        )

        active_structures = Structure.objects.filter(
            structure_type__slug__in=['af-signprot-refined-cem', 'af-signprot-refined-xray']
        ).values_list(
            'pdb_code__index',
            'pdb_data__pdb',
            'preferred_chain'
        )

        # Process active structures
        for pdb_code, pdb_data, preferred_chain in active_structures:
            output_filename = f"{output_db_dir}/{pdb_code}_ref_info.pdb"
            print(f"Processing structure {pdb_code} with preferred chain {preferred_chain}")
            process_structure(pdb_code, pdb_data, output_filename, preferred_chain)

        # Process inactive structures
        for pdb_code, pdb_data in inactive_structures:
            output_filename = f"{output_db_dir}/{pdb_code}_ref_info.pdb"
            print(f"Processing structure {pdb_code}")
            process_structure(pdb_code, pdb_data, output_filename)

class Command(BaseBuild):
    help = "Assigns generic numbers to structures and extracts residues based on B-factors, for Foldseek databases"

    def add_arguments(self, parser):
        parser.add_argument(
            '--structure_type',
            nargs='*',
            choices=['raw', 'ref', 'af'],
            help='Specify which structure types to process. If omitted, all types are processed.',
        )

    def handle(self, *args, **options):
        # Define the output directory
        output_dir = os.path.join(settings.DATA_DIR, 'structure_data')

        # Get the structure types to process
        structure_types = options['structure_type'] or ['raw', 'ref', 'af']

        for structure_type in structure_types:
            print(f"Processing structure type: {structure_type}")
            process_structures(structure_type, output_dir)
