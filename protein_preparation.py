#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/env python
# coding: utf-8

# ## Protein preparation

# Please install the following packages:
# 
# - PDBFixer
# - Biopython
# - OpenMM (and OpenMMForceFields)
# - OpenBabel

# In[ ]:


import os
import config
from Bio.PDB import PDBParser, Select, PDBIO
from pdbfixer import PDBFixer
from openmm.app import PDBFile, ForceField, Simulation
from openmm import VerletIntegrator
import openmm.unit as unit
import openbabel.pybel as pybel

def select_chain(pdb_file):
    """Extract individual protein chains and return chain A."""

    print("\n=== Selecting chain A if protein contains multiple chains ===")

    class ChainSelector(Select):

        def __init__(self, target_chain):
            self.target_chain = target_chain

        def accept_chain(self, chain):
            return chain.id == self.target_chain

    # Load structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    # Save each chain
    io = PDBIO()

    chain_a_file = None

    for model in structure:
        for chain in model:

            chain_id = chain.id

            output_file = os.path.join(
                config.PROTEIN_DIRECTORY,
                f"{config.PDB_ID}_{chain_id}.pdb"
            )

            io.set_structure(structure)
            io.save(
                output_file,
                ChainSelector(chain_id)
            )

            print(
                f"\n=== {config.PDB_ID}_{chain_id}.pdb "
                f"has been extracted. ==="
            )

            if chain_id == "A":
                chain_a_file = output_file

    if chain_a_file is None:
        raise ValueError("Chain A was not found in the PDB structure.")

    print(
        "\n=== Chain A selected for further processing! ==="
    )

    return chain_a_file
    
def fix_protein(pdb_file):
    """Fix missing residues/atoms and minimize the protein."""

    print("\n=== Starting PDBFixer... ===")

    # Create PDBFixer using the input file
    fixer = PDBFixer(filename=pdb_file)

    # Find and fix missing/nonstandard residues
    fixer.findMissingResidues()

    fixer.findNonstandardResidues()
    print("Nonstandard residues:")
    print(fixer.nonstandardResidues)

    fixer.replaceNonstandardResidues()

    # Remove heterogens but keep protein
    fixer.removeHeterogens(keepWater=False)

    # Find and add missing atoms
    fixer.findMissingAtoms()

    print("Missing atoms:")
    print(fixer.missingAtoms)

    print("Missing terminal atoms:")
    print(fixer.missingTerminals)

    fixer.addMissingAtoms()

    # Add hydrogens at physiological pH
    fixer.addMissingHydrogens(pH=7.4)

    print(
        "\n=== Loading force field "
        "(amber14-all.xml, amber14/tip3p.xml)... ==="
    )

    forcefield = ForceField(
        "amber14-all.xml",
        "amber14/tip3p.xml"
    )

    # Create OpenMM system
    system = forcefield.createSystem(
        fixer.topology,
        ignoreExternalBonds=False
    )

    print("\n=== Force field loaded ===")

    # Create simulation
    print("\n=== Creating simulation for minimization ===")

    integrator = VerletIntegrator(
        0.001 * unit.picoseconds
    )

    platform = None

    simulation = Simulation(
        fixer.topology,
        system,
        integrator,
        platform
    )

    simulation.context.setPositions(
        fixer.positions
    )

    print("\n=== Minimizing energy... ===")

    simulation.minimizeEnergy()

    minimized_positions = (
        simulation.context
        .getState(getPositions=True)
        .getPositions()
    )

    return fixer, minimized_positions


def save_minimized_structure(fixer, minimized_positions):
    """Save minimized protein and generate PDBQT."""

    fixed_pdb = os.path.join(
        config.PROTEIN_DIRECTORY,
        f"{config.PDB_ID}_A_fixed.pdb"
    )

    with open(fixed_pdb, "w") as output:

        PDBFile.writeFile(
            fixer.topology,
            minimized_positions,
            output
        )

    print(
        f"\n=== Minimization complete. "
        f"Minimized structure saved to {fixed_pdb} ==="
    )

    # --------------------------------------------------------
    # Generate PDBQT
    # --------------------------------------------------------

    print("\n=== Generating PDBQT file ===")

    receptor_pdbqt_path = os.path.join(
        config.PROTEIN_DIRECTORY,
        f"{config.PDB_ID}_A.pdbqt"
    )

    if os.path.exists(receptor_pdbqt_path):

        print(
            f"\n=== {receptor_pdbqt_path} already exists. "
            f"Skipping PDBQT generation. ==="
        )

    else:

        print(
            f"\n=== Generating {receptor_pdbqt_path}... ==="
        )

        mol = next(
            pybel.readfile(
                "pdb",
                fixed_pdb
            )
        )

        mol.write(
            "pdbqt",
            receptor_pdbqt_path,
            overwrite=True
        )

        print(
            f"\n=== {config.PDB_ID}_A.pdbqt "
            f"has been generated and saved ==="
        )

    print("\n=== Fixing protein complete ===")

