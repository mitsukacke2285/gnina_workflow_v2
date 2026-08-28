# config.py

# ============================================================
# PROJECT INFORMATION
# ============================================================

PROJECT_NAME = "GNINA docking workflow"
VERSION = "0.1.0"


# ============================================================
# FILE PATHS
# ============================================================

PROTEIN_DIRECTORY = "molecular_docking/protein_files"
LIGAND_DIRECTORY = "molecular_docking/ligand_structures"
PDB_ID =  input(
        "Enter desired PDB code to be downloaded from RCSB: "
    ).strip().upper()
MOLECULAR_DOCKING_DIRECTORY = "molecular_docking/docking_results"