# Semi-automated workflow for docking small molecules
This simple workflow for docking small molecule ligands into receptors is ideal for hit-to-lead or lead optimization. Its core modules are:

1. Protein preparation with PDBFixer, rdkit, MDAnalysis, Biopython
2. Ligand preparation with Scrubber (molscrub)
4. Docking with Gnina (based on Autodock Vina and Smina)
5. Visual analysis with py3Dmol

This workflow can be used by beginner and advanced users alike. The only input needed is the PDB id. The ligand(s) will be automatically identified after entering the PDB id in the ligand preparation module.

Scope and limitations:
1. rigid and flexible docking (later: free energy perturbation)
2. small molecule ligand docking
3. docking into known binding site of co-crystallized ligand
4. docking into unknown binding site(s)

# Requirements/packages needed to be installed
Biopython, Gnina, MDAnalysis, Numpy, OpenBabel, OpenMM, Os, PDBFixer, Requests, Rdkit utils, Scrubber, Subprocess, a list of SMILES of compounds to be prepared for docking as a csv file
NOTE: This workflow was run in WSL (Windows Subsystem for linux).

# Installation
When running the scripts, the packages should be installed automatically since the code is implemented in the script. 
Should installation fail somehow, the packages can be installed manually via pip install or conda install -c conda-forge.
When you are ready, open WSL.

## Get gnina
Type in the following command: 

wget https://github.com/gnina/gnina/releases/download/v1.3/gnina.fix

## Make gnina executable
Type in the following commands:

mv gnina.fix gnina

chmod +x gnina

## Install packages
Type in the following commands: 

pip install biopython mdanalysis numpy openmm os pdbfixer requests useful_rdkit_utils scrubber subprocess

conda install -c conda-forge openbabel 

# Usage
## Protein preparation
Type in the following command: 

python Protein_preparation.py

Packages will be installed; script will prompt you to enter a PDB ID for your target protein. Folders will be created (default is molecular_docking/protein_files). The target protein will be prepared automatically and a pdbqt file will be generated at the end of the script.

### Before running the next script
Create a folder named "molecular_docking/ligand_structures" by typing 

mkdir -p molecular_docking/ligand_structures.

Save a csv file containing SMILES of compounds that you want to prepare for docking into the folder you just created.

## Ligand preparation
type in the following command: 

python Ligand_preparation.py

Packages will be installed; script will prompt you to enter a PDB ID for your target protein. Folders will be created (default is molecular_docking/ligand_structures). The script will identify all ligands bound to the PDB structure and prompts you to select a ligand you are interested in. This ligand will be saved as a ligand_id variable and will stick with the rest of the script. Several ligand files (sdf, pdb) will be generated along the way. The final output is:
1. a ligand that is pose corrected against an ideal ligand (downloaded automatically from https://www.rcsb.org/ and saved as {ligand_id}_ideal.sdf) and prepared by scrubber (aka molscrub) and saved as {ligand_id}.sdf.

2. an sdf file containing several scrubbed and pose corrected compounds obtained from the csv file.

## Docking with Gnina
type in the following command: 

python Docking_with_Gnina.py

Several packages will be instgalled, gnina will be downloaded and installed, and a folder named molecular_docking/docking_results will be created. The script will prompt you to enter the PDB ID and ligand ID used in the previous scripts to set the variables for the rest of the script. The you will prompted to select the docking mode by typing in the desired letter (selections a-d). After some time when the script finishes, the output will be a docked ligand sdf file (name output depends on the docking mode selected earlier). The results can be analyzed in an external program or script.
