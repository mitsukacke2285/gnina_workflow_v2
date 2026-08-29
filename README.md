# Simplified Docking with Gnina (SDoG) - a fully automated Gnina workflow

```
 / \__
(    @\___
 /         O
/   (_____/
/_____/   U

```

A semi-automated workflow for docking small-molecule ligands into protein targets. The workflow is designed to support **hit-to-lead** and **lead-optimization** projects while remaining accessible to both beginner and advanced users.

The workflow takes a **PDB ID** as the primary target input and automates protein preparation, identification and preparation of a co-crystallized ligand, ligand preparation, and molecular docking with **GNINA**.

---

# Workflow Overview

The workflow consists of three main stages:

1. **Protein preparation**
2. **Ligand preparation**
3. **Molecular docking with GNINA**

### 1. Protein preparation

The workflow downloads the selected protein structure from the **RCSB Protein Data Bank** and prepares it for docking using:

* PDBFixer
* RDKit
* MDAnalysis
* Biopython
* OpenMM

The final output is a docking-ready **PDBQT** protein structure.

### 2. Ligand preparation

The workflow identifies small-molecule ligands present in the selected PDB structure and allows the user to select the ligand of interest.

The selected co-crystallized ligand is used as a reference for binding-site and pose preparation.

User-provided compounds are supplied as SMILES in a CSV file and prepared using **Scrubber (molscrub)**.

### 3. Molecular docking

Prepared ligands are docked into the prepared protein using **GNINA**, a docking program based on AutoDock Vina and Smina that also provides CNN-based scoring capabilities.

---

## Workflow

```text
                         PDB ID
                           │
                           ▼
                  Download PDB Structure
                           │
                           ▼
                  Protein Preparation
            PDBFixer / RDKit / MDAnalysis
                    Biopython / OpenMM
                           │
                           ▼
                 Prepared Protein PDBQT
                           │
              ┌────────────┴────────────┐
              │                         │
              ▼                         ▼
      Co-crystallized Ligand       User Compounds
              │                         │
              ▼                         ▼
       Reference Ligand          SMILES CSV Input
         Preparation                    │
              │                         ▼
              │                  Scrubber / molscrub
              │                         │
              └────────────┬────────────┘
                           ▼
                   Prepared Ligands
                           │
                           ▼
                    GNINA Docking
                           │
                           ▼
                    Docked SDF Files
                           │
                           ▼
                  Post-Docking Analysis
```

---

# Scope and Limitations

The current workflow is intended for:

* Small-molecule docking
* Hit-to-lead projects
* Lead-optimization projects
* Rigid docking
* Flexible docking
* Docking into known binding sites
* Docking using a co-crystallized ligand as a reference
* Docking into unknown or predicted binding sites when appropriate docking coordinates are provided

### Planned future development

Potential future extensions include:

* Free-energy perturbation (FEP)
* Binding-site prediction
* Automated docking-score analysis
* Compound ranking and prioritization
* Protein-ligand interaction analysis
* Automated post-docking analysis

> **Important:** Molecular docking provides computational predictions of ligand binding poses and scores. Docking results should not be interpreted as experimental measurements of binding affinity or biological activity. Experimental validation is required.

---

# Requirements

This workflow was developed and tested under **WSL (Windows Subsystem for Linux)**.

## Software requirements

The workflow requires:

* Python 3
* Conda
* GNINA
* Open Babel

The Python dependencies are provided in:

```text
requirements.txt
```

### Python packages

The `requirements.txt` file contains:

* **Biopython** — biological structure and sequence handling
* **MDAnalysis** — molecular structure and trajectory analysis
* **NumPy** — numerical computing
* **RDKit** — cheminformatics and molecular manipulation
* **OpenMM** — molecular simulation and energy minimization
* **PDBFixer** — protein structure preparation and repair
* **Scrubber** — ligand preparation and standardization
* **useful_rdkit_utils** — RDKit utility functions
* **Requests** — downloading files from web resources

### Dependencies not included in `requirements.txt`

Some software used by the workflow is not installed through `pip`.

#### GNINA

GNINA is a standalone executable and must be downloaded separately.

#### Open Babel

Open Babel is a cheminformatics command-line application and is installed separately through Conda.

#### Python standard-library modules

The workflow also uses Python standard-library modules such as:

* `os`
* `subprocess`

These are included with Python and **do not need to be installed separately**.

---

# Installation

## 1. Open WSL or Ubuntu

The workflow was developed for **Windows Subsystem for Linux (WSL)**.

Open your WSL or Ubuntu terminal before proceeding.

A dedicated Conda environment is recommended to prevent dependency conflicts with other Python projects.

---

## 2. Create a Conda environment

Create a new environment:

```bash
conda create -n molecular_docking python=3.11
```

Activate the environment:

```bash
conda activate molecular_docking
```

Verify the Python version:

```bash
python --version
```

---

## 3. Install Python dependencies

Navigate to the directory containing the workflow and `requirements.txt`.

For example:

```bash
cd /path/to/your/workflow
```

Install all Python dependencies using:

```bash
pip install -r requirements.txt
```

This command reads the `requirements.txt` file and automatically installs all packages listed in it.

You can verify the installation with:

```bash
pip list
```

### If installation fails

Some packages may occasionally fail to install through `pip` because of Python versions, operating-system differences, or dependency conflicts.

If this happens, the affected package can generally be installed through Conda:

```bash
conda install -c conda-forge <package-name>
```

---

# 4. Install Open Babel

Open Babel is **not included in `requirements.txt`** because it is installed separately from the Python dependencies.

Install it using Conda:

```bash
conda install -c conda-forge openbabel
```

Verify the installation:

```bash
obabel -V
```

A version number should be displayed if the installation was successful.

---

# 5. Download GNINA

GNINA is also **not included in `requirements.txt`** because it is a standalone executable rather than a Python package.

Download GNINA:

```bash
wget https://github.com/gnina/gnina/releases/download/v1.3/gnina.fix
```

Rename the executable:

```bash
mv gnina.fix gnina
```

Make GNINA executable:

```bash
chmod +x gnina
```

Test GNINA:

```bash
./gnina --help
```

If the GNINA help information is displayed, the executable is working correctly.

### Optional: Add GNINA to your PATH

If you want to run GNINA from any directory, move it into `/usr/local/bin`:

```bash
sudo mv gnina /usr/local/bin/gnina
```

Then verify:

```bash
gnina --help
```

---

# Complete Installation

For a fresh WSL/Conda environment, the complete installation can be performed with:

```bash
# Create Conda environment
conda create -n molecular_docking python=3.11

# Activate environment
conda activate molecular_docking

# Navigate to workflow
cd /path/to/your/workflow

# Install Python dependencies
pip install -r requirements.txt

# Install Open Babel
conda install -c conda-forge openbabel

# Download GNINA
wget https://github.com/gnina/gnina/releases/download/v1.3/gnina.fix

# Rename GNINA
mv gnina.fix gnina

# Make GNINA executable
chmod +x gnina

# Test GNINA
./gnina --help

# Test Open Babel
obabel -V
```

Once these steps have completed successfully, the workflow is ready to run.

---

# Installation Verification

You can test several of the major Python dependencies individually:

```bash
python -c "import Bio; print('Biopython OK')"
python -c "import MDAnalysis; print('MDAnalysis OK')"
python -c "import numpy; print('NumPy OK')"
python -c "from rdkit import Chem; print('RDKit OK')"
python -c "import openmm; print('OpenMM OK')"
python -c "import pdbfixer; print('PDBFixer OK')"
```

Test Open Babel:

```bash
obabel -V
```

Test GNINA:

```bash
./gnina --help
```

If these commands execute successfully, the core software environment is installed.

---

# Usage

The workflow is divided into three scripts:

```text
protein_preparation.py
ligand_extraction_and_preparation.py
in_silico_screening_and_reporting.py
```

The scripts should generally be executed in this order:

```text
Protein Preparation
        ↓
Ligand Preparation
        ↓
GNINA Docking
```

---

# 1. Protein Preparation

Run:

```bash
python Protein_preparation.py
```

The script will prompt you to enter the **PDB ID** of your target protein.

The workflow will then:

1. Download the selected PDB structure.
2. Create the required output directories.
3. Identify and prepare the protein structure.
4. Perform the required structure cleanup and preparation.
5. Generate a docking-ready PDBQT file.

By default, protein-related files are stored under:

```text
molecular_docking/
└── protein_files/
```

The prepared protein PDBQT file generated during this step will be used during docking.

---

# 2. Prepare the Ligand Input File

Before running the ligand-preparation script, create the ligand directory:

```bash
mkdir -p molecular_docking/ligand_structures
```

Place your compound CSV file inside this directory.

The CSV file should contain the **SMILES strings of the compounds that you want to dock**.

For example:

```text
SMILES
CCOc1ccc2nc(S(N)(=O)=O)sc2c1
CC1=CC=C(C=C1)C(=O)NCC2CC2
```

The exact column name and CSV format should match the requirements of `Ligand_preparation.py`.

---

# 3. Ligand Preparation

Run:

```bash
python Ligand_preparation.py
```

The script will prompt you to enter the **PDB ID** of the target protein.

The workflow will:

1. Identify ligands present in the selected PDB structure.
2. Display the available ligands.
3. Prompt you to select the ligand of interest.
4. Store the selected ligand as the `ligand_id` variable.
5. Retrieve the corresponding ideal ligand structure from the RCSB PDB.
6. Correct the experimentally observed ligand pose against the ideal ligand.
7. Prepare the reference ligand using Scrubber (`molscrub`).
8. Read the user-provided compound CSV.
9. Scrub and prepare the user-provided compounds.
10. Generate SDF files suitable for docking.

---

## Reference Ligand

The selected co-crystallized ligand is retrieved from the RCSB PDB and saved as:

```text
{ligand_id}_ideal.sdf
```

The prepared reference ligand is subsequently saved as:

```text
{ligand_id}.sdf
```

This ligand can be used as the reference structure for the docking workflow.

---

## Prepared Compound Library

The compounds supplied in the input CSV file are processed using Scrubber and prepared for docking.

The resulting SDF file contains the prepared compounds and is stored under:

```text
molecular_docking/
└── ligand_structures/
```

Several intermediate files may also be generated during ligand preparation.

---

# 4. Docking with GNINA

Once protein and ligand preparation are complete, run:

```bash
python Docking_with_Gnina.py
```

The script will:

1. Check the required dependencies.
2. Locate or download GNINA as required by the implementation.
3. Create the docking-results directory.
4. Prompt you for the **PDB ID** used in the previous steps.
5. Prompt you for the **ligand ID** selected during ligand preparation.
6. Set the required docking variables.
7. Prompt you to select a docking mode.
8. Run GNINA using the prepared protein and ligand structures.

The available docking modes are selected using the provided options (`a–d`).

---

# Docking Output

The docking results are stored under:

```text
molecular_docking/
└── docking_results/
```

The exact output filename depends on the docking mode selected.

The final output is a **docked SDF file** containing the predicted ligand poses and associated docking information.

The results can be analyzed using external molecular-visualization or cheminformatics software.

---

# Directory Structure

After running the workflow, the project may have a structure similar to:

```text
project/
│
├── Protein_preparation.py
├── Ligand_preparation.py
├── Docking_with_Gnina.py
├── requirements.txt
│
└── molecular_docking/
    │
    ├── protein_files/
    │   ├── original_structure.pdb
    │   └── prepared_protein.pdbqt
    │
    ├── ligand_structures/
    │   ├── input_compounds.csv
    │   ├── {ligand_id}_ideal.sdf
    │   ├── {ligand_id}.sdf
    │   └── prepared_compounds.sdf
    │
    └── docking_results/
        └── docked_compounds.sdf
```

The exact filenames and number of intermediate files may vary depending on the selected PDB structure, ligand, and docking mode.

---

# Input and Output Summary

| Stage                 | Input                               | Output                                            |
| --------------------- | ----------------------------------- | ------------------------------------------------- |
| Protein preparation   | PDB ID                              | Prepared protein PDBQT                            |
| Ligand preparation    | PDB ID + ligand ID + compound CSV   | Prepared reference ligand + prepared compound SDF |
| GNINA docking         | Protein PDBQT + prepared ligand SDF | Docked ligand SDF                                 |
| Post-docking analysis | Docked SDF                          | Docking scores, poses, and interaction analysis   |

---

# Typical Workflow

Once the software environment has been installed, a typical workflow is:

### Step 1 — Prepare the protein

```bash
python Protein_preparation.py
```

Enter the desired PDB ID when prompted.

### Step 2 — Prepare the ligands

Make sure your compound CSV is located in:

```text
molecular_docking/ligand_structures/
```

Then run:

```bash
python Ligand_preparation.py
```

Enter the same PDB ID and select the desired co-crystallized ligand.

### Step 3 — Run docking

```bash
python Docking_with_Gnina.py
```

Enter the PDB ID and ligand ID when prompted, then select the desired docking mode.

### Overall workflow

```text
PDB ID
  ↓
Protein preparation
  ↓
Prepared protein PDBQT
  ↓
Ligand identification
  ↓
Reference ligand preparation
  +
User compound preparation
  ↓
Prepared ligand SDF
  ↓
GNINA
  ↓
Docked compound SDF
  ↓
Post-docking analysis
```

---

# Notes and Best Practices

## Co-crystallized ligands

The workflow is particularly useful for PDB structures containing a known co-crystallized small-molecule ligand.

The selected ligand provides an experimentally observed reference for the binding site and ligand pose.

## Unknown binding sites

Docking into an unknown binding site is possible, but the appropriate docking coordinates must be supplied by the user or obtained from an upstream binding-site prediction method.

The workflow itself does not establish whether a predicted binding site is biologically valid.

## Protein preparation

Protein preparation is an important part of the docking workflow. The quality of the input structure, protonation states, missing residues, alternate conformations, cofactors, waters, and other structural features can influence docking results.

Users should therefore inspect the prepared protein structure when applying the workflow to a new target.

## Docking parameters

Docking results can depend strongly on the selected GNINA parameters, including:

* Search space
* Exhaustiveness
* Number of poses
* Scoring function
* Receptor flexibility
* Ligand flexibility
* CPU/GPU availability

For reproducible experiments, these parameters should be recorded alongside the docking results.

---

# Reproducibility

For reproducible computational experiments, it is recommended to record:

* PDB ID
* Selected ligand ID
* Input compound CSV
* Python version
* Conda environment
* Package versions
* GNINA version
* Open Babel version
* Docking mode
* Docking parameters
* Protein preparation settings

Keeping the original input and intermediate structures is also recommended.

---

# Troubleshooting

## `ModuleNotFoundError`

If Python reports:

```text
ModuleNotFoundError: No module named '...'
```

make sure the Conda environment is activated:

```bash
conda activate molecular_docking
```

Then reinstall the dependencies:

```bash
pip install -r requirements.txt
```

---

## GNINA cannot be executed

If you receive an error such as:

```text
Permission denied
```

make GNINA executable:

```bash
chmod +x gnina
```

Then test:

```bash
./gnina --help
```

If `gnina` cannot be found after moving it to `/usr/local/bin`, check that the executable is present:

```bash
which gnina
```

---

## Open Babel cannot be found

Test:

```bash
obabel -V
```

If the command is not found, install Open Babel:

```bash
conda install -c conda-forge openbabel
```

---

# Disclaimer

This workflow is intended as a **computational research and molecular-modeling tool**.

Docking scores and predicted binding poses are computational hypotheses and should not be considered experimental measurements of binding affinity, potency, selectivity, or biological activity.

Docking results should therefore be interpreted together with appropriate structural, biochemical, pharmacological, and experimental data.

---

# Future Development

Potential future additions to the workflow include:

* Automated binding-site detection
* Improved protein preparation
* Multiple receptor conformations
* Ensemble docking
* Automated docking-score ranking
* Ligand interaction analysis
* Pose clustering
* Rescoring
* Molecular dynamics simulations
* Free-energy perturbation (FEP)
* Automated reporting and visualization

The long-term goal is to develop the workflow into a more comprehensive **semi-automated computational drug-discovery platform** for small-molecule hit identification and lead optimization.
