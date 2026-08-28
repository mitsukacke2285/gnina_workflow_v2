#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/env python
# coding: utf-8

"""
Gnina Molecular Docking Workflow

Docking modes:
    a = Single-ligand redocking + RMSD calculation
    b = Batch docking
    c = Flexible docking
    d = Unknown-site docking

The script:
    1. Checks for Gnina and downloads it if necessary.
    2. Reads PDB and ligand IDs from params.json.
    3. Runs the selected docking protocol.
    4. Extracts docking scores into a CSV report.
    5. Calculates RMSD for single-ligand redocking.
"""

import json
import os
import subprocess

import pandas as pd
import useful_rdkit_utils as uru

from rdkit import Chem, RDLogger
from rdkit.Chem import PandasTools
from rdkit.rdBase import BlockLogs

import config


# ============================================================================
# 1. Gnina installation
# ============================================================================

GNINA_FILE = "gnina"

GNINA_URL = (
    "https://github.com/gnina/gnina/releases/download/v1.3/gnina"
)


def install_gnina():
    """Download Gnina if it is not already available."""

    if os.path.exists(GNINA_FILE):
        print(f"✅ {GNINA_FILE} already exists, skipping download.")
        return

    print(f"⬇️ {GNINA_FILE} not found, downloading...")

    subprocess.run(
        [
            "wget",
            GNINA_URL,
            "-O",
            GNINA_FILE,
        ],
        check=True,
    )

    subprocess.run(
        ["chmod", "+x", GNINA_FILE],
        check=True,
    )

    print("✅ Gnina download complete.")


# ============================================================================
# 2. Working directories and docking parameters
# ============================================================================

protein_directory = config.PROTEIN_DIRECTORY
ligand_directory = config.LIGAND_DIRECTORY
docking_results_directory = config.MOLECULAR_DOCKING_DIRECTORY

os.makedirs(docking_results_directory, exist_ok=True)

pdb_id = config.PDB_ID

# ============================================================================
# 3. Docking modes
# ============================================================================

def single_ligand_docking(
    pdb_id,
    ligand_id,
    exhaustiveness,
    gpu_flag,
    cnn_flag,
):
    """
    Perform single-ligand docking / redocking.
    """

    receptor = os.path.join(
        protein_directory,
        f"{pdb_id}_A.pdbqt",
    )

    ligand = os.path.join(
        ligand_directory,
        f"{ligand_id}.sdf",
    )

    reference_ligand = os.path.join(
        ligand_directory,
        f"{ligand_id}_corrected_pose.sdf",
    )

    output = os.path.join(
        docking_results_directory,
        f"{ligand_id}_docked_{pdb_id}.sdf",
    )

    command = [
        f"./{GNINA_FILE}",
        "-r",
        receptor,
        "-l",
        ligand,
        "--autobox_ligand",
        reference_ligand,
        "-o",
        output,
        "--seed",
        "0",
        "--exhaustiveness",
        str(exhaustiveness),
        *gpu_flag,
        *cnn_flag,
    ]

    print("\nRunning Gnina:")
    print(" ".join(command))

    subprocess.run(command, check=True)

    return output


def batch_docking(
    pdb_id,
    ligand_id,
    exhaustiveness,
    gpu_flag,
    cnn_flag,
):
    """
    Perform batch docking of multiple ligands.
    """

    receptor = os.path.join(
        protein_directory,
        f"{pdb_id}_A.pdbqt",
    )

    ligand_file = os.path.join(
        ligand_directory,
        "ligands_to_dock.sdf",
    )

    reference_ligand = os.path.join(
        ligand_directory,
        f"{ligand_id}_corrected_pose.sdf",
    )

    output = os.path.join(
        docking_results_directory,
        f"multiple_ligands_docked_{pdb_id}.sdf",
    )

    command = [
        f"./{GNINA_FILE}",
        "-r",
        receptor,
        "-l",
        ligand_file,
        "--autobox_ligand",
        reference_ligand,
        "-o",
        output,
        "--seed",
        "0",
        "--exhaustiveness",
        str(exhaustiveness),
        *gpu_flag,
        *cnn_flag,
    ]

    print("\nRunning Gnina:")
    print(" ".join(command))

    subprocess.run(command, check=True)

    return output


def flexible_docking(
    pdb_id,
    ligand_id,
    gpu_flag,
    cnn_flag,
):
    """
    Perform flexible docking using maximum exhaustiveness (64).
    """

    receptor = os.path.join(
        protein_directory,
        f"{pdb_id}_A.pdbqt",
    )

    ligand = os.path.join(
        ligand_directory,
        f"{ligand_id}.sdf",
    )

    reference_ligand = os.path.join(
        ligand_directory,
        f"{ligand_id}_corrected_pose.sdf",
    )

    output = os.path.join(
        docking_results_directory,
        f"{ligand_id}_flex_{pdb_id}.sdf",
    )

    command = [
        f"./{GNINA_FILE}",
        "-r",
        receptor,
        "-l",
        ligand,
        "--autobox_ligand",
        reference_ligand,
        "-o",
        output,
        "--flexdist_ligand",
        reference_ligand,
        "--flexdist",
        "3.59",
        "--seed",
        "0",
        "--exhaustiveness",
        "64",
        *gpu_flag,
        *cnn_flag,
    ]

    print("\nRunning Gnina:")
    print(" ".join(command))

    subprocess.run(command, check=True)

    return output


def unknown_site_docking(
    pdb_id,
    ligand_id,
    exhaustiveness,
    gpu_flag,
    cnn_flag,
):
    """
    Perform docking over the whole protein.
    """

    receptor = os.path.join(
        protein_directory,
        f"{pdb_id}_A.pdbqt",
    )

    ligand = os.path.join(
        ligand_directory,
        f"{ligand_id}.sdf",
    )

    output = os.path.join(
        docking_results_directory,
        f"{ligand_id}_docked_whole_{pdb_id}.sdf",
    )

    command = [
        f"./{GNINA_FILE}",
        "-r",
        receptor,
        "-l",
        ligand,
        "--autobox_ligand",
        receptor,
        "-o",
        output,
        "--seed",
        "0",
        "--exhaustiveness",
        str(exhaustiveness),
        *gpu_flag,
        *cnn_flag,
    ]

    print("\nRunning Gnina:")
    print(" ".join(command))

    subprocess.run(command, check=True)

    return output


# ============================================================================
# 4. Docking report
# ============================================================================

def report(docking_results, pdb_id):
    """
    Convert Gnina SDF output into a CSV report.

    Parameters
    ----------
    docking_results : str or list[str]
        Gnina SDF output file(s).

    pdb_id : str
        Protein/PDB identifier.

    Returns
    -------
    str
        Path to the generated CSV file.
    """

    score_columns = [
        "minimizedAffinity",
        "CNNscore",
        "CNNaffinity",
        "CNN_VS",
        "CNNaffinity_variance",
    ]

    # Normalize input
    if isinstance(docking_results, str):
        sdf_paths = [docking_results]
    else:
        sdf_paths = docking_results

    df_list = []

    for filename in sdf_paths:

        if not os.path.exists(filename):
            raise FileNotFoundError(
                f"Docking result file not found: {filename}"
            )

        print(f"Reading docking results: {filename}")

        with BlockLogs():
            df = PandasTools.LoadSDF(
                filename,
                molColName="ROMol",
            )

        if df.empty:
            print(f"⚠️ No molecules found in {filename}")
            continue

        df_list.append(df)

    if not df_list:
        raise ValueError(
            "No valid docking results were found."
        )

    combo_df = pd.concat(
        df_list,
        ignore_index=True,
    )

    # ------------------------------------------------------------------------
    # Convert Gnina scores to numeric
    # ------------------------------------------------------------------------

    for col in score_columns:

        if col in combo_df.columns:

            combo_df[col] = pd.to_numeric(
                combo_df[col],
                errors="coerce",
            )

            # Remove completely empty columns
            if combo_df[col].isna().all():
                combo_df.drop(
                    columns=[col],
                    inplace=True,
                )

    # ------------------------------------------------------------------------
    # Extract SMILES
    # ------------------------------------------------------------------------

    if "ROMol" in combo_df.columns:

        combo_df["SMILES"] = combo_df["ROMol"].apply(
            lambda mol: (
                Chem.MolToSmiles(mol)
                if mol is not None
                else None
            )
        )

    # ------------------------------------------------------------------------
    # Remove RDKit molecule column from CSV
    # ------------------------------------------------------------------------

    if "ROMol" in combo_df.columns:
        combo_df.drop(
            columns=["ROMol"],
            inplace=True,
        )

    # ------------------------------------------------------------------------
    # Sort by docking affinity
    # ------------------------------------------------------------------------

    if "minimizedAffinity" in combo_df.columns:

        combo_df.sort_values(
            by="minimizedAffinity",
            ascending=True,
            inplace=True,
        )

        combo_df.reset_index(
            drop=True,
            inplace=True,
        )

        combo_df.insert(
            0,
            "Rank",
            range(1, len(combo_df) + 1),
        )

    # ------------------------------------------------------------------------
    # Save CSV
    # ------------------------------------------------------------------------

    output_csv = os.path.join(
        docking_results_directory,
        f"docking_results_{pdb_id}.csv",
    )

    combo_df.to_csv(
        output_csv,
        index=False,
    )

    print(
        f"\n✅ Docking results saved to:\n"
        f"{output_csv}"
    )

    # Display first five results
    print("\nTop docking results:")

    print(
        combo_df.head().to_string(
            index=False
        )
    )

    return output_csv


# ============================================================================
# 5. RMSD calculation
# ============================================================================

def rmsd_calculation(
    ligand_directory,
    docking_results_directory,
    ligand_id,
    pdb_id,
):
    """
    Calculate MCS RMSD between the cognate ligand and docked poses.
    """

    print(
        "\n=== Running MCS RMSD calculation ==="
    )

    cognate_path = os.path.join(
        ligand_directory,
        f"{ligand_id}_corrected_pose.sdf",
    )

    docking_path = os.path.join(
        docking_results_directory,
        f"{ligand_id}_docked_{pdb_id}.sdf",
    )

    if not os.path.exists(cognate_path):
        raise FileNotFoundError(
            f"Cognate ligand not found: {cognate_path}"
        )

    if not os.path.exists(docking_path):
        raise FileNotFoundError(
            f"Docking result not found: {docking_path}"
        )

    cognate = Chem.MolFromMolFile(
        cognate_path
    )

    if cognate is None:
        raise ValueError(
            f"Could not read cognate ligand: {cognate_path}"
        )

    poses = Chem.SDMolSupplier(
        docking_path
    )

    log_path = os.path.join(
        docking_results_directory,
        f"{ligand_id}_{pdb_id}_rmsd.log",
    )

    with open(
        log_path,
        "w",
    ) as log_file:

        log_file.write(
            "Pose_Index\tNum_Matches\tRMSD\n"
        )

        for i, pose in enumerate(poses):

            if pose is None:
                continue

            RDLogger.DisableLog(
                "rdApp.warning"
            )

            n_match, rmsd = uru.mcs_rmsd(
                cognate,
                pose,
            )

            line = (
                f"{i}\t"
                f"{n_match}\t"
                f"{rmsd:.2f}\n"
            )

            print(line.strip())

            log_file.write(line)

    print(
        f"\n=== RMSD results saved to ===\n"
        f"{log_path}"
    )

    return log_path


# ============================================================================
# 6. Docking workflow
# ============================================================================

def docking_main(
    selection,
    pdb_id,
    ligand_id,
    exhaustiveness,
    gpu_flag,
    cnn_flag,
):
    """
    Run the selected docking protocol.
    """

    # ------------------------------------------------------------------------
    # Single ligand / redocking
    # ------------------------------------------------------------------------

    if selection == "a":

        print(
            "\n=== Single ligand docking ==="
        )

        output = single_ligand_docking(
            pdb_id,
            ligand_id,
            exhaustiveness,
            gpu_flag,
            cnn_flag,
        )

        rmsd_calculation(
            ligand_directory,
            docking_results_directory,
            ligand_id,
            pdb_id,
        )

        report(
            output,
            pdb_id,
        )

    # ------------------------------------------------------------------------
    # Batch docking
    # ------------------------------------------------------------------------

    elif selection == "b":

        print(
            "\n=== Batch docking ==="
        )

        output = batch_docking(
            pdb_id,
            ligand_id,
            exhaustiveness,
            gpu_flag,
            cnn_flag,
        )

        report(
            output,
            pdb_id,
        )

    # ------------------------------------------------------------------------
    # Flexible docking
    # ------------------------------------------------------------------------

    elif selection == "c":

        print(
            "\n=== Flexible docking ==="
        )

        output = flexible_docking(
            pdb_id,
            ligand_id,
            gpu_flag,
            cnn_flag,
        )

        report(
            output,
            pdb_id,
        )

    # ------------------------------------------------------------------------
    # Unknown-site docking
    # ------------------------------------------------------------------------

    elif selection == "d":

        print(
            "\n=== Unknown-site docking ==="
        )

        output = unknown_site_docking(
            pdb_id,
            ligand_id,
            exhaustiveness,
            gpu_flag,
            cnn_flag,
        )

        report(
            output,
            pdb_id,
        )

    else:

        raise ValueError(
            "Invalid docking selection."
        )


# ============================================================================
# 7. User input
# ============================================================================

def get_user_settings():
    """
    Collect docking parameters from the user.
    """

    print(
        """
=== Welcome to Gnina ===

Select docking mode:

    a = Single ligand docking
        RMSD and CNN scores will be calculated

    b = Batch docking

    c = Flexible docking
        Exhaustiveness fixed at 64

    d = Docking on unknown site
"""
    )

    selection = input(
        "Select docking mode [a/b/c/d]: "
    ).strip().lower()

    if selection not in {"a", "b", "c", "d"}:
        raise ValueError(
            "Invalid docking mode. "
            "Please select a, b, c, or d."
        )

    # ------------------------------------------------------------------------
    # Exhaustiveness
    # ------------------------------------------------------------------------

    allowed_exhaustiveness = [
        8,
        16,
        24,
        32,
        40,
        48,
        56,
        64,
    ]

    # Flexible docking always uses 64
    if selection == "c":

        exhaustiveness = 64

        print(
            "\nFlexible docking selected."
        )

        print(
            "Exhaustiveness automatically set to 64."
        )

    else:

        print(
            "\nAllowed exhaustiveness levels:"
        )

        print(
            ", ".join(
                map(
                    str,
                    allowed_exhaustiveness,
                )
            )
        )

        try:
            exhaustiveness = int(
                input(
                    "\nDefine exhaustiveness: "
                )
            )

        except ValueError:
            raise ValueError(
                "Exhaustiveness must be an integer."
            )

        if exhaustiveness not in allowed_exhaustiveness:
            raise ValueError(
                "Invalid exhaustiveness level. "
                "Please select one of: "
                + ", ".join(
                    map(
                        str,
                        allowed_exhaustiveness,
                    )
                )
            )

    # ------------------------------------------------------------------------
    # GPU
    # ------------------------------------------------------------------------

    ask_gpu = input(
        "\nRun with GPU? [y/n]: "
    ).strip().lower()

    if ask_gpu == "y":

        gpu_flag = []

    elif ask_gpu == "n":

        gpu_flag = ["--no_gpu"]

    else:

        raise ValueError(
            "Invalid GPU selection. "
            "Please enter y or n."
        )

    # ------------------------------------------------------------------------
    # CNN
    # ------------------------------------------------------------------------

    ask_cnn = input(
        "\nRun with CNN scoring? [y/n]: "
    ).strip().lower()

    if ask_cnn == "y":

        cnn_flag = []

    elif ask_cnn == "n":

        cnn_flag = [
            "--cnn_scoring=none"
        ]

    else:

        raise ValueError(
            "Invalid CNN selection. "
            "Please enter y or n."
        )

    return (
        selection,
        exhaustiveness,
        gpu_flag,
        cnn_flag,
    )


# ============================================================================
# 8. Main
# ============================================================================

def gnina(ligand_id):

    # ------------------------------------------------------------------------
    # Install/check Gnina
    # ------------------------------------------------------------------------

    install_gnina()

    # ------------------------------------------------------------------------
    # User settings
    # ------------------------------------------------------------------------

    (
        selection,
        exhaustiveness,
        gpu_flag,
        cnn_flag,
    ) = get_user_settings()

    print(f"\nSelected docking mode: {selection}")
    print(f"Exhaustiveness: {exhaustiveness}")
    print(
        f"GPU: {'enabled' if not gpu_flag else 'disabled'}"
    )
    print(
        f"CNN scoring: "
        f"{'enabled' if not cnn_flag else 'disabled'}"
    )

    # ------------------------------------------------------------------------
    # Run docking
    # ------------------------------------------------------------------------

    docking_main(
        selection=selection,
        pdb_id=pdb_id,
        ligand_id=ligand_id,
        exhaustiveness=exhaustiveness,
        gpu_flag=gpu_flag,
        cnn_flag=cnn_flag,
    )
