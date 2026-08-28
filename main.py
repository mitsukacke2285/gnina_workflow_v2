#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import sys

if "download_pdb" in sys.modules:
    del sys.modules["download_pdb"]

import config

from download_pdb import download_pdb_file
from protein_preparation import select_chain, fix_protein, save_minimized_structure
from ligand_extraction_and_preparation import select_ligand_from_pdb, download_ideal_ligand, fix_and_align, scrubbing_ligands
from in_silico_screening_and_reporting import gnina 
# from compound_profiling import create_admet_model, profile_file
# from filter import combined_filter, select_compounds
# from file_io import load_csv, save_csv, save_sdf


pdb_id = config.PDB_ID
protein_directory = config.PROTEIN_DIRECTORY
ligand_directory = config.LIGAND_DIRECTORY


def main():

    # ========================================================
    # 1. Download PDB
    # ========================================================

    print("Download PDB file...")  

    download_pdb_file(
        pdb_id
    )
    
    print(
        f"{pdb_id}.pdb has been saved to: "
        f"{protein_directory}"
    )

    # ========================================================
    # 2. Protein preparation
    # ========================================================

    select_chain(f"{protein_directory}/{pdb_id}.pdb")
    fixer, minimized_positions = fix_protein(f"{protein_directory}/{pdb_id}_A.pdb")
    save_minimized_structure(fixer, minimized_positions)

    # ========================================================
    # 3. Ligand extraction and preparation
    # ========================================================

    ligand_id, pose_mol = select_ligand_from_pdb()

    ideal_ligand, ideal_mol = download_ideal_ligand(ligand_id)

    corrected_pose_with_H = fix_and_align(
        ideal_mol,
        pose_mol,
        ligand_id
    )

    scrubbing_ligands(ligand_id)


    # ========================================================
    # 4. In silico screening
    # ========================================================

    gnina(ligand_id)
    
if __name__ == "__main__":
    main()


# In[ ]:


# In[ ]:




