#!/usr/bin/env python
# coding: utf-8

# In[2]:


import os
import urllib.request
import config


def download_pdb_file(pdb_id):

    # Ensure output directory exists
    os.makedirs(config.PROTEIN_DIRECTORY, exist_ok=True)

    # Define output file
    pdb_file = os.path.join(
        config.PROTEIN_DIRECTORY,
        f"{pdb_id}.pdb"
    )

    # RCSB download URL
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"

    try:
        urllib.request.urlretrieve(url, pdb_file)
        print(f"Downloaded: {pdb_file}")
        return pdb_file

    except Exception as e:
        print(f"Download error: {e}")
        return None

