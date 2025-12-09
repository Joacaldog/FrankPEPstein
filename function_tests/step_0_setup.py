#@title 0.2 Install Dependencies & Setup Tools
#@markdown This cell clones the repository and installs dependencies in the base environment.

import os
import sys
import subprocess

# --- 1. Clone Repository ---
if not os.path.exists("FrankPEPstein"):
    print("Cloning repository...")
    !git clone https://github.com/Joacaldog/FrankPEPstein.git
else:
    print("Repository already exists.")

# --- 2. Install Dependencies (Base Environment) ---
print("Installing dependencies in base environment (matches Colab kernel)...")
# Using 'mamba install' installs into the active environment (base)
# This prevents Bio.PDB.ccealign errors due to binary incompatibility between envs
!mamba install -q -y -c conda-forge -c salilab openbabel biopython fpocket joblib tqdm py3dmol vina python=3.10 salilab::modeller

# --- 3. Run Notebook Setup Utils ---
print(f"Dependencies installed.")

repo_path = os.path.abspath("FrankPEPstein")
if repo_path not in sys.path:
    sys.path.append(repo_path)
from scripts import notebook_utils

# DRIVE CONFIGURATION: Enter your File IDs here
drive_ids = {
    "adfr_id": "1gmRj8mva84-JB7UXUcQfB3Ziw_nwwdox",       # ADFRsuite_x86_64Linux_1.0.tar.gz
    "db_id": "1a4GoZ1ZT-DNYMyvVtKJukNdF6TAaLJU5",    # minipockets_..._curated-db.tar.gz
    "dict_id": "1nrwSUof0lox9fp8Ow5EICIN9u0lglu7U"      # reduce_wwPDB_het_dict.tar.gz
}

print("Setting up external tools...")
notebook_utils.setup_external_tools(drive_ids)

print("Configuring Modeller...")
notebook_utils.configure_modeller()

print("Setup Complete!")
