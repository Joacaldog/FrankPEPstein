#@title 0.2 Install Dependencies & Setup Tools
#@markdown This cell clones the repository and installs dependencies in a dedicated environment.

import os
import sys
import subprocess

# --- 1. Clone Repository ---
if not os.path.exists("FrankPEPstein"):
    print("Cloning repository...")
    subprocess.run("git clone https://github.com/Joacaldog/FrankPEPstein.git", shell=True, check=True)
else:
    print("Repository already exists.")

# --- 1.5 Install CondaColab ---
try:
    import condacolab
    print("condacolab already installed.")
except ImportError:
    print("Installing CondaColab...")
    subprocess.run("pip install -q condacolab", shell=True, check=True)
    import condacolab
    condacolab.install()
    print("Please restart the runtime if asked, then run this cell again.")

# --- 2. Create Conda Environment 'FrankPEPstein' ---
env_path = "/usr/local/envs/FrankPEPstein"
if os.path.exists(env_path):
    print(f"Environment 'FrankPEPstein' already exists at {env_path}. Skipping creation.")
else:
    print("Creating 'FrankPEPstein' environment with Python 3.10 (this may take a few minutes)...")
    # Create environment with all dependencies including Modeller
    subprocess.run("mamba create -n FrankPEPstein -q -y -c conda-forge -c salilab openbabel biopython fpocket joblib tqdm py3dmol vina python=3.10 salilab::modeller", shell=True, check=True)

# --- 3. Configure Path for Colab Usage ---
# Since Colab runs on the 'base' kernel, we need to manually add the new env to paths
site_packages = f"{env_path}/lib/python3.10/site-packages"

if site_packages not in sys.path:
    sys.path.append(site_packages)

# Add binary path for tools like fpocket, obabel, etc.
os.environ['PATH'] = f"{env_path}/bin:" + os.environ['PATH']

print(f"Environment 'FrankPEPstein' created and configured.")

# --- 4. Run Notebook Setup Utils ---
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
