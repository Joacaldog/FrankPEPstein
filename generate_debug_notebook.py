
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

nb = new_notebook()

# Read local notebook_utils.py for patching
with open('scripts/notebook_utils.py', 'r') as f:
    local_utils_content = f.read()

# Title and Introduction
nb.cells.append(new_markdown_cell("""
# FrankPEPstein: Interactive Peptide Fragment Design

This notebook implements the FrankPEPstein pipeline for designing peptide fragments binding to a specific protein pocket.

**Steps:**
1.  Setup & Dependencies
2.  Input Data (Receptor)
3.  Pocket Selection
4.  Fragment Generation & Ranking
"""))


# Step 0: Dependency Installation
nb.cells.append(new_markdown_cell("## 0. Setup & Dependencies"))

# 0.1 CondaColab (Restarts Kernel)
nb.cells.append(new_code_cell("""
#@title 0.1 Install CondaColab (Running this will restart the session)
import sys
import os

# Helper to suppress output
class SuppressStdout:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')
    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

try:
    import condacolab
    with SuppressStdout():
        condacolab.check()
except ImportError:
    print("Installing CondaColab...")
    !pip install -q condacolab
    import condacolab
    with SuppressStdout():
        condacolab.install()
print("CondaColab installed.")
"""))

# Main Setup
setup_code_src = """
#@title Install Dependencies & Setup Tools
#@markdown This cell clones the repository and creates the 'FrankPEPstein' environment with Python 3.10.

import os
import sys
import subprocess

# --- 1. Clone Repository ---
if not os.path.exists("FrankPEPstein"):
    print("Cloning repository...")
    !git clone https://github.com/Joacaldog/FrankPEPstein.git
else:
    print("Repository already exists.")

# --- 2. Create Conda Environment 'FrankPEPstein' ---
print("Creating 'FrankPEPstein' environment with Python 3.10...")
# Create environment with all dependencies including Modeller
!mamba create -n FrankPEPstein -q -y -c conda-forge -c salilab openbabel biopython fpocket joblib tqdm py3dmol vina python=3.10 salilab::modeller

# --- 3. Configure Path for Colab Usage ---
# Since Colab runs on the 'base' kernel, we need to manually add the new env to paths
env_path = "/usr/local/envs/FrankPEPstein"
site_packages = f"{env_path}/lib/python3.10/site-packages"

if site_packages not in sys.path:
    sys.path.append(site_packages)

# Add binary path for tools like fpocket, obabel, etc.
os.environ['PATH'] = f"{env_path}/bin:" + os.environ['PATH']

print(f"Environment 'FrankPEPstein' created and configured.")

# --- PATCH: Update notebook_utils.py with local changes ---
# This ensures we use the latest path logic without needing a git push
patched_utils_content = r'''<<NOTEBOOK_UTILS_CONTENT>>'''

os.makedirs("FrankPEPstein/scripts", exist_ok=True)
with open("FrankPEPstein/scripts/notebook_utils.py", "w") as f:
    f.write(patched_utils_content)
print("Patched notebook_utils.py with latest local version.")

# --- 4. Setup External Tools & Config ---
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
"""

# Inject local content
setup_code_final = setup_code_src.replace("<<NOTEBOOK_UTILS_CONTENT>>", local_utils_content)
nb.cells.append(new_code_cell(setup_code_final))

# Step 1-3: Input, Detection, and Selection (Combined)
nb.cells.append(new_markdown_cell("## 1. Input & Pocket Selection"))
nb.cells.append(new_code_cell("""
#@title Select Receptor & Detect Pockets
#@markdown **Instructions:**
#@markdown 1. Upload your Receptor PDB.
#@markdown 2. Choose Mode: **Auto Detect** (runs fpocket) or **Manual Upload** (upload your specific pocket PDB).
#@markdown 3. Select the pocket from the dropdown to visualize.

import os
import subprocess
import py3Dmol
import ipywidgets as widgets
from google.colab import files
from IPython.display import display

# --- configuration ---
detection_mode = "Auto Detect" #@param ["Auto Detect", "Manual Upload"]

# Global variables for next steps
receptor_filename = None
pockets_dir = "pockets_upload" # Default for manual
final_pockets_list = []

# --- 1. Upload Receptor ---
print(f"--- Upload Receptor PDB ({detection_mode}) ---")
uploaded_r = files.upload()

if not uploaded_r:
    print("No receptor file uploaded.")
else:
    receptor_filename = list(uploaded_r.keys())[0]
    print(f"Receptor: {receptor_filename}")

    # --- 2. Pocket Handling ---
    if detection_mode == "Auto Detect":
        print(f"\\nRunning fpocket on {receptor_filename}...")
        try:
            # Fix: Quotes for filenames with spaces
            subprocess.run(f"fpocket -f '{receptor_filename}'", shell=True, check=True)
            
            # Robust folder finding
            base_name = os.path.splitext(receptor_filename)[0]
            possible_folders = [f"{receptor_filename}_out", f"{base_name}_out"]
            output_folder = next((f for f in possible_folders if os.path.exists(f)), None)

            if output_folder:
                pockets_dir = os.path.join(output_folder, "pockets")
                if os.path.exists(pockets_dir):
                    final_pockets_list = [f for f in os.listdir(pockets_dir) if f.endswith(".pdb")]
                    print(f"Auto-detection finished. Found {len(final_pockets_list)} pockets.")
                else:
                    print("Warning: pockets subdirectory not found.")
            else:
                print("Error: fpocket output not found.")
                
        except subprocess.CalledProcessError:
             print("Error running fpocket.")

    elif detection_mode == "Manual Upload":
        print(f"\\n--- Upload Manual Pocket PDB ---")
        os.makedirs(pockets_dir, exist_ok=True)
        uploaded_p = files.upload()
        if uploaded_p:
            for p_file in uploaded_p.keys():
                # Move to a pockets folder to keep structure consistent
                os.rename(p_file, os.path.join(pockets_dir, p_file))
                final_pockets_list.append(p_file)
            print(f"Manual upload finished. Available pockets: {len(final_pockets_list)}")

    # --- 3. Visualization & Selection ---
    if final_pockets_list:
        print("\\n--- Pocket Selection & Visualization ---")
        
        pocket_dropdown = widgets.Dropdown(
            options=sorted(final_pockets_list),
            description='Select Pocket:',
            disabled=False,
        )

        def view_pocket(pocket_file):
            view = py3Dmol.view(width=800, height=600)
            
            # 1. Receptor Surface (White, Transparent)
            with open(receptor_filename, 'r') as f:
                view.addModel(f.read(), "pdb")
            view.setStyle({}) 
            view.addSurface(py3Dmol.SES, {'opacity': 0.8, 'color': 'white'})
            
            # 2. Selected Pocket (Red Spheres)
            full_path = os.path.join(pockets_dir, pocket_file)
            if os.path.exists(full_path):
                with open(full_path, 'r') as f:
                    view.addModel(f.read(), "pdb")
                view.setStyle({'model': -1}, {'sphere': {'color': 'red', 'opacity': 0.7}})
            else:
                print(f"Error: Could not find {full_path}")

            view.zoomTo()
            view.show()
            
        display(widgets.interactive(view_pocket, pocket_file=pocket_dropdown))
    else:
        print("No pockets available to select.")
""", metadata={"cellView": "form"}))
# 
# # Step 4: Pocket Extraction
nb.cells.append(new_markdown_cell("## 4. Pocket Extraction & Superposer Prep"))
nb.cells.append(new_code_cell("""
#@title 4. Pocket Extraction & Box Generation
#@markdown This step extracts the selected pocket and calculates the grid box center and size.

import os
from Bio.PDB import PDBParser, PDBIO, Select

# --- Configuration ---
try:
    selected_pocket_file = pocket_dropdown.value
    print(f"Selected Pocket File: {selected_pocket_file}")
    
    pocket_path = os.path.join(pockets_dir, selected_pocket_file)
    print(f"Pocket Path: {pocket_path}")
    
except NameError:
    print("Error: 'pocket_dropdown' or 'pockets_dir' not defined. Did you run the previous cell?")

# --- Helper Functions ---
def get_box_center_size(pdb_file, buffer=0.0):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("pocket", pdb_file)
    coords = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    coords.append(atom.get_coord())
    
    if not coords:
        return None, None

    min_coord = [min([c[i] for c in coords]) for i in range(3)]
    max_coord = [max([c[i] for c in coords]) for i in range(3)]
    
    center = [(min_coord[i] + max_coord[i]) / 2 for i in range(3)]
    size = [(max_coord[i] - min_coord[i]) + buffer for i in range(3)]
    
    return center, size

# --- Main Extraction Logic ---
if os.path.exists(pocket_path):
    print("Calculating box parameters...")
    center, size = get_box_center_size(pocket_path, buffer=10.0) 
    
    if center:
        center_str = f"{center[0]:.3f} {center[1]:.3f} {center[2]:.3f}"
        size_str = f"{size[0]:.3f} {size[1]:.3f} {size[2]:.3f}"
        
        print("-" * 30)
        print(f"Box Center: {center_str}")
        print(f"Box Size:   {size_str}")
        print("-" * 30)
        
        box_center = center
        box_size = size
        print("Pocket parameters ready for FrankPEPstein.")
        
    else:
        print("Error: Could not calculate coordinates from pocket file.")
else:
    print(f"Error: Pocket file not found at {pocket_path}")
"""))
# 
# # Step 5: Superposer
# nb.cells.append(new_markdown_cell("## 5. Superposer (Fragment Generation)"))
# nb.cells.append(new_code_cell("""
# # Run Superposer
# # Note: pepbdb_folder must be set by user
# pepbdb_folder = "/content/pepbdb" # Placeholder
# print("Running Superposer...")
# 
# # cmd = f"python scripts/superposerV5.2_leave1out.py -T {receptor_filename} -d {pepbdb_folder} ..."
# # !{cmd}
# print("Superposer command placeholder (Requires PepBDB and adjust arguments)")
# """))
# 
# # Step 6: FrankVINA (FNKPSTN)
# nb.cells.append(new_markdown_cell("## 6. FrankVINA (FNKPSTN)"))
# nb.cells.append(new_code_cell("""
# # Run FrankVINA FNKPSTN
# print("Running FrankVINA FNKPSTN...")
# # !python scripts/frankVINA_FNKPSTN.py receptor.pdb 4
# """))
# 
# # Step 7: Patch Clustering
# nb.cells.append(new_markdown_cell("## 7. Patch Clustering"))
# nb.cells.append(new_code_cell("""
# # Run Patch Clustering
# print("Running Patch Clustering...")
# # !python scripts/patch_clustering_V8.7.py -w 6 -t 4
# """))
# 
# # Step 8: FrankVINA (V3)
# nb.cells.append(new_markdown_cell("## 8. FrankVINA (V3)"))
# nb.cells.append(new_code_cell("""
# # Run FrankVINA V3
# print("Running FrankVINA V3...")
# # !python scripts/frankVINA_V3.py receptor.pdb 4 10
# """))
# 
# # Step 9: Download
# nb.cells.append(new_markdown_cell("## 9. Download Results"))
# nb.cells.append(new_code_cell("""
# # Zip and Download
# !zip -r frankpepstein_results.zip results_folder/
# files.download('frankpepstein_results.zip')
# """))

# Save Notebook
with open('FrankPEPstein_DEBUG.ipynb', 'w') as f:
    nbf.write(nb, f)

print("Notebook generated: FrankPEPstein_DEBUG.ipynb")

