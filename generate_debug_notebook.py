
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

nb = new_notebook()

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

# Optional: Factory Reset
nb.cells.append(new_code_cell("""
#@title Factory Reset (Uncomment and run to wipe runtime)
# from google.colab import runtime
# runtime.unassign()
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

# 0.2 Main Setup
nb.cells.append(new_code_cell("""
#@title 0.2 Install Dependencies & Setup Tools
#@markdown This cell clones the repository, installs bio-dependencies, and sets up external tools.

import os
import sys
import subprocess

# --- 1. Clone Repository ---
if not os.path.exists("FrankPEPstein"):
    print("Cloning repository...")
    !git clone https://github.com/Joacaldog/FrankPEPstein.git
else:
    print("Repository already exists.")

# --- 2. Install Core Bio-Dependencies ---
print("Installing bio-dependencies (this may take a few minutes)...")
!mamba install -q -c conda-forge -c salilab openbabel biopython fpocket joblib tqdm py3dmol vina python=3.10 salilab::modeller

# --- 3. Setup External Tools & Config ---
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
"""))

# Step 1: File Upload
nb.cells.append(new_markdown_cell("## 1. Input Data: Upload Receptor"))
nb.cells.append(new_code_cell("""
from google.colab import files
import os

print("Upload your Receptor PDB file (e.g., receptor.pdb)")
uploaded = files.upload()
receptor_filename = list(uploaded.keys())[0]
print(f"Receptor uploaded: {receptor_filename}")
"""))


# Step 2: Pocket Detection
# nb.cells.append(new_markdown_cell("## 2. Pocket Detection (fpocket)"))
# nb.cells.append(new_code_cell("""
# import subprocess
# 
# # Run fpocket
# print(f"Running fpocket on {receptor_filename}...")
# subprocess.run(f"fpocket -f {receptor_filename}", shell=True, check=True)
# 
# output_folder = f"{receptor_filename}_out"
# if os.path.exists(output_folder):
#     print("fpocket finished successfully.")
#     pockets_dir = os.path.join(output_folder, "pockets")
#     pockets = [f for f in os.listdir(pockets_dir) if f.endswith(".pdb")]
#     print(f"Found {len(pockets)} pockets.")
# else:
#     print("Error: fpocket output folder not found.")
# """))
# 
# # Step 3: Visualization & Selection
# nb.cells.append(new_markdown_cell("## 3. Pocket Selection"))
# nb.cells.append(new_code_cell("""
# import py3Dmol
# import ipywidgets as widgets
# from IPython.display import display
# 
# # Simple Pocket Selector Widget
# pocket_dropdown = widgets.Dropdown(
#     options=sorted(pockets),
#     description='Select Pocket:',
#     disabled=False,
# )
# 
# def view_pocket(pocket_file):
#     view = py3Dmol.view(width=800, height=600)
#     
#     # Load Receptor
#     with open(receptor_filename, 'r') as f:
#         view.addModel(f.read(), "pdb")
#     view.setStyle({'cartoon': {'color': 'white'}})
#     
#     # Load Selected Pocket (Red spheres)
#     pocket_path = os.path.join(pockets_dir, pocket_file)
#     with open(pocket_path, 'r') as f:
#         view.addModel(f.read(), "pdb")
#     view.setStyle({'model': -1}, {'sphere': {'color': 'red', 'opacity': 0.7}})
#     
#     view.zoomTo()
#     view.show()
# 
# # Link widget to view
# widgets.interactive(view_pocket, pocket_file=pocket_dropdown)
# """))
# 
# # Step 4: Pocket Extraction
# nb.cells.append(new_markdown_cell("## 4. Pocket Extraction & Superposer Prep"))
# nb.cells.append(new_code_cell("""
# # Extract selected pocket coordinates and size
# selected_pocket_file = pocket_dropdown.value
# print(f"Selected: {selected_pocket_file}")
# 
# # Logic to parse extraction (Simulated based on known file format)
# # In reality, we might read the PDB and calculate min/max coordinates
# # For now, we will assume generic values or implement a simple parser
# 
# def get_box_center_size(pdb_file):
#     # This is a placeholder for the actual logic to get box center/size
#     # You would use BioPython here
#     return (0, 0, 0), (20, 20, 20)
# 
# center, size = get_box_center_size(os.path.join(pockets_dir, selected_pocket_file))
# print(f"Box Center: {center}, Size: {size}")
# """))
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

