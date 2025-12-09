
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

# Step 0: Dependency Installation
nb.cells.append(new_markdown_cell("## 0. Setup & Dependencies"))
nb.cells.append(new_code_cell("""
# Install CondaColab
!pip install -q condacolab
import condacolab
condacolab.install()
"""))

nb.cells.append(new_code_cell("""
# Install Core Bio-Dependencies via Mamba
import sys
!mamba install -q -c conda-forge openbabel biopython fpocket joblib tqdm py3dmol
"""))

nb.cells.append(new_code_cell("""
# Setup External Tools (Click & ADFR)
import os

# Create folders
os.makedirs("utilities", exist_ok=True)
os.makedirs("scripts", exist_ok=True)

# User must ensure these files are present or uploaded
if os.path.exists("Click.tar.gz") and os.path.getsize("Click.tar.gz") > 0:
    !tar -xzf Click.tar.gz -C utilities/
    # Add Click to PATH (Placeholder - depends on actual folder structure inside tar)
    # os.environ['PATH'] += ":/content/utilities/click/bin" 
    print("Click installed (supposedly)")
else:
    print("WARNING: Click.tar.gz missing or empty. Please upload the 'click' tool tarball.")

if os.path.exists("ADFRsuite_x86_64Linux_1.0.tar.gz"):
    !tar -xzf ADFRsuite_x86_64Linux_1.0.tar.gz -C utilities/
    # Setup ADFR path
    # os.environ['PATH'] += ":/content/utilities/ADFRsuite/bin"
    print("ADFR Suite installed")
else:
    print("WARNING: ADFRsuite_x86_64Linux_1.0.tar.gz missing. Please upload it.")

# Setup Modeller License
license_key = input("Enter your MODELLER license key (or press Enter if configured in file): ")
if license_key:
    # Update config file logic here if needed
    print(f"Key {license_key} received (Configuration logic to be added)")
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
nb.cells.append(new_markdown_cell("## 2. Pocket Detection (fpocket)"))
nb.cells.append(new_code_cell("""
import subprocess

# Run fpocket
print(f"Running fpocket on {receptor_filename}...")
subprocess.run(f"fpocket -f {receptor_filename}", shell=True, check=True)

output_folder = f"{receptor_filename}_out"
if os.path.exists(output_folder):
    print("fpocket finished successfully.")
    pockets_dir = os.path.join(output_folder, "pockets")
    pockets = [f for f in os.listdir(pockets_dir) if f.endswith(".pdb")]
    print(f"Found {len(pockets)} pockets.")
else:
    print("Error: fpocket output folder not found.")
"""))

# Step 3: Visualization & Selection
nb.cells.append(new_markdown_cell("## 3. Pocket Selection"))
nb.cells.append(new_code_cell("""
import py3Dmol
import ipywidgets as widgets
from IPython.display import display

# Simple Pocket Selector Widget
pocket_dropdown = widgets.Dropdown(
    options=sorted(pockets),
    description='Select Pocket:',
    disabled=False,
)

def view_pocket(pocket_file):
    view = py3Dmol.view(width=800, height=600)
    
    # Load Receptor
    with open(receptor_filename, 'r') as f:
        view.addModel(f.read(), "pdb")
    view.setStyle({'cartoon': {'color': 'white'}})
    
    # Load Selected Pocket (Red spheres)
    pocket_path = os.path.join(pockets_dir, pocket_file)
    with open(pocket_path, 'r') as f:
        view.addModel(f.read(), "pdb")
    view.setStyle({'model': -1}, {'sphere': {'color': 'red', 'opacity': 0.7}})
    
    view.zoomTo()
    view.show()

# Link widget to view
widgets.interactive(view_pocket, pocket_file=pocket_dropdown)
"""))

# Step 4: Pocket Extraction
nb.cells.append(new_markdown_cell("## 4. Pocket Extraction & Superposer Prep"))
nb.cells.append(new_code_cell("""
# Extract selected pocket coordinates and size
selected_pocket_file = pocket_dropdown.value
print(f"Selected: {selected_pocket_file}")

# Logic to parse extraction (Simulated based on known file format)
# In reality, we might read the PDB and calculate min/max coordinates
# For now, we will assume generic values or implement a simple parser

def get_box_center_size(pdb_file):
    # This is a placeholder for the actual logic to get box center/size
    # You would use BioPython here
    return (0, 0, 0), (20, 20, 20)

center, size = get_box_center_size(os.path.join(pockets_dir, selected_pocket_file))
print(f"Box Center: {center}, Size: {size}")
"""))

# Step 5: Superposer
nb.cells.append(new_markdown_cell("## 5. Superposer (Fragment Generation)"))
nb.cells.append(new_code_cell("""
# Run Superposer
# Note: pepbdb_folder must be set by user
pepbdb_folder = "/content/pepbdb" # Placeholder
print("Running Superposer...")

# cmd = f"python scripts/superposerV5.2_leave1out.py -T {receptor_filename} -d {pepbdb_folder} ..."
# !{cmd}
print("Superposer command placeholder (Requires PepBDB and adjust arguments)")
"""))

# Step 6: FrankVINA (FNKPSTN)
nb.cells.append(new_markdown_cell("## 6. FrankVINA (FNKPSTN)"))
nb.cells.append(new_code_cell("""
# Run FrankVINA FNKPSTN
print("Running FrankVINA FNKPSTN...")
# !python scripts/frankVINA_FNKPSTN.py receptor.pdb 4
"""))

# Step 7: Patch Clustering
nb.cells.append(new_markdown_cell("## 7. Patch Clustering"))
nb.cells.append(new_code_cell("""
# Run Patch Clustering
print("Running Patch Clustering...")
# !python scripts/patch_clustering_V8.7.py -w 6 -t 4
"""))

# Step 8: FrankVINA (V3)
nb.cells.append(new_markdown_cell("## 8. FrankVINA (V3)"))
nb.cells.append(new_code_cell("""
# Run FrankVINA V3
print("Running FrankVINA V3...")
# !python scripts/frankVINA_V3.py receptor.pdb 4 10
"""))

# Step 9: Download
nb.cells.append(new_markdown_cell("## 9. Download Results"))
nb.cells.append(new_code_cell("""
# Zip and Download
!zip -r frankpepstein_results.zip results_folder/
files.download('frankpepstein_results.zip')
"""))

# Save Notebook
with open('FrankPEPstein.ipynb', 'w') as f:
    nbf.write(nb, f)

print("Notebook generated: FrankPEPstein.ipynb")
