import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
import os

def read_file_content(path):
    if os.path.exists(path):
        with open(path, 'r') as f:
            return f.read()
    else:
        print(f"Warning: File not found {path}")
        return f"# Error: File {path} not found"

nb = new_notebook()

# --- Title Cell ---
nb.cells.append(new_markdown_cell("""
# **FrankPEPstein: De Novo Peptide Generation**
### *Fragment-based generation of high-affinity peptides for protein pockets.*

**Steps:**
1.  **Setup**: Install dependencies (ADFR, OpenBabel, etc.) and download the database.
2.  **Input**: Upload your receptor and select a pocket.
3.  **Run**: Configure parameters and launch the generation pipeline.
"""))

# --- Step 0: Setup ---
content_0 = read_file_content("functions/step_0_setup.py")
nb.cells.append(new_code_cell(content_0))

# --- Step 1: Input & Pocket ---
content_1 = read_file_content("functions/step_1_pocket_processing.py")
nb.cells.append(new_code_cell(content_1))

# --- Step 2: Pipeline Execution ---
content_2 = read_file_content("functions/step_2_FrankPEPstein.py")
nb.cells.append(new_code_cell(content_2))

# --- Save ---
output_filename = 'FrankPEPstein.ipynb'
with open(output_filename, 'w') as f:
    nbf.write(nb, f)

print(f"✅ Generated {output_filename} successfully.")
print("You can now upload this notebook to Google Colab.")
