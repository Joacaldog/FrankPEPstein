
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
import os

nb = new_notebook()

# --- Load Cell Contents from Function Tests ---
# We assume the user wants the exact content of these files in the notebook cells.

def read_file_content(path):
    with open(path, 'r') as f:
        return f.read()

# Cell 1: Setup
setup_content = read_file_content('function_tests/step_0_setup.py')

# Cell 2: Combined Workflow (Input -> Selection -> Extraction)
combined_workflow_content = read_file_content('function_tests/step_1_4_workflow.py')

# Cell 3: FrankPEPstein Pipeline
pipeline_content = read_file_content('function_tests/step_5_8_pipeline.py')

# --- Build Notebook ---

# Title
nb.cells.append(new_markdown_cell("""
# FrankPEPstein: Incremental Debug Notebook

This notebook tests the pipeline using a consolidated workflow:
1.  **Setup**: Installs dependencies in a dedicated environment.
2.  **Workflow**: Uploads Receptor -> Detects/Uploads Pockets -> Selects Pocket -> Calculates Box.
"""))

# Cell 1: Setup
nb.cells.append(new_code_cell(setup_content))

# Cell 2: Combined Workflow
nb['cells'].append(nbf.v4.new_code_cell(combined_workflow_content))
    
# Cell 3: FrankPEPstein Pipeline
if pipeline_content:
    nb['cells'].append(nbf.v4.new_code_cell(pipeline_content))
else:
    print("Warning: Pipeline content empty or file not found")

output_filename = 'FrankPEPstein_DEBUG.ipynb' # Define output_filename
with open(output_filename, 'w') as f:
    nbf.write(nb, f)

print(f"Notebook generated: {output_filename}")
