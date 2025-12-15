#@title 1. Pocket Selection & Gridbox Generation
#@markdown **Instructions:**
#@markdown 1. Select **Detection Mode**.
#@markdown 2. If **Auto Detect**, run fpocket and select a predicted pocket.
#@markdown 3. If **Manual Upload**, upload your pre-defined pocket PDB.
#@markdown 4. Use the controls to adjust the gridbox (Cyan Box) if needed.
#@markdown 5. Click **Confirm & Extract**.

import os
import shutil
import subprocess
import sys
import glob
import re
import ipywidgets as widgets
from IPython.display import display, clear_output
import py3Dmol

# --- Configuration ---
# Ensure directories exist
initial_path = os.getcwd()
pockets_dir = os.path.join(initial_path, "pockets")
fpocket_storage_dir = os.path.join(initial_path, "fpocket_pockets")
os.makedirs(pockets_dir, exist_ok=True)
os.makedirs(fpocket_storage_dir, exist_ok=True)

# State Variables
receptor_filename = "receptor.pdb"
if not os.path.exists(receptor_filename):
    # Try finding it in pipeline state or current dir
    pass 

pipeline_state_file = "pipeline_state.json"
import json
def save_pipeline_state(data):
    current = {}
    if os.path.exists(pipeline_state_file):
        try:
            with open(pipeline_state_file, 'r') as f: current = json.load(f)
        except: pass
    current.update(data)
    with open(pipeline_state_file, 'w') as f: json.dump(current, f)

# --- Part 1: Receptor Input ---
receptor_upload_widget = widgets.FileUpload(description="Upload Receptor (pdb)", accept=".pdb", multiple=False, layout=widgets.Layout(width='300px'))
receptor_status = widgets.Output()

def handle_receptor_upload(change):
    receptor_status.clear_output()
    with receptor_status:
        if not receptor_upload_widget.value: return
        # ipywidgets 7/8 compat
        upl_file = list(receptor_upload_widget.value.values())[0] if isinstance(receptor_upload_widget.value, dict) else receptor_upload_widget.value[0]
        content = upl_file['content']
        with open(receptor_filename, "wb") as f:
            f.write(content)
        print("✅ Receptor uploaded successfully.")

receptor_upload_widget.observe(handle_receptor_upload, names='value')

print("Step 1.1: Load Receptor")
display(widgets.VBox([receptor_upload_widget, receptor_status]))

with receptor_status:
    if os.path.exists(receptor_filename):
        print(f"✅ Receptor file present: {receptor_filename}")
    else:
        print("Waiting for receptor upload...")

# --- Part 2: Detection & Upload Logic ---

mode_selector = widgets.ToggleButtons(
    options=['Auto Detect', 'Manual Upload'],
    description='Mode:',
    disabled=False,
    button_style='',
)

# Auto Detect Widgets
detect_btn = widgets.Button(description="Run fpocket", button_style='info')
pocket_dropdown = widgets.Dropdown(description="Select Pocket:", options=[], disabled=True)

# Manual Upload Widgets
upload_btn = widgets.FileUpload(description="Upload Pocket PDB", accept=".pdb", multiple=False, layout=widgets.Layout(width='300px'))

# Containers
log_output_1 = widgets.Output()

def run_fpocket(b):
    log_output_1.clear_output()
    with log_output_1:
        if not os.path.exists(receptor_filename):
            print("❌ Receptor not found! Please run setup/upload first.")
            return
            
        print("Running fpocket... (this may take a minute)")
        # Clean previous
        subprocess.run(f"rm -rf {fpocket_storage_dir}/*", shell=True)
        
        cmd = f"fpocket -f {receptor_filename}"
        subprocess.run(cmd, shell=True)
        
        # Check results
        # fpocket output: receptor_out/pockets/pocketX_atm.pdb
        out_dir = receptor_filename.replace(".pdb", "_out")
        pockets_found = glob.glob(os.path.join(out_dir, "pockets", "pocket*_atm.pdb"))
        
        if pockets_found:
            print(f"✅ Found {len(pockets_found)} pockets.")
            # Copy to storage
            options = []
            for p in pockets_found:
                basename = os.path.basename(p)
                dest = os.path.join(fpocket_storage_dir, basename)
                shutil.copy(p, dest)
                options.append(basename)
            
            # Sort naturally
            options.sort(key=lambda x: int(re.search(r'\d+', x).group()) if re.search(r'\d+', x) else 0)
            
            pocket_dropdown.options = options
            pocket_dropdown.disabled = False
            if options: pocket_dropdown.value = options[0]
            
            # Trigger Box Init
            initialize_ui(None)
        else:
            print("❌ No pockets found.")

def handle_upload(change):
    log_output_1.clear_output()
    with log_output_1:
        if not upload_btn.value: return
        
        # Get file
        # ipywidgets 7 vs 8 compat logic
        upl_file = list(upload_btn.value.values())[0] if isinstance(upload_btn.value, dict) else upload_btn.value[0]
        
        content = upl_file['content']
        name = upl_file['name']
        
        dest = os.path.join(fpocket_storage_dir, name)
        with open(dest, "wb") as f:
            f.write(content)
            
        print(f"✅ Uploaded {name}")
        
        pocket_dropdown.options = [name]
        pocket_dropdown.value = name
        pocket_dropdown.disabled = False
        
        # Trigger Box Init
        initialize_ui(None)

detect_btn.on_click(run_fpocket)
upload_btn.observe(handle_upload, names='value')

def on_pocket_change(change):
    if change['type'] == 'change' and change['name'] == 'value':
        initialize_ui(None)

pocket_dropdown.observe(on_pocket_change)

# Selection UI
selection_ui = widgets.VBox([
    mode_selector,
    widgets.HBox([detect_btn, pocket_dropdown]),
    upload_btn,
    log_output_1
])

# Visibiity logic
def update_mode_ui(change):
    if mode_selector.value == 'Auto Detect':
        detect_btn.layout.display = 'block'
        pocket_dropdown.layout.display = 'block'
        upload_btn.layout.display = 'none'
    else:
        detect_btn.layout.display = 'none'
        # pocket_dropdown.layout.display = 'none' # Keep dropdown to show selected uploaded file?
        # Usually manual upload implies ONE file. But let's keep it for consistency.
        upload_btn.layout.display = 'block'

mode_selector.observe(update_mode_ui, names='value')
update_mode_ui(None) # Init

print("Select Detection Mode:")
display(selection_ui)


# --- Part 2: Box Calculation & Extraction (Isolated) ---
# This matches the previous updated logic

# Helper Functions (Subprocess)
def run_processing_isolated(receptor_path, pocket_path, output_pocket_path, mode="extract", buffer=0.0):
    """
    Runs extraction/processing and box calculation in isolated environment.
    """
    script_content = f"""
import sys
import os
from Bio import PDB
from Bio.PDB import PDBParser, PDBIO, Select, NeighborSearch

def process_and_box(receptor_file, pocket_file, output_file, mode, buffer_val):
    try:
        parser = PDBParser(QUIET=True)
        
        # 1. Load Pocket
        pocket_struct = parser.get_structure("pocket", pocket_file)
        
        # Determine atoms for Box Calculation
        atoms_for_box = []
        residues_for_saving = [] 
        
        if mode == 'extract':
             # Fpocket mode: Load Receptor, Find Neighbors (5A)
             receptor_struct = parser.get_structure("receptor", receptor_file)
             pocket_atoms = [atom for atom in pocket_struct.get_atoms()]
             if not pocket_atoms:
                 print("ERROR: No atoms in pocket file")
                 return
                 
             receptor_atoms = list(receptor_struct.get_atoms())
             ns = NeighborSearch(receptor_atoms)
             selected_residues = set()
             for p_atom in pocket_atoms:
                 nearby = ns.search(p_atom.get_coord(), 5.0, level='R')
                 for res in nearby:
                     selected_residues.add((res.parent.id, res.id))
            
             class PocketSelect(Select):
                 def accept_residue(self, residue):
                     return (residue.parent.id, residue.id) in selected_residues
                     
             io = PDBIO()
             io.set_structure(receptor_struct)
             temp_out = output_file + ".tmp"
             io.save(temp_out, PocketSelect())
             
             saved_struct = parser.get_structure("saved", temp_out)
             for model in saved_struct:
                 for chain in model:
                     chain.id = 'p'
                     for residue in chain:
                         for atom in residue:
                             atoms_for_box.append(atom)
             io.set_structure(saved_struct)
             io.save(output_file)
             os.remove(temp_out)
             
        elif mode == 'direct':
             for model in pocket_struct:
                 for chain in model:
                     chain.id = 'p'
                     for residue in chain:
                         for atom in residue:
                             atoms_for_box.append(atom)
             io = PDBIO()
             io.set_structure(pocket_struct)
             io.save(output_file)
             
        if not atoms_for_box:
            print("ERROR: No atoms/residues for box calculation")
            return

        coords = [a.get_coord() for a in atoms_for_box]
        min_coord = [min([c[i] for c in coords]) for i in range(3)]
        max_coord = [max([c[i] for c in coords]) for i in range(3)]
        
        center = [(min_coord[i] + max_coord[i]) / 2 for i in range(3)]
        size = [(max_coord[i] - min_coord[i]) + float(buffer_val) for i in range(3)]
        
        print(f"CENTER:{{center[0]}},{{center[1]}},{{center[2]}}")
        print(f"SIZE:{{size[0]}},{{size[1]}},{{size[2]}}")
        print("SUCCESS")

    except Exception as e:
        print(f"ERROR:{{e}}")

if __name__ == "__main__":
    process_and_box(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], float(sys.argv[5]))
"""
    script_name = "process_box_isolated.py"
    with open(script_name, "w") as f:
        f.write(script_content)
        
    python_exe = "/usr/local/envs/FrankPEPstein/bin/python"
    if not os.path.exists(python_exe): python_exe = sys.executable

    try:
        result = subprocess.run(
            [python_exe, script_name, receptor_path, pocket_path, output_pocket_path, mode, str(buffer)],
            capture_output=True, text=True, check=True
        )
        
        center, size = None, None
        success = False
        
        for line in result.stdout.splitlines():
            if line.startswith("CENTER:"):
                center = [float(x) for x in line.split(":")[1].split(",")]
            elif line.startswith("SIZE:"):
                size = [float(x) for x in line.split(":")[1].split(",")]
            elif line.startswith("SUCCESS"):
                success = True
            
        return center, size, success
    except subprocess.CalledProcessError as e:
        return None, None, False

def get_pocket_volume_hull(pocket_pdb_path):
    """Calculates Convex Hull volume of C-alpha atoms for visualization."""
    script_content = f"""
import sys
import numpy as np
from Bio.PDB import PDBParser
from scipy.spatial import ConvexHull

def calc_vol(pdb_file):
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("struct", pdb_file)
        # Use all atoms
        points = []
        for atom in structure.get_atoms():
             points.append(atom.get_coord())
             
        if len(points) < 4:
            print("VOL:0.0")
            return
            
        hull = ConvexHull(points)
        print(f"VOL:{{hull.volume}}")
    except Exception as e:
        print("VOL:0.0")

if __name__ == "__main__":
    calc_vol(sys.argv[1])
"""
    script_name = "calc_volume_isolated.py"
    with open(script_name, "w") as f:
        f.write(script_content)
        
    python_exe = "/usr/local/envs/FrankPEPstein/bin/python"
    if not os.path.exists(python_exe): python_exe = sys.executable

    try:
        result = subprocess.run(
            [python_exe, script_name, pocket_pdb_path],
            capture_output=True, text=True, check=True
        )
        for line in result.stdout.splitlines():
            if line.startswith("VOL:"):
                return float(line.split(":")[1])
    except:
        pass
    return 0.0

# --- Custom Widget Helpers (Box Adjustment) ---

def create_control_group(label, initial_val, step, color_hex):
    # Style buttons
    btn_layout = widgets.Layout(width='40px')
    
    minus_button = widgets.Button(description='\u25c0', layout=btn_layout)
    plus_button = widgets.Button(description='\u25b6', layout=btn_layout)
    
    val_widget = widgets.FloatText(value=initial_val, step=step, description='', layout=widgets.Layout(width='80px'))
    
    def on_minus(b):
        val_widget.value -= step
        update_visual(None)
        
    def on_plus(b):
        val_widget.value += step
        update_visual(None)
            
    minus_button.on_click(on_minus)
    plus_button.on_click(on_plus)
    val_widget.observe(lambda c: update_visual(None) if c['type'] == 'change' and c['name'] == 'value' else None)

    label_html = widgets.HTML(f"<b style='color:{color_hex}; font-size:14px; margin-right:5px;'>{label}</b>")
    
    box = widgets.HBox([label_html, minus_button, val_widget, plus_button], layout=widgets.Layout(align_items='center'))
    return val_widget, box

controls = {}
viz_output = widgets.Output() # Holds the 3D widget
output_log = widgets.Output() # Holds text logs

def initialize_ui(b):
    output_log.clear_output()
    viz_output.clear_output()
    
    if not pocket_dropdown.value:
        # Don't print error if just init, wait for selection
        return

    # 1. Calc Initial
    selected_pocket = pocket_dropdown.value
    src_pocket_path = os.path.join(fpocket_storage_dir, selected_pocket)
    temp_pocket_path = os.path.join(pockets_dir, "temp_calc.pdb")
    
    # Determine mode from toggle
    mode = "extract" if mode_selector.value == "Auto Detect" else "direct"
    
    with output_log: print(f"Processing {selected_pocket}...")
    center, size, success = run_processing_isolated(receptor_filename, src_pocket_path, temp_pocket_path, mode=mode)
    
    if not success:
        with output_log: print("Failed calculating defaults.")
        center = [0.0, 0.0, 0.0]
        size = [20.0, 20.0, 20.0]

    cx, cy, cz = center
    sx, sy, sz = size

    # 2. Build Control Rows (Reset controls)
    cx_w, cx_box = create_control_group("Center X", cx, 0.5, "#FF0000")
    sx_w, sx_box = create_control_group("Size X  ", sx, 1.0, "#FF0000")
    
    cy_w, cy_box = create_control_group("Center Y", cy, 0.5, "#00AA00")
    sy_w, sy_box = create_control_group("Size Y  ", sy, 1.0, "#00AA00")
    
    cz_w, cz_box = create_control_group("Center Z", cz, 0.5, "#0000FF")
    sz_w, sz_box = create_control_group("Size Z  ", sz, 1.0, "#0000FF")

    controls['cx'] = cx_w; controls['sx'] = sx_w
    controls['cy'] = cy_w; controls['sy'] = sy_w
    controls['cz'] = cz_w; controls['sz'] = sz_w
    
    # Header
    header = widgets.HTML("<h3>Manual Gridbox Adjustment</h3>")
    
    row_x = widgets.HBox([cx_box, widgets.HTML("&nbsp;&nbsp;|&nbsp;&nbsp;"), sx_box])
    row_y = widgets.HBox([cy_box, widgets.HTML("&nbsp;&nbsp;|&nbsp;&nbsp;"), sy_box])
    row_z = widgets.HBox([cz_box, widgets.HTML("&nbsp;&nbsp;|&nbsp;&nbsp;"), sz_box])
    
    ui_container = widgets.VBox([
        header,
        widgets.HTML("<hr style='border-top: 1px solid #ccc;'>"),
        row_x,
        row_y,
        row_z,
        widgets.HTML("<hr style='border-top: 1px solid #ccc;'>"),
        confirm_btn
    ])
    
    # Display UI
    with output_log:
        clear_output()
        display(ui_container)
    
    # Trigger first viz
    update_visual(None, create_new=True)


def update_visual(b, create_new=False):
    viz_output.clear_output(wait=True)
    with viz_output:
        if 'cx' not in controls: return
        
        # Get Values
        cx = controls['cx'].value; sx = controls['sx'].value
        cy = controls['cy'].value; sy = controls['sy'].value
        cz = controls['cz'].value; sz = controls['sz'].value
        
        # Setup View
        view = py3Dmol.view(width=800, height=600)
        
        # 0. Receptor Surface (Background)
        # It's helpful to see the receptor too
        if os.path.exists(receptor_filename):
            with open(receptor_filename, 'r') as f:
                view.addModel(f.read(), "pdb")
            view.setStyle({'model': -1}, {}) # Hide atoms
            view.addSurface(py3Dmol.SES, {'opacity': 0.3, 'color': 'gray'}, {'model': -1})

        # 1. Expanded Pocket (Cartoon) - The one with 3A buffer
        temp_pocket_path = os.path.join(pockets_dir, "temp_calc.pdb")
        if os.path.exists(temp_pocket_path):
            with open(temp_pocket_path, 'r') as f:
                view.addModel(f.read(), "pdb")
            # Cartoon Representation
            view.setStyle({'model': -1}, {'cartoon': {'color': 'cyan', 'opacity': 0.8}})
        
        # 2. Original Fpocket (White Surface)
        selected_pocket = pocket_dropdown.value
        src_pocket_path = os.path.join(fpocket_storage_dir, selected_pocket)
        if os.path.exists(src_pocket_path):
             with open(src_pocket_path, 'r') as f:
                view.addModel(f.read(), "pdb")
             # White Surface
             view.addSurface(py3Dmol.SES, {'opacity': 0.5, 'color': 'white'}, {'model': -1})
             # Hide atoms of this one, just surface
             view.setStyle({'model': -1}, {}) 
             
             # Calculate Volume
             vol = get_pocket_volume_hull(src_pocket_path)
             print(f"Original Pocket Volume (Convex Hull): {vol:.2f} A^3")
             
             # 3. Volume Representation (Red Surface)
             view.addSurface(py3Dmol.SAS, {'opacity': 0.3, 'color': 'red'}, {'model': -1})

        # 4. Gridbox
        view.addBox({
            'center': {'x': cx, 'y': cy, 'z': cz},
            'dimensions': {'w': sx, 'h': sy, 'd': sz},
            'color': 'cyan',
            'opacity': 0.2,
            'wireframe': False
        })
        view.addBox({
            'center': {'x': cx, 'y': cy, 'z': cz},
            'dimensions': {'w': sx, 'h': sy, 'd': sz},
            'color': 'black',
            'wireframe': True
        })

        if create_new:
            view.zoomTo()
            
        view.show()

def finalize_process(b):
    output_log.clear_output()
    with output_log:
        if 'cx' not in controls: return
        print("Finalizing extraction with custom box...")
        
        final_center = [controls['cx'].value, controls['cy'].value, controls['cz'].value]
        final_size   = [controls['sx'].value, controls['sy'].value, controls['sz'].value]
        
        selected_pocket = pocket_dropdown.value
        src_pocket_path = os.path.join(fpocket_storage_dir, selected_pocket)
        final_pocket_name = "pocket.pdb"
        final_pocket_path = os.path.join(pockets_dir, final_pocket_name)
        
        mode = "extract" if mode_selector.value == "Auto Detect" else "direct"
        
        # Determine specific 3A logic? 
        # The logic is embedded in 'extract' mode defaults in process_and_box
        
        _, _, success = run_processing_isolated(
            receptor_filename, src_pocket_path, final_pocket_path, mode=mode, buffer=0.0
        )
        
        if success:
            print(f"\u2705 Box Center: {final_center}")
            print(f"\u2705 Box Size:   {final_size}")
            
            # Save State
            global box_center, box_size, extracted_pocket_path
            box_center = final_center
            box_size = final_size
            extracted_pocket_path = os.path.abspath(final_pocket_path)
            
            save_pipeline_state({
                "box_center": box_center,
                "box_size": box_size,
                "extracted_pocket_path": extracted_pocket_path
            })
            
            # Save to root for visibility/Step 2
            root_pocket = os.path.join(initial_path, "pocket.pdb")
            try:
                shutil.copy(final_pocket_path, root_pocket)
                print(f"\u2705 Ready! (Copied to {root_pocket})")
            except: pass
        else:
            print("Error saving final pocket.")

# Main Buttons for Box Section
init_btn = widgets.Button(description='Reset/Refresh View', button_style='primary', icon='refresh', layout=widgets.Layout(width='200px'))
confirm_btn = widgets.Button(description='Confirm & Extract', button_style='success', icon='check', layout=widgets.Layout(width='100%'))

# init_btn.on_click(initialize_ui) # initialize_ui is triggered by dropdown change now
confirm_btn.on_click(finalize_process)

print("\n--- Gridbox Adjustment ---")
display(widgets.VBox([viz_output, output_log]))
