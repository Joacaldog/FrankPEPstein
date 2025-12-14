#@title 1. Input & Pocket Selection
#@markdown **Instructions:**
#@markdown 1. Upload your Receptor PDB.
#@markdown 2. Choose Mode: **Auto Detect** (runs fpocket) or **Manual Upload** (upload your specific pocket PDB).
#@markdown 3. Select the pocket from the dropdown to visualize.

import os
import sys
import subprocess
try:
    import py3Dmol
except ImportError:
    # Try adding FrankPEPstein env to path
    env_path = "/usr/local/envs/FrankPEPstein"
    site_packages = f"{env_path}/lib/python3.10/site-packages"
    if os.path.exists(site_packages):
        if site_packages not in sys.path:
            sys.path.append(site_packages)
        # Add binary path too
        if f"{env_path}/bin" not in os.environ['PATH']:
            os.environ['PATH'] = f"{env_path}/bin:" + os.environ['PATH']
    
    # Retry import
    try:
        import py3Dmol
    except ImportError:
        print("py3Dmol not found. Installing...")
        subprocess.run("pip install -q py3dmol", shell=True, check=True)
        import py3Dmol

import ipywidgets as widgets
from google.colab import files
from IPython.display import display
import shutil
import json
import re

# --- configuration ---
detection_mode = "Auto Detect" #@param ["Auto Detect", "Manual Upload"]

# Global variables
receptor_filename = None
initial_path = os.getcwd() # Main Directory
# Refactor: Use FrankPEPstein_run as centralized storage for execution
pockets_dir = os.path.join(initial_path, "FrankPEPstein_run") 
final_pockets_list = []

# Ensure pockets dir exists
if not os.path.exists(pockets_dir):
    os.makedirs(pockets_dir)

# --- Persistence Logic ---
def save_pipeline_state(updates):
    state_file = "pipeline_state.json"
    current_state = {}
    if os.path.exists(state_file):
        try:
            with open(state_file, "r") as f:
                current_state = json.load(f)
        except:
            pass
    current_state.update(updates)
    with open(state_file, "w") as f:
        json.dump(current_state, f, indent=4)
    print(f"State saved to {state_file}")

# --- 1. Upload Receptor ---
print(f"--- Upload Receptor PDB ({detection_mode}) ---")
uploaded_r = files.upload()

if not uploaded_r:
    print("No receptor file uploaded.")
else:
    raw_filename = list(uploaded_r.keys())[0]
        
    # Check for Colab duplicate naming (e.g. receptor(1).pdb)
    match = re.search(r'^(.*?)\s?\(\d+\)(\.[^.]*)?$', raw_filename)
    if match:
        clean_name = match.group(1) + (match.group(2) if match.group(2) else "")
        print(f"Detected duplicate upload: {raw_filename} -> overwriting {clean_name}")
        
        if os.path.exists(clean_name):
            os.remove(clean_name)
        os.rename(raw_filename, clean_name)
        receptor_filename = os.path.abspath(clean_name)
    else:
        receptor_filename = os.path.abspath(raw_filename)
        
    print(f"Receptor: {receptor_filename}")
    save_pipeline_state({"receptor_filename": receptor_filename})

    # --- 2. Pocket Handling ---

    # Determine fpocket path
    fpocket_bin = "fpocket"
    if shutil.which(fpocket_bin) is None:
        # Try specific env path
        env_fpocket = "/usr/local/envs/FrankPEPstein/bin/fpocket"
        if os.path.exists(env_fpocket):
            fpocket_bin = env_fpocket
        else:
            print("⚠️ fpocket executable not found in PATH or FrankPEPstein env.")
            # We let it fail in subprocess if still not found, but this warning helps.

    if detection_mode == "Auto Detect":
        try:
            print(f"Running fpocket on {receptor_filename}")
            # Capture output for debugging
            # Using -m to filter small pockets as requested - REMOVED due to user report of bugs
            result = subprocess.run(f"{fpocket_bin} -f '{receptor_filename}'", shell=True, capture_output=True, text=True)

            
            if result.returncode != 0:
                print("❌ Error running fpocket.")
                print(f"Exit Code: {result.returncode}")
                print(f"STDERR:\n{result.stderr}")
                print(f"STDOUT:\n{result.stdout}")
                # We can try to look at why.
            else:
                # Success logic
                pass 
                
            # Check for output ONLY if successful or to diagnose
            base_name_no_ext = os.path.splitext(os.path.basename(receptor_filename))[0]
            base_dir = os.path.dirname(receptor_filename)
            
            folder_name_1 = f"{os.path.basename(receptor_filename)}_out"
            folder_name_2 = f"{base_name_no_ext}_out"
            
            possible_folders = [
                os.path.join(base_dir, folder_name_1),
                os.path.join(base_dir, folder_name_2)
            ]
            
            output_folder = next((f for f in possible_folders if os.path.exists(f)), None)

            if output_folder:
                fpocket_pockets_dir = os.path.join(output_folder, "pockets")
                if os.path.exists(fpocket_pockets_dir):
                    # Move/Copy relevant pockets to our centralized dir
                    found_pockets = [f for f in os.listdir(fpocket_pockets_dir) if f.endswith(".pdb")]
                    for p in found_pockets:
                        src = os.path.join(fpocket_pockets_dir, p)
                        dst = os.path.join(pockets_dir, p)
                        shutil.copy(src, dst)
                        final_pockets_list.append(p)
                        
                    print(f"Auto-detection finished. Found {len(final_pockets_list)} pockets.")
                    if not final_pockets_list:
                        print(f"⚠️ No pockets found!)")
                else:
                    print(f"Warning: pockets subdirectory not found in {output_folder}")
            else:
                 if result.returncode == 0:
                     print("Error: fpocket finished but output folder not found.")
                
        except Exception as e:
             print(f"Unexpected error running fpocket: {e}")

    elif detection_mode == "Manual Upload":
        print(f"\n--- Upload Manual Pocket PDB ---")
        uploaded_p = files.upload()
        if uploaded_p:
            for p_file in uploaded_p.keys():
                match = re.search(r'^(.*?)\s?\(\d+\)(\.[^.]*)?$', p_file)
                if match:
                    clean_name = match.group(1) + (match.group(2) if match.group(2) else "")
                    print(f"Detected duplicate upload: {p_file} -> overwriting {clean_name}")
                else:
                    clean_name = p_file

                target_path = os.path.join(pockets_dir, clean_name)
                if os.path.exists(target_path):
                    os.remove(target_path)
                
                os.rename(p_file, target_path)
                
                if clean_name not in final_pockets_list:
                    final_pockets_list.append(clean_name)
                    
            print(f"Manual upload finished. Available pockets: {len(final_pockets_list)}")

    # --- 3. Visualization & Selection ---
    if final_pockets_list:
        print("\n--- Pocket Selection & Visualization ---")
        
        pocket_dropdown = widgets.Dropdown(
            options=sorted(final_pockets_list),
            description='Select Pocket:',
            disabled=False,
        )

        def view_pockets(selected_pocket_file):
            view = py3Dmol.view(width=800, height=600)
            
            # Receptor
            with open(receptor_filename, 'r') as f:
                view.addModel(f.read(), "pdb")
            view.setStyle({}) 
            view.addSurface(py3Dmol.SES, {'opacity': 0.3, 'color': 'white'})
            
            colors = ['#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#FF00FF', '#00FFFF', '#FFA500', '#800080', '#008000', '#800000']
            
            for i, p_file in enumerate(sorted(final_pockets_list)):
                full_path = os.path.join(pockets_dir, p_file)
                if os.path.exists(full_path):
                    with open(full_path, 'r') as f:
                        view.addModel(f.read(), "pdb")
                    
                    match = re.search(r'(\d+)', p_file)
                    label_text = match.group(1) if match else p_file
                    
                    is_selected = (p_file == selected_pocket_file)
                    
                    if is_selected:
                         color = 'red'
                         opacity = 1.0
                         label_style = {'fontSize': 18, 'fontColor': 'red', 'backgroundColor': 'white', 'backgroundOpacity': 0.8, 'border': '2px solid red'}
                    else:
                        color = colors[i % len(colors)]
                        opacity = 0.6
                        label_style = {'fontSize': 12, 'fontColor': 'black', 'backgroundColor': 'white', 'backgroundOpacity': 0.5}

                    view.setStyle({'model': -1}, {'sphere': {'color': color, 'opacity': opacity}})
                    view.addLabel(label_text, label_style, {'model': -1})

            view.zoomTo()
            view.show()
            
        display(widgets.interactive(view_pockets, selected_pocket_file=pocket_dropdown))
    else:
        print("No pockets available to select.")

#@title 4. Pocket Extraction, Manual Adjustment & Box Generation
#@markdown This step calculates the initial box, allows manual adjustment with +/- buttons, and then extracts the final pocket.

import os
import subprocess
import sys
import ipywidgets as widgets
from IPython.display import display, clear_output

# --- Helper Functions (Subprocess) ---
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

# --- Custom Widget Helpers ---

def create_control_group(label, initial_val, step, color_hex):
    """
    Creates a [-] [Value] [+] control group.
    Returns: value_widget, minus_btn, plus_btn, container
    """
    # Style buttons
    btn_layout = widgets.Layout(width='30px')
    # Using 'info' or 'primary' doesn't give custom colors easily in standard buttons without style attributes
    # We will use the container border/label to indicate color.
    
    minus_button = widgets.Button(description='-', layout=btn_layout)
    plus_button = widgets.Button(description='+', layout=btn_layout)
    
    # Custom coloring for buttons is tricky in basic ipywidgets without custom CSS.
    # We'll use the button_style for generic look, but wrap in a colored box.
    # User requested ARROWS because text was invisible.
    
    minus_button = widgets.Button(icon='arrow-left', layout=btn_layout)
    plus_button = widgets.Button(icon='arrow-right', layout=btn_layout)
    
    minus_button.style.button_color = '#e0e0e0'
    plus_button.style.button_color = '#e0e0e0'
    
    def on_minus(b):
        val_widget.value -= step
        update_visual(None)
        
    def on_plus(b):
        val_widget.value += step
        update_visual(None)
            
    minus_button.on_click(on_minus)
    plus_button.on_click(on_plus)
    val_widget.observe(lambda c: update_visual(None) if c['type'] == 'change' and c['name'] == 'value' else None)

    # Label styling
    # We use HTML to color the label text
    label_html = widgets.HTML(f"<b style='color:{color_hex}; font-size:14px; margin-right:5px;'>{label}</b>")
    
    box = widgets.HBox([label_html, minus_button, val_widget, plus_button], layout=widgets.Layout(align_items='center'))
    return val_widget, box

# --- UI State Containers ---
# We store widgets here to access them in callbacks
controls = {}

# --- Logic ---
viz_output = widgets.Output()
output_log = widgets.Output()

def initialize_ui(b):
    output_log.clear_output()
    viz_output.clear_output()
    
    if 'pocket_dropdown' not in globals() or not pocket_dropdown.value:
        with output_log: print("Please select a pocket first.")
        return

    # 1. Calc Initial
    selected_pocket = pocket_dropdown.value
    src_pocket_path = os.path.join(pockets_dir, selected_pocket)
    temp_pocket_path = os.path.join(pockets_dir, "temp_calc.pdb")
    mode = "extract" if detection_mode == "Auto Detect" else "direct"
    
    with output_log: print("Checking pocket parameters...")
    center, size, success = run_processing_isolated(receptor_filename, src_pocket_path, temp_pocket_path, mode=mode)
    
    if not success:
        with output_log: print("Failed calculating defaults.")
        center = [0.0, 0.0, 0.0]
        size = [20.0, 20.0, 20.0]

    cx, cy, cz = center
    sx, sy, sz = size

    # 2. Build Control Rows (Reset controls)
    # X Axis (Red)
    cx_w, cx_box = create_control_group("Center X", cx, 0.5, "#FF0000")
    sx_w, sx_box = create_control_group("Size X  ", sx, 1.0, "#FF0000")
    
    # Y Axis (Green)
    cy_w, cy_box = create_control_group("Center Y", cy, 0.5, "#00AA00")
    sy_w, sy_box = create_control_group("Size Y  ", sy, 1.0, "#00AA00")
    
    # Z Axis (Blue)
    cz_w, cz_box = create_control_group("Center Z", cz, 0.5, "#0000FF")
    sz_w, sz_box = create_control_group("Size Z  ", sz, 1.0, "#0000FF")

    # Store for global access
    controls['cx'] = cx_w; controls['sx'] = sx_w
    controls['cy'] = cy_w; controls['sy'] = sy_w
    controls['cz'] = cz_w; controls['sz'] = sz_w
    
    # Layout
    # Organize by Function (Center Col, Size Col) or by Axis (Row X, Row Y)?
    # User asked for color coding. Row by Axis is cleanest for color grouping.
    
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
    update_visual(None)


def update_visual(b):
    viz_output.clear_output(wait=True)
    with viz_output:
        if 'cx' not in controls: return
        
        # Get Current Values
        cx = controls['cx'].value; sx = controls['sx'].value
        cy = controls['cy'].value; sy = controls['sy'].value
        cz = controls['cz'].value; sz = controls['sz'].value
        
        # Load View
        view = py3Dmol.view(width=800, height=600)
        
        # 1. ADD POCKET ONLY
        # We need the path. Using the one we calculated/detected
        selected_pocket = pocket_dropdown.value
        # If 'Auto', we might be using the raw one. If 'Manual', same.
        # But wait, we want to see the pocket relative to the box.
        # Ideally we use the 'temp_calc.pdb' if it was valid, or the source.
        # Let's use source for visualization to avoid confusion
        src_pocket_path = os.path.join(pockets_dir, selected_pocket)
        if os.path.exists(src_pocket_path):
             with open(src_pocket_path, 'r') as f:
                view.addModel(f.read(), "pdb")
             # Show Atoms nicely
             view.setStyle({'stick':{'colorscheme':'greenCarbon'}})
             view.addSurface(py3Dmol.SES, {'opacity': 0.6, 'color': 'white'})
        
        # 2. ADD GRIDBOX
        # User requested: "caras con colores pero transparente"
        view.addBox({
            'center': {'x': cx, 'y': cy, 'z': cz},
            'dimensions': {'w': sx, 'h': sy, 'd': sz},
            'color': 'cyan',   # Single color for the box itself
            'opacity': 0.4,    # Transparent
            'wireframe': False # Solid faces
        })
        # Add wireframe on top for definition?
        view.addBox({
            'center': {'x': cx, 'y': cy, 'z': cz},
            'dimensions': {'w': sx, 'h': sy, 'd': sz},
            'color': 'black',
            'wireframe': True
        })

        view.zoomTo()
        view.show()

def finalize_process(b):
    output_log.clear_output()
    with output_log:
        if 'cx' not in controls: return
        print("Finalizing extraction with custom box...")
        
        # 1. Use user adjusted values
        final_center = [controls['cx'].value, controls['cy'].value, controls['cz'].value]
        final_size   = [controls['sx'].value, controls['sy'].value, controls['sz'].value]
        
        selected_pocket = pocket_dropdown.value
        src_pocket_path = os.path.join(pockets_dir, selected_pocket)
        final_pocket_name = "pocket.pdb"
        final_pocket_path = os.path.join(pockets_dir, final_pocket_name)
        
        mode = "extract" if detection_mode == "Auto Detect" else "direct"
        
        # Run extraction logic again to ensure clean PDB
        _, _, success = run_processing_isolated(
            receptor_filename, src_pocket_path, final_pocket_path, mode=mode, buffer=0.0
        )
        
        if success:
            print(f"✅ Box Center: {final_center}")
            print(f"✅ Box Size:   {final_size}")
            
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
            
            root_pocket = os.path.join(initial_path, "pocket.pdb")
            try:
                import shutil
                shutil.copy(final_pocket_path, root_pocket)
                print(f"✅ Ready! (Copied to {root_pocket})")
            except: pass
        else:
            print("Error saving final pocket.")

# Main Buttons
init_btn = widgets.Button(description='Start Box Adjustment', button_style='primary', icon='edit', layout=widgets.Layout(width='200px'))
confirm_btn = widgets.Button(description='Confirm & Extract', button_style='success', icon='check', layout=widgets.Layout(width='100%'))

init_btn.on_click(initialize_ui)
confirm_btn.on_click(finalize_process)

print("\n--- 4. Pocket Processing & Box Calculation ---")
display(widgets.VBox([init_btn, viz_output, output_log]))
