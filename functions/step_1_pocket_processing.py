#@title 1. Pocket Selection & Gridbox Generation
#@markdown **Instructions:**
#@markdown 1. Upload **Receptor**.
#@markdown 2. Select **Detection Mode**.
#@markdown 3. **Auto**: Run fpocket, select candidate (Red). **Manual**: Upload pocket (Red).
#@markdown 4. Click **Confirm Selection & Calculate Box** to generate buffer and gridbox.
#@markdown 5. Adjust Gridbox (RGB Wireframe) and click **Finalize Extraction**.

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
initial_path = os.getcwd()
pockets_dir = os.path.join(initial_path, "pockets")
fpocket_storage_dir = os.path.join(initial_path, "fpocket_pockets")
os.makedirs(pockets_dir, exist_ok=True)
os.makedirs(fpocket_storage_dir, exist_ok=True)

receptor_filename = "receptor.pdb"
pipeline_state_file = "pipeline_state.json"
fpocket_exe_path = "/usr/local/envs/FrankPEPstein/bin/fpocket" # Forced path

import json
def save_pipeline_state(data):
    current = {}
    if os.path.exists(pipeline_state_file):
        try:
            with open(pipeline_state_file, 'r') as f: current = json.load(f)
        except: pass
    current.update(data)
    with open(pipeline_state_file, 'w') as f: json.dump(current, f)

# --- Helpers ---

def get_uploaded_file_data(widget_value):
    """Safe extraction for ipywidgets 7 and 8 with sanitization"""
    if not widget_value: return None, None
    
    name, content = None, None
    
    # Version 7: value is dict {filename: {content:..., metadata:...}}
    if isinstance(widget_value, dict):
        filename = list(widget_value.keys())[0]
        file_data = widget_value[filename]
        content = file_data['content']
        if 'metadata' in file_data and 'name' in file_data['metadata']:
             name = file_data['metadata']['name']
        else:
             name = filename
        
    # Version 8: value is list [{'name':..., 'content':...}]
    elif isinstance(widget_value, list):
        file_data = widget_value[0]
        name = file_data['name']
        content = file_data['content']
        
    # Sanitize Name (remove " (1)", etc.)
    if name:
        name = re.sub(r'\s*\(\d+\)', '', name)
        
    return name, content

def get_pocket_center_and_vol(pdb_file):
    """Returns center [x,y,z] and volume for a pocket PDB."""
    script_content = f"""
import sys
import numpy as np
from Bio.PDB import PDBParser
from scipy.spatial import ConvexHull

def analyze(pdb_file):
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("struct", pdb_file)
        points = []
        for atom in structure.get_atoms():
             points.append(atom.get_coord())
             
        if not points:
            print("CENTER:0,0,0")
            print("VOL:0.0")
            return

        coords = np.array(points)
        center = np.mean(coords, axis=0)
        print(f"CENTER:{{center[0]}},{{center[1]}},{{center[2]}}")
            
        if len(points) < 4:
            print("VOL:0.0")
        else:
            hull = ConvexHull(points)
            print(f"VOL:{{hull.volume}}")
            
    except Exception as e:
        print("CENTER:0,0,0")
        print("VOL:0.0")

if __name__ == "__main__":
    analyze(sys.argv[1])
"""
    script_name = "calc_metrics_isolated.py"
    with open(script_name, "w") as f: f.write(script_content)
        
    python_exe = "/usr/local/envs/FrankPEPstein/bin/python"
    if not os.path.exists(python_exe): python_exe = sys.executable

    center = [0,0,0]
    vol = 0.0
    try:
        result = subprocess.run([python_exe, script_name, pdb_file], capture_output=True, text=True)
        for line in result.stdout.splitlines():
            if line.startswith("CENTER:"):
                center = [float(x) for x in line.split(":")[1].split(",")]
            elif line.startswith("VOL:"):
                vol = float(line.split(":")[1])
    except: pass
    return center, vol

# --- UI Components ---

# 1. Receptor Output
receptor_upload_widget = widgets.FileUpload(description="Upload Receptor (pdb)", accept=".pdb", multiple=False, layout=widgets.Layout(width='300px'))
receptor_status = widgets.Output()

# 2. Mode Selection
mode_selector = widgets.ToggleButtons(options=['Auto Detect', 'Manual Upload'], description='Mode:', button_style='')

# 3. Auto Mode Widgets
run_fpocket_btn = widgets.Button(description="Run fpocket", button_style='info', icon='play')
pocket_dropdown = widgets.Dropdown(description="Select Pocket:", options=[], disabled=True)

# 4. Manual Mode Widgets
manual_upload_btn = widgets.FileUpload(description="Upload Pocket PDB", accept=".pdb", multiple=False, layout=widgets.Layout(width='300px'))

# 5. Shared Action
calc_box_btn = widgets.Button(description="Confirm Selection & Calculate Box", button_style='warning', icon='calculator', layout=widgets.Layout(width='100%'), disabled=True)

# 6. Gridbox Controls (Starts Hidden)
viz_output = widgets.Output()
log_output = widgets.Output()

# --- Logic: Part 1 (Selection) ---

def handle_receptor_upload(change):
    receptor_status.clear_output()
    with receptor_status:
        if not receptor_upload_widget.value: return
        name, content = get_uploaded_file_data(receptor_upload_widget.value)
        if content:
            with open(receptor_filename, "wb") as f: f.write(content)
            print("✅ Receptor uploaded successfully.")
            update_viz_phase1()

def run_fpocket_action(b):
    log_output.clear_output()
    with log_output:
        if not os.path.exists(receptor_filename):
            print("❌ Receptor not found!")
            return
        print("Running fpocket (forced conda path)...")
        subprocess.run(f"rm -rf {fpocket_storage_dir}/*", shell=True)
        
        # Force execution
        cmd = [fpocket_exe_path, "-f", receptor_filename]
        if not os.path.exists(fpocket_exe_path):
             print(f"⚠️ {fpocket_exe_path} not found. Trying 'fpocket' from PATH.")
             cmd = ["fpocket", "-f", receptor_filename]
             
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
             print(f"❌ fpocket failed:\n{res.stderr}")
             return

        # Find results
        base = os.path.splitext(receptor_filename)[0]
        out_dir = base + "_out"
        pockets = glob.glob(os.path.join(out_dir, "pockets", "pocket*_atm.pdb"))
        
        if pockets:
            print(f"✅ Found {len(pockets)} pockets.")
            options = []
            for p in pockets:
                basename = os.path.basename(p)
                dest = os.path.join(fpocket_storage_dir, basename)
                shutil.copy(p, dest)
                options.append(basename)
            
            # Sort by number for labels
            options.sort(key=lambda x: int(re.search(r'\d+', x).group()) if re.search(r'\d+', x) else 0)
            
            pocket_dropdown.options = options
            pocket_dropdown.disabled = False
            if options: pocket_dropdown.value = options[0]
            calc_box_btn.disabled = False
            update_viz_phase1()
        else:
            print("❌ No pockets found.")

def handle_manual_upload(change):
    log_output.clear_output()
    with log_output:
        if not manual_upload_btn.value: return
        name, content = get_uploaded_file_data(manual_upload_btn.value)
        if not name or not content: return
        
        dest = os.path.join(fpocket_storage_dir, name)
        with open(dest, "wb") as f: f.write(content)
        print(f"✅ Uploaded {name}")
        
        pocket_dropdown.options = [name] # Use dropdown to hold state even if hidden
        pocket_dropdown.value = name
        calc_box_btn.disabled = False
        update_viz_phase1()

# Visualization Phase 1: Selection
def update_viz_phase1(change=None):
    viz_output.clear_output(wait=True)
    with viz_output:
        view = py3Dmol.view(width=800, height=600)
        
        # 1. Receptor: White Surface
        if os.path.exists(receptor_filename):
            with open(receptor_filename, 'r') as f: view.addModel(f.read(), "pdb")
            view.setStyle({'model': -1}, {})
            view.addSurface(py3Dmol.SES, {'opacity': 0.7, 'color': 'white'}, {'model': -1})

        # 2. Pockets
        if mode_selector.value == 'Auto Detect':
            # Rainbow Colors for Non-Selected
            colors = ['#FF00FF', '#00FFFF', '#FFFF00', '#00FF00', '#0000FF', '#FFA500']
            
            all_pockets = pocket_dropdown.options
            selected = pocket_dropdown.value
            
            for i, p_name in enumerate(all_pockets):
                path = os.path.join(fpocket_storage_dir, p_name)
                is_selected = (p_name == selected)
                
                with open(path, 'r') as f: view.addModel(f.read(), "pdb")
                model_idx = -1
                view.setStyle({'model': -1}, {}) # Hide atoms
                
                # Center for Label
                center, _ = get_pocket_center_and_vol(path)
                # Label ID (extract number)
                pid = re.search(r'\d+', p_name)
                label_text = pid.group() if pid else str(i+1)
                
                view.addLabel(label_text, {'position': {'x':center[0], 'y':center[1], 'z':center[2]}, 
                                           'backgroundColor': 'black', 'fontColor': 'white'})

                if is_selected:
                    # Red Surface, 1.0 Opacity (No transparency)
                    view.addSurface(py3Dmol.SES, {'opacity': 1.0, 'color': 'red'}, {'model': model_idx})
                else:
                    # Rainbow Surface, 0.9 Opacity (10% transp)
                    color = colors[i % len(colors)]
                    view.addSurface(py3Dmol.SES, {'opacity': 0.9, 'color': color}, {'model': model_idx})
                    
        else: # Manual
            # Just show upload if exists
            if pocket_dropdown.value:
                path = os.path.join(fpocket_storage_dir, pocket_dropdown.value)
                if os.path.exists(path):
                    with open(path, 'r') as f: view.addModel(f.read(), "pdb")
                    view.setStyle({'model': -1}, {}) 
                    # Red Surface, 1.0 Opacity
                    view.addSurface(py3Dmol.SES, {'opacity': 1.0, 'color': 'red'}, {'model': -1})

        view.zoomTo()
        view.show()

# --- Logic: Part 2 (Calculation & Gridbox) ---

# Gridbox Controls
cx_w = widgets.FloatText(description='Center X')
cy_w = widgets.FloatText(description='Center Y')
cz_w = widgets.FloatText(description='Center Z')
sx_w = widgets.FloatText(description='Size X')
sy_w = widgets.FloatText(description='Size Y')
sz_w = widgets.FloatText(description='Size Z')
controls = {'cx': cx_w, 'cy': cy_w, 'cz': cz_w, 'sx': sx_w, 'sy': sy_w, 'sz': sz_w}
finalize_btn = widgets.Button(description="Finalize Extraction", button_style='success', icon='check', layout=widgets.Layout(width='100%'))
gridbox_ui = widgets.VBox([
    widgets.HTML("<h3>Gridbox Adjustment</h3>"),
    widgets.HBox([cx_w, sx_w]), widgets.HBox([cy_w, sy_w]), widgets.HBox([cz_w, sz_w]),
    finalize_btn
])
gridbox_ui.layout.display = 'none' # Hidden initially

def calculate_box_and_transition(b):
    log_output.clear_output()
    gridbox_ui.layout.display = 'none' # Hide while working
    
    with log_output:
        if not pocket_dropdown.value: return
        print("Calculating 3A buffer and Gridbox...")
        
        # Determine paths
        selected_pocket = pocket_dropdown.value
        src_path = os.path.join(fpocket_storage_dir, selected_pocket)
        # Temp output for buffered pocket
        temp_pocket_path = os.path.join(pockets_dir, "temp_calc.pdb")
        
        # ALWAYS EXTRACT (to get buffer)
        # We need the extraction script inline here or reuse helper?
        # Reusing the existing extraction script logic is best, but defining it here for self-containedness
        extract_script = f"""
import sys, os
from Bio.PDB import PDBParser, PDBIO, Select, NeighborSearch
def run(receptor, pocket, out):
    try:
        parser = PDBParser(QUIET=True)
        rec = parser.get_structure("r", receptor)
        poc = parser.get_structure("p", pocket)
        
        # 1. Select neighbors (3A buffer logic? user said "buffer de 3A" but typically script used 5A in NeighborSearch? using 5A to be safe/standard or user defined)
        # Original script used 5.0 for neighbors. User asked for "buffer de 3A".
        # Let's use 3.0 if specifically requested, or sticking to 5.0 if that was 'standard'.
        # User prompt: "Boton de seleccionar pocket. se añade un buffer de 3A." -> Use 3.0
        
        atoms_p = list(poc.get_atoms())
        atoms_r = list(rec.get_atoms())
        ns = NeighborSearch(atoms_r)
        
        selected = set()
        for a in atoms_p:
            nearby = ns.search(a.get_coord(), 3.0, level='R') # 3A Buffer
            for res in nearby: selected.add((res.parent.id, res.id))
            
        class SelectP(Select):
            def accept_residue(self, residue): return (residue.parent.id, residue.id) in selected
            
        io = PDBIO()
        io.set_structure(rec)
        io.save(out, SelectP())
        
        # 2. Convert to chain P 
        st = parser.get_structure("clean", out)
        coords = []
        for model in st:
            for chain in model:
                chain.id = 'p'
                for residue in chain:
                    for atom in residue:
                        coords.append(atom.get_coord())
        io.set_structure(st)
        io.save(out)
        
        # 3. Print Box Stats
        if coords:
            min_c = [min([c[i] for c in coords]) for i in range(3)]
            max_c = [max([c[i] for c in coords]) for i in range(3)]
            center = [(min_c[i] + max_c[i])/2 for i in range(3)]
            size = [(max_c[i] - min_c[i]) for i in range(3)] # Box fits exactly? Or + buffer? Default usually fits. User can adjust.
            print(f"BOX:{{center[0]}},{{center[1]}},{{center[2]}}|{{size[0]}},{{size[1]}},{{size[2]}}")
            
    except Exception as e: print(e)

if __name__ == "__main__": run(sys.argv[1], sys.argv[2], sys.argv[3])
"""
        script_path = "extract_temp.py"
        with open(script_path, "w") as f: f.write(extract_script)
        
        python_exe = "/usr/local/envs/FrankPEPstein/bin/python"
        if not os.path.exists(python_exe): python_exe = sys.executable
        
        res = subprocess.run([python_exe, script_path, receptor_filename, src_path, temp_pocket_path], capture_output=True, text=True)
        
        # Parse Box defaults
        center = [0,0,0]; size = [20,20,20]
        for line in res.stdout.splitlines():
            if line.startswith("BOX:"):
                parts = line.split("|")
                center = [float(x) for x in parts[0].split(":")[1].split(",")]
                size = [float(x) for x in parts[1].split(",")]
        
        # Set Controls
        cx_w.value = center[0]; cy_w.value = center[1]; cz_w.value = center[2]
        sx_w.value = size[0]; sy_w.value = size[1]; sz_w.value = size[2]
        
        # Show Controls & Switch Viz
        gridbox_ui.layout.display = 'block'
        print("✅ Processing complete. Adjust Gridbox.")
        update_viz_phase2()

# Visualization Phase 2: Gridbox
def update_viz_phase2(change=None):
    viz_output.clear_output(wait=True)
    with viz_output:
        view = py3Dmol.view(width=800, height=600)
        
        # RECEPTOR REMOVED (User Request)
        
        temp_pocket_path = os.path.join(pockets_dir, "temp_calc.pdb")
        selected_pocket = pocket_dropdown.value
        src_path = os.path.join(fpocket_storage_dir, selected_pocket)
        
        mode = mode_selector.value
        
        # 1. Visuals based on Mode
        if mode == 'Auto Detect':
            # Auto: Buffered=Cartoon(0.7), Original=Surface(1.0), Volume=Spheres(1.0)
            if os.path.exists(temp_pocket_path):
                with open(temp_pocket_path, 'r') as f: view.addModel(f.read(), "pdb")
                view.setStyle({'model': -1}, {'cartoon': {'color': 'cyan', 'opacity': 0.7}})
            
            if os.path.exists(src_path):
                with open(src_path, 'r') as f: view.addModel(f.read(), "pdb")
                view.setStyle({'model': -1}, {})
                view.addSurface(py3Dmol.SES, {'opacity': 0.9, 'color': 'white'}, {'model': -1})
                # Volume (using SAS roughly or spheres?)
                # user said "volumen calculado dentro de la gridbox en rojo esferas" -> Sphere style 
                # Ideally we calculate spheres, but for now we can use VDW/Sphere representation of atoms
                view.addStyle({'model': -1}, {'sphere': {'color': 'red', 'opacity': 0.9, 'radius': 1.0}}) 
                
        else: # Manual
            # Manual: Buffered=Surface(0.9), Volume=Spheres(0.9)
            if os.path.exists(temp_pocket_path):
                with open(temp_pocket_path, 'r') as f: view.addModel(f.read(), "pdb")
                view.setStyle({'model': -1}, {})
                view.addSurface(py3Dmol.SES, {'opacity': 0.9, 'color': 'white'}, {'model': -1})
                
            # Volume? If manual upload is just the atoms, show them as spheres too?
            # Or assume the "buffer" is the surface and the "upload" is the volume?
            # Usually Manual Upload = Pocket Definition (Volume).
            if os.path.exists(src_path):
                 with open(src_path, 'r') as f: view.addModel(f.read(), "pdb")
                 view.setStyle({'model': -1}, {'sphere': {'color': 'red', 'opacity': 0.9, 'radius': 1.0}})

        # 2. Gridbox (RGB Cylinders)
        cx, cy, cz = cx_w.value, cy_w.value, cz_w.value
        sx, sy, sz = sx_w.value, sy_w.value, sz_w.value
        
        x1, x2 = cx - sx/2, cx + sx/2
        y1, y2 = cy - sy/2, cy + sy/2
        z1, z2 = cz - sz/2, cz + sz/2
        r = 0.2
        
        # X (Red)
        for y, z in [(y1,z1), (y2,z1), (y1,z2), (y2,z2)]:
            view.addCylinder({'start':{'x':x1,'y':y,'z':z}, 'end':{'x':x2,'y':y,'z':z}, 'radius':r, 'color':'red'})
        # Y (Green)
        for x, z in [(x1,z1), (x2,z1), (x1,z2), (x2,z2)]:
            view.addCylinder({'start':{'x':x,'y':y1,'z':z}, 'end':{'x':x,'y':y2,'z':z}, 'radius':r, 'color':'green'})
        # Z (Blue)
        for x, y in [(x1,y1), (x2,y1), (x1,y2), (x2,y2)]:
            view.addCylinder({'start':{'x':x,'y':y,'z':z1}, 'end':{'x':x,'y':y,'z':z2}, 'radius':r, 'color':'blue'})

        view.zoomTo()
        view.show()

def finalize_process(b):
    log_output.clear_output()
    with log_output:
        print("Finalizing...")
        final_pocket_path = os.path.join(pockets_dir, "pocket.pdb")
        temp_pocket_path = os.path.join(pockets_dir, "temp_calc.pdb")
        if os.path.exists(temp_pocket_path):
            shutil.copy(temp_pocket_path, final_pocket_path)
            shutil.copy(final_pocket_path, os.path.join(initial_path, "pocket.pdb"))
            
            save_pipeline_state({
                "box_center": [cx_w.value, cy_w.value, cz_w.value],
                "box_size": [sx_w.value, sy_w.value, sz_w.value],
                "extracted_pocket_path": os.path.abspath(final_pocket_path)
            })
            print("✅ Pocket and Gridbox Saved!")
        else:
            print("❌ Error: No processed pocket found.")

# Link Events
receptor_upload_widget.observe(handle_receptor_upload, names='value')
run_fpocket_btn.on_click(run_fpocket_action)
manual_upload_btn.observe(handle_manual_upload, names='value')
pocket_dropdown.observe(lambda c: update_viz_phase1() if c['type'] == 'change' and c['name'] == 'value' else None)
mode_selector.observe(lambda c: update_viz_phase1() if c['type'] == 'change' and c['name'] == 'value' else None)

calc_box_btn.on_click(calculate_box_and_transition)
finalize_btn.on_click(finalize_process)

# Live Gridbox Updates
for w in controls.values():
    w.observe(lambda c: update_viz_phase2() if c['type'] == 'change' and c['name'] == 'value' else None)

# Update Visibility
def update_ui_mode(c):
    is_auto = (mode_selector.value == 'Auto Detect')
    run_fpocket_btn.layout.display = 'block' if is_auto else 'none'
    pocket_dropdown.layout.display = 'block' if is_auto else 'none'
    manual_upload_btn.layout.display = 'none' if is_auto else 'block'
    update_viz_phase1()

mode_selector.observe(update_ui_mode, names='value')

print("Step 1.1: Load Receptor")
display(widgets.VBox([receptor_upload_widget, receptor_status]))

print("\nStep 1.2: Pocket Selection")
display(widgets.VBox([mode_selector, widgets.HBox([run_fpocket_btn, pocket_dropdown]), manual_upload_btn]))

print("\nStep 1.3: Calculation & Viz")
display(calc_box_btn)
display(widgets.VBox([viz_output, gridbox_ui, log_output]))

# Init
update_ui_mode(None)
