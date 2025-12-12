#@title 2. Structure-Guided Peptide Generation
#@markdown **Instructions:**
#@markdown 1. Configure parameters (optional).
#@markdown 2. Click **Run Generation** to start the FrankPEPstein pipeline.
#@markdown 3. Use **Stop & Reset** to cancel and clean up if needed.

import os
import sys
import json
import time
import threading
import subprocess
import glob
import ipywidgets as widgets
from IPython.display import display
import multiprocessing

# --- Dependency Check ---
try:
    import py3Dmol
except ImportError:
    subprocess.run("pip install -q py3dmol", shell=True, check=True)
    import py3Dmol

# --- Configuration & State ---
initial_path = os.getcwd()
repo_folder = os.path.join(initial_path, "FrankPEPstein")
state_file = "pipeline_state.json"

# Fixed DB Paths
minipockets_folder = os.path.join(initial_path, "DB/minipockets_surface80_winsize3_size3_curated-db")
db_folder = os.path.join(initial_path, "DB/filtered_DB_P5-15_R30_id10")

# Executables permission fix
def fix_permissions():
    executables = [
        f"{repo_folder}/utilities/vina_1.2.4_linux_x86_64",
        f"{initial_path}/utilities/click/click"
    ]
    for exe in executables:
        if os.path.exists(exe):
            os.chmod(exe, 0o755)

fix_permissions()

# --- State Loading ---
pipeline_state = {}
if os.path.exists(state_file):
    try:
        with open(state_file, "r") as f:
            pipeline_state = json.load(f)
    except:
        pass

receptor_path = pipeline_state.get("receptor_filename", None)
# extracted_pocket_path = pipeline_state.get("extracted_pocket_path", None)
# Enforce standard location as per Step 1 refactor
pockets_dir = os.path.join(initial_path, "pockets")
standard_pocket_path = os.path.join(pockets_dir, "pocket.pdb")

if os.path.exists(standard_pocket_path):
    extracted_pocket_path = standard_pocket_path
else:
    extracted_pocket_path = pipeline_state.get("extracted_pocket_path", None)

box_center = pipeline_state.get("box_center", None)
box_size = pipeline_state.get("box_size", None)

if not receptor_path or not extracted_pocket_path:
    print("⚠️ Warning: Receptor or Extracted Pocket not found in state. Please run Step 1 first.")

# --- Widgets ---
style = {'description_width': 'initial'}

w_pep_size = widgets.IntSlider(value=8, min=5, max=15, description='Peptide Size:', style=style)
w_threads = widgets.IntText(value=multiprocessing.cpu_count(), description='Threads:', style=style)
w_candidates = widgets.IntText(value=10, description='Candidates:', style=style)

btn_run = widgets.Button(description='Run Generation', button_style='success', icon='play')
btn_stop = widgets.Button(description='Stop & Reset', button_style='danger', icon='stop')

out_log = widgets.Output(layout={'border': '1px solid black', 'height': '200px', 'overflow_y': 'scroll'})
out_vis = widgets.Output()

# --- Execution Logic ---
process_handle = None
stop_event = threading.Event()

def log(msg):
    with out_log:
        print(msg)

def run_pipeline(pep_size, threads, candidates):
    global process_handle
    stop_event.clear()
    
    # Define Output Paths (for monitoring)
    run_folder_name = "FrankPEPstein_run"
    super_out_name = f"superpockets_residuesAligned3_RMSD0.1"
    output_superposer_path = os.path.join(initial_path, run_folder_name, super_out_name)
    
    # Construct CLI Command for run_FrankPEPstein.py
    script_path = os.path.join(repo_folder, "scripts/run_FrankPEPstein.py")
    
    cmd_list = [
        sys.executable, script_path,
        "-w", str(pep_size),
        "-t", str(threads),
        "-c", str(candidates)
    ]
    
    # Append Gridbox coordinates from Step 1
    if box_center and box_size:
        cmd_list.extend([
            "-xc", str(box_center[0]),
            "-yc", str(box_center[1]),
            "-zc", str(box_center[2]),
            "-xs", str(box_size[0]),
            "-ys", str(box_size[1]),
            "-zs", str(box_size[2])
        ])
    else:
        log("⚠️ Gridbox definition missing in state. Execution will fail.")
        return

    log(f"--- Starting Pipeline ---")
    log(f"Command: {' '.join(cmd_list)}")
    
    try:
        process_handle = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=initial_path)
        
        # Monitor Loop
        while True:
            ret = process_handle.poll()
            if stop_event.is_set():
                process_handle.terminate()
                log("🛑 Process stopped by user.")
                return
            
            if ret is not None:
                # Finished
                break
                
            # Optional: Read output line by line to show progress in log?
            # For now, simplistic monitoring.
            time.sleep(1)
            
        if ret != 0:
            stderr = process_handle.stderr.read()
            log(f"❌ Pipeline failed: {stderr}")
            return
            
        log("✅ Pipeline Finished Successfully.")
        
    except Exception as e:
        log(f"❌ Error executing pipeline: {e}")
        return

def on_run_click(b):
    out_log.clear_output()
    out_vis.clear_output()
    btn_run.disabled = True
    
    # Start Visualization Thread
    threading.Thread(target=viz_loop, daemon=True).start()
    
    # Initial Visualization (Static) - Immediate Feedback
    with out_vis:
        out_vis.clear_output(wait=True)
        view = py3Dmol.view(width=800, height=600)
        
        # Receptor
        if receptor_path and os.path.exists(receptor_path):
            with open(receptor_path, 'r') as f:
                view.addModel(f.read(), "pdb")
            view.setStyle({'model': -1}, {})
            view.addSurface(py3Dmol.SES, {'opacity': 0.9, 'color': 'white'})
        
        # Pocket (Orange)
        if extracted_pocket_path and os.path.exists(extracted_pocket_path):
             with open(extracted_pocket_path, 'r') as f:
                view.addModel(f.read(), "pdb")
             view.setStyle({'model': -1}, {'sphere': {'color': 'orange', 'opacity': 0.5}})
        
        # GridBox (Red)
        if box_center and box_size:
             cx, cy, cz = box_center
             sx, sy, sz = box_size
             view.addBox({
                 'center': {'x': cx, 'y': cy, 'z': cz},
                 'dimensions': {'w': sx, 'h': sy, 'd': sz},
                 'color': 'red',
                 'opacity': 0.5
             })
             
        view.zoomTo()
        view.show()
        print("Scanning for fragments...")
        
    # Run Pipeline in Thread
    def target():
        run_pipeline(w_pep_size.value, w_threads.value, w_candidates.value)
        btn_run.disabled = False
        
    threading.Thread(target=target).start()

def on_stop_click(b):
    stop_event.set()
    if process_handle:
        process_handle.terminate()
    
    # Cleanup logic requested by user
    run_folder_name = "FrankPEPstein_run"
    output_superposer_path = os.path.join(initial_path, run_folder_name, f"superpockets_residuesAligned3_RMSD0.1")
    
    if os.path.exists(output_superposer_path):
        log(f"Cleaning up {output_superposer_path}...")
        os.system(f"rm -rf {output_superposer_path}")
        
    btn_run.disabled = False
    log("Stopped and Reset.")

btn_run.on_click(on_run_click)
btn_stop.on_click(on_stop_click)

# --- Visualization Logic ---
def viz_loop():
    # Monitor output folder
    run_folder_name = "FrankPEPstein_run"
    super_out_name = f"superpockets_residuesAligned3_RMSD0.1"
    target_dir = os.path.join(initial_path, run_folder_name, super_out_name)
    
    last_count = 0
    
    # Wait for directory to exist
    while not os.path.exists(target_dir):
        if stop_event.is_set() or (btn_run.disabled == False): return # Stop if finished or cancelled
        time.sleep(1)
        
    log("Visualization started monitoring...")
    
    while True:
        if stop_event.is_set() or (btn_run.disabled == False): break
        
        files = glob.glob(os.path.join(target_dir, "patch_file_*.pdb"))
        count = len(files)
        
        # Update every 10 new fragments
        if count >= last_count + 10:
            with out_vis:
                out_vis.clear_output(wait=True)
                
                view = py3Dmol.view(width=800, height=600)
                
                # Receptor
                if receptor_path and os.path.exists(receptor_path):
                    with open(receptor_path, 'r') as f:
                        view.addModel(f.read(), "pdb")
                    view.setStyle({'model': -1}, {})
                    view.addSurface(py3Dmol.SES, {'opacity': 0.9, 'color': 'white'}) # Opacity 0.9 as requested
                
                # Pocket (Orange) if extracted exists
                # Or use the one from state? extracted_pocket_path
                if extracted_pocket_path and os.path.exists(extracted_pocket_path):
                     with open(extracted_pocket_path, 'r') as f:
                        view.addModel(f.read(), "pdb")
                     view.setStyle({'model': -1}, {'sphere': {'color': 'orange', 'opacity': 0.5}})
                
                # GridBox (Red)
                if box_center and box_size:
                     # Add box representation using shape
                     # py3Dmol shape spec: {vertex: [...], ...}
                     # Easier: add box wireframe
                     
                     # Calculate corners
                     cx, cy, cz = box_center
                     sx, sy, sz = box_size
                     
                     min_x, max_x = cx - sx/2, cx + sx/2
                     min_y, max_y = cy - sy/2, cy + sy/2
                     min_z, max_z = cz - sz/2, cz + sz/2
                     
                     view.addBox({
                         'center': {'x': cx, 'y': cy, 'z': cz},
                         'dimensions': {'w': sx, 'h': sy, 'd': sz},
                         'color': 'red',
                         'opacity': 0.5
                     })

                # Fragments
                # Load a subset if too many? or all?
                # Loading hundreds of PDBs into py3Dmol via string might slow down browser.
                # Let's load the latest 50 or so.
                sorted_files = sorted(files, key=os.path.getmtime, reverse=True)[:50]
                
                for pf in sorted_files:
                    with open(pf, 'r') as f:
                        view.addModel(f.read(), "pdb")
                    view.setStyle({'model': -1}, {'stick': {'colorscheme': 'greenCarbon'}})
                
                view.zoomTo()
                view.show()
                print(f"Fragments found: {count} (Showing last {len(sorted_files)})")
                
            last_count = count
            
        time.sleep(2)

# --- Layout ---
ui = widgets.VBox([
    widgets.HBox([w_pep_size, w_threads, w_candidates]),
    widgets.HBox([btn_run, btn_stop]),
    out_log,
    out_vis
])

display(ui)
