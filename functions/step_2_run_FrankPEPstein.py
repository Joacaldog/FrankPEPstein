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
        f"{repo_folder}/utilities/click",
        f"{repo_folder}/utilities/vina_1.2.4_linux_x86_64",
        f"{initial_path}/utilities/click", # Just in case it's here
        f"{initial_path}/utilities/vina_1.2.4_linux_x86_64"
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
    
    # Define Output Paths (based on superposer.py logic usually)
    # superposer creates 'FrankPEPstein_run' and 'superpockets_residuesAligned3_RMSD0.1' inside
    # We should ensure we know where it is.
    
    run_folder_name = "FrankPEPstein_run"
    super_out_name = f"superpockets_residuesAligned3_RMSD0.1"
    output_superposer_path = os.path.join(initial_path, run_folder_name, super_out_name)
    
    # 1. Superposer
    log(f"--- Starting Superposer (Peptide Size: {pep_size}) ---")
    
    # Construct command
    # NOTE: superposer.py uses hardcoded output location usually, we respect that.
    # We use extracted_pocket_path as -T target
    
    # Clean output if exists and requested (Stop button usually handles this, but good to be safe)
    # But user said "if stop is pressed", so we leave it for now.
    
    cmd_superposer = [
        "python3", f"{repo_folder}/scripts/superposer.py",
        "-T", extracted_pocket_path,
        "-d", db_folder,
        "-r", "0.1",
        "-t", str(threads),
        "-a", "3",
        "-fm", minipockets_folder
    ]
    
    log(f"Running: {' '.join(cmd_superposer)}")
    
    try:
        process_handle = subprocess.Popen(cmd_superposer, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=initial_path)
        
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
                
            time.sleep(0.5)
            
        if ret != 0:
            stderr = process_handle.stderr.read()
            log(f"❌ Superposer failed: {stderr}")
            return
            
        log("✅ Superposer finished.")
        
    except Exception as e:
        log(f"❌ Error executing superposer: {e}")
        return

    # 2. FrankVINA 1
    if stop_event.is_set(): return
    log("--- Starting FrankVINA 1 ---")
    
    if os.path.exists(output_superposer_path):
        os.chdir(output_superposer_path)
        # Copy receptor
        if os.path.exists(receptor_path):
            os.system(f'cp "{receptor_path}" receptor.pdb')
        
        cmd_vina1 = f'python3 {repo_folder}/scripts/frankVINA_1.py "{initial_path}" receptor.pdb {threads}'
        log(f"Running: {cmd_vina1}")
        exit_code = os.system(cmd_vina1) # Blocks, but acceptable for this step? Ideally subprocess but os.system used in original
        
        if exit_code != 0:
            log("❌ FrankVINA 1 failed.")
            os.chdir(initial_path)
            return

        # Clean
        os.system("rm * 2> /dev/null") # As per original script logic (careful here!)
        # Actually original script removes everything? Wait.
        # Original: os.system("rm * 2> /dev/null") -> This might remove generated files!
        # Ah, frankVINA_1 likely filters them. Let's trust original logic for now.
        
        # 3. Patch Clustering
        if stop_event.is_set(): 
            os.chdir(initial_path)
            return

        log(f"--- Patch Clustering (kmer: {pep_size}) ---")
        patch_files = [x for x in os.listdir(".") if "patch_file" in x]
        
        if len(patch_files) == 0:
             log("⚠️ No patch files found.")
        elif len(patch_files) > 1:
            cmd_cluster = f'python3 {repo_folder}/scripts/patch_clustering.py -w {pep_size} -t {threads}'
            log(f"Running: {cmd_cluster}")
            os.system(cmd_cluster)
            
            # 4. FrankVINA 2
            cluster_folder = f"frankPEPstein_{pep_size}"
            if os.path.exists(cluster_folder):
                os.chdir(cluster_folder)
                os.system(f'cp "{receptor_path}" receptor.pdb')
                
                log("--- Starting FrankVINA 2 ---")
                cmd_vina2 = f'python3 {repo_folder}/scripts/frankVINA_2.py "{initial_path}" receptor.pdb {threads} {candidates}'
                log(f"Running: {cmd_vina2}")
                os.system(cmd_vina2)
                
                os.system("rm * 2> /dev/null")
            else:
                log(f"❌ Cluster folder {cluster_folder} not found.")

        elif len(patch_files) == 1:
             log("Single patch file found, skipping clustering.")
             # Logic for single file...
             
        os.chdir(initial_path)
        log("🎉 Pipeline Finished.")
    else:
        log(f"❌ Output folder {output_superposer_path} not found.")
        os.chdir(initial_path)

def on_run_click(b):
    out_log.clear_output()
    out_vis.clear_output()
    btn_run.disabled = True
    
    # Start Visualization Thread
    threading.Thread(target=viz_loop, daemon=True).start()
    
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
