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

# --- Progress Widgets ---
prog_style = {'description_width': 'initial'}
layout_full = widgets.Layout(width='98%')

# 1. Superposition / Fragment Finding
pb_super = widgets.IntProgress(value=0, min=0, max=100, description='1. Superposition:', style=prog_style, layout=layout_full)
lbl_super = widgets.Label(value="Waiting...")

# 2. Fragment Filtering (FrankVINA 1)
pb_filter1 = widgets.IntProgress(value=0, min=0, max=100, description='2. Frag. Filter:', style=prog_style, layout=layout_full)
lbl_filter1 = widgets.Label(value="Waiting...")

# 3. Clustering
pb_cluster = widgets.IntProgress(value=0, min=0, max=100, description='3. Clustering:', style=prog_style, layout=layout_full)
lbl_cluster = widgets.Label(value="Waiting...")

# 4. Peptide Filtering (FrankVINA 2)
pb_filter2 = widgets.IntProgress(value=0, min=0, max=100, description='4. Pep. Filter:', style=prog_style, layout=layout_full)
lbl_filter2 = widgets.Label(value="Waiting...")

ui_progress = widgets.VBox([
    widgets.HBox([pb_super, lbl_super]),
    widgets.HBox([pb_filter1, lbl_filter1]),
    widgets.HBox([pb_cluster, lbl_cluster]),
    widgets.HBox([pb_filter2, lbl_filter2])
])

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

    log("--- Starting Pipeline ---")
    log(f"Command: {' '.join(cmd_list)}")
    
    # reset bars
    pb_super.value = 0; lbl_super.value = "Waiting..."
    pb_filter1.value = 0; lbl_filter1.value = "Waiting..."
    pb_cluster.value = 0; lbl_cluster.value = "Waiting..."
    pb_filter2.value = 0; lbl_filter2.value = "Waiting..."
    
    try:
        # Use Popen with PIPE for stdout to read line by line
        process_handle = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, cwd=initial_path, bufsize=1, universal_newlines=True)
        
        # Monitor Loop
        for line in iter(process_handle.stdout.readline, ''):
            if stop_event.is_set():
                process_handle.terminate()
                log("🛑 Process stopped by user.")
                return
            
            line_str = line.strip()
            if not line_str: continue
            
            # log(line_str) # Optional: too verbose for main log? Maybe just errors or key lines? 
            # Let's log it to out_log so user sees details if they scroll
            # But maybe filter tqdm bars?
            if "Minimizing" not in line_str and "%" not in line_str:
                 log(line_str)
            
            # --- Progress Logic ---
            if "Running Superposer" in line_str:
                pb_super.value = 10; lbl_super.value = "Running..."
                pb_filter1.value = 0; lbl_filter1.value = "Waiting..."
                
            elif "Running FrankVINA 1" in line_str:
                pb_super.value = 100; lbl_super.value = "Done"
                pb_filter1.value = 10; lbl_filter1.value = "Running..."
                
            elif "checking for patches" in line_str.lower():
                pb_filter1.value = 90; lbl_filter1.value = "Checking..."
                
            elif "Running patch_clustering" in line_str:
                pb_filter1.value = 100; lbl_filter1.value = "Done"
                pb_cluster.value = 10; lbl_cluster.value = "Running..."
                
            elif "Running FrankVINA 2" in line_str:
                pb_cluster.value = 100; lbl_cluster.value = "Done"
                pb_filter2.value = 10; lbl_filter2.value = "Running..."
            
            elif "Minimizing complexes" in line_str:
                pb_filter2.value = 50; lbl_filter2.value = "Minimizing..."
                
            elif "Converting top candidates" in line_str:
                pb_filter2.value = 90; lbl_filter2.value = "Finalizing..."
            
            # TQDM parsing (Simplistic)
            # If we see tqdm output we could update the active bar value?
            # It's hard to parse \r lines reliably without logic. 
            # We'll rely on stage headers for big jumps.
            
        process_handle.wait()
        ret = process_handle.returncode
            
        if ret != 0:
            log(f"❌ Pipeline failed with code {ret}")
            return
            
        # Final set
        pb_filter2.value = 100; lbl_filter2.value = "Done"
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
    temp_folder_path = os.path.join(initial_path, run_folder_name, f"temp_folder_residuesAligned3_RMSD0.1")
    
    if os.path.exists(output_superposer_path):
        log(f"Cleaning up {output_superposer_path}...")
        os.system(f"rm -rf {output_superposer_path}")
        
    if os.path.exists(temp_folder_path):
        log(f"Cleaning up {temp_folder_path}...")
        os.system(f"rm -rf {temp_folder_path}")
        
    btn_run.disabled = False
    log("Stopped and Reset.")
btn_run.on_click(on_run_click)
btn_stop.on_click(on_stop_click)

# --- Visualization Logic ---
def viz_loop():
    try:
        # Monitor output folder
        run_folder_name = "FrankPEPstein_run"
        super_out_name = "superpockets_residuesAligned3_RMSD0.1"
        target_dir = os.path.join(initial_path, run_folder_name, super_out_name)
        all_patches_dir = os.path.join(target_dir, "all_patches_found") 
        
        last_count = 0
        showing_final = False
        
        log("Visualization monitoring started...")
        
        # Initial Base View (Immediate)
        with out_vis:
            out_vis.clear_output(wait=True)
            view = py3Dmol.view(width=800, height=600)
            
            if receptor_path and os.path.exists(receptor_path):
                with open(receptor_path, 'r') as f: view.addModel(f.read(), "pdb")
                view.setStyle({'model': -1}, {})
                view.addSurface(py3Dmol.SES, {'opacity': 0.9, 'color': 'white'})
            
            if extracted_pocket_path and os.path.exists(extracted_pocket_path):
                 with open(extracted_pocket_path, 'r') as f: view.addModel(f.read(), "pdb")
                 view.setStyle({'model': -1}, {'sphere': {'color': 'orange', 'opacity': 0.5}})
                 
            if box_center and box_size:
                 cx, cy, cz = box_center; sx, sy, sz = box_size
                 view.addBox({'center': {'x': cx, 'y': cy, 'z': cz},'dimensions': {'w': sx, 'h': sy, 'd': sz},'color': 'red','opacity': 0.5})
                 
            view.zoomTo()
            # USE display() explicitly for thread-safe widget updating
            display(view)
        
        while True:
            if stop_event.is_set(): break
            
            # Final Check if pipeline finished
            if btn_run.disabled == False:
                if not showing_final and os.path.exists(target_dir):
                    final_folders = glob.glob(os.path.join(target_dir, "top_*_peps"))
                    if final_folders:
                         final_folder = final_folders[0]
                         final_peps = glob.glob(os.path.join(final_folder, "*.pdb"))
                         
                         if final_peps:
                             log(f"Pipeline Finished. Found {len(final_peps)} Final Peptides.")
                             
                             with out_vis:
                                out_vis.clear_output(wait=True)
                                view = py3Dmol.view(width=800, height=600)
                                
                                if receptor_path and os.path.exists(receptor_path):
                                    with open(receptor_path, 'r') as f: view.addModel(f.read(), "pdb")
                                    view.setStyle({'model': -1}, {})
                                    view.addSurface(py3Dmol.SES, {'opacity': 0.9, 'color': 'white'})
                                    
                                if extracted_pocket_path and os.path.exists(extracted_pocket_path):
                                     with open(extracted_pocket_path, 'r') as f: view.addModel(f.read(), "pdb")
                                     view.setStyle({'model': -1}, {'sphere': {'color': 'orange', 'opacity': 0.3}})
                                
                                for pep in final_peps:
                                    with open(pep, 'r') as f: view.addModel(f.read(), "pdb")
                                    view.setStyle({'model': -1}, {'stick': {'colorscheme': 'blueCarbon'}})
                                    
                                view.zoomTo()
                                display(view)
                             showing_final = True
                break
            
            # Live Fragment Monitoring
            if os.path.exists(target_dir):
                files = glob.glob(os.path.join(target_dir, "patch_file_*.pdb"))
                count = len(files)
                
                if count >= last_count + 5 and count > 0:
                     with out_vis:
                        out_vis.clear_output(wait=True)
                        view = py3Dmol.view(width=800, height=600)
                        
                        if receptor_path and os.path.exists(receptor_path):
                            with open(receptor_path, 'r') as f: view.addModel(f.read(), "pdb")
                            view.setStyle({'model': -1}, {})
                            view.addSurface(py3Dmol.SES, {'opacity': 0.9, 'color': 'white'})
                            
                        if extracted_pocket_path and os.path.exists(extracted_pocket_path):
                            with open(extracted_pocket_path, 'r') as f: view.addModel(f.read(), "pdb")
                            view.setStyle({'model': -1}, {'sphere': {'color': 'orange', 'opacity': 0.5}})
                            
                        if box_center and box_size:
                             cx, cy, cz = box_center; sx, sy, sz = box_size
                             view.addBox({'center': {'x': cx, 'y': cy, 'z': cz},'dimensions': {'w': sx, 'h': sy, 'd': sz},'color': 'red','opacity': 0.5})
    
                        sorted_files = sorted(files, key=os.path.getmtime, reverse=True)[:50]
                        for pf in sorted_files:
                            with open(pf, 'r') as f: view.addModel(f.read(), "pdb")
                            view.setStyle({'model': -1}, {'stick': {'colorscheme': 'greenCarbon'}})
                            
                        view.zoomTo()
                        display(view)
                     
                     last_count = count
            
            time.sleep(2)
            
    except Exception as e:
        log(f"❌ Visualization Thread Error: {e}")
        import traceback
        traceback.print_exc()

# --- Layout ---
ui = widgets.VBox([
    widgets.HBox([w_pep_size, w_threads, w_candidates]),
    widgets.HBox([btn_run, btn_stop]),
    ui_progress,
    out_log,
    out_vis
])

display(ui)
