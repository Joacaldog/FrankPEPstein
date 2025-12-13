#@title 2. Structure-Guided Peptide Generation
#@markdown **Instructions:**
#@markdown 1. Configure parameters.
#@markdown 2. Run this cell to start generation.
#@markdown 3. The 3D view will update every 30 seconds if new fragments are found.

import os
import sys
import json
import time
import glob
import subprocess
import threading
import multiprocessing
import ipywidgets as widgets
from IPython.display import display

# --- Dependency Check ---
try:
    import py3Dmol
except ImportError:
    subprocess.run("pip install -q py3dmol", shell=True)
    import py3Dmol

# --- Parameters ---
peptide_size = 8 #@param {type:"slider", min:5, max:15, step:1}
threads = 0 #@param {type:"integer"}
if threads <= 0:
    threads = multiprocessing.cpu_count()
candidates = 10 #@param {type:"integer"}

# --- Configuration & State ---
initial_path = os.getcwd()
repo_folder = os.path.join(initial_path, "FrankPEPstein")
state_file = "pipeline_state.json"

# Fix permissions
def fix_permissions():
    executables = [
        f"{repo_folder}/utilities/vina_1.2.4_linux_x86_64",
        f"{initial_path}/utilities/click/click"
    ]
    for exe in executables:
        if os.path.exists(exe):
            os.chmod(exe, 0o755)

fix_permissions()

# Load State
pipeline_state = {}
if os.path.exists(state_file):
    try:
        with open(state_file, "r") as f:
            pipeline_state = json.load(f)
    except:
        pass

receptor_path = pipeline_state.get("receptor_filename", None)
pockets_dir = os.path.join(initial_path, "pockets")
standard_pocket_path = os.path.join(pockets_dir, "pocket.pdb")

if os.path.exists(standard_pocket_path):
    extracted_pocket_path = standard_pocket_path
else:
    extracted_pocket_path = pipeline_state.get("extracted_pocket_path", None)

box_center = pipeline_state.get("box_center", [0.0, 0.0, 0.0])
box_size = pipeline_state.get("box_size", [20.0, 20.0, 20.0])

if not box_center or not box_size:
    # Fallback to defaults if missing (shouldn't happen if Step 1 ran)
    box_center = [0.0, 0.0, 0.0]
    box_size = [20.0, 20.0, 20.0]

# --- UI Layout ---
# Using HTML widget for robust threaded updates
out_vis = widgets.HTML(layout={'border': '1px solid #ddd', 'height': '600px', 'width': '100%'})

# --- Logic ---

def generate_view_html(extra_pdbs=None):
    try:
        view = py3Dmol.view(width=800, height=600)
        
        # 1. Receptor (0.85 opacity as requested)
        if receptor_path and os.path.exists(receptor_path):
            with open(receptor_path, 'r') as f:
                view.addModel(f.read(), "pdb")
            view.setStyle({'model': -1}, {})
            view.addSurface(py3Dmol.SES, {'opacity': 0.85, 'color': 'white'})
            
        # 2. Pocket
        if extracted_pocket_path and os.path.exists(extracted_pocket_path):
            with open(extracted_pocket_path, 'r') as f:
                view.addModel(f.read(), "pdb")
            view.setStyle({'model': -1}, {'sphere': {'color': 'orange', 'opacity': 0.6}})

        # 3. Gridbox (Wireframe Edges)
        if box_center and box_size:
            cx, cy, cz = box_center
            sx, sy, sz = box_size
            
            # Calculate corners
            min_x, max_x = cx - sx/2, cx + sx/2
            min_y, max_y = cy - sy/2, cy + sy/2
            min_z, max_z = cz - sz/2, cz + sz/2
            
            p1 = {'x':min_x, 'y':min_y, 'z':min_z}
            p2 = {'x':max_x, 'y':min_y, 'z':min_z}
            p3 = {'x':max_x, 'y':max_y, 'z':min_z}
            p4 = {'x':min_x, 'y':max_y, 'z':min_z}
            
            p5 = {'x':min_x, 'y':min_y, 'z':max_z}
            p6 = {'x':max_x, 'y':min_y, 'z':max_z}
            p7 = {'x':max_x, 'y':max_y, 'z':max_z}
            p8 = {'x':min_x, 'y':max_y, 'z':max_z}
            
            def add_line_edge(start, end):
                view.addLine({'start': start, 'end': end, 'color': 'red', 'linewidth': 5})

            # Bottom Face
            add_line_edge(p1, p2); add_line_edge(p2, p3); add_line_edge(p3, p4); add_line_edge(p4, p1)
            # Top Face
            add_line_edge(p5, p6); add_line_edge(p6, p7); add_line_edge(p7, p8); add_line_edge(p8, p5)
            # Verticals
            add_line_edge(p1, p5); add_line_edge(p2, p6); add_line_edge(p3, p7); add_line_edge(p4, p8)

        # 4. Extra Fragments (Live Updates)
        if extra_pdbs:
            for pdb_file in extra_pdbs:
                 if os.path.exists(pdb_file):
                     with open(pdb_file, 'r') as f:
                         view.addModel(f.read(), "pdb")
                     view.setStyle({'model': -1}, {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.15}})

        view.zoomTo()
        return view._make_html()
        
    except Exception as e:
        return f"<b>Viz Error:</b> {e}"

# Initial view
initial_html = generate_view_html()
out_vis.value = initial_html

# Display UI immediately
display(out_vis)


# --- Threading ---
stop_event = threading.Event()

def monitor_fragments():
    run_folder_name = "FrankPEPstein_run"
    fragments_dir = os.path.join(initial_path, run_folder_name, "superpockets_residuesAligned3_RMSD0.1")
    
    last_count = 0
    
    while not stop_event.is_set():
        if os.path.exists(fragments_dir):
            files = glob.glob(os.path.join(fragments_dir, "patch_file_*.pdb"))
            current_count = len(files)
            
            if current_count > last_count:
                print(f"[Monitor] Updates found! Refreshing 3D View... (Total fragments: {current_count})")
                
                # Sort by modification time to show newest first, limit to 50
                files.sort(key=os.path.getmtime, reverse=True)
                new_html = generate_view_html(files[:50])
                out_vis.value = new_html
                
                last_count = current_count
        
        # Check stop event every 1s, but wait 30s total interval
        for _ in range(30):
            if stop_event.is_set(): break
            time.sleep(1)


def run_step_2():
    # Input Validation
    if not receptor_path or not extracted_pocket_path:
        print("❌ Error: Receptor or Pocket not found. Please run Step 1 successfully.")
        return
    if not box_center or not box_size:
        print("❌ Error: Pocket Gridbox not defined. Please run Step 1 successfully.")
        return

    print(f"--- Starting FrankPEPstein Generation ---")
    print(f"Peptide Size: {peptide_size}")
    print(f"Threads: {threads}")
    print(f"Candidates: {candidates}")
    
    # Determine Python Executable
    # Try to find the Conda Env python first
    conda_python = "/usr/local/envs/FrankPEPstein/bin/python"
    if os.path.exists(conda_python):
        python_exe = conda_python
    else:
        python_exe = sys.executable

    script_path = os.path.join(repo_folder, "scripts/run_FrankPEPstein.py")
    
    cmd_list = [
        python_exe, "-u", script_path,
        "-w", str(peptide_size),
        "-t", str(threads),
        "-c", str(candidates),
        "-xc", str(box_center[0]),
        "-yc", str(box_center[1]),
        "-zc", str(box_center[2]),
        "-xs", str(box_size[0]),
        "-ys", str(box_size[1]),
        "-zs", str(box_size[2])
    ]
    
    global process
    stop_event.clear()
    
    # Start Monitor Thread
    t = threading.Thread(target=monitor_fragments, daemon=True)
    t.start()
    
    try:
        process = subprocess.Popen(
            cmd_list, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.STDOUT, 
            text=True, 
            bufsize=0, # Unbuffered
            universal_newlines=True
        )
        
        # Stream output character by character to handle \r (tqdm) correctly
        while True:
            char = process.stdout.read(1)
            if not char and process.poll() is not None:
                break
            if char:
                print(char, end='')
            
        process.wait()
        
        stop_event.set() # Stop monitor
        
        if process.returncode == 0:
            print("\n✅ Pipeline Finished Successfully.")
            # Final Viz Update
            run_folder_name = "FrankPEPstein_run"
            fragments_dir = os.path.join(initial_path, run_folder_name, "superpockets_residuesAligned3_RMSD0.1")
            if os.path.exists(fragments_dir):
                files = glob.glob(os.path.join(fragments_dir, "patch_file_*.pdb"))
                files.sort(key=os.path.getmtime, reverse=True)
                out_vis.value = generate_view_html(files[:50])
                print(f"Final visualization updated with {len(files)} fragments.")

        else:
            print(f"\n❌ Pipeline failed with exit code {process.returncode}")
        
    except KeyboardInterrupt:
        print("\n🛑 Pipeline interrupted by user.")
        stop_event.set()
        if 'process' in locals():
            process.terminate()
        cleanup()
    except Exception as e:
        stop_event.set()
        print(f"\n❌ Execution Error: {e}")
        
def cleanup():
    run_folder_name = "FrankPEPstein_run"
    output_superposer_path = os.path.join(initial_path, run_folder_name, f"superpockets_residuesAligned3_RMSD0.1")
    temp_folder_path = os.path.join(initial_path, run_folder_name, f"temp_folder_residuesAligned3_RMSD0.1")
    
    if os.path.exists(output_superposer_path):
        subprocess.run(f"rm -rf {output_superposer_path}", shell=True)
        print(f"Removed {output_superposer_path}")
        
    if os.path.exists(temp_folder_path):
        subprocess.run(f"rm -rf {temp_folder_path}", shell=True)
        print(f"Removed {temp_folder_path}")
    print("Cleanup complete.")

if __name__ == "__main__":
    run_step_2()
