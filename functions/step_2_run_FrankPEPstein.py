#@title 2. Structure-Guided Peptide Generation
#@markdown **Instructions:**
#@markdown 1. Adjust the Gridbox coordinates if needed.
#@markdown 2. Click **Update View** to see the new box.
#@markdown 3. Click **Run Generation** to start.

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

# --- Parameters (Defaults) ---
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

box_center_init = pipeline_state.get("box_center", [0.0, 0.0, 0.0])
box_size_init = pipeline_state.get("box_size", [20.0, 20.0, 20.0])

if not box_center_init: box_center_init = [0.0, 0.0, 0.0]
if not box_size_init: box_size_init = [20.0, 20.0, 20.0]

# --- Widgets ---
style = {'description_width': 'initial'}
layout_half = widgets.Layout(width='48%')

# Box Center Widgets
w_xc = widgets.FloatText(value=box_center_init[0], description='Center X:', style=style, layout=layout_half)
w_yc = widgets.FloatText(value=box_center_init[1], description='Center Y:', style=style, layout=layout_half)
w_zc = widgets.FloatText(value=box_center_init[2], description='Center Z:', style=style, layout=layout_half)

# Box Size Widgets
w_xs = widgets.FloatText(value=box_size_init[0], description='Size X:', style=style, layout=layout_half)
w_ys = widgets.FloatText(value=box_size_init[1], description='Size Y:', style=style, layout=layout_half)
w_zs = widgets.FloatText(value=box_size_init[2], description='Size Z:', style=style, layout=layout_half)

btn_update = widgets.Button(description='Update View', button_style='info', icon='refresh')
btn_run = widgets.Button(description='Run Generation', button_style='success', icon='play')
btn_stop = widgets.Button(description='Stop & Reset', button_style='danger', icon='stop')

# Use HTML widget for robust threaded updates
out_vis = widgets.HTML(layout={'border': '1px solid #ddd', 'height': '600px', 'width': '100%'})
out_log = widgets.Output(layout={'border': '1px solid #ccc', 'height': '300px', 'overflow_y': 'scroll'})

# Grouping
box_ui = widgets.VBox([
    widgets.Label("Gridbox Parameters:"),
    widgets.HBox([w_xc, w_xs]),
    widgets.HBox([w_yc, w_ys]),
    widgets.HBox([w_zc, w_zs]),
    btn_update
])

main_ui = widgets.VBox([
    box_ui,
    widgets.HBox([btn_run, btn_stop]),
    out_vis,
    out_log
])

# --- Logic ---

def draw_view(extra_pdbs=None):
    try:
        view = py3Dmol.view(width=800, height=600)
        
        # 1. Receptor
        if receptor_path and os.path.exists(receptor_path):
            with open(receptor_path, 'r') as f:
                view.addModel(f.read(), "pdb")
            view.setStyle({'model': -1}, {})
            view.addSurface(py3Dmol.SES, {'opacity': 0.9, 'color': 'white'})
            
        # 2. Pocket
        if extracted_pocket_path and os.path.exists(extracted_pocket_path):
            with open(extracted_pocket_path, 'r') as f:
                view.addModel(f.read(), "pdb")
            view.setStyle({'model': -1}, {'sphere': {'color': 'orange', 'opacity': 0.6}})

        # 3. Gridbox (From Widgets)
        try:
            cx, cy, cz = float(w_xc.value), float(w_yc.value), float(w_zc.value)
            sx, sy, sz = float(w_xs.value), float(w_ys.value), float(w_zs.value)
            
            view.addBox({
                'center': {'x': cx, 'y': cy, 'z': cz},
                'dimensions': {'w': sx, 'h': sy, 'd': sz},
                'color': 'red',
                'opacity': 0.5
            })
        except ValueError:
            pass # Handle transient empty widget values

        # 4. Extra Fragments (Live Updates)
        if extra_pdbs:
            for pdb_file in extra_pdbs:
                 if os.path.exists(pdb_file):
                     with open(pdb_file, 'r') as f:
                         view.addModel(f.read(), "pdb")
                     view.setStyle({'model': -1}, {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.15}})

        view.zoomTo()
        
        # Set HTML content directly - Thread Safe for Widget property
        out_vis.value = view._make_html()
        
    except Exception as e:
        # Avoid printing to stdout during viz loop to prevent interfering with tqdm
        pass

def on_update_click(b):
    draw_view()

# --- Visualization Monitor Thread ---
stop_event = threading.Event()

def monitor_fragments():
    run_folder_name = "FrankPEPstein_run"
    fragments_dir = os.path.join(initial_path, run_folder_name, "superpockets_residuesAligned3_RMSD0.1")
    
    last_count = 0
    
    while not stop_event.is_set():
        if os.path.exists(fragments_dir):
            # Look for patch files
            files = glob.glob(os.path.join(fragments_dir, "patch_file_*.pdb"))
            current_count = len(files)
            
            if current_count > last_count:
                # Update View with new files (Limit to last 50 to avoid lag)
                # Sort by modification time to show newest
                files.sort(key=os.path.getmtime, reverse=True)
                draw_view(files[:50])
                last_count = current_count
                
        time.sleep(2)


def on_run_click(b):
    out_log.clear_output()
    btn_run.disabled = True
    btn_update.disabled = True
    stop_event.clear()
    
    # Start Monitor Thread
    t = threading.Thread(target=monitor_fragments, daemon=True)
    t.start()
    
    with out_log:
        # Input Validation
        if not receptor_path or not extracted_pocket_path:
            print("❌ Error: Receptor or Pocket not found. Please run Step 1 successfully.")
            btn_run.disabled = False
            return

        box_center = [w_xc.value, w_yc.value, w_zc.value]
        box_size = [w_xs.value, w_ys.value, w_zs.value]

        print(f"--- Starting FrankPEPstein Generation ---")
        print(f"Peptide Size: {peptide_size}")
        print(f"Threads: {threads}")
        print(f"Candidates: {candidates}")
        print(f"Gridbox Center: {box_center}")
        print(f"Gridbox Size: {box_size}")
        
        script_path = os.path.join(repo_folder, "scripts/run_FrankPEPstein.py")
        
        cmd_list = [
            sys.executable, "-u", script_path,
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
        try:
            process = subprocess.Popen(
                cmd_list, 
                stdout=subprocess.PIPE, 
                stderr=subprocess.STDOUT, 
                text=True, 
                bufsize=1, 
                universal_newlines=True
            )
            
            # Streaming loop
            for line in iter(process.stdout.readline, ''):
                print(line, end='') 
                
            process.wait()
            
            # Stop monitoring
            stop_event.set()
            
            if process.returncode == 0:
                print("\n✅ Pipeline Finished Successfully.")
                # Final Viz Check
                monitor_fragments() # One last update
            else:
                print(f"\n❌ Pipeline failed with exit code {process.returncode}")
            
        except KeyboardInterrupt:
            stop_event.set()
            print("\n🛑 Pipeline interrupted.")
        except Exception as e:
            stop_event.set()
            print(f"\n❌ Execution Error: {e}")
        finally:
            stop_event.set()
            btn_run.disabled = False
            btn_update.disabled = False

def on_stop_click(b):
    stop_event.set()
    global process
    if 'process' in globals() and process:
         process.terminate()
         
    with out_log:
        print("\n🛑 Stopped by user. Cleaning up...")
        cleanup()
    
    btn_run.disabled = False
    btn_update.disabled = False

def cleanup():
    run_folder_name = "FrankPEPstein_run"
    output_superposer_path = os.path.join(initial_path, run_folder_name, f"superpockets_residuesAligned3_RMSD0.1")
    temp_folder_path = os.path.join(initial_path, run_folder_name, f"temp_folder_residuesAligned3_RMSD0.1")
    
    if os.path.exists(output_superposer_path):
        subprocess.run(f"rm -rf {output_superposer_path}", shell=True)
        # print(f"Removed {output_superposer_path}")
        
    if os.path.exists(temp_folder_path):
        subprocess.run(f"rm -rf {temp_folder_path}", shell=True)
        # print(f"Removed {temp_folder_path}")

btn_update.on_click(on_update_click)
btn_run.on_click(on_run_click)
btn_stop.on_click(on_stop_click)

# Initial Draw
display(main_ui)
draw_view()
