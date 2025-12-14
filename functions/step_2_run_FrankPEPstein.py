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

# Import Viz Utils
try:
    from functions import viz_utils
except ImportError:
    # Fallback if running directly or path issues
    sys.path.append(os.path.dirname(__file__))
    import viz_utils

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
modeller_key = 'MODELIRANJE'


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

def ensure_modeller_config(key='MODELIRANJE'):
    """Finds and fixes Modeller config.py with the provided key."""
    try:
        import modeller
        modeller_path = os.path.dirname(modeller.__file__)
        config_path = os.path.join(modeller_path, "config.py")
    except ImportError:
        # Fallback search if module not importable in this context
        config_path = None
        possible_paths = [
             "/usr/local/envs/FrankPEPstein/lib/modeller-*/modlib/modeller/config.py",
             f"{sys.prefix}/lib/modeller-*/modlib/modeller/config.py"
        ]
        for pattern in possible_paths:
            found = glob.glob(pattern)
            if found:
                config_path = found[0]
                break
    
    if config_path and os.path.exists(config_path):
        try:
            with open(config_path, 'w') as f:
                f.write(f"license = '{key}'\n")
                f.write("install_dir = r'/usr/local/envs/FrankPEPstein/lib/modeller-10.8'\n") # Basic config
            return True, config_path
        except Exception as e:
            return False, str(e)
    return False, "Not found"

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
standard_pocket_path = os.path.join(initial_path, "pocket.pdb")

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
# Widgets
# REPLACED: HTML with Image for static Matplotlib rendering
out_vis = widgets.Image(
    layout=widgets.Layout(border='1px solid #ddd', height='500px', width='600px')
)

progress_bar = widgets.FloatProgress(
    value=0.0,
    min=0.0,
    max=100.0,
    description='Progress:',
    bar_style='info',
    style={'bar_color': '#4287f5'},
    layout=widgets.Layout(width='100%')
)

status_label = widgets.Label(
    value="Ready to start...",
    layout=widgets.Layout(width='100%')
)

log_output = widgets.Output(
    layout={'border': '1px solid #ccc', 'height': '200px', 'overflow_y': 'scroll'}
)

# Container
ui_container = widgets.VBox([
    widgets.HBox([out_vis], layout=widgets.Layout(justify_content='center')),
    widgets.HBox([progress_bar]),
    status_label,
    log_output
])

# --- Logic ---

def update_static_viz(extra_pdbs=None, title="Pipeline Running..."):
    try:
        if extracted_pocket_path and os.path.exists(extracted_pocket_path):
            img_bytes = viz_utils.render_static_view(
                extracted_pocket_path, 
                extra_pdbs if extra_pdbs else [],
                title=title
            )
            if img_bytes:
                out_vis.value = img_bytes
    except Exception as e:
        # Fail silently visually, log if needed
        pass

# Initial view
update_static_viz(title="Ready")
display(ui_container)

# --- Threading & Execution ---
stop_event = threading.Event()
pipeline_phase = "Initializing"

def monitor_fragments():
    run_folder_name = "FrankPEPstein_run"
    fragments_dir = os.path.join(initial_path, run_folder_name, "superpockets_residuesAligned3_RMSD0.1")
    
    last_count = 0
    
    while not stop_event.is_set():
        if os.path.exists(fragments_dir):
            files = glob.glob(os.path.join(fragments_dir, "patch_file_*.pdb"))
            current_count = len(files)
            
            # We update periodically regardless of new files to ensure phase title is current
            # Sort by modification time to show newest first, limit to more if static (it's fast)
            files.sort(key=os.path.getmtime, reverse=True)
            
            # Pass top 100 fragments
            update_static_viz(files[:100], title=f"{pipeline_phase} (Fragments: {current_count})")
            
            if current_count > last_count:
                last_count = current_count
        
        # Check stop event every 1s, but wait 10s total interval (faster update for static image)
        for _ in range(10):
            if stop_event.is_set(): break
            time.sleep(1)

import re

def run_step_2():
    global pipeline_phase
    # Input Validation
    # We clear the log output for a new run
    log_output.clear_output()
    progress_bar.value = 0
    progress_bar.bar_style = 'info'
    status_label.value = "Initializing..."
    pipeline_phase = "Initializing"
    
    # 0. Fix Modeller License
    with log_output:
        print("Checking Modeller License...")
        success, msg = ensure_modeller_config(modeller_key)
        if success:
             print(f"✅ Modeller license configured in {msg}")
        else:
             print(f"⚠️ Warning: Could not configure Modeller license: {msg}")

    if not receptor_path or not extracted_pocket_path:
        with log_output:
            print("❌ Error: Receptor or Pocket not found. Please run Step 1 successfully.")
        return
    if not box_center or not box_size:
        with log_output:
            print("❌ Error: Pocket Gridbox not defined. Please run Step 1 successfully.")
        return
    
    with log_output:
        print(f"--- Starting FrankPEPstein Generation ---")
        print(f"Peptide Size: {peptide_size}")
        print(f"Threads: {threads}")
        
    # Determine Python Executable
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
    pipeline_phase = "Scanning minipockets" # Default start
    
    # Start Monitor Thread
    t = threading.Thread(target=monitor_fragments, daemon=True)
    t.start()
    
    try:
        process = subprocess.Popen(
            cmd_list, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.STDOUT, 
            text=True, 
            bufsize=1, # Line buffered
            universal_newlines=True
        )
        
        # Regex for tqdm: "  10%|#         | 10/100 [00:01<00:09,  9.15it/s]"
        # We look for a percentage pattern e.g. " 10%" or "100%"
        tqdm_pattern = re.compile(r'(\d+)%\|.*\| (\d+)/(\d+) \[(.*)\]')
        
        while True:
            line = process.stdout.readline()
            if not line and process.poll() is not None:
                break
            if line:
                clean_line = line.strip()
                
                # Check for Phase Markers
                if "--- Running Superposer ---" in clean_line:
                    pipeline_phase = "Scanning minipockets"
                    status_label.value = pipeline_phase
                elif "--- Running FrankVINA 1 ---" in clean_line:
                    pipeline_phase = "Selecting fragment candidates"
                    status_label.value = pipeline_phase
                    progress_bar.bar_style = 'warning' # Change color to indicate change
                elif "--- Checking for patches ---" in clean_line:
                    pipeline_phase = "Clustering fragments and obtaining combinations of peptides"
                    status_label.value = pipeline_phase
                    progress_bar.bar_style = 'info'
                elif "--- Running FrankVINA 2 ---" in clean_line:
                    pipeline_phase = "Refining peptide candidates and selecting"
                    status_label.value = pipeline_phase
                    progress_bar.bar_style = 'success'

                # Check for progress bar
                match = tqdm_pattern.search(line)
                if match:
                    # Update Widget
                    pct = int(match.group(1))
                    current = match.group(2)
                    total = match.group(3)
                    timing = match.group(4)
                    
                    progress_bar.value = pct
                    status_label.value = f"{pipeline_phase}: {pct}% ({current}/{total}) - {timing}"
                else:
                    # Normal Log
                    if clean_line:
                        with log_output:
                            print(clean_line)
            
        process.wait()
        stop_event.set() # Stop monitor
        
        if process.returncode == 0:
            with log_output:
                print("\n✅ Pipeline Finished Successfully.")
            progress_bar.value = 100
            progress_bar.bar_style = 'success'
            status_label.value = "Completed Successfully"
            
            # Final Viz Update
            pipeline_phase = "Completed"
            monitor_fragments() # One last update calls update_static_viz

        else:
            with log_output:
                print(f"\n❌ Pipeline failed with exit code {process.returncode}")
            progress_bar.bar_style = 'danger'
            status_label.value = "Failed"
        
    except KeyboardInterrupt:
        with log_output:
            print("\n🛑 Pipeline interrupted by user.")
        stop_event.set()
        if 'process' in locals():
            process.terminate()
        cleanup()
    except Exception as e:
        stop_event.set()
        with log_output:
            print(f"\n❌ Execution Error: {e}")
            
def cleanup():
    run_folder_name = "FrankPEPstein_run"
    output_superposer_path = os.path.join(initial_path, run_folder_name, f"superpockets_residuesAligned3_RMSD0.1")
    temp_folder_path = os.path.join(initial_path, run_folder_name, f"temp_folder_residuesAligned3_RMSD0.1")
    
    with log_output:
        if os.path.exists(output_superposer_path):
            subprocess.run(f"rm -rf {output_superposer_path}", shell=True)
            print(f"Removed {output_superposer_path}")
            
        if os.path.exists(temp_folder_path):
            subprocess.run(f"rm -rf {temp_folder_path}", shell=True)
            print(f"Removed {temp_folder_path}")
        print("Cleanup complete.")

if __name__ == "__main__":
    run_step_2()
