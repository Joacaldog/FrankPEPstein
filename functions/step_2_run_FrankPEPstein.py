#@title 2. Structure-Guided Peptide Generation
#@markdown **Instructions:**
#@markdown 1. Configure parameters.
#@markdown 2. Run this cell to start generation.
#@markdown 3. To Stop: Click the "Stop" button on the cell implementation (Interrupt Execution). Cleanup will run automatically.

import os
import sys
import json
import subprocess
import multiprocessing

# --- Parameters ---
peptide_size = 8 #@param {type:"slider", min:5, max:15, step:1}
threads = multiprocessing.cpu_count() #@param {type:"integer"}
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
# Enforce standard pocket
pockets_dir = os.path.join(initial_path, "pockets")
standard_pocket_path = os.path.join(pockets_dir, "pocket.pdb")

if os.path.exists(standard_pocket_path):
    extracted_pocket_path = standard_pocket_path
else:
    extracted_pocket_path = pipeline_state.get("extracted_pocket_path", None)

box_center = pipeline_state.get("box_center", None)
box_size = pipeline_state.get("box_size", None)

# --- Execution ---

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
    
    script_path = os.path.join(repo_folder, "scripts/run_FrankPEPstein.py")
    
    cmd_list = [
        sys.executable, script_path,
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
    
    try:
        # Run subprocess allowing stdout to stream to cell output (inherits streams)
        # This preserves tqdm and prints from the script
        subprocess.run(cmd_list, check=True)
        print("\n✅ Pipeline Finished Successfully.")
        
    except KeyboardInterrupt:
        print("\n🛑 Pipeline interrupted by user.")
        cleanup()
    except subprocess.CalledProcessError as e:
        print(f"\n❌ Pipeline failed with exit code {e.returncode}")
        
def cleanup():
    print("Cleaning up temporary files...")
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
