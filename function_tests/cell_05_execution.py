import os
import sys
import shutil

# Add repo root
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import argparse

def run_execution(center_x, center_y, center_z, size_x, size_y, size_z):
    print("--- Cell 05: Pipeline Execution ---")
    
    repo_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    scripts_dir = os.path.join(repo_dir, "scripts")
    db_dir = os.path.join(repo_dir, "DB")
    
    run_dir = "run_frankpepstein"
    os.makedirs(run_dir, exist_ok=True)
    os.chdir(run_dir)
    print(f"Working in: {os.getcwd()}")
    
    # Check inputs (should be copied by previous steps or manual copy for test)
    if not os.path.exists("pocket.pdb"):
         # Try to find from work_dir if running test sequence?
         # Or just assume user copied it?
         # For test automation, let's look for it in ../work_dir/...
         # Hardcoding for 1tig test flow:
         source_pocket = "../work_dir/1tig_out/pockets/pocket1_atm.pdb"
         if os.path.exists(source_pocket):
             shutil.copy(source_pocket, "pocket.pdb")
             shutil.copy("../work_dir/1tig.pdb", "receptor.pdb")
         else:
             print("Warning: pocket.pdb not found in run dir and could not auto-locate.")

    superposer = os.path.join(scripts_dir, "superposerV5.2_leave1out.py")
    if os.path.exists(superposer):
        print(f"Found Superposer script: {superposer}")
        # DB Structure is flat, minipockets are directly in the top folder
        pepbdb_path = f"{db_dir}/minipockets_surface80_winsize3_size3_curated-db/"
        minipockets_folder = f"{db_dir}/minipockets_surface80_winsize3_size3_curated-db/"
        
        # Construct CMD
        # Note: center and size components passed individually
        cmd = [
            "python", superposer,
            "-T", "pocket.pdb",
            "-d", pepbdb_path,
            "-r", "0.1",
            "-t", "4",  # Increase threads slightly for speed if possible
            "-a", "3",
            "-fm", minipockets_folder,
            "-x_center", str(center_x), "-y_center", str(center_y), "-z_center", str(center_z),
            "-x_size", str(size_x), "-y_size", str(size_y), "-z_size", str(size_z)
        ]
        
        print("Executing Superposer...")
        print(" ".join(cmd))
        try:
             import subprocess
             subprocess.run(cmd, check=True)
             print("Superposer execution completed successfully.")
        except subprocess.CalledProcessError as e:
             print(f"Superposer failed with error: {e}")

    else:
        print("Superposer script NOT found.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cx", type=float, default=0)
    parser.add_argument("--cy", type=float, default=0)
    parser.add_argument("--cz", type=float, default=0)
    parser.add_argument("--sx", type=float, default=20)
    parser.add_argument("--sy", type=float, default=20)
    parser.add_argument("--sz", type=float, default=20)
    args = parser.parse_args()
    
    run_execution(args.cx, args.cy, args.cz, args.sx, args.sy, args.sz)
