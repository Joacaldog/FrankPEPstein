import os
import sys
import shutil

# Add repo root
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

def run_execution(center=[0,0,0], size=[20,20,20]):
    print("--- Cell 05: Pipeline Execution ---")
    
    repo_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    scripts_dir = os.path.join(repo_dir, "scripts")
    db_dir = os.path.join(repo_dir, "DB")
    
    run_dir = "run_frankpepstein"
    os.makedirs(run_dir, exist_ok=True)
    os.chdir(run_dir)
    print(f"Working in: {os.getcwd()}")
    
    # Dry Run check for scripts
    superposer = os.path.join(scripts_dir, "superposerV5.2_leave1out.py")
    if os.path.exists(superposer):
        print(f"Found Superposer script: {superposer}")
        # Construct command (Printing mainly, as we don't have inputs here without previous steps)
        pepbdb_path = f"{db_dir}/minipockets_surface80_winsize3_size3_curated-db/"
        minipockets_folder = f"{db_dir}/minipockets_surface80_winsize3_size3_curated-db/minipockets/"
        
        cmd = f"python {superposer} -T pocket.pdb -d {pepbdb_path} -r 0.1 -t 2 -a 3 -fm {minipockets_folder} -x_center {center[0]} -y_center {center[1]} -z_center {center[2]} -x_size {size[0]} -y_size {size[1]} -z_size {size[2]}"
        print("Superposer Command (Ready to run):")
        print(cmd)
    else:
        print("Superposer script NOT found.")

if __name__ == "__main__":
    run_execution()
