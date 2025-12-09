import os
import sys

# Add repo root to path so we can import scripts
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from scripts import notebook_utils

def run_patching(pockets_dir="work_dir/receptor.pdb_out/pockets", selected_pocket="pocket1_atm.pdb"):
    print("--- Cell 04: Preparation & Patching ---")
    
    pocket_abs_path = os.path.abspath(os.path.join(pockets_dir, selected_pocket))
    if not os.path.exists(pocket_abs_path):
        print(f"Error: Pocket {pocket_abs_path} not found.")
        return

    print(f"Selected Pocket: {pocket_abs_path}")
    
    # Calculate Box
    center, size = notebook_utils.get_pocket_box(pocket_abs_path)
    if center:
        print(f"Box Center: {center}")
        print(f"Box Size: {size}")
    else:
        print("Failed to calculate box.")

    # Patch Scripts
    repo_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    scripts_dir = os.path.join(repo_dir, "scripts")
    utilities_dir = os.path.join(repo_dir, "utilities")
    db_dir = os.path.join(repo_dir, "DB")
    
    path_replacements = {
        "/home/jgutierrez/scripts/": f"{scripts_dir}/",
        "/home/jgutierrez/utilities/./vina_1.2.4_linux_x86_64": "vina", 
        "/home/jgutierrez/scripts/reduce_wwPDB_het_dict.txt": f"{db_dir}/reduce_wwPDB_het_dict.txt",
        "/home/jgutierrez/FrankPEPstein/filtered_DB_P5-15_R30_id10": f"{db_dir}/minipockets_surface80_winsize3_size3_curated-db"
    }
    
    print("Patching scripts...")
    count = notebook_utils.patch_scripts(scripts_dir, path_replacements)
    print(f"Patched {count} scripts.")
    
    return center, size

if __name__ == "__main__":
    run_patching()
