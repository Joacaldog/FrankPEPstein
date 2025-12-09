import subprocess
import os
import glob
import sys

def run_pocket_detection(receptor_path="work_dir/receptor.pdb"):
    print("--- Cell 03: Pocket Detection ---")
    
    if not os.path.exists(receptor_path):
        print(f"Error: Receptor {receptor_path} not found. Run cell_02 first.")
        return

    print(f"Running fpocket on {receptor_path}...")
    try:
        subprocess.run(f"fpocket -f {receptor_path}", shell=True, check=True)
    except subprocess.CalledProcessError:
        print("Error running fpocket. Is it installed?")
        return

    output_folder = f"{receptor_path}_out"
    pockets_dir = os.path.join(output_folder, "pockets")
    
    if os.path.exists(pockets_dir):
        pockets = [f for f in os.listdir(pockets_dir) if f.endswith(".pdb")]
        print(f"Found {len(pockets)} pockets:")
        for p in sorted(pockets):
            print(f" - {p}")
        return pockets_dir
    else:
        print("No pockets directory generated.")
        return None

if __name__ == "__main__":
    if len(sys.argv) > 1:
        run_pocket_detection(sys.argv[1])
    else:
        run_pocket_detection()
