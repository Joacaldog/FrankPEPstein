import os
import shutil
import sys

def run_input(receptor_path="receptor.pdb"):
    print("--- Cell 02: Input Data ---")
    
    if not os.path.exists(receptor_path):
        print(f"Error: Receptor file '{receptor_path}' not found.")
        # Create a dummy receptor for testing if not exists? 
        # Better to warn user.
        return None

    print(f"Receptor found: {receptor_path}")
    
    # Move to working directory
    os.makedirs("work_dir", exist_ok=True)
    dest_path = f"work_dir/{os.path.basename(receptor_path)}"
    shutil.copy(receptor_path, dest_path)
    
    print(f"Receptor copied to: {dest_path}")
    return dest_path

if __name__ == "__main__":
    # Allow passing file from command line
    if len(sys.argv) > 1:
        run_input(sys.argv[1])
    else:
        run_input()
