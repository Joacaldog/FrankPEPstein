import os
import sys
import shutil
import subprocess
import json

# --- CONFIGURATION ---
BASE_DIR = os.getcwd() # Run from /home/joacaldo/Onedrive/FrankPEPstein
FRANK_PYTHON = "/home/joacaldo/miniforge3/envs/FrankPEPstein/bin/python"
CLUSTERING_PYTHON = "/home/joacaldo/miniforge3/envs/FrankPEPstein/bin/python" # Assuming same for now
SCRIPTS_DIR = os.path.join(BASE_DIR, "scripts")
DB_PATH = os.path.join(BASE_DIR, "DB", "minipockets_surface80_winsize3_size3_curated-db")

# Test Data
RECEPTOR_FILE = "2AYG.pdb"
POCKET_FILE = "2AYG_out/pockets/pocket9_atm.pdb"
BOX_CENTER = [-3.833, 29.811, -9.029]
BOX_SIZE = [15.231, 9.840, 12.541] # Added buffer to size

# Parameters
PEP_LENGTH = 8
N_PEPS = 10
THREADS = 8

def run_pipeline():
    print(f"Base Directory: {BASE_DIR}")
    
    # 2. Setup Run Directory
    run_dir = os.path.join(BASE_DIR, "FrankPEPstein_Run_Local")
    if os.path.exists(run_dir):
        shutil.rmtree(run_dir)
    os.makedirs(run_dir)
    print(f"Run Directory: {run_dir}")

    # Copy Receptor
    target_receptor = "receptor.pdb"
    shutil.copy(RECEPTOR_FILE, os.path.join(run_dir, target_receptor))
    
    # Prepare Pocket (Chain p)
    target_pocket_file = "target_pocket.pdb"
    print(f"Preparing pocket from: {POCKET_FILE}")
    
    # Isolate BioPython usage
    prep_script_content = f"""
import sys
from Bio import PDB
import os

def prepare_pocket(receptor_path, pocket_atm_path, output_path):
    try:
        parser = PDB.PDBParser(QUIET=True)
        
        # 1. Parse Input Structures
        print(f"Reading receptor from {{receptor_path}}...")
        receptor_struct = parser.get_structure("receptor", receptor_path)
        
        print(f"Reading pocket atoms from {{pocket_atm_path}}...")
        pocket_struct = parser.get_structure("pocket_atm", pocket_atm_path)
        
        # 2. Collect Pocket Atoms
        pocket_atoms = []
        for model in pocket_struct:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        pocket_atoms.append(atom)
        
        print(f"Pocket defined by {{len(pocket_atoms)}} atoms.")

        # 3. Neighbor Search (10.0 Angstrom expansion)
        print("Running NeighborSearch (10.0 A)...")
        # Collect all receptor atoms for search
        receptor_atoms = list(receptor_struct.get_atoms())
        ns = PDB.NeighborSearch(receptor_atoms)
        
        selected_residues = set()
        
        for p_atom in pocket_atoms:
            # search returns atoms/residues/chains/models/level
            # 'R' for residues
            nearby_residues = ns.search(p_atom.get_coord(), 10.0, level='R')
            for res in nearby_residues:
                # Store unique identifier (ChainID, ResID)
                # Ensure we capture chain ID correctly
                selected_residues.add((res.parent.id, res.id))
                
        print(f"Selected {{len(selected_residues)}} residues for new pocket.")

        # 4. Extract and Save
        class PocketSelect(PDB.Select):
            def accept_residue(self, residue):
                return (residue.parent.id, residue.id) in selected_residues

        # Save temporarily
        temp_out = output_path + ".temp"
        io = PDB.PDBIO()
        io.set_structure(receptor_struct)
        io.save(temp_out, PocketSelect())
        
        # 5. Reload and Rename Chain to 'p'
        temp_struct = parser.get_structure("temp", temp_out)
        for model in temp_struct:
            for chain in model:
                chain.id = 'p'
        
        io.set_structure(temp_struct)
        io.save(output_path)
        os.remove(temp_out)
        
        print(f"SUCCESS: Saved {{output_path}}")

    except Exception as e:
        print(f"ERROR: {{e}}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    # Ensure paths are absolute for safety
    # We use f-string interpolation to inject the absolute path calculated by the host script
    rec = '{os.path.join(run_dir, target_receptor)}'
    poc = '{os.path.abspath(POCKET_FILE)}'
    out = '{os.path.join(run_dir, "pocket.pdb")}'
    
    prepare_pocket(rec, poc, out)
"""
    prep_script_path = os.path.join(run_dir, "prep_pocket.py")
    with open(prep_script_path, "w") as f:
        f.write(prep_script_content)
        
    res = subprocess.run([FRANK_PYTHON, prep_script_path], capture_output=True, text=True)
    print(res.stdout)
    if "SUCCESS" not in res.stdout:
        print(f"❌ Error preparing pocket: {res.stderr}")
        return

    print("✅ Pocket Prepared as 'pocket.pdb'.")
    print("🛑 HALTING for Manual Verification.")
    print(f"Please check: {os.path.join(run_dir, 'pocket.pdb')}")
    return # STOP HERE

    # Debug DB
    if not os.path.exists(DB_PATH):

        print(f"❌ DB not found: {DB_PATH}")
        return
    print(f"DB exists. Checking files...")
    # try:
    #     print(len(os.listdir(DB_PATH)))
    # except:
    #     pass

    # --- SUPERPOSER ---
    print("\n🚀 Starting Step 1: Superposer...")
    superposer_script = os.path.join(SCRIPTS_DIR, "superposerV5.2_leave1out.py")
    
    cmd_superposer = [
        FRANK_PYTHON, superposer_script,
        "-T", target_pocket_file, 
        "-d", DB_PATH,
        "-a", "3", 
        "-r", "0.1",
        "-x_center", str(BOX_CENTER[0]),
        "-y_center", str(BOX_CENTER[1]),
        "-z_center", str(BOX_CENTER[2]),
        "-x_size", str(BOX_SIZE[0]),
        "-y_size", str(BOX_SIZE[1]),
        "-z_size", str(BOX_SIZE[2]),
        "-t", str(THREADS),
        "-fm", DB_PATH
    ]

    # Environment
    click_dir = os.path.join(BASE_DIR, "utilities", "Click")
    # ADFR Suite installed in base dir now
    adfr_bin = os.path.join(BASE_DIR, "ADFRsuite_x86_64Linux_1.0", "bin")
    
    conda_bin = os.path.dirname(FRANK_PYTHON)
    env = os.environ.copy()
    env["PATH"] = f"{click_dir}:{adfr_bin}:{conda_bin}:{env.get('PATH', '')}"
    print(f"PATH modified to include: {click_dir} and {adfr_bin}")
    
    # Run
    process = subprocess.Popen(
        cmd_superposer,
        cwd=run_dir,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        env=env,
        bufsize=1
    )
    
    for line in iter(process.stdout.readline, ''):
        print(line, end='')
        
    process.stdout.close()
    ret = process.wait()
    
    if ret != 0:
        print(f"❌ Failed with {ret}")
    else:
        print("✅ Success")

if __name__ == "__main__":
    run_pipeline()
