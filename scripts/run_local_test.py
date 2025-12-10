
import os
import sys
import argparse
import subprocess
import shutil

def setup_test_environment(repo_root, test_dir):
    """Sets up the isolated test environment with symlinks and data copying."""
    print(f"🔧 Setting up test environment in {test_dir}...")
    
    if os.path.exists(test_dir):
        print(f"Cleaning existing test directory {test_dir}...")
        # Be careful not to delete repo_root if user points there!
        if os.path.abspath(test_dir) == os.path.abspath(repo_root):
            print("❌ Cannot use repo root as test directory!")
            sys.exit(1)
        shutil.rmtree(test_dir)
    
    os.makedirs(test_dir)
    
    # Symlink heavy directories
    for folder in ["utilities", "DB", "scripts"]:
        src = os.path.join(repo_root, folder)
        dst = os.path.join(test_dir, folder)
        if os.path.exists(src):
            os.symlink(src, dst)
            print(f"🔗 Linked {folder}")
        else:
            print(f"⚠️ Warning: Source {folder} not found in repo root.")
            
    # Setup Data
    complex_src_root = os.path.join(repo_root, "complex_test")
    test_complex_name = "7ugb_A_I"
    src_target = os.path.join(complex_src_root, test_complex_name)
    
    complex_dst_root = os.path.join(test_dir, "complex_test")
    os.makedirs(complex_dst_root, exist_ok=True)
    
    dst_target = os.path.join(complex_dst_root, test_complex_name)
    
    if os.path.exists(src_target):
        if os.path.exists(dst_target):
             shutil.rmtree(dst_target)
        shutil.copytree(src_target, dst_target)
        print(f"📂 Copied {test_complex_name} to test dir")
    else:
        print(f"❌ {test_complex_name} not found in {complex_src_root}")
        
    return complex_dst_root

def run_pipeline(test_dir, repo_root):
    """Executes the pipeline scripts within the test directory context."""
    os.chdir(test_dir)
    print(f"🚀 Running Pipeline in {os.getcwd()}")
    
    # We must ensure we use the scripts from the test_dir (symlinked) OR absolute paths.
    # We'll use absolute paths to the *symlinked* scripts to ensure relative imports inside them work?
    # No, python resolves symlinks. But let's assume `run_superposer` uses `config.py` which uses REPO_ROOT.
    # `config.REPO_ROOT` is hardcoded absolute path to the REAL repo. 
    # That is GOOD. We want to use the installed tools/DBs.
    # But `run_superposer` writes outputs to CWD (or target dir).
    
    scripts_dir = os.path.join(test_dir, "scripts")
    complex_dir = os.path.join(test_dir, "complex_test")
    run_superposer = os.path.join(scripts_dir, "run_superposer_cluster_biolipcurated.py")
    run_patchclust = os.path.join(scripts_dir, "run_patchclust_cluster_biolipcurated.py")
    
    # Step 0: Prep Pocket (Ensure Chain 'p')
    print("\n▶️  Stage 0: Pocket Preparation (Chain renaming)")
    # Iterate over complexes to prep specific pockets?
    # run_superposer iterates over folders.
    # We should probably do this prep FOR EACH complex in `run_superposer`?
    # Or `run_local_test` preps the test data?
    # `run_local_test` sets up `test_dir/complex_test/7ugb_A_I`.
    # I will prep the pocket inside that folder.
    
    target_7ugb = os.path.join(complex_dir, "7ugb_A_I")
    pocket_pdb = os.path.join(target_7ugb, "pocket.pdb")
    
    if os.path.exists(pocket_pdb):
        prep_script_content = """
import sys
from Bio import PDB
import os

def prepare_pocket(input_path, output_path):
    try:
        parser = PDB.PDBParser(QUIET=True)
        struct = parser.get_structure("pocket", input_path)
        
        # Rename all chains to 'p'
        for model in struct:
            for chain in model:
                chain.id = 'p'
        
        io = PDB.PDBIO()
        io.set_structure(struct)
        io.save(output_path)
        print("SUCCESS")
    except Exception as e:
        print(f"ERROR: {e}")

if __name__ == "__main__":
    prepare_pocket(sys.argv[1], sys.argv[2])
"""
        prep_script_path = os.path.join(test_dir, "prep_pocket.py")
        with open(prep_script_path, "w") as f:
            f.write(prep_script_content)
            
        print(f"Running pocket prep on {pocket_pdb}...")
        subprocess.run([sys.executable, prep_script_path, pocket_pdb, pocket_pdb], check=True)
        print("✅ Pocket chain set to 'p'")
        
    
    # Step 1: Superposer + Vina I
    print("\n▶️  Stage 1: Superposer & Initial Docking")
    
    # Define paths to the DB components within the test directory
    minipockets_db_path = os.path.join(test_dir, "DB", "minipockets_surface80_winsize3_size3_curated-db")
    db_path = os.path.join(test_dir, "DB", "biolip_curated_pockets")
    
    # Ensure submodule exists in test DB (if symlinked DB didn't contain it yet)
    os.makedirs(db_path, exist_ok=True)

    # Augment DB with target for self-matching
    db_target = os.path.join(db_path, "7ugb_A_I")
    if not os.path.exists(db_target):
        print(f"⚠️ Augmenting DB subfolder with target {target_7ugb}...")
        os.symlink(target_7ugb, db_target)

    # 1. Run Superposer
    print("\n🚀 Step 1: Running Superposer...")
    run_superposer_script = os.path.join(scripts_dir, "run_superposer_cluster_biolipcurated.py")
    
    # Ensure config can be imported by adding scripts to path
    env = os.environ.copy()
    env["PYTHONPATH"] = scripts_dir + os.pathsep + env.get("PYTHONPATH", "")

    # Arguments: minipockets_db target_dir -d db_path
    cmd = [
        sys.executable, run_superposer_script,
        minipockets_db_path,
        "--target_dir", complex_dir,
        "-d", db_path
    ]
    
    subprocess.run(cmd, env=env, check=True)
    
    # 2. Run Patch Clustering (Original)
    print("\n🚀 Step 2: Running Patch Clustering...")
    run_patchclust_script = os.path.join(scripts_dir, "run_patchclust_cluster_biolipcurated.py")
    
    cmd_patch = [
        sys.executable, run_patchclust_script,
        "--target_dir", complex_dir
    ]
    subprocess.run(cmd_patch, env=env, check=True)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo_root", default=os.getcwd(), help="Path to main repository")
    parser.add_argument("--test_dir", default="/home/joacaldo/FrankPEPstein_test", help="Path to isolated test environment")
    args = parser.parse_args()
    
    repo_root = os.path.abspath(args.repo_root)
    test_dir = os.path.abspath(args.test_dir)
    
    setup_test_environment(repo_root, test_dir)
    run_pipeline(test_dir, repo_root)

if __name__ == "__main__":
    main()
