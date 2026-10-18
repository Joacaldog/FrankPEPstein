import os
import sys
import argparse
import glob
import subprocess
import shutil

def main():
    parser = argparse.ArgumentParser(description="Run FrankPEPstein Pipeline")
    
    # Pipeline Parameters
    parser.add_argument("-w", "--pep_size", type=int, default=8, help="Peptide size (kmer)")
    parser.add_argument("-t", "--threads", type=int, default=36, help="Number of threads")
    parser.add_argument("-c", "--candidates", type=int, default=10, help="Number of candidates for Step 2")
    parser.add_argument("-s", "--sampling", type=int, default=500, help="Subsampling limit for Vina screening")
    parser.add_argument("-rmsd", "--rmsd_allowed", type=float, default=0.5, help="RMSD cutoff for superposer")
    parser.add_argument("-seed", "--seed", type=int, default=0, help="Random seed for the sampling steps, for reproducible runs")
    
    # Gridbox Parameters (Required for superposer)
    parser.add_argument("-xc", "--x_center", type=float, required=True, help="Resulting box center X")
    parser.add_argument("-yc", "--y_center", type=float, required=True, help="Resulting box center Y")
    parser.add_argument("-zc", "--z_center", type=float, required=True, help="Resulting box center Z")
    parser.add_argument("-xs", "--x_size", type=float, required=True, help="Resulting box size X")
    parser.add_argument("-ys", "--y_size", type=float, required=True, help="Resulting box size Y")
    parser.add_argument("-zs", "--z_size", type=float, required=True, help="Resulting box size Z")
    
    args = parser.parse_args()
    
    pep_size = args.pep_size
    threads = args.threads
    candidates_number = args.candidates
    
    initial_path = os.getcwd()
    # Derive the repo root from this file rather than cwd (see run_local.py for the rationale).
    # initial_path stays cwd: it locates DB/, utilities/ and the run folder, which are data, not code.
    repo_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    # DB Paths
    minipockets_folder = os.path.join(initial_path, "DB/minipockets_surface80_winsize3_size3")
    db_folder = os.path.join(initial_path, "DB/filtered_DB_P5-15_R30_id10")
    
    # Pocket Path (Now in FrankPEPstein_run due to refactor)
    # pockets_folder was removed.
    run_folder_init = os.path.join(initial_path, "FrankPEPstein_run")
    pocket_pdb = os.path.join(run_folder_init, "pocket.pdb")
    
    if not os.path.exists(pocket_pdb):
        print(f"Error: Pocket file not found at {pocket_pdb}")
        sys.exit(1)
        
    # Output Folders
    run_folder = os.path.join(initial_path, "FrankPEPstein_run")
    if not os.path.exists(run_folder):
        os.makedirs(run_folder)
        
    # Switch to run_folder for superposer execution
    os.chdir(run_folder)
    
    output_superposer_path = os.path.join(run_folder, f"superpockets_residuesAligned3_RMSD{args.rmsd_allowed}")
    
    # 1. Superposer
    print(f"--- Running Superposer ---")
    
    # Construct superposer command with explicit gridbox args
    # superposer.py expects: -x_center ... -x_size ...
    superposer_cmd_list = [
        sys.executable, f"{repo_folder}/scripts/superposer.py",
        "-i", initial_path,
        "-T", "pocket.pdb", # Simply filename, we are in the dir
        "-d", db_folder,
        "-t", str(threads),
        "-fm", minipockets_folder,
        "-x_center", str(args.x_center),
        "-y_center", str(args.y_center),
        "-z_center", str(args.z_center),
        "-x_size", str(args.x_size),
        "-y_size", str(args.y_size),
        "-z_size", str(args.z_size),
        "-rmsd", str(args.rmsd_allowed)
    ]
    
    superposer_cmd_str = " ".join(superposer_cmd_list)
    print(f"CMD: {superposer_cmd_str}")
    
    subprocess.run(superposer_cmd_list, check=True)

    # 2. FrankVINA 1
    if os.path.exists(output_superposer_path):
        os.chdir(output_superposer_path)
        shutil.copy(os.path.join(initial_path, "receptor.pdb"), ".") # Replaced os.system('cp ...') with shutil.copy
        
        # DIAGNOSTIC STAGE. FrankVINA 1 docks and scores every patch and writes the ten best to
        # top_10_patches/, but nothing reads that folder: patch_clustering below globs every
        # patch_file_*.pdb from this directory, unfiltered. The ranking is a report for the user,
        # not a filter. Do not wire it in as one without checking that the obabel pdbqt->pdb
        # round-trip preserves residue names, which the rest of the pipeline depends on.
        print(f"--- Running FrankVINA 1 (diagnostic scoring; does not filter the patch set) ---")
        cmd_vina1 = [sys.executable, f'{repo_folder}/scripts/frankVINA_1.py', initial_path, str(threads)]
        print(f"CMD: {' '.join(cmd_vina1)}")
        subprocess.run(cmd_vina1, check=True)
        
        # NOTE: do NOT add a cleanup sweep of this directory here. An earlier `rm *` lived at this
        # point and never actually deleted anything (no shell, so the `*` was a literal filename).
        # The patch discovery immediately below reads patch_file_*.pdb from this very directory, so a
        # working glob would wipe the superposer's entire output and leave the run silently empty.
        # The only leftovers are H_receptor.pdb and receptor.pdbqt, which are overwritten on re-run.

        # 3. Patch Clustering & FrankVINA 2
        print(f"--- Checking for patches ---")
        patch_files = [x for x in os.listdir(".") if "patch_file" in x]
        
        if len(patch_files) == 0:
            print("No patch files in folder")
        elif len(patch_files) > 1:
            print(f"Running patch_clustering with kmer: {pep_size} ")
            cmd_clustering = [sys.executable, f'{repo_folder}/scripts/patch_clustering.py', '-w', str(pep_size), '-t', str(threads), '-c', str(args.sampling), '-seed', str(args.seed)]
            print(f"CMD: {' '.join(cmd_clustering)}")
            # patch_clustering creates its output directory at import time, so os.path.exists()
            # below cannot tell success from failure. Check the exit status instead.
            clustering_rc = subprocess.run(cmd_clustering, check=False).returncode

            cluster_dir = f"frankPEPstein_{pep_size}"
            if clustering_rc != 0:
                print(f"Error: patch_clustering exited with code {clustering_rc}; skipping FrankVINA 2.")
            elif os.path.exists(cluster_dir):
                os.chdir(cluster_dir)
                shutil.copy(pocket_pdb, ".") # Modified: Copy pocket instead of receptor for VINA 2 speedup
                # Also stage the full receptor: frankVINA_2 screens against the pocket shell for
                # speed, then rescores the surviving candidates against the whole protein.
                full_receptor = os.path.join(initial_path, "receptor.pdb")
                if os.path.exists(full_receptor):
                    shutil.copy(full_receptor, ".")
                else:
                    print(f"Warning: {full_receptor} not found; skipping the full-receptor rescore.")

                print(f"--- Running FrankVINA 2 ---")
                # List form, not a shell string: initial_path is an arbitrary filesystem path.
                cmd_vina2 = [sys.executable, f'{repo_folder}/scripts/frankVINA_2.py',
                             initial_path, str(threads), str(candidates_number),
                             str(args.sampling), str(args.seed)]
                print(f"CMD: {' '.join(cmd_vina2)}")
                vina2_rc = subprocess.run(cmd_vina2, check=False).returncode

                if vina2_rc != 0:
                    # Do NOT clean up after a failure: the fragments are the only record of what
                    # stage 4 was given, and a bare `rm *` here used to erase them while the run
                    # still reported success.
                    print(f"Error: FrankVINA 2 exited with code {vina2_rc}; leaving {cluster_dir} intact for inspection.")
                else:
                    # Narrow cleanup of the staged inputs only, never a bare glob.
                    for leftover in glob.glob("frag_*.pdb"):
                        os.remove(leftover)
                    for staged in ("pocket.pdb", "receptor.pdb"):
                        if os.path.exists(staged):
                            os.remove(staged)

                # Return to the superposer output dir: the patch-file move below assumes it.
                os.chdir(output_superposer_path)
            else:
                 print(f"Error: {cluster_dir} not created.")

        elif len(patch_files) == 1:
             print("Only one patch file in folder")
             pass
             
        # [ADDED] Move used patch files to 'all_patches_found' folder for organization
        all_patches_dir = "all_patches_found"
        if os.path.exists(output_superposer_path): # Ensure we are in the right context/dir or check files locally
           # We are already in output_superposer_path due to os.chdir above line 83
           if not os.path.exists(all_patches_dir):
               os.makedirs(all_patches_dir)
           
           print("Moving patch files to 'all_patches_found'...")
           for patch_file in glob.glob("patch_file_*.pdb"):
               os.replace(patch_file, os.path.join(all_patches_dir, patch_file))

    else:
        print(f"Error: {output_superposer_path} not found.")
        sys.exit(1)

if __name__ == "__main__":
    main()
