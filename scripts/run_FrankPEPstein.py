import os
import sys
import argparse
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
    repo_folder = os.path.join(initial_path, "FrankPEPstein")
    
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
        
        print(f"--- Running FrankVINA 1 ---")
        cmd_vina1 = [sys.executable, f'{repo_folder}/scripts/frankVINA_1.py', initial_path, str(threads)]
        print(f"CMD: {' '.join(cmd_vina1)}")
        subprocess.run(cmd_vina1, check=True)
        
        # os.system("rm * 2> /dev/null") # Cleanup (careful, this might delete logs/pdbs if script failed, but following original logic)
        # Replaced with subprocess.run
        subprocess.run(["rm", "*"], check=False, stderr=subprocess.DEVNULL) # check=False as original had 2> /dev/null
        
        # 3. Patch Clustering & FrankVINA 2
        print(f"--- Checking for patches ---")
        patch_files = [x for x in os.listdir(".") if "patch_file" in x]
        
        if len(patch_files) == 0:
            print("No patch files in folder")
        elif len(patch_files) > 1:
            print(f"Running patch_clustering with kmer: {pep_size} ")
            cmd_clustering = [sys.executable, f'{repo_folder}/scripts/patch_clustering.py', '-w', str(pep_size), '-t', str(threads), '-c', '2000']
            print(f"CMD: {' '.join(cmd_clustering)}")
            subprocess.run(cmd_clustering, check=True)
            
            cluster_dir = f"frankPEPstein_{pep_size}"
            if os.path.exists(cluster_dir):
                os.chdir(cluster_dir)
                shutil.copy(pocket_pdb, ".") # Modified: Copy pocket instead of receptor for VINA 2 speedup
                
                print(f"--- Running FrankVINA 2 ---")
                cmd_vina2 = f'{sys.executable} {repo_folder}/scripts/frankVINA_2.py {initial_path} {threads} {candidates_number} {args.sampling}'
                print(f"CMD: {cmd_vina2}")
                os.system(cmd_vina2)
                
                os.system("rm * 2> /dev/null")
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
           os.system(f"mv patch_file_*.pdb {all_patches_dir} 2> /dev/null")

    else:
        print(f"Error: {output_superposer_path} not found.")
        sys.exit(1)

if __name__ == "__main__":
    main()
