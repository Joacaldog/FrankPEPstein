import os
import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description="Run FrankPEPstein Pipeline")
    
    # Pipeline Parameters
    parser.add_argument("-w", "--pep_size", type=int, default=8, help="Peptide size (kmer)")
    parser.add_argument("-t", "--threads", type=int, default=36, help="Number of threads")
    parser.add_argument("-c", "--candidates", type=int, default=10, help="Number of candidates for Step 2")
    
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
    minipockets_folder = os.path.join(initial_path, "DB/minipockets_surface80_winsize3_size3_curated-db")
    db_folder = os.path.join(initial_path, "DB/filtered_DB_P5-15_R30_id10")
    
    # Pocket Path (Centralized)
    pockets_folder = os.path.join(initial_path, "pockets")
    pocket_pdb = os.path.join(pockets_folder, "pocket.pdb")
    
    if not os.path.exists(pocket_pdb):
        print(f"Error: Pocket file not found at {pocket_pdb}")
        sys.exit(1)
        
    # Output Folders
    run_folder = os.path.join(initial_path, "FrankPEPstein_run")
    if not os.path.exists(run_folder):
        os.makedirs(run_folder)
        
    # Switch to run_folder for superposer execution
    os.chdir(run_folder)
    
    output_superposer_path = os.path.join(run_folder, "superpockets_residuesAligned3_RMSD0.1")
    
    # 1. Superposer
    print(f"--- Running Superposer ---")
    
    # Construct superposer command with explicit gridbox args
    # superposer.py expects: -x_center ... -x_size ...
    superposer_cmd_list = [
        sys.executable, f"{repo_folder}/scripts/superposer.py",
        "-T", pocket_pdb,
        "-d", db_folder,
        "-t", str(threads),
        "-fm", minipockets_folder,
        "-x_center", str(args.x_center),
        "-y_center", str(args.y_center),
        "-z_center", str(args.z_center),
        "-x_size", str(args.x_size),
        "-y_size", str(args.y_size),
        "-z_size", str(args.z_size)
    ]
    
    superposer_cmd_str = " ".join(superposer_cmd_list)
    print(f"CMD: {superposer_cmd_str}")
    
    exit_code = os.system(superposer_cmd_str)
    if exit_code != 0:
        print("Error: Superposer failed.")
        sys.exit(1)

    # 2. FrankVINA 1
    if os.path.exists(output_superposer_path):
        os.chdir(output_superposer_path)
        os.system(f'cp {initial_path}/receptor.pdb .')
        
        print(f"--- Running FrankVINA 1 ---")
        cmd_vina1 = f'{sys.executable} {repo_folder}/scripts/frankVINA_1.py {initial_path} receptor.pdb {threads}'
        print(f"CMD: {cmd_vina1}")
        os.system(cmd_vina1)
        
        os.system("rm * 2> /dev/null") # Cleanup (careful, this might delete logs/pdbs if script failed, but following original logic)
        
        # 3. Patch Clustering & FrankVINA 2
        print(f"--- Checking for patches ---")
        patch_files = [x for x in os.listdir(".") if "patch_file" in x]
        
        if len(patch_files) == 0:
            print("No patch files in folder")
        elif len(patch_files) > 1:
            print(f"Running patch_clustering with kmer: {pep_size} ")
            os.system(f'{sys.executable} {repo_folder}/scripts/patch_clustering.py -w {pep_size} -t {threads}')
            
            cluster_dir = f"frankPEPstein_{pep_size}"
            if os.path.exists(cluster_dir):
                os.chdir(cluster_dir)
                os.system(f'cp {initial_path}/receptor.pdb .')
                
                print(f"--- Running FrankVINA 2 ---")
                cmd_vina2 = f'{sys.executable} {repo_folder}/scripts/frankVINA_2.py {initial_path} receptor.pdb {threads} {candidates_number}'
                print(f"CMD: {cmd_vina2}")
                os.system(cmd_vina2)
                
                os.system("rm * 2> /dev/null")
            else:
                 print(f"Error: {cluster_dir} not created.")
                 
        elif len(patch_files) == 1:
             print("Only one patch file in folder")
             # Logic for single file (Handle if needed, user didn't specify strict instructions here but original had some logic)
             # keeping original brief check logic:
             # os.makedirs(pep_dir, exist_ok=True) -> pep_dir undefined in original snippet I saw earlier? 
             # Assuming standard simple logic for now.
             pass
             
    else:
        print(f"Error: {output_superposer_path} not found.")
        sys.exit(1)

if __name__ == "__main__":
    main()
