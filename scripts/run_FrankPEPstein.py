import os
import sys
initial_path = os.getcwd()
repo_folder = os.path.join(initial_path, "FrankPEPstein")
minipockets_folder = os.path.join(initial_path, "DB/minipockets_surface80_winsize3_size3_curated-db")
db_folder = os.path.join(initial_path, "DB/filtered_DB_P5-15_R30_id10")
fpocket_run_folder = os.path.join(initial_path, "FPocket_run")
run_folder = os.path.join(initial_path, "FrankPEPstein_run")
output_superposer_path = os.path.join(run_folder, "superpockets_residuesAligned3_RMSD0.1")
superposer_cmd = f'python3 {repo_folder}/scripts/superposer.py -T {fpocket_folder}/pocket.pdb -d {db_folder} -r 0.1 -t {threads} -a 3 -fm {minipockets_folder}'
os.system(superposer_cmd)
os.chdir(output_superposer_path)
os.system(f'cp {initial_path}/receptor.pdb .')
print(f'python3 {repo_folder}/scripts/frankVINA_1.py {initial_path} receptor.pdb {threads}')
os.system(f'python3 {repo_folder}/scripts/frankVINA_1.py {initial_path} receptor.pdb {threads}')
os.system("rm * 2> /dev/null")
pep_size = 8
threads = 36
candidates_number = 10
os.chdir(output_superposer_path)
if len([x for x in os.listdir(".") if "patch_file" in x]) == 0:
    print("No patch files in folder")
elif len([x for x in os.listdir(".") if "patch_file" in x]) > 1:
    print(f"Running patch_clustering with kmer: {pep_size} ")
    os.system(f'python3 {repo_folder}/scripts/patch_clustering.py -w {pep_size} -t {threads}')
    os.chdir(f"frankPEPstein_{pep_size}")
    os.system(f'cp {initial_path}/receptor.pdb .')
    print(f'python3 {repo_folder}/scripts/frankVINA_2.py {initial_path} receptor.pdb {threads} {candidates_number}')
    os.system(f'python3 {repo_folder}/scripts/frankVINA_2.py {initial_path} receptor.pdb {threads} {candidates_number}')
    os.system("rm * 2> /dev/null")
elif len([x for x in os.listdir(".") if "patch_file" in x]) == 1:
    print("Only one patch file in folder")
    os.makedirs(pep_dir, exist_ok=True)
    os.system(f"cp patch_file*.pdb {pep_dir}frag_1.pdb")
elif os.path.exists(pep_tsv):
    print("already processed")
