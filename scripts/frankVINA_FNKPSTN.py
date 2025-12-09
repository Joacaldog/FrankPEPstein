import os
import sys
from tqdm import tqdm
import fnmatch
from operator import itemgetter
from modeller import *
from modeller.scripts import complete_pdb
from modeller.optimizers import ConjugateGradients, MolecularDynamics
from modeller.automodel import autosched
from modeller import log
from joblib import Parallel, delayed
import math
import multiprocessing


tqdm._instances.clear()

receptor_file_chain = sys.argv[1]
receptor_file = receptor_file_chain
threads = sys.argv[2]

def scoring_filter():
    # file_list = []
    energy_patch_list = []
    for file in os.listdir("."):
        if fnmatch.fnmatch(file, '*.log'):
            with open(file) as f:
                for line in f.readlines():
                    line = line.strip()
                    if "Estimated" in line:
                        energy_value = line.split(":")[1].split(" ")[1]
                        energy_patch = float(energy_value), file
                        energy_patch_list.append(energy_patch)

    sorted_list = sorted(energy_patch_list,key=itemgetter(0))
    if len(sorted_list) >= 1:
        folder_output3 = f"../top_10_patches"
        if not os.path.exists(folder_output3):
            os.makedirs(folder_output3)
        with open(f"{folder_output3}/top10_patches.tsv", "w") as outfile:
            outfile.write("AffinityBindingPred\tPEP\n")
            for selected in sorted_list[:10]:
                energy_value = selected[0]
                outfile.write(f"{energy_value}\t{selected[1].replace('.log', '')}\n")
                patch_file = selected[1].replace(".log", ".pdb")
                cmd_cp3 = ("cp {}_out.pdbqt {} 2> /dev/null").format(patch_file.replace(".pdb", ""), folder_output3)
                os.system(cmd_cp3)
                # os.system(f"cp {patch_file} {folder_output3} 2> /dev/null")
    #             final_file = f'{patch_file.replace(".pdb", "")}_out.pdbqt'
    #             file_list.append(final_file)
    # return file_list


def main():
    os.system(f"reduce -Quiet -DB \"/mnt/c/Users/Joacaldo/OneDrive - Universidad Católica de Chile/FrankPEPstein/scripts/reduce_wwPDB_het_dict.txt\" {receptor_file} 1> H_{receptor_file} 2> /dev/null")
    os.system(f"sed -i '/END/d' H_{receptor_file}")
    os.system(f'prepare_receptor -r H_{receptor_file} -o {receptor_file}qt 1> /dev/null 2> /dev/null')
    if not os.path.exists("temp_folder"):
        os.makedirs("temp_folder")
    def vina_scorer(file):
        if fnmatch.fnmatch(file, 'patch_file*.pdb'):
            os.system(f"reduce -Quiet -DB \"/mnt/c/Users/Joacaldo/OneDrive - Universidad Católica de Chile/FrankPEPstein/scripts/reduce_wwPDB_het_dict.txt\" {file} 1> H_{file} 2> /dev/null")
            # os.system(f'prepare_ligand -A bonds,bonds_hydrogens,hydrogens -g -l H_{file} -o {file.replace(".pdb", ".pdbqt")} 1> /dev/null 2> /dev/null')
            os.system(f'prepare_ligand -l H_{file} -o {file.replace(".pdb", ".pdbqt")} 1> /dev/null 2> /dev/null')
            os.system(f'vina --verbosity 0 --autobox --local_only --receptor {receptor_file}qt --ligand {file}qt > {file.replace(".pdb", "")}.log')
            os.system(f'mv {file.replace(".pdb", "_out.pdbqt")} {file.replace(".pdb", "")}.log temp_folder')
            os.system(f'rm H_{file} {file}qt 2> /dev/null')


    Parallel(n_jobs=int(threads))(delayed(vina_scorer)(file) for file in tqdm(os.listdir("."), total=len(os.listdir(".")), 
                                                                                        desc=f"filtering peps by energy"))
    os.chdir("temp_folder")
    scoring_filter()
    os.chdir("../")
    os.system("rm -r temp_folder")
    os.chdir("top_10_patches")
    for pep_pdbqt in os.listdir("."):
        if fnmatch.fnmatch(pep_pdbqt, '*.pdbqt'):
            base = pep_pdbqt.replace("_out.pdbqt", "")
            os.system(f"obabel -ipdbqt {pep_pdbqt} -o pdb -O {base}.pdb")
            os.system(f'rm {pep_pdbqt} 1> /dev/null 2> /dev/null')

if __name__ == '__main__':
    main()
