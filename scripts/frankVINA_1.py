import os
import sys
import shutil
import subprocess
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
# Configuration Variables

initial_path = sys.argv[1]
receptor_file = "receptor.pdb"
threads = sys.argv[2]
# [ADDED] ADFR Configuration
ADFR_DIR = os.path.join(initial_path, "utilities/ADFRsuite_x86_64Linux_1.0")
ADFR_BIN = os.path.join(ADFR_DIR, "bin")
ADFR_LIB = os.path.join(ADFR_DIR, "lib")

# Prepend to PATH so tools can find their python2.7 and libraries
os.environ["PATH"] = ADFR_BIN + os.pathsep + os.environ["PATH"]
# Set LD_LIBRARY_PATH so binaries find libpython2.7.so.1.0
os.environ["LD_LIBRARY_PATH"] = ADFR_LIB + os.pathsep + os.environ.get("LD_LIBRARY_PATH", "")

# Configuration Variables
# Using absolute paths via ADFR_BIN for clarity
REDUCE_PATH = os.path.join(ADFR_BIN, "reduce")
REDUCE_DB_PATH = f"{initial_path}/DB/reduce_wwPDB_het_dict.txt"
# Using the binaries directly now that they are in PATH context
PREPARE_RECEPTOR_PATH = os.path.join(ADFR_BIN, "prepare_receptor")
PREPARE_LIGAND_PATH = os.path.join(ADFR_BIN, "prepare_ligand")
VINA_PATH = f"{initial_path}/utilities/vina_1.2.4_linux_x86_64"
OBABEL_PATH = os.path.join(ADFR_BIN, "obabel")



tqdm._instances.clear()


def run_tool(command, stdout=None, stderr=None):
    """Run an external tool without a shell; like os.system, a failure is not fatal here."""
    if isinstance(stdout, str):
        with open(stdout, "wb") as handle:
            return subprocess.run(command, stdout=handle, stderr=stderr, check=False)
    return subprocess.run(command, stdout=stdout, stderr=stderr, check=False)


def move_into(folder, *paths):
    """mv paths folder, overwriting, skipping any path the tool did not produce."""
    for path in paths:
        try:
            os.replace(path, os.path.join(folder, os.path.basename(path)))
        except OSError:
            pass


def remove_files(*paths):
    for path in paths:
        try:
            os.remove(path)
        except OSError:
            pass


def strip_end_records(path):
    """sed -i '/END/d'"""
    if os.path.exists(path):
        with open(path, "rb") as handle:
            lines = handle.readlines()
        with open(path, "wb") as handle:
            handle.writelines(line for line in lines if b"END" not in line)



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
                        # Keep only favourable scores, as stage 2 does. Without this the top-10 is
                        # filled with the least-bad positive (i.e. clashing) poses when nothing binds.
                        if float(energy_value) < 0:
                            energy_patch_list.append((float(energy_value), file))

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
                try:
                    shutil.copy(patch_file.replace(".pdb", "_out.pdbqt"), folder_output3)
                except OSError:
                    pass


def main():
    run_tool([REDUCE_PATH, "-Quiet", "-DB", REDUCE_DB_PATH, receptor_file],
             stdout=f"H_{receptor_file}", stderr=subprocess.DEVNULL)
    strip_end_records(f"H_{receptor_file}")
    run_tool([PREPARE_RECEPTOR_PATH, "-r", f"H_{receptor_file}", "-o", f"{receptor_file}qt"])
    if not os.path.exists("temp_folder"):
        os.makedirs("temp_folder")
    def vina_scorer(file):
        if fnmatch.fnmatch(file, 'patch_file*.pdb'):
            base = file.replace(".pdb", "")
            run_tool([REDUCE_PATH, "-Quiet", "-DB", REDUCE_DB_PATH, file],
                     stdout=f"H_{file}", stderr=subprocess.DEVNULL)
            run_tool([PREPARE_LIGAND_PATH, "-l", f"H_{file}", "-o", f"{base}.pdbqt"],
                     stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            run_tool([VINA_PATH, "--verbosity", "0", "--autobox", "--local_only",
                      "--receptor", f"{receptor_file}qt", "--ligand", f"{file}qt"],
                     stdout=f"{base}.log")
            move_into("temp_folder", f"{base}_out.pdbqt", f"{base}.log")
            remove_files(f"H_{file}", f"{file}qt")


    Parallel(n_jobs=int(threads))(delayed(vina_scorer)(file) for file in tqdm(os.listdir("."), total=len(os.listdir(".")), 
                                                                                        desc=f"filtering peps by energy"))
    os.chdir("temp_folder")
    scoring_filter()
    os.chdir("../")
    shutil.rmtree("temp_folder", ignore_errors=True)
    if os.path.exists("top_10_patches"):
        os.chdir("top_10_patches")
        for pep_pdbqt in os.listdir("."):
            if fnmatch.fnmatch(pep_pdbqt, '*.pdbqt'):
                base = pep_pdbqt.replace("_out.pdbqt", "")
                run_tool([OBABEL_PATH, "-ipdbqt", pep_pdbqt, "-o", "pdb", "-O", f"{base}.pdb"])
                remove_files(pep_pdbqt)
    else:
        print("Warning: 'top_10_patches' folder not found. No patches passed energy filter.")

if __name__ == '__main__':
    main()
