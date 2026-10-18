#GRIDBOX ADDED

import os
import sys
import shutil
import subprocess
from tqdm import tqdm
import fnmatch
from Bio.PDB import *
import Bio.PDB
pdbio = Bio.PDB.PDBIO()
parser = Bio.PDB.MMCIFParser()
import glob
import warnings
from joblib import Parallel, delayed
import numpy as np
import argparse

program_description = "Select and generate fragment of peptides that could eventually bind to target receptor based on minipocket alignments"
parser = argparse.ArgumentParser(description=program_description)

# [ADDED] Initial Path Argument
parser.add_argument("-i", "--initial_path", type=str,
                    help="Absolute path to the main directory (containing DB, utilities, etc.)", required=True)

parser.add_argument("-T", "--target_receptor", type=str,
                    help="target receptor to be scanned for posibles fragments of peptides to bind, MUST be on the same folder", required=True)
parser.add_argument("-d", "--pepbdb_folder", type=str,
                    help="absolute path to database folder", required=True)
parser.add_argument("-x_center", "--x_center", type=float,
                    help="x coordinate of center for gridbox", required=True)
parser.add_argument("-y_center", "--y_center", type=float,
                    help="y coordinate of center for gridbox", required=True)
parser.add_argument("-z_center", "--z_center", type=float,
                    help="z coordinate of center for gridbox", required=True)

parser.add_argument("-x_size", "--x_size", type=float,
                    help="x size for gridbox", required=True)
parser.add_argument("-y_size", "--y_size", type=float,
                    help="y size for gridbox", required=True)
parser.add_argument("-z_size", "--z_size", type=float,
                    help="z size for gridbox", required=True)

parser.add_argument("-t", "--threads", type=int,
                    help="Number of threads", required=True)
parser.add_argument("-fm", "--folder_minipockets", type=str,
                    help="folder containing minipockets", required=True)
parser.add_argument("-rmsd", "--rmsd_allowed", type=float, default=0.5,
                    help="RMSD cutoff for superposition")

args = parser.parse_args()

# --- [ADDED] Path Configuration ---
initial_path = args.initial_path
CLICK_DIR = os.path.join(initial_path, "utilities/click")
CLICK_PATH = os.path.join(CLICK_DIR, "click")
if not os.path.exists(CLICK_PATH):
     CLICK_PATH = os.path.join(initial_path, "utilities/click")

PARAMETERS_INP_PATH = os.path.join(CLICK_DIR, "Parameters.inp")

print(f"DEBUG: Using CLICK_PATH={CLICK_PATH}")
print(f"DEBUG: Using PARAMETERS_INP_PATH={PARAMETERS_INP_PATH}")

# Fail loudly here rather than silently producing zero patches. Every CLICK invocation below
# sends stdout and stderr to /dev/null inside a worker with a bare except, so a missing or
# non-executable binary is otherwise indistinguishable from "the database had no good hits".
# Note the fallback above can leave CLICK_PATH pointing at the *directory*, which is executable.
if not (os.path.isfile(CLICK_PATH) and os.access(CLICK_PATH, os.X_OK)):
    sys.exit(f"ERROR: CLICK binary not found or not executable at {CLICK_PATH}. "
             f"Re-run the setup step so utilities/click/click exists.")
if not os.path.isfile(PARAMETERS_INP_PATH):
    sys.exit(f"ERROR: CLICK Parameters.inp not found at {PARAMETERS_INP_PATH}.")
# ----------------------------------

fau_file = args.target_receptor
pepbdb_folder = args.pepbdb_folder
cutoff = 3
threads = args.threads
folder_minipockets = args.folder_minipockets
RMSD_allowed = args.rmsd_allowed
x_center = args.x_center
y_center = args.y_center
z_center = args.z_center
x_size = args.x_size
y_size = args.y_size
z_size = args.z_size

folder_temp = (f"temp_folder_residuesAligned{cutoff}_RMSD{RMSD_allowed}")
# Wipe rather than reuse: this is scratch, and leftovers from an aborted run get appended to
# (see the super_file merge below) and silently corrupt the current run's superpositions.
# Must happen before the two copies into the temp dir below.
if os.path.exists(folder_temp):
    shutil.rmtree(folder_temp, ignore_errors=True)
os.makedirs(folder_temp, exist_ok=True)

folder_output = (f"superpockets_residuesAligned{cutoff}_RMSD{RMSD_allowed}")
if not os.path.exists(folder_output):
    os.makedirs(folder_output)

working_directory = os.getcwd()

# Copy Target to Temp (shutil, not a shell: these paths are user-supplied and may contain spaces)
shutil.copy(fau_file, folder_temp)
# Copy Parameters (using configured path)
shutil.copy(PARAMETERS_INP_PATH, folder_temp)


def remove_files(*paths):
    """rm -f: delete whatever exists and ignore the rest."""
    for path in paths:
        try:
            os.remove(path)
        except OSError:
            pass


def read_lines(path):
    """The file's lines as bytes, or none at all if CLICK did not write it."""
    try:
        with open(path, "rb") as handle:
            return handle.readlines()
    except OSError:
        return []


def strip_end_records(path):
    """sed -i '/END/d': CLICK's outputs end in END, which would cut a concatenation short."""
    if os.path.exists(path):
        lines = read_lines(path)
        with open(path, "wb") as handle:
            handle.writelines(line for line in lines if b"END" not in line)


def concatenate(destination, *sources):
    """cat sources > destination. A missing source contributes nothing, as with cat."""
    with open(destination, "wb") as handle:
        for source in sources:
            handle.writelines(read_lines(source))


def run_click_binary(*click_args):
    """Invoke CLICK without a shell, so paths containing spaces still work."""
    subprocess.run([CLICK_PATH, *click_args],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)

def run_click(file):
    try:
        warnings.simplefilter('ignore')

        class SelectFragOnly(Select):
            """
            Select pocket residues only
            """

            def accept_residue(self, residue):
                if (residue in fragment):
                    return True  # Accept residue
                else:
                    return False  # Discard residue

        class SelectFragOnly2(Select):
            """
            Select pocket residues only
            """

            def accept_residue(self, residue):
                if (residue in new_peptide):
                    return True  # Accept residue
                else:
                    return False  # Discard residue

        if fnmatch.fnmatch(file, 'minipocket_*.pdb'):
            folder_file = folder_minipockets + "/" + file
            folder = "_".join(file.split("_")[1:4])
            # print(folder.split("_")[0], working_directory)
            if folder.split("_")[0] not in working_directory:
                peptide_chain = folder.split("_")[-1]
                Faufile_noExtension = fau_file.replace(".pdb", "")
                file_noExtension = file.replace(".pdb", "")
                patch = file_noExtension.split("_")[-1].split("-")
                patch_length = len(patch)
                log_file1 = (f"{file_noExtension}-{Faufile_noExtension}.pdb.1.clique")
                out_file1 = (f"{file_noExtension}-{Faufile_noExtension}.1.pdb")
                sup_needed1 = (f"{Faufile_noExtension}-{file_noExtension}.1.pdb")
                # print("-------", folder, peptide_chain, patch, patch_length, "-------")
                shutil.copy(folder_file, ".")

                run_click_binary(file, f"{Faufile_noExtension}.pdb")

                pass_criteria = 0
                if os.path.exists(log_file1):
                    with open(log_file1) as f:
                        for line in f.readlines():
                            try:
                                line = line.strip()
                                if "Overlap" in line:
                                    overlap = float(
                                        line.split("=")[1].split(" ")[-1])
                                    if overlap < 100:
                                            remove_files(log_file1, out_file1, sup_needed1, file)
                                            break
                                    if overlap >= 100:
                                        pass_criteria += 1
                                if "RMSD" in line:
                                    RMSD_score = float(line.split("=")[1].split(" ")[-1])
                                    if RMSD_score > RMSD_allowed:
                                        remove_files(log_file1, out_file1, sup_needed1, file)
                                        break
                                    if RMSD_score <= RMSD_allowed:
                                        pass_criteria += 1
                                if "The number of matched atoms" in line:
                                    matched_Ca = float(line.split("=")[1].split(" ")[-1])
                                    if matched_Ca < cutoff:
                                        remove_files(log_file1, out_file1, sup_needed1, file)
                                        break
                                    if matched_Ca >= cutoff:
                                        pass_criteria += 1
                            except:
                                pass

            if pass_criteria >= 3:
                try:
                    name_outfile1_supneeded1 = f'merge_{sup_needed1.replace(".1.pdb", "")}'
                    strip_end_records(out_file1)
                    concatenate(f"{name_outfile1_supneeded1}.pdb", out_file1, sup_needed1)
                    
                    run_click_binary(f"{name_outfile1_supneeded1}.pdb", file)

                    log_file_tmp = (f"{name_outfile1_supneeded1}-{file_noExtension}.pdb.1.clique")
                    out_file1_tmp = (f"{file_noExtension}-{name_outfile1_supneeded1}.1.pdb")
                    sup_needed1_tmp = (f"{name_outfile1_supneeded1}-{file_noExtension}.1.pdb")
                    leftovers = (log_file1, out_file1_tmp, log_file_tmp, out_file1,
                                 sup_needed1_tmp, f"{name_outfile1_supneeded1}.pdb")
                    # Write to a temporary and only then replace sup_needed1. The old
                    # `grep ... > sup_needed1` truncated the file before grep ran, so a failed
                    # CLICK merge left a 0-byte pocket behind and out_file1 (the only other copy)
                    # was deleted on the next line. Falling back to the stale sup_needed1 is not an
                    # option either: it is in the target frame, while what follows needs the
                    # donor-complex frame. Skip the minipocket instead.
                    if not os.path.exists(sup_needed1_tmp) or os.path.getsize(sup_needed1_tmp) == 0:
                        remove_files(*leftovers, sup_needed1, file)
                        return
                    # grep ' p': keep the pocket chain. grep newline-terminates its last line.
                    with open(f"{sup_needed1}.new", "wb") as handle:
                        handle.writelines(line if line.endswith(b"\n") else line + b"\n"
                                          for line in read_lines(sup_needed1_tmp) if b" p" in line)
                    if os.path.getsize(f"{sup_needed1}.new") == 0:
                        remove_files(*leftovers, f"{sup_needed1}.new", sup_needed1, file)
                        return
                    os.replace(f"{sup_needed1}.new", sup_needed1)
                    remove_files(*leftovers)
                except:
                    pass

                complex_file = (f"{pepbdb_folder}/{folder}/peptide_complex.pdb")
                # print("complex_file", complex_file)
                if not os.path.exists(complex_file):
                    return # Skip this minipocket if DB entry is missing
                
                fragment = []
                allowed_coord_X_s = x_center - (x_size/2)
                allowed_coord_X_e = x_center + (x_size/2)
                allowed_coord_Y_s = y_center - (y_size/2)
                allowed_coord_Y_e = y_center + (y_size/2)
                allowed_coord_Z_s = z_center - (z_size/2)
                allowed_coord_Z_e = z_center + (z_size/2)

                for res in patch:
                    res_id = res[3:]
                    res_name = res[:3]
                    try:
                        parser = PDBParser() # Re-instantiate to avoid concurrency issues?
                        structure = parser.get_structure(file, complex_file)
                        io = PDBIO()
                        io.set_structure(structure)
                        for model in structure:
                            for chain in model:
                                if chain.get_id() == peptide_chain:
                                    chain.id = "x"
                                    for residue in chain:
                                        if str(residue.get_id()[1]) == res_id:
                                            if str(residue.get_resname()) == res_name:
                                                fragment.append(residue)
                    except Exception as e:
                        print(f"DEBUG: Error parsing complex file {complex_file}: {e}")
                
                if len(fragment) == patch_length:
                    sup_needed2 = None
                    try:
                        patch_file = (f'patch_file_{folder}_{"-".join(patch)}.pdb')
                        io.save(patch_file, select=SelectFragOnly())
                        minipocket_file_noExtension = file_noExtension
                        strip_end_records(sup_needed1)
                        super_file = (f"super_{minipocket_file_noExtension}.pdb")
                        # Truncate, do not append: `>>` reused a stale super_file left behind by an
                        # aborted run, silently merging two different superpositions.
                        concatenate(super_file, sup_needed1, patch_file)
                        
                        run_click_binary(super_file, f"{Faufile_noExtension}.pdb")

                        out_file2 = (f'{Faufile_noExtension}-{super_file.replace(".pdb", "")}.1.pdb')
                        sup_needed2 = (f'{super_file.replace(".pdb", "")}-{Faufile_noExtension}.1.pdb')
                        log_file2 = (f'{super_file.replace(".pdb", "")}-{Faufile_noExtension}.pdb.1.clique')
                        strip_end_records(sup_needed2)
                        remove_files(out_file2, log_file2, sup_needed1, super_file, patch_file, file)
                    except Exception as e:
                        print(f"DEBUG: superposition step failed for {file}: {e}")

                    if sup_needed2 is None:
                        # The try above aborted before naming the superposed file. Reading it below
                        # would raise NameError, and so would the handler that reports it.
                        return

                    fragment = []
                    pocket = []
                    fragment_output_lumen = []
                    try:
                        parser = PDBParser()
                        structure = parser.get_structure(
                            sup_needed2.replace(".pdb", ""), sup_needed2)
                        io = PDBIO()
                        io.set_structure(structure)
                        for model in structure:
                            for chain in model:
                                if chain.get_id() == "x":
                                    for residue in chain:
                                        fragment.append(residue)
                                if chain.get_id() == "p":
                                    for residue in chain:
                                        pocket.append(residue)
                    except Exception as e:
                        print(f"DEBUG: Error loading final superposed file {sup_needed2}: {e}")

                    remove_files(sup_needed2)
                    new_peptide = []
                    if not pocket:
                        # No chain 'p' was recovered from the superposed file, so there is nothing to
                        # check clashes against. The test below is vacuously true on an empty distance
                        # list, which would wave every residue through unchecked; bail out instead.
                        return
                    for residue_pep in fragment:
                        dist_pocket_list = []
                        for residue_poc in pocket:
                            for atoms_poc in residue_poc:
                                for atoms_pep in residue_pep:
                                    coord_res_pep = atoms_pep.get_coord()
                                    coord_res_poc = atoms_poc.get_coord()
                                    dist_pocket_pep = np.linalg.norm(
                                        coord_res_pep - coord_res_poc)
                                    dist_pocket_list.append(
                                        dist_pocket_pep)
                        dist_pocket_pass = [
                            dist_pocket for dist_pocket in dist_pocket_list if dist_pocket >= 1.8]
                        if len(dist_pocket_pass) == len(dist_pocket_list):  # DISTANCIA
                            new_peptide.append(residue_pep)

                    # Debug Distance Filter
                    # if len(new_peptide) == 0:
                         # print(f"DEBUG: {file} failed distance filter (Clash check).")

                    for residue_pep2 in new_peptide:
                        try:
                            coord_res_pep = residue_pep2["CA"].get_coord()
                            x_coord_res_pep = coord_res_pep[0]
                            y_coord_res_pep = coord_res_pep[1]
                            z_coord_res_pep = coord_res_pep[2]
                            if x_coord_res_pep > allowed_coord_X_s and x_coord_res_pep < allowed_coord_X_e:
                                if y_coord_res_pep > allowed_coord_Y_s and y_coord_res_pep < allowed_coord_Y_e:
                                    if z_coord_res_pep > allowed_coord_Z_s and z_coord_res_pep < allowed_coord_Z_e:
                                        res_nameid = residue_pep2.get_resname() + str(residue_pep2.get_id()[1])
                                        fragment_output_lumen.append(res_nameid)
                                    # else: print(f"DEBUG: Z out of bounds {z_coord_res_pep}")
                                # else: print(f"DEBUG: Y out of bounds {y_coord_res_pep}")
                            # else: print(f"DEBUG: X out of bounds {x_coord_res_pep}")
                        except:
                            pass
                    fragment_output_lumen = set(fragment_output_lumen)
                    error_lumen = len(fragment) - len(fragment_output_lumen)
                    
                    if error_lumen <= 0 and len(new_peptide) > 0:
                        patch = [res.get_resname() + str(res.get_id()[1])
                                for res in new_peptide]
                        # Include the source DB entry in the name. Without it, every database hit
                        # whose surviving residues share names and numbers wrote to the same path
                        # and the last writer won, silently discarding distinct superpositions.
                        # The patch_file_ prefix is load-bearing: four downstream globs match on it.
                        patch_file2 = ("../{out_folder}/patch_file_{folder}_{patch}.pdb").format(
                            out_folder=folder_output, folder=folder, patch="-".join(patch))

                        io.save(patch_file2, select=SelectFragOnly2())

    except:
        pass

os.chdir(folder_temp)
Parallel(n_jobs=threads)(delayed(run_click)(file) for file in tqdm(os.listdir(folder_minipockets), total=len(
    os.listdir(folder_minipockets)), desc="aligning minipockets to target_pocket and defining patches"))
os.chdir("../")
shutil.rmtree(folder_temp, ignore_errors=True)
