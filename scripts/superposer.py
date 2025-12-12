#GRIDBOX ADDED

import os
import sys
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

# Configuration Variables
PARAMETERS_INP_PATH = "~/FrankPEPstein/utilities/Parameters.inp"
CLICK_PATH = "~/FrankPEPstein/utilities/click"

program_description = "Select and generate fragment of peptides that could eventually bind to target receptor based on minipocket alignments"
parser = argparse.ArgumentParser(description=program_description)
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

args = parser.parse_args()
fau_file = args.target_receptor
pepbdb_folder = args.pepbdb_folder
cutoff = 3
threads = args.threads
folder_minipockets = args.folder_minipockets
RMSD_allowed = 0.1
x_center = args.x_center
y_center = args.y_center
z_center = args.z_center
x_size = args.x_size
y_size = args.y_size
z_size = args.z_size

folder_temp = (f"temp_folder_residuesAligned{cutoff}_RMSD{RMSD_allowed}")
if not os.path.exists(folder_temp):
    os.makedirs(folder_temp)

folder_output = (f"superpockets_residuesAligned{cutoff}_RMSD{RMSD_allowed}")
if not os.path.exists(folder_output):
    os.makedirs(folder_output)

working_directory = os.getcwd()

os.system(f'cp {fau_file} {folder_temp}')
os.system(f"cp {PARAMETERS_INP_PATH} {folder_temp}")

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
                os.system(f"cp {folder_file} .")
                cmd = (f"{CLICK_PATH} {file} {Faufile_noExtension}.pdb 1> /dev/null 2> /dev/null")
                os.system(cmd)
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
                                            os.system(f'rm {log_file1} {out_file1} {sup_needed1} {file} 2> /dev/null')
                                            break
                                    if overlap >= 100:
                                        pass_criteria += 1
                                if "RMSD" in line:
                                    RMSD_score = float(line.split("=")[1].split(" ")[-1])
                                    if RMSD_score > RMSD_allowed:
                                        os.system(f'rm {log_file1} {out_file1} {sup_needed1} {file} 2> /dev/null')
                                        break
                                    if RMSD_score <= RMSD_allowed:
                                        pass_criteria += 1
                                if "The number of matched atoms" in line:
                                    matched_Ca = float(line.split("=")[1].split(" ")[-1])
                                    if matched_Ca < cutoff:
                                        os.system(f'rm {log_file1} {out_file1} {sup_needed1} {file} 2> /dev/null')
                                        break
                                    if matched_Ca >= cutoff:
                                        pass_criteria += 1
                            except:
                                pass

            if pass_criteria >= 3:
                try:
                    name_outfile1_supneeded1 = f'merge_{sup_needed1.replace(".1.pdb", "")}'
                    os.system(f"sed -i '/END/d' {out_file1}")
                    os.system(f'cat {out_file1} {sup_needed1} > {name_outfile1_supneeded1}.pdb')
                    cmd = (f"{CLICK_PATH} {name_outfile1_supneeded1}.pdb {file} 1> /dev/null 2> /dev/null")
                    os.system(cmd)
                    log_file_tmp = (f"{name_outfile1_supneeded1}-{file_noExtension}.pdb.1.clique")
                    out_file1_tmp = (f"{file_noExtension}-{name_outfile1_supneeded1}.1.pdb")
                    sup_needed1_tmp = (f"{name_outfile1_supneeded1}-{file_noExtension}.1.pdb")
                    os.system(f"cat {sup_needed1_tmp} | grep ' p' > {sup_needed1}")
                    os.system(f'rm {log_file1} {out_file1_tmp} {log_file_tmp} {out_file1} {sup_needed1_tmp} {name_outfile1_supneeded1}.pdb 2> /dev/null')
                except:
                    pass
                complex_file = (f"{pepbdb_folder}/{folder}/peptide_complex.pdb")
                # print("complex_file", complex_file)
                fragment = []
                allowed_coord_X_s = x_center - (x_size/2 * 0.385)
                allowed_coord_X_e = x_center + (x_size/2 * 0.385)
                allowed_coord_Y_s = y_center - (y_size/2 * 0.385)
                allowed_coord_Y_e = y_center + (y_size/2 * 0.385)
                allowed_coord_Z_s = z_center - (z_size/2 * 0.385)
                allowed_coord_Z_e = z_center + (z_size/2 * 0.385)

                for res in patch:
                    res_id = res[3:]
                    res_name = res[:3]
                    parser = PDBParser()
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
                if len(fragment) == patch_length:
                    try:
                        patch_file = (f'patch_file_{folder}_{"-".join(patch)}.pdb')
                        io.save(patch_file, select=SelectFragOnly())
                        minipocket_file_noExtension = file_noExtension
                        os.system(f"sed -i '/END/d' {sup_needed1}")
                        super_file = (f"super_{minipocket_file_noExtension}.pdb")
                        cmd_merge_peptide_pocket = f"cat {sup_needed1} {patch_file} >> {super_file}"
                        os.system(cmd_merge_peptide_pocket)
                        cmd_super = (f'{CLICK_PATH} {super_file} {Faufile_noExtension}.pdb 1> /dev/null')
                        os.system(cmd_super)
                        out_file2 = (f'{Faufile_noExtension}-{super_file.replace(".pdb", "")}.1.pdb')
                        sup_needed2 = (f'{super_file.replace(".pdb", "")}-{Faufile_noExtension}.1.pdb')
                        log_file2 = (f'{super_file.replace(".pdb", "")}-{Faufile_noExtension}.pdb.1.clique')
                        os.system(f"sed -i '/END/d' {sup_needed2}")
                        os.system(f"rm {out_file2} {log_file2} {sup_needed1} {super_file} {patch_file} {file}")
                    except:
                        pass
                    
                    fragment = []
                    pocket = []
                    fragment_output_lumen = []
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
                    try:
                        os.system(f"rm {sup_needed2}")
                    except:
                        pass
                    new_peptide = []
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
                        except:
                            pass
                    fragment_output_lumen = set(fragment_output_lumen)
                    error_lumen = len(fragment) - len(fragment_output_lumen)
                    if error_lumen <= 0 and len(new_peptide) > 0:
                        patch = [res.get_resname() + str(res.get_id()[1])
                                for res in new_peptide]
                        patch_file2 = ("../{out_folder}/patch_file_{patch}.pdb").format(
                            out_folder=folder_output, patch="-".join(patch))
                        io.save(patch_file2, select=SelectFragOnly2())
    except:
        pass

os.chdir(folder_temp)
Parallel(n_jobs=threads)(delayed(run_click)(file) for file in tqdm(os.listdir(folder_minipockets), total=len(
    os.listdir(folder_minipockets)), desc="aligning minipockets to target_pocket and defining patches"))
os.chdir("../")
os.system("rm -r " + folder_temp)
