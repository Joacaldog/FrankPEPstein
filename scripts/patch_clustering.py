import os
import time
from tqdm import tqdm
import fnmatch
from Bio.PDB import *
import Bio.PDB
pdbio = Bio.PDB.PDBIO()
parser = Bio.PDB.MMCIFParser()
import glob
import warnings
import itertools
from Bio.PDB.StructureBuilder import *
import copy
import numpy as np
warnings.simplefilter('ignore')
import argparse
from joblib import Parallel, delayed
from string import ascii_lowercase
from itertools import product
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster, complete
from collections import Counter, OrderedDict
import statistics
import math
from operator import itemgetter
from sklearn.linear_model import LinearRegression
import random

# Configuration Variables
# MAX_COMBINATIONS moved to argparse
OUTLIER_THRESHOLD = 3.5

tqdm._instances.clear()

program_description = "Constructs peptides from patches"
parser = argparse.ArgumentParser(description=program_description)
parser.add_argument("-w", "--winsize", type=int,
                    help="Peptide lenght", required=True)
parser.add_argument("-t", "--threads", type=int,
                    help="Number of threads to utilize", required=True)
parser.add_argument("-c", "--max_combinations", type=int,
                    help="Max number of combinations to sample", required=False, default=100)
args = parser.parse_args()
winsize = args.winsize
threads = args.threads
MAX_COMBINATIONS = args.max_combinations

d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    
initial_path = os.getcwd()
folder_output = (f"frankPEPstein_{winsize}")
if not os.path.exists(folder_output):
    os.makedirs(folder_output)

def save_pdb(peptide_ordered):
    id_pep = "pep"
    builder = StructureBuilder()
    builder.init_structure(id_pep)
    builder.init_seg(" ")
    builder.init_model(0)
    builder.init_chain("x")
    final_peptide2 = copy.deepcopy(peptide_ordered)
    j = 0
    for residue in final_peptide2:
        j += 1
        res_name = residue.get_resname()
        builder.init_residue(res_name, " ", j, " ")
        res = builder.residue
        for atom in residue.get_atoms():
            res.add(atom)
    structure2 = builder.get_structure()
    io = PDBIO()
    io.set_structure(structure2)
    peptide_name = "pep_" + \
        "".join([d3to1.get(res.get_resname())
                 for res in final_peptide2])
    peptide_file = folder_output + "/" + peptide_name + ".pdb"
    io.save(peptide_file, preserve_atom_numbering=False)
    # time.sleep(0.5)
    return peptide_file


def pdb_parser(peptide_file):
    pdb_file_no_extension = peptide_file.replace(".pdb", "")
    cmd_remove_END = (
        f'cat {peptide_file} 2> /dev/null | grep -v "END" > {pdb_file_no_extension}_noEND.pdb')

    os.system(cmd_remove_END)

    with open(pdb_file_no_extension + "_connects.txt", "w") as outfile:
        bond_list = []
        with open(peptide_file) as f:
            for line in f.readlines():
                atom_number = line[6:11].strip()
                atom_name = line[12:16].strip()
                if atom_name == "N" or atom_name == "C":
                    if atom_number != "1":
                        bond_list.append(atom_number)
        bond_list = bond_list[:-1]
        number_of_bonds = len(bond_list)
        for i in range(0, number_of_bonds, 2):
            first_atom = bond_list[i].rjust(5)
            second_atom = bond_list[i + 1].rjust(5)
            output = ("CONECT{}{}\n").format(
                first_atom, second_atom)
            outfile.write(output)
        outfile.write("END")

    os.system(("rm -f {}.pdb 2> /dev/null").format(pdb_file_no_extension))
    merge_files = ("cat {f1}_noEND.pdb {f1}_connects.txt > {f1}.pdb 2> /dev/null").format(
        f1=pdb_file_no_extension)
    os.system(merge_files)
    os.system(("rm -f {}_noEND.pdb 2> /dev/null").format(pdb_file_no_extension))
    os.system(("rm -f {}_connects.txt 2> /dev/null").format(pdb_file_no_extension))
    return pdb_file_no_extension
    
def kmerizator(pdb_file_no_extension, k_mers_list):
    file = pdb_file_no_extension
    class SelectFragOnly(Select):
        """
        Select fragment residues only
        """

        def accept_residue(self, residue):
            if (residue in fragment):
                return True  # Accept residue
            else:
                return False  # Discard residue
            
    # Parse filename command line arguments with DOS/Unix-save wildcard treatment

    filenames = []
    file = f"{file}"
    filenames.append(file)
    # ​
    # Loop over all PDB files
    for filename in filenames:
        try:
            # ​
            # Split filename into /path/basename.ext
            (path, ext) = os.path.splitext(filename)
            basename = os.path.basename(path)
        # ​
        # Open file
            if ext == '.gz':
                fp = gzip.open(filename)
            else:
                fp = open(filename, 'r')
        # ​
        # Parse PDB
            parser = PDBParser()
            structure = parser.get_structure(basename, fp)

            # Prepare writer
            io = PDBIO()
            io.set_structure(structure)
        # Iterate over PDB structure
            for model in structure:
                for chain in model:
                    residues = [residue for residue in chain]  # Get residues as list
                    # Calculate last index permitting full windows size
                    if len(residues) >= winsize:
                        endidx = len(residues) - winsize + 1
                        # print(endidx)
                        # Iterate over residues
                        for residx, residue in enumerate(residues[:endidx]):
                            fragment = []  # Empty protein fragment (peptide)
                            # Iterate over protein fragment of size winsize
                            for residx2, residue in enumerate(residues[residx:residx + winsize]):
                                fragment.append(residue)  # Add residue to protein fragment
                            # Write protein fragment to disk (in PDB format)
                            kmer_name = "".join([d3to1.get(res.get_resname()) for res in fragment])
                            if kmer_name not in k_mers_list:
                                k_mers_list.append(kmer_name)
                                io.save(f'{folder_output}/frag_{kmer_name}.pdb', select=SelectFragOnly())
                                pdb_parser(f'{folder_output}/frag_{kmer_name}.pdb')
                            # os.system(f'rm sequence_analysis_{kmer_name}.txt 2> /dev/null')


                    if len(residues) < winsize:
                        fragment = residues
                        kmer_name = "".join([d3to1.get(res.get_resname()) for res in fragment])
                        io.save(f'{folder_output}/frag_{kmer_name}.pdb', select=SelectFragOnly())
                        pdb_parser(f'{folder_output}/frag_{kmer_name}.pdb')
                        # os.system(f'rm sequence_analysis_{kmer_name}.txt 2> /dev/null')
                        

        except:
            continue
  

def find_outliers_linear_trend(coordinates_list, threshold_multiplier=OUTLIER_THRESHOLD):
    # Separar las coordenadas en X, Y, y Z
    x_coordinates = np.array([coord[0][0] for coord in coordinates_list])
    y_coordinates = np.array([coord[0][1] for coord in coordinates_list])
    z_coordinates = np.array([coord[0][2] for coord in coordinates_list])

    # Ajustar una regresión lineal a las coordenadas
    def fit_linear_model(x, y):
        model = LinearRegression()
        model.fit(x.reshape(-1, 1), y)
        return model

    x_model = fit_linear_model(np.arange(len(x_coordinates)), x_coordinates)
    y_model = fit_linear_model(np.arange(len(y_coordinates)), y_coordinates)
    z_model = fit_linear_model(np.arange(len(z_coordinates)), z_coordinates)

    # Calcular los residuos de cada coordenada con respecto a su modelo lineal
    x_residuals = x_coordinates - x_model.predict(np.arange(len(x_coordinates)).reshape(-1, 1))
    y_residuals = y_coordinates - y_model.predict(np.arange(len(y_coordinates)).reshape(-1, 1))
    z_residuals = z_coordinates - z_model.predict(np.arange(len(z_coordinates)).reshape(-1, 1))

    # Calcular la mediana de los residuos para cada dimensión
    x_median_residual = np.median(x_residuals)
    y_median_residual = np.median(y_residuals)
    z_median_residual = np.median(z_residuals)

    # Identificar los outliers basados en los residuos
    outliers = []
    for i in range(len(coordinates_list)):
        coord, name = coordinates_list[i]
        is_outlier = (
            abs(x_residuals[i] - x_median_residual) > threshold_multiplier * np.std(x_residuals) or
            abs(y_residuals[i] - y_median_residual) > threshold_multiplier * np.std(y_residuals) or
            abs(z_residuals[i] - z_median_residual) > threshold_multiplier * np.std(z_residuals)
        )
        if is_outlier:
            outliers.append((coord, name))

    return outliers


def delete_outsider_frag():
    coordinates_3d = []
    for file in tqdm(os.listdir("."), total=len(os.listdir(".")), desc="loading structures", position=0, leave=True):
        warnings.simplefilter('ignore')
        if fnmatch.fnmatch(file, 'patch_file_*.pdb'):
            name_file = file.replace(".pdb", "")
            parser = PDBParser()
            structure = parser.get_structure(name_file, file)
            for model in structure:
                for chain in model:
                    ca_list = []
                    for residue in chain:
                        ca_coord = residue["CA"].coord
                        ca_list.append(ca_coord)
                    array = list(np.average(np.array(ca_list), axis=0))
                    coords_array_name = (array, file)
                    coordinates_3d.append(coords_array_name)
    return coordinates_3d

if len([x for x in os.listdir(".") if "patch_file" in x]) > 2:
    coordinates_3d = delete_outsider_frag()
    outliers_3d = find_outliers_linear_trend(coordinates_3d)
    if len(outliers_3d) > 0:
        if not os.path.exists("outlier_folder"):
            os.makedirs("outlier_folder")
        for outlier in outliers_3d:
            file = outlier[1]
            # print(f"mv {file} outlier_folder")
            os.system(f"mv {file} outlier_folder")

def combinator():
    res_dict = {}
    res_dict_pull = {}
    final_peptides_list = [] #ACA LISTA FINAL DE PEPTIDOS ORDENADOS
    if len([x for x in os.listdir(".") if "patch_file" in x]) == 0:
        print("No patches files in folder")
    if len([x for x in os.listdir(".") if "patch_file" in x]) <= 1 and len([x for x in os.listdir(".") if "patch_file" in x]) > 0:
        cmd_cp = (f"cp patch_file* {folder_output}")
        os.system(cmd_cp)
    if len([x for x in os.listdir(".") if "patch_file" in x]) > 1:
        for file in tqdm(os.listdir("."), total=len(os.listdir(".")), desc="loading structures", position=0, leave=True):
            try:
                warnings.simplefilter('ignore')
                if fnmatch.fnmatch(file, 'patch_file_*.pdb'):
                    name_file = file.replace(".pdb", "")
                    parser = PDBParser()
                    structure = parser.get_structure(name_file, file)
                    for model in structure:
                        for chain in model:
                            for residue in chain:
                                resnameid_coord = residue.get_resname(
                                ) + str(residue.get_id()[1]) + "_" + "_".join([str(tuple(residue["N"].coord)), str(tuple(residue["CA"].coord)), str(tuple(residue["C"].coord))])
                                if res_dict.get(resnameid_coord) is None:
                                    res_dict[resnameid_coord] = []
                                res_dict[resnameid_coord].append(residue)
                                if res_dict_pull.get(resnameid_coord) is None:
                                    res_dict_pull[resnameid_coord] = residue
            except:
                continue

        duplicate_dict = {}
        for residue_list in res_dict.values():
            for residue in residue_list:
                resname = residue.get_resname()
                if duplicate_dict.get(resname) is None:
                    duplicate_dict[resname] = []
                duplicate_dict[resname].append(residue)
        reformed_residues = []
        for resname, residues_list in duplicate_dict.items():
            res_cluster_dict = {}
            if len(residues_list) > 1:
                all_coords_array = [(np.average(np.array(
                    [residues["N"].coord, residues["CA"].coord, residues["C"].coord]), axis=0)) for residues in residues_list]
                dist_matrix = pdist(all_coords_array, metric='euclidean')
                complete_dist_matrix = complete(dist_matrix)
                clusters = fcluster(complete_dist_matrix,
                                    1, criterion='distance')
                clustered_orderer_res = [(y, "C" + str(x)) for x, y in sorted(
                    zip(clusters, residues_list))]
                for residue, order in clustered_orderer_res:
                    if res_cluster_dict.get(order) is None:
                        res_cluster_dict[order] = []
                    res_cluster_dict[order].append(residue)
                for cluster_residues in res_cluster_dict.values():
                    cluster_residues_str = [residue.get_resname() + str(residue.get_id()[1]) + "_" + "_".join([str(tuple(
                        residue["N"].coord)), str(tuple(residue["CA"].coord)), str(tuple(residue["C"].coord))]) for residue in cluster_residues]
                    counter = Counter(cluster_residues_str)
                    cluster_residues = res_dict_pull.get(
                        list(counter.keys())[0])
                    reformed_residues.append(cluster_residues)
            if len(residues_list) == 1:
                reformed_residues.append(residues_list[0])

        res_cluster_dict = {}
        all_coords_array = [(np.average(np.array([residues["N"].coord, residues["CA"].coord,residues["C"].coord]), axis=0)) for residues in reformed_residues]
        dist_matrix = pdist(all_coords_array, metric='euclidean')
        complete_dist_matrix = complete(dist_matrix)
        clusters = fcluster(complete_dist_matrix,
                            3, criterion='distance')
        clustered_orderer_res = [(y, "C" + str(x)) for x, y in sorted(
            zip(clusters, reformed_residues))]
        for residue, order in clustered_orderer_res:
            if res_cluster_dict.get(order) is None:
                res_cluster_dict[order] = []
            res_cluster_dict[order].append(residue)

        equivalentRes_dup = []
        single_res = []
        for order, equivalent_residues in res_cluster_dict.items():
            if len(equivalent_residues) > 1:
                equivalentRes_dup.append(equivalent_residues)
            if len(equivalent_residues) == 1:
                residue = equivalent_residues[0]
                single_res.append(residue)

        dup_comb_list = []
        print("initiating combinations...")
        combination_res_list = product(*equivalentRes_dup)
        for combination in combination_res_list:
            combination = list(set(combination))
            dup_comb_list.append(combination)

        combination_res_list = list(product([single_res], dup_comb_list))
        print("combinations finished...")
        final_names_list = []
        len_comb = len(combination_res_list)
        if len_comb >= MAX_COMBINATIONS:
            combination_res_list = random.sample(combination_res_list, MAX_COMBINATIONS)
        for final_comb in tqdm(combination_res_list, total=len(combination_res_list), desc="Working on each combination...", position=0, leave=True):
            res_dist_dict = {}
            res_dist_list = []
            final_peptide = list(itertools.chain.from_iterable(
                [final_comb[0], final_comb[1]]))
            nameresid_list = list(
                set([res.get_resname(
                ) + str(res.get_id()[1]) + "_" + "_".join([str(tuple(res["N"].coord)), str(tuple(res["CA"].coord)), str(tuple(res["C"].coord))]) for res in final_peptide]))
            final_peptide = [res_dict_pull.get(nameresid)
                                for nameresid in nameresid_list]
            combination_NC_list = list(itertools.permutations(
                final_peptide, 2))
            for CA_prox in combination_NC_list:
                residue1 = CA_prox[0]
                residue2 = CA_prox[1]
                residue1_CA_coord = (np.average(np.array([residue1["N"].coord, residue1["CA"].coord,
                                                            residue1["C"].coord]), axis=0))
                residue2_CA_coord = (np.average(np.array([residue2["N"].coord, residue2["CA"].coord,
                                                            residue2["C"].coord]), axis=0))
                dist_CA = np.linalg.norm(
                    residue1_CA_coord - residue2_CA_coord)
                if res_dist_dict.get(residue1) is None:
                    res_dist_dict[residue1] = []
                res_dist_dict[residue1].append(dist_CA)

            for res, dist_list in res_dist_dict.items():
                dist_max = max(dist_list)
                # dist_min = min(dist_list)
                # st_dist = np.std([dist_max,dist_min])
                tuple_res = res, dist_max
                res_dist_list.append(tuple_res)

            sorted_list = sorted(res_dist_list,key=itemgetter(1))
            start_res = [sorted_list[-1][0]]
            end_res = sorted_list[-2][0]
            mid_residues = [x[0] for x in sorted_list[:-2]]
            prod_res = product(start_res, mid_residues)
            peptide_ordered = []
            dist_bond_pass = 0
            prod_rest_list = []
            for res_pair in prod_res:
                start_res = res_pair[0]
                next_res = res_pair[1]
                dist_res = next_res["CA"] - start_res["CA"]
                triple_data = start_res, next_res, dist_res
                prod_rest_list.append(triple_data)
            sorted_list = sorted(prod_rest_list,key=itemgetter(2))
            start_res = sorted_list[0][0]
            next_res = sorted_list[0][1]
            peptide_ordered.append(start_res)
            peptide_ordered.append(next_res)
            mid_residues.remove(next_res)
            dist_bond1 = next_res["N"] - start_res["C"]
            dist_bond2 = start_res["N"] - next_res["C"]
            if dist_bond1 < dist_bond2:
                dist_bond_pass += 1

                    
            for i in range(1, len(mid_residues)):
                prod_rest_list2 = []  
                prev_res = [peptide_ordered[-1]]
                prod_res = product(prev_res, mid_residues)
                for res_pair in prod_res:
                    start_res = res_pair[0]
                    next_res = res_pair[1]
                    dist_res = next_res["CA"] - start_res["CA"]
                    triple_data = start_res, next_res, dist_res
                    prod_rest_list2.append(triple_data)
                sorted_list = sorted(prod_rest_list2,key=itemgetter(2))
                next_res = sorted_list[0][1]
                peptide_ordered.append(next_res)
                mid_residues.remove(next_res)
            peptide_ordered.append(end_res)

            # if len(peptide_ordered) >= 5:
            if dist_bond_pass == 1:
                peptide_name = "".join([d3to1.get(res.get_resname()) for res in peptide_ordered])
                if peptide_name not in final_names_list:
                    final_names_list.append(peptide_name)
                    final_peptides_list.append(peptide_ordered)
            if dist_bond_pass == 0:
                peptide_ordered = peptide_ordered[::-1]
                peptide_name = "".join([d3to1.get(res.get_resname()) for res in peptide_ordered])
                if peptide_name not in final_names_list:
                    final_names_list.append(peptide_name)
                    final_peptides_list.append(peptide_ordered)

    return final_peptides_list
    
def main():
    final_peptides_list = combinator()
    if len(final_peptides_list) > 0:
        def frag_min(peptide_ordered, k_mers_list):
            peptide_file = save_pdb(peptide_ordered)
            kmerizator(peptide_file, k_mers_list)
            os.system(f"rm -f {peptide_file} 2> /dev/null")

        k_mers_list = []
        Parallel(n_jobs=threads, backend="threading")(delayed(frag_min)(peptide_ordered, k_mers_list) for peptide_ordered in tqdm(final_peptides_list, total=len(final_peptides_list), 
                                                                                         desc=f"generating peptides of length {winsize}", position=0, leave=True))
    os.chdir(folder_output)
    os.system("rm *noEND.pdb *connects.txt 2> /dev/null")
        
if __name__ == '__main__':
    main()
