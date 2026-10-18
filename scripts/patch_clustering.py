import os
import sys
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
parser.add_argument("-seed", "--seed", type=int,
                    help="Seed for the combination sampling, for reproducible runs", required=False, default=0)
args = parser.parse_args()
winsize = args.winsize
threads = args.threads
MAX_COMBINATIONS = args.max_combinations
random.seed(args.seed)

d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

_warned_resnames = set()


def usable_residue(res):
    """A residue this pipeline can name with one letter and place by its backbone.

    Anything else (modified residues, ligands, waters, caps) is filtered out at ingest rather
    than defaulted to 'X'. An 'X' would flow into the peptide name, and step 3 recovers the
    final sequence from the filename with a regex over the 20 standard letters, so a single X
    silently truncates the reported sequence instead of failing.
    """
    return res.get_resname() in d3to1 and all(atom in res for atom in ("N", "CA", "C"))


def warn_unusable(res):
    name = res.get_resname()
    if name not in _warned_resnames:
        _warned_resnames.add(name)
        print(f"Warning: skipping residue type {name!r} (not a standard residue, or missing N/CA/C).")


def sequence_of(residues):
    """One-letter sequence, or None if any residue is not representable."""
    letters = [d3to1.get(res.get_resname()) for res in residues]
    return None if any(letter is None for letter in letters) else "".join(letters)


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
    # Copy atoms only. A deepcopy of a residue follows its parent links and duplicates the whole
    # patch structure; the copy is only there so res.add() re-parents a copy, not the shared atom.
    j = 0
    for residue in peptide_ordered:
        j += 1
        res_name = residue.get_resname()
        builder.init_residue(res_name, " ", j, " ")
        res = builder.residue
        for atom in residue.get_atoms():
            res.add(atom.copy())
    structure2 = builder.get_structure()
    io = PDBIO()
    io.set_structure(structure2)
    sequence = sequence_of(peptide_ordered)
    if sequence is None:
        return None
    peptide_name = "pep_" + sequence
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
                            kmer_name = sequence_of(fragment)
                            if kmer_name is not None and kmer_name not in k_mers_list:
                                k_mers_list.append(kmer_name)
                                io.save(f'{folder_output}/frag_{kmer_name}.pdb', select=SelectFragOnly())
                                pdb_parser(f'{folder_output}/frag_{kmer_name}.pdb')
                            # os.system(f'rm sequence_analysis_{kmer_name}.txt 2> /dev/null')


                    if len(residues) < winsize:
                        fragment = residues
                        kmer_name = sequence_of(fragment)
                        if kmer_name is None:
                            continue
                        io.save(f'{folder_output}/frag_{kmer_name}.pdb', select=SelectFragOnly())
                        pdb_parser(f'{folder_output}/frag_{kmer_name}.pdb')
                        # os.system(f'rm sequence_analysis_{kmer_name}.txt 2> /dev/null')
                        

        except:
            continue
  

def remove_residue(residues, target):
    """Remove `target` itself from `residues`.

    list.remove() matches by ==, and Biopython residue equality compares full_id without the
    structure id, so a residue from another patch file with the same chain and number counts
    as equal. remove() would then drop that other residue, leaving `target` duplicated in the
    assembled peptide and the other one lost.
    """
    del residues[next(i for i, res in enumerate(residues) if res is target)]


def find_outliers_linear_trend(coordinates_list, threshold_multiplier=OUTLIER_THRESHOLD):
    """Flag patch centroids that sit far from the bulk of the patch cloud.

    This used to fit a linear regression of each coordinate against the patch's *index* in the
    list, i.e. against os.listdir order. There is no independent variable here, so the "trend"
    was an artefact of filesystem ordering, and sorting the listing would not help: patch
    filenames encode residue numbers, so alphabetical order correlates with position along the
    donor chain and therefore with position in space. The statistic is now computed from the
    point cloud itself, using a median-absolute-deviation scale with a 1 A floor so it cannot
    fire on a near-degenerate set of almost-identical patches.
    """
    centroids = np.array([coord for coord, _ in coordinates_list])
    distances = np.linalg.norm(centroids - np.median(centroids, axis=0), axis=1)
    median_distance = np.median(distances)
    scale = max(1.4826 * np.median(np.abs(distances - median_distance)), 1.0)
    cut = median_distance + threshold_multiplier * scale
    return [coordinates_list[i] for i in range(len(distances)) if distances[i] > cut]


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
                        if not usable_residue(residue):
                            continue
                        ca_list.append(residue["CA"].coord)
                    if not ca_list:
                        # A patch of only caps/waters would otherwise give a NaN centroid.
                        continue
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
                                # Load-bearing filter: nothing unusable enters res_dict, so it can
                                # never reach reformed_residues, peptide_ordered, the peptide name,
                                # or a frag_*.pdb that frankVINA_2 would have to minimise.
                                if not usable_residue(residue):
                                    warn_unusable(residue)
                                    continue
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
        print("initiating combinations...")
        
        # Calculate total possible combinations
        total_stats = 1
        for group in equivalentRes_dup:
            total_stats *= len(group)
            
        dup_comb_list = []
        
        if total_stats <= MAX_COMBINATIONS:
            # Small enough to generate all
            combination_res_list = product(*equivalentRes_dup)
            for combination in combination_res_list:
                combination = list(set(combination))
                dup_comb_list.append(combination)
        else:
            # Too many: draw MAX_COMBINATIONS DISTINCT combinations. Independent per-group draws
            # sample with replacement, so asking for N combinations silently returned noticeably
            # fewer unique ones. Index into the cartesian product instead.
            print(f"Total combinations ({total_stats}) exceeds limit ({MAX_COMBINATIONS}). Random sampling...")
            group_sizes = [len(group) for group in equivalentRes_dup]

            def decode_combination(index):
                combo = []
                for group, size in zip(equivalentRes_dup, group_sizes):
                    combo.append(group[index % size])
                    index //= size
                return combo

            if total_stats <= sys.maxsize:
                for index in random.sample(range(total_stats), MAX_COMBINATIONS):
                    dup_comb_list.append(list(set(decode_combination(index))))
            else:
                # Astronomically many combinations: collisions are negligible, reject the few.
                seen_indices = set()
                attempts = 0
                while len(dup_comb_list) < MAX_COMBINATIONS and attempts < 20 * MAX_COMBINATIONS:
                    attempts += 1
                    index = tuple(random.randrange(size) for size in group_sizes)
                    if index in seen_indices:
                        continue
                    seen_indices.add(index)
                    dup_comb_list.append(list(set(
                        group[i] for group, i in zip(equivalentRes_dup, index))))

        # Skip the original product line as we built dup_comb_list manually
        combination_res_list = list(product([single_res], dup_comb_list))
        print("combinations finished...")
        final_names_list = []
        len_comb = len(combination_res_list)
        # Random sample again if strictly needed technically (though we just curbed it), 
        # but logic dictates we might have up to MAX_COMBINATIONS now.
        if len_comb >= MAX_COMBINATIONS:
            combination_res_list = random.sample(combination_res_list, MAX_COMBINATIONS)
        short_combinations = 0
        # Backbone centroid of every residue, computed once. The pairwise loop below used to
        # recompute both centroids for every ordered pair of every combination. Keyed by id(),
        # because Biopython residue equality ignores the structure id, so residues from two
        # patch files with the same chain and number compare equal.
        centroids = {id(res): np.average(np.array([res["N"].coord, res["CA"].coord, res["C"].coord]), axis=0)
                     for res in res_dict_pull.values()}
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
                residue1_CA_coord = centroids[id(residue1)]
                residue2_CA_coord = centroids[id(residue2)]
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
            if len(sorted_list) < 3:
                # The 3 A clustering collapsed this combination to fewer than three distinct
                # positions, so there is no start/mid/end split to make. Indexing [-2] here used
                # to raise IndexError and destroy every peptide assembled so far.
                short_combinations += 1
                continue
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
            remove_residue(mid_residues, next_res)
            dist_bond1 = next_res["N"] - start_res["C"]
            dist_bond2 = start_res["N"] - next_res["C"]
            if dist_bond1 < dist_bond2:
                dist_bond_pass += 1

                    
            # Consume mid_residues until empty. The previous `range(1, len(mid_residues))` bound
            # discounted the residue already removed above, so the loop ran one iteration short and
            # every assembled peptide silently lost its last middle residue.
            while mid_residues:
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
                remove_residue(mid_residues, next_res)
            peptide_ordered.append(end_res)

            # if len(peptide_ordered) >= 5:
            if dist_bond_pass == 0:
                peptide_ordered = peptide_ordered[::-1]
            peptide_name = sequence_of(peptide_ordered)
            if peptide_name is not None and peptide_name not in final_names_list:
                final_names_list.append(peptide_name)
                final_peptides_list.append(peptide_ordered)

        if short_combinations:
            print(f"Warning: skipped {short_combinations} combination(s) with fewer than 3 distinct "
                  f"residue positions - the 3 A clustering collapsed them. Patches may be too few "
                  f"or too similar.")

    return final_peptides_list
    
def main():
    final_peptides_list = combinator()
    if len(final_peptides_list) > 0:
        def frag_min(peptide_ordered, k_mers_list):
            peptide_file = save_pdb(peptide_ordered)
            if peptide_file is None:
                return
            kmerizator(peptide_file, k_mers_list)
            os.system(f"rm -f {peptide_file} 2> /dev/null")

        k_mers_list = []
        Parallel(n_jobs=threads, backend="threading")(delayed(frag_min)(peptide_ordered, k_mers_list) for peptide_ordered in tqdm(final_peptides_list, total=len(final_peptides_list), 
                                                                                         desc=f"generating peptides of length {winsize}", position=0, leave=True))
    os.chdir(folder_output)
    os.system("rm *noEND.pdb *connects.txt 2> /dev/null")
        
if __name__ == '__main__':
    main()
