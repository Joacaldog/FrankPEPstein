import os
import warnings
from Bio import PDB

initial_path = os.getcwd()

def count_residues_in_chain(pdb_file_path):
    warnings.simplefilter('ignore')
    structure = PDB.PDBParser().get_structure("pdb_structure", pdb_file_path)
    
    # Initialize a set to store unique residue numbers
    unique_residues = set()
    
    for model in structure:
        for chain in model:
            for residue in chain:
                # Only consider standard amino acid residues (not water, ions, etc.)
                if PDB.is_aa(residue):
                    unique_residues.add(residue.get_id()[1])
    
    return len(unique_residues)

for folder in next(os.walk('.'))[1]:
    os.chdir(initial_path)
    os.chdir(folder)
    pdb = folder.split("_")[0]
    peptide_path = base_path = os.path.join(initial_path, folder, "peptide.pdb")
    receptor_path = base_path = os.path.join(initial_path, folder, "receptor.pdb")
    if not os.path.isfile(peptide_path):
        print(f"No peptide.pdb file in folder {folder}")
        continue
    if not os.path.isfile(receptor_path):
        print(f"No receptor.pdb file in folder {folder}")
        continue

    pep_size = count_residues_in_chain(f"peptide.pdb")
    base_path = os.path.join(initial_path, folder, "superpockets_residuesAligned3_RMSD0.1/top_10_patches")
    pep_tsv = f"{base_path}/frankPEPstein_{pep_size}/top_1_peps/top1_peps.tsv"
    pep_dir = f"{base_path}/frankPEPstein_{pep_size}/top_1_peps/"
    scripts_dir = os.path.dirname(os.path.realpath(__file__))

    if os.path.isdir(base_path):
        if not os.path.exists(pep_tsv):
            print("Running complex --->", folder)
            os.chdir("superpockets_residuesAligned3_RMSD0.1/top_10_patches/")
            if len([x for x in os.listdir(".") if "patch_file" in x]) == 0:
                print("No patches files in folder")
            if len([x for x in os.listdir(".") if "patch_file" in x]) > 1:
                print(f"Running patch_clustering with kmer: {pep_size} ")
                os.system(f'python3 "{scripts_dir}/patch_clustering_V8.7.py" -w {pep_size} -t 36')
                os.chdir(f"frankPEPstein_{pep_size}")
                os.system(f'cp ../../../receptor.pdb .')
                # print(f'python3 "{scripts_dir}/frankVINA_V3.py" receptor.pdb 36 1')
                os.system(f'python3 "{scripts_dir}/frankVINA_V3.py" receptor.pdb 36 1')
                os.system("rm * 2> /dev/null")
            if len([x for x in os.listdir(".") if "patch_file" in x]) == 1:
                print("Only one patch file in folder")
                os.makedirs(pep_dir, exist_ok=True)
                os.system(f"cp patch_file*.pdb {pep_dir}frag_1.pdb")
        if os.path.exists(pep_tsv):
            print(folder, "already processed")
    if not os.path.isdir(base_path):
        print(folder, "no top_10_patches")
