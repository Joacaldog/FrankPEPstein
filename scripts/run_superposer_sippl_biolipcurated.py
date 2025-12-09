import os
import sys

initial_path = os.getcwd()
minipockets_folder = sys.argv[1]
cmd1 = f'python3 ~/scripts/superposerV5.2_leave1out.py -T pocket.pdb -d /mnt/c/Users/Joacaldo/OneDrive - Universidad Católica de Chile/FrankPEPstein/DB/minipockets_surface80_winsize3_size3_curated-db -r 0.1 -t 36 -a 3 -fm {minipockets_folder}'

for folder in next(os.walk('.'))[1]:
    os.chdir(initial_path)
    os.chdir(folder)
    pdb = folder
    peptide_path = base_path = os.path.join(initial_path, folder, "peptide.pdb")
    receptor_path = base_path = os.path.join(initial_path, folder, "receptor.pdb")
    base_path = os.path.join(initial_path, folder, "superpockets_residuesAligned3_RMSD0.1/top_10_patches/")
    if not os.path.isdir(base_path):
        if not os.path.isfile(peptide_path):
            print(f"No peptide.pdb file in folder {folder}")
            continue
        if not os.path.isfile(receptor_path):
            print(f"No receptor.pdb file in folder {folder}")
            continue
        with open(f"{pdb}.gpf") as f:
            for line in f.readlines():
                if "npts" in line:
                    size = line.split()
                    x_size = size[1]
                    y_size = size[2]
                    z_size = size[3]
                if "gridcenter" in line:
                    gridcenter = line.split()
                    x_center = gridcenter[1]
                    y_center = gridcenter[2]
                    z_center = gridcenter[3]
            superposer_cmd = f'{cmd1} -x_center {x_center} -y_center {y_center} -z_center {z_center} -x_size {x_size} -y_size {y_size} -z_size {z_size}'
            print("Running with complex -->", folder)
            os.system(superposer_cmd)
            os.chdir("superpockets_residuesAligned3_RMSD0.1")
            os.system(f'cp ../receptor.pdb .')
            os.system(f'python3 ~/scripts/frankVINA_FNKPSTN.py receptor.pdb 36')
            os.system("rm * 2> /dev/null")
    if os.path.isdir(base_path):
        print(folder, "already processed")