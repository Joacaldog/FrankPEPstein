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
import shutil
import subprocess
import random  # Para hacer la muestra aleatoria

frank_folder_init = os.getcwd()
tqdm._instances.clear()

receptor_file_chain = sys.argv[1]
receptor_file = receptor_file_chain
threads = sys.argv[2]
selected_peps = int(sys.argv[3])

MAX_PEPTIDES = 100  # Límite de muestreo

def run_cmd(cmd):
    subprocess.run(cmd, shell=True, check=False,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def refine2(atmsel):
    md = MolecularDynamics(md_time_step=4, md_return='FINAL')
    md.optimize(atmsel, temperature=300, max_iterations=100)

def optimize2(atmsel, sched):
    for step in sched:
        step.optimize(atmsel, max_iterations=50)
    refine2(atmsel)
    cg = ConjugateGradients()
    cg.optimize(atmsel, max_iterations=100)

def minimization(pdb_file_no_extension, mode):
    log.level(output=0, notes=0, warnings=0, errors=0, memory=0)
    env = Environ()
    env.io.atom_files_directory = ['./']
    env.edat.dynamic_sphere = True
    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')

    code = pdb_file_no_extension
    if mode == "prot":
        mdl2 = complete_pdb(env, code, transfer_res_num=True)
        if len(mdl2.chains) > 1:
            pep = mdl2.chains[1]
            first_num_pep = pep.residues[0].num
            pep_chain_name = str(pep.name)
            pep_chain_length = pep.residues[-1].num
            pep_sel = selection(mdl2.residue_range(f'{first_num_pep}:{pep_chain_name}',
                                                   f'{pep_chain_length}:{pep_chain_name}'))
            atmsel = selection(pep_sel)
            mdl2.restraints.make(atmsel, restraint_type='stereo', spline_on_site=False)
            mdl2.env.edat.nonbonded_sel_atoms = 2
            sched = autosched.loop.make_for_model(mdl2)
            mdl2.env.edat.nonbonded_sel_atoms = 1
            optimize2(atmsel, sched)
            atmsel.energy()
            mdl2.write(file=code+'_min.pdb')
            run_cmd(f"rm {code}.pdb 2> /dev/null")

    if mode == "pep":
        mdl1 = complete_pdb(env, code, transfer_res_num=True)
        if len(mdl1.chains) == 1:
            pep = mdl1.chains[0]
            first_num_pep = pep.residues[0].num
            pep_chain_name = str(pep.name)
            pep_chain_length = pep.residues[-1].num
            pep_sel = selection(mdl1.residue_range(f'{first_num_pep}:{pep_chain_name}',
                                                   f'{pep_chain_length}:{pep_chain_name}'))
            atmsel = selection(pep_sel)
            mdl1.restraints.make(atmsel, restraint_type='stereo', spline_on_site=False)
            atmsel.energy()
            cg = ConjugateGradients()
            md = MolecularDynamics()
            cg.optimize(atmsel, max_iterations=60)
            md.optimize(atmsel, temperature=300, max_iterations=30)
            atmsel.energy()
            mdl1.write(file='min_'+code+'.pdb')
            run_cmd(f"rm {code}.pdb 2> /dev/null")

def scoring_filter():
    energy_patch_list = []
    for file in os.listdir("."):
        if fnmatch.fnmatch(file, '*.log'):
            with open(file) as f:
                for line in f:
                    line = line.strip()
                    if "Estimated" in line:
                        energy_value = line.split(":")[1].split(" ")[1]
                        if float(energy_value) < 0:
                            energy_patch_list.append((float(energy_value), file))

    sorted_list = sorted(energy_patch_list, key=itemgetter(0))
    if len(sorted_list) >= 1:
        folder_output3 = f"../top_{selected_peps}_peps"
        if not os.path.exists(folder_output3):
            os.makedirs(folder_output3)
        with open(f"{folder_output3}/top{selected_peps}_peps.tsv", "w") as outfile:
            outfile.write("AffinityBindingPred\tPEP\n")
            for selected in sorted_list[:selected_peps]:
                energy_value = selected[0]
                outfile.write(f"{energy_value}\t{selected[1].replace('.log', '')}\n")
                patch_file = selected[1].replace(".log", ".pdb")
                run_cmd(f"mv MinPEP_{patch_file.replace('.pdb', '')}_out.pdbqt {folder_output3} 2> /dev/null")

def main():
    # 1) Crear temp_folder
    if not os.path.exists("temp_folder"):
        os.makedirs("temp_folder")

    # reduce receptor
    run_cmd(f"reduce -Quiet -DB \"/mnt/c/Users/Joacaldo/OneDrive - Universidad Católica de Chile/FrankPEPstein/scripts/reduce_wwPDB_het_dict.txt\" {receptor_file} 1> temp_folder/H_{receptor_file} 2> /dev/null")

    os.chdir("temp_folder")
    run_cmd(f"sed -i '/END/d' H_{receptor_file} ; prepare_receptor -r H_{receptor_file} -o MinREC_{receptor_file}qt 1> /dev/null 2> /dev/null")
    os.chdir(frank_folder_init)

    if not os.path.exists("results_folder"):
        os.makedirs("results_folder")

    def vina_scorer(file):
        os.chdir(frank_folder_init)
        if "noEND" not in file:
            run_cmd(f"mv {file} temp_folder 2> /dev/null")
            os.chdir("temp_folder")
            minimization(file.replace(".pdb", ""), "pep")
            min_file = f"min_{file.replace('.pdb','')}.pdb"
            run_cmd(f"cat H_{receptor_file} {min_file} > complex_{min_file} 2> /dev/null")
            complex_file = f"complex_{min_file.replace('.pdb','')}"
            minimization(complex_file, "prot")
            complex_min_file = f'{complex_file}_min.pdb'
            run_cmd(
                f'cat {complex_min_file} | grep " x " | grep -v "TER" 1> MinPEP_{min_file} 2> /dev/null ; '
                f'reduce -Quiet -DB /mnt/c/Users/Joacaldo/OneDrive - Universidad Católica de Chile/FrankPEPstein/scripts/reduce_wwPDB_het_dict.txt MinPEP_{min_file} 1> H_MinPEP_{min_file} 2> /dev/null ; '
                f"sed -i '/END/d' H_MinPEP_{min_file} ; "
                f"prepare_ligand -l H_MinPEP_{min_file} -o MinPEP_{min_file.replace('.pdb', '.pdbqt')} 1> /dev/null 2> /dev/null"
            )
            log_file = f"{min_file.replace('.pdb','')}.log"
            cmd_vina = (
                f"vina --verbosity 0 --autobox --local_only "
                f"--receptor MinREC_{receptor_file}qt --ligand MinPEP_{min_file.replace('.pdb', '.pdbqt')} > {log_file}"
            )
            run_cmd(cmd_vina)
            out_pdbqt = f"MinPEP_{min_file.replace('.pdb','')}_out.pdbqt"
            run_cmd(f"mv {out_pdbqt} {log_file} ../results_folder")
            run_cmd(
                f"rm {complex_min_file} MinREC_{receptor_file} complex_{min_file} MinPEP_{min_file}qt "
                f"H_MinPEP_{min_file} MinREC_{min_file}qt {complex_file}_min.pdb MinPEP_{min_file} "
                f"{min_file} min_{file.replace('.pdb','')} 2> /dev/null"
            )

    # Filtrar los 'frag*.pdb'
    all_files = os.listdir(".")
    frag_files = [f for f in all_files if fnmatch.fnmatch(f, 'frag*.pdb')]

    # Si hay más de 1000, hacemos una muestra aleatoria
    if len(frag_files) > MAX_PEPTIDES:
        frag_files = random.sample(frag_files, MAX_PEPTIDES)
    
    # Resto de archivos se ignoran, solo procesamos la muestra

    # Paralelizar con joblib usando la lista reducida
    Parallel(n_jobs=int(threads))(
        delayed(vina_scorer)(file)
        for file in tqdm(frag_files, total=len(frag_files), desc="Minimizing complexes")
    )

    os.chdir("results_folder")
    scoring_filter()
    os.chdir(frank_folder_init)
    run_cmd("rm -r results_folder")
    run_cmd("rm -r temp_folder")

    # Convertir a .pdb
    current_dir = os.getcwd()
    top_dir = f"top_{selected_peps}_peps"
    full_path_top = os.path.join(current_dir, top_dir)
    if os.path.exists(full_path_top):
        os.chdir(top_dir)
        for pep_pdbqt in os.listdir("."):
            if fnmatch.fnmatch(pep_pdbqt, '*.pdbqt'):
                base = pep_pdbqt.replace(".pdbqt", "")
                run_cmd(f"obabel -ipdbqt {pep_pdbqt} -o pdb -O {base}.pdb ; rm {pep_pdbqt}")
    else:
        print(f"No final candidates found in {current_dir}.")
def main_wrapper():
    main()

if __name__ == '__main__':
    main_wrapper()
