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
from collections import Counter

initial_path = sys.argv[1]
# The per-peptide MODELLER refinement and the screening dock run against the pocket shell only,
# which is what makes this stage affordable. That also means the screening score ignores the rest
# of the protein, so the survivors are rescored against the full receptor at the end (see
# rescore_against_full_receptor).
receptor_file = "pocket.pdb"
full_receptor_file = "receptor.pdb"
threads = sys.argv[2]
selected_peps = int(sys.argv[3])
try:
    MAX_PEPTIDES = int(sys.argv[4])
except IndexError:
    MAX_PEPTIDES = 500 # Default if not provided
try:
    SEED = int(sys.argv[5])
except IndexError:
    SEED = 0
random.seed(SEED)
ADFR_DIR = os.path.join(initial_path, "utilities/ADFRsuite_x86_64Linux_1.0")
ADFR_BIN = os.path.join(ADFR_DIR, "bin")
ADFR_LIB = os.path.join(ADFR_DIR, "lib")

# Prepend to PATH so tools can find their python2.7 and libraries
os.environ["PATH"] = ADFR_BIN + os.pathsep + os.environ["PATH"]
# Set LD_LIBRARY_PATH so binaries find libpython2.7.so.1.0
os.environ["LD_LIBRARY_PATH"] = ADFR_LIB + os.pathsep + os.environ.get("LD_LIBRARY_PATH", "")

# Configuration Variables
REDUCE_PATH = os.path.join(ADFR_BIN, "reduce")
REDUCE_DB_PATH = f"{initial_path}/DB/reduce_wwPDB_het_dict.txt"
PREPARE_RECEPTOR_PATH = os.path.join(ADFR_BIN, "prepare_receptor")
PREPARE_LIGAND_PATH = os.path.join(ADFR_BIN, "prepare_ligand")
VINA_PATH = f"{initial_path}/utilities/vina_1.2.4_linux_x86_64"
OBABEL_PATH = os.path.join(ADFR_BIN, "obabel")

frank_folder_init = os.getcwd()
tqdm._instances.clear()
HAS_FULL_RECEPTOR = False

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
        if len(mdl2.chains) <= 1:
            # No peptide chain to refine. Returning silently here used to leave the expected
            # _min.pdb unwritten, and every consumer downstream is a shell command with stderr
            # discarded, so the peptide just vanished from the results with no explanation.
            return None, f"prot: expected more than 1 chain, got {len(mdl2.chains)}"
        if True:
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
            return code + '_min.pdb', None

    if mode == "pep":
        mdl1 = complete_pdb(env, code, transfer_res_num=True)
        if len(mdl1.chains) != 1:
            return None, f"pep: expected exactly 1 chain, got {len(mdl1.chains)}"
        if True:
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
            return 'min_' + code + '.pdb', None

    return None, f"unknown minimization mode {mode!r}"

def scoring_filter(n_submitted=None):
    energy_patch_list = []
    n_rejected = 0
    n_failed = 0
    for file in os.listdir("."):
        if fnmatch.fnmatch(file, '*.log'):
            scored = False
            with open(file) as f:
                for line in f:
                    line = line.strip()
                    if "Estimated" in line:
                        energy_value = line.split(":")[1].split(" ")[1]
                        scored = True
                        if float(energy_value) < 0:
                            energy_patch_list.append((float(energy_value), file))
                        else:
                            n_rejected += 1
                        # One score per log: a log with two Estimated lines used to occupy two
                        # slots in the top-N.
                        break
            if not scored:
                n_failed += 1

    # A crashed Vina run and a peptide that simply did not bind used to be indistinguishable:
    # both just never appeared. Report the breakdown so an empty result set is diagnosable.
    n_scored = len(energy_patch_list)
    summary = (f"Scored: {n_scored} | Non-binding (>=0): {n_rejected} | "
               f"Failed (log without an Estimated line): {n_failed}")
    if n_submitted is not None:
        summary += f" | No log at all: {n_submitted - n_scored - n_rejected - n_failed} | Submitted: {n_submitted}"
    print(summary)

    sorted_list = sorted(energy_patch_list, key=itemgetter(0))
    if len(sorted_list) >= 1:
        folder_output3 = f"../top_{selected_peps}_peps"
        if not os.path.exists(folder_output3):
            os.makedirs(folder_output3)
        with open(f"{folder_output3}/top{selected_peps}_peps.tsv", "w") as outfile:
            outfile.write("AffinityBindingPred\tPEP\n")
            written = 0
            missing = 0
            # Check the pose exists BEFORE writing its row, and keep going down the ranking
            # instead of leaving a short table. Previously the row was written first and the
            # pose moved with `mv ... 2> /dev/null`, so a missing pose produced a TSV entry
            # with no structure behind it and silently consumed one of the N slots.
            for energy_value, log_name in sorted_list:
                if written >= selected_peps:
                    break
                base = log_name.replace(".log", "")
                pose = f"MinPEP_{base}_out.pdbqt"
                if not os.path.exists(pose):
                    missing += 1
                    continue
                outfile.write(f"{energy_value}\t{base}\n")
                shutil.move(pose, os.path.join(folder_output3, pose))
                written += 1
        if missing:
            print(f"Note: skipped {missing} ranked candidate(s) whose pose file was missing.")
        if written < selected_peps:
            print(f"Note: only {written} of the {selected_peps} requested candidates had both a "
                  f"score and a pose.")

def rescore_against_full_receptor(top_dir):
    """Rescore the selected poses against the whole protein, not just the pocket shell.

    --score_only, not --local_only: the pose was already optimised inside the pocket; what this
    measures is whether it also clashes with the parts of the receptor the screening dock never
    saw. A pose that looks good against the shell and bad here is sitting in the protein body.
    """
    receptor_qt = os.path.join(frank_folder_init, "temp_folder", f"MinREC_{full_receptor_file}qt")
    tsv = os.path.join(top_dir, f"top{selected_peps}_peps.tsv")
    if not os.path.exists(tsv):
        return
    with open(tsv) as handle:
        rows = [line.rstrip("\n").split("\t") for line in handle if line.strip()]
    if len(rows) < 2:
        return

    rescored = []
    for energy_value, base in rows[1:]:
        pose = os.path.join(top_dir, f"MinPEP_{base}_out.pdbqt")
        score = ""
        if os.path.exists(pose):
            result = subprocess.run(
                [VINA_PATH, "--verbosity", "0", "--score_only",
                 "--receptor", receptor_qt, "--ligand", pose],
                check=False, capture_output=True, text=True)
            for line in result.stdout.splitlines():
                if "Estimated" in line:
                    score = line.split(":")[1].split()[0]
                    break
        rescored.append((energy_value, base, score))

    with open(tsv, "w") as handle:
        handle.write("AffinityBindingPred_pocket\tPEP\tAffinityBindingPred_fullReceptor\n")
        for energy_value, base, score in rescored:
            handle.write(f"{energy_value}\t{base}\t{score}\n")
    print(f"Rescored {sum(1 for r in rescored if r[2])} of {len(rescored)} candidates "
          f"against the full receptor.")


def main():
    # 1) Crear temp_folder
    if not os.path.exists("temp_folder"):
        os.makedirs("temp_folder")

    # reduce receptor
    run_cmd(f"{REDUCE_PATH} -Quiet -DB {REDUCE_DB_PATH} {receptor_file} 1> temp_folder/H_{receptor_file} 2> /dev/null")

    os.chdir("temp_folder")
    run_cmd(f"sed -i '/END/d' H_{receptor_file} ; {PREPARE_RECEPTOR_PATH} -r H_{receptor_file} -o MinREC_{receptor_file}qt 1> /dev/null 2> /dev/null")
    os.chdir(frank_folder_init)

    # Prepare the full receptor for the final rescore. Not fatal if absent: the pocket-only
    # ranking still works, but say so rather than quietly reporting shell-only affinities.
    global HAS_FULL_RECEPTOR
    HAS_FULL_RECEPTOR = False
    if os.path.exists(full_receptor_file):
        run_cmd(f"{REDUCE_PATH} -Quiet -DB {REDUCE_DB_PATH} {full_receptor_file} "
                f"1> temp_folder/H_{full_receptor_file} 2> /dev/null")
        os.chdir("temp_folder")
        run_cmd(f"sed -i '/END/d' H_{full_receptor_file}")
        prep = subprocess.run(
            [PREPARE_RECEPTOR_PATH, "-r", f"H_{full_receptor_file}",
             "-o", f"MinREC_{full_receptor_file}qt"],
            check=False, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
        os.chdir(frank_folder_init)
        HAS_FULL_RECEPTOR = os.path.exists(
            os.path.join("temp_folder", f"MinREC_{full_receptor_file}qt"))
        if not HAS_FULL_RECEPTOR:
            print(f"WARNING: could not prepare the full receptor for rescoring. "
                  f"prepare_receptor said: {(prep.stderr or '').strip()[:300]}")
    else:
        print(f"WARNING: {full_receptor_file} not found next to the fragments. Reported affinities "
              f"will come from the pocket shell alone and ignore the rest of the protein.")

    if not os.path.exists("results_folder"):
        os.makedirs("results_folder")

    def vina_scorer(file):
        # joblib's SequentialBackend (n_jobs=1, e.g. -t 1 or a 1-vCPU host) runs this in the calling
        # process, so the chdir below would leak into main() and break every path after the loop.
        # Exceptions are caught per peptide: a raise here propagates through joblib and cancels
        # every remaining task, turning "lose 1 of 500 peptides" into "lose the whole screen".
        try:
            return _vina_scorer(file)
        except Exception as exc:
            return (file, f"exception: {exc}")
        finally:
            os.chdir(frank_folder_init)

    def _vina_scorer(file):
        os.chdir(frank_folder_init)
        if "noEND" not in file:
            run_cmd(f"mv {file} temp_folder 2> /dev/null")
            os.chdir("temp_folder")
            _, why = minimization(file.replace(".pdb", ""), "pep")
            min_file = f"min_{file.replace('.pdb','')}.pdb"
            if why or not os.path.exists(min_file):
                return (file, why or "peptide minimization produced no output")
            run_cmd(f"cat H_{receptor_file} {min_file} > complex_{min_file} 2> /dev/null")
            complex_file = f"complex_{min_file.replace('.pdb','')}"
            _, why = minimization(complex_file, "prot")
            complex_min_file = f'{complex_file}_min.pdb'
            if why or not os.path.exists(complex_min_file):
                return (file, why or "complex minimization produced no output")
            run_cmd(
                f'cat {complex_min_file} | grep " x " | grep -v "TER" 1> MinPEP_{min_file} 2> /dev/null ; '
                f'{REDUCE_PATH} -Quiet -DB {REDUCE_DB_PATH} MinPEP_{min_file} 1> H_MinPEP_{min_file} 2> /dev/null ; '
                f"sed -i '/END/d' H_MinPEP_{min_file} ; "
                f"{PREPARE_LIGAND_PATH} -l H_MinPEP_{min_file} -o MinPEP_{min_file.replace('.pdb', '.pdbqt')} 1> /dev/null 2> /dev/null"
            )
            log_file = f"{min_file.replace('.pdb','')}.log"
            cmd_vina = (
                f"{VINA_PATH} --verbosity 0 --autobox --local_only "
                f"--receptor MinREC_{receptor_file}qt --ligand MinPEP_{min_file.replace('.pdb', '.pdbqt')} > {log_file}"
            )
            run_cmd(cmd_vina)
            out_pdbqt = f"MinPEP_{min_file.replace('.pdb','')}_out.pdbqt"
            if not os.path.exists(out_pdbqt):
                return (file, "vina produced no output pose")
            run_cmd(f"mv {out_pdbqt} {log_file} ../results_folder")
            # NOTE: `MinREC_{receptor_file}` below names MinREC_pocket.pdb, which is never
            # created, so it is a harmless no-op. Do NOT "correct" it to MinREC_pocket.pdbqt:
            # that IS the shared prepared receptor, built once and used by every worker, and
            # the first peptide to finish would delete it for all the others.
            run_cmd(
                f"rm {complex_min_file} MinREC_{receptor_file} complex_{min_file} MinPEP_{min_file}qt "
                f"H_MinPEP_{min_file} MinREC_{min_file}qt {complex_file}_min.pdb MinPEP_{min_file} "
                f"{min_file} min_{file.replace('.pdb','')} 2> /dev/null"
            )

    # Filtrar los 'frag*.pdb'
    all_files = os.listdir(".")
    # Sorted so that the seeded sample below does not depend on os.listdir order.
    frag_files = sorted(f for f in all_files if fnmatch.fnmatch(f, 'frag*.pdb'))

    total_combinations = len(frag_files)
    print(f"Total Combinations Found: {total_combinations}")

    # Si hay más de MAX_PEPTIDES, hacemos una muestra aleatoria
    if len(frag_files) > MAX_PEPTIDES:
        print(f"Subsampling to {MAX_PEPTIDES} for processing speed...")
        frag_files = random.sample(frag_files, MAX_PEPTIDES)
    
    # Resto de archivos se ignoran, solo procesamos la muestra

    # Paralelizar con joblib usando la lista reducida
    statuses = Parallel(n_jobs=int(threads))(
        delayed(vina_scorer)(file)
        for file in tqdm(frag_files, total=len(frag_files), desc="Minimizing complexes")
    )

    skipped = [st for st in statuses if st]
    if skipped:
        print(f"{len(skipped)} of {len(frag_files)} peptides were dropped before scoring:")
        reasons = Counter(reason for _, reason in skipped)
        for reason, count in reasons.most_common():
            print(f"  {count:>5}  {reason}")

    print("Minimization complete. Ranking results...")
    os.chdir(os.path.join(frank_folder_init, "results_folder"))
    scoring_filter(len(frag_files))
    os.chdir(frank_folder_init)

    # Must run before temp_folder is removed: the prepared full receptor lives there.
    top_dir_path = os.path.join(frank_folder_init, f"top_{selected_peps}_peps")
    if HAS_FULL_RECEPTOR and os.path.isdir(top_dir_path):
        print("Rescoring top candidates against the full receptor...")
        rescore_against_full_receptor(top_dir_path)

    print("Cleaning up temporary files...")
    run_cmd("rm -r results_folder")
    run_cmd("rm -r temp_folder")

    # Convertir a .pdb
    print("Converting top candidates to PDB...")
    current_dir = os.getcwd()
    top_dir = f"top_{selected_peps}_peps"
    full_path_top = os.path.join(current_dir, top_dir)
    if os.path.exists(full_path_top):
        os.chdir(top_dir)
        for pep_pdbqt in os.listdir("."):
            if fnmatch.fnmatch(pep_pdbqt, '*.pdbqt'):
                base = pep_pdbqt.replace(".pdbqt", "")
                run_cmd(f"{OBABEL_PATH} -ipdbqt {pep_pdbqt} -o pdb -O {base}.pdb ; rm {pep_pdbqt}")
        print(f"Done! Results in {top_dir}")
    else:
        print(f"No final candidates found in {current_dir}.")
def main_wrapper():
    main()

if __name__ == '__main__':
    main_wrapper()
