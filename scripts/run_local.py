import os
import sys
import argparse
import shutil
import subprocess
import multiprocessing
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Chain import Chain

def merge_to_chain_p(model):
    """Collapse every chain of a model into a single chain 'p' with unique residue numbers.

    Renaming each chain to 'p' independently produced several blocks all labelled 'p' with
    colliding residue numbers; Biopython keeps only the first residue of each duplicate
    (chain, resseq) when the file is read back, so on a multi-chain receptor roughly half the
    pocket silently disappeared. The single chain id 'p' is load-bearing further down the
    pipeline (superposer.py selects the pocket with `grep ' p'`), so the ids are merged rather
    than diversified.
    """
    chains = list(model)
    if len(chains) == 1:
        chains[0].id = 'p'
        return chains[0]
    merged = Chain('p')
    index = 0
    for chain in chains:
        for residue in list(chain):
            chain.detach_child(residue.id)
            index += 1
            residue.id = (residue.id[0], index, ' ')
            merged.add(residue)
        model.detach_child(chain.id)
    model.add(merged)
    return merged


def calculate_box_and_save(input_pocket, output_pocket):
    """
    Reads input pocket, merges all chains into chain 'p', saves to output_pocket,
    and calculates the gridbox as the tight bounding box of the pocket atoms.
    """
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("pocket", input_pocket)

        atoms = []
        for model in structure:
            merged = merge_to_chain_p(model)
            for residue in merged:
                for atom in residue:
                    atoms.append(atom)

        if not atoms:
            print("Error: No atoms found in pocket.")
            return None, None
            
        # Save processed pocket
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_pocket)
        
        # Calculate Box
        coords = [a.get_coord() for a in atoms]
        min_coord = [min([c[i] for c in coords]) for i in range(3)]
        max_coord = [max([c[i] for c in coords]) for i in range(3)]
        
        # Deliberately no padding: the box is the tight bounding box of the pocket atoms.
        # In the Colab UI the Size controls let the user add clearance interactively.
        buffer = 0.0
        
        center = [(min_coord[i] + max_coord[i]) / 2 for i in range(3)]
        size = [(max_coord[i] - min_coord[i]) + buffer for i in range(3)]
        
        return center, size
        
    except Exception as e:
        print(f"Error processing pocket: {e}")
        return None, None

def main():
    parser = argparse.ArgumentParser(description="Run FrankPEPstein Locally")
    parser.add_argument("-r", "--receptor", required=True, help="Path to Receptor PDB")
    parser.add_argument("-p", "--pocket", required=True, help="Path to Pocket PDB")
    parser.add_argument("-w", "--pep_size", type=int, default=8, help="Peptide size (kmer)")
    parser.add_argument("-t", "--threads", type=int, default=multiprocessing.cpu_count(), help="Number of threads (default: all cpus)")
    parser.add_argument("-c", "--candidates", type=int, default=10, help="Number of candidates")
    parser.add_argument("-s", "--sampling", type=int, default=500, help="Subsampling limit for Vina screening")
    parser.add_argument("-rmsd", "--rmsd_allowed", type=float, default=0.5, help="RMSD cutoff for superposer")
    parser.add_argument("-seed", "--seed", type=int, default=0, help="Random seed for the sampling steps, for reproducible runs")

    args = parser.parse_args()
    
    # Setup Paths
    current_dir = os.getcwd()
    # Derive the repo root from this file, not from cwd: setup_local.py requires cwd to be
    # the repo root, while the Colab layout has the repo as a subdirectory of cwd. Deriving
    # from __file__ satisfies both.
    repo_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    scripts_folder = os.path.join(repo_folder, "scripts")
    run_script = os.path.join(scripts_folder, "run_FrankPEPstein.py")
    
    # 1. Prepare Receptor
    dest_receptor = os.path.join(current_dir, "receptor.pdb")
    if os.path.abspath(args.receptor) != dest_receptor:
        print(f"Copying receptor to {dest_receptor}...")
        shutil.copy(args.receptor, dest_receptor)
        
    # 2. Prepare Pocket
    # User Request: "dejar el pocket.pdb inicial junto a receptor.pdb en el main directory"
    pocket_basename = os.path.basename(args.pocket)
    dest_initial_pocket = os.path.join(current_dir, pocket_basename)
    
    if os.path.abspath(args.pocket) != dest_initial_pocket:
        print(f"Copying initial pocket to {dest_initial_pocket}...")
        shutil.copy(args.pocket, dest_initial_pocket)
        
    # Process Pocket -> FrankPEPstein_run/pocket.pdb
    pockets_dir = os.path.join(current_dir, "FrankPEPstein_run")
    os.makedirs(pockets_dir, exist_ok=True)
    dest_processed_pocket = os.path.join(pockets_dir, "pocket.pdb")
    dest_receptor_run = os.path.join(pockets_dir, "receptor.pdb") # [ADDED] Save receptor to run folder
    shutil.copy(dest_receptor, dest_receptor_run)
    
    print(f"Processing pocket to {dest_processed_pocket}...")
    center, size = calculate_box_and_save(dest_initial_pocket, dest_processed_pocket)
    
    # User Requirement: The pocket in the root should also be 'p'
    # Overwrite the initial copy with the processed one
    if center and size:
        print(f"Updating {dest_initial_pocket} with processed chain 'p'...")
        shutil.copy(dest_processed_pocket, dest_initial_pocket)
    
    if not center or not size:
        print("Failed to process pocket.")
        sys.exit(1)
        
    print(f"Gridbox Calculated:")
    print(f"  Center: {center}")
    print(f"  Size:   {size}")
    
    # 3. Execution
    cmd = [
        sys.executable, run_script,
        "-w", str(args.pep_size),
        "-t", str(args.threads),
        "-c", str(args.candidates),
        "-s", str(args.sampling),
        "-rmsd", str(args.rmsd_allowed),
        "-seed", str(args.seed),
        "-xc", str(center[0]),
        "-yc", str(center[1]),
        "-zc", str(center[2]),
        "-xs", str(size[0]),
        "-ys", str(size[1]),
        "-zs", str(size[2])
    ]
    
    print("-" * 40)
    print(f"Running Pipeline: {' '.join(cmd)}")
    print("-" * 40)
    
    try:
        subprocess.run(cmd, check=True)
        print("✅ Local run completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"❌ Execution failed with code {e.returncode}")
        sys.exit(e.returncode)

if __name__ == "__main__":
    main()
