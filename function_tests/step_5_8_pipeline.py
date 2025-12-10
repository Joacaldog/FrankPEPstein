#@title 3. FrankPEPstein Pipeline
#@markdown **Instructions:**
#@markdown 1. Configure the pipeline parameters (Length, Number of Peptides).
#@markdown 2. Click **Run Pipeline** to start the search and generation process.

import os
import sys
import subprocess
import shutil
import ipywidgets as widgets
from IPython.display import display

# --- Configuration Widgets ---
style = {'description_width': 'initial'}

length_slider = widgets.IntSlider(
    value=10,
    min=3,
    max=30,
    step=1,
    description='Target Peptide Length:',
    style=style,
    layout=widgets.Layout(width='50%')
)

num_peptides_slider = widgets.IntSlider(
    value=10,
    min=1,
    max=100,
    step=1,
    description='Number of Output Peptides:',
    style=style,
    layout=widgets.Layout(width='50%')
)

threads_slider = widgets.IntSlider(
    value=2,
    min=1,
    max=8,
    step=1,
    description='CPU Threads:',
    style=style,
    layout=widgets.Layout(width='50%')
)

# Output area
run_output = widgets.Output()

# --- Execution Logic ---
import json

def run_frankpepstein_pipeline(b):
    run_output.clear_output()
    with run_output:
        # 0. Load State (Persistence Check)
        global box_center, box_size, extracted_pocket_path, receptor_filename
        
        state_file = "pipeline_state.json"
        if os.path.exists(state_file):
            try:
                with open(state_file, "r") as f:
                    state = json.load(f)
                    
                # Restore globals if missing
                if 'box_center' not in globals() and 'box_center' in state:
                    box_center = state['box_center']
                    print("✅ Restored box_center from state file.")
                if 'box_size' not in globals() and 'box_size' in state:
                    box_size = state['box_size']
                    print("✅ Restored box_size from state file.")
                if 'extracted_pocket_path' not in globals() and 'extracted_pocket_path' in state:
                    extracted_pocket_path = state['extracted_pocket_path']
                    print("✅ Restored extracted_pocket_path from state file.")
                if 'receptor_filename' not in globals() and 'receptor_filename' in state:
                    receptor_filename = state['receptor_filename']
                    print("✅ Restored receptor_filename from state file.")
            except Exception as e:
                print(f"⚠️ Warning: Could not load state file: {e}")

        # 1. Validate Prereqs
        if 'box_center' not in globals() or 'box_size' not in globals():
             print("❌ Error: Box parameters (center/size) not defined. Please run Step 4 first.")
             return
        if 'receptor_filename' not in globals() or not receptor_filename:
             print("❌ Error: Receptor filename not defined. Please run Step 1 first.")
             return
             
        # Check environment
        frank_python = "/usr/local/envs/FrankPEPstein/bin/python"
        # Paths to scripts (Assuming repo structure)
        # Force absolute base path
        if os.path.exists("/content/FrankPEPstein"):
            base_dir = "/content/FrankPEPstein"
        else:
            base_dir = os.path.abspath(os.getcwd())
            
        print(f"Base Directory: {base_dir}")

        scripts_dir = os.path.join(base_dir, "scripts")
        db_path = os.path.join(base_dir, "DB", "minipockets_surface80_winsize3_size3_curated-db")
        
        # 2. Setup Run Directory (Absolute Path)
        run_dir = os.path.join(base_dir, "FrankPEPstein_Run")
        
        if os.path.exists(run_dir):
            shutil.rmtree(run_dir) # Clean start
        os.makedirs(run_dir)
        
        # Copy Receptor (Full)
        target_receptor_full = "receptor.pdb"
        try:
            # Source file check
            if os.path.isabs(receptor_filename) and os.path.exists(receptor_filename):
                src_receptor = receptor_filename
            elif os.path.exists(os.path.join(base_dir, receptor_filename)):
                src_receptor = os.path.join(base_dir, receptor_filename)
            elif os.path.exists(os.path.join("/content", receptor_filename)):
                src_receptor = os.path.join("/content", receptor_filename)
            else:
                src_receptor = os.path.abspath(receptor_filename)

            if not os.path.exists(src_receptor):
                 print(f"❌ Error: Could not find receptor file '{receptor_filename}'.")
                 return
                 
            shutil.copy(src_receptor, os.path.join(run_dir, target_receptor_full))
        except FileNotFoundError:
             print(f"❌ Error copying receptor.")
             return

        # Prepare Pocket for Superposer (Chain 'p')
        target_pocket_file = "target_pocket.pdb"
        if 'extracted_pocket_path' not in globals() or not extracted_pocket_path:
            print("❌ Error: Pocket path not defined. Please run extraction (Step 3/4) first.")
            return

        # Prepare Pocket for Superposer (Chain 'p')
        target_pocket_file = "target_pocket.pdb"
        if 'extracted_pocket_path' not in globals() or not extracted_pocket_path:
            print("❌ Error: Pocket path not defined. Please run extraction (Step 3/4) first.")
            return

        print(f"Preparing pocket from: {extracted_pocket_path}")
        
        # Isolate BioPython usage to avoid Colab kernel issues
        prep_script_content = f"""
import sys
from Bio import PDB
import os

def prepare_pocket(input_path, output_path):
    try:
        parser = PDB.PDBParser(QUIET=True)
        struct = parser.get_structure("pocket", input_path)
        
        # Rename all chains to 'p'
        for model in struct:
            for chain in model:
                chain.id = 'p'
        
        io = PDB.PDBIO()
        io.set_structure(struct)
        io.save(output_path)
        print("SUCCESS")
    except Exception as e:
        print(f"ERROR: {{e}}")

if __name__ == "__main__":
    prepare_pocket('{extracted_pocket_path}', '{os.path.join(run_dir, target_pocket_file)}')
"""
        prep_script_path = os.path.join(run_dir, "prep_pocket.py")
        with open(prep_script_path, "w") as f:
            f.write(prep_script_content)
            
        try:
            res = subprocess.run([frank_python, prep_script_path], capture_output=True, text=True)
            if "SUCCESS" in res.stdout:
                print(f"✅ Created {target_pocket_file} with chain 'p' for Superposer.")
            else:
                 print(f"❌ Error preparing pocket PDB: {res.stdout} {res.stderr}")
                 return
        except Exception as e:
            print(f"❌ Error executing isolated prep script: {e}")
            return
            

        # 3. Define Parameters
        pep_length = length_slider.value
        n_peps = num_peptides_slider.value
        n_threads = threads_slider.value
        
        if not os.path.exists(db_path):
             print(f"❌ Error: Database not found at {db_path}")
             return

        print("\n--- Pipeline Configuration ---")
        print(f"Target Length: {pep_length}")
        print(f"Output Count : {n_peps}")
        print(f"Threads      : {n_threads}")

        # --- A. SUPERPOSER (Fragment Scanning) ---
        print("\n🚀 Starting Step 1: Fragment Scanning (Superposer)...")
        
        superposer_script = os.path.join(scripts_dir, "superposerV5.2_leave1out.py")
        
        # Uses target_pocket_file (chain p) as -T
        cmd_superposer = [
            frank_python, superposer_script,
            "-T", target_pocket_file, 
            "-d", db_path,
            "-a", "3", 
            "-r", "0.1",
            "-x_center", str(box_center[0]),
            "-y_center", str(box_center[1]),
            "-z_center", str(box_center[2]),
            "-x_size", str(box_size[0]),
            "-y_size", str(box_size[1]),
            "-z_size", str(box_size[2]),
            "-t", str(n_threads),
            "-fm", db_path
        ]
        
        try:
            # We run this INSIDE the run_dir to keep outputs contained
            subprocess.run(cmd_superposer, cwd=run_dir, check=True)
            print("✅ Superposer finished.")
        except subprocess.CalledProcessError as e:
            print(f"❌ Error running Superposer: {e}")
            return

        # --- B. FrankVINA Part 1 (Filtering Patches) ---
        print("\n🚀 Starting Step 2: Patch Filtering (FrankVINA I)...")
        
        # Output of superposer is in 'superpockets_residuesAligned3_RMSD0.1' inside run_dir
        super_out_dir = os.path.join(run_dir, "superpockets_residuesAligned3_RMSD0.1")
        
        if not os.path.exists(super_out_dir):
            print("❌ Error: Superposer output directory not found. Did it find any fragments?")
            return
            
        vina_script_1 = os.path.join(scripts_dir, "frankVINA_FNKPSTN.py")
        
        # Need to copy receptor there? run_superposer did it.
        # superposerV5.2 seems to put files in subfolders inside the output dir?
        # Wait, superposerV5.2 iterates over minipockets and puts output in `superpockets...`
        # Let's inspect what frankVINA_FNKPSTN expects. It expects to run IN the folder with patches.
        # But superposer puts all matches in that folder?
        # `os.system(f"cp \"{folder_file}\" .")` -> it copies minipockets to CWD (which is run_dir).
        # Actually superposerV5.2 creates `superpockets...` and puts outputs there?
        # "folder_output = ...; if not os.path.exists...makedirs"
        
        # The logic in `run_superposer...py` was:
        # os.chdir("superpockets_residuesAligned3_RMSD0.1")
        # cp ../receptor.pdb .
        # run frankVINA_FNKPSTN.py receptor.pdb threads
        
        # So we replicate that:
        try:
            # Copy receptor into the superpockets folder
            shutil.copy(os.path.join(run_dir, target_receptor), os.path.join(super_out_dir, target_receptor))
            
            cmd_vina1 = [
                frank_python, vina_script_1,
                target_receptor,
                str(n_threads)
            ]
            
            subprocess.run(cmd_vina1, cwd=super_out_dir, check=True)
            print("✅ Patch Filtering finished.")
            
        except Exception as e:
            print(f"❌ Error in Vina Step 1: {e}")
            return
            
        # --- C. Patch Clustering (Assembly) ---
        print(f"\n🚀 Starting Step 3: Peptide Assembly (PatchClustering) for {pep_length}-mers...")
        
        # Input for this step is in `superpockets.../top_10_patches`
        # Created by frankVINA_FNKPSTN
        patches_dir = os.path.join(super_out_dir, "top_10_patches")
        
        if not os.path.exists(patches_dir):
             print("❌ Error: 'top_10_patches' folder not found. No patches passed filtering?")
             return
             
        # Script expects to run inside that folder
        patch_clust_script = os.path.join(scripts_dir, "patch_clustering_V8.7.py")
        
        cmd_clust = [
            frank_python, patch_clust_script,
            "-w", str(pep_length),
            "-t", str(n_threads)
        ]
        
        try:
            subprocess.run(cmd_clust, cwd=patches_dir, check=True)
            print("✅ Patch Clustering finished.")
        except subprocess.CalledProcessError as e:
             print(f"❌ Error in Patch Clustering: {e}")
             return

        # --- D. FrankVINA Part 2 (Final Scoring) ---
        print("\n🚀 Starting Step 4: Final Refinement & Scoring (FrankVINA II)...")
        
        # Output of clustering is in `frankPEPstein_{winsize}` inside patches_dir
        final_dir = os.path.join(patches_dir, f"frankPEPstein_{pep_length}")
        
        if not os.path.exists(final_dir):
             print(f"❌ Error: Output directory '{final_dir}' not found. No peptides assembled?")
             return
             
        vina_script_2 = os.path.join(scripts_dir, "frankVINA_V3.py")
        
        try:
            # Copy receptor to final dir
            shutil.copy(os.path.join(super_out_dir, target_receptor), os.path.join(final_dir, target_receptor))
            
            cmd_vina2 = [
                frank_python, vina_script_2,
                target_receptor,
                str(n_threads),
                str(n_peps) # Number of top peptides to keep
            ]
            
            subprocess.run(cmd_vina2, cwd=final_dir, check=True)
            print("✅ Final Scoring finished.")
            
            # --- Results ---
            results_tsv = os.path.join(final_dir, f"top_{n_peps}_peps", f"top{n_peps}_peps.tsv")
            if os.path.exists(results_tsv):
                print(f"\n🎉 Success! Top peptides saved in: {results_tsv}")
                # Optional: specific display code for results
            else:
                print("⚠️ Warning: Pipeline finished but results TSV not found.")
                
        except Exception as e:
            print(f"❌ Error in Vina Step 2: {e}")
            return


# Draw UI
run_btn = widgets.Button(
    description='Run FrankPEPstein Pipeline',
    disabled=False,
    button_style='danger', # 'success', 'info', 'warning', 'danger' or ''
    tooltip='Start the magic',
    icon='rocket'
)
run_btn.on_click(run_frankpepstein_pipeline)

display(widgets.VBox([
    widgets.HBox([length_slider, num_peptides_slider]),
    threads_slider,
    run_btn,
    run_output
]))
