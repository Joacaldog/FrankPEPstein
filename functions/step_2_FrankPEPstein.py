#@title 2. FrankPEPstein Pipeline
#@markdown **Instructions:**
#@markdown 1. Configure parameters.
#@markdown 2. Run to start the pipeline. Real-time visualization will appear below.

import os
import sys
import subprocess
import shutil
import time
import threading
import ipywidgets as widgets
from IPython.display import display, Javascript

# --- Configuration Widgets ---
style = {'description_width': 'initial'}

length_slider = widgets.IntSlider(
    value=10, min=3, max=30, step=1,
    description='Target Peptide Length:', style=style, layout=widgets.Layout(width='50%')
)
num_peptides_slider = widgets.IntSlider(
    value=10, min=1, max=100, step=1,
    description='Number of Output Peptides:', style=style, layout=widgets.Layout(width='50%')
)
threads_slider = widgets.IntSlider(
    value=os.cpu_count() or 2, min=1, max=os.cpu_count() or 8, step=1,
    description='CPU Threads:', style=style, layout=widgets.Layout(width='50%')
)

run_output = widgets.Output()

# --- Visualization Helper ---
def update_viz_js(mol_view, fragment_pdbs):
    """
    Generates JS to update the viewer with new fragments without full reload.
    """
    # This is tricky in Colab. Usually we re-render the 3D mol.
    # For simplicity, we might just re-render the view every X seconds if new files found.
    pass

# --- Pipeline Logic ---
import json

def run_frankpepstein_pipeline(b):
    run_output.clear_output()
    with run_output:
        # Load State
        state_file = "pipeline_state.json"
        global box_center, box_size, extracted_pocket_path, receptor_filename
        
        if os.path.exists(state_file):
            try:
                with open(state_file, "r") as f:
                    state = json.load(f)
                if 'box_center' in state: box_center = state['box_center']
                if 'box_size' in state: box_size = state['box_size']
                if 'extracted_pocket_path' in state: extracted_pocket_path = state['extracted_pocket_path']
                if 'receptor_filename' in state: receptor_filename = state['receptor_filename']
            except: pass

        if 'box_center' not in globals() or 'extracted_pocket_path' not in globals():
            print("Error: Missing pipeline state (Pocket/Box). Run Step 1/4 first.")
            return

        # Setup Paths
        base_dir = os.path.abspath("FrankPEPstein") if os.path.exists("FrankPEPstein") else os.path.abspath(".")
        scripts_dir = os.path.join(base_dir, "scripts")
        # DB Path: filtered_DB_P5-15_R30_id10 (extracted from ...optim.tar.gz)
        db_path = os.path.join(base_dir, "DB", "filtered_DB_P5-15_R30_id10") 
        complexes_path = db_path # User says this is the path for -d too? Or inside?
        # User said: "la carpeta que se pone como -d es filtered_DB_P5-15_R30_id10"
        # Superposer requires "-d" (database path) and "-fm" (folder minipockets? or also DB?)
        # Step 2 script sets -d to complexes_path and -fm to db_path. 
        # I will set both to this directory.
        complexes_path = db_path
        
        # Fallback for complexes if not extracted yet (User mentioned parallel download)
        # But for now we assume it exists
        
        run_dir = os.path.join(base_dir, "FrankPEPstein_Run")
        if os.path.exists(run_dir): shutil.rmtree(run_dir)
        os.makedirs(run_dir)
        
        print(f"Initializing Run in {run_dir}...")
        
        # 0. Copy Receptor
        target_receptor = "receptor.pdb"
        try:
            shutil.copy(receptor_filename, os.path.join(run_dir, target_receptor))
        except:
            print("Error copying receptor.")
            return

        # 1. Prepare Pocket (Rename Chain 'p')
        print("Preparing Pocket...")
        prep_script = os.path.join(run_dir, "prep_pocket.py")
        target_pocket = "target_pocket.pdb"
        
        with open(prep_script, "w") as f:
            f.write(f"""
import sys
from Bio import PDB
try:
    parser = PDB.PDBParser(QUIET=True)
    s = parser.get_structure("p", "{extracted_pocket_path}")
    for m in s:
        for c in m: c.id = 'p'
    io = PDB.PDBIO()
    io.set_structure(s)
    io.save("{os.path.join(run_dir, target_pocket)}")
    print("SUCCESS")
except Exception as e: print(e)
""")
        
        # Run with env python
        frank_python = "/usr/local/envs/FrankPEPstein/bin/python"
        subprocess.run([frank_python, prep_script], check=True)
        
        if not os.path.exists(os.path.join(run_dir, target_pocket)):
            print("Pocket preparation failed.")
            return

        # 2. Superposer (Step 1)
        print("\nStep 1: Fragment Scanning (Superposer)...")
        # Ensure 'click' and 'adfr' in PATH
        env = os.environ.copy()
        
        # Find paths dynamically
        utilities_dir = os.path.join(base_dir, "utilities")
        
        # Determine Minipockets DB (starts with minipockets...)
        db_root = os.path.dirname(db_path) # parent of filtered_DB
        minipockets_path = db_path # Fallback
        try:
            candidates = [d for d in os.listdir(db_root) if d.startswith("minipockets")]
            if candidates:
                minipockets_path = os.path.join(db_root, candidates[0])
        except:
             pass

        # Add Click
        click_bin = os.path.join(utilities_dir, "Click/bin")
        if os.path.exists(click_bin): env["PATH"] = f"{click_bin}:{env['PATH']}"
        
        # Add ADFR
        adfr_bin = os.path.join(utilities_dir, "ADFRsuite_x86_64Linux_1.0/bin") # Common struct
        if os.path.exists(adfr_bin): env["PATH"] = f"{adfr_bin}:{env['PATH']}"

        superposer_script = os.path.join(scripts_dir, "superposerV5.2_leave1out_cluster.py")
        
        cmd_super = [
            frank_python, superposer_script,
            "-T", target_pocket,
            "-d", complexes_path,
            "-a", "3", "-r", "0.1",
            "-x_center", str(box_center[0]), "-y_center", str(box_center[1]), "-z_center", str(box_center[2]),
            "-x_size", str(box_size[0]), "-y_size", str(box_size[1]), "-z_size", str(box_size[2]),
            "-t", str(threads_slider.value),
            "-fm", minipockets_path
        ]
        
        # Start Visualization Thread
        # We will monitor the output folder of superposer
        # output folder format: superpockets_residuesAligned3_RMSD0.1
        super_out_dir = os.path.join(run_dir, "superpockets_residuesAligned3_RMSD0.1")
        
        # viz_monitor writes to global viz_output checking run folder
        
        stop_event = threading.Event()
        
        def viz_monitor():
            import py3Dmol
            import time
            from IPython.display import display, clear_output
            
            # Initial Wait
            time.sleep(2)
            
            seen_fragments = 0
            
            while not stop_event.is_set():
                # 1. Check for fragments
                fragments = []
                if os.path.exists(super_out_dir):
                    fragments = [f for f in os.listdir(super_out_dir) if f.startswith("patch_file")]
                
                n_frags = len(fragments)
                
                # Update every 10s or if first time
                # We force update for feedback
                
                with viz_output:
                    clear_output(wait=True)
                    
                    print(f"Visualization Update: {n_frags} fragments found.")
                    
                    view = py3Dmol.view(width=800, height=600)
                    
                    # A. Receptor
                    view.addModel(open(receptor_filename, 'r').read(), "pdb")
                    view.setStyle({'model': 0}, {}) 
                    view.addSurface(py3Dmol.SES, {'opacity': 0.3, 'color': 'white'}, {'model': 0})
                    
                    # B. Box (Center/Size)
                    # Convert center/size to min/max
                    # center is [x, y, z], size is [dx, dy, dz]
                    b_min = [box_center[i] - box_size[i]/2 for i in range(3)]
                    b_max = [box_center[i] + box_size[i]/2 for i in range(3)]
                    
                    view.addBox({
                        'center': {'x': box_center[0], 'y': box_center[1], 'z': box_center[2]},
                        'dimensions': {'w': box_size[0], 'h': box_size[1], 'd': box_size[2]},
                        'color': 'cyan',
                        'opacity': 0.5,
                        'wireframe': True
                    })
                    
                    # C. Fragments (Sample if too many?)
                    # Loading all might be heavy if thousands. Load last 50?
                    # or random sample?
                    # Let's show up to 50 random ones to keep it responsive
                    import random
                    viz_frags = fragments
                    if len(viz_frags) > 50:
                         viz_frags = random.sample(fragments, 50)
                         
                    for frag_file in viz_frags:
                        frag_path = os.path.join(super_out_dir, frag_file)
                        try:
                            view.addModel(open(frag_path, 'r').read(), "pdb")
                            # Style fragments as sticks, colorful?
                            # Last added model is -1
                            view.setStyle({'model': -1}, {'stick': {'colorscheme': 'greenCarbon'}})
                        except: pass
                        
                    view.zoomTo()
                    view.show()
                
                time.sleep(10)
                
        viz_thread = threading.Thread(target=viz_monitor)
        viz_thread.start()
        
        try:
            # Run Superposer
            # Check DB exist
            if not os.path.exists(complexes_path):
                 print(f"Warning: Complex DB not found at {complexes_path}")
            
            process = subprocess.Popen(cmd_super, cwd=run_dir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env)
            
            for line in iter(process.stdout.readline, ''):
                print(line, end='')
                
            process.wait()
            
            stop_event.set()
            viz_thread.join()
            
            print("\nSuperposer Completed.")
            
        except Exception as e:
            stop_event.set()
            print(f"Error: {e}")
            return

        # 3. Patch Clustering (Step 2 & 3)
        # Note: The original pipeline split this into FrankVINA I (filtering) -> PatchClust -> FrankVINA II (scoring)
        # We will follow that.
        
        pep_len = length_slider.value
        n_peps = num_peptides_slider.value
        
        # A. Filter (FrankVINA I)
        print("\nStep 2: Patch Filtering...")
        if not os.path.exists(super_out_dir):
            print("No Superposer output.")
            return

        # Copy receptor to output dir for Vina
        shutil.copy(os.path.join(run_dir, target_receptor), os.path.join(super_out_dir, target_receptor))
        
        vina_script_1 = os.path.join(scripts_dir, "frankVINA_FNKPSTN_cluster.py")
        subprocess.run([frank_python, vina_script_1, target_receptor, str(threads_slider.value)], cwd=super_out_dir, check=True, env=env)
        
        # B. Cluster
        print("\nStep 3: Clustering...")
        patches_dir = os.path.join(super_out_dir, "top_10_patches")
        if not os.path.exists(patches_dir):
            print("No qualified patches found.")
            return
            
        clust_script = os.path.join(scripts_dir, "patch_clustering_V8.7_cluster.py")
        subprocess.run([frank_python, clust_script, "-w", str(pep_len), "-t", str(threads_slider.value)], cwd=patches_dir, check=True, env=env)
        
        # C. Score (FrankVINA II)
        print("\nStep 4: Final Scoring...")
        final_dir = os.path.join(patches_dir, f"frankPEPstein_{pep_len}")
        if not os.path.exists(final_dir):
             print("Clustering failed.")
             return
        
        shutil.copy(os.path.join(run_dir, target_receptor), os.path.join(final_dir, target_receptor))
        vina_script_2 = os.path.join(scripts_dir, "frankVINA_V3_cluster.py")
        subprocess.run([frank_python, vina_script_2, target_receptor, str(threads_slider.value), str(n_peps)], cwd=final_dir, check=True, env=env)
        
        print("\nPipeline Finished!")
        # Show Top Results?
        results_file = os.path.join(final_dir, f"top_{n_peps}_peps", f"top{n_peps}_peps.tsv")
        if os.path.exists(results_file):
             print(f"Results saved to: {results_file}")
             # Display top 5
             with open(results_file) as f:
                 print("\nTop Candidates:")
                 print(f.read())

run_btn = widgets.Button(description='Run Pipeline', button_style='danger', icon='rocket')
run_btn.on_click(run_frankpepstein_pipeline)


viz_output = widgets.Output()

display(widgets.VBox([
    widgets.HBox([length_slider, num_peptides_slider]),
    threads_slider,
    run_btn,
    run_output,
    viz_output
]))

