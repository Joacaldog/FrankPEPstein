
import os
import time
import subprocess
import shutil
import py3Dmol
from IPython.display import display
import ipywidgets as widgets

def run_superposer_with_viz(cmd_string, run_dir, receptor_path, box_center=None, box_size=None, output_dir_name="superpockets_residuesAligned3_RMSD0.1"):
    """
    Executes the Superposer command while visualizing found fragments in real-time.
    
    Args:
        cmd_string (str): The full command to run superposer.
        run_dir (str): Directory where the run takes place.
        receptor_path (str): Path to the receptor PDB file.
        box_center (list/tuple): [x, y, z] coordinates of the box center.
        box_size (list/tuple): [x, y, z] dimensions of the box.
        output_dir_name (str): Name of the folder where fragments are saved.
    """
    
    print(f"Working Directory: {run_dir}")
    print(f"Command: {cmd_string}")
    
    # Store current dir to restore later
    original_cwd = os.getcwd()
    
    try:
        if not os.path.exists(run_dir):
            os.makedirs(run_dir)
        
        os.chdir(run_dir)
        
        # Initialize Viewer
        view = py3Dmol.view(width=800, height=600)
        
        # Load Receptor
        if os.path.exists(receptor_path):
            with open(receptor_path, 'r') as f:
                view.addModel(f.read(), "pdb")
            view.setStyle({'cartoon': {'color': 'white', 'opacity': 0.8}})
        else:
            print(f"Warning: Receptor not found at {receptor_path}")

        # Visualize Search Box
        if box_center is not None and box_size is not None:
            # py3Dmol expects dicts for center and dimensions
            # Dimensions in py3Dmol addBox are typically w, h, d
            try:
                view.addBox({
                    'center': {'x': box_center[0], 'y': box_center[1], 'z': box_center[2]},
                    'dimensions': {'w': box_size[0], 'h': box_size[1], 'd': box_size[2]},
                    'color': 'cyan',
                    'opacity': 0.3
                })
                # Add wireframe box for better visibility of edges
                view.addBox({
                    'center': {'x': box_center[0], 'y': box_center[1], 'z': box_center[2]},
                    'dimensions': {'w': box_size[0], 'h': box_size[1], 'd': box_size[2]},
                    'color': 'blue',
                    'wireframe': True
                })
            except Exception as e:
                print(f"Error adding box visualization: {e}")

        view.zoomTo()
        view.show()
        display_handle = display(view, display_id=True)
        
        # Output widget for logs
        log_output = widgets.Output()
        display(log_output)
        
        # Start Process
        # We use shell=True if cmd_string is a string, or split it. 
        # Splitting is safer for Popen with shell=False.
        process = subprocess.Popen(
            cmd_string, 
            shell=True,
            stdout=subprocess.PIPE, 
            stderr=subprocess.STDOUT, 
            text=True,
            bufsize=1
        )
        
        seen_fragments = set()
        out_path = os.path.abspath(output_dir_name)
        
        start_time = time.time()
        
        # Monitoring Loop
        while process.poll() is None:
            # Check for new fragments
            if os.path.exists(out_path):
                current_files = set([f for f in os.listdir(out_path) if f.endswith(".pdb")])
                new_files = current_files - seen_fragments
                
                if new_files:
                    with log_output:
                        print(f"[{time.strftime('%H:%M:%S')}] Found {len(new_files)} new fragments!")
                    
                    for frag_file in new_files:
                        full_frag_path = os.path.join(out_path, frag_file)
                        try:
                            with open(full_frag_path, 'r') as f:
                                # Add fragment
                                view.addModel(f.read(), "pdb")
                                # Style last added model (fragment)
                                view.setStyle({'model': -1}, {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.2}})
                        except Exception as e:
                            with log_output:
                                print(f"Error loading {frag_file}: {e}")
                    
                    seen_fragments.update(new_files)
                    display_handle.update(view)
            
            # Read stdout specific lines if needed? 
            # Or just let it run. Reading stdout might block if not careful.
            # Using read1 or similar is complex in blocking loop. 
            # For now, let's just loop.
            
            time.sleep(2)
        
        # Final check
        if process.returncode == 0:
            print("Superposer completed successfully.")
        else:
            print(f"Superposer failed with return code {process.returncode}")
            print(process.stdout.read())

    except Exception as e:
        print(f"An error occurred: {e}")
    finally:
        os.chdir(original_cwd)

