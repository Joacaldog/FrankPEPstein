#@title 1. Input & Pocket Selection
#@markdown **Instructions:**
#@markdown 1. Upload your Receptor PDB.
#@markdown 2. Choose Mode: **Auto Detect** (runs fpocket) or **Manual Upload** (upload your specific pocket PDB).
#@markdown 3. Select the pocket from the dropdown to visualize.

import os
import subprocess
import py3Dmol
import ipywidgets as widgets
from google.colab import files
from IPython.display import display

# --- configuration ---
detection_mode = "Auto Detect" #@param ["Auto Detect", "Manual Upload"]

# Global variables for next steps
receptor_filename = None
pockets_dir = "pockets_upload" # Default for manual
final_pockets_list = []

# --- 1. Upload Receptor ---
print(f"--- Upload Receptor PDB ({detection_mode}) ---")
uploaded_r = files.upload()

if not uploaded_r:
    print("No receptor file uploaded.")
else:
    receptor_filename = list(uploaded_r.keys())[0]
    print(f"Receptor: {receptor_filename}")

    # --- 2. Pocket Handling ---
    if detection_mode == "Auto Detect":
        print(f"\nRunning fpocket on {receptor_filename}...")
        try:
            # Fix: Quotes for filenames with spaces
            subprocess.run(f"fpocket -f '{receptor_filename}'", shell=True, check=True)
            
            # Robust folder finding
            base_name = os.path.splitext(receptor_filename)[0]
            possible_folders = [f"{receptor_filename}_out", f"{base_name}_out"]
            output_folder = next((f for f in possible_folders if os.path.exists(f)), None)

            if output_folder:
                pockets_dir = os.path.join(output_folder, "pockets")
                if os.path.exists(pockets_dir):
                    final_pockets_list = [f for f in os.listdir(pockets_dir) if f.endswith(".pdb")]
                    print(f"Auto-detection finished. Found {len(final_pockets_list)} pockets.")
                else:
                    print("Warning: pockets subdirectory not found.")
            else:
                print("Error: fpocket output not found.")
                
        except subprocess.CalledProcessError:
             print("Error running fpocket.")

    elif detection_mode == "Manual Upload":
        print(f"\n--- Upload Manual Pocket PDB ---")
        os.makedirs(pockets_dir, exist_ok=True)
        uploaded_p = files.upload()
        if uploaded_p:
            for p_file in uploaded_p.keys():
                # Move to a pockets folder to keep structure consistent
                os.rename(p_file, os.path.join(pockets_dir, p_file))
                final_pockets_list.append(p_file)
            print(f"Manual upload finished. Available pockets: {len(final_pockets_list)}")

    # --- 3. Visualization & Selection ---
    if final_pockets_list:
        print("\n--- Pocket Selection & Visualization ---")
        
        pocket_dropdown = widgets.Dropdown(
            options=sorted(final_pockets_list),
            description='Select Pocket:',
            disabled=False,
        )

        def view_pocket(pocket_file):
            view = py3Dmol.view(width=800, height=600)
            
            # 1. Receptor Surface (White, Transparent)
            with open(receptor_filename, 'r') as f:
                view.addModel(f.read(), "pdb")
            view.setStyle({}) 
            view.addSurface(py3Dmol.SES, {'opacity': 0.8, 'color': 'white'})
            
            # 2. Selected Pocket (Red Spheres)
            # Ensure we look in the correct dir (either fpocket out or manual upload dir)
            full_path = os.path.join(pockets_dir, pocket_file)
            if os.path.exists(full_path):
                with open(full_path, 'r') as f:
                    view.addModel(f.read(), "pdb")
                view.setStyle({'model': -1}, {'sphere': {'color': 'red', 'opacity': 0.7}})
            else:
                print(f"Error: Could not find {full_path}")

            view.zoomTo()
            view.show()
            
            # Set a global var for the *path* so next cell can find it easily? 
            # Actually next cell will read `pocket_dropdown.value` and `pockets_dir`
            
        display(widgets.interactive(view_pocket, pocket_file=pocket_dropdown))
    else:
        print("No pockets available to select.")
