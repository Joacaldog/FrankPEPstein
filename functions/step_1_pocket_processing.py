#@title 1. Input & Pocket Selection
#@markdown **Instructions:**
#@markdown 1. Upload your Receptor PDB.
#@markdown 2. Choose Mode: **Auto Detect** (runs fpocket) or **Manual Upload** (upload your specific pocket PDB).
#@markdown 3. Select the pocket from the dropdown to visualize.

import os
import sys
import subprocess
try:
    import py3Dmol
except ImportError:
    # Try adding FrankPEPstein env to path
    env_path = "/usr/local/envs/FrankPEPstein"
    site_packages = f"{env_path}/lib/python3.10/site-packages"
    if os.path.exists(site_packages):
        if site_packages not in sys.path:
            sys.path.append(site_packages)
        # Add binary path too
        if f"{env_path}/bin" not in os.environ['PATH']:
            os.environ['PATH'] = f"{env_path}/bin:" + os.environ['PATH']
    
    # Retry import
    try:
        import py3Dmol
    except ImportError:
        print("py3Dmol not found. Installing...")
        subprocess.run("pip install -q py3dmol", shell=True, check=True)
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

# --- Persistence Logic ---
import json

def save_pipeline_state(updates):
    state_file = "pipeline_state.json"
    current_state = {}
    if os.path.exists(state_file):
        try:
            with open(state_file, "r") as f:
                current_state = json.load(f)
        except:
            pass
    current_state.update(updates)
    with open(state_file, "w") as f:
        json.dump(current_state, f, indent=4)
    print(f"State saved to {state_file}")

# --- 1. Upload Receptor ---
print(f"--- Upload Receptor PDB ({detection_mode}) ---")
uploaded_r = files.upload()

import re

if not uploaded_r:
    print("No receptor file uploaded.")
else:
    raw_filename = list(uploaded_r.keys())[0]
        
    # Check for Colab duplicate naming (e.g. receptor(1).pdb)
    match = re.search(r'^(.*?)\s?\(\d+\)(\.[^.]*)?$', raw_filename)
    if match:
        clean_name = match.group(1) + (match.group(2) if match.group(2) else "")
        print(f"Detected duplicate upload: {raw_filename} -> overwriting {clean_name}")
        
        if os.path.exists(clean_name):
            os.remove(clean_name)
        os.rename(raw_filename, clean_name)
        receptor_filename = os.path.abspath(clean_name)
    else:
        receptor_filename = os.path.abspath(raw_filename)
        
    print(f"Receptor: {receptor_filename}")
    save_pipeline_state({"receptor_filename": receptor_filename})

    # --- 2. Pocket Handling ---
    if detection_mode == "Auto Detect":
        print(f"\nRunning fpocket on {receptor_filename}...")
        try:
            # Fix: Quotes for filenames with spaces
            subprocess.run(f"fpocket -f '{receptor_filename}'", shell=True, check=True)
            
            # Robust folder finding
            # fpocket creates output based on the filename in the SAME directory
            # receptor_filename is absolute, so output should be absolute too?
            # fpocket output format: /path/to/file_out/
            
            base_name_no_ext = os.path.splitext(os.path.basename(receptor_filename))[0]
            base_dir = os.path.dirname(receptor_filename)
            
            # Possible output folder names
            folder_name_1 = f"{os.path.basename(receptor_filename)}_out"
            folder_name_2 = f"{base_name_no_ext}_out"
            
            # Check in the same directory as the receptor
            possible_folders = [
                os.path.join(base_dir, folder_name_1),
                os.path.join(base_dir, folder_name_2)
            ]
            
            output_folder = next((f for f in possible_folders if os.path.exists(f)), None)

            if output_folder:
                pockets_dir = os.path.join(output_folder, "pockets")
                if os.path.exists(pockets_dir):
                    final_pockets_list = [f for f in os.listdir(pockets_dir) if f.endswith(".pdb")]
                    print(f"Auto-detection finished. Found {len(final_pockets_list)} pockets.")
                    pockets_dir = os.path.abspath(pockets_dir) # Ensure absolute
                else:
                    print(f"Warning: pockets subdirectory not found in {output_folder}")
            else:
                print("Error: fpocket output not found. Checked:", possible_folders)
                
        except subprocess.CalledProcessError:
             print("Error running fpocket.")

    elif detection_mode == "Manual Upload":
        print(f"\n--- Upload Manual Pocket PDB ---")
        os.makedirs(pockets_dir, exist_ok=True)
        uploaded_p = files.upload()
        import re
        if uploaded_p:
            for p_file in uploaded_p.keys():
                # Colab renames duplicate uploads to filename(1).ext. 
                # User wants to overwrite instead.
                
                # Check for pattern like "name(1).pdb" or "name (1).pdb"
                # Regex matches: (any content) optional space (digits) (extension)
                match = re.search(r'^(.*?)\s?\(\d+\)(\.[^.]*)?$', p_file)
                
                if match:
                    clean_name = match.group(1) + (match.group(2) if match.group(2) else "")
                    print(f"Detected duplicate upload: {p_file} -> overwriting {clean_name}")
                else:
                    clean_name = p_file

                target_path = os.path.join(pockets_dir, clean_name)
                
                # If target exists, log that we are replacing it
                if os.path.exists(target_path):
                    print(f"Replacing existing file: {clean_name}")
                    os.remove(target_path)
                
                # Move (rename) the uploaded file to the target path
                # Note: 'p_file' is in CWD (content/), target is in pockets_dir
                os.rename(p_file, target_path)
                
                if clean_name not in final_pockets_list:
                    final_pockets_list.append(clean_name)
                    
            print(f"Manual upload finished. Available pockets: {len(final_pockets_list)}")

    # --- 3. Visualization & Selection ---
    if final_pockets_list:
        print("\n--- Pocket Selection & Visualization ---")
        print("Displaying all detected pockets. Select one below for extraction.")
        
        pocket_dropdown = widgets.Dropdown(
            options=sorted(final_pockets_list),
            description='Select Pocket:',
            disabled=False,
        )

        def view_pockets(selected_pocket_file):
            view = py3Dmol.view(width=800, height=600)
            
            # 1. Receptor Surface (White, Transparent)
            with open(receptor_filename, 'r') as f:
                view.addModel(f.read(), "pdb")
            view.setStyle({}) 
            view.addSurface(py3Dmol.SES, {'opacity': 0.3, 'color': 'white'})
            
            # 2. Add ALL pockets with distinct colors
            colors = ['#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#FF00FF', '#00FFFF', '#FFA500', '#800080', '#008000', '#800000']
            
            # Regex to find number in filename (e.g. pocket5_atm.pdb -> 5, or just 5.pdb -> 5)
            import re
            
            for i, p_file in enumerate(sorted(final_pockets_list)):
                full_path = os.path.join(pockets_dir, p_file)
                if os.path.exists(full_path):
                    with open(full_path, 'r') as f:
                        view.addModel(f.read(), "pdb")
                    
                    # Try to extract a short label (number)
                    # Common patterns: "pocket5_atm.pdb", "pocket5.pdb", "5.pdb"
                    match = re.search(r'(\d+)', p_file)
                    label_text = match.group(1) if match else p_file
                    
                    # Determine styling
                    is_selected = (p_file == selected_pocket_file)
                    
                    if is_selected:
                         color = 'red'
                         opacity = 1.0
                         label_style = {'fontSize': 18, 'fontColor': 'red', 'backgroundColor': 'white', 'backgroundOpacity': 0.8, 'border': '2px solid red'}
                    else:
                        color = colors[i % len(colors)]
                        opacity = 0.6
                        label_style = {'fontSize': 12, 'fontColor': 'black', 'backgroundColor': 'white', 'backgroundOpacity': 0.5}

                    view.setStyle({'model': -1}, {'sphere': {'color': color, 'opacity': opacity}})
                    
                    # Add 3D Label (Number only)
                    view.addLabel(label_text, label_style, {'model': -1})

            view.zoomTo()
            view.show()
            
        display(widgets.interactive(view_pockets, selected_pocket_file=pocket_dropdown))
    else:
        print("No pockets available to select.")
#@title 4. Pocket Extraction & Box Generation
#@markdown This step extracts the selected pocket (residues within 5Å of fpocket spheres) and calculates the grid box.

import os
import subprocess # Added subprocess import here for consistency
import sys # Added sys import for python_exe fallback

# --- Helper Functions (Subprocess) ---
def run_extraction_and_box_isolated(receptor_path, pocket_path, output_pocket_path, buffer=5.0):
    """
    Runs extraction (NeighborSearch) and box calculation in isolated environment.
    """
    
    script_content = """
import sys
import os
from Bio import PDB
from Bio.PDB import PDBParser, PDBIO, Select, NeighborSearch

def extract_and_box(receptor_file, pocket_atm_file, output_file, buffer_val):
    try:
        parser = PDBParser(QUIET=True)
        
        # 1. Load Structures
        # print(f"Loading receptor: {receptor_file}")
        receptor_struct = parser.get_structure("receptor", receptor_file)
        
        # print(f"Loading pocket atoms: {pocket_atm_file}")
        pocket_struct = parser.get_structure("pocket_atm", pocket_atm_file)
        
        # 2. Collect Pocket Atoms (e.g. Alpha Spheres)
        pocket_atoms = []
        for model in pocket_struct:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        pocket_atoms.append(atom)
        
        if not pocket_atoms:
            print("ERROR: No atoms in pocket file")
            return

        # 3. Neighbor Search (5.0 Angstrom)
        # Find PROTEIN residues near the pocket spheres
        receptor_atoms = list(receptor_struct.get_atoms())
        ns = NeighborSearch(receptor_atoms)
        
        selected_residues = set()
        
        for p_atom in pocket_atoms:
            # Level 'R' returns residues
            nearby_residues = ns.search(p_atom.get_coord(), 5.0, level='R')
            for res in nearby_residues:
                selected_residues.add((res.parent.id, res.id))
                
        # 4. Save Extracted Pocket (Chain 'p')
        class PocketSelect(Select):
            def accept_residue(self, residue):
                return (residue.parent.id, residue.id) in selected_residues

        io = PDBIO()
        io.set_structure(receptor_struct)
        
        # We need to save, but also RENAME chain to 'p'
        # To do this cleanly without modifying original struct in memory too much:
        # Save to temp, reload, rename, save final.
        
        temp_out = output_file + ".tmp"
        io.save(temp_out, PocketSelect())
        
        # Reload temp to verify and calculate box
        extracted_struct = parser.get_structure("extracted", temp_out)
        
        all_coords = []
        for model in extracted_struct:
            for chain in model:
                chain.id = 'p' # RENAME TO p
                for residue in chain:
                    for atom in residue:
                        all_coords.append(atom.get_coord())
        
        # Save Final
        io.set_structure(extracted_struct)
        io.save(output_file)
        os.remove(temp_out)
        
        # 5. Box Calculation
        if not all_coords:
            print("ERROR: No residues selected")
            return

        min_coord = [min([c[i] for c in all_coords]) for i in range(3)]
        max_coord = [max([c[i] for c in all_coords]) for i in range(3)]
        
        # Basic Center
        center = [(min_coord[i] + max_coord[i]) / 2 for i in range(3)]
        
        # Size + Buffer
        # buffer_val is passed (usually 10 to add padding)
        # User requested trick? Standard padding is usually sufficient if 5A extraction was done.
        size = [(max_coord[i] - min_coord[i]) + float(buffer_val) for i in range(3)]
        
        print(f"CENTER:{center[0]},{center[1]},{center[2]}")
        print(f"SIZE:{size[0]},{size[1]},{size[2]}")
        print("SUCCESS")

    except Exception as e:
        print(f"ERROR:{e}")

if __name__ == "__main__":
    extract_and_box(sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]))

"""
    script_name = "extract_box_isolated.py"
    with open(script_name, "w") as f:
        f.write(script_content)
        
    python_exe = "/usr/local/envs/FrankPEPstein/bin/python"
    # Fallback to sys.executable if env not found (local testing)
    if not os.path.exists(python_exe): python_exe = sys.executable

    try:
        result = subprocess.run(
            [python_exe, script_name, receptor_path, pocket_path, output_pocket_path, str(buffer)],
            capture_output=True, text=True, check=True
        )
        
        center, size = None, None
        success = False
        
        for line in result.stdout.splitlines():
            if line.startswith("CENTER:"):
                center = [float(x) for x in line.split(":")[1].split(",")]
            elif line.startswith("SIZE:"):
                size = [float(x) for x in line.split(":")[1].split(",")]
            elif line.startswith("SUCCESS"):
                success = True
            elif line.startswith("ERROR"):
                print(f"Script Error: {line}")
                
        return center, size, success

    except subprocess.CalledProcessError as e:
        print(f"Execution Error: {e.stderr}")
        return None, None, False

# --- GUI ---
output_log = widgets.Output()

def extract_and_calculate_box(b):
    output_log.clear_output()
    with output_log:
        if 'pocket_dropdown' not in globals() or not pocket_dropdown.value:
            print("No pocket selected.")
            return
        
        if 'receptor_filename' not in globals() or not receptor_filename:
             print("Receptor not defined.")
             return

        selected_pocket = pocket_dropdown.value
        pocket_path = os.path.join(pockets_dir, selected_pocket)
        
        # Define output path for the extracted (real) pocket
        # stored in the same dir as receptor usually, or pockets dir
        extracted_filename = "pocket.pdb"
        extracted_path = os.path.join(pockets_dir, extracted_filename)
        
        print(f"Processing {selected_pocket}...")
        print(f"Receptor: {receptor_filename}")
        print("Running NeighborSearch (5Å) + Box Calculation...")
        
        center, size, success = run_extraction_and_box_isolated(
            receptor_filename, pocket_path, extracted_path, buffer=10.0 # Buffer for BOX size (padding), extraction radius is 5.0 inside script
        )
        
        if success and center:
            center_str = f"{center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f}"
            size_str = f"{size[0]:.3f}, {size[1]:.3f}, {size[2]:.3f}"
            
            print("-" * 30)
            print(f"Box Center: {center_str}")
            print(f"Box Size:   {size_str}")
            print("-" * 30)
            print(f"✅ Created Extracted Pocket: {extracted_path} (pocket.pdb)")
            
            global box_center, box_size, extracted_pocket_path
            box_center = center
            box_size = size
            extracted_pocket_path = os.path.abspath(extracted_path)
            
            save_pipeline_state({
                "box_center": box_center,
                "box_size": box_size,
                "extracted_pocket_path": extracted_pocket_path
            })
        else:
            print("Extraction Failed.")

extract_btn = widgets.Button(
    description='Extract Pocket & Calculate Box',
    button_style='success',
    icon='box',
    layout=widgets.Layout(width='50%')
)
extract_btn.on_click(extract_and_calculate_box)

print("\n--- 4. Extraction Control ---")
display(extract_btn, output_log)

