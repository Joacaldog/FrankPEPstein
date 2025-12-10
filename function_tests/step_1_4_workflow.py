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
            
            for i, p_file in enumerate(sorted(final_pockets_list)):
                full_path = os.path.join(pockets_dir, p_file)
                if os.path.exists(full_path):
                    with open(full_path, 'r') as f:
                        view.addModel(f.read(), "pdb")
                    
                    # Highlight selected pocket
                    if p_file == selected_pocket_file:
                         # Selected: Solid, bright, larger spheres? or just distinct standard color
                         # Let's make the selected one FLASH or be very obvious.
                         # Maybe just opaque vs transparent?
                         view.setStyle({'model': -1}, {'sphere': {'color': 'red', 'opacity': 1.0, 'radius': 1.0}})
                    else:
                        # Others: Different colors, slightly transparent
                        color = colors[i % len(colors)]
                        view.setStyle({'model': -1}, {'sphere': {'color': color, 'opacity': 0.6}})
                        
                    # Add label?
                    # view.addLabel(p_file, {'fontSize': 12, 'fontColor': 'black', 'backgroundColor': 'white'}, {'model': -1})

            view.zoomTo()
            view.show()
            
        display(widgets.interactive(view_pockets, selected_pocket_file=pocket_dropdown))
    else:
        print("No pockets available to select.")
#@title 4. Pocket Extraction & Box Generation
#@markdown This step extracts the selected pocket and calculates the grid box center and size.

import os
# --- Helper Functions (Subprocess) ---
def run_box_calculation_isolated(pdb_file, buffer=10.0):
    """
    Runs the box calculation in the FrankPEPstein environment (isolated)
    to avoid Bio.PDB binary incompatibility with the Colab kernel.
    """
    
    # 1. Create the script
    script_content = """
import sys
import os
from Bio.PDB import PDBParser

def get_box_center_size(pdb_file, buffer):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("pocket", pdb_file)
    coords = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    coords.append(atom.get_coord())
    
    if not coords:
        return None, None

    min_coord = [min([c[i] for c in coords]) for i in range(3)]
    max_coord = [max([c[i] for c in coords]) for i in range(3)]
    
    center = [(min_coord[i] + max_coord[i]) / 2 for i in range(3)]
    size = [(max_coord[i] - min_coord[i]) + float(buffer) for i in range(3)]
    
    return center, size

if __name__ == "__main__":
    pdb_path = sys.argv[1]
    buf = sys.argv[2]
    
    try:
        center, size = get_box_center_size(pdb_path, buf)
        if center:
            print(f"CENTER:{center[0]},{center[1]},{center[2]}")
            print(f"SIZE:{size[0]},{size[1]},{size[2]}")
        else:
            print("ERROR:No coordinates")
    except Exception as e:
        print(f"ERROR:{e}")
"""
    
    script_name = "calculate_box_isolated.py"
    with open(script_name, "w") as f:
        f.write(script_content)
        
    # 2. Run with isolated python
    # We assume the environment is at /usr/local/envs/FrankPEPstein
    python_exe = "/usr/local/envs/FrankPEPstein/bin/python"
    
    if not os.path.exists(python_exe):
        print(f"Error: Python executable not found at {python_exe}. Is the environment created?")
        return None, None
        
    try:
        result = subprocess.run(
            [python_exe, script_name, pdb_file, str(buffer)], 
            capture_output=True, 
            text=True, 
            check=True
        )
        
        # 3. Parse Output
        center = None
        size = None
        
        for line in result.stdout.splitlines():
            if line.startswith("CENTER:"):
                parts = line.strip().split(":")[1].split(",")
                center = [float(p) for p in parts]
            elif line.startswith("SIZE:"):
                parts = line.strip().split(":")[1].split(",")
                size = [float(p) for p in parts]
            elif line.startswith("ERROR:"):
                print(f"Script Error: {line}")
                
        return center, size
        
    except subprocess.CalledProcessError as e:
        print(f"Error executing isolated script: {e}")
        print(f"Stderr: {e.stderr}")
        return None, None

# --- GUI & Interaction Logic ---

# Output widget to capture print statements from the callback
output_log = widgets.Output()

def extract_and_calculate_box(b):
    """
    Callback function triggered by the 'Extract Box' button.
    It reads the current dropdown value and runs the isolated box calculation.
    """
    output_log.clear_output()
    with output_log:
        # Check if dropdown exists and has a value
        try:
            if 'pocket_dropdown' not in globals() or not pocket_dropdown.value:
                print("Error: No pocket selected. Please upload a receptor and detect pockets first.")
                return
            
            selected_pocket_file = pocket_dropdown.value
            print(f"Selected Pocket: {selected_pocket_file}")
            
            # pockets_dir must be defined globally from the upload step
            if 'pockets_dir' not in globals():
                 print("Error: pockets_dir not defined.")
                 return

            pocket_path = os.path.join(pockets_dir, selected_pocket_file)
            
            if os.path.exists(pocket_path):
                print(f"Path: {pocket_path}")
                print("Calculating box parameters...")
                
                # Run isolated calculation
                center, size = run_box_calculation_isolated(pocket_path, buffer=10.0)
                
                if center:
                    center_str = f"{center[0]:.3f} {center[1]:.3f} {center[2]:.3f}"
                    size_str = f"{size[0]:.3f} {size[1]:.3f} {size[2]:.3f}"
                    
                    print("-" * 30)
                    print(f"Box Center: {center_str}")
                    print(f"Box Size:   {size_str}")
                    print("-" * 30)
                    
                    # Store these in global namespace so other cells can access them if needed
                    global box_center, box_size
                    box_center = center
                    box_size = size
                    
                    print("✅ Pocket parameters ready for FrankPEPstein!")
                else:
                    print("❌ Error: Could not calculate coordinates.")
            else:
                print(f"❌ Error: File not found {pocket_path}")
        except Exception as e:
            print(f"An unexpected error occurred: {e}")

# Create Button
extract_btn = widgets.Button(
    description='Extract Pocket & Calculate Box',
    disabled=False,
    button_style='success', # 'success', 'info', 'warning', 'danger' or ''
    tooltip='Click to extract the selected pocket',
    icon='box'
)

extract_btn.on_click(extract_and_calculate_box)

# Display GUI
print("\n--- 4. Extraction Control ---")
display(extract_btn, output_log)

