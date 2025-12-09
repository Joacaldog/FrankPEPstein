#@title 4. Pocket Extraction & Box Generation
#@markdown This step extracts the selected pocket and calculates the grid box center and size.

import os
from Bio.PDB import PDBParser, PDBIO, Select

# --- Configuration ---
# Read variables from previous step (Assuming they exist in global scope)
# In Colab, variables persist across cells.
try:
    selected_pocket_file = pocket_dropdown.value
    print(f"Selected Pocket File: {selected_pocket_file}")
    
    # pockets_dir was defined in Step 1-3
    pocket_path = os.path.join(pockets_dir, selected_pocket_file)
    print(f"Pocket Path: {pocket_path}")
    
except NameError:
    print("Error: 'pocket_dropdown' or 'pockets_dir' not defined. Did you run the previous cell?")
    # Fallback for testing/debugging if not running sequentially
    # pocket_path = "pockets/pocket1.pdb" 

# --- Helper Functions ---
def get_box_center_size(pdb_file, buffer=0.0):
    """
    Calculates center and size of the box from PDB atoms.
    Buffer adds padding to the box size (total).
    """
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
    # Size is the dimension length
    size = [(max_coord[i] - min_coord[i]) + buffer for i in range(3)]
    
    return center, size

# --- Main Extraction Logic ---
if os.path.exists(pocket_path):
    print("Calculating box parameters...")
    # Buffer can be adjusted. Usually 4-10 Angstroms padding is standard for Vina
    center, size = get_box_center_size(pocket_path, buffer=10.0) 
    
    if center:
        # Round for cleaner output
        center_str = f"{center[0]:.3f} {center[1]:.3f} {center[2]:.3f}"
        size_str = f"{size[0]:.3f} {size[1]:.3f} {size[2]:.3f}"
        
        print("-" * 30)
        print(f"Box Center: {center_str}")
        print(f"Box Size:   {size_str}")
        print("-" * 30)
        
        # Save as global variables for next steps (FrankPEPstein)
        box_center = center
        box_size = size
        
        # Optional: Save a 'clean' pocket PDB if needed?
        # fpocket output is usually already clean atoms of the pocket.
        # But we might want to ensure it is standard PDB.
        
        print("Pocket parameters ready for FrankPEPstein.")
        
    else:
        print("Error: Could not calculate coordinates from pocket file.")
else:
    print(f"Error: Pocket file not found at {pocket_path}")
