import os
import sys
import glob
import shutil

def configure_modeller(license_key='MODELIRANJE', repo_dir='FrankPEPstein'):
    """
    Configures Modeller by locating the config.py file in the installation
    and replacing the license key placeholder with the provided key.
    """
    # Template location in the repo
    template_config = os.path.join(repo_dir, "utilities/config.py")
    
    # Search for config.py in multiple possible locations
    possible_paths = [
        f"{sys.prefix}/lib/modeller-*/modlib/modeller/config.py", # Standard standalone install
        f"{sys.prefix}/lib/python*/site-packages/modeller/config.py" # Site-packages install
    ]
    
    dest_config_paths = []
    for pattern in possible_paths:
        dest_config_paths.extend(glob.glob(pattern))
    
    if dest_config_paths and os.path.exists(template_config):
        dest_config = dest_config_paths[0]
        print(f"Found modeller config at: {dest_config}")
        print(f"Using template {template_config} to update {dest_config}")
        
        with open(template_config, 'r') as f:
            content = f.read()
        
        # Replace placeholder 'MODELIRANJE' with actual key
        new_content = content.replace("'MODELIRANJE'", f"'{license_key}'")
        
        with open(dest_config, 'w') as f:
            f.write(new_content)
        print("Modeller configured successfully.")
        return True
    else:
        print("Error: Modeller config destination or template not found.")
        return False

def get_pocket_box(pdb_file):
    """
    Calculates the center and size of a box surrounding the atoms in the given PDB file.
    Adds a buffer of 5.0 units to the size.
    """
    import Bio.PDB
    parser = Bio.PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("pocket", pdb_file)
    coords = []
    for atom in structure.get_atoms():
        coords.append(atom.get_coord())
    
    if not coords:
        return None, None

    min_coord = [min([c[i] for c in coords]) for i in range(3)]
    max_coord = [max([c[i] for c in coords]) for i in range(3)]
    
    center = [(min_coord[i] + max_coord[i]) / 2 for i in range(3)]
    size = [(max_coord[i] - min_coord[i]) + 5.0 for i in range(3)] # Add buffer
    return center, size

def patch_scripts(scripts_dir, path_replacements):
    """
    Iterates through .py files in scripts_dir and applies string replacements.
    """
    print("Patching scripts...")
    count = 0
    for script_name in os.listdir(scripts_dir):
        if script_name.endswith(".py"):
            full_path = os.path.join(scripts_dir, script_name)
            with open(full_path, 'r') as f:
                content = f.read()
            
            original_content = content
            for old, new in path_replacements.items():
                content = content.replace(old, new)
            
            # Additional patches for command calls
            content = content.replace("vina_1.2.4_linux_x86_64", "vina") 
            
            if content != original_content:
                with open(full_path, 'w') as f:
                    f.write(content)
                print(f"Patched {script_name}")
                count += 1
    return count
