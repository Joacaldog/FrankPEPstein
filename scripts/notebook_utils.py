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
    
    # Try using python import to find the location
    dest_config = None
    try:
        import modeller
        modeller_path = os.path.dirname(modeller.__file__)
        candidate = os.path.join(modeller_path, "config.py")
        if os.path.exists(candidate):
            dest_config = candidate
    except Exception:
        # Modeller raises an error on import if not configured, which is expected.
        # We just want to find where it is installed.
        pass

    # Fallback to search if import finding failed
    if not dest_config:
        possible_paths = [
            f"{sys.prefix}/lib/modeller-*/modlib/modeller/config.py", # Standard standalone install
            f"{sys.prefix}/lib/python*/site-packages/modeller/config.py", # Site-packages install
            f"{sys.prefix}/pkgs/modeller-*/lib/modeller-*/modlib/modeller/config.py", # Conda pkgs cache structure
            "/usr/local/envs/FrankPEPstein/lib/modeller-*/modlib/modeller/config.py" # Targeted Conda Environment
        ]
        
        dest_config_paths = []
        for pattern in possible_paths:
            found = glob.glob(pattern)
            dest_config_paths.extend(found)
        
        if dest_config_paths:
            dest_config = dest_config_paths[0]

    
    if dest_config and os.path.exists(template_config):
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
        print(f"Error: Modeller config destination ({dest_config}) or template ({template_config}) not found.")
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
            content = content.replace("vina", "vina") 
            
            if content != original_content:
                with open(full_path, 'w') as f:
                    f.write(content)
                print(f"Patched {script_name}")
                count += 1
    return count

def setup_external_tools(drive_ids=None):
    """
    Sets up external tools (ADFR, Click, DB).
    If drive_ids is provided, downloads missing files from Google Drive.
    """
    import subprocess
    
    # Ensure gdown is installed
    try:
        import gdown
    except ImportError:
        print("Installing gdown...")
        subprocess.run([sys.executable, "-m", "pip", "install", "gdown"], check=True)
        import gdown

    base_dir = "FrankPEPstein"
    # Adjust base_dir if we are running from root vs inside scripts?
    # The utils assume repo_dir='FrankPEPstein' usually implies subfolder.
    # But if looking for "utilities", it usually expects to find them relative to CWD?
    # Let's check config.
    # In notebook setup: repo_path = os.path.abspath("FrankPEPstein").
    # If we run cell_01_setup.py from FrankPEPstein root, base_dir "FrankPEPstein" might be wrong if we are IN it?
    # notebook_utils.py:
    #   configure_modeller default repo_dir='FrankPEPstein'.
    #   BUT when running locally in the repo, 'FrankPEPstein' folder DOES NOT EXIST inside 'FrankPEPstein'.
    #   The repo IS the cwd.
    #   When cloning in Colab: cwd is /content/, repo is /content/FrankPEPstein.
    #   So 'FrankPEPstein/utilities' is correct there.
    #   But LOCALLY, if I am in ~/FrankPEPstein/, 'FrankPEPstein/utilities' does not exist. 'utilities' exists.
    
    # I need to handle this path difference!
    
    if os.path.exists("utilities"):
        # We are likely INSIDE the repo root (Local execution)
        base_dir = "."
    elif os.path.exists("FrankPEPstein/utilities"):
        # We are likely in parent dir (Colab default)
        base_dir = "FrankPEPstein"
    else:
        # Fallback or create?
        base_dir = "FrankPEPstein" # Default to colab behavior for safety, or create it.

    utilities_dir = os.path.join(base_dir, "utilities")
    db_dir = os.path.join(base_dir, "DB")
    
    os.makedirs(utilities_dir, exist_ok=True)
    os.makedirs(db_dir, exist_ok=True)

    # File definitions
    # --- Bundle Download Logic ---
    # Downloads everything in two main packages if IDs are provided
    bundles = {
        "utilities_pkg": {
            "path": os.path.join(base_dir, "utilities.tar.gz"),
            "id_key": "utilities_pkg_id",
            "extract_to": os.path.join(base_dir, "utilities"),
            "desc": "Utilities Bundle"
        },
        "db_pkg": {
            "path": os.path.join(base_dir, "DB.tar.gz"),
            "id_key": "db_pkg_id",
            "extract_to": os.path.join(base_dir, "DB"),
            "desc": "Database Bundle"
        }
    }

    if drive_ids is None:
        drive_ids = {}

    # Download and extract bundles first
    for name, info in bundles.items():
        bundle_id = drive_ids.get(info["id_key"])
        if bundle_id and not os.path.exists(info["extract_to"]): # Only if dir doesn't exist? Or check manifest?
             # Actually, we should check if the CONTENT exists, but downloading bundle is safer if unsure.
             # Simple check: if tarball doesn't exist, download.
             if not os.path.exists(info["path"]):
                 print(f"Downloading {info['desc']}...")
                 url = f'https://drive.google.com/uc?id={bundle_id}'
                 gdown.download(url, info["path"], quiet=False)
             
             # Extract
             if os.path.exists(info["path"]):
                 print(f"Extracting {info['desc']}...")
                 os.makedirs(info["extract_to"], exist_ok=True)
                 # strip-components=0 because we tarred content of utilies into utilities.tar.gz?
                 # I tarred with -C utilities, so it contains "ADFR..." at root.
                 subprocess.run(f"tar -xzf {info['path']} -C {info['extract_to']}", shell=True, check=True)

    # --- Individual Tool Verification ---
    # Even after bundle extraction, we run this to ensure paths are set and bins are executable.
    # It also handles legacy cases (individual IDs provided).
    files = {
        "adfr": {
            "path": os.path.join(utilities_dir, "ADFRsuite_x86_64Linux_1.0.tar.gz"),
            "id_key": "adfr_id", # Fallback key
            "extract_cmd": f"tar -xzf {{}} -C {utilities_dir}",
            "bin_path": os.path.join(os.path.abspath(utilities_dir), "ADFRsuite_x86_64Linux_1.0/bin") 
        },
        "click": {
            "path": os.path.join(utilities_dir, "Click.tar.gz"),
            "id_key": "click_id",
            "extract_cmd": f"tar -xzf {{}} -C {utilities_dir}",
            "bin_path": os.path.join(os.path.abspath(utilities_dir), "Click/bin")
        },
        "db": {
            "path": os.path.join(db_dir, "minipockets_surface80_winsize3_size3_curated-db.tar.gz"),
            "id_key": "db_id",
            "extract_cmd": f"tar -xzf {{}} -C {db_dir}"
        },
        "dict": {
            "path": os.path.join(db_dir, "reduce_wwPDB_het_dict.tar.gz"), 
            "id_key": "dict_id",
            "extract_cmd": f"tar -xzf {{}} -C {db_dir}"
        }
    }

    for name, info in files.items():
        if not os.path.exists(info["path"]):
            # Check if we have an ID to download (fallback)
            file_id = drive_ids.get(info["id_key"])
            if file_id:
                print(f"Downloading {name} (Fallback)...")
                url = f'https://drive.google.com/uc?id={file_id}'
                gdown.download(url, info["path"], quiet=False)
        
        # Extract if exists
        should_extract = False
        if os.path.exists(info["path"]):
             should_extract = True
             if "bin_path" in info and os.path.exists(info["bin_path"]):
                 should_extract = False
             # For DB files that don't have bin_path, we might re-extract unnecessarily?
             # Check if destination exists
             if name == "db" and os.path.exists(os.path.join(db_dir, "minipockets_surface80_winsize3_size3_curated-db")):
                 should_extract = False
             if name == "dict" and os.path.exists(os.path.join(db_dir, "reduce_wwPDB_het_dict.txt")):
                 should_extract = False
        
        if should_extract:
            print(f"Extracting {name}...")
            subprocess.run(info["extract_cmd"].format(info["path"]), shell=True, check=True)
            
        # Add to PATH if needed
        if "bin_path" in info and os.path.exists(info["bin_path"]):
            if info['bin_path'] not in os.environ['PATH']:
                os.environ['PATH'] += f":{info['bin_path']}"
                print(f"Added {name} to PATH: {info['bin_path']}")
            if name == "click":
                    subprocess.run(f"chmod +x {info['bin_path']}/click", shell=True)
    
    # Handle dictionary txt
    dict_txt = os.path.join(db_dir, "reduce_wwPDB_het_dict.txt")
    if os.path.exists(dict_txt):
        print("Dictionary txt found.")
    else:
        print("WARNING: reduce_wwPDB_het_dict.txt not found (maybe inside another folder after extraction?)")

