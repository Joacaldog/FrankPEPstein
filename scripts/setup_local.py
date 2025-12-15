#!/usr/bin/env python3
import os
import sys
import shutil
import subprocess
import argparse

# Colors for output
GREEN = "\033[92m"
RED = "\033[91m"
YELLOW = "\033[93m"
RESET = "\033[0m"

def log(msg, color=RESET):
    print(f"{color}{msg}{RESET}")

def check_command(cmd):
    return shutil.which(cmd) is not None

def get_conda_cmd():
    if check_command("mamba"):
        return "mamba"
    elif check_command("conda"):
        return "conda"
    return None

def main():
    parser = argparse.ArgumentParser(description="FrankPEPstein Local Setup")
    # Key is now hardcoded default
    modeller_key = "MODELIRANJE"
    
    args = parser.parse_args()
    
    repo_root = os.getcwd()
    if not os.path.exists(os.path.join(repo_root, "scripts", "setup_local.py")):
        log("❌ Please run this script from the repository root.", RED)
        sys.exit(1)

    log(f"--- FrankPEPstein Local Setup ---", GREEN)
    
    # 1. Check Conda
    conda_cmd = get_conda_cmd()
    if not conda_cmd:
        log("❌ Conda or Mamba not found. Please install Mambaforge or Miniconda.", RED)
        sys.exit(1)
    log(f"✅ Found Package Manager: {conda_cmd}", GREEN)

    # 2. Create Environment
    env_name = "FrankPEPstein"
    check_env = subprocess.run(f"{conda_cmd} env list", shell=True, capture_output=True, text=True)
    
    if env_name in check_env.stdout:
        log(f"⚠️ Environment '{env_name}' already exists. Skipping creation.", YELLOW)
        log("To recreate, delete it first: conda env remove -n FrankPEPstein", YELLOW)
    else:
        log(f"Creating Conda Environment '{env_name}'... (This may take a while)", parse_color=GREEN)
        cmd = (
            f"{conda_cmd} create -n {env_name} -y -c conda-forge -c salilab "
            "openbabel biopython fpocket joblib tqdm py3dmol vina pigz scipy scikit-learn matplotlib python=3.10 salilab::modeller"
        )
        log(f"Running: {cmd}")
        ret = subprocess.run(cmd, shell=True)
        if ret.returncode != 0:
            log("❌ Environment creation failed.", RED)
            sys.exit(1)
        log("✅ Environment Created.", GREEN)

    # 3. Download Databases & Tools
    files_id = "1M30wmaf6vaXJl1kmj-0cD5yhBYDCx_xw"
    archive_name = "files.tar.gz"
    
    if os.path.exists("DB") and os.path.exists("utilities"):
        log("✅ DB and utilities folders already exist. Skipping download.", GREEN)
    else:
        log("Downloading Databases & Tools...", GREEN)
        
        # Check if archive exists
        if not os.path.exists(archive_name):
            # Try gdown first
            if check_command("gdown"):
                subprocess.run(f"gdown --id {files_id} -O {archive_name}", shell=True)
            else:
                log("⚠️ gdown not found. Attempting to install in current python...", YELLOW)
                subprocess.run(f"{sys.executable} -m pip install gdown", shell=True)
                subprocess.run(f"gdown --id {files_id} -O {archive_name}", shell=True)
        
        if os.path.exists(archive_name):
            log("Extracting archive...", GREEN)
            # Try pigz for speed
            if check_command("pigz"):
                subprocess.run(f"tar -I pigz -xf {archive_name}", shell=True)
            else:
                subprocess.run(f"tar -xf {archive_name}", shell=True)
            log("✅ Extraction complete.", GREEN)
        else:
            log("❌ Failed to download files.tar.gz", RED)
            # Don't exit, might be partially setup manually? 
    
    # 4. Utilities Setup
    log("Configuring Utilities...", GREEN)
    utilities_dir = os.path.join(repo_root, "utilities")
    
    if os.path.exists(utilities_dir):
        # Rename ADFR
        adfr_extracted = os.path.join(utilities_dir, "ADFR")
        adfr_target = os.path.join(utilities_dir, "ADFRsuite_x86_64Linux_1.0")
        if os.path.exists(adfr_extracted) and not os.path.exists(adfr_target):
             log(f"Renaming ADFR folder...", YELLOW)
             os.rename(adfr_extracted, adfr_target)
             
        # [ADDED] Patch ADFR scripts with current path (fix hardcoded /home/joacaldo paths)
        if os.path.exists(adfr_target):
             log("Patching ADFR scripts paths...", GREEN)
             bin_dir = os.path.join(adfr_target, "bin")
             abs_adfr_path = os.path.abspath(adfr_target)
             for fname in os.listdir(bin_dir):
                 fpath = os.path.join(bin_dir, fname)
                 if os.path.isfile(fpath) and not os.path.islink(fpath):
                     try:
                         # Read (ignore errors for binaries)
                         with open(fpath, 'rb') as f:
                             content_bytes = f.read()
                         
                         # Check if text file (shebang or ADS_ROOT)
                         try:
                             content = content_bytes.decode('utf-8')
                             if "ADS_ROOT=" in content:
                                 import re
                                 # Replace ADS_ROOT="..." with correct path
                                 new_content = re.sub(r'ADS_ROOT="[^"]+"', f'ADS_ROOT="{abs_adfr_path}"', content)
                                 with open(fpath, 'w') as f:
                                     f.write(new_content)
                                 log(f"Patched {fname}", YELLOW)
                         except UnicodeDecodeError:
                             pass # Binary file
                     except Exception as e:
                         pass
             
        # Permissions
        log("Fixing permissions...", GREEN)
        subprocess.run(f"chmod -R +x {utilities_dir}", shell=True)
        
        # Config Local 'scripts/notebook_utils.py' is NOT needed locally 
        # because Step 0 generated it for Colab environment patching.
        # But we do need Modeller Config
        
        # 5. Configure Modeller
        log("Configuring Modeller License...", GREEN)
        # We must run the configuration script INSIDE the conda environment
        # to find the correct config.py file in that env.
        
        config_script = os.path.join(repo_root, "scripts", "configure_modeller.py")
        
        # Using 'conda run' to execute inside the env
        modeller_cmd = f"{conda_cmd} run -n {env_name} python {config_script}"
        # configure_modeller.py currently doesn't take args, it has a default. 
        # But wait, looking at configure_modeller.py source...
        # It calls `configure_modeller(license_key='MODELIRANJE')` in main block.
        # It DOES NOT parse argparse.
        # I should simply edit the file OR rely on user editing it?
        # Better: Pass it via environment variable or just assume user changes it later?
        # Setup said "replicate logic". In colab it patches with 'MODELIRANJE'.
        # I will assume the script uses the default 'MODELIRANJE' and user changes it later if they have a real key.
        # Or I can update configure_modeller to take an arg.
        
        # For now, let's just run it as is.
        modeller_cmd = f"{conda_cmd} run -n {env_name} python {config_script}"
        
        ret = subprocess.run(modeller_cmd, shell=True)
        if ret.returncode == 0:
            log("✅ Modeller Configured.", GREEN)
        else:
            log("⚠️ Modeller configuration threw an error (maybe env not ready?). Check manually.", YELLOW)

    log(f"\n{'='*30}", GREEN)
    log("✅ Setup Finished!", GREEN)
    log(f"To run the pipeline:", GREEN)
    log(f"1. conda activate {env_name}", GREEN)
    log(f"2. python scripts/run_local.py ...", GREEN)
    log(f"{'='*30}\n", GREEN)

if __name__ == "__main__":
    main()
