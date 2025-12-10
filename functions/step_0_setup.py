
#@title 0. Install CondaColab & Setup Tools (~2 min)
import sys
import os
import subprocess
from IPython.display import clear_output

# Helper to suppress output
class SuppressStdout:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')
    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

def run_setup():
    # Install tqdm first if missing (fast)
    try:
        from tqdm.notebook import tqdm
    except ImportError:
        subprocess.run("pip install -q tqdm", shell=True)
        from tqdm.notebook import tqdm

    print("Setting up FrankPEPstein environment...")
    
    steps = [
        ("Installing CondaColab", "condacolab"),
        ("Cloning Repository", "git"),
        ("Creating Conda Environment (Slow)", "env"),
        ("Configuring Notebook Utils", "patch"),
        ("Setting up External Tools (Parallel DB Download)", "tools"),
        ("Configuring Modeller", "modeller")
    ]
    
    with tqdm(total=len(steps)) as pbar:
        # 1. CondaColab
        pbar.set_description(steps[0][0])
        try:
            with SuppressStdout():
                import condacolab
                condacolab.check()
        except ImportError:
            with SuppressStdout():
                subprocess.run("pip install -q condacolab", shell=True, check=True)
                import condacolab
                condacolab.install()
        pbar.update(1)

        # 2. Git Clone
        pbar.set_description(steps[1][0])
        with SuppressStdout():
            if not os.path.exists("FrankPEPstein"):
                subprocess.run("git clone https://github.com/Joacaldog/FrankPEPstein.git", shell=True, check=True)
        pbar.update(1)

        # 3. Create Environment
        pbar.set_description(steps[2][0])
        env_path = "/usr/local/envs/FrankPEPstein"
        if not os.path.exists(env_path):
             # Added 'pigz' to the package list
             subprocess.run("mamba create -n FrankPEPstein -q -y -c conda-forge -c salilab openbabel biopython fpocket joblib tqdm py3dmol vina pigz python=3.10 salilab::modeller > /dev/null 2>&1", shell=True, check=True)
        
        # Configure Path
        site_packages = f"{env_path}/lib/python3.10/site-packages"
        if site_packages not in sys.path:
            sys.path.append(site_packages)
        os.environ['PATH'] = f"{env_path}/bin:" + os.environ['PATH']
        pbar.update(1)

        # 4. Patch Utils
        pbar.set_description(steps[3][0])
        # We write a clean utils file
        patched_utils_content = r'''import os
import sys
import glob
import shutil
import subprocess

def configure_modeller(license_key='MODELIRANJE', repo_dir='FrankPEPstein'):
    template_config = os.path.join(repo_dir, "utilities/config.py")
    dest_config = None
    try:
        import modeller
        modeller_path = os.path.dirname(modeller.__file__)
        candidate = os.path.join(modeller_path, "config.py")
        if os.path.exists(candidate):
            dest_config = candidate
    except Exception:
        pass

    if not dest_config:
        possible_paths = [
            f"{sys.prefix}/lib/modeller-*/modlib/modeller/config.py",
            f"{sys.prefix}/lib/python*/site-packages/modeller/config.py",
            "/usr/local/envs/FrankPEPstein/lib/modeller-*/modlib/modeller/config.py"
        ]
        for pattern in possible_paths:
            found = glob.glob(pattern)
            if found:
                dest_config = found[0]
                break

    if dest_config and os.path.exists(template_config):
        with open(template_config, 'r') as f:
            content = f.read()
        new_content = content.replace("'MODELIRANJE'", f"'{license_key}'")
        with open(dest_config, 'w') as f:
            f.write(new_content)
        return True
    return False

def setup_external_tools(drive_ids=None):
    if drive_ids is None: drive_ids = {}
    
    # Install gdown if needed
    try: import gdown
    except ImportError: subprocess.run([sys.executable, "-m", "pip", "install", "gdown"], check=True); import gdown

    # Determine Base Dir
    base_dir = "FrankPEPstein" if os.path.exists("FrankPEPstein") else "."
    utilities_dir = os.path.join(base_dir, "utilities")
    db_dir = os.path.join(base_dir, "DB")
    os.makedirs(utilities_dir, exist_ok=True)
    os.makedirs(db_dir, exist_ok=True)


    # 1. Database (Parallel Download & Decompress)
    db_id = drive_ids.get("db_id")
    # Determine pigz availability
    use_pigz = True if subprocess.run("which pigz", shell=True).returncode == 0 else False
    extract_flag = "-I pigz -xf" if use_pigz else "-xzf"

    db_tar = os.path.join(base_dir, "filtered_DB_P5-15_R30_id10_optim.tar.gz")
    
    # Check if Main DB is already populated (simple check)
    if not os.path.exists(os.path.join(db_dir, "filtered_DB_P5-15_R30_id10")):
        if db_id:
             if not os.path.exists(db_tar):
                 print(f"Downloading DB (ID: {db_id})...")
                 # We assume single big tarball for now based on user context
                 gdown.download(f'https://drive.google.com/uc?id={db_id}', db_tar, quiet=True)
             
             if os.path.exists(db_tar):
                 print(f"Extracting DB using {extract_flag}...")
                 subprocess.run(f"tar {extract_flag} {db_tar} -C {db_dir} > /dev/null 2>&1", shell=True, check=True)
                 # Optional: delete tar to save space?
                 # os.remove(db_tar)

    # 1.5 Minipockets Database
    mini_id = drive_ids.get("minipockets_id")
    mini_tar = os.path.join(base_dir, "minipockets.tar.gz")
    # We check for a known folder inside, e.g. "minipockets_surface80_winsize3_size3_curated"
    # Or just check if the tar was extracted.
    # User said "-fm" points to it.
    # Let's assume the tar extracts to a folder inside DB.
    
    if mini_id:
        # Check standard folder name (guess or generic)
        # If user didn't specify name, we rely on tar content.
        # But we need a check to avoid redownload.
        # Let's check if 'minipockets*' exists in DB dir
        existing_mini = glob.glob(os.path.join(db_dir, "minipockets*"))
        if not existing_mini:
             if not os.path.exists(mini_tar):
                  print(f"Downloading Minipockets DB (ID: {mini_id})...")
                  gdown.download(f'https://drive.google.com/uc?id={mini_id}', mini_tar, quiet=True)
             
             if os.path.exists(mini_tar):
                  print(f"Extracting Minipockets DB...")
                  subprocess.run(f"tar {extract_flag} {mini_tar} -C {db_dir} > /dev/null 2>&1", shell=True, check=True)

    # 2. Utilities (ADFR)
    # Usually passed as zip or tar
    # User might have "utilities_pkg_id"
    util_id = drive_ids.get("utilities_pkg_id")
    util_tar = os.path.join(base_dir, "utilities.tar.gz")
    
    if util_id and not os.path.exists(os.path.join(utilities_dir, "ADFRsuite_x86_64Linux_1.0")):
         if not os.path.exists(util_tar):
              gdown.download(f'https://drive.google.com/uc?id={util_id}', util_tar, quiet=True)
         subprocess.run(f"tar {extract_flag} {util_tar} -C {utilities_dir} > /dev/null 2>&1", shell=True, check=True)

    # ADFR Installation / License Check
    # Look for install script
    adfr_root = os.path.join(utilities_dir, "ADFRsuite_x86_64Linux_1.0")
    install_sh = os.path.join(utilities_dir, "install.sh") # Common location? Or inside extracted folder?
    
    # Sometimes it extracts AS "ADFRsuite..." directly.
    # If the user provides a zip of the INSTALLED directory, no installation needed.
    # Checks bin
    adfr_bin = os.path.join(adfr_root, "bin")
    if os.path.exists(adfr_bin):
        # Add to PATH
        if adfr_bin not in os.environ['PATH']:
            os.environ['PATH'] += f":{adfr_bin}"
            print(f"Added ADFR bin to PATH: {adfr_bin}")
    else:
        # If bin missing, maybe we need to run install?
        # Check for an install script in utilities
        candidates = glob.glob(os.path.join(utilities_dir, "*install*.sh"))
        if candidates:
             installer = candidates[0]
             print(f"Running ADFR installer: {installer}")
             # Pipe 'yes' to accept license
             # usage: yes | ./install.sh -d target_dir
             # or simply ./install.sh and expect yes
             # We assume it asks "Do you accept... (y/n)"
             subprocess.run(f"yes | bash {installer}", shell=True) # Try bash execution
             
             # Recheck bin
             if os.path.exists(adfr_bin):
                  os.environ['PATH'] += f":{adfr_bin}"

    # Click
    click_bin = os.path.join(utilities_dir, "Click/bin") 
    if os.path.exists(click_bin):
        if click_bin not in os.environ['PATH']:
             os.environ['PATH'] += f":{click_bin}"
        subprocess.run(f"chmod +x {click_bin}/click", shell=True)

'''
        os.makedirs("FrankPEPstein/scripts", exist_ok=True)
        with open("FrankPEPstein/scripts/notebook_utils.py", "w") as f:
            f.write(patched_utils_content)
        pbar.update(1)

        # 5. External Tools Setup
        pbar.set_description(steps[4][0])
        repo_path = os.path.abspath("FrankPEPstein")
        if repo_path not in sys.path:
            sys.path.append(repo_path)
        from scripts import notebook_utils
        
        # Correct Drive IDs
        drive_ids = {
            "db_id": "13a6M_UVham9SiBCE6PCQQi6CnUvGNvNO", # Provided by user
            "minipockets_id": "1a4GoZ1ZT-DNYMyvVtKJukNdF6TAaLJU5", # New Minipockets DB
            "utilities_pkg_id": "1gmRj8mva84-JB7UXUcQfB3Ziw_nwwdox", # Kept existing
            # If explicit ADFR installer needed:
            # "adfr_installer_id": "..."
        }
        
        with SuppressStdout():
             notebook_utils.setup_external_tools(drive_ids)
        pbar.update(1)

        # 6. Configure Modeller
        pbar.set_description(steps[5][0])
        with SuppressStdout():
            notebook_utils.configure_modeller()
        pbar.update(1)
        
    clear_output()
    print("✅ Setup Ready!")

if __name__ == "__main__":
    run_setup()

