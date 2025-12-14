
#@title 0. Install CondaColab & Setup Tools (~3 min)
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
        ("Setting up External Tools (DB Download)", "tools"),
        ("Configuring Modeller", "modeller")
    ]
    
    with tqdm(total=len(steps)) as pbar:
        # 1. CondaColab
        pbar.set_description(steps[0][0])
        try:
            with SuppressStdout():
                import condacolab
                condacolab.check()
        except (ImportError, AssertionError):
            print("Installing CondaColab... (Kernel will restart and colab will say it crashes, you dont need to do anything)")
            print("Installing python dependencies...")
            subprocess.run("pip install -q py3dmol logomaker", shell=True, check=True)
            subprocess.run("pip install -q biopython", shell=True, check=True) # Ensure biopython is there too
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



def setup_external_tools(files_id=None):
    # Install gdown if needed
    try: import gdown
    except ImportError: subprocess.run([sys.executable, "-m", "pip", "install", "gdown"], check=True); import gdown

    # Determine Base Dir
    base_dir = "." # Always download to current directory (root of colab)
    
    # 1. Download Single Archive
    if files_id:
        archive_path = os.path.join(base_dir, "files.tar.gz")
        # Check if already extracted (simple check for DB and utilities)
        if not (os.path.exists(os.path.join(base_dir, "DB")) and os.path.exists(os.path.join(base_dir, "utilities"))):
            if not os.path.exists(archive_path):
                 print(f"Downloading files.tar.gz (ID: {files_id})...")
                 gdown.download(f'https://drive.google.com/uc?id={files_id}', archive_path, quiet=True)
            
            if os.path.exists(archive_path):
                 print("Extracting files.tar.gz (Parallel with pigz)...")
                 # pigz -d -c files.tar.gz | tar xf -
                 # or tar -I pigz -xf ...
                 subprocess.run(f"tar -I pigz -xf {archive_path} -C {base_dir}", shell=True, check=True)
                 # os.remove(archive_path) # Optional clean up
    
    # 2. Permissions & Final Setup
    utilities_dir = os.path.join(base_dir, "utilities")
    if os.path.exists(utilities_dir):
        # [ADDED] Rename extracted ADFR folder to Suite name required by tools
        adfr_extracted = os.path.join(utilities_dir, "ADFR")
        adfr_target = os.path.join(utilities_dir, "ADFRsuite_x86_64Linux_1.0")
        if os.path.exists(adfr_extracted) and not os.path.exists(adfr_target):
             print(f"Renaming {adfr_extracted} to {adfr_target}...")
             os.rename(adfr_extracted, adfr_target)
        
        # [ADDED] Patch ADFR scripts with current path (fix hardcoded /home/joacaldo paths)
        if os.path.exists(adfr_target):
             print("Patching ADFR scripts paths...")
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
                         except UnicodeDecodeError:
                             pass # Binary file
                     except Exception as e:
                         pass
             
         # if os.path.exists(utilities_dir): # REMOVED redundant line causing indent error
        print("Fixing permissions...")
        subprocess.run(f"chmod -R +x {utilities_dir}", shell=True)
        
        # Add paths
        adfr_bin = os.path.join(utilities_dir, "ADFR/bin")
        if os.path.exists(adfr_bin) and adfr_bin not in os.environ['PATH']:
             os.environ['PATH'] += f":{adfr_bin}"
             
        click_bin = os.path.join(utilities_dir, "click") # assuming it extracts directly or in bin
        # Adjust check if structure is different (e.g. utilities/Click/bin)
        # User said "utilities que quedaran en ~/. Ahora ADFR al extraerse queda listo... lo mismo con click"
        # We assume standard structure.
        

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
        files_id = "1M30wmaf6vaXJl1kmj-0cD5yhBYDCx_xw"
        
        with SuppressStdout():
             notebook_utils.setup_external_tools(files_id)
        pbar.update(1)

        # 6. Configure Modeller
        pbar.set_description(steps[5][0])
        with SuppressStdout():
            notebook_utils.configure_modeller()
        pbar.update(1)
        

    # 3. Verify Executables
    print(f"\n{'='*20}")
    print("Verifying Executables...")
    
    base_dir = "FrankPEPstein" if os.path.exists("FrankPEPstein") else "."
    utilities_dir = os.path.join(base_dir, "utilities")
    
    # Check Click
    click_bin = os.path.join(utilities_dir, "Click", "bin", "click")
    if not os.path.exists(click_bin):
         # Try logic from superposer logic/config
         click_bin = os.path.join(utilities_dir, "Click", "click")
         
    if os.path.exists(click_bin):
         if os.access(click_bin, os.X_OK):
             print(f"✅ Click is executable.")
         else:
             print(f"⚠️ Click found but NOT executable. Fixing...")
             subprocess.run(f"chmod +x {click_bin}", shell=True)
             if os.access(click_bin, os.X_OK):
                 print(f"✅ Click fixed.")
             else:
                 print(f"❌ Failed to fix Click permissions.")
    else:
         print(f"❌ Click binary not found at {click_bin}")

    # Check Vina
    vina_bin = os.path.join(utilities_dir, "vina_1.2.4_linux_x86_64")
    if os.path.exists(vina_bin):
         if os.access(vina_bin, os.X_OK):
             print(f"✅ Vina is executable.")
         else:
             print(f"⚠️ Vina found but NOT executable. Fixing...")
             subprocess.run(f"chmod +x {vina_bin}", shell=True)
             if os.access(vina_bin, os.X_OK):
                 print(f"✅ Vina fixed.")
             else:
                 print(f"❌ Failed to fix Vina permissions.")
    else:
         print(f"❌ Vina binary not found at {vina_bin}")
    print(f"{'='*20}\n")

    clear_output()
    print("✅ Setup Ready!")

if __name__ == "__main__":
    run_setup()

