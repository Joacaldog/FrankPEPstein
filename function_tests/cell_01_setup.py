import sys
import os

# Add repo root to path so we can import scripts
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from scripts import notebook_utils

def run_setup():
    print("--- Cell 01: Setup & Dependencies ---")
    
    # IDs provided by user
    drive_ids = {
        "adfr_id": "1gmRj8mva84-JB7UXUcQfB3Ziw_nwwdox",
        "db_id": "1a4GoZ1ZT-DNYMyvVtKJukNdF6TAaLJU5", 
        "dict_id": "1nrwSUof0lox9fp8Ow5EICIN9u0lglu7U"
    }
    
    print("Setting up external tools...")
    notebook_utils.setup_external_tools(drive_ids)
    
    print("Configuring Modeller...")
    notebook_utils.configure_modeller()
    print("Setup Complete.")

if __name__ == "__main__":
    run_setup()
