#@title 4. Download Results
#@markdown **Instructions:**
#@markdown Click the button below to download a ZIP archive containing:
#@markdown - Candidate Peptide PDBs
#@markdown - Sequence Alignment (Fasta)
#@markdown - Sequence Motif Plot (Logo)

import os
import shutil
import glob
from google.colab import files
import ipywidgets as widgets
from IPython.display import display
from datetime import datetime

# Logic to find the target folder (same as Step 3)
initial_path = os.getcwd()
run_base = os.path.join(initial_path, "FrankPEPstein_run")


def find_target_folder():
    """Resolve the results folder at click time, not at cell-execution time.

    Binding this once when the cell runs meant the button kept pointing at whatever existed
    then: it zipped a previous run's results, or refused to do anything, after a later run had
    already produced new ones.
    """
    # Recursive glob to match: FrankPEPstein_run/**/top_*_peps
    candidate_folders = glob.glob(os.path.join(run_base, "**", "top_*_peps"), recursive=True)
    if not candidate_folders:
        return None
    return sorted(candidate_folders, key=os.path.getmtime, reverse=True)[0]


out_log = widgets.Output()

def download_results(b):
    out_log.clear_output()
    with out_log:
        target_folder = find_target_folder()
        if not target_folder or not os.path.exists(target_folder):
            print("❌ No results found to download.")
            return

        stamp = datetime.fromtimestamp(os.path.getmtime(target_folder)).strftime("%Y-%m-%d %H:%M:%S")
        print(f"Compressing results from: {target_folder}  (last modified {stamp})")
        
        # Timestamped filename
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        zip_filename = f"FrankPEPstein_Results_{timestamp}"
        zip_path = os.path.join(initial_path, zip_filename) # shutil.make_archive adds .zip extension automatically
        
        try:
            # Create ZIP
            shutil.make_archive(zip_path, 'zip', target_folder)
            final_zip = zip_path + ".zip"
            
            print(f"✅ Created archive: {final_zip}")
            print("Downloading...")
            
            files.download(final_zip)
            
        except Exception as e:
            print(f"Error during download: {e}")

btn_download = widgets.Button(
    description='Download Results (ZIP)',
    button_style='info',
    icon='download',
    layout=widgets.Layout(width='50%')
)
btn_download.on_click(download_results)

display(btn_download, out_log)
