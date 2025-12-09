import os
import shutil

def run_results():
    print("--- Cell 06: Results & Download ---")
    
    results_dir = "results_folder" # Placeholder name from notebook
    # Check if results exist (mocking check)
    # real path would be run_frankpepstein/results or similar?
    # Notebook implies 'results_folder/' exists in CWD.
    
    if not os.path.exists(results_dir):
        print(f"Warning: {results_dir} not found (Pipeline hasn't run or failed).")
        print("Creating dummy results for test...")
        os.makedirs(results_dir, exist_ok=True)
        with open(f"{results_dir}/result_summary.txt", "w") as f:
            f.write("Dummy results")
            
    print(f"Visualizing results in {results_dir}...")
    files = os.listdir(results_dir)
    print(f"Generated {len(files)} result files.")
    
    archive_name = "frankpepstein_results"
    shutil.make_archive(archive_name, 'zip', results_dir)
    print(f"Results compressed to: {os.path.abspath(archive_name)}.zip")

if __name__ == "__main__":
    run_results()
