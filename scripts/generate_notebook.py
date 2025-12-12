import os
import json
import glob
import re

def create_notebook():
    functions_dir = "functions"
    output_notebook = "FrankPEPstein.ipynb"
    
    # Notebook Structure (nbformat 4)
    notebook = {
        "cells": [],
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3"
            },
            "language_info": {
                "codemirror_mode": {
                    "name": "ipython",
                    "version": 3
                },
                "file_extension": ".py",
                "mimetype": "text/x-python",
                "name": "python",
                "nbconvert_exporter": "python",
                "pygments_lexer": "ipython3",
                "version": "3.10.12"
            },
            "colab": {
                "provenance": []
            }
        },
        "nbformat": 4,
        "nbformat_minor": 0
    }

    # Get python files sorted
    # We want step_0, step_1, step_2 in order.
    # Glob might return random order, so we sort.
    files = sorted(glob.glob(os.path.join(functions_dir, "step_*.py")))

    if not files:
        print(f"No files found in {functions_dir}")
        return

    for script_path in files:
        print(f"Processing {script_path}...")
        with open(script_path, "r") as f:
            content = f.read()
        
        # Extract Title from first line if it looks like #@title ...
        # (Colab form cell)
        # We put the whole content in one cell.
        
        cell = {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {
                "id": os.path.basename(script_path).replace(".py", "") 
            },
            "outputs": [],
            "source": content.splitlines(keepends=True)
        }
        
        # If the file has #@title, we might want to ensure it is detectable by Colab
        # It usually is just the first line.
        
        notebook["cells"].append(cell)

    with open(output_notebook, "w") as f:
        json.dump(notebook, f, indent=2)
    
    print(f"✅ Generated {output_notebook} with {len(files)} cells.")

if __name__ == "__main__":
    # Ensure we are in repo root or relative correct path
    if not os.path.exists("functions"):
        # Try to find it relative to script location
        script_dir = os.path.dirname(os.path.abspath(__file__))
        repo_root = os.path.dirname(script_dir)
        if os.path.exists(os.path.join(repo_root, "functions")):
            os.chdir(repo_root)
            print(f"Changed directory to {repo_root}")
    
    create_notebook()
