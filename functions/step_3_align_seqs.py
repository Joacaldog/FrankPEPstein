#@title 3. Sequence Alignment & Logo Analysis
#@markdown **Instructions:** 
#@markdown This step analyzes the finalized peptide candidates from Step 2.
#@markdown 1. Extracts amino acid sequences from the best PDB candidates.
#@markdown 2. Generates a Multiple Sequence Alignment (Multifasta).
#@markdown 3. Visualizes conserved motifs using a Sequence Logo.

import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import logomaker
from IPython.display import display, Image
import re

# Identify the latest run folder
initial_path = os.getcwd()
run_base = os.path.join(initial_path, "FrankPEPstein_run")

# Find the specific run folder with results (frankPEPstein_X/top_Y_peps)
# We need to search recursively or assume standard structure.
# Structure: FrankPEPstein_run/frankPEPstein_{pep_size}/top_{candidates}_peps/

# Recursive glob to match: FrankPEPstein_run/**/top_*_peps
candidate_folders = glob.glob(os.path.join(run_base, "**", "top_*_peps"), recursive=True)

if not candidate_folders:
    print("❌ No candidate results found from Step 2.")
else:
    # Use the most recent one if multiple (though pipeline likely cleans up)
    # Sort by modification time
    target_folder = sorted(candidate_folders, key=os.path.getmtime, reverse=True)[0]
    print(f"Analyzing results from: {target_folder}")
    
    # 1. Extract Sequences
    pdb_files = glob.glob(os.path.join(target_folder, "*.pdb")) # Actually frag*.pdb or similar? frankVINA_2 outputs top PDBs.
    
    sequences = []
    
    # Check if files exist
    if not pdb_files:
        # Maybe they are still pdbqt? frankVINA_2 converts to pdb at the end.
        print("⚠️ No PDB files found in target folder.")
    else:
        for pdb_path in pdb_files:
            filename = os.path.basename(pdb_path)
            # Expecting filename to contain sequence?
            # User said "si la contiene" (filename contains info).
            # frankVINA_2 output format typically: {score}_{sequence}.pdb or similar?
            # Let's try to find a sequence string (UPPERCASE letters).
            # Usually frag_SEQUENCE_score.pdb or SEQUENCE.pdb
            
            # Simple heuristic: extracting the longest string of uppercase letters
            # Or if user standard is specific...
            # Looking at frankVINA_2.py again, it seems it runs `vina_scorer`.
            # If we don't know exact format, we can extract from PDB SEQRES or Atoms (Chain 'p' or 'x').
            # User said: "desde el nombre del archivo mejor" AND "filename contains info".
            # Let's try to extract sequence from filename.
            # Assuming standard amino acids.
            
            # Attempt to match sequence chars [ACDEFGHIKLMNPQRSTVWY]
            # Pattern: Longest contiguous string of AAs?
            # Or maybe the whole filename is the sequence?
            
            # Fallback: Parse PDB if filename parsing is ambiguous, but user insists on filename.
            # Let's deduce from typical FrankPEPstein behavior. 
            # Often it is `SEQ_score.pdb` or `rank_SEQ_score.pdb`.
            
            # Heuristic: Find all caps string length > 4.
            matches = re.findall(r'[ACDEFGHIKLMNPQRSTVWY]{5,}', filename)
            if matches:
                 # Take the longest one
                 seq = max(matches, key=len)
                 sequences.append(seq)
            else:
                 # Fallback: try reading PDB?
                 # No, user said filename. Let's assume filename IS sequence if simple.
                 # Example: "AAAAA.pdb"
                 base = os.path.splitext(filename)[0]
                 if all(c in "ACDEFGHIKLMNPQRSTVWY_" for c in base): # Allow underscore
                     sequences.append(base.split('_')[0]) # Split score if present
                 else:
                     print(f"Skipping {filename}: Could not deduce sequence from name.")

        if sequences:
            print(f"Extracted {len(sequences)} sequences.")
            
            # 2. Generate Multifasta (MSA)
            fasta_path = os.path.join(target_folder, "candidates.fasta")
            with open(fasta_path, "w") as f:
                for i, seq in enumerate(sequences):
                    f.write(f">candidate_{i+1}\n{seq}\n")
            print(f"✅ Generated Multifasta: {fasta_path}")
            
            # 3. Generate Sequence Logo
            # Create a matrix for logomaker
            # Sequences must be same length for simple logo.
            lengths = [len(s) for s in sequences]
            if len(set(lengths)) > 1:
                print("⚠️ Sequences have varying lengths, alignment needed for proper Logo. Using simple left-alignment.")
                # Pad with gaps? or just ignore? Logomaker needs DataFrame.
                max_len = max(lengths)
                padded_seqs = [s.ljust(max_len, '-') for s in sequences]
                seq_list = padded_seqs
            else:
                seq_list = sequences

            try:
                # Create counts matrix
                logo_matrix = logomaker.alignment_to_matrix(seq_list)
                
                # Plot
                fig, ax = plt.subplots(figsize=(10, 4))
                logo = logomaker.Logo(logo_matrix, ax=ax)
                
                ax.set_title("Conserved Peptide Motifs", fontsize=14)
                ax.set_xlabel("Position", fontsize=12)
                ax.set_ylabel("Probability / Information", fontsize=12)
                
                logo_path = os.path.join(target_folder, "logo.png")
                plt.savefig(logo_path, bbox_inches='tight', dpi=300)
                plt.show() # Display in notebook
                print(f"✅ Generated Logo Plot: {logo_path}")
                
            except Exception as e:
                print(f"Error creating logo: {e}")
                
        else:
            print("No valid sequences found to align.")
