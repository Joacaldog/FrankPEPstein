#@title 3. Sequence Alignment & Logo Analysis
#@markdown **Instructions:** 
#@markdown This step analyzes the finalized peptide candidates from Step 2.
#@markdown 1. Extracts amino acid sequences from the best PDB candidates.
#@markdown 2. Generates a Multiple Sequence Alignment (Multifasta).
#@markdown 3. Visualizes conserved motifs using a Sequence Logo.

import os
import glob
import datetime
from collections import Counter
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
    # Pick the most recent, but say so out loud. A run that produced no candidates leaves the
    # previous run's folder in place, and silently analysing it looks exactly like success.
    ranked = sorted(candidate_folders, key=os.path.getmtime, reverse=True)
    target_folder = ranked[0]
    stamp = datetime.datetime.fromtimestamp(os.path.getmtime(target_folder)).strftime("%Y-%m-%d %H:%M:%S")
    print(f"Analyzing results from: {target_folder}  (last modified {stamp})")
    if len(ranked) > 1:
        print(f"⚠️ {len(ranked)} result folders exist. If the run you just launched failed, these "
              f"may be from an earlier run — check the timestamp above. Others found:")
        for other in ranked[1:]:
            other_stamp = datetime.datetime.fromtimestamp(os.path.getmtime(other)).strftime("%Y-%m-%d %H:%M:%S")
            print(f"    {other}  ({other_stamp})")
    
    # 1. Extract Sequences
    pdb_files = glob.glob(os.path.join(target_folder, "*.pdb")) # Actually frag*.pdb or similar? frankVINA_2 outputs top PDBs.
    
    sequences = []
    
    # Check if files exist
    if not pdb_files:
        # Maybe they are still pdbqt? frankVINA_2 converts to pdb at the end.
        print("⚠️ No PDB files found in target folder.")
    else:
        # The peptide sequence is embedded in the filename by patch_clustering's kmerizator
        # (frag_<SEQ>.pdb), which frankVINA_2 then prefixes and suffixes into
        # MinPEP_min_frag_<SEQ>_out.pdb. Anchor on that structure instead of scanning for the
        # longest run of amino-acid letters: the generic scan silently dropped peptides shorter
        # than five residues, and its fallback branch could never fire because the prefixes
        # contain lowercase characters.
        AA = "ACDEFGHIKLMNPQRSTVWY"
        patterns = [
            rf'frag_([{AA}]+)_out',     # MinPEP_min_frag_SEQ_out.pdb
            rf'frag_([{AA}]+)\.',       # frag_SEQ.pdb
            rf'^([{AA}]+)\.',           # SEQ.pdb
        ]
        for pdb_path in sorted(pdb_files):
            filename = os.path.basename(pdb_path)
            seq = None
            for pattern in patterns:
                found = re.search(pattern, filename)
                if found:
                    seq = found.group(1)
                    break
            if seq:
                sequences.append(seq)
            else:
                print(f"Skipping {filename}: could not read a sequence from the name.")

        if sequences:
            print(f"Extracted {len(sequences)} sequences.")
            
            # 2. Generate Multifasta (MSA)
            fasta_path = os.path.join(target_folder, "candidates.fasta")
            with open(fasta_path, "w") as f:
                for i, seq in enumerate(sequences):
                    f.write(f">candidate_{i+1}\n{seq}\n")
            print(f"✅ Generated Multifasta: {fasta_path}")
            
            # 3. Generate Sequence Logo
            # Only equal-length sequences go into the logo. These peptides are independently
            # docked, not homologous, so there is nothing to align: right-padding shorter ones
            # with gaps pairs non-corresponding positions, and because logomaker ignores gap
            # characters the padded columns are built from fewer observations. The result is
            # letters that shrink toward the C-terminus, which reads as loss of conservation
            # when it is only loss of sample size.
            lengths = [len(s) for s in sequences]
            modal_length = Counter(lengths).most_common(1)[0][0]
            seq_list = [s for s in sequences if len(s) == modal_length]
            if len(seq_list) < len(sequences):
                excluded = len(sequences) - len(seq_list)
                print(f"⚠️ Logo built from the {len(seq_list)} sequences of length {modal_length}; "
                      f"excluded {excluded} of other lengths (they remain in the FASTA). "
                      f"These candidates are independently docked, not homologous, so padding "
                      f"them to a common length would not be an alignment.")

            try:
                # information, not raw counts: the y-axis says bits and must mean bits.
                logo_matrix = logomaker.alignment_to_matrix(seq_list, to_type='information')

                # Plot
                fig, ax = plt.subplots(figsize=(10, 4))
                logo = logomaker.Logo(logo_matrix, ax=ax)

                ax.set_title(f"Residue frequencies across {len(seq_list)} candidates "
                             f"(length {modal_length})", fontsize=14)
                ax.set_xlabel("Position", fontsize=12)
                ax.set_ylabel("Information (bits)", fontsize=12)

                logo_path = os.path.join(target_folder, "logo.png")
                plt.savefig(logo_path, bbox_inches='tight', dpi=300)
                plt.show() # Display in notebook
                print(f"✅ Generated Logo Plot: {logo_path}")
                if len(seq_list) < 20:
                    print(f"⚠️ Only {len(seq_list)} sequences: treat apparent conservation with "
                          f"caution, a small sample looks conserved by chance.")

            except Exception as e:
                print(f"Error creating logo: {e}")
                
        else:
            print("No valid sequences found to align.")
