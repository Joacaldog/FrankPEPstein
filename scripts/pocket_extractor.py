#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Added arg parse

import argparse
import numpy as np
import pandas as pd
import os
import sys
import fnmatch
from collections import Counter
from tqdm import tqdm
from Bio.PDB import *
import glob
import itertools
import warnings
from joblib import Parallel, delayed

warnings.simplefilter('ignore')

program_description = "generates pocket files"
parser = argparse.ArgumentParser(description=program_description)
parser.add_argument("-s", "--minimum_surface", type=int,
                    help="Minimum contact surface for receptor residues", required=True)
parser.add_argument("-t", "--threads", type=int,
                    help="Number of threads", required=True)

args = parser.parse_args()
above = args.minimum_surface
threads = args.threads


def get_res_contact_area(above, folder):
    try:
        warnings.simplefilter('ignore')

        class SelectFragOnly(Select):
            def accept_residue(self, residue):
                if (residue in recRes_pocket):
                    return True  # Accept residue
                else:
                    return False  # Discard residue

        class SelectFragOnly2(Select):
            def accept_residue(self, residue):
                if (residue in ligRes_lumen):
                    return True  # Accept residue
                else:
                    return False  # Discard residue

        pocket_list = []
        folder_file = f'{folder}/{folder}.matrix.{"".join(sorted(folder.split("_")[1:]))}.by_res.tsv' 
        with open(folder_file) as f:
            l = f.readline()
            columns = l.strip().split("\t")
        df = pd.read_csv(folder_file, sep="\t", index_col=0,
                            skipinitialspace=True, skip_blank_lines=True, usecols=columns)
        columns = df.columns
        df = df.drop(
            [col for col in columns if df[col].sum() <= above], axis=1)
        columns = df.columns
        indexes = df.index

        for i in range(len(indexes)):
            for c in range(len(columns)):
                area_of_res = df.iloc[i, c]
                peptide_res = df.index[i]
                receptor_residue = df.columns[c]
                chain_id = folder.split("_")[-1]
                chain = "/" + chain_id + "/"
                if area_of_res >= above:
                    if chain in peptide_res:
                        pocket_list.append(receptor_residue)

        folder_file = folder + "/peptide_complex.pdb"
        recRes_pocket = []
        for recRes in pocket_list:
            parser = PDBParser()
            structure = parser.get_structure(
                folder + "1", folder_file)
            io = PDBIO()
            io.set_structure(structure)
            chain_receptor = recRes.split("/")[1]
            res_id = recRes.split("/")[2]
            res_name = recRes.split("/")[0]
            for model in structure:
                for chain in model:
                    if chain.get_id() == chain_receptor:
                        chain.id = "p"
                        for residue in chain:
                            if str(residue.get_id()[1]) == res_id:
                                recRes_pocket.append(residue)
            outfile_name = ("{f}/pocket.pdb").format(f=folder)
            io.save(outfile_name, select=SelectFragOnly())
    except Exception as e:
        print(f"Error in {folder}: {e}")
Parallel(n_jobs=threads)(delayed(get_res_contact_area)(above, folder) for folder in tqdm(next(os.walk('.'))[1], total=len(
    next(os.walk('.'))[1]), desc="processing matrix files and generating lumen and pocket files"))
