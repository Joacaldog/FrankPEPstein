
import os
import math
import numpy as np
import matplotlib.pyplot as plt
from io import BytesIO

def get_atom_coords(pdb_file, atom_type=None):
    """
    Parses a PDB file and returns a list of (x, y, z) coordinates.
    If atom_type is specified (e.g., 'CA'), only those atoms are returned.
    """
    coords = []
    if not os.path.exists(pdb_file):
        return np.array(coords)
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Check atom type if requested
                name = line[12:16].strip()
                if atom_type and name != atom_type:
                    continue
                
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
                except ValueError:
                    pass
    return np.array(coords)

def align_to_principal_axes(coords):
    """
    Aligns coordinates to their principal axes using SVD.
    Returns aligned coordinates and the rotation matrix.
    """
    if coords.shape[0] < 3:
        return coords, np.eye(3)
        
    # Center
    center = np.mean(coords, axis=0)
    centered = coords - center
    
    # SVD
    # U, S, Vh = np.linalg.svd(centered)
    # Vh contains the rotation (principal axes)
    # Rotated coordinates = centered @ Vh.T
    
    # Covariance matrix
    cov = np.cov(centered, rowvar=False)
    evals, evecs = np.linalg.eigh(cov)
    
    # Sort eigenvalues/vectors (largest to smallest)
    idx = evals.argsort()[::-1]
    evecs = evecs[:, idx]
    
    # Rotate
    aligned = np.dot(centered, evecs)
    
    # Force X to be the long axis (already done by sort order usually)
    return aligned, evecs, center

def render_static_view(pocket_path, fragments_paths, title="Processing..."):
    """
    Generates a static matplotlib image (buffer) of the pocket and fragments.
    Aligns the view so the pocket is horizontal.
    """
    
    # 1. Load Pocket Coords
    pocket_coords = get_atom_coords(pocket_path)
    if len(pocket_coords) == 0:
        return None
        
    # 2. Align Pocket (Principal Axis to X)
    p_aligned, rot_matrix, center = align_to_principal_axes(pocket_coords)
    
    # 3. Load Fragments and Transform with same matrix
    frag_coords_list = []
    for fp in fragments_paths:
        fc = get_atom_coords(fp)
        if len(fc) > 0:
            # Center then Rotate
            fc_centered = fc - center
            fc_aligned = np.dot(fc_centered, rot_matrix)
            frag_coords_list.append(fc_aligned)
            
    # 4. Plot
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    
    # Plot Pocket (White/Ghostly Surface)
    # We use Z (depth) to control alpha or sorting if we wanted, 
    # but for simple "ghost" we use low alpha and large marker.
    # X axis is p_aligned[:, 0], Y is p_aligned[:, 1]
    
    # Scatter Pocket
    ax.scatter(
        p_aligned[:, 0], 
        p_aligned[:, 1], 
        c='orange', # User requested orange surface in Step 3.5, but for this "scanning" step 
                    # white/ghost was mentioned in previous context. Let's use light gray/orange mix.
                    # User request: "el pocket se vea como superficie en blanco" (the pocket should look like a white surface)
        edgecolors='none', 
        alpha=0.1, 
        s=100,
        label='Pocket'
    )
    
    # Plot Fragments (Green)
    for i, fc in enumerate(frag_coords_list):
        ax.scatter(
            fc[:, 0], 
            fc[:, 1], 
            c='#00FF00', # Green
            s=20,
            alpha=0.8,
            marker='.'
        )
        # Connect atoms in fragment with lines? Maybe too messy. 
        # Fragments are small peptides. Just dots is safer for speed.
        
    ax.set_aspect('equal')
    ax.set_axis_off()
    ax.set_title(title, color='black')
    
    # Dark background like py3dmol default? Or white?
    # User screenshot showed dark background. Matplotlib usually white.
    # Let's check prompt. "check that the pocket surface looks like a white surface" 
    # implies dark background might be better to see "white".
    
    # Set dark background
    fig.patch.set_facecolor('#1a1a1a')
    ax.set_facecolor('#1a1a1a')
    
    # If pocket is white, we change color above
    ax.collections[0].set_color('white') 
    
    # Save to buffer
    buf = BytesIO()
    plt.savefig(buf, format='png', facecolor='#1a1a1a') # Save with background
    plt.close(fig)
    buf.seek(0)
    return buf.read()
