
import os
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
from io import BytesIO

# --- Constants ---
ATOM_RADII = {'C': 1.7, 'N': 1.55, 'O': 1.52, 'S': 1.8, 'H': 1.2, 'P': 1.8}
ATOM_COLORS = {'C': '#00FF00', 'N': 'blue', 'O': 'red', 'S': 'yellow', 'P': 'orange', 'H': 'white'} # Carbon Green as requested

def get_atom_data(pdb_file, atom_type=None):
    """
    Parses PDB and returns list of dicts: {'x', 'y', 'z', 'element', 'name', 'resid'}
    """
    atoms = []
    if not os.path.exists(pdb_file):
        return atoms
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                name = line[12:16].strip()
                if atom_type and name != atom_type:
                    continue
                
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    element = line[76:78].strip()
                    if not element:
                        element = name[0] # Fallback
                    
                    resid = line[22:26].strip()
                    
                    atoms.append({
                        'x': x, 'y': y, 'z': z,
                        'element': element,
                        'name': name,
                        'resid': resid
                    })
                except ValueError:
                    pass
    return atoms

def transform_atoms(atoms, center, rot_matrix):
    """Applies centering and rotation to atom list."""
    coords = np.array([[a['x'], a['y'], a['z']] for a in atoms])
    if len(coords) == 0: return []
    
    centered = coords - center
    aligned = np.dot(centered, rot_matrix)
    
    new_atoms = []
    for i, atom in enumerate(atoms):
        new_a = atom.copy()
        new_a['x'], new_a['y'], new_a['z'] = aligned[i]
        new_atoms.append(new_a)
    return new_atoms

def align_perspective(pocket_atoms):
    """Calculates rotation matrix based on pocket PCA + Isometric tilt."""
    coords = np.array([[a['x'], a['y'], a['z']] for a in pocket_atoms])
    if len(coords) < 3:
        return np.eye(3), np.mean(coords, axis=0) if len(coords)>0 else [0,0,0]
        
    center = np.mean(coords, axis=0)
    centered = coords - center
    cov = np.cov(centered, rowvar=False)
    evals, evecs = np.linalg.eigh(cov)
    # PCA aligns with principal axes. Often flat.
    idx = evals.argsort()[::-1]
    evecs = evecs[:, idx]
    
    # --- ADD ISOMETRIC-LIKE OFFSET ---
    # We apply a rigid rotation AFTER aligning to PCA to show depth.
    # Rotation X (45 deg) * Rotation Y (30 deg)
    
    # Rx (45)
    theta_x = np.radians(45)
    rx = np.array([
        [1, 0, 0],
        [0, np.cos(theta_x), -np.sin(theta_x)],
        [0, np.sin(theta_x), np.cos(theta_x)]
    ])
    
    # Ry (30)
    theta_y = np.radians(30)
    ry = np.array([
        [np.cos(theta_y), 0, np.sin(theta_y)],
        [0, 1, 0],
        [-np.sin(theta_y), 0, np.cos(theta_y)]
    ])
    
    iso_rot = np.dot(rx, ry)
    
    # Final Matrix = PCA * Isometric
    # Note: evecs columns are axes. Mapping world to aligned space.
    # Transpose of evecs rotates vector TO principal frame.
    # Then we apply iso_rot.
    # Combined = Iso * (PCA_Transpose) 
    # But here we return a matrix 'rot_matrix' such that: v_new = dot(v_old - c, rot_matrix)
    # So we want M such that v_new = (v_old - c) @ M
    
    # evecs maps local -> world. evecs.T maps world -> local.
    # If we want v_local = (v-c) @ evecs (if row vectors?)
    # Numpy eigh returns column eigenvectors. 
    # So v_world = evecs @ v_local. v_local = evecs.T @ v_world.
    # If using row vectors: v_local^T = v_world^T @ evecs
    
    # We want to rotate the cloud. 
    # v_aligned = (v-c) @ evecs   (This aligns PC1 to X, PC2 to Y, PC3 to Z approximately)
    # Then v_iso = v_aligned @ iso_rot.T
    
    final_matrix = np.dot(evecs, iso_rot.T)
    
    return final_matrix, center

def render_static_view(receptor_path, pocket_path, box_center, box_size, fragments_paths, title="Processing..."):
    # 1. Load Data
    pocket_atoms = get_atom_data(pocket_path)
    if not pocket_atoms: return None
    
    # receptor_atoms = get_atom_data(receptor_path) # Disabled for speed
    receptor_atoms = []
    
    # 2. Calculate Alignment (Focus on Pocket)
    rot_matrix, center = align_perspective(pocket_atoms)
    
    # 3. Transform All
    p_aligned = transform_atoms(pocket_atoms, center, rot_matrix)
    r_aligned = transform_atoms(receptor_atoms, center, rot_matrix)
    
    frag_aligned_list = []
    for fp in fragments_paths:
        f_atoms = get_atom_data(fp)
        frag_aligned_list.append(transform_atoms(f_atoms, center, rot_matrix))
        
    # Transform Gridbox Corners
    cx, cy, cz = box_center
    sx, sy, sz = box_size
    corners = [
        [cx-sx/2, cy-sy/2, cz-sz/2], [cx+sx/2, cy-sy/2, cz-sz/2],
        [cx-sx/2, cy+sy/2, cz-sz/2], [cx+sx/2, cy+sy/2, cz-sz/2],
        [cx-sx/2, cy-sy/2, cz+sz/2], [cx+sx/2, cy-sy/2, cz+sz/2],
        [cx-sx/2, cy+sy/2, cz+sz/2], [cx+sx/2, cy+sy/2, cz+sz/2]
    ]
    corners_aligned = np.dot(np.array(corners) - center, rot_matrix)
    
    # 4. Plot
    fig = plt.figure(figsize=(10, 8)) # Larger
    ax = fig.add_subplot(111)
    
    # Z-sorting for "3D-like" 2D plot
    # We collect all render operations and sort by Z depth
    render_queue = []
    
    # A. Gridbox Lines (Draw first/behind or make explicit?)
    # Lines don't occlude well in simple 2D scatter, check Z.
    # We will draw gridbox on top usually or behind? Pymol draws on top often.
    # Let's add segments to queue.
    edges = [
        (0,1), (1,3), (3,2), (2,0), # Back face? depends on rot.
        (4,5), (5,7), (7,6), (6,4), # Front face
        (0,4), (1,5), (2,6), (3,7)  # Connecting
    ]
    
    for s, e in edges:
        p1, p2 = corners_aligned[s], corners_aligned[e]
        z_avg = (p1[2] + p2[2])/2
        render_queue.append({
            'type': 'line', 'z': z_avg,
            'x': [p1[0], p2[0]], 'y': [p1[1], p2[1]],
            'color': 'red', 'width': 3
        })
        
    # B. Pocket (Surfaces)
    # Receptor commented out for speed
    # for atom in r_aligned: ...
        
    for atom in p_aligned:
        render_queue.append({
            'type': 'atom', 'z': atom['z'],
            'x': atom['x'], 'y': atom['y'],
            'r': 130, # Slightly larger/different
            'c': 'white', 'alpha': 0.8
        })

    # C. Fragments (Sticks)
    # Need connectivity. Simple distance check < 1.6A
    for frag in frag_aligned_list:
        # Atoms
        for i, atom in enumerate(frag):
            elem = atom['element']
            color = ATOM_COLORS.get(elem, 'gray')
            if elem == 'C': color = '#00FF00' # Explicit Green Carbon
            
            render_queue.append({
                'type': 'atom', 'z': atom['z'],
                'x': atom['x'], 'y': atom['y'],
                'r': 40, # Smaller than surface
                'c': color, 'alpha': 1.0
            })
            
            # Bonds (Forward check)
            for j in range(i+1, len(frag)):
                atom2 = frag[j]
                # Distance
                d2 = (atom['x']-atom2['x'])**2 + (atom['y']-atom2['y'])**2 + (atom['z']-atom2['z'])**2
                if d2 < 2.6: # 1.6**2 ~ 2.56
                    z_avg = (atom['z'] + atom2['z'])/2
                    render_queue.append({
                        'type': 'line', 'z': z_avg,
                        'x': [atom['x'], atom2['x']], 'y': [atom['y'], atom2['y']],
                        'color': 'white', 'width': 2 # Bond color usually gray/white
                    })

    # 5. Execute Sort and Draw
    render_queue.sort(key=lambda item: item['z'])
    
    for item in render_queue:
        if item['type'] == 'atom':
            ax.scatter(item['x'], item['y'], s=item['r'], c=item['c'], alpha=item['alpha'], edgecolors='none')
        elif item['type'] == 'line':
            ax.plot(item['x'], item['y'], c=item['color'], linewidth=item['width'], alpha=0.8)

    # 6. Zoom / Limits
    # Focus on Gridbox Size + Padding
    max_dim = max(sx, sy, sz)
    limit = (max_dim / 2) * 1.5 # 1.5x zoom padding
    
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    
    ax.set_aspect('equal')
    ax.set_axis_off()
    ax.set_title(title, color='white')
    
    # Dark Background
    fig.patch.set_facecolor('black')
    ax.set_facecolor('black')
    
    buf = BytesIO()
    plt.savefig(buf, format='png', facecolor='black')
    plt.close(fig)
    buf.seek(0)
    return buf.read()
