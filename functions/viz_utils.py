
import os
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import ConvexHull, Delaunay
from io import BytesIO

# --- Constants ---
ATOM_RADII = {'C': 1.7, 'N': 1.55, 'O': 1.52, 'S': 1.8, 'H': 1.2, 'P': 1.8}
ATOM_COLORS = {'C': '#00FF00', 'N': 'blue', 'O': 'red', 'S': 'yellow', 'P': 'orange', 'H': 'white'} 

def get_atom_data(pdb_file, atom_type=None):
    """
    Parses PDB and returns list of dicts: {'x', 'y', 'z', 'element', 'name', 'resid', 'chain'}
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
                    
                    resid = int(line[22:26].strip())
                    chain = line[21:22].strip()
                    
                    atoms.append({
                        'x': x, 'y': y, 'z': z,
                        'element': element,
                        'name': name,
                        'resid': resid,
                        'chain': chain
                    })
                except ValueError:
                    pass
    return atoms

def get_points_in_hull(hull, density=1.0):
    """
    Generates a grid of points inside a ConvexHull.
    """
    # Bounding box
    min_x, min_y, min_z = hull.points[hull.vertices].min(axis=0)
    max_x, max_y, max_z = hull.points[hull.vertices].max(axis=0)
    
    # Grid
    grid_x, grid_y, grid_z = np.mgrid[min_x:max_x:density, min_y:max_y:density, min_z:max_z:density]
    grid_points = np.vstack([grid_x.ravel(), grid_y.ravel(), grid_z.ravel()]).T
    
    # Check identifying using Delaunay
    # Using scipy.spatial.Delaunay to check if points are in hull (simplex check)
    # This can be slow for large grids.
    # Alternatives: Just use the hull equations?
    
    delaunay = Delaunay(hull.points[hull.vertices])
    valid = delaunay.find_simplex(grid_points) >= 0
    
    return grid_points[valid]

def render_static_view(receptor_path, pocket_path, box_center, box_size, fragments_paths, title="Processing..."):
    # print("Loading visualization...")
    
    # 1. Load Data
    # Pocket: Use CA for hull calc? Or all atoms for volume?
    # User said: "between residues of fpocket file" -> better use all atoms for volume definition
    pocket_atoms = get_atom_data(pocket_path)
    if not pocket_atoms: return None
    
    pocket_coords = np.array([[a['x'], a['y'], a['z']] for a in pocket_atoms])
    
    
    # 2. Setup Plot
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_facecolor('black')
    fig.patch.set_facecolor('black')

    # 5. Pocket Volume (Red Spheres)
    # Hull of pocket atoms
    if len(pocket_coords) > 4:
        try:
            hull_p = ConvexHull(pocket_coords)
            # Generate internal points
            # Density approx 1.5A for speed vs look
            internal_pts = get_points_in_hull(hull_p, density=2.0)
            
            if len(internal_pts) > 0:
                ax.scatter(internal_pts[:,0], internal_pts[:,1], internal_pts[:,2], 
                           color='red', s=20, alpha=0.3, edgecolors='none', label='Volume')
        except: pass

    # 6. Peptide Candidates
    for fp in fragments_paths:
        f_atoms = get_atom_data(fp)
        f_coords = np.array([[a['x'], a['y'], a['z']] for a in f_atoms])
        if len(f_coords) > 0:
            # Sticks
            # Simple connectivity based on distance (neighbor graph)
            # Matplotlib 3D lines
            # Just plotting atoms for now? Or basic bond logic
            
            # Draw atoms
            colors = [ATOM_COLORS.get(a['element'], 'gray') for a in f_atoms]
            ax.scatter(f_coords[:,0], f_coords[:,1], f_coords[:,2], c=colors, s=30, alpha=1.0)
            
            # Draw bonds (simple loop)
            for i in range(len(f_atoms)):
                for j in range(i+1, len(f_atoms)):
                    p1 = f_coords[i]
                    p2 = f_coords[j]
                    dist = np.linalg.norm(p1-p2)
                    if dist < 1.8: # Bond length
                        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], color='white', linewidth=1)

    # 7. Gridbox (Red Corners/Edges)
    cx, cy, cz = box_center
    sx, sy, sz = box_size
    
    # Draw corners/edges? Use wireframe box
    # Plot transparent box?
    # Just red edges
    
    def plot_cube(center, size, ax):
        ox, oy, oz = center
        l, w, h = size
        
        x = [ox-l/2, ox+l/2]
        y = [oy-w/2, oy+w/2]
        z = [oz-h/2, oz+h/2]
        
        # Draw edges
        # Bottom
        ax.plot([x[0], x[1]], [y[0], y[0]], [z[0], z[0]], color='white', linewidth=2)
        ax.plot([x[0], x[1]], [y[1], y[1]], [z[0], z[0]], color='white', linewidth=2)
        ax.plot([x[0], x[0]], [y[0], y[1]], [z[0], z[0]], color='white', linewidth=2)
        ax.plot([x[1], x[1]], [y[0], y[1]], [z[0], z[0]], color='white', linewidth=2)
        # Top
        ax.plot([x[0], x[1]], [y[0], y[0]], [z[1], z[1]], color='white', linewidth=2)
        ax.plot([x[0], x[1]], [y[1], y[1]], [z[1], z[1]], color='white', linewidth=2)
        ax.plot([x[0], x[0]], [y[0], y[1]], [z[1], z[1]], color='white', linewidth=2)
        ax.plot([x[1], x[1]], [y[0], y[1]], [z[1], z[1]], color='white', linewidth=2)
        # Verticals
        ax.plot([x[0], x[0]], [y[0], y[0]], [z[0], z[1]], color='white', linewidth=2)
        ax.plot([x[1], x[1]], [y[0], y[0]], [z[0], z[1]], color='white', linewidth=2)
        ax.plot([x[0], x[0]], [y[1], y[1]], [z[0], z[1]], color='white', linewidth=2)
        ax.plot([x[1], x[1]], [y[1], y[1]], [z[0], z[1]], color='white', linewidth=2)

    plot_cube(box_center, box_size, ax)

    # Camera
    ax.set_axis_off()
    ax.set_title(title, color='white')
    
    # Auto-center view on pocket
    max_range = max(box_size) * 1.5
    ax.set_xlim(cx - max_range/2, cx + max_range/2)
    ax.set_ylim(cy - max_range/2, cy + max_range/2)
    ax.set_zlim(cz - max_range/2, cz + max_range/2)
    
    buf = BytesIO()
    plt.savefig(buf, format='png', facecolor='black')
    plt.close(fig)
    buf.seek(0)
    return buf.read()
