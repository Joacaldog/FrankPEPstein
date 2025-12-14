
import os
import glob
import json
import ipywidgets as widgets
from IPython.display import display

def visualize_candidates():
    # Dependency Check inside function to ensure environment
    try:
        import py3Dmol
    except ImportError:
        print("py3Dmol not installed.")
        return

    initial_path = os.getcwd()
    
    # 1. Load State
    state_file = "pipeline_state.json"
    receptor_path = None
    extracted_pocket_path = None
    box_center = None
    box_size = None
    
    if os.path.exists(state_file):
        try:
            with open(state_file, "r") as f:
                state = json.load(f)
                receptor_path = state.get("receptor_filename")
                box_center = state.get("box_center")
                box_size = state.get("box_size")
        except:
            pass
            
    # Pocket Path
    standard_pocket_path = os.path.join(initial_path, "pocket.pdb")
    if os.path.exists(standard_pocket_path):
        extracted_pocket_path = standard_pocket_path
    
    # 2. Find Candidates
    run_base = os.path.join(initial_path, "FrankPEPstein_run")
    candidate_folders = glob.glob(os.path.join(run_base, "frankPEPstein_*", "top_*_peps"))
    
    if not candidate_folders:
        print("❌ No candidate results found. Run fragments generation first.")
        return
        
    target_folder = sorted(candidate_folders, key=os.path.getmtime, reverse=True)[0]
    print(f"Visualizing candidates from: {target_folder}")
    
    pdb_files = glob.glob(os.path.join(target_folder, "*.pdb"))
    if not pdb_files:
        print("❌ No PDB files found in target folder.")
        return

    # 3. Render
    try:
        view = py3Dmol.view(width=1000, height=800, js='https://3dmol.org/build/3Dmol.js')
        
        # Receptor: White Surface + Cartoon
        if receptor_path and os.path.exists(receptor_path):
            with open(receptor_path, 'r') as f:
                view.addModel(f.read(), "pdb")
            # Cartoon for structure
            view.setStyle({'model': -1}, {'cartoon': {'color': 'white', 'opacity': 0.4}})
            # Surface for volume
            view.addSurface(py3Dmol.SES, {'opacity': 0.3, 'color': 'white'})

        # Pocket: Orange Surface
        if extracted_pocket_path and os.path.exists(extracted_pocket_path):
            with open(extracted_pocket_path, 'r') as f:
                view.addModel(f.read(), "pdb")
            # Orange surface as requested
            view.addSurface(py3Dmol.SES, {'opacity': 0.6, 'color': 'orange'})
            # Also show atoms slightly to define center
            view.setStyle({'model': -1}, {'sphere': {'radius': 0.5, 'color': 'orange', 'opacity': 0.0}}) # Hidden atoms, just surface

        # Gridbox: Thick Red Bars
        if box_center and box_size:
            cx, cy, cz = box_center
            sx, sy, sz = box_size
            
            min_x, max_x = cx - sx/2, cx + sx/2
            min_y, max_y = cy - sy/2, cy + sy/2
            min_z, max_z = cz - sz/2, cz + sz/2
            
            def draw_edge(p1, p2):
                view.addLine({
                    'start': p1, 'end': p2,
                    'color': 'red', 'linewidth': 10 # Thick
                })
            
            # Corners
            c000 = {'x':min_x, 'y':min_y, 'z':min_z}
            c100 = {'x':max_x, 'y':min_y, 'z':min_z}
            c010 = {'x':min_x, 'y':max_y, 'z':min_z}
            c110 = {'x':max_x, 'y':max_y, 'z':min_z}
            c001 = {'x':min_x, 'y':min_y, 'z':max_z}
            c101 = {'x':max_x, 'y':min_y, 'z':max_z}
            c011 = {'x':min_x, 'y':max_y, 'z':max_z}
            c111 = {'x':max_x, 'y':max_y, 'z':max_z}
            
            # Edges
            draw_edge(c000, c100); draw_edge(c100, c110); draw_edge(c110, c010); draw_edge(c010, c000) # Bottom
            draw_edge(c001, c101); draw_edge(c101, c111); draw_edge(c111, c011); draw_edge(c011, c001) # Top
            draw_edge(c000, c001); draw_edge(c100, c101); draw_edge(c110, c111); draw_edge(c010, c011) # Sides

        # Peptides: Sticks
        for pdb_file in pdb_files:
            with open(pdb_file, 'r') as f:
                view.addModel(f.read(), "pdb")
            # Green/Multicolor sticks
            view.setStyle({'model': -1}, {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.2}})

        view.zoomTo()
        return view.show()
        
    except Exception as e:
        print(f"Visualization Error: {e}")
