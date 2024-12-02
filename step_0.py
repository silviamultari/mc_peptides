import numpy as np
import subprocess
import os
from functions.simulate_complex import simulate_complex
from functions.score_trajectory import score_complex

def parse_pdb(file_path):
    coordinates = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                coo_line = line.split()
                x = float(coo_line[6])
                y = float(coo_line[7])
                z = float(coo_line[8])
                coordinates.append((x, y, z))
    coordinates = np.array(coordinates)
    min_coord = np.min(coordinates, axis=0)
    max_coord = np.max(coordinates, axis=0)
    center = (min_coord + max_coord) / 2
    return center

def change_pdb(ref_path, new_path):
    # Definizione del centro geometrico di box e ligando
    box_center = parse_pdb(ref_path)
    print("Center of the box set at: ", box_center)
    lig_center = parse_pdb(new_path)
    print("Center of the ligand set at: ", lig_center)
    # Calcolo del vettore di traslazione per ogni punto
    v_trans = lig_center - box_center
    print("Translation vector: ", v_trans)
    dx, dy, dz = v_trans
    
    with open(new_path, 'r') as file:
        lines = file.readlines()
    with open(new_path, 'w') as file:
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                coo_line = line.split()
                x = float(coo_line[6]) - dx
                y = float(coo_line[7]) - dy
                z = float(coo_line[8]) - dz

                # Reformat the line with the new coordinates
                new_line = f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
                file.write(new_line)
            else:
                file.write(line)
                
change_pdb("ALLIANCE/inputs/ref/pdb/lig_1AKJ.pdb", "ALLIANCE/inputs/ref/pdb/cyc_1AKJ.pdb")

simulate_complex("ALLIANCE/inputs/ref/pdb/rece_1AKJ.pdb", "ALLIANCE/inputs/ref/pdb/cyc_1AKJ.pdb",  
                            output_base="ALLIANCE/inputs/ref/pdb/ref_system", 
                            steps=250000, solvate=True)

def run_analysis(base):
    command = ['python', 'ALLIANCE/functions/analyse.py', '-p', f'{base}_minimised.pdb', '-t', f'{base}_traj.dcd', '-o', f'{base}_reimaged', '-r']
    subprocess.call(command)
                
def reimage_trajectory(base):
    command = ['mdconvert', f'{base}_reimaged.dcd', '-o', f'{base}_traj_reimaged.pdb', '-t', f'{base}_reimaged.pdb']
    subprocess.call(command)
    
    
run_analysis(base="ALLIANCE/inputs/ref/pdb/ref_system")
reimage_trajectory(base="ALLIANCE/inputs/ref/pdb/ref_system")

score_complex(score_list=["bach","pisa","zrank","irad","firedock","bmf-bluues"],
                        pdb_name=f'ref_system_traj_reimaged',src_route="ALLIANCE/functions",
                        local_path=os.getcwd(),traj_path=f'ALLIANCE/inputs/ref/pdb',
                        chain_join=["A", "B"],binder="C")

