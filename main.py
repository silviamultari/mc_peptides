import argparse
import pandas as pd
import random
from numpy.random import default_rng
from pathlib import Path
import os
import numpy as np
import traceback

from functions.objects import PDBConverter, mutation, replace_ter_with_conect
from functions.vinadock import vina_dock
from functions.ligand_in_water import run_simulation
from functions.clustering import get_most_probable_conformation
from functions.simulate_complex import simulate_complex
from functions.analysis import run_analysis, reimage_trajectory, scores_extraction, find_negatives
from functions.score_trajectory import score_complex
from functions.af2 import alphafold2

from functions.checkpoints import save_checkpoint
import pickle

# Flags configuration
parser = argparse.ArgumentParser(description="Run mutation-docking loop with specified paths and parameters.")
parser.add_argument("--base_dir", required=True, help="Base directory path")  # ALLIANCE
parser.add_argument("--input_rece_file", required=True, help="Receptor file name")  # rece_1AKJ
parser.add_argument("--ref_lig_file", required=True, help="Reference ligand file name")  # lig_1AKJ
parser.add_argument("--out_base", type=str, default="system", help="Name of the outputs for MD simulation")
parser.add_argument("--max_template_date", default="2021-11-01", help="AlphaFold2 max template date")
parser.add_argument("--model_preset", default="monomer", help="AlphaFold2 model preset")
parser.add_argument("--steps", type=int, default=5000000, help="Number of simulation steps")
parser.add_argument("--step_size", type=float, default=0.002, help="Size of each simulation step")
parser.add_argument("--consensus_threshold", type=int, default=3, help="Consensus threshold for acceptance")
parser.add_argument("--temperature", type=float, default=0.8, help="Metropolis temperature")
parser.add_argument("--ref_seq", type=str, default="CAAAAAAAAAAC", help="Reference sequence")  # E.g., CAADQTQDTEAAC
parser.add_argument("--keep_pos", type=int, nargs='+', required=True, help="Positions to keep unchanged")  # E.g., 0 5 6 7 12
parser.add_argument("--iter", type=int, default=500, help="Number of iterations of the loop")

args = parser.parse_args()


# Folders and paths configuration
BASE_DIR = Path(args.base_dir)
REF_INPUT_FOLDER = BASE_DIR / "inputs" / "ref" / "pdb"
REF_OUTPUT_FOLDER = BASE_DIR / "inputs" / "ref" / "pdbqt"
LIG_INPUT_FOLDER = BASE_DIR / "inputs" / "ligands" / "pdb"
LIG_OUTPUT_FOLDER = BASE_DIR / "inputs" / "ligands" / "pdbqt"
FASTA_DIR = BASE_DIR / "inputs" / "ligands" / "fasta"
MINIMIZED_OUTPUT_PATH = BASE_DIR / "docking" / "minimized_outputs"
DOCKING_OUT_PATH = BASE_DIR / "docking" / "outputs"
DOCKING_PDB_PATH = DOCKING_OUT_PATH / "pdb"
OUTPUT_PDB = BASE_DIR / "output" / "complexes"
CHECKPOINT_DIR = BASE_DIR / "checkpoints"

directories = [
    REF_INPUT_FOLDER, 
    REF_OUTPUT_FOLDER, 
    LIG_INPUT_FOLDER, 
    LIG_OUTPUT_FOLDER, 
    FASTA_DIR, 
    MINIMIZED_OUTPUT_PATH, 
    DOCKING_OUT_PATH, 
    DOCKING_PDB_PATH,
    OUTPUT_PDB,
    CHECKPOINT_DIR
]

for directory in directories:
    if not directory.exists():
        directory.mkdir(parents=True)
        print(f"New directory: {directory}")

# Ref definition
ref_converter = PDBConverter(REF_INPUT_FOLDER, REF_OUTPUT_FOLDER)
ref_pdbqt = ref_converter.pdb2pdbqt_rece(args.base_dir)

# Variables initialization
rng = default_rng()
xrange = (0.0, 1.0)
T = args.temperature
consensus_threshold = args.consensus_threshold
ref_seq = args.ref_seq
keep_pos = args.keep_pos  
amino_acids = ['A', 'R', 'N', 'D', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
input_rece_file = args.input_rece_file
ref_lig_file = args.ref_lig_file
receptor_path = f"{REF_OUTPUT_FOLDER}/{input_rece_file}.pdbqt"
ref_ligand_path = f"{REF_OUTPUT_FOLDER}/{ref_lig_file}.pdbqt"

old_pos = "ref"
old_i = "seq"

out_base = args.out_base

# Checkpoints creation
checkpoint_file_base = CHECKPOINT_DIR
latest_checkpoint_file = max(CHECKPOINT_DIR.glob("checkpoint_*.pkl"), default=None, key=os.path.getctime)

if latest_checkpoint_file:
    with open(latest_checkpoint_file, 'rb') as f:
        latest_checkpoint = pickle.load(f)
    ref_seq = latest_checkpoint.get('ref_seq', ref_seq)
    old_pos = latest_checkpoint.get('old_pos', old_pos)
    old_i = latest_checkpoint.get('old_i', old_i)
    start_iteration = latest_checkpoint.get('iteration', 0)
    print(f"Restarting from iteration {start_iteration} with ref_seq {ref_seq}")
else:
    start_iteration = 0
    print(f"Restarting from iteration {start_iteration} with ref_seq {ref_seq}")

# Main mutation loop
for i in range(start_iteration,500):
    valid_indices = [j for j in range(len(ref_seq)) if j not in keep_pos]
    pos = random.choice(valid_indices)
    mutate_seq = mutation(ref_seq, pos, amino_acids)
    fasta_string = f">sequence_{pos}_{i}\n{mutate_seq}"
    fasta_path = FASTA_DIR / f"seq_{pos}_{i}.fasta"
            
    with open(fasta_path, 'w') as file:
        file.write(fasta_string)
    
    if not os.path.exists(f'{LIG_INPUT_FOLDER}/seq_{pos}_{i}'):
        os.makedirs(f'{LIG_INPUT_FOLDER}/seq_{pos}_{i}')

    alphafold2(fasta_path, args.max_template_date, LIG_INPUT_FOLDER)
    
    try:
        run_simulation(f"{LIG_INPUT_FOLDER}/seq_{pos}_{i}/ranked_0.pdb", 
                       f"{LIG_INPUT_FOLDER}/seq_{pos}_{i}/traj_mut_{pos}_{i}.pdb", 
                       steps=20000, report_interval=300)
    except Exception as e:
        print(f"Error during pos={pos}, i={i} execution: {e}")
        traceback.print_exc()
        save_checkpoint(i, checkpoint_file_base, ref_seq, old_pos, old_i, interval=10)
        continue

    representative_frame, largest_cluster_label, cluster_counts = get_most_probable_conformation(
        f"{LIG_INPUT_FOLDER}/seq_{pos}_{i}/traj_mut_{pos}_{i}.pdb", pos, i, 
        output_path=f"{LIG_INPUT_FOLDER}/seq_{pos}_{i}", subsample_rate=100)
    
    replace_ter_with_conect(f"{LIG_INPUT_FOLDER}/seq_{pos}_{i}/rep_mut_{pos}_{i}.pdb")
    
    lig_converter = PDBConverter(f"{LIG_INPUT_FOLDER}/seq_{pos}_{i}", f"{LIG_OUTPUT_FOLDER}/seq_{pos}_{i}")
    new_ligand = lig_converter.pdb2pdbqt_lig(args.base_dir, pos, i)
    new_ligand_path = f"{LIG_OUTPUT_FOLDER}/seq_{pos}_{i}/rep_mut_{pos}_{i}.pdbqt"
    
    try:
        best_docking_score, mean_docking_score = vina_dock(
            receptor_path, ref_ligand_path, new_ligand_path, 
            f"{MINIMIZED_OUTPUT_PATH}/min_{pos}_{i}.pdbqt", 
            f"{DOCKING_OUT_PATH}/out_{pos}_{i}.pdbqt", 
            box_size=[20, 20, 20], exhaustiveness=32, n_poses=10
        )
    except Exception as e:
        print(f"Error during pos={pos}, i={i} execution: {e}")
        traceback.print_exc()
        save_checkpoint(i, checkpoint_file_base, ref_seq, old_pos, old_i, interval=10)
        continue

    if best_docking_score < -5:
        if not os.path.exists(f'{OUTPUT_PDB}/seq_{pos}_{i}'):
            os.makedirs(f'{OUTPUT_PDB}/seq_{pos}_{i}')
            
        pdb_converter = PDBConverter(DOCKING_OUT_PATH, DOCKING_OUT_PATH / "pdb")
        pdb_ligand = pdb_converter.pdbqt2pdb(args.base_dir, pos, i)
        
        try:
            simulate_complex(receptor_path, f"{DOCKING_OUT_PATH}/pdb/out_{pos}_{i}.pdb",  
                             output_base=f'{OUTPUT_PDB}/seq_{pos}_{i}/system_{pos}_{i}', 
                             steps=args.steps, solvate=True)
        except Exception as e:
            print(f"Error during pos={pos}, i={i} execution: {e}")
            traceback.print_exc()
            save_checkpoint(i, checkpoint_file_base, ref_seq, old_pos, old_i, interval=10)
            continue
        
        run_analysis(BASE_DIR, out_base, pos, i)
        reimage_trajectory(BASE_DIR, out_base, pos, i)
        
        replace_ter_with_conect(f"{OUTPUT_PDB}/seq_{pos}_{i}/system_{pos}_{i}_traj_reimaged.pdb")
        replace_ter_with_conect(f"{OUTPUT_PDB}/seq_{pos}_{i}/system_{pos}_{i}_reimaged.pdb")
        
        score_complex(
            score_list=["bach", "pisa", "zrank", "irad", "firedock", "bmf-bluues"],
            pdb_name=f'system_{pos}_{i}_traj_reimaged',
            src_route=f"{BASE_DIR}/functions",
            local_path=os.getcwd(),
            traj_path=f'{OUTPUT_PDB}/seq_{pos}_{i}',
            chain_join=["A", "B"],
            binder="C"
        )
        
        new_scores = scores_extraction(out_base, pos, i)
        old_scores = scores_extraction(out_base, old_pos, old_i)
        
        delta_scores, n_neg = find_negatives(new_scores, old_scores)
        
        Pacc = min(1, np.exp(-delta_scores[0] / T))
        x = rng.uniform(*xrange)
        if n_neg >= consensus_threshold or (n_neg < consensus_threshold and x <= Pacc):
            ref_seq = mutate_seq
            old_pos = pos
            old_i = i

        # Checkpoint creation every 10 iteration
        save_checkpoint(i, checkpoint_file_base, ref_seq, old_pos, old_i, interval=10)
