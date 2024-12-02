import os
from glob import glob
import subprocess

        
class PDBConverter:
    def __init__(self, input_folder, output_folder):
        self.input_folder = input_folder
        self.output_folder = output_folder

    def pdb2pdbqt_rece(self, base):
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

        pdb_files = glob(os.path.join(self.input_folder, '*.pdb'))

        for pdb_file in pdb_files:
            pdbqt_file = os.path.join(self.output_folder, os.path.basename(pdb_file).replace('.pdb', '.pdbqt'))
            command = [f'{base}functions/ADFRsuite-1.1dev/ADFRsuite_x86_64Linux_1.1dev/docking/bin/prepare_receptor',  '-r' , pdb_file, '-o', pdbqt_file]
            subprocess.call(command)

    # def pdb2pdbqt_lig(self, pos, i):
    #     if not os.path.exists(self.output_folder):
    #         os.makedirs(self.output_folder)

    #     pdb_files = glob(os.path.join(self.input_folder, '*.pdb'))

    #     for pdb_file in pdb_files:
    #         pdbqt_file = os.path.join(self.output_folder, os.path.basename(pdb_file).replace(f'mut_{pos}_{i}.pdb', f'mut_{pos}_{i}.pdbqt'))
    #         command = ['ALLIANCE/functions/ADFRsuite-1.1dev/ADFRsuite_x86_64Linux_1.1dev/docking/bin/obabel', '-ipdb', pdb_file, '-opdbqt', pdbqt_file, '--gen3d', '-h', '-m']
    #         # command = ['ALLIANCE/functions/ADFRsuite-1.1dev/ADFRsuite_x86_64Linux_1.1dev/docking/bin/prepare_ligand',  '-l' , pdb_file, '-o', pdbqt_file]
    #         subprocess.call(command)
    
    def pdb2pdbqt_lig(self, base, pos, i):
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

        # Prendi un singolo file pdb basato su 'pos' e 'i'
        pdb_file = os.path.join(self.input_folder, f'rep_mut_{pos}_{i}.pdb')
        
        # Verifica che il file esista prima di continuare
        if os.path.exists(pdb_file):
            pdbqt_file = os.path.join(self.output_folder, f'rep_mut_{pos}_{i}.pdbqt')
            # command = ['ALLIANCE/functions/ADFRsuite-1.1dev/ADFRsuite_x86_64Linux_1.1dev/docking/bin/obabel', pdb_file, '-O', pdbqt_file] # '--gen3d', '-h', '-m']
            command = [f'{base}functions/ADFRsuite-1.1dev/ADFRsuite_x86_64Linux_1.1dev/docking/bin/prepare_ligand',  '-l' , pdb_file, '-o', pdbqt_file]
            subprocess.call(command)
        else:
            print(f"Il file {pdb_file} non esiste.")

            
    def pdbqt2pdb(self, base, pos, i):
        """
        Converts a PDBQT file to PDB format.

        Args:
            input_file (str): Path to the input PDBQT file.
        """
        pdbqt_files = glob(os.path.join(self.input_folder, f'out_{pos}_{i}.pdbqt'))
        for pdbqt_file in pdbqt_files:
            
            with open(pdbqt_file, "r") as f:
                lines = f.readlines()
        
            # Write back the modified contents
            with open(pdbqt_file, "w") as f:
                in_first_model = False
                in_section_to_remove = False
            
                for line in lines:
                    if line.startswith("MODEL"):
                        model_number = int(line.strip()[6:])
                        if model_number == 1:
                            in_first_model = True
                        else:
                            # Stop writing from the second "MODEL" onwards
                            in_first_model = False
                            continue
                    
                    if line.startswith("ENDROOT"):
                        in_section_to_remove = True
                    
                    if not in_section_to_remove and in_first_model:
                        f.write(line)
                    
                    if line.startswith("ENDMDL"):
                        in_section_to_remove = False
   
            pdb_file = os.path.join(self.output_folder, os.path.basename(pdbqt_file).replace('.pdbqt', '.pdb'))
            command = [f'{base}functions/ADFRsuite-1.1dev/ADFRsuite_x86_64Linux_1.1dev/docking/bin/obabel', '-ipdbqt', pdbqt_file, '-O', pdb_file]
            subprocess.call(command)
    


import random

def mutation(sequence, position, amino_acids):
    """
    Sostituisce un amminoacido in una sequenza con uno casuale dalla lista di amminoacidi.

    Args:
        sequence (str): La sequenza originale di amminoacidi.
        position (int): La posizione dell'amminoacido da sostituire (0-indice).
        amino_acids (list): Lista di amminoacidi possibili.

    Returns:
        str: La sequenza mutata con l'amminoacido sostituito.
    """
    if position < 0 or position >= len(sequence):
        raise ValueError("La posizione deve essere valida per la sequenza.")

    # Scegli un amminoacido casuale dalla lista
    new_amino_acid = random.choice(amino_acids)

    # Crea una nuova sequenza con l'amminoacido mutato
    mutated_sequence = sequence[:position] + new_amino_acid + sequence[position + 1:]

    return mutated_sequence



def replace_ter_with_conect(file_path):
    with open(file_path, 'r') as file:
        pdb_data = file.readlines()

    # Find the indexes for the SG atoms
    sg_indexes = []
    for line in pdb_data:
        if line.startswith("ATOM") and " SG " in line:
            atom_index = int(line[6:11].strip())
            sg_indexes.append(atom_index)

    if len(sg_indexes) >= 2:
        # Replace the TER line with the new CONECT lines
        for i, line in enumerate(pdb_data):
            if line.startswith("END"):
                pdb_data[i] = f"CONECT{sg_indexes[0]:>5}{sg_indexes[1]:>5}\n"
                pdb_data.insert(i + 1, f"CONECT{sg_indexes[1]:>5}{sg_indexes[0]:>5}\n")
                pdb_data.insert(i + 2, "END")
                break

    # Overwrite the file with the modified data
    with open(file_path, 'w') as file:
        file.writelines(pdb_data)
