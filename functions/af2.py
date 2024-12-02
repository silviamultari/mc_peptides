import subprocess
def alphafold2(fasta_path, max_template_date, LIG_INPUT_FOLDER, num_predict="1",
               model="monomer_casp14", run_relax="False", data_dir="/mnt/data/alphafold_data"):
    command = [
    "bash", "/opt/miniconda3/envs/alphafold-2.3.1/run_alphafold.sh", "-d", f"{data_dir}", "-o", 
    f"{LIG_INPUT_FOLDER}", "-f", f"{fasta_path}", "-t", f"{max_template_date}", 
    "-m", model, "-r", run_relax, "-l", num_predict]
    
    result = subprocess.run(command, capture_output=True, text=True)
    print("Standard Output:\n", result.stdout)
    print("Standard Error:\n", result.stderr)
