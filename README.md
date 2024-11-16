# NAME: Monte Carlo Algorithm for Cyclic Peptide Optimization

Here we introduce **NAME**, an open-source Monte Carlo algorithm tailored to optimize cyclic peptide binders.  
The workflow begins with introducing random mutations to the peptide sequence, evaluated by sampling their bound conformations through molecular dynamics simulations.  
Interaction strengths between the peptide and target protein are assessed using diverse scoring metrics, and decisions to accept or discard mutations are guided by a consensus strategy that integrates these scores.  
By iteratively refining sequences, the method efficiently uncovers new candidates with enhanced binding affinity to their targets.

## Dependencies

To run the code, install the packages listed in `requirements.txt`.  
It is also necessary to download **ADFRsuite-1.1dev**: the download file is in the `functions` folder, where the suite must be installed.

## Usage

### Reference system definition
First, you have to define a reference system, consisting of the target receptor, its known ligand, and the cyclic peptide you want to optimize.
```bash
usage: step_0.py [-h] --ligand_pdb LIGAND_PDB --cyclic_pdb CYCLIC_PDB --receptor_pdb RECEPTOR_PDB [--steps STEPS] --src_route SRC_ROUTE

Run simulation steps and analysis.

optional arguments:
  -h, --help            show this help message and exit
  --ligand_pdb LIGAND_PDB
                        Path to the reference ligand PDB file. WARNING: it has to be inputs/ref/pdb/filename.pdb
  --cyclic_pdb CYCLIC_PDB
                        Path to the starting cyclic peptide PDB file. WARNING: it has to be inputs/ref/pdb/filename.pdb
  --receptor_pdb RECEPTOR_PDB
                        Path to the receptor PDB file. WARNING: it has to be inputs/ref/pdb/filename.pdb
  --steps STEPS         Number of simulation steps.
  --src_route SRC_ROUTE
                        Path to functions directory
```

### Main algorithm run
Then, you can run the Monte Carlo algorithm.

```bash
usage: main.py [-h] --base_dir BASE_DIR --input_rece_file INPUT_RECE_FILE
                 --ref_lig_file REF_LIG_FILE [--max_template_date MAX_TEMPLATE_DATE]
                 [--model_preset MODEL_PRESET] [--steps STEPS]
                 [--consensus_threshold CONSENSUS_THRESHOLD]
                 [--temperature TEMPERATURE] --ref_seq REF_SEQ
                 --keep_pos KEEP_POS [KEEP_POS ...]

Run the main mutation loop with specified paths and parameters.

optional arguments:
  -h, --help            show this help message and exit
  --base_dir BASE_DIR   Base directory path
  --input_rece_file INPUT_RECE_FILE
                        Receptor file name
  --ref_lig_file REF_LIG_FILE
                        Reference ligand file name
  --max_template_date MAX_TEMPLATE_DATE
                        AlphaFold2 max template date (default: 2021-11-01)
  --model_preset MODEL_PRESET
                        AlphaFold2 model preset (default: monomer)
  --steps STEPS         Number of simulation steps (default: 5000000)
  --consensus_threshold CONSENSUS_THRESHOLD
                        Consensus threshold for acceptance (default: 3)
  --temperature TEMPERATURE
                        Metropolis temperature (default: 0.8)
  --ref_seq REF_SEQ     Reference sequence (default: CAAAAAAAAAAAC)
  --keep_pos KEEP_POS [KEEP_POS ...]
                        Positions to keep unchanged
  --iter ITER           Number of iterations of the loop (default: 500)
```




