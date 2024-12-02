import sys, time
from openmmforcefields.generators import SystemGenerator
import openmm
from openmm import app, unit, LangevinIntegrator
from openmm.app import PDBFile, Simulation, Modeller, DCDReporter, StateDataReporter
from pdbfixer import PDBFixer
from openmm.app import PDBxFile
import functions.utils as utils

def fix_pdb(pdb_path, output_path):
    """
    Fixes a PDB file using PDBFixer, adding missing atoms and hydrogens.

    Args:
        pdb_path (str): Path to the input PDB file.
        output_path (str): Path to save the fixed PDB file.
        is_peptide (bool): If True, treat the molecule as a peptide; otherwise, treat as a protein.
    """
    fixer = PDBFixer(filename=pdb_path)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.addMissingHydrogens(7.4)
    
    # Save the fixed PDB to the specified output path
    with open(output_path, 'w') as out_file:
        PDBFile.writeFile(fixer.topology, fixer.positions, out_file)
    print(f'Fixed PDB saved to {output_path}')
    

def replace_last_chain_identifier(pdb_file, new_chain_id):
    with open(pdb_file, 'r') as file:
        lines = file.readlines()

    # Find all unique chain identifiers
    chain_ids = set()
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            chain_ids.add(line[21])

    # Determine the last chain ID
    last_chain_id = None
    for line in reversed(lines):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            last_chain_id = line[21]
            break
    atoms = range(3585,3920)
    # Replace the last chain ID if it is 'A'
    if last_chain_id == 'A':
        lines = [line if not (line.startswith("ATOM") or line.startswith("HETATM")) or line[21] != 'A' or int(line[7:11]) not in atoms else line[:21] + new_chain_id + line[22:] for line in lines]

    # Write the modified lines to a new file
    with open(f"{pdb_file}", 'w') as file:
        file.writelines(lines)
        
        

def simulate_complex(protein_pdb_path, ligand_pdb_path, output_base='output', steps=5000, step_size=0.002,
                     friction_coeff=1, interval=1000, temperature=300, solvate=False, padding=10,
                     water_model="tip3p", positive_ion="Na+", negative_ion="Cl-", ionic_strength=0,
                     no_neutralize=False, equilibration_steps=200, protein_force_field='amber/ff14SB.xml',
                     water_force_field='amber/tip3p_standard.xml'):
    """
    Run an MD simulation for a protein-peptide complex, optionally adding a solvent box.

    Args:
        protein_pdb_path (str): Path to the protein PDB file.
        ligand_pdb_path (str): Path to the peptide ligand PDB file.
        output_base (str): Base name for output files.
        steps (int): Number of simulation steps.
        step_size (float): Step size in picoseconds.
        friction_coeff (float): Friction coefficient in 1/ps.
        interval (int): Reporting interval.
        temperature (float): Simulation temperature in Kelvin.
        solvate (bool): Whether to add a solvent box.
        padding (float): Padding for the solvent box in Ångströms.
        water_model (str): Water model for solvation.
        positive_ion (str): Positive ion for solvation.
        negative_ion (str): Negative ion for solvation.
        ionic_strength (float): Ionic strength for solvation in molar.
        no_neutralize (bool): Whether to skip neutralizing the system with ions.
        equilibration_steps (int): Number of equilibration steps.
        protein_force_field (str): Force field file for the protein and peptide.
        water_force_field (str): Force field file for water.
    """
    t0 = time.time()

    output_complex = output_base + '_complex.pdb'
    output_traj_dcd = output_base + '_traj.dcd'
    output_min = output_base + '_minimised.pdb'
    fixed_protein_pdb = output_base + '_fixed_protein.pdb'
    fixed_ligand_pdb = output_base + '_fixed_ligand.pdb'
    temperature_unit = temperature * unit.kelvin

    print('Processing', protein_pdb_path, 'and', ligand_pdb_path, 'with', steps, 'steps generating outputs',
          output_complex, output_min, output_traj_dcd)

    # Fix the protein PDB file
    print('Fixing protein PDB...')
    fix_pdb(protein_pdb_path, fixed_protein_pdb)

    # Fix the ligand PDB file
    print('Fixing peptide ligand PDB...')
    fix_pdb(ligand_pdb_path, fixed_ligand_pdb)

    # get the chosen or fastest platform
    platform = utils.get_platform()

    print('Reading ligand')
    ligand_pdb = PDBFile(fixed_ligand_pdb)
    ligand_topology = ligand_pdb.topology

    print('Preparing system')
    # Initialize a SystemGenerator with the protein force field for both the protein and peptide
    forcefield_kwargs = {'constraints': app.HBonds, 'rigidWater': True, 'removeCMMotion': False, 'hydrogenMass': 4 * unit.amu}
    system_generator = SystemGenerator(
        forcefields=[protein_force_field, water_force_field],
        forcefield_kwargs=forcefield_kwargs)

    # Use Modeller to combine the protein and ligand into a complex
    print('Reading protein')
    protein_pdb = PDBFile(fixed_protein_pdb)

    print('Preparing complex')
    modeller = Modeller(protein_pdb.topology, protein_pdb.positions)
    print('System has %d atoms' % modeller.topology.getNumAtoms())

    print('Adding peptide ligand...')
    modeller.add(ligand_topology, ligand_pdb.getPositions())
    print('System has %d atoms' % modeller.topology.getNumAtoms())

    # Solvate
    if solvate:
        print('Adding solvent...')
        modeller.addSolvent(system_generator.forcefield, model=water_model, padding=padding * unit.angstroms,
                            positiveIon=positive_ion, negativeIon=negative_ion,
                            ionicStrength=ionic_strength * unit.molar, neutralize=not no_neutralize)
        print('System has %d atoms' % modeller.topology.getNumAtoms())

    with open(output_complex, 'w') as outfile:
        PDBFile.writeFile(modeller.topology, modeller.positions, outfile)

    # Create the system using the SystemGenerator
    system = system_generator.create_system(modeller.topology)

    friction_coeff_unit = friction_coeff / unit.picosecond
    step_size_unit = step_size * unit.picoseconds
    duration = (step_size_unit * steps).value_in_unit(unit.nanoseconds)
    print('Simulating for {} ns'.format(duration))

    integrator = LangevinIntegrator(temperature_unit, friction_coeff_unit, step_size_unit)
    if solvate:
        system.addForce(openmm.MonteCarloBarostat(1 * unit.atmospheres, temperature_unit, 25))

    if system.usesPeriodicBoundaryConditions():
        print('Default Periodic box: {}'.format(system.getDefaultPeriodicBoxVectors()))
    else:
        print('No Periodic Box')

    simulation = Simulation(modeller.topology, system, integrator, platform=platform)
    context = simulation.context
    context.setPositions(modeller.positions)

    state = simulation.context.getState(getEnergy=True)
    print(f"Initial potential energy: {state.getPotentialEnergy()}")
    print('Minimising ...')
    simulation.minimizeEnergy()

    # Write out the minimised PDB.
    with open(output_min, 'w') as outfile:
        PDBFile.writeFile(modeller.topology, context.getState(getPositions=True, enforcePeriodicBox=True).getPositions(), file=outfile, keepIds=True)
    replace_last_chain_identifier(output_min, "C")

    # Equilibrate
    simulation.context.setVelocitiesToTemperature(temperature_unit)
    print('Equilibrating ...')
    simulation.step(equilibration_steps)

    # Run the simulation.
    simulation.reporters.append(DCDReporter(output_traj_dcd, interval, enforcePeriodicBox=True))
    simulation.reporters.append(StateDataReporter(sys.stdout, interval * 5, step=True, potentialEnergy=True, temperature=True))
    
    print('Starting simulation with', steps, 'steps ...')
    t1 = time.time()
    simulation.step(steps)
    t2 = time.time()
    print('Simulation complete in {} mins at {}. Total wall clock time was {} mins'.format(
        round((t2 - t1) / 60, 3), temperature_unit, round((t2 - t0) / 60, 3)))
    print('Simulation time was', round(duration, 3), 'ns')


# Example usage:
# simulate_complex('openmm_simulations/ALLIANCE/reference/rece_1AKJ.pdb', 'openmm_simulations/ALLIANCE/reference/lig_1AKJ.pdb', steps=10000, solvate=True)