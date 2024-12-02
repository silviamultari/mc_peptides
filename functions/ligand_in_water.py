from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
from pdbfixer import PDBFixer
from io import StringIO
from simtk.openmm.app import PDBFile, ForceField, Simulation, PDBReporter, StateDataReporter
from simtk.openmm import LangevinIntegrator
from simtk.unit import kelvin, picosecond, nanometer
from simtk.openmm.app import NoCutoff, AllBonds
from pdbfixer import PDBFixer
from io import StringIO
import os

def run_simulation(input_pdb, output_pdb, steps, report_interval, topology_file='ref_topology.pdb'):
    # Load the PDB file and fix it using PDBFixer
    pdb_fixer = PDBFixer(filename=input_pdb)
    pdb_fixer.findMissingResidues()
    pdb_fixer.findMissingAtoms()
    pdb_fixer.addMissingAtoms()
    pdb_fixer.addMissingHydrogens(7.4)

    # Save the fixed PDB to a string buffer
    fixed_pdb = StringIO()
    PDBFile.writeFile(pdb_fixer.topology, pdb_fixer.positions, fixed_pdb, keepIds=True)
    fixed_pdb.seek(0)

    # Load the fixed PDB into an OpenMM PDBFile object
    pdb = PDBFile(fixed_pdb)

    # Define the force field
    forcefield = ForceField('amber99sb.xml', 'tip3p.xml')

    # Create the system
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff, nonbondedCutoff=1*nanometer, constraints=AllBonds)

    # Create an integrator
    integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.0025*picoseconds)

    # Create the simulation object
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)

    # Minimize energy
    simulation.minimizeEnergy()

    # Save the topology file
    with open(topology_file, 'w') as f:
        PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), f)

    # Set up reporters
    simulation.reporters.append(PDBReporter(output_pdb, report_interval))
    simulation.reporters.append(StateDataReporter(stdout, report_interval, step=True, potentialEnergy=True, temperature=True))

    # Run the simulation
    simulation.step(steps)


