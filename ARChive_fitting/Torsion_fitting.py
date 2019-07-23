# Fitting QCArchive TorsionDrives using QUBEKit
from QUBEKit.parametrisation import AnteChamber
from QUBEKit.utils.helpers import collect_archive_tdrive
from QUBEKit.dihedrals import TorsionOptimiser
from QUBEKit.ligand import Ligand

import qcportal as ptl

print('This script will optimise the torsion parameters of c1c[cH:1][c:2](cc1)[C:3](=[O:4])O using data from QCArchive')
# Open up an instance of the client
client = ptl.FractalClient()
client
# Gather the torsiondrive data set
ds = client.get_collection("TorsionDriveDataset", "OpenFF Fragmenter Phenyl Benchmark")

# Access the B3LYP-D3 data
ds.query("B3LYP-D3")

print('Getting TorsionDrive record')
# Extract the tdrive record to find the initial input molecule
td = ds.df.loc["c1c[cH:1][c:2](cc1)[C:3](=[O:4])O", "B3LYP-D3"]

# Get the molecule id that can be used to make a qubekit molecule
mol_id = td.initial_molecule
molecule = client.query_molecules(id=mol_id)[0]

print('Creating QUBEKit ligand molecule')
# Instance QUBEKit ligand from a json dict
mol = Ligand(molecule.json_dict(), 'torsion_example')

# Now we should reset the basis and theory to match
mol.basis = 'dzvp'
mol.theory = 'b3lyp-d3bj'

print('Parameterising using Anetchamber')
# Now check the structure of the molecule and assign initial parameters using antechamber
mol.write_pdb()
AnteChamber(mol)

# Get the central bond of the scanned dihedral
dihedral = td.keywords.dihedrals[0]
scan = (dihedral[1], dihedral[2])
# The ligand needs to know which torsion of the two rotatable options it has found is being fit
mol.scan_order = [scan]  # This is the central bond identified by td.keywords.dihedrals

# Now we need to extract all of the optimised geometries, we have a helper function that stores them into the
# required internal data structures
# Set the qm_scan data to a dict with keys coresponding to the central dihedral bond
mol.qm_scans = {}

print('Collecting torsiondrive results from the archive')
# Pass the torsiondrive record and the client instance to the helper function to get ordered lists of energies and geometries
mol.qm_scans[mol.scan_order[0]] = collect_archive_tdrive(td, client)

# Now try and do a quick single point fitting of the molecule, turn of the refinement method
mol.refinement_method = None
# Import the optimiser set the constraints to None and run the optimisation.
mol.constraints_file = None
optimiser = TorsionOptimiser(mol)
optimiser.refinement = None

print('Running optimisation')
optimiser.run()
