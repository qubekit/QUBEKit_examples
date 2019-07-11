import qcportal as ptl
import numpy as np

from QUBEKit.ligand import Ligand
from QUBEKit.mod_seminario import ModSeminario
from QUBEKit.parametrisation.base_parametrisation import Parametrisation


# Cache VEHICLe dataset
client = ptl.FractalClient()
client.list_collections('OptimizationDataset')
ds = client.get_collection('OptimizationDataset', 'OpenFF VEHICLe Set 1')
ds.df.head()
ds.list_specifications()

# Abstracted out for later use in a loop
smiles = 'c1cc[nH]c1'

# Specific ID of given smiles string
opt_record = ds.get_entry(f'{smiles}-0').object_map['default']

# Optimisation of molecule at ID: opt_record
optimisation = client.query_procedures(id=opt_record)[0]

# Extract hessian
hessian = client.query_results(molecule=optimisation.final_molecule, driver="hessian")[0].return_result
# Reshape hessian
hessian = np.array(hessian).reshape(int(len(hessian) ** 0.5), -1) * 627.509391 / (0.529 ** 2)

# Extract optimised structure
opt_struct = client.query_procedures(id=opt_record)[0].get_final_molecule()

# Initialise Ligand object using the smiles string
mol = Ligand(opt_struct.json_dict(), name='initial_test')

# Set the qm coords to the input coords from qcengine
mol.coords['qm'] = mol.coords['input']

# Insert hessian and optimised coordinates
mol.hessian = hessian
mol.parameter_engine = 'none'

# Create empty parameter dicts
Parametrisation(mol).gather_parameters()

# Get Mod Sem angle and bond params
ModSeminario(mol).modified_seminario_method()

# Write out final xml
mol.write_parameters()
