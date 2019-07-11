import qcportal as ptl
import numpy as np

from QUBEKit.ligand import Ligand
from QUBEKit.engines.rdkit import RDKit


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
opt_struct = client.query_results(molecule=optimisation.final_molecule, driver='gradient')[0].return_result
# Reshape optimised structure
opt_struct = np.array(opt_struct).reshape((len(opt_struct) // 3, 3)) * 0.529

# Initialise Ligand object using the smiles string
mol = Ligand(RDKit().smiles_to_pdb(smiles))

# Insert hessian and optimised coordinates
mol.hessian = hessian
mol.coords['qm'] = opt_struct

mol.write_pdb(input_type='qm', name='test')
