from QUBEKit.ligand import Ligand
from QUBEKit.mod_seminario import ModSeminario
from QUBEKit.parametrisation.base_parametrisation import Parametrisation
from QUBEKit.utils import constants

import numpy as np
import qcportal as ptl


# Cache VEHICLe dataset
client = ptl.FractalClient()
client.list_collections('OptimizationDataset')
ds = client.get_collection('OptimizationDataset', 'OpenFF VEHICLe Set 1')
ds.df.head()
ds.list_specifications()

# Cache records
records = ds.data.records

for item in records:
    smiles = item.strip('\n')

    # Specific ID of given smiles string
    try:
        opt_record = ds.get_entry(smiles).object_map['default']
    except KeyError:
        continue

    # Optimisation of molecule at ID: opt_record
    optimisation = client.query_procedures(id=opt_record)[0]

    # Extract hessian
    try:
        hessian = client.query_results(molecule=optimisation.final_molecule, driver="hessian")[0].return_result
    except IndexError:
        # Molecule has been optimised but no hessian has been calculated yet
        continue

    # Reshape hessian
    conversion = constants.HA_TO_KCAL_P_MOL / (constants.BOHR_TO_ANGS ** 2)
    hessian = np.array(hessian).reshape(int(len(hessian) ** 0.5), -1) * conversion

    # Extract optimised structure
    opt_struct = client.query_procedures(id=opt_record)[0].get_final_molecule()

    # Initialise Ligand object using the json dict from qcengine
    mol = Ligand(opt_struct.json_dict(), name='initial_test')

    # Set the qm coords to the input coords from qcengine
    mol.coords['qm'] = mol.coords['input']

    # Insert hessian and optimised coordinates
    mol.hessian = hessian
    mol.parameter_engine = 'none'

    # Create empty parameter dicts
    Parametrisation(mol).gather_parameters()

    print(item)

    with open('Modified_Seminario_Bonds.txt', 'a+') as bonds_file, \
            open('Modified_Seminario_Angles.txt', 'a+') as angles_file:
        bonds_file.write(f'\n\n{smiles}\n\n')
        angles_file.write(f'\n\n{smiles}\n\n')

    # Get Mod Sem angle and bond params
    ModSeminario(mol).modified_seminario_method()
