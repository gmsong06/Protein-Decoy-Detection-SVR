import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from Bio.PDB import PDBParser
import numpy as np
from scipy.spatial import KDTree

def get_contacts(mol, cutoff=5.0):
    conf = mol.GetConformer()
    num_atoms = mol.GetNumAtoms()

    # List to store indices and positions of non-hydrogen atoms
    non_hydrogen_indices = []
    non_hydrogen_positions = []

    # List to store indices and positions of hydrogen atoms
    hydrogen_indices = []
    hydrogen_positions = []

    # Separate non-hydrogen and hydrogen atoms
    for i in range(num_atoms):
        atom_i = mol.GetAtomWithIdx(i)
        pos_i = conf.GetAtomPosition(i)
        
        if atom_i.GetSymbol() in {'H', 'h'}:
            hydrogen_indices.append(i)
            hydrogen_positions.append(pos_i)
        else:
            non_hydrogen_indices.append(i)
            non_hydrogen_positions.append(pos_i)
    
    non_hydrogen_positions = np.array(non_hydrogen_positions)
    tree = KDTree(non_hydrogen_positions)
    contacts = set()

    # Check contacts between non-hydrogen atoms
    for i, pos_i in zip(non_hydrogen_indices, non_hydrogen_positions):
        indices = tree.query_ball_point(pos_i, cutoff)
        for idx in indices:
            j = non_hydrogen_indices[idx]
            if i < j:
                contacts.add((i, j))
    
    
    return list(contacts)


def get_chain_info(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_file)

    chains = {}
    for model in structure:
        for chain in model:
            chain_id = chain.id
            chains[chain_id] = []

            for residue in chain:
                for atom in residue:
                    chains[chain_id].append(atom.coord)
    
    return chains

def test_mmff94_charges(pdb_file):
    try:
        mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
        if mol is None:
            raise ValueError("Molecule could not be read from the PDB file.")
    except Exception as e:
        print(f"Error reading PDB file: {e}", flush=True)
        return

    mol = Chem.AddHs(mol)
    print(f"Number of atoms after adding Hs: {mol.GetNumAtoms()}", flush=True)

    optimize_status = AllChem.MMFFOptimizeMolecule(mol, maxIters=200, nonBondedThresh=100.0)
    if optimize_status == 0:
        print("Molecule optimization successful", flush=True)
    else:
        print("Molecule optimization failed", flush=True)

    try:
        mmff_props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94')
        if mmff_props is None:
            raise ValueError("MMFF properties could not be retrieved.")
    except Exception as e:
        print(f"Error getting MMFF properties: {e}", flush=True)
        return

    chains = get_chain_info(pdb_file)
    
    contacts = get_contacts(mol, cutoff=5.0)
    print(f"Number of contacts: {len(contacts)}", flush=True)

    contact_charges = []
    electric_potential = 0
    conf = mol.GetConformer()
    
    print("Contacts (indexes) and Partial Charges of Atoms in Contact:", flush=True)
    for (i, j) in contacts:
        charge_i = mmff_props.GetMMFFPartialCharge(i)
        charge_j = mmff_props.GetMMFFPartialCharge(j)
        pos_i = conf.GetAtomPosition(i)
        pos_j = conf.GetAtomPosition(j)
        distance = pos_i.Distance(pos_j)
        electric_potential += (charge_i * charge_j) / distance
        contact_charges.append((i, j, charge_i, charge_j, distance))
    
    for (i, j, charge_i, charge_j, distance) in contact_charges:
        print(f"Atoms in contact: Index i = {i}, Index j = {j}, Atom i Partial Charge: {charge_i}, Atom j Partial Charge: {charge_j}, Distance: {distance}", flush=True)

    print(f"Total Electric Potential: {electric_potential}", flush=True)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python test_mmff94_charges.py <pdb_file>", flush=True)
        sys.exit(1)

    pdb_file = sys.argv[1]
    test_mmff94_charges(pdb_file)



