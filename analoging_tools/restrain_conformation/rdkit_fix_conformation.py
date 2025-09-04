from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms
from rdkit.Chem import rdFMCS 
import argparse
import os
import subprocess

# Brandon L

def write_mol2(mol, path):
    from rdkit.Chem import rdmolfiles
    writer = Chem.SDWriter(path)
    writer.write(mol)
    writer.close()

def read_smi_file(filename):
    smiles_list = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            smiles = parts[0]
            name = parts[1] if len(parts) > 1 else ""
            smiles_list.append((smiles, name))
    return smiles_list

def find_custom_pattern(mol, smarts):
    patt = Chem.MolFromSmarts(smarts)
    matches = mol.GetSubstructMatches(patt)
    if not matches:
        return None
    return matches[0]

def constrained_embed_from_reference(mol, ref_mol, smarts):
    ref_match = find_custom_pattern(ref_mol, smarts)
    query_match = find_custom_pattern(mol, smarts)
    if ref_match is None or query_match is None:
        print("SMARTS did not match in reference or target molecule.")
        return None, None
    ref_conf = ref_mol.GetConformer()
    coordMap = {}
    for qi, ri in zip(query_match, ref_match):  # All motif atoms
        coordMap[qi] = ref_conf.GetAtomPosition(ri)
    try:
        cid = AllChem.EmbedMolecule(mol, coordMap=coordMap, useRandomCoords=True, maxAttempts=1000)
        if cid != 0:
            print("Embedding failed for motif-constrained molecule.")
            return None, None
    except Exception as e:
        print(f"Embedding failed: {e}")
        return None, None
    return mol, query_match

'''
def set_dihedral_by_indices(mol, atom_indices, atom_order, angle=0):
    conf = mol.GetConformer()
    idxs = [atom_indices[i] for i in atom_order]
    rdMolTransforms.SetDihedralDeg(conf, *idxs, angle)
'''

def process_smiles_with_full_constraint(smiles_and_names, ref_mol, smarts, verbose=True):
    results = []
    for smi, name in smiles_and_names:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            if verbose:
                print(f"Invalid SMILES: {smi} ({name})")
            results.append(None)
            continue
        mol = Chem.AddHs(mol)
        mol, query_match = constrained_embed_from_reference(mol, ref_mol, smarts)
        if mol is None:
            results.append(None)
            continue
        # MMFF optimize with strong constraints on all matched atoms
        mp = AllChem.MMFFGetMoleculeProperties(mol)
        ff = AllChem.MMFFGetMoleculeForceField(mol, mp)
        for i in query_match:
            ff.MMFFAddPositionConstraint(i, 0, 1.e4)
        ff.Minimize(maxIts=100000)
        if verbose:
            print(f"Processed: {smi} ({name})")
        if name:
            mol.SetProp("_Name", name)
        results.append(mol)
    return results

def load_reference(ref_file):
    ext = os.path.splitext(ref_file)[1].lower()
    if ext == ".mol2":
        mol = Chem.MolFromMol2File(ref_file, removeHs=False)
    elif ext in [".pdb", ".ent"]:
        mol = Chem.MolFromPDBFile(ref_file, removeHs=False)
    elif ext in [".mol", ".sdf"]:
        mol = Chem.MolFromMolFile(ref_file, removeHs=False)
    else:
        raise ValueError(f"Unknown file extension: {ext}")
    if mol is None or mol.GetNumConformers() == 0:
        raise ValueError("Reference structure could not be loaded or has no 3D coordinates.")
    return mol


def main():
    parser = argparse.ArgumentParser(
        description='Generate 3D conformers from smiles with substructure constrained to conformation of refernce ligand'
    )

    parser.add_argument('smiles', type=str, help='*smi files containing smiles of close analogues to reference ligand')
    parser.add_argument('reference_ligand', type=str, help='3D structure of reference ligand')
    parser.add_argument('smarts', type=str, help='SMARTS pattern used to constrain conformation of analogues to reference ligand conformation')

    args = parser.parse_args()

    smi_filename = args.smiles
    ref_filename = args.reference_ligand
    motif_smarts = args.smarts

    smiles_and_names = read_smi_file(smi_filename)
    ref_mol = load_reference(ref_filename)

    mols = process_smiles_with_full_constraint(
        smiles_and_names,
        ref_mol,
        motif_smarts,
        verbose=True
    )
    output_dir = './fix_confirmation_output'
    output_name = 'fix_confirmation_output'
    os.makedirs(output_dir, exist_ok=True)

    for idx, m in enumerate(mols):
        if m is None:
            continue
        name = m.GetProp("_Name") if m.HasProp("_Name") and m.GetProp("_Name") else f"mol_{idx+1}"
        # Sanitize name for files
        name = "".join([c if c.isalnum() or c in "-_." else "_" for c in name])
        mol2_path = os.path.join(output_dir, f"{name}.sdf")
        write_mol2(m, mol2_path)

    full_cmd = (
    "source /lustre/fs6/lyu_lab/scratch/lyu_soft//openbabel/current/env.sh && "
    "obabel *.sdf -omol2 -m && "
    "rm *.sdf"
    )

    subprocess.run(
        full_cmd,
        shell=True,
        executable="/bin/bash",
        cwd= os.path.join(os.getcwd(), output_name),
        check=True
    )

if __name__ == '__main__':
    main()





    
    



    

