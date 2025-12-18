from multiprocessing import Pool, cpu_count

from rdkit import Chem


def converter(smiles):
    try:
        inchi = Chem.MolToInchi(Chem.MolFromSmiles(smiles))
    except:
        inchi = smiles

    return inchi


def smiles_to_inchi(smiles_list):
    with Pool(cpu_count()) as p:
        inchi_list = p.map(converter, smiles_list)

    return inchi_list
