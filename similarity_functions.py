import csv
import operator
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import Descriptors
from rdkit.Chem.Fingerprints import FingerprintMols
from sanifix import fix_mol


def similar_FDA_molecules(m):
    """
    Search known FDA-approved molecules for those similar to the provided SMILES string.
    """
    query_fp = FingerprintMols.FingerprintMol(m)
    fda_dict = {}
    csv_reader = csv.reader(open("./resources/FDAdrugs_name_smiles.csv"))
    for row in csv_reader:
        fda_dict[row[1]] = row[0]
    fda_list = list(fda_dict.keys())
    fda_ms = [Chem.MolFromSmiles(x) for x in fda_list]
    fda_fps = [FingerprintMols.FingerprintMol(x) for x in fda_ms]
    similarity_dict = {}
    for idx, fda_fp in enumerate(fda_fps):
        similarity = DataStructs.FingerprintSimilarity(fda_fp, query_fp)
        similarity_dict[fda_list[idx]] = similarity
    most_similar_smiles = max(similarity_dict.items(), key=operator.itemgetter(1))[0]
    most_similar_name = fda_dict[most_similar_smiles]
    return_dict = {'FDA smiles': most_similar_smiles, 'name': most_similar_name}
    return return_dict


