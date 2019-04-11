from rdkit.Chem import Descriptors

def get_rdkit_descriptors(m):
    """
    Return RDKit descriptors, from molecule object
    """

    return_dict = {}
    return_dict["ExactMolWt"] = Descriptors.ExactMolWt(m)
    return_dict["FpDensityMorgan1"] = Descriptors.FpDensityMorgan1(m)
    return_dict["FpDensityMorgan2"] = Descriptors.FpDensityMorgan2(m)
    return_dict["FpDensityMorgan3"] = Descriptors.FpDensityMorgan3(m)
    return_dict["HeavyAtomMolWt"] = Descriptors.HeavyAtomMolWt(m)
    return_dict["MaxAbsPartialCharge"] = Descriptors.MaxAbsPartialCharge(m)
    return_dict["MaxPartialCharge"] = Descriptors.MaxPartialCharge(m)
    return_dict["MinAbsPartialCharge"] = Descriptors.MinAbsPartialCharge(m)
    return_dict["MinPartialCharge"] = Descriptors.MinPartialCharge(m)
    return_dict["MolWt"] = Descriptors.MolWt(m)
    return_dict["NumRadicalElectrons"] = Descriptors.NumRadicalElectrons(m)
    return_dict["NumValenceElectrons"] = Descriptors.NumValenceElectrons(m)
    return return_dict
