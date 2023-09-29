from Bio.PDB import PDBList, PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.PDB import DSSP
import freesasa
import re

def get_pdb_files(ids):
    pdbl = PDBList()
    for id in ids:
        pdbl.retrieve_pdb_file(id, file_format="pdb")

def get_pdb_sequences(ids):
    parser = PDBParser()
    ppb = PPBuilder(radius=1.8) #CaPP default value
    result = {key: "" for key in ids}
    for id in ids:
        structure = parser.get_structure(id, f"{id[1:3]}/pdb{id}.ent")
        
        pps = ppb.build_peptides(structure)
        regex = re.compile(r"start=(\d+) end=(\d+)")
        polypeptides = [(int(start), int(end)) for start, end in regex.findall(pps.__repr__())]
        i = 0
        for pp in pps:
            if polypeptides[0][0] > 0:
                result[id] += "?" * polypeptides[0][0]
            result[id] += str(pp.get_sequence())
            if i < len(polypeptides) - 1:
                gap_size = polypeptides[i + 1][0] - polypeptides[i][1] - 1
                if gap_size > 0:
                    result[id] += "?" * gap_size
            i += 1

    return result

ids = ["4esq","2h34"]
get_pdb_files(ids)
result = get_pdb_sequences(ids)
print(result, sep="\n")

features = {key: [] for key in ids}
for res in result:
    analysis = ProteinAnalysis(res.value)
    features[res.key].molecular_weight = analysis.molecular_weight()
    features[res.key].isoelectric_point = analysis.isoelectric_point()
    features[res.key].aromaticity = analysis.aromaticity()
    features[res.key].instability_index = analysis.instability_index()
    features[res.key].secondary_structure_fraction = analysis.secondary_structure_fraction()
    

# Surface area can be calculated using third-party libraries like freesasa or PySASA

# import freesasa

# def surface_area(structure):
#     structure = freesasa.Structure('data/raw/1tup.pdb')
#     result = freesasa.calc(structure)
#     return result.totalArea()

# # Example usage:
# sa = surface_area(structure)
# print(f"Surface Area: {sa}")

# from Bio.PDB import DSSP

# def solvent_accessibility(structure):
#     model = structure[0]
#     dssp = DSSP(model, 'data/raw/1tup.pdb')  # Assuming 1tup.pdb is the protein structure file
#     solvent_accessibility = [(residue_id, rel_ASA) for residue_id, (res, ss, rel_ASA, phi, psi) in dssp]
#     return solvent_accessibility

# # Example usage:
# structure = parser.get_structure('Protein_Name', 'data/raw/1tup.pdb')
# sa = solvent_accessibility(structure)
# print(f"Solvent Accessibility: {sa}")
