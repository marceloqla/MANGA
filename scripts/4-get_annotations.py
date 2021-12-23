import os
import sys
import re
import math
import json
import textwrap
import itertools
import Bio
import pickle
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist
import urllib.request
from bs4 import BeautifulSoup
import time
from Bio.PDB.PDBList import PDBList
from Bio.Data import SCOPData
from Bio.PDB.DSSP import DSSP
from prody import parsePDB, ANM, calcPerturbResponse, calcMechStiff, sliceAtomicData, calcSqFlucts#, calcDynamicCouplingIndex, 
import numpy as np
import networkx as nx

ONE_AA_TYPES = [
    "A","C","D","E","F",
    "G","H","I","K","L",
    "M","N","P","Q","R",
    "S","T","V","W","Y"
]

THREE_TO_ONE_AA_TYPES_CANONICAL = {
    "ALA": "A", "MET": "M",
    "CYS": "C", "ASN": "N",
    "ASP": "D", "PRO": "P",
    "GLU": "E", "GLN": "Q",
    "PHE": "F", "ARG": "R",
    "GLY": "G", "SER": "S",
    "HIS": "H", "THR": "T",
    "ILE": "I", "VAL": "V",
    "LYS": "K", "TRP": "W",
    "LEU": "L", "TYR": "Y"
}

THREE_AA_TYPES_CANONICAL_SMALL = list(THREE_TO_ONE_AA_TYPES_CANONICAL.keys())
THREE_AA_TYPES_CANONICAL_SMALL = [aa[0]+aa[1:].lower() for aa in THREE_AA_TYPES_CANONICAL_SMALL]
THREE_TO_ONE_AA_TYPES_CANONICAL_SMALL = dict((aathree[0]+aathree[1:].lower(), aaone) for aathree, aaone in THREE_TO_ONE_AA_TYPES_CANONICAL.items())

ONE_TO_THREE1 = {
    "A":"Ala", "M":"Met",
    "C":"Cys", "N":"Asn",
    "D":"Asp", "P":"Pro",
    "E":"Glu", "Q":"Gln",
    "F":"Phe", "R":"Arg",
    "G":"Gly", "S":"Ser",
    "H":"His", "T":"Thr",
    "I":"Ile", "V":"Val",
    "K":"Lys", "W":"Trp",
    "L":"Leu", "Y":"Tyr"
}

class JSONObject:
    def toJSON(self):
        return json.dumps(self, default=lambda o: o.__dict__, 
            sort_keys=True, indent=4)
    def __str__(self):
        # for k, v in self.__dict__:
        #     return f"{k}, {v}"
        return str(self.__dict__)

def pickle_writer(fpath, fobj):
    fopenp = open(fpath, 'wb')
    pickle.dump(fobj, fopenp, protocol=3)
    fopenp.close()

def pickle_loader(fpath):
    fopenp = open(fpath, 'rb')
    fobj = pickle.load(fopenp)
    fopenp.close()
    return fobj

def make_output_dir(outputdir):
    if outputdir[len(outputdir)-1] == '/':#Remove /
        outputdir = outputdir[:-1]
    if not os.path.exists(outputdir):
        # print(f"{outputdir} is being created.")
        os.makedirs(outputdir)
    # print(f"{outputdir} folder existance checked.")
    return outputdir

def get_bad_ligand_list(ions=False):
    bleicher = ["SO4","GOL","PO4","EDO","ACT","BME","TRS","MPD","ACY","MES","CIT","FMT","EPE","ACE","PG4","DMS","IPA","NO3","CO3","CAC","BOG","IPH","IMD","AZI","FLC","DTT","1PE","PGA","NH4","EOH","P6G","WO4","SIN","BCT","BTB","PO3","SO3","P33","OAA","15P","PE4","NCO","1PG","12P","MOH","CH3","PE5","PEU","MO6","P4G","PIG","P2K"]
    extension_mqla = ["UNX", "OLA", "OLB", "OLC", "TLA"]
    all_bads = bleicher.copy()
    all_bads.extend(extension_mqla)
    if ions:
        # extension_ions = ["CA", "NA", "ZN"]
        extension_ions = ["NA"]
        all_bads.extend(extension_ions)
    return all_bads

def parseQualifier(new_qualifier):
    if not new_qualifier:
        return tuple()
    if "/" in new_qualifier[0]:
        qualifier_type = new_qualifier[0].split("/")[1].split("=")[0]
        qualifier_split = f"/{qualifier_type}="
        first_line = new_qualifier[0].split(qualifier_split)[1].rstrip()
    else:
        qualifier_type = "note"
        if "{" in qualifier_type:
            qualifier_type = "evidence"
        first_line = new_qualifier[0].split()[1:]
    qualifier_lines_text = ""
    if len(new_qualifier) > 1:
        qualifier_lines = [" ".join(line.split()[1:]) for line in new_qualifier[1:]]
        qualifier_lines_text = " ".join([line.rstrip() for line in qualifier_lines])
    qualifier_lines_text = first_line + " " + qualifier_lines_text
    return (qualifier_type, qualifier_lines_text)

def parseQualifiers(qualifier_lines):
    feature_qualifiers = []
    new_qualifier = None
    for feature_line in qualifier_lines:
        if feature_line.split()[1][0] == "/":
            if new_qualifier:
                feature_qualifiers.append(parseQualifier(new_qualifier))
            new_qualifier = []
        new_qualifier.append(feature_line)
    feature_qualifiers.append(parseQualifier(new_qualifier))
    return feature_qualifiers

def retrieve_uniprot_features(uniprot_code, uniprot_file):
    feature_type_alt_list = [
        "SITE",
        "MOD_RES",
        "CROSSLNK",
        "CARBOHYD",
        "LIPID",
        "MOTIF",
    ]
    feature_type_keep_list = [ #https://www.uniprot.org/docs/userman.htm
        "DNA_BIND",
        "BINDING",
        "DISULFID",
        "ACT_SITE",
        "NP_BIND",
    ]
    content = ""
    try:
        f = open(uniprot_file, "r")
        content = f.read()
        # Do something with the file
        f.close()
    except FileNotFoundError:
        download_link = f"https://www.uniprot.org/uniprot/{uniprot_code}.txt"
        handle = urllib.request.urlopen(download_link)
        content = handle.read().decode("utf-8")
        f = open(uniprot_file, "w")
        f.write(""+content)
        f.close()
        # print("content")
        # print(content)
    ft_lines = [line for line in content.split('\n') if line[:2] == "FT"]
    each_features = []
    current_feature = None
    for i_ft_line, ft_line in enumerate(ft_lines):
        if ft_line[5] != " ":
            if current_feature: 
                each_features.append(current_feature)
            current_feature = []
        current_feature.append(ft_line)
    each_features.append(current_feature)

    each_features_transformed = []
    alt_features_transformed = []

    for feature in each_features:
        feature_type = feature[0].split()[1]
        feature_resids = []
        feature_qualifiers = []
        if feature_type in feature_type_keep_list:
            qualifier_lines = []
            if feature_type == "BINDING":
                if len(feature[0].split()) == 3:
                    # feature_resids = re.findall('\d+', ft_resids_line)
                    feature_resids.append(feature[0].split()[-1])
                elif len(feature[0].split()) == 4:
                    feature_resids.append(feature[0].split()[2])
                    feature_resids.append(feature[0].split()[3])
                elif len(feature[0].split()) == 5:
                    feature_resids.append(feature[0].split()[2])
                    feature_resids.append(feature[0].split()[3])
                    qualifier_lines.append(feature[0].split()[-1])
                else:
                    feature_resids = []
                    # print(feature)
                    # quit()
            elif feature_type == "DNA_BIND":
                if len(feature[0].split()) == 3:
                    if ".." in feature[0].split()[-1]:
                        feature_resids = re.findall('\d+', feature[0].split()[-1])
                        if len(feature_resids) > 1:
                            feature_resids = list(range(int(feature_resids[0]),int(feature_resids[1])+1))
                    else:
                        feature_resids = []
                        # print(feature)
                        # quit()
                    # feature_resids.append(feature[0].split()[-1])
                elif len(feature[0].split()) == 4:
                    feature_resids.append(feature[0].split()[2])
                    feature_resids.append(feature[0].split()[3])
                    # feature_resids = [re.search('\d+', ftres).group() for ftres in feature_resids]
                    feature_resids = [re.search('\d+', ftres) for ftres in feature_resids]
                    feature_resids = [ftres.group() for ftres in feature_resids if ftres]
                    if len(feature_resids) > 1:
                        feature_resids = list(range(int(feature_resids[0]),int(feature_resids[1])+1))
                elif len(feature[0].split()) == 5:
                    feature_resids.append(feature[0].split()[2])
                    feature_resids.append(feature[0].split()[3])
                    # feature_resids = [re.search('\d+', ftres).group() for ftres in feature_resids]
                    feature_resids = [re.search('\d+', ftres) for ftres in feature_resids]
                    feature_resids = [ftres.group() for ftres in feature_resids if ftres]
                    if len(feature_resids) > 1:
                        feature_resids = list(range(int(feature_resids[0]),int(feature_resids[1])+1))
                    qualifier_lines.append(feature[0].split()[-1])
                else:
                    feature_resids = []
                    # print(feature)
                    # quit()
            elif feature_type == "DISULFID":
                if len(feature[0].split()) == 3:
                    if ".." in feature[0].split()[-1]:
                        feature_resids = re.findall('\d+', feature[0].split()[-1])
                    else:
                        feature_resids = []
                        # print(feature)
                        # quit()
                    # feature_resids.append(feature[0].split()[-1])
                elif len(feature[0].split()) == 4:
                    feature_resids.append(feature[0].split()[2])
                    feature_resids.append(feature[0].split()[3])
                elif len(feature[0].split()) == 5:
                    feature_resids.append(feature[0].split()[2])
                    feature_resids.append(feature[0].split()[3])
                    qualifier_lines.append(feature[0].split()[-1])
                else:
                    feature_resids = []
                    # print(feature)
                    # quit()
            elif feature_type == "ACT_SITE":
                if len(feature[0].split()) == 3:
                    feature_resids.append(feature[0].split()[-1])
                elif len(feature[0].split()) == 4:
                    feature_resids.append(feature[0].split()[2])
                    feature_resids.append(feature[0].split()[3])
                elif len(feature[0].split()) == 5:
                    feature_resids.append(feature[0].split()[2])
                    feature_resids.append(feature[0].split()[3])
                    qualifier_lines.append(feature[0].split()[-1])
                else:
                    feature_resids = []
                    # print(feature)
                    # quit()
            elif feature_type == "NP_BIND":
                if len(feature[0].split()) == 3:
                    if ".." in feature[0].split()[-1]:
                        feature_resids = re.findall('\d+', feature[0].split()[-1])
                        if len(feature_resids) > 1:
                            feature_resids = list(range(int(feature_resids[0]),int(feature_resids[1])+1))
                    else:
                        feature_resids = []
                        # print(feature)
                        # quit()
                elif len(feature[0].split()) == 4:
                    feature_resids.append(feature[0].split()[2])
                    feature_resids.append(feature[0].split()[3])
                    # feature_resids = [re.search('\d+', ftres).group() for ftres in feature_resids]
                    feature_resids = [re.search('\d+', ftres) for ftres in feature_resids]
                    feature_resids = [ftres.group() for ftres in feature_resids if ftres]
                    if len(feature_resids) > 1:
                        feature_resids = list(range(int(feature_resids[0]),int(feature_resids[1])+1))
                elif len(feature[0].split()) == 5:
                    feature_resids.append(feature[0].split()[2])
                    feature_resids.append(feature[0].split()[3])
                    # feature_resids = [re.search('\d+', ftres).group() for ftres in feature_resids]
                    feature_resids = [re.search('\d+', ftres) for ftres in feature_resids]
                    feature_resids = [ftres.group() for ftres in feature_resids if ftres]
                    if len(feature_resids) > 1:
                        feature_resids = list(range(int(feature_resids[0]),int(feature_resids[1])+1))
                    qualifier_lines.append(feature[0].split()[-1])
                else:
                    feature_resids = []
                    # print(feature)
                    # quit()
            if len(feature) > 1:
                qualifier_lines.extend(feature[1:])
            feature_qualifiers = parseQualifiers(qualifier_lines)
            # feature_resids = [re.search('\d+', str(ftres)).group() for ftres in feature_resids]
            feature_resids = [re.search('\d+', str(ftres)) for ftres in feature_resids]
            feature_resids = [ftres.group() for ftres in feature_resids if ftres]
            feature_resids = [int(ftres) for ftres in feature_resids if ftres]
            feature_resids = list(set(feature_resids))
            if len(feature_resids) > 0:
                each_features_transformed.append((feature_type,feature_resids,feature_qualifiers))
        if feature_type in feature_type_alt_list:
            qualifier_lines = []
            #TODO: fillup uniprot parsing and test
            if feature_type == "SITE":
                if len(feature[0].split()) == 3:
                    feature_resids = re.findall('\d+', feature[0].split()[-1])
                elif len(feature[0].split()) == 4:
                    feature_resids.append(feature[0].split()[2])
                    feature_resids.append(feature[0].split()[3])
                elif len(feature[0].split()) == 5:
                    feature_resids.append(feature[0].split()[2])
                    feature_resids.append(feature[0].split()[3])
                    qualifier_lines.append(feature[0].split()[-1])
                else:
                    feature_resids = []
                    # print(feature)
                    # quit()
            if feature_type == "MOD_RES":
                if len(feature[0].split()) == 3:
                    feature_resids.append(feature[0].split()[-1])
                elif len(feature[0].split()) == 4:
                    feature_resids.append(feature[0].split()[2])
                    feature_resids.append(feature[0].split()[3])
                elif len(feature[0].split()) == 5:
                    feature_resids.append(feature[0].split()[2])
                    feature_resids.append(feature[0].split()[3])
                    qualifier_lines.append(feature[0].split()[-1])
                else:
                    feature_resids = []
                    # print(feature)
                    # quit()
            if feature_type == "CROSSLNK":
                if len(feature[0].split()) == 3:
                    feature_resids = re.findall('\d+', feature[0].split()[-1])
                elif len(feature[0].split()) == 4:
                    feature_resids.append(feature[0].split()[2])
                    feature_resids.append(feature[0].split()[3])
                elif len(feature[0].split()) == 5:
                    feature_resids.append(feature[0].split()[2])
                    feature_resids.append(feature[0].split()[3])
                    qualifier_lines.append(feature[0].split()[-1])
                else:
                    feature_resids = []
                    # print(feature)
                    # quit()
            if feature_type == "CARBOHYD":
                if len(feature[0].split()) == 3:
                    feature_resids.append(feature[0].split()[-1])
                elif len(feature[0].split()) == 4:
                    feature_resids.append(feature[0].split()[2])
                    feature_resids.append(feature[0].split()[3])
                elif len(feature[0].split()) == 5:
                    feature_resids.append(feature[0].split()[2])
                    feature_resids.append(feature[0].split()[3])
                    qualifier_lines.append(feature[0].split()[-1])
                else:
                    feature_resids = []
                    # print(feature)
                    # quit()
            if feature_type == "LIPID":
                if len(feature[0].split()) == 3:
                    feature_resids.append(feature[0].split()[-1])
                elif len(feature[0].split()) == 4:
                    feature_resids.append(feature[0].split()[2])
                    feature_resids.append(feature[0].split()[3])
                elif len(feature[0].split()) == 5:
                    feature_resids.append(feature[0].split()[2])
                    feature_resids.append(feature[0].split()[3])
                    qualifier_lines.append(feature[0].split()[-1])
                else:
                    feature_resids = []
                    # print(feature)
                    # quit()
            if feature_type == "MOTIF":
                if len(feature[0].split()) == 3:
                    if ".." in feature[0].split()[-1]:
                        feature_resids = re.findall('\d+', feature[0].split()[-1])
                        if len(feature_resids) > 1:
                            feature_resids = list(range(int(feature_resids[0]),int(feature_resids[1])+1))
                    else:
                        feature_resids = []
                        # print(feature)
                        # quit()
                elif len(feature[0].split()) == 4:
                    feature_resids.append(feature[0].split()[2])
                    feature_resids.append(feature[0].split()[3])
                    # feature_resids = [re.search('\d+', ftres).group() for ftres in feature_resids]
                    feature_resids = [re.search('\d+', ftres) for ftres in feature_resids]
                    feature_resids = [ftres.group() for ftres in feature_resids if ftres]
                    if len(feature_resids) > 1:
                        feature_resids = list(range(int(feature_resids[0]),int(feature_resids[1])+1))
                elif len(feature[0].split()) == 5:
                    feature_resids.append(feature[0].split()[2])
                    feature_resids.append(feature[0].split()[3])
                    # feature_resids = [re.search('\d+', ftres).group() for ftres in feature_resids]
                    feature_resids = [re.search('\d+', ftres) for ftres in feature_resids]
                    feature_resids = [ftres.group() for ftres in feature_resids if ftres]
                    if len(feature_resids) > 1:
                        feature_resids = list(range(int(feature_resids[0]),int(feature_resids[1])+1))
                    qualifier_lines.append(feature[0].split()[-1])
                else:
                    feature_resids = []
                    # print(feature)
                    # quit()
            if len(feature) > 1:
                qualifier_lines.extend(feature[1:])
            feature_qualifiers = parseQualifiers(qualifier_lines)
            # feature_resids = [re.search('\d+', str(ftres)).group() for ftres in feature_resids]
            feature_resids = [re.search('\d+', str(ftres)) for ftres in feature_resids]
            feature_resids = [ftres.group() for ftres in feature_resids if ftres]
            feature_resids = [int(ftres) for ftres in feature_resids if ftres]
            feature_resids = list(set(feature_resids))
            if len(feature_resids) > 0:
                alt_features_transformed.append((feature_type,feature_resids,feature_qualifiers))
    return each_features_transformed, alt_features_transformed

def get_sequence_from_uniprot_code(uniprot_code):
    seq = ""
    download_link = f"https://www.uniprot.org/uniprot/{uniprot_code}.fasta"
    handle = urllib.request.urlopen(download_link)
    content = handle.read().decode("utf-8")
    seq = "".join(content.split("\n")[1:])
    make_output_dir("./seqfasta/")
    with open(f"./seqfasta/{uniprot_code}.fasta", "w") as fasta_handle:
        fasta_string = '>'+ uniprot_code+'\n'+'\n'.join(textwrap.wrap(seq, 60))+''
        fasta_handle.write(fasta_string)
    return seq

def retrieve_pdb_seq(pdb_code, selchain="", download_dir="./PDB", verbose=False):
    pdbl = PDBList(verbose=verbose)
    pdbl.retrieve_pdb_file(pdb_code,pdir=download_dir,file_format="pdb")
    
    pdb_file_name = "pdb" + pdb_code.lower() + ".ent"
    structure = None
    try:
        structure = Bio.PDB.PDBParser().get_structure(pdb_code, download_dir + "/" + pdb_file_name)
    except FileNotFoundError:
        pdbl.retrieve_pdb_file(pdb_code,pdir=download_dir,file_format="mmCif")
        pdb_file_name = pdb_code.lower() + ".cif"
        structure = Bio.PDB.MMCIFParser().get_structure(pdb_code, download_dir + "/" + pdb_file_name)
    modelo = structure[0]
    chain_to_seqs = {}
    for chain in modelo:
        chain_to_seqs[chain.id] = []

        pdbsequence = ""
        pdbsequencenum = []
        pdbsequencefull = []
        
        for residue in chain:
            if Bio.PDB.Polypeptide.is_aa(residue) == True:
                number = residue.get_id()
                oneletter = ""
                try:
                    oneletter = Bio.PDB.Polypeptide.three_to_one(residue.get_resname())
                except KeyError:
                    oneletter = SCOPData.protein_letters_3to1[residue.get_resname()]
                try:
                    num = str(number[1]) + str((number[2].rstrip( ))[0])
                except IndexError:
                    num = str(number[1])
                if len(oneletter) > 1:
                    letters = list(oneletter)
                    for letter in letters:
                        pdbsequence += letter
                        pdbsequencenum.append(residue.get_resname() + num)
                        pdbsequencefull.append(residue.get_id())
                    continue
                else:
                    pdbsequence += oneletter
                    pdbsequencenum.append(residue.get_resname() + num)
                    pdbsequencefull.append(residue.get_id())
        chain_to_seqs[chain.id].append(pdbsequence)
        chain_to_seqs[chain.id].append(pdbsequencenum)
        chain_to_seqs[chain.id].append(pdbsequencefull)
    if selchain not in chain_to_seqs.keys():
        return None
    if selchain:
        return chain_to_seqs[selchain][0], chain_to_seqs[selchain][1], chain_to_seqs[selchain][2]
    return chain_to_seqs

def min_dist_lists(atlist1, atlist2):
    distance_pairs = list(itertools.product(atlist1, atlist2))
    min_distance = min([atoms[0]-atoms[1] for atoms in distance_pairs])
    return min_distance

def get_min_dist_matrix(pdb_code, selchain, download_dir="./PDB", verbose=False):
    pdbl = PDBList(verbose=verbose)
    pdbl.retrieve_pdb_file(pdb_code,pdir=download_dir,file_format="pdb")

    make_output_dir(f"./MAPPING_JOBS/{pdb_code}:{selchain}/distmatrix")
    min_dist_matrix_f = f"./MAPPING_JOBS/{pdb_code}:{selchain}/distmatrix/{pdb_code}_{selchain}_mindist.pickle"
    if os.path.exists(min_dist_matrix_f) == False:
        pdb_file_name = "pdb" + pdb_code.lower() + ".ent"
        structure = None
        try:
            structure = Bio.PDB.PDBParser().get_structure(pdb_code, download_dir + "/" + pdb_file_name)
        except FileNotFoundError:
            pdbl.retrieve_pdb_file(pdb_code,pdir=download_dir,file_format="mmCif")
            pdb_file_name = pdb_code.lower() + ".cif"
            structure = Bio.PDB.MMCIFParser().get_structure(pdb_code, download_dir + "/" + pdb_file_name)
        modelo = structure[0]
        chain = modelo[selchain]

        residue_list = []

        residue_index = 1
        for residue in chain:
            if Bio.PDB.Polypeptide.is_aa(residue) == False:
                continue
            current_residue = residue
            number = residue.get_id()
            oneletter = ""
            try:
                oneletter = Bio.PDB.Polypeptide.three_to_one(residue.get_resname())
            except KeyError:
                oneletter = SCOPData.protein_letters_3to1[residue.get_resname()]
            try:
                num = str(number[1]) + str((number[2].rstrip( ))[0])
            except IndexError:
                num = str(number[1])
            residue_list.append((residue, residue.get_resname(), num, oneletter))
        # product
        residue_pairs = itertools.combinations(residue_list, r=2)
        dict_of_distances = {}
        for res_pair in residue_pairs:
            res_1_tuple = res_pair[0]
            res_1_obj = res_1_tuple[0]
            res_1_name = res_1_tuple[1]+res_1_tuple[2]
            res_2_tuple = res_pair[1]
            res_2_obj = res_2_tuple[0]
            res_2_name = res_2_tuple[1]+res_2_tuple[2]
            min_distance = min_dist_lists(res_1_obj.get_atoms(), res_2_obj.get_atoms())
            res_key = "-".join(sorted([res_1_name, res_2_name]))
            dict_of_distances[res_key] = min_distance
        pickle_writer(min_dist_matrix_f, dict_of_distances)
    else:
        dict_of_distances = pickle_loader(min_dist_matrix_f)
    return dict_of_distances

def calc_dssp(pdb_code, selchain, pdbsequencenum, pdbsequencefull, dssp_command, download_dir="./PDB", verbose=False):

    pdbl = PDBList(verbose=verbose)
    pdbl.retrieve_pdb_file(pdb_code,pdir=download_dir,file_format="pdb")
    
    pdb_file_name = "pdb" + pdb_code.lower() + ".ent"
    structure = None
    try:
        structure = Bio.PDB.PDBParser().get_structure(pdb_code, download_dir+"/" + pdb_file_name)
    except FileNotFoundError:
        pdbl.retrieve_pdb_file(pdb_code,pdir=download_dir,file_format="mmCif")
        pdb_file_name = pdb_code.lower() + ".cif"
        structure = Bio.PDB.MMCIFParser().get_structure(pdb_code, download_dir + "/" + pdb_file_name)
    modelo = structure[0]
    chain = modelo[selchain]

    make_output_dir(f"./MAPPING_JOBS/{pdb_code}:{selchain}/DSSP")
    dssp_file = f"./MAPPING_JOBS/{pdb_code}:{selchain}/DSSP/" + pdb_code.lower() + "_dssp.pickle"
    dssp = None
    # dssp_dict = None
    if os.path.exists(dssp_file) == False:
        dssp = DSSP(modelo, download_dir + '/'+pdb_file_name, dssp=dssp_command, acc_array="Wilke")
        # dssp_dict = dssp.property_dict
        pickle_writer(dssp_file, dssp)
    else:
        dssp = pickle_loader(dssp_file)
    dssp_info = []
    dssp_info_dict = {}
    for i in range(0, len(pdbsequencenum)):
        chain_res = pdbsequencenum[i]
        residue_key = pdbsequencefull[i]
        if (chain_res[0:3] in THREE_TO_ONE_AA_TYPES_CANONICAL.keys()):
            try:
                dssp_res = dssp[selchain, residue_key]
                dssp_info.append({"name": chain_res, "sec": str(dssp_res[2]), "phi":str(dssp_res[4]), "psi":str(dssp_res[5]), "depth":str(dssp_res[3])})
                dssp_info_dict[chain_res] = {"sec": str(dssp_res[2]), "phi":str(dssp_res[4]), "psi":str(dssp_res[5]), "depth":str(dssp_res[3])}
            except KeyError:
                if verbose:
                    print(chain_res + " not found in DSSP.")
                pass
    # print(dssp_info)
    return dssp_info_dict

def sw_align_seqs(seq1, seq2, gap_open=-12, gap_extend=-4, used_matrix=matlist.blosum100):
    """
    Align two sequences by Smith-Waterman algorithm, parameters allow controling gap opening and extension and blosum matrix used for scoring
    """
    seq1 = seq1.replace('U', 'X')
    seq2 = seq2.replace('U', 'X')
    alns = pairwise2.align.localds(seq1, seq2, used_matrix, gap_open, gap_extend)
    if len(alns) < 1:
        return [0, None]
    top_aln = alns[0]
    aln_1, aln_2, score, begin, end = top_aln
    aln_str_formatted = format_alignment(*top_aln, full_sequences=True)
    if len(aln_str_formatted) < 5:
        return [0, None]
    aln_f1, match_f, aln_f2, score_str, nothing = aln_str_formatted.split("\n")
    return [1, [aln_f1, aln_f2, begin, end, score]]

def get_index_if_matched(query_aln, match_aln, aaitype=""):
    query_ipos, query_iaa, match_ipos, match_iaa = 0,0,0,0
    aaindex_list = []
    for query_char, match_char in zip(query_aln, match_aln):
        if query_char not in ["-", "."]:
            if match_char not in ["-", "."]: # matched
                if aaitype == "query":
                    aaindex_list.append(query_iaa)
                elif aaitype == "match":
                    aaindex_list.append(match_iaa)
                else:
                    aaindex_list.append(match_ipos)
            query_iaa += 1 # match iaa represents aminoacid 0-index of each aa in query_aln
        if match_char not in ["-", "."]:
            match_iaa += 1 # match iaa represents aminoacid 0-index of each aa in matched_aln
        match_ipos += 1 # match and query ipos are equal
        query_ipos += 1 # they represent the local aln position index
    return aaindex_list

def map_comms(current_comms, backup_comms, uniprot_aaindex_list, pdb_aaindex_list, pdbsequencenum, verbose=False):
    found_comm_resids = [pos-1 for pos in current_comms if pos-1 in uniprot_aaindex_list]
    found_comm_resids_resids_indexed = [uniprot_aaindex_list.index(pos) for pos in found_comm_resids]
    found_comm_resids_resids_pdbs = [pdbsequencenum[pdb_aaindex_list[ind]] for ind in found_comm_resids_resids_indexed]
    if len(found_comm_resids_resids_pdbs) == 0 and len(backup_comms) > 0:
        if verbose:
            print("using backup comms")
        found_comm_resids = [pos-1 for pos in backup_comms if pos-1 in uniprot_aaindex_list]
        found_comm_resids_resids_indexed = [uniprot_aaindex_list.index(pos) for pos in found_comm_resids]
        found_comm_resids_resids_pdbs = [pdbsequencenum[pdb_aaindex_list[ind]] for ind in found_comm_resids_resids_indexed]
    return found_comm_resids_resids_pdbs

def map_uniprot_features(uniprot_features, uniprot_aaindex_list, pdb_aaindex_list, pdbsequencenum, pdb_mapped_features):
    for feature in uniprot_features:
        feature_resids = [pos-1 for pos in feature[1] if pos-1 in uniprot_aaindex_list]
        feature_resids_indexed = [uniprot_aaindex_list.index(pos) for pos in feature_resids]
        feature_resids_pdbs = [pdbsequencenum[pdb_aaindex_list[ind]] for ind in feature_resids_indexed]
        if feature[0] == "DISULFID":
            for feature_resids_pdb in feature_resids_pdbs:
                if feature_resids_pdb[:3] != "CYS":
                    feature_resids_pdbs = []
        if len(feature_resids_pdbs) > 0:
            pdb_mapped_features.append((feature[0],feature_resids_pdbs,feature[2]))
    return pdb_mapped_features

def retrieve_content_from_url(url_text, identifier, content_type, decoding=True, max_tries=5, verbose=False):
    make_output_dir(f"./pdbsum_data/{content_type}")
    pdbsum_data_file = f"./pdbsum_data/{content_type}/{identifier}.txt"
    content = False
    if os.path.exists(pdbsum_data_file) == False:
        base_url = "https://www.ebi.ac.uk"
        tries = 1
        content = None
        while (tries <= max_tries and not content):
            try:
                handle = urllib.request.urlopen(f"{base_url}{url_text}")
                content = handle.read().decode("iso-8859-1")
            except Exception as e:
                content = None
                time.sleep(5)
            if content:
                break
            tries += 1
        if not content:
            if verbose:
                print("404:")
                print(f"{base_url}{url_text}")
            pass
            # raise Exception("No content for file")
        # except UnicodeDecodeError:
            # content = handle.read().decode("iso-8859-1")
        with open(pdbsum_data_file, "w") as pdbsum_data_file_handle:
            pdbsum_data_file_handle.write(content)
    else:
        with open(pdbsum_data_file, "r") as pdbsum_data_file_handle:
            content = pdbsum_data_file_handle.read()
    if content == False:
        # raise Exception("No content found!!")
        if verbose:
            print("No content found!!")
        pass
    return content

def explore_residue_mappings(pdb_code, tag_example_dict, allowed_chains=[], allowed_ligands=[]):
    pdb_code = pdb_code.lower()
    pdbsum_info_url = f"/thornton-srv/databases/cgi-bin/pdbsum/GetPage.pl?pdbcode={pdb_code}&template=main.html"
    content = retrieve_content_from_url(pdbsum_info_url, pdb_code, "summary")
    
    soup = BeautifulSoup(content, 'html.parser')
    filtered_tags = []
    
    allowed_tags = ["/thornton-srv/databases/cgi-bin/pdbsum/"]
    required_tag = f".pl?pdbcode={pdb_code.lower()}"

    chains_urls = {}
    for tag in soup.find_all('a'):
        href_text = tag.get('href')
        mouseover_text = tag.get("onmouseover")
        if mouseover_text:
            if "Go to protein page for this chain" in mouseover_text:
                chain = href_text.split("chain=")[1]
                chains_urls[chain] = href_text
    for selchain, chain_url in chains_urls.items():
        # print(selchain, chain_url)
        tag_example_dict = explore_chain_information(pdb_code, selchain, chain_url, tag_example_dict)
    return tag_example_dict

def explore_chain_information(pdb_code, selchain, chain_url, tag_example_dict):

    chain_content = retrieve_content_from_url(chain_url, f"{pdb_code}_{selchain}", "summary_chain")
    detailed_soup = BeautifulSoup(chain_content, 'html.parser')
    disulphide_urls = {}
    catalytic_urls = {}
    activesite_urls = {}

    required_tag = f".pl?pdbcode={pdb_code.lower()}"
    
    for tag in detailed_soup.find_all('a'):
        href_text = tag.get('href')
        if required_tag in href_text and "o=" in href_text:
            get_text = href_text.split("o=")[1].split("&")[0]
            if get_text not in tag_example_dict.keys():
                tag_example_dict[get_text] = href_text
    return tag_example_dict

def retrieve_residue_mappings(pdb_code, allowed_chains=[], allowed_ligands=[], verbose=False):
    pdb_code = pdb_code.lower()
    pdbsum_info_url = f"/thornton-srv/databases/cgi-bin/pdbsum/GetPage.pl?pdbcode={pdb_code}&template=main.html"
    content = retrieve_content_from_url(pdbsum_info_url, pdb_code, "summary")
    
    soup = BeautifulSoup(content, 'html.parser')
    filtered_tags = []
    
    allowed_tags = ["/thornton-srv/databases/cgi-bin/pdbsum/"]
    required_tag = f".pl?pdbcode={pdb_code.lower()}"

    chains_urls = {}
    for tag in soup.find_all('a'):
        href_text = tag.get('href')
        mouseover_text = tag.get("onmouseover")
        if mouseover_text:
            if "Go to protein page for this chain" in mouseover_text:
                chain = href_text.split("chain=")[1]
                chains_urls[chain] = href_text
    disulphides_per_chain = {}
    catalytic_per_chain = {}
    activesite_per_chain = {}
    for selchain, chain_url in chains_urls.items():
        # print(selchain, chain_url)
        results = get_chain_information(pdb_code, selchain, chain_url)
        if len(results[0]) > 0:
            disulphides_per_chain[selchain] = results[0]
        if len(results[1]) > 0:
            catalytic_per_chain[selchain] = results[1]
        if len(results[2]) > 0:
            activesite_per_chain[selchain] = results[2]
    if verbose:
        print(f"Done retrieving disulphides, catalytic and active site information for each {pdb_code} chain")
    
    # LIGAND INFO PARSING
    pdbsum_ligand_url = f"/thornton-srv/databases/cgi-bin/pdbsum/GetPage.pl?pdbcode={pdb_code}&template=ligands.html&c=999"
    ligand_content = retrieve_content_from_url(pdbsum_ligand_url, pdb_code, "ligand")
    soup = BeautifulSoup(ligand_content, 'html.parser')

    ligand_urls = {}

    for tag in soup.find_all('a'):
        href_text = tag.get('href')
        if required_tag in href_text and "template=ligands.html" in href_text:
            tag_img = tag.find("img")
            if tag_img:
                continue
            lig_id = tag.get_text().strip()
            lig_idx = href_text.split("l=")[1].split('">')[0]
            if href_text in list(ligand_urls.values()):
                continue
            ligand_urls[f"{lig_id}_{lig_idx}"] = href_text
    interacting_residues_per_ligand = {}
    for lig_fid, lig_url in ligand_urls.items():
        results = get_ligand_information(pdb_code, lig_fid, lig_url)
        if len(results[1]) > 0:
            interacting_residues_per_ligand[results[0]] = results[1]
        # results = get_ligand_information(pdb_code, lig_fid, lig_url)
        # if len(results) > 0:
        #     interacting_residues_per_ligand[lig_fid] = results
    if verbose:
        print(f"Done retrieving interacting residues for each {pdb_code} ligand/metal")

    # DNA INFO PARSING
    pdbsum_dna_url = f"/thornton-srv/databases/cgi-bin/pdbsum/GetPage.pl?pdbcode={pdb_code}&template=dna.html&c=999"
    dna_content = retrieve_content_from_url(pdbsum_dna_url, pdb_code, "dna")
    
    dna_urls = []
    required_dna = f"GetPS.pl?pdbcode={pdb_code}"
    soup = BeautifulSoup(dna_content, 'html.parser')
    # http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetPage.pl?pdbcode=6kiv&template=dna.html
    for tag in soup.find_all('a'):
        href_text = tag.get('href')
        if required_dna in href_text and "psfile=nucplot.ps" in href_text:
            if "&pdf=YES" in href_text:
                continue
            if href_text not in dna_urls:
                dna_urls.append(href_text)
    dna_interacting_residues = []
    for idna, dna_url in enumerate(dna_urls):
        dna_interacting_residues.extend(get_dna_information(pdb_code,idna+1, dna_url))
    dna_interacting_residues = list(set(dna_interacting_residues))
    if verbose:
        print(f"Done retrieving interacting residues for each {pdb_code} DNA chain")

    # INTERFACE INFO PARSING
    pdbsum_interface_url =f"/thornton-srv/databases/cgi-bin/pdbsum/GetPage.pl?pdbcode={pdb_code}&template=interfaces.html&c=999"
    interface_content = retrieve_content_from_url(pdbsum_interface_url, pdb_code, "interface")
    soup = BeautifulSoup(interface_content, 'html.parser')

    interface_urls = {}
    unique_interfaces = []
    for tag in soup.find_all('a'):
        href_text = tag.get('href')
        if required_tag in href_text and "template=interfaces.html" in href_text and "RESIDUE" in href_text and href_text not in unique_interfaces:
            unique_interfaces.append(href_text)
            img_tags = tag.find_all('img')
            if len(img_tags) != 2:
                # raise Exception("More than 2 chains in pairwise interface")
                if verbose:
                    print("More than 2 chains in pairwise interface")
                continue
            from_img_tag = img_tags[0]
            from_img_tag_name = from_img_tag.get('src').split("/gif/chain")[1].split(".jpg")[0]
            to_img_tag = img_tags[1]
            to_img_tag_name = to_img_tag.get('src').split("/gif/chain")[1].split(".jpg")[0]
            interface_name = f"{from_img_tag_name}-{to_img_tag_name}"
            if href_text in interface_urls.values():
                continue
            interface_urls[interface_name] = href_text
    interacting_residues_per_interface = {}
    for interface_name,interface_url in interface_urls.items():
        results = get_interface_information(pdb_code, interface_name, interface_url)
        if len(results) > 0:
            interacting_residues_per_interface[interface_name] = results
    if verbose:
        print(f"Done retrieving interacting residues for each {pdb_code} interface")

    pdbsum_mapped_info = JSONObject()
    pdbsum_mapped_info.disulphides_per_chain = disulphides_per_chain
    pdbsum_mapped_info.catalytic_per_chain = catalytic_per_chain
    pdbsum_mapped_info.activesite_per_chain = activesite_per_chain
    pdbsum_mapped_info.interacting_residues_per_ligand = interacting_residues_per_ligand
    pdbsum_mapped_info.interacting_residues_per_interface = interacting_residues_per_interface
    pdbsum_mapped_info.dna_interacting_residues = dna_interacting_residues

    return pdbsum_mapped_info

def get_chain_information(pdb_code, selchain, chain_url):
    chain_content = retrieve_content_from_url(chain_url, f"{pdb_code}_{selchain}", "summary_chain")
    detailed_soup = BeautifulSoup(chain_content, 'html.parser')
    disulphide_urls = {}
    catalytic_urls = {}
    activesite_urls = {}

    required_tag = f".pl?pdbcode={pdb_code.lower()}"
    
    for tag in detailed_soup.find_all('a'):
        href_text = tag.get('href')
        if required_tag in href_text and "o=DISULPHIDE" in href_text:
            tag_name = tag.get_text().strip()
            disulphide_urls[tag_name] = href_text
        elif required_tag in href_text and "o=CSA" in href_text:
            tag_name = tag.get_text().strip()
            catalytic_urls[tag_name] = href_text
        elif required_tag in href_text and "o=SITE" in href_text:
            tag_name = tag.get_text().strip()
            activesite_urls[tag_name] = href_text
    
    disulphide_list = []
    for k,disulphide_url in disulphide_urls.items():
        disulphide_list.extend(get_disulphide(pdb_code, selchain, disulphide_url))
    disulphide_list = list(set(disulphide_list))

    catalytic_list = []
    for k, catalytic_url in catalytic_urls.items():
        catalytic_list.extend(get_catalytic(pdb_code, selchain, catalytic_url))
    catalytic_list = list(set(catalytic_list))

    activesite_list = []
    for k,activesite_url in activesite_urls.items():
        activesite_list.extend(get_activesite(pdb_code, selchain, activesite_url))
    activesite_list = list(set(activesite_list))

    return disulphide_list, catalytic_list, activesite_list

def get_disulphide(pdb_code, selchain, disulphide_url):
    disulphide_chain_content = retrieve_content_from_url(disulphide_url, f"{pdb_code}_{selchain}_disulphide", "summary_chain/disulphide")
    detailed_soup2 = BeautifulSoup(disulphide_chain_content, 'html.parser')

    disulphide_pairs = []
    for tag in detailed_soup2.find_all('table'):
        table_content = tag.decode_contents()
        if "cysteine" in table_content and "table" not in table_content:
            rows = tag.find_all('tr')
            valid_rows = 0
            for row in rows:
                row_content = row.decode_contents()
                cols = row.find_all("td")
                if len(cols) == 15 and "cysteine" not in row_content:
                    chain_from, cys_from_num = cols[0].find("font").get_text().split()
                    chain_to, cys_to_num = cols[2].find("font").get_text().split()
                    if chain_from == selchain and chain_to == selchain:
                        disulphide_pairs.append((cys_from_num,cys_to_num))
                    valid_rows += 1
    return disulphide_pairs

def get_catalytic(pdb_code, selchain, catalytic_url):
    catalytic_content = retrieve_content_from_url(catalytic_url, f"{pdb_code}_{selchain}_catalytic", "summary_chain/catalytic")
    detailed_soup2 = BeautifulSoup(catalytic_content, 'html.parser')

    catalytic_residues = []
    for tag in detailed_soup2.find_all('b'):
        b_content = tag.decode_contents()
        if "Catalytic residues" in b_content:
            resid_tag = tag.find("font")
            if resid_tag:
                residues = resid_tag.get_text().split("-")
                residues = [res.strip() for res in residues]
                chains = [res.split('(')[1].split(')')[0] for res in residues]
                residues = [res.split('(')[0] for res in residues]
                residues = [res for i, res in enumerate(residues) if chains[i] == selchain]
                residues = [(res[:3].upper(), res[3:]) for res in residues]
                # print("residues")
                # print(residues)
                # print("chains")
                # print(chains)
                catalytic_residues.extend(residues)
    # quit()
    return catalytic_residues

def get_activesite(pdb_code, selchain, activesite_url):
    activesite_content = retrieve_content_from_url(activesite_url, f"{pdb_code}_{selchain}_activesite", "summary_chain/activesite")
    detailed_soup2 = BeautifulSoup(activesite_content, 'html.parser')

    activesite_residues = []
    for tag in detailed_soup2.find_all('td'):
        tag_content = tag.decode_contents()
        if "Active site residue(s):-" in tag_content and "View site in RasMol" not in tag_content:
            # print("tag")
            # print(tag)
            pre = tag.find("pre")
            data_list = pre.get_text().strip().split()
            for i in range(0,len(data_list),3):
                restype = data_list[i]
                chain = data_list[i+1]
                resnumber = data_list[i+2]
                if chain == selchain:
                    activesite_residues.append((restype.upper(), resnumber))
    # print("activesite_residues")
    # print(activesite_residues)
    # quit()
    return activesite_residues

def get_ligand_information(pdb_code, lig_fid, lig_url, verbose=False):
    lig_content = retrieve_content_from_url(lig_url, f"{pdb_code}_{lig_fid}_lig", "ligand/lig")
    detailed_soup = BeautifulSoup(lig_content, 'html.parser')

    ligand_interacting_residues = []

    txt_file_url = False
    required_lig = f"GetLigInt.pl?pdb={pdb_code}"
    for tag in detailed_soup.find_all('a'):
        href_text = tag.get('href')
        if required_lig in href_text:
            txt_file_url = href_text
            break
    if txt_file_url == False:
        if verbose:
            print("lig_url")
            print(lig_url)
            print("txt_file_url")
            print(txt_file_url)
            print("Url 404 ligand")
        # raise Exception("Url 404 ligand")
    lig_content_txt = retrieve_content_from_url(txt_file_url, f"{pdb_code}_{lig_fid}_lig_txt", "ligand/lig")
    ligand_interacting_residues = contact_parse(lig_content_txt, name=True)
    return ligand_interacting_residues

def get_interface_information(pdb_code, interface_name, interface_url, verbose=False):
    interface_content = retrieve_content_from_url(interface_url, f"{pdb_code}_{interface_name}_interface", "interface/each")
    detailed_soup = BeautifulSoup(interface_content, 'html.parser')

    interface_interacting_residues = []

    txt_file_url = False
    required_interface = f"GetIface.pl?pdb={pdb_code}"
    for tag in detailed_soup.find_all('a'):
        href_text = tag.get('href')
        if required_interface in href_text:
            txt_file_url = href_text
            break
    if txt_file_url == False:
        # raise Exception("Url 404 interface")
        if verbose:
            print("Url 404 interface")
    interface_content_txt = retrieve_content_from_url(txt_file_url, f"{pdb_code}_{interface_name}_interface_txt", "interface/each")
    interface_interacting_residues = contact_parse(interface_content_txt)
    # print("interface_interacting_residues")
    # print(interface_interacting_residues)
    return interface_interacting_residues

def contact_parse(contact_text, name=False):
    partner = False
    interaction_list = []
    interaction_classes = ["Hydrogen bonds", "Non-bonded contacts", "Disulphide bonds", "Salt bridges"]
    interaction_types = ["hydrogen","non_bonded","disulphide","salt_bridge"]
    for index_class, interaction_class in enumerate(interaction_classes):
        interaction_type = interaction_types[index_class]
        if interaction_class not in contact_text:
            continue
        interaction_block = contact_text.split(f"{interaction_class}")[1]
        class_number = 1
        interaction_block_lines = interaction_block.split("\n")
        parse_start = "None"
        for line in interaction_block_lines:
            line = line.strip()
            if not line:
                continue
            line_list = line.split()
            if line_list[0] == f"{class_number}.":
                parse_start = "Running"

                restype = line_list[3]
                resnum = line_list[4]
                chain = line_list[5]

                restype2 = line_list[9]
                resnum2 = line_list[10]
                chain2 = line_list[11]
                partner = restype2+"__"+resnum2+"("+chain2+")"
                if restype in THREE_TO_ONE_AA_TYPES_CANONICAL.keys():
                    interaction_list.append((restype, resnum, chain, interaction_type))
                if restype2 in THREE_TO_ONE_AA_TYPES_CANONICAL.keys():
                    interaction_list.append((restype2, resnum2, chain2, interaction_type))
                class_number += 1
            elif parse_start == "Running":
                # parse_start = "Stopped"
                break
    interaction_list = list(set(interaction_list))
    if not name:
        return interaction_list
    # print(partner, interaction_list)
    return [partner, interaction_list]

def get_dna_information(pdb_code, dna_header, dna_url, verbose=False):
    ps_content = retrieve_content_from_url(dna_url, f"{pdb_code}_{dna_header}", "dna/each", decoding=True)
    lines_of_pscontent = ps_content.split('\n')
    lines_with_residues = [line for line in lines_of_pscontent if any(aa in line for aa in THREE_AA_TYPES_CANONICAL_SMALL)]
    lines_with_residues = [line for line in lines_with_residues if '%%' not in line]
    current_align = False
    dna_interacting_residues = []
    for line in lines_with_residues:
        resid = None
        if "Lalign" in line:
            current_align = "Lalign"
            resid = line.split(" ) ")[0].split()
            if len(resid) > 1:
                resid = resid[1]
            else:
                resid = resid[0][1:]
            resname = resid[:3].upper()
            chain = resid.split("(")[1][0]
            resnum = resid.split("(")[0][3:]
            dna_interacting_residues.append((resname, resnum, chain))
        elif "Ralign" in line:
            current_align = "Ralign"
            resid = line.split(") ")[0].split()[0][1:]
            resname = resid[:3].upper()
            chain = resid.split("(")[1][0]
            resnum = resid.split("(")[0][3:]
            dna_interacting_residues.append((resname, resnum, chain))
        elif "Print" in line:
            if not current_align:
                # raise Exception("No align in PS file")
                if verbose:
                    print("No align in PS file")
            elif current_align == "Lalign":
                resid = line.split(" ) ")[0].split()
                if len(resid) > 1:
                    resid = resid[1]
                else:
                    resid = resid[0][1:]
                resname = resid[:3].upper()
                chain = resid.split("(")[1][0]
                resnum = resid.split("(")[0][3:]
                dna_interacting_residues.append((resname, resnum, chain))
            elif current_align == "Ralign":
                resid = line.split(") ")[0].split()[0][1:]
                resname = resid[:3].upper()
                chain = resid.split("(")[1][0]
                resnum = resid.split("(")[0][3:]
                dna_interacting_residues.append((resname, resnum, chain))
        else:
            # raise Exception(f"in PS file: Cant interpret line:\n{line}")
            if verbose:
                print("in PS file: Cant interpret line:\n{line}")
            pass
    dna_interacting_residues = list(set(dna_interacting_residues))
    return dna_interacting_residues

def merge_information_act_site(selchain, residue_information, pdb_mapped_features):
    act_features = [feature for feature in pdb_mapped_features if feature[0] == "ACT_SITE"]
    for feature in act_features:
        if selchain not in residue_information.activesite_per_chain.keys():
            residue_information.activesite_per_chain[selchain] = []
        residue_information.activesite_per_chain[selchain].extend([(respos_f[:3], respos_f[3:]) for respos_f in feature[1]])
        residue_information.activesite_per_chain[selchain] = list(set(residue_information.activesite_per_chain[selchain]))
    return residue_information

def merge_information_ligands(selchain, residue_information, pdb_mapped_features, bad_ligand_list):
    ligand_features = [feature for feature in pdb_mapped_features if feature[0] == "BINDING" or feature[0] == "NP_BIND"]
    for lf_i, feature in enumerate(ligand_features):
        # lig_id = "uniprot_lig" + str(lf_i)
        lig_id = "UnP" + str(lf_i)
        lig_name = None
        # residue_information.interacting_residues_per_ligand
        if len(feature[2]) > 0:
            for qualifier in feature[2]:
                if 'note' in qualifier[0]:
                    lig_name = re.search('[A-Z]{2,3}', qualifier[1])
                    if lig_name:
                        lig_name = lig_name.group()
                        lig_id = lig_name + "__" + str(lf_i)
                        break
        if lig_name not in bad_ligand_list:
            residue_information.interacting_residues_per_ligand[lig_id] = [(respos_f[:3], respos_f[3:], selchain) for respos_f in feature[1]]
    return residue_information

def merge_information_additional(selchain, residue_information, pdb_mapped_features_additional):
    for feature in pdb_mapped_features_additional:
        if selchain not in residue_information.additional_info.keys():
            residue_information.additional_info[selchain] = []
        residue_information.additional_info[selchain].extend([(respos_f[:3], respos_f[3:]) for respos_f in feature[1]])
        residue_information.additional_info[selchain] = list(set(residue_information.additional_info[selchain]))
    return residue_information

def merge_information_disulphides(selchain, residue_information, pdb_mapped_features, disulphide_list):
    residue_information.intrachain_disulphide_list = disulphide_list
    if selchain in residue_information.disulphides_per_chain.keys():
        residue_information.intrachain_disulphide_list.extend(residue_information.disulphides_per_chain[selchain])
    dis_features = [feature for feature in pdb_mapped_features if feature[0] == "DISULFID"]
    for feature in dis_features:
        # dis_tuple = tuple([respos_f[3:] for respos_f in feature[1]])
        # dis_tuple = [feature[1][0], feature[1][1]]
        dis_tuple = feature[1].copy()
        if len(dis_tuple) == 1:
            dis_tuple = tuple([dis_tuple[0][:3], dis_tuple[0][3:]])
        elif len(dis_tuple) == 2:
            dis_tuple = tuple(sorted([tuple([dis_tuple[0][:3], dis_tuple[0][3:]]),tuple([dis_tuple[1][:3], dis_tuple[1][3:]])]))
        # else:
        #     print("error")
        #     quit()
        residue_information.intrachain_disulphide_list.append(dis_tuple)
    residue_information.intrachain_disulphide_list = list(set(residue_information.intrachain_disulphide_list))
    return residue_information

def filter_for_chain(selchain, residue_information, bad_ligand_list):
    residue_information.ilv_for_chain = [v[0]+v[1] for v in residue_information.ilvclusters]
    residue_information.ilv_contact_for_chain = [v[0]+v[1] for v in residue_information.ilvexpanded]
    
    active_site_for_chain = []
    if selchain in residue_information.catalytic_per_chain.keys():
        active_site_for_chain.extend([v[0]+v[1] for v in residue_information.catalytic_per_chain[selchain]])
    if selchain in residue_information.activesite_per_chain.keys():
        active_site_for_chain.extend([v[0]+v[1] for v in residue_information.activesite_per_chain[selchain]])
    residue_information.active_site_for_chain = active_site_for_chain

    ligand_keys = residue_information.interacting_residues_per_ligand.keys()

    uniprot_ligand_keys = [k for k in ligand_keys if "UnP" in k]
    filter_ligand_keys = [k for k in ligand_keys if k not in uniprot_ligand_keys]

    all_keys_filtered = [k for k in filter_ligand_keys if k.split("__")[0] not in bad_ligand_list]
    all_keys_filtered.extend(uniprot_ligand_keys)
    
    all_ligand_binding_for_chain = []
    for k in all_keys_filtered:
        values_filtered = [v[0]+v[1] for v in residue_information.interacting_residues_per_ligand[k] if v[2] == selchain]
        all_ligand_binding_for_chain.extend(values_filtered)
    residue_information.all_ligand_binding_for_chain = all_ligand_binding_for_chain

    dna_binding_for_chain = []
    dna_binding_for_chain.extend([v[0]+v[1] for v in residue_information.dna_interacting_residues if v[2] == selchain])
    residue_information.dna_binding_for_chain = dna_binding_for_chain

    interfacing_for_chain = []
    interfacing_for_chain_all = []
    for k,v in residue_information.interacting_residues_per_interface.items():
        if selchain in k:
            interfacing_for_chain.extend([x[0]+x[1] for x in v if x[3] != "non_bonded"])
            interfacing_for_chain_all.extend([x[0]+x[1] for x in v])
    residue_information.interfacing_for_chain = interfacing_for_chain
    residue_information.interfacing_for_chain_all = interfacing_for_chain_all
    return residue_information

def get_interactions_charged_residues_intrachain(pdb_code, selchain, download_dir="./PDB", salt_bridge_cutoff=6.0,verbose=False):
    # conditions for charged: Pos and Neg resids
    # This includes:
    # Amide of First residue, Arg and Lys Side-Chain Nitrogens
    # Carboxyl of Last residue, Glu and Asp Side-Chain Oxygens
    positive_first = ["N"]
    negative_last = ["O", "OXT"]
    negative_dict = {
        "ASP": ["OD1", "OD2"],
        "GLU": ["OE1", "OE2"]
    }
    positive_dict = {
        "ARG":["NE","NH1","NH2"],
        "LYS":["NZ"]
    }
    pdbl = PDBList(verbose=verbose)
    pdbl.retrieve_pdb_file(pdb_code,pdir=download_dir,file_format="pdb")
    
    pdb_file_name = "pdb" + pdb_code.lower() + ".ent"
    structure = None
    try:
        structure = Bio.PDB.PDBParser().get_structure(pdb_code, download_dir + "/" + pdb_file_name)
    except FileNotFoundError:
        pdbl.retrieve_pdb_file(pdb_code,pdir=download_dir,file_format="mmCif")
        pdb_file_name = pdb_code.lower() + ".cif"
        structure = Bio.PDB.MMCIFParser().get_structure(pdb_code, download_dir + "/" + pdb_file_name)
    modelo = structure[0]
    chain = modelo[selchain]

    positives_atoms_list = []
    negatives_atoms_list = []

    positives_residues_list = []
    negatives_residues_list = []

    residue_index = 1
    current_residue = None
    current_resname = None
    current_resnumb = None
    for residue in chain:
        if Bio.PDB.Polypeptide.is_aa(residue) == False:
            continue
        current_residue = residue
        number = residue.get_id()
        oneletter = ""
        try:
            oneletter = Bio.PDB.Polypeptide.three_to_one(residue.get_resname())
        except KeyError:
            oneletter = SCOPData.protein_letters_3to1[residue.get_resname()]
        try:
            num = str(number[1]) + str((number[2].rstrip( ))[0])
        except IndexError:
            num = str(number[1])
        current_resname = residue.get_resname()
        current_resnumb = num
        if (residue_index == 1):
            positive_ats = [at for at in residue if at.get_name() in positive_first]
            positives_atoms_list.extend(positive_ats)
            positives_residues_list.extend([(current_resname,current_resnumb) for at in positive_ats])
        if residue.get_resname() in positive_dict.keys():
            positive_list = positive_dict[residue.get_resname()]
            positive_ats = [at for at in residue if at.get_name() in positive_list]
            positives_atoms_list.extend(positive_ats)
            positives_residues_list.extend([(current_resname,current_resnumb) for at in positive_ats])
        if residue.get_resname() in negative_dict.keys():
            negative_list = negative_dict[residue.get_resname()]
            negative_ats = [at for at in residue if at.get_name() in negative_list]
            negatives_atoms_list.extend(negative_ats)
            negatives_residues_list.extend([(current_resname,current_resnumb) for at in negative_ats])
        residue_index += 1
    negative_ats = [at for at in current_residue if at.get_name() in negative_last]
    negatives_atoms_list.extend(negative_ats)
    negatives_residues_list.extend([(current_resname,current_resnumb) for at in negative_ats])

    salt_bridge_pairs = []
    for inegpos, negat in enumerate(negatives_atoms_list):
        negative_residue = negatives_residues_list[inegpos]
        for ipospos, posat in enumerate(positives_atoms_list):
            positive_residue = positives_residues_list[ipospos]
            if abs(negat - posat) <= salt_bridge_cutoff:
                salt_bridge_pairs.append((negative_residue,positive_residue))
    salt_bridge_pairs = list(set(salt_bridge_pairs))
    return salt_bridge_pairs

def get_interactions_disulphide_intrachain(pdb_code, selchain, download_dir="./PDB", disulphide_cutoff=2.2,verbose=False):
    # conditions for disulphide: Cys residues
    # Sulphurs within 2.2 Angstroms

    pdbl = PDBList(verbose=verbose)
    pdbl.retrieve_pdb_file(pdb_code,pdir=download_dir,file_format="pdb")
    
    pdb_file_name = "pdb" + pdb_code.lower() + ".ent"
    structure = None
    try:
        structure = Bio.PDB.PDBParser().get_structure(pdb_code, download_dir+"/" + pdb_file_name)
    except FileNotFoundError:
        pdbl.retrieve_pdb_file(pdb_code,pdir=download_dir,file_format="mmCif")
        pdb_file_name = pdb_code.lower() + ".cif"
        structure = Bio.PDB.MMCIFParser().get_structure(pdb_code, download_dir+"/" + pdb_file_name)
    modelo = structure[0]
    chain = modelo[selchain]

    cysteine_atom_list = []
    cysteine_residue_list = []

    residue_index = 1
    current_residue = None
    for residue in chain:
        if Bio.PDB.Polypeptide.is_aa(residue) == False:
            continue
        current_residue = residue
        number = residue.get_id()
        oneletter = ""
        try:
            oneletter = Bio.PDB.Polypeptide.three_to_one(residue.get_resname())
        except KeyError:
            oneletter = SCOPData.protein_letters_3to1[residue.get_resname()]
        try:
            num = str(number[1]) + str((number[2].rstrip( ))[0])
        except IndexError:
            num = str(number[1])
        current_resname = residue.get_resname()
        current_resnumb = num
        if residue.get_resname() == "CYS":
            sulphur_list = [at for at in residue if at.get_name() == "SG"]
            cysteine_atom_list.extend(sulphur_list)
            cysteine_residue_list.extend([(current_resname,current_resnumb) for at in sulphur_list])

    disulphide_pairs = []
    for icys1 in range(0, len(cysteine_atom_list)-1):
        cys_at_1 = cysteine_atom_list[icys1]
        cys_res_1 = cysteine_residue_list[icys1]
        for icys2 in range(icys1+1, len(cysteine_atom_list)):
            cys_at_2 = cysteine_atom_list[icys2]
            cys_res_2 = cysteine_residue_list[icys2]
            if (abs(cys_at_1 - cys_at_2) <= disulphide_cutoff):
                disulphide_pairs.append((cys_res_1,cys_res_2))
    return disulphide_pairs

def get_ilv_clusters(pdb_code, selchain, download_dir="./PDB", neighbour_max=6.56, verbose=False, phe=False):

    pdbl = PDBList(verbose=verbose)
    pdbl.retrieve_pdb_file(pdb_code,pdir=download_dir,file_format="pdb")
    
    pdb_file_name = "pdb" + pdb_code.lower() + ".ent"
    structure = None
    try:
        structure = Bio.PDB.PDBParser().get_structure(pdb_code, download_dir+"/" + pdb_file_name)
    except FileNotFoundError:
        pdbl.retrieve_pdb_file(pdb_code,pdir=download_dir,file_format="mmCif")
        pdb_file_name = pdb_code.lower() + ".cif"
        structure = Bio.PDB.MMCIFParser().get_structure(pdb_code, download_dir+"/" + pdb_file_name)
    modelo = structure[0]
    chain = modelo[selchain]

    residue_index = 1
    current_residue = None

    # main_chain_at_names = []
    side_chains  = {
        "ILE":["CB","CD1","CG1","CG2"],
        "LEU":["CB","CD1","CD2","CG"],
        "VAL":["CB","CG1","CG2"],
    }
    if phe ==  True:
        side_chains['PHE'] = ["CB","CD1","CD2","CE1","CE2","CG","CZ"]

    ilv_atom_residue_list = []
    ilv_at_dict = {}
    for residue in chain:
        if Bio.PDB.Polypeptide.is_aa(residue) == False:
            continue
        current_residue = residue
        number = residue.get_id()
        oneletter = ""
        try:
            oneletter = Bio.PDB.Polypeptide.three_to_one(residue.get_resname())
        except KeyError:
            oneletter = SCOPData.protein_letters_3to1[residue.get_resname()]
        try:
            num = str(number[1]) + str((number[2].rstrip( ))[0])
        except IndexError:
            num = str(number[1])
        current_resname = residue.get_resname()
        current_resnumb = num
        if residue.get_resname() in side_chains.keys():
            filter_ilv = side_chains[current_resname]
            ilv_atom_list = [at for at in residue if at.get_name() in filter_ilv]
            ilv_atom_residue_list.append(((current_resname,current_resnumb), ilv_atom_list))
            ilv_at_dict[(current_resname,current_resnumb)] = ilv_atom_list

    # get neighbours within neighbour_max
    ilv_residues_distance_pairs = list(itertools.product(ilv_atom_residue_list, ilv_atom_residue_list))
    ilv_residues_distance_pairs_filter_same = [pair for pair in ilv_residues_distance_pairs if pair[0][0] != pair[1][0]]
    ilv_residues_distance_pairs_filter_dist = [pair for pair in ilv_residues_distance_pairs_filter_same if min_dist_lists(pair[0][1],pair[1][1]) <= neighbour_max]

    # ilv_residues_networks = 
    # ilv_edges = [(pair[0][0][0]+pair[0][0][1],pair[1][0][0]+pair[1][0][1]) for pair in ilv_residues_distance_pairs_filter_dist]
    ilv_edges = [(pair[0][0],pair[1][0]) for pair in ilv_residues_distance_pairs_filter_dist]
    ilv_resids = [pair[0][0] for pair in ilv_residues_distance_pairs_filter_dist]
    ilv_resids.extend([pair[1][0] for pair in ilv_residues_distance_pairs_filter_dist])
    ilv_resids = list(set(ilv_resids))

    G = nx.Graph()
    G.add_edges_from(ilv_edges)
    remove = [node for node,degree in dict(G.degree()).items() if degree < 2]
    G.remove_nodes_from(remove)
    return ilv_resids

def hydrophobic_contacting_ilv_cluster(pdb_code, selchain, ilv_resids, download_dir="./PDB", hydrophobic_cutoff=5.0,verbose=False):
    hydrophobic_strict = {
        "ALA": ["CB"], "CYS":["CB"],
        "ILE":["CB","CD1","CG1","CG2"],
        "LEU":["CB","CD1","CD2","CG"],
        "MET":["CB","CE","CG","SD"],
        "PHE":["CB","CD1","CD2","CE1","CE2","CG","CZ"],
        "TRP":["CB","CD1","CD2","CE2","CE3","CG","CH2","CZ2","CZ3"],
        "VAL":["CB","CG1","CG2"],
    }

    pdbl = PDBList(verbose=verbose)
    pdbl.retrieve_pdb_file(pdb_code,pdir=download_dir,file_format="pdb")
    
    pdb_file_name = "pdb" + pdb_code.lower() + ".ent"
    structure = None
    try:
        structure = Bio.PDB.PDBParser().get_structure(pdb_code, download_dir+"/" + pdb_file_name)
    except FileNotFoundError:
        pdbl.retrieve_pdb_file(pdb_code,pdir=download_dir,file_format="mmCif")
        pdb_file_name = pdb_code.lower() + ".cif"
        structure = Bio.PDB.MMCIFParser().get_structure(pdb_code, download_dir+"/" + pdb_file_name)
    modelo = structure[0]
    chain = modelo[selchain]

    hydrophobic_residue_list = []
    for residue in chain:
        if Bio.PDB.Polypeptide.is_aa(residue) == False:
            continue
        current_residue = residue
        number = residue.get_id()
        oneletter = ""
        try:
            oneletter = Bio.PDB.Polypeptide.three_to_one(residue.get_resname())
        except KeyError:
            oneletter = SCOPData.protein_letters_3to1[residue.get_resname()]
        try:
            num = str(number[1]) + str((number[2].rstrip( ))[0])
        except IndexError:
            num = str(number[1])
        current_resname = residue.get_resname()
        current_resnumb = num
        if residue.get_resname() in hydrophobic_strict.keys():
            hydrophobic_list_filter = hydrophobic_strict[current_resname]
            hydrophobic_atom_list = [at for at in residue if at.get_name() in hydrophobic_list_filter]
            hydrophobic_residue_list.append(((current_resname,current_resnumb),hydrophobic_atom_list))
    hydrophobic_residue_list_filtered = [restuple for restuple in hydrophobic_residue_list if restuple[0] not in ilv_resids]
    ilv_residue_list_filtered = [restuple for restuple in hydrophobic_residue_list if restuple[0] in ilv_resids]

    expanded_ilv_residues_distance_pairs = list(itertools.product(ilv_residue_list_filtered, hydrophobic_residue_list_filtered))
    expanded_ilv_residues_distance_pairs_filter_same = [pair for pair in expanded_ilv_residues_distance_pairs if pair[0][0] != pair[1][0]]
    expanded_ilv_residues_distance_pairs_filter_dist = [pair for pair in expanded_ilv_residues_distance_pairs_filter_same if min_dist_lists(pair[0][1],pair[1][1]) <= hydrophobic_cutoff]   
    
    expanded_ilv_residues = [pair[0][0] for pair in expanded_ilv_residues_distance_pairs_filter_dist]
    expanded_ilv_residues.extend([pair[1][0] for pair in expanded_ilv_residues_distance_pairs_filter_dist])
    expanded_ilv_residues = [resp for resp in expanded_ilv_residues if resp not in ilv_resids]
    expanded_ilv_residues = list(set(expanded_ilv_residues))
    return expanded_ilv_residues

def calcDynamicFlexibilityIndex(prs_matrix, atoms, select):
    profiles = sliceAtomicData(prs_matrix, atoms, select, axis=0)
    return np.sum(profiles, axis=1)/np.sum(prs_matrix)

def calc_prody_anm_get_features(pdb_code, selchain):
    make_output_dir(f"./MAPPING_JOBS/{pdb_code}:{selchain}/anm")
    anm_file = f"./MAPPING_JOBS/{pdb_code}:{selchain}/anm/{pdb_code}_{selchain}_anm.pickle"
    prody_model_ca = parsePDB(f"{pdb_code.lower()}{selchain}", subset='ca')
    if os.path.exists(anm_file) == False:
        anm = ANM()
        anm.buildHessian(prody_model_ca, cutoff=13.0)
        anm.calcModes('all')
        pickle_writer(anm_file, anm)
    else:
        anm = pickle_loader(anm_file)
    names = list(prody_model_ca.getResnames())
    names = [str(name) for name in names]
    nums = list(prody_model_ca.getResnums())
    nums = [str(num) for num in nums]
    resnums = [name+num for name, num in zip(names, nums)]
    resnum_characteristics = dict((aapos, {}) for aapos in resnums)
    prs_mat, effectiveness, sensitivity = calcPerturbResponse(anm)
    stiffness = calcMechStiff(anm, prody_model_ca)
    meanStiff = np.array([np.mean(stiffness, axis=0)])
    meanStiff = list(meanStiff[0])
    dfis = calcDynamicFlexibilityIndex(prs_mat, prody_model_ca, "all")
    sq_fluctuations = calcSqFlucts(anm)
    ires = 0
    for effect, sensiv, stiff, dfi, sq_flutc in zip(effectiveness, sensitivity, meanStiff, dfis, sq_fluctuations):
        aapos = resnums[ires]
        resnum_characteristics[aapos]["effect"] = effect
        resnum_characteristics[aapos]["sensv"] = sensiv
        resnum_characteristics[aapos]["stiff"] = stiff
        resnum_characteristics[aapos]["dfi"] = dfi
        resnum_characteristics[aapos]["fluc"] = sq_flutc
        ires += 1
    # print(resnum_characteristics)
    # quit()
    os.remove(f"{pdb_code.lower()}.pdb.gz")
    return resnum_characteristics


def main_run_get_annotations(arguments, verbose=False):
    DSSP_COMMAND = 'EDIT ME'
    #uniprot_code, pdb_code, selchain, map_file, comms_file
    if len(arguments) != 6:
        return([0, "Missing arguments", None])

    #current_comms, backup_comms,
    uniprot_code = arguments[1]
    pdb_code = arguments[2]
    selchain = arguments[3]
    if os.path.exists(f"{uniprot_code}_{pdb_code}:{selchain}_annomap.pickle") == False:
        map_file = arguments[4]
        comm_file = arguments[5]
        map_loaded = pickle_loader(map_file)
        pdb_chain_list = [pdb_c[0] for pdb_c in map_loaded[1]]
        if f"{pdb_code}:{selchain}" not in pdb_chain_list:
            return([0, "Unmapped PDB", None])
        comm_loaded = pickle_loader(comm_file)
        current_comms = comm_loaded[1]
        backup_comms = comm_loaded[2]
        if verbose:
            print("uniprot_code")
            print(uniprot_code)
            print("pdb_code")
            print(pdb_code)
            print("selchain")
            print(selchain)
            print("current_comms")
            print(current_comms)
            print("backup_comms")
            print(backup_comms)
        # make_output_dir(f"./features/{uniprot_code}/")
        make_output_dir(f"./MAPPING_JOBS/{pdb_code}:{selchain}/UNIPROT_FEATURES/{uniprot_code}")
        make_output_dir(f"./MAPPING_JOBS/{pdb_code}:{selchain}/MAPPED_INFO")
        uniprot_features, uniprot_features_additional = retrieve_uniprot_features(uniprot_code, f"./MAPPING_JOBS/{pdb_code}:{selchain}/UNIPROT_FEATURES/{uniprot_code}/{uniprot_code}.txt")
        if verbose:
            print("uniprot_features")
            print(uniprot_features)
            print("uniprot_features_additional")
            print(uniprot_features_additional)
        uniprot_sequence = get_sequence_from_uniprot_code(uniprot_code)

        pdbsequence, pdbsequencenum, pdbsequencefull = retrieve_pdb_seq(pdb_code, selchain)
        converted_pdbsequence = [ONE_TO_THREE1[aa1].upper()+pdbsequencenum[i][3:] for i, aa1 in enumerate(pdbsequence)]
        if verbose:
            print("pdbsequencefull")
            print(pdbsequencefull)
            print("converted_pdbsequence")
            print(converted_pdbsequence)
        distance_matrix = get_min_dist_matrix(pdb_code, selchain)
        dssp_info_dict = calc_dssp(pdb_code, selchain, pdbsequencenum,pdbsequencefull, DSSP_COMMAND)
        if verbose:
            print("dssp_info_dict")
            print(dssp_info_dict)

        aln_results = sw_align_seqs(uniprot_sequence, pdbsequence)
        uniprot_aln, pdbchain_aln, begin, end, score = aln_results[1]
        
        uniprot_aaindex_list = get_index_if_matched(uniprot_aln, pdbchain_aln, "query")
        pdb_aaindex_list = get_index_if_matched(uniprot_aln, pdbchain_aln, "match")
        pos_index_list = get_index_if_matched(uniprot_aln, pdbchain_aln, "")

        found_comm_resids_in_pdb = map_comms(current_comms, backup_comms, uniprot_aaindex_list, pdb_aaindex_list, pdbsequencenum, verbose=False)
        if verbose:
            print("found_comm_resids_in_pdb")
            print(found_comm_resids_in_pdb)
        pdb_mapped_features = []
        pdb_mapped_features = map_uniprot_features(uniprot_features, uniprot_aaindex_list, pdb_aaindex_list, pdbsequencenum, pdb_mapped_features)
        pdb_mapped_features_additional = map_uniprot_features(uniprot_features_additional, uniprot_aaindex_list, pdb_aaindex_list, pdbsequencenum, pdb_mapped_features)
        if verbose:
            print("pdb_mapped_features")
            print(pdb_mapped_features)
            print("pdb_mapped_features_additional")
            print(pdb_mapped_features_additional)
        residue_information = retrieve_residue_mappings(pdb_code)

        residue_information.pdbsequencenum = pdbsequencenum
        residue_information.pdbsequencefull = pdbsequencefull
        residue_information.converted_pdbsequence = converted_pdbsequence
        residue_information.distance_matrix = distance_matrix
        residue_information.dssp_info = dssp_info_dict
        residue_information.mapped_communities = found_comm_resids_in_pdb
        residue_information.additional_info = uniprot_features_additional

        disulphide_list = get_interactions_disulphide_intrachain(pdb_code, selchain)
        residue_information = merge_information_disulphides(selchain, residue_information, pdb_mapped_features, disulphide_list)
        residue_information = merge_information_act_site(selchain, residue_information, pdb_mapped_features)
        bad_ligand_list = get_bad_ligand_list()
        residue_information = merge_information_ligands(selchain, residue_information, pdb_mapped_features, bad_ligand_list)

        # residue_information = merge_information_additional(selchain, residue_information, pdb_mapped_features_additional)

        residue_information.salt_bridges = get_interactions_charged_residues_intrachain(pdb_code, selchain)
        residue_information.ilvclusters = get_ilv_clusters(pdb_code, selchain)
        residue_information.ilvexpanded = hydrophobic_contacting_ilv_cluster(pdb_code, selchain, residue_information.ilvclusters)

        residue_information = filter_for_chain(selchain, residue_information, bad_ligand_list)
        anm_characteristics = calc_prody_anm_get_features(pdb_code, selchain)
        residue_information.anm_characteristics = anm_characteristics
        if verbose:
            for prop, value in vars(residue_information).items():
                print(prop)
                # if prop not in ["distance_matrix", "dssp_info"]:
                # if prop in ["ilv_for_chain", "active_site_for_chain", "all_ligand_binding_for_chain", "dna_binding_for_chain", "interfacing_for_chain", "interfacing_for_chain_all"]:
                # if prop in ["anm_characteristics"]:
                    # print(prop, ":", value)
        results = [residue_information]
        return [1, "Annotations mapped", results]
    else:
        return [1, "Annotations mapped", pickle_loader(f"{uniprot_code}_{pdb_code}:{selchain}_annomap.pickle")]

if __name__ == "__main__":
    arguments = sys.argv
    status, status_string, results = main_run_get_annotations(arguments)
    if results[0] == 0:
        raise Exception(results[1])
    pickle_writer(f"{arguments[1]}_{arguments[2]}:{arguments[3]}_annomap.pickle", results)
    # disulphides_per_chain
    # catalytic_per_chain
    # activesite_per_chain
    # interacting_residues_per_ligand
    # interacting_residues_per_interface
    # dna_interacting_residues
    # pdbsequencenum
    # converted_pdbsequence
    # distance_matrix
    # dssp_info
    # intrachain_disulphide_list
    # salt_bridges
    # ilvclusters
    # ilvexpanded
    # ilv_for_chain
    # active_site_for_chain
    # all_ligand_binding_for_chain
    # dna_binding_for_chain
    # interfacing_for_chain
    # interfacing_for_chain_all
    # anm_characteristics
