import sys
import re
import requests
import pickle
import os
import urllib.request
import textwrap
import Bio
from Bio.PDB.PDBList import PDBList
from Bio.Data import SCOPData
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist

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

def validate_uniprot_code(uniprot_code):
    is_valid_code = False
    if re.match(r'[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', uniprot_code):
        is_valid_code = True
    return is_valid_code

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

def get_pdb_chain_from_uniprot_code(uniprot_code):
    url_uniprot_code = f"https://www.uniprot.org/uniprot/{uniprot_code}.txt"
    uniprot_txt_data = requests.get(url_uniprot_code)
    uniprot_txt_data = uniprot_txt_data.text
    uniprot_txt_data_split = uniprot_txt_data.split("\n")
    pdb_lines = [line for line in uniprot_txt_data_split if line.startswith("DR   PDB;")]
    
    pdb_codes = [line.split(";")[1].replace(" ", "").upper() for line in pdb_lines]
    dettypes = [line.split(";")[2] for line in pdb_lines]
    resolutions = [line.split(";")[3].replace(" ", "").replace("A", "") for line in pdb_lines]
    chains = [line.split(";")[-1].split("=")[0].replace(" ", "") for line in pdb_lines]

    # index_out = [i for i,e in enumerate(dettypes) if dettype in e]
    index_out = list(range(len(pdb_codes)))
    pdb_chain_list = []
    for i in index_out:
        pdb_code = pdb_codes[i]
        dettype = dettypes[i]
        resolution = resolutions[i]
        chain = chains[i]
        if "/" not in chain:
            pdb_chain_list.append((f"{pdb_code}:{chain}", pdb_code, chain, dettype, resolution))
        else:
            chain_list = chain.split("/")
            for chain in chain_list:
                pdb_chain_list.append((f"{pdb_code}:{chain}", pdb_code, chain, dettype, resolution))
    return pdb_chain_list

def get_sequences_from_pdb_chain_list(pdb_chain_list, download_dir="./PDB", verbose=False):
    pdb_chain_seq_dict = {}
    for pdb_chain_tuple in pdb_chain_list:
        pdb_chain = pdb_chain_tuple[0]
        pdb_code = pdb_chain.split(":")[0]
        selchain = pdb_chain.split(":")[1]

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
                # break
            chain_to_seqs[chain.id].append(pdbsequence)
            chain_to_seqs[chain.id].append(pdbsequencenum)
            chain_to_seqs[chain.id].append(pdbsequencefull)
        if selchain in chain_to_seqs.keys():
            pdb_chain_seq_dict[pdb_chain] = chain_to_seqs[selchain]
    return pdb_chain_seq_dict

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

def map_uniprot_pdb_residues(uniprot_sequence, pdb_chain_seq_dict):
    uniprot_pdb_chain_mappings_dict = {}
    uniprot_pdb_chain_mutations_dict = {}
    for pdb_chain, pdb_sequence_tuple in pdb_chain_seq_dict.items():
        pdb_sequence_singleaa = pdb_sequence_tuple[0]
        pdb_sequence_fullids = pdb_sequence_tuple[2]

        done, results = sw_align_seqs(uniprot_sequence, pdb_sequence_singleaa)
        mappings_dict = {}
        mutations_dict = {}
        if done == 1:
            uniprot_aln, pdbchain_aln, begin, end, score  = results

            uniprot_aaindex_list = get_index_if_matched(uniprot_aln, pdbchain_aln, "query")
            pdb_aaindex_list = get_index_if_matched(uniprot_aln, pdbchain_aln, "match")
            pos_index_list = get_index_if_matched(uniprot_aln, pdbchain_aln, "")

            if len(pdb_sequence_singleaa) != len(pdb_sequence_fullids):
                #mapping errors
                pass
            for i, current_aa_index in enumerate(pos_index_list):
                aaindex_uniprot = uniprot_aaindex_list[i]
                aaindex_pdb = pdb_aaindex_list[i]
                if uniprot_aln[current_aa_index] != pdbchain_aln[current_aa_index]:
                    # mut_list.append((pdb_sequence_fullids[aaindex_pdb, aaindex_pdb, aaindex_uniprot))
                    mutations_dict[aaindex_uniprot] = [pdb_sequence_fullids[aaindex_pdb], aaindex_pdb]
                mappings_dict[aaindex_uniprot] = [pdb_sequence_fullids[aaindex_pdb], aaindex_pdb]
        uniprot_pdb_chain_mutations_dict[pdb_chain] = mutations_dict
        uniprot_pdb_chain_mappings_dict[pdb_chain] = mappings_dict
    return uniprot_pdb_chain_mappings_dict, uniprot_pdb_chain_mutations_dict

def main_run_get_seqstr(arguments, verbose=False):
    #UniProtCode
    uniprot_code = None
    if len(arguments) < 2 or len(arguments) > 3:
        return([0, "No UniProt code", None])
    uniprot_code = arguments[1]
    if not uniprot_code:
        return([0, "Invalid UniProt code", None])
    is_valid_code = validate_uniprot_code(uniprot_code)
    if not is_valid_code:
        return([0, "Invalid UniProt code", None])
    uniprot_sequence = get_sequence_from_uniprot_code(uniprot_code)
    pdb_chain_list = get_pdb_chain_from_uniprot_code(uniprot_code)
    pdb_chain_seq_dict = get_sequences_from_pdb_chain_list(pdb_chain_list)
    uniprot_pdb_chain_mappings_dict, uniprot_pdb_chain_mutations_dict = map_uniprot_pdb_residues(uniprot_sequence, pdb_chain_seq_dict)
    if verbose:
        print("uniprot_sequence")
        print(uniprot_sequence)
        print("pdb_chain_list")
        print(pdb_chain_list)
        print("pdb_chain_seq_dict")
        print(pdb_chain_seq_dict)
        print("uniprot_pdb_chain_mappings_dict")
        print(uniprot_pdb_chain_mappings_dict)
        print("uniprot_pdb_chain_mutations_dict")
        print(uniprot_pdb_chain_mutations_dict)
    results = [uniprot_sequence,pdb_chain_list, pdb_chain_seq_dict, uniprot_pdb_chain_mappings_dict, uniprot_pdb_chain_mutations_dict]
    return(1, "UniProt validated and parsed", results)

if __name__ == "__main__":
    arguments = sys.argv
    status, status_string, results = main_run_get_seqstr(arguments)
    if results[0] == 0:
        raise Exception(results[1])
    pickle_writer(f"{arguments[1]}_seqstr.pickle", results)
