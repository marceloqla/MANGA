import sys
import re
import pickle
import os
from subprocess import Popen, PIPE
from hmmer3_phmmer import main_run_phmmer
import urllib.parse
import urllib.request
import textwrap

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

def validate_email(email):
    email_regex = r'\b[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Z|a-z]{2,}\b'
    is_valid_email = False
    if re.fullmatch(email_regex, email):
        is_valid_email = True
    return is_valid_email

def get_single_header(seqfile):
    seqheader = ""
    with open(seqfile, "r") as seqfile_handle:
        seqheader = seqfile_handle.readlines()[0]
        seqheader = seqheader.rstrip().split(">")[1]
    return seqheader

def run_phmmer_api(email, seqfile, outfile, db="uniprotrefprot", ince=0.0001, nhits=5000, quiet=False):
    quiet_str = ""
    if quiet:
        quiet_str = " --quiet"
    arguments_phmmer = f"--email {email} --database {db} --incE {ince} --nhits {nhits} --sequence {seqfile} --outfile {outfile}{quiet_str}"
    main_run_phmmer(arguments_phmmer.split())

def parse_phmmer_results(idfile, verbose=False):
    ids_list = []
    with open(idfile, "r") as uniprot_result_read:
        lines = uniprot_result_read.readlines()
        # ids_list = [line.strip()[2:] for line in lines]
        start_idx = 2
        if lines[0][0] == ">":
            start_idx = 1
        elif lines[0][2] ==  ":":
            start_idx = 3
        ids_list = [line.strip()[start_idx:] for line in lines]
    if verbose:
        print("Number of homologs:")
        print(len(ids_list))
        print("List of homologs:")
        print(ids_list)
    return ids_list

def get_protein_sequences(jobid, list_index, uniprot_list):
    """Retrieves the sequences from the UniProt database based on the list of
    UniProt ids.
    In general, 
        1. Compose your query here with the advanced search tool:
    https://www.uniprot.org/uniprot/?query=id%3Ap40925+OR+id%3Ap40926+OR+id%3Ao43175&sort=score
        2. Replace `&sort=score` with `&format=fasta`
        3. Edit this function as necessary
    Returns:
        protein_dict (dict): the updated dictionary
    """
    make_output_dir("./buffers/")
    uniprot_file = f"./buffers/{jobid}_{list_index}.fasta"
    fasta = None
    if os.path.exists(uniprot_file) == False:
        list_length = len(uniprot_list)
        uniprot_list = ['id%3A'+id for id in uniprot_list]
        line = '+OR+'.join(uniprot_list)
        url = f'https://www.uniprot.org/uniprot/?query={line}&format=fasta'
        with urllib.request.urlopen(url) as f:
            fasta = f.read().decode('utf-8').strip()
        if not fasta:
            raise Exception("could not download list")
        with open(uniprot_file, "w") as uniprot_file_h:
            uniprot_file_h.write(fasta)
    else:
        with open(uniprot_file, "r") as uniprot_file_h:
            fasta = uniprot_file_h.read()
    return fasta

def divide_chunks(l, n):
    # looping till length l
    for i in range(0, len(l), n): 
        yield l[i:i + n]

def retrieve_homolog_sequences(jobid, homolog_list, phmmer_rawfile, seqfile, verbose=False):
    uniprot_list = list(divide_chunks(homolog_list, 300))
    full_fasta = ""
    for i,l in enumerate(uniprot_list):
        if verbose:
            print(f"Downloading list: {i+1} of {len(uniprot_list)}...")
        # l = [a[1:] for a in l]
        fasta = get_protein_sequences(jobid, i+1, l)
        full_fasta += fasta
    with open(phmmer_rawfile, "w") as full_fasta_file_h:
        full_fasta_file_h.write(full_fasta)
    seqheader = ""
    seq = ""
    with open(seqfile, "r") as seqfile_read:
        lines = seqfile_read.readlines()
        seqheader = lines[0]
        seqheader = seqheader.rstrip().split(">")[1]
        seq = "".join([line.rstrip() for line in lines[1:]])
    clean_alignment(phmmer_rawfile, phmmer_rawfile)
    nmsa = {}
    nmsa[seqheader] = seq
    h, msa = readFastaMSA(phmmer_rawfile, mode="dict")
    for k, v in msa.items():
        if k not in nmsa:
            nmsa[k] = v
    writeFasta(nmsa, phmmer_rawfile)

def align_by_mafft_linsi(inputfilepath, outfilepath, mafft_path, verbose=False):
    # os.system(f"{MAFFT_PATH} --maxiterate 1000 --localpair  {inputfilepath} > {outfilepath}")
    p = Popen(f"{mafft_path} --maxiterate 1000 --localpair  {inputfilepath} > {outfilepath}", shell=True, stdout=PIPE, stderr=PIPE, stdin=PIPE)
    stdout, stderr = p.communicate()
    if verbose:
        print(f"{stdout.decode('UTF-8')}")
        print(f"{stderr.decode('UTF-8')}")

def align_by_prank(inputfilepath, outfilepath, prank_path, verbose=False):
    # os.system(f"{PRANK_PATH} -d={inputfilepath} -o={outfilepath}")
    p = Popen(f"{prank_path} -d={inputfilepath} -o={outfilepath}", shell=True, stdout=PIPE, stderr=PIPE, stdin=PIPE)
    stdout, stderr = p.communicate()
    if verbose:
        print(f"{stdout.decode('UTF-8')}")
        print(f"{stderr.decode('UTF-8')}")

def align_by_mafft_auto(inputfilepath, outfilepath, mafft_path, verbose=False):
    # os.system(f"{MAFFT_PATH} --auto  {inputfilepath} > {outfilepath}")
    p = Popen(f"{mafft_path} --auto  {inputfilepath} > {outfilepath}", shell=True, stdout=PIPE, stderr=PIPE, stdin=PIPE)
    stdout, stderr = p.communicate()
    if verbose:
        print(f"{stdout.decode('UTF-8')}")
        print(f"{stderr.decode('UTF-8')}")

def align_homologs(phmmer_rawfile, phmmer_msafile, mafft_path, prank_path, program="mafftauto", verbose=False):
    if program == "mafftlinsi":
        align_by_mafft_linsi(phmmer_rawfile, phmmer_msafile, mafft_path, verbose=verbose)
    elif program == "prank":
        align_by_prank(phmmer_rawfile, phmmer_msafile, prank_path, verbose=verbose)
    elif program == "mafftauto":
        align_by_mafft_auto(phmmer_rawfile, phmmer_msafile, mafft_path, verbose=verbose)

def readFastaMSA(inputfile, mode="list", clean=False, replacement="-"):
    """
    Read Fasta formatted MSA file, outputs dict with header : sequence
    """
    if mode == "list":
        header_list = []
        msa = []
    elif mode == "dict":
        msa = {}
    fr = open(inputfile,'r')
    header = None
    sequence = ""
    for line in fr:
        line = line.strip()
        if line.startswith('>'):
            if header is not None:
                if clean:
                    sequence = sequence.upper()
                    sequence = sequence.replace(".", "-")
                    sequence = sequence.replace("B",replacement)
                    sequence = sequence.replace("J",replacement)
                    sequence = sequence.replace("O",replacement)
                    sequence = sequence.replace("U",replacement)
                    sequence = sequence.replace("X",replacement)
                    sequence = sequence.replace("Z",replacement)
                if mode == "list":
                    header_list.append(header)
                    msa.append(sequence)
                elif mode == "dict":
                    msa[header] = sequence
                header = line[1:]
                sequence = ""
            else:
                header = line[1:]
        else:
            sequence += line
    if header is not None:
        if clean:
            sequence = sequence.upper()
            sequence = sequence.replace(".", "-")
            sequence = sequence.replace("B",replacement)
            sequence = sequence.replace("J",replacement)
            sequence = sequence.replace("O",replacement)
            sequence = sequence.replace("U",replacement)
            sequence = sequence.replace("X",replacement)
            sequence = sequence.replace("Z",replacement)
        if mode == "list":
            header_list.append(header)
            msa.append(sequence)
        elif mode == "dict":
            msa[header] = sequence
    fr.close()
    if mode == "dict":
        header_list = list(msa.keys())
    return header_list, msa

def writeFasta(msa_dict, outname):
    fo = open(outname, "w+")
    z = 1
    msa_seq_num = len(msa_dict.keys())
    for header_f, sequence in msa_dict.items():
        fo.write('>'+header_f+'\n'+'\n'.join(textwrap.wrap(sequence, 60))+'')
        if z < msa_seq_num:
            fo.write('\n')
        z += 1
    fo.close()

def check_alignment(phmmer_msafile, minseq=0, mincol=0):
    header_list, msa = readFastaMSA(phmmer_msafile)
    checked = False
    if len(header_list) > minseq and len(msa[0]) > mincol:
        checked = True
    return checked

def clean_alignment(phmmer_msafile, phmmer_msafile_clean, replacement="-"):
    header_list, msa = readFastaMSA(phmmer_msafile, mode="dict", clean=True, replacement=replacement)
    writeFasta(msa, phmmer_msafile_clean)

def writeUnalignedFasta(msa_dict, fname="unaligned.fa",outputdir="/", check=False):
    """
    From global MSA var, writes unaligned fasta file with name defined in fname variable (default: 'unaligned.fa')
    Input variable check verifies for existance of the fname variable named file before creating a new one
    """
    unal_file = outputdir + "/" + fname
    if check and os.path.exists(unal_file):
        return unal_file
    fw = open(unal_file,'w')
    for seqname,sequence in msa_dict.items():
        fw.write('>' + seqname + "\n")
        fw.write(sequence.replace('.','').replace('-','') + "\n")
    fw.close()
    return unal_file

def filter_msa_by_maxid_cdhit(msa_dict, alnprefix, cd_hit_path, keep=None, maxid=None, outputdir="."):
    """
    Returns new multiple sequence alignment filtered by the cd-hit program
    Cd-hit path should be defined as CD_HIT_PATH

    Input variable check verifies for the 'cluster' file existance in the output directory before running cd-hit
    Before running Cd-hit, an unaligned fasta file is created (see above function)
    Actual sequences are retrieved from MSA by the 'cluster' file indexation
    """
    if not maxid:
        maxid = 0.7
    unal_file = writeUnalignedFasta(msa_dict, outputdir=outputdir)
    make_output_dir(outputdir + "/CD-HIT_CLUSTER")
    out_file = outputdir + "/CD-HIT_CLUSTER/" + alnprefix + "_" + str(maxid)
    n = 5
    if maxid > 0.7:
        n = 5
    elif maxid > 0.6:
        n = 4
    elif maxid > 0.5:
        n = 3
    else:
        n = 2

    if os.path.exists(out_file) == False:
        os.system(cd_hit_path + ' -i ' + unal_file + ' -o ' + out_file + ' -c ' + str(maxid) + ' -n ' + str(n) + ' -M 3000 -T 2')

    new_msa = {}
    if keep:
        new_msa[keep] = msa_dict[keep]
    fr = open(out_file)
    for line in fr:
        line = line.strip()
        if len(line) > 1:
            if line[0] == '>':
                seqname = line[1:]
                sequence = msa_dict[seqname]
                new_msa[seqname] = sequence
    fr.close()
    return new_msa

def maxid_filter_keep_main(phmmer_msafile_clean, phmmer_msafile_maxid, header_seq, cd_hit_path, min_fseqs=500):
    header_list, msa = readFastaMSA(phmmer_msafile_clean, mode="dict", clean=False)
    if len(header_list) > min_fseqs:
        new_msa = filter_msa_by_maxid_cdhit(msa, header_seq, cd_hit_path, keep=header_seq, maxid=0.7)
        if len(new_msa.keys()) > min_fseqs:
            msa = new_msa
    writeFasta(msa, phmmer_msafile_maxid)

def filter_for_header(msa_dict, header):
    alnseq = msa_dict[header]
    to_keep_columns = [ipos for ipos, pos in enumerate(alnseq) if pos not in ["-", "."]]
    print("generated to keep columns")
    new_msa_dict = {}
    v = 1
    for e_header, e_seq in msa_dict.items():
        print(f"filtering: {v} of {len(msa_dict.keys())}")
        e_seq_filter = [pos for ipos, pos in enumerate(list(e_seq)) if ipos in to_keep_columns]
        e_seq_filter = "".join(e_seq_filter)
        new_msa_dict[e_header] = e_seq_filter
        v += 1
    return new_msa_dict

def seq_column_filtering(phmmer_msafile_maxid, phmmer_msafile_refseq, header_seq):
    header_list, msa = readFastaMSA(phmmer_msafile_maxid, mode="dict", clean=False)
    msa = filter_for_header(msa, header_seq)
    writeFasta(msa, phmmer_msafile_refseq)

def main_run_get_homologs(arguments, verbose=False):
    MAFFT_PATH = "EDIT ME"
    PRANK_PATH = "EDIT ME"
    CD_HIT_PATH = 'EDIT ME'
    #Email SeqFastaFile JobId
    if len(arguments) != 4:
        return([0, "Missing arguments", None])
    email = arguments[1]
    is_valid_email = validate_email(email)
    if not is_valid_email:
        return([0, "Invalid email", None])
    seqfile = arguments[2]
    jobid = arguments[3]
    header_seq = get_single_header(seqfile)
    if verbose:
        print("arguments")
        print(arguments)
        print("email")
        print(email)
        print("seqfile")
        print(seqfile)
        print("header_seq")
        print(header_seq)
    make_output_dir(f"./PHMMER_JOBS/{jobid}/homologs/")
    make_output_dir(f"./PHMMER_JOBS/{jobid}/msa/")
    phmmer_outfile = f"./PHMMER_JOBS/{jobid}/homologs/{jobid}"
    phmmer_idfile = f"./PHMMER_JOBS/{jobid}/homologs/{jobid}.ids.txt"
    if os.path.exists(phmmer_idfile) == False:
        notverbose = verbose and False
        run_phmmer_api(email, seqfile, phmmer_outfile, quiet=notverbose)
    homolog_list = parse_phmmer_results(phmmer_idfile, verbose=verbose)
    phmmer_rawfile = f"./PHMMER_JOBS/{jobid}/msa/seqs_raw.fasta"
    retrieve_homolog_sequences(jobid, homolog_list, phmmer_rawfile, seqfile, verbose=verbose)
    # clean_insert_alignment(phmmer_rawfile, phmmer_rawfile)
    # quit()
    method = "mafftauto"
    phmmer_msafile = f"./PHMMER_JOBS/{jobid}/msa/seqs_aln_{method}.fasta"
    if os.path.exists(phmmer_msafile) == False:
        align_homologs(phmmer_rawfile, phmmer_msafile, MAFFT_PATH, PRANK_PATH, program=method, verbose=verbose)
    has_seqs = check_alignment(phmmer_msafile)
    if has_seqs == False:
        os.remove(phmmer_msafile)
    phmmer_msafile_clean = f"./PHMMER_JOBS/{jobid}/msa/seqs_aln_{method}_cln.fasta"
    clean_alignment(phmmer_msafile, phmmer_msafile_clean)
    phmmer_msafile_maxid = f"./PHMMER_JOBS/{jobid}/msa/seqs_aln_{method}_cln_maxid.fasta"
    maxid_filter_keep_main(phmmer_msafile_clean, phmmer_msafile_maxid, header_seq, CD_HIT_PATH)
    phmmer_msafile_refseq = f"./PHMMER_JOBS/{jobid}/msa/seqs_aln_{method}_cln_maxid_refseq.fasta"
    seq_column_filtering(phmmer_msafile_maxid, phmmer_msafile_refseq, header_seq)
    return([1, "Homologs aligned, MSA cleaned", [phmmer_msafile_refseq]])

if __name__ == "__main__":
    arguments = sys.argv
    status, status_string, results = main_run_get_homologs(arguments)
    if results[0] == 0:
        raise Exception(results[1])
