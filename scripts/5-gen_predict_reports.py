#TODO
import sys
import pickle
import os
import subprocess
import matplotlib.pyplot as plt
import shap
import Bio
import numpy as np
from Bio.PDB.PDBList import PDBList
from sklearn.ensemble import RandomForestRegressor
import shutil

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

def get_pdb_path(pdb_chain, download_dir="./PDB"):
    pdb_code = pdb_chain.split(":")[0]
    pdb_file_name = "pdb" + pdb_code.lower() + ".ent"
    structure = None
    try:
        structure = Bio.PDB.PDBParser().get_structure(pdb_code, download_dir + "/" + pdb_file_name)
    except FileNotFoundError:
        pdbl.retrieve_pdb_file(pdb_code,pdir=download_dir,file_format="mmCif")
        pdb_file_name = pdb_code.lower() + ".cif"
    return download_dir + "/" + pdb_file_name

def parse_edit_file(to_parse_file):
    residues_to_add = None
    residues_to_remove = None
    with open(to_parse_file, "r") as to_parse_handle:
        all_lines = to_parse_handle.read()
        lines_for_add = all_lines.split("#ADD")[1].split("#REMOVE")[0]
        residues_to_add = [res.rstrip() for res in lines_for_add]
        lines_for_remove = all_lines.split("#REMOVE")[1]
        residues_to_remove = [res.rstrip() for res in lines_for_remove]
    return residues_to_add, residues_to_remove

def add_functional_mappings(residue_information, functional_mappings):
    residue_information.mapped_annotations = []
    for functional_mapping in functional_mappings:
        to_add_list = getattr(residue_information, functional_mapping)
        residue_information.mapped_annotations.extend(to_add_list)
    residue_information.mapped_annotations = list(set(residue_information.mapped_annotations))
    return residue_information

def filter_residue_info(comm_add, comm_remove, residue_information, mapping_name):
    old_mapping = getattr(residue_information, mapping_name)
    new_mapping = [aapos for aapos in old_mapping if aapos not in comm_remove]
    new_mapping.extend(comm_add)
    new_mapping = list(set(new_mapping))
    setattr(residue_information, mapping_name)
    return residue_information

def generate_distance_to_sites(residue_information, mapping_name, to_map_name, pseudoval=0.1):
    site_mindist = {}
    sites_list = getattr(residue_information, mapping_name)
    if len(sites_list) == 0 and mapping_name == "mapped_communities":
        for pdbres in residue_information.pdbsequencenum:
            site_mindist[pdbres] = 20.0
            setattr(residue_information, to_map_name, site_mindist)
            return residue_information
    sites_keys = [k for k in residue_information.distance_matrix.keys() if k.split("-")[0] in sites_list or k.split("-")[1] in sites_list]
    for pdbres in residue_information.pdbsequencenum:
        if pdbres in sites_list:
            site_mindist[pdbres] = pseudoval
            continue
        pdbres_keys = [(k, residue_information.distance_matrix[k]) for k in sites_keys if pdbres in k]
        pdbres_keys = sorted(pdbres_keys, key= lambda x: x[1])
        if len(pdbres_keys) > 0:
            site_mindist[pdbres] = pdbres_keys[0][1]
    setattr(residue_information, to_map_name, site_mindist)
    return residue_information

def generate_mutation_list(uniprot_sequence, uniprot_pdb_chain_mappings_dict, uniprot_pdb_chain_mutations_dict, residue_information, conservation_dict, residue_sequence_mapping):
    mutation_list = [
        "A","C","D","E","F",
        "G","H","I","K","L",
        "M","N","P","Q","R",
        "S","T","V","W","Y"
    ]
    
    mqla_matrix = {'A': {'C': -0.21585091784198246, 'S': -0.0860983981693364, 'V': -0.1524179104477612, 'M': -0.2624752021974669, 'L': -0.33493062433826293, 'T': -0.14255375866662112, 'F': -0.3656030175210902, 'H': -0.3209988729123262, 'I': -0.3602419142782014, 'Y': -0.3603582089552239, 'Q': -0.2593250444049734, 'G': -0.10834547612342599, 'N': -0.3355509320411594, 'R': -0.34586716015609287, 'E': -0.34718592157342276, 'W': -0.4259880970547841, 'K': -0.3756966570210174, 'D': -0.3174327754757212, 'P': -0.3744044588785911}, 'C': {'S': -0.22942322629022321, 'V': -0.2353670293348024, 'A': -0.1994029850746269, 'M': -0.28057702802428486, 'L': -0.3042091836734694, 'T': -0.3543500572375953, 'F': -0.30699447804364965, 'H': -0.4670008354218879, 'I': -0.3667839476415198, 'Y': -0.28048222188310695, 'Q': -0.4419349916068976, 'G': -0.311965811965812, 'N': -0.4646567164179105, 'R': -0.4483418367346939, 'E': -0.5, 'W': -0.38137313432835823, 'K': -0.5988734064630892, 'D': -0.4918504337853231, 'P': -0.5471698113207547}, 'D': {'C': -0.207134003580774, 'S': -0.16905034324942791, 'V': -0.2345555517606476, 'A': -0.1588075880758808, 'M': -0.2547609075094569, 'L': -0.2947101163651884, 'T': -0.20353238015138772, 'F': -0.2766292040763302, 'H': -0.1755956479823716, 'I': -0.3441931255095638, 'Y': -0.2395762352158401, 'Q': -0.24187487846540612, 'G': -0.1454925373134328, 'N': -0.10076653480782456, 'R': -0.2497897392767031, 'E': -0.10604599917366753, 'W': -0.36934254317064585, 'K': -0.2674443888139305, 'P': -0.4030104442667533}, 'E': {'C': -0.1560354691075515, 'S': -0.12506738544474394, 'V': -0.14702517162471396, 'A': -0.10834732183732644, 'M': -0.20050125313283207, 'L': -0.21641682324718453, 'T': -0.13373134328358216, 'F': -0.2515820895522388, 'H': -0.13488520408163274, 'I': -0.22823880597014928, 'Y': -0.1788785660553866, 'Q': -0.093948491611478, 'G': -0.1010042265787766, 'N': -0.11894616019318249, 'R': -0.16688552731856682, 'W': -0.22624123747175712, 'K': -0.1793141440572924, 'D': -0.1064690026954178, 'P': -0.36813930348258705}, 'F': {'C': -0.24130567106623613, 'S': -0.33806146572104023, 'V': -0.25216714824926656, 'A': -0.385957424894914, 'M': -0.3339867091635971, 'L': -0.14372957127987412, 'T': -0.43813297072987467, 'H': -0.38692047572521776, 'I': -0.26070789147500345, 'Y': -0.17090961098398172, 'Q': -0.4598476675716504, 'G': -0.48396029618427017, 'N': -0.488955223880597, 'R': -0.5563346228239845, 'E': -0.5286128490767588, 'W': -0.261864794750496, 'K': -0.5136578666259729, 'D': -0.5710050400237178, 'P': -0.5702307294596756}, 'G': {'C': -0.2321083172147003, 'S': -0.20162512050681725, 'V': -0.33355795148247985, 'A': -0.22029020799890714, 'M': -0.3482388059701493, 'L': -0.3940416365993802, 'T': -0.3604179104477612, 'F': -0.3913092704126325, 'H': -0.2859922178988326, 'I': -0.48868143985516116, 'Y': -0.40896531653749635, 'Q': -0.31105752769174055, 'N': -0.3042985074626866, 'R': -0.23002130529252168, 'E': -0.30236670037932273, 'W': -0.37949529720397546, 'K': -0.4021059056920495, 'D': -0.3078805970149254, 'P': -0.4729131695406684}, 'H': {'C': -0.23733325791410081, 'S': -0.18163004361516943, 'V': -0.24841791044776124, 'A': -0.2615643397813289, 'M': -0.25403986694142955, 'L': -0.2440774804905239, 'T': -0.22807412216390138, 'F': -0.20727172565762286, 'I': -0.2932103016113483, 'Y': -0.1431405463146651, 'Q': -0.1119680484781711, 'G': -0.317507706411296, 'N': -0.1860697940503433, 'R': -0.16337313432835826, 'E': -0.22897938206541424, 'W': -0.2740669329293486, 'K': -0.3399979157052608, 'D': -0.19893448082993015, 'P': -0.37045148247978443}, 'I': {'C': -0.25011445139630706, 'S': -0.3119769079132639, 'V': -0.07703970557744144, 'A': -0.3661852875821965, 'M': -0.15698206278754906, 'L': -0.11671641791044772, 'T': -0.19764179104477614, 'F': -0.22829238516709904, 'H': -0.3844022169437846, 'Y': -0.3618667823124713, 'Q': -0.4544005102040817, 'G': -0.4593134328358209, 'N': -0.3781431334622825, 'R': -0.4562885159143506, 'E': -0.5160994964138563, 'W': -0.4593031262911445, 'K': -0.4719975583702122, 'D': -0.5484510911033115, 'P': -0.5659701492537313}, 'K': {'C': -0.22326247505431632, 'S': -0.19207601556353382, 'V': -0.21058640492753922, 'A': -0.17520215633423175, 'M': -0.14645308924485126, 'L': -0.21861787926410398, 'T': -0.14639856769040085, 'F': -0.23507757033920584, 'H': -0.1810481264556113, 'I': -0.24690026954177896, 'Y': -0.2276773712637779, 'Q': -0.1202971684889263, 'G': -0.23247761194029853, 'N': -0.0924115135656246, 'R': -0.06513432835820898, 'E': -0.14395060571959328, 'W': -0.31478509398248267, 'D': -0.310179104477612, 'P': -0.43610752331187397}, 'L': {'C': -0.16544944731628453, 'S': -0.3225441683682329, 'V': -0.14202334630350189, 'A': -0.2555318174881734, 'M': -0.07578904956371812, 'T': -0.2947252444566864, 'F': -0.1583237986270023, 'H': -0.3625820235006868, 'I': -0.14202172096908938, 'Y': -0.3259575766824356, 'Q': -0.33665753780630847, 'G': -0.44085383931218497, 'N': -0.45464852126965993, 'R': -0.4252032520325204, 'E': -0.4699372058979361, 'W': -0.3558352402745996, 'K': -0.4699656135826349, 'D': -0.5779587391877039, 'P': -0.46841305907514164}, 'M': {'C': -0.17059948979591846, 'S': -0.32288694295742193, 'V': -0.088177108038281, 'A': -0.2498278474039389, 'L': -0.06794029850746275, 'T': -0.10591448332013226, 'F': -0.21698906644238855, 'H': -0.28317092433998803, 'I': -0.11935531531579802, 'Y': -0.33080842859110315, 'Q': -0.26147959183673475, 'G': -0.3163476105219667, 'N': -0.34260154738878146, 'R': -0.3081031588585381, 'E': -0.44436554173041953, 'W': -0.407942606959462, 'K': -0.25619954901158715, 'D': -0.4034983694040913, 'P': -0.500104278223262}, 'N': {'C': -0.23041857455252984, 'S': -0.08274872448979598, 'V': -0.25322388059701495, 'A': -0.12104848123984951, 'M': -0.14878647932352013, 'L': -0.2893981607978414, 'T': -0.17056879217738602, 'F': -0.25606469002695426, 'H': -0.09661989795918367, 'I': -0.3038149015287151, 'Y': -0.16817910447761195, 'Q': -0.13410465422035234, 'G': -0.1541110039939403, 'R': -0.16663019384715305, 'E': -0.2126913265306123, 'W': -0.30699447804364965, 'K': -0.09048410497316273, 'D': -0.09255850644228236, 'P': -0.410249614404967}, 'P': {'C': -0.2418478260869565, 'S': -0.15617848970251713, 'V': -0.2072368421052632, 'A': -0.16013673415724428, 'M': -0.28119211071156147, 'L': -0.22149577425693395, 'T': -0.186751136207134, 'F': -0.32322388059701496, 'H': -0.18306404634509396, 'I': -0.2621241872590706, 'Y': -0.31761517615176155, 'Q': -0.20280353918963784, 'G': -0.23551212938005395, 'N': -0.2769647696476965, 'R': -0.242219173673412, 'E': -0.2955522388059702, 'W': -0.3951658234963462, 'K': -0.25795424664738364, 'D': -0.32532413874515514}, 'Q': {'C': -0.11431343511346065, 'S': -0.07301492537313438, 'V': -0.15280597014925373, 'A': -0.07666999392714263, 'M': -0.10724657176613228, 'L': -0.084453722009112, 'T': -0.09085492749731472, 'F': -0.11874626865671645, 'H': -0.0323582089552239, 'I': -0.14920071047957373, 'Y': -0.15039250791901948, 'G': -0.13150920195769122, 'N': -0.14317602040816335, 'R': -0.06618943230277191, 'E': -0.07653731343283587, 'W': -0.22068077803203656, 'K': -0.0949213401158211, 'D': -0.16238805970149262, 'P': -0.302129897449382}, 'R': {'C': -0.2340561224489797, 'S': -0.16908286281092633, 'V': -0.2789799344239899, 'A': -0.20647031893789114, 'M': -0.2155351319300326, 'L': -0.22690992476968258, 'T': -0.23049320111035126, 'F': -0.31438096889601347, 'H': -0.202156478404559, 'I': -0.25176119402985075, 'Y': -0.2678788054390451, 'Q': -0.22267569402089443, 'G': -0.2382331643736423, 'N': -0.2690149253731344, 'E': -0.33007398116793585, 'W': -0.2793574835718362, 'K': -0.09375, 'D': -0.4104691075514874, 'P': -0.4679919137466307}, 'S': {'C': -0.13592742773564148, 'V': -0.17367867473047593, 'A': -0.0935354691075515, 'M': -0.19565217391304346, 'L': -0.16013673415724428, 'T': -0.09149559587986139, 'F': -0.22023880597014928, 'H': -0.18620024789973832, 'I': -0.17579734472760578, 'Y': -0.2082646436012159, 'Q': -0.16175797344727616, 'G': -0.14002103602419141, 'N': -0.14087528604118993, 'R': -0.17659701492537316, 'E': -0.14817910447761196, 'W': -0.2592712574027608, 'K': -0.23219942156727724, 'D': -0.21767734553775744, 'P': -0.28043166033800937}, 'T': {'C': -0.23555491990846678, 'S': -0.10627057651779792, 'V': -0.17530226442159913, 'A': -0.11490514905149052, 'M': -0.20154142552543408, 'L': -0.2434679334916864, 'F': -0.2762985074626866, 'H': -0.18889082683733016, 'I': -0.16880341880341876, 'Y': -0.34056603773584904, 'Q': -0.16242616903427565, 'G': -0.25244456686406835, 'N': -0.20427553444180524, 'R': -0.2577313432835821, 'E': -0.2424199084668193, 'W': -0.3287164179104478, 'K': -0.2358208955223881, 'D': -0.3542840201600949, 'P': -0.2866469297612517}, 'V': {'C': -0.1179701492537314, 'S': -0.27301068891955627, 'A': -0.13133512043028595, 'M': -0.1706254948535233, 'L': -0.1294135643293829, 'T': -0.11889158870198913, 'F': -0.2405874539675676, 'H': -0.42255455516557305, 'I': -0.11295522388059707, 'Y': -0.36152036389545494, 'Q': -0.3922193877551021, 'G': -0.3141492537313433, 'N': -0.4095744680851065, 'R': -0.4942332850407191, 'E': -0.3839312184998517, 'W': -0.4123086734693878, 'K': -0.5253868471953579, 'D': -0.4704980600674815, 'P': -0.5268694650193901}, 'W': {'C': -0.2125398822277556, 'S': -0.36409214092140924, 'V': -0.4506775067750677, 'A': -0.4653116531165312, 'M': -0.374141876430206, 'L': -0.3478100442136945, 'T': -0.4515369493993448, 'F': -0.32463682182033793, 'H': -0.3124217665121555, 'I': -0.4326530612244898, 'Y': -0.2841084049603704, 'Q': -0.5103159587111105, 'G': -0.4156189555125726, 'N': -0.5476115075469593, 'R': -0.41876208897485495, 'E': -0.4671403197158082, 'K': -0.5328687877991058, 'D': -0.6012451823302698, 'P': -0.5048780487804878}, 'Y': {'C': -0.30134770889487866, 'S': -0.3238211506180376, 'V': -0.2621316399857143, 'A': -0.40091322723119105, 'M': -0.24099196845156715, 'L': -0.25795424664738364, 'T': -0.3842276904832493, 'F': -0.07468052079322302, 'H': -0.20328358208955227, 'I': -0.3082557607202808, 'Q': -0.40827709681862, 'G': -0.4902905935821512, 'N': -0.3136673584346279, 'R': -0.45106006552602845, 'E': -0.5094813212034093, 'W': -0.28758208955223885, 'K': -0.474056603773585, 'D': -0.463611859838275, 'P': -0.6088955223880598}}

    position_names = []
    mutation_names = []
    mutation_names_pdb = []
    X_predict = []
    # 'Shannon-Valdar (Min)',
    # 'Escore Matriz',
    # 'Distância Sítios Estruturais', 'Distância Comunidades',
    # 'Resistência ANM', 'Sensibilidade ANM',
    # 'Maestro ddG',
    for unipos, pdb_id in uniprot_pdb_chain_mappings_dict.items():
        pdbpos = residue_sequence_mapping[pdb_id[0]]
        if unipos in uniprot_pdb_chain_mutations_dict.keys():
            continue
        to_add = True
        restype_one = uniprot_sequence[int(unipos)]
        conservation = conservation_dict[int(unipos)]
        distance_structural = residue_information.distance_to_communities[pdbpos]
        distance_communities = residue_information.distance_to_annotations[pdbpos]
        resistence = residue_information.anm_characteristics[pdbpos]["stiff"]
        sensitivity = residue_information.anm_characteristics[pdbpos]["sensv"]
        if pdbpos not in position_names:
            position_names.append(pdbpos)
        for mut in mutation_list:
            mut_name = f"{restype_one}{unipos}{mut}"
            mut_name_pdb = f"{pdbpos}{mut}"
            if restype_one == mut or to_add == "None":
                mutation_names.append(mut_name)
                mutation_names_pdb.append(mut_name_pdb)
                X_predict.append("None")
                continue
            matriz_escore = mqla_matrix[restype_one][mut]
            X_array = [conservation, matriz_escore, distance_structural, distance_communities, resistence, sensitivity]
            mutation_names.append(mut_name)
            mutation_names_pdb.append(mut_name_pdb)
            X_predict.append(X_array)
    return position_names, mutation_names, mutation_names_pdb, X_predict

def run_maestro_predictions(maestro_command_line, maestro_config_xml, pdb_path, pdb_chain, selchain, mutation_names_pdb, position_names, X_predict, verbose=False):
    mutation_list = [
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
    X_predict_new = []
    imuts = 0
    maestro_predictions = []
    #Generate mutations file
    # for mut_name, X_array in zip(mutation_names, X_predict):
    #Run Maestro
    #Get predictions in order of mutation_names
    wtpos_muts_dict = {}
    if os.path.exists(f"./PREDICT_JOBS/MAESTRO/{pdb_chain}/ddg.pickle") == False:
        for i_pos, position_name in enumerate(position_names):
            wt_res = THREE_TO_ONE_AA_TYPES_CANONICAL[position_name[:3]]
            wt_pos = wt_res + position_name[3:]
            mutation_list_filter = [mut for mut in mutation_list if mut != wt_res]
            mutations = ",".join(mutation_list_filter)
            command_line = f"{maestro_command_line} {maestro_config_xml} {pdb_path} --evalmut='{wt_pos}.{selchain}" + "{" + f"{mutations}"+"}"+ "'"
            if verbose:
                print(wt_res)
                print(mutations)
                print(command_line)
            proc = subprocess.Popen(command_line, stdout=subprocess.PIPE, shell=True)
            (out, err) = proc.communicate()
            if err:
                if verbose:
                    print(f"MAESTRO error for position: {position_name}")
                continue
            output_lines = str(out.decode("utf-8")).split('\n')
            output_lines_results = [line for line in output_lines if "wildtype" not in line]
            output_lines_results = [line for line in output_lines_results if len(line) > 3]
            output_lines_results = [line for line in output_lines_results if line[0] != "#"]

            mut_to_ddg = {}
            for line in output_lines_results:
                mut_to_ddg[line.split()[2].split("{")[1][:-1]] = line.split()[4]
            if verbose:
                print(mut_to_ddg)
            wtpos_muts_dict[position_name] = mut_to_ddg
        make_output_dir(f"./PREDICT_JOBS/MAESTRO/{pdb_chain}/")
        pickle_writer(f"./PREDICT_JOBS/MAESTRO/{pdb_chain}/ddg.pickle", wtpos_muts_dict)
    else:
        print("loading previous maestro run...")
        wtpos_muts_dict = pickle_loader(f"./PREDICT_JOBS/MAESTRO/{pdb_chain}/ddg.pickle")
    #Push predictions to array
    for mut_name, X_array in zip(mutation_names_pdb, X_predict):
        position_name = mut_name[:-1]
        mut = mut_name[-1]
        # maestro_ddg = maestro_predictions[imuts]
        if position_name not in wtpos_muts_dict.keys():
            maestro_ddg = None
        else:
            if mut not in wtpos_muts_dict[position_name].keys():
                maestro_ddg = None
            else:
                maestro_ddg = wtpos_muts_dict[position_name][mut]
        if maestro_ddg:
            X_predict_new.append(X_array)
            X_predict_new[-1].append(float(maestro_ddg))
        else:
            X_predict_new.append("None")
        imuts += 1
    return X_predict_new

def evaluate_mutation_list(uniprot_code, pdb_chain, mutation_names, X_predict):
    column_names = [
        "Conservation", "Empirical Matrix", "Functional Distance", "Comm Distance", "Stiffness (ANM)", "Sensibility (ANM)", "ddG (MAESTRO)"
    ]
    rforest = pickle_loader("all_rf_model.pickle")
    # shap_exp = pickle_loader("all_rf_model_shap.pickle")
    valid_X_predict = np.array([np.array(data) for data in X_predict if data != "None"])

    # shap_name = f"./{uniprot_code}_{pdb_chain}_shap.pickle"
    shap_explainer = shap.TreeExplainer(rforest, valid_X_predict)
    # pickle_writer(shap_name, shap_explainer)
    # shap_explainer = pickle_loader(shap_name)

    shap_values = shap_explainer(valid_X_predict)
    shap_values.feature_names = column_names
    manga_predictions = {}
    mutation_effect_dict = {}
    results_dir = f"./RESULTS/{uniprot_code}_{pdb_chain}/"
    # make_output_dir(f"{results_dir}/SHAP/")
    make_output_dir(f"{results_dir}/")
    valid_i = 0
    for mut_name, x_array in zip(mutation_names, X_predict):
        pos_name = mut_name[:-1]
        mut = mut_name[-1]
        if pos_name not in manga_predictions.keys():
            manga_predictions[pos_name] = []
        predict = 'None'
        mut_shap = 'None'
        if x_array != 'None':
            predict = rforest.predict(np.array([x_array]))
            predict = predict[0]
            mut_shap = {"base":shap_values[valid_i].base_values, "values":list(shap_values[valid_i].values), "names":column_names}
            valid_i += 1
        manga_predictions[pos_name].append([mut, predict, mut_shap])
        mutation_effect_dict[mut_name] = [mut, predict, mut_shap]
    print(mutation_effect_dict.keys())
    print(manga_predictions.keys())
    results_file_name = f"{results_dir}/mutations_info.js"
    #TODO: export libs and templates to folder
    with open(results_file_name, "w") as results_file_handle:
        results_file_handle.write(f"var prediction_per_pos = {manga_predictions};")
        results_file_handle.write("\n")
        results_file_handle.write(f"var prediction_per_mut = {mutation_effect_dict};")
    # for valid_mut_name, valid_x_array in zip(valid_mutation_names, valid_X_predict):


def generate_reports():
    pass

def main_run_prediction(arguments, verbose=True):
    MAESTRO_COMMAND_LINE = "/home/marcelo/Documentos/Programas/MAESTRO/MAESTRO_linux_x64/./maestro"
    MAESTRO_CONFIG_XML = "/home/marcelo/Documentos/Programas/MAESTRO/MAESTRO_linux_x64/config.xml"
    allowed_mappings = ["active_site_for_chain", "all_ligand_binding_for_chain", "interfacing_for_chain", "dna_binding_for_chain"]
    functional_mappings = [
        "active_site_for_chain", #MANGA validated with this feature
        "all_ligand_binding_for_chain", #MANGA validated with this feature
        "interfacing_for_chain", #MANGA validated with this feature
        "dna_binding_for_chain", #MANGA validated with this feature
        # "interfacing_for_chain_all", #NON validated. Includes nonbonded ligplot terms
    ]
    if len(arguments) < 4:
        return([0, "Missing arguments", None])
    elif len(arguments) > 7:
        return([0, "Wrong usage", None])
    seqstrmap_file = arguments[1]
    seqstrmap = pickle_loader(seqstrmap_file)
    conscevomap_file = arguments[2]
    conscevomap = pickle_loader(conscevomap_file)
    annomap_file = arguments[3]
    annomap = pickle_loader(annomap_file)
    residue_information = annomap[0]
    uniprot_residue_mapping = {}
    residue_sequence_mapping = {}
    for res_id, pdbres in zip(residue_information.pdbsequencefull, residue_information.pdbsequencenum):
        residue_sequence_mapping[res_id] = pdbres
    
    uniprot_code = annomap_file.split("_")[0]
    pdb_chain = annomap_file.split("_")[1]
    pdb_keys = [pdbc[0] for pdbc in seqstrmap[1]]
    if pdb_chain not in pdb_keys:
        return([0, "Unmapped pdb", None])
    # print("pdb_chain")
    # print(pdb_chain)
    # print("seqstrmap")
    # print(seqstrmap[0])
    # print(seqstrmap[1])
    pdb_chain_seq_dict = seqstrmap[2][pdb_chain]
    uniprot_pdb_chain_mappings_dict = seqstrmap[3][pdb_chain]
    for unipos, pdb_id in uniprot_pdb_chain_mappings_dict.items():
        pdbpos = residue_sequence_mapping[pdb_id[0]]
        uniprot_residue_mapping[pdbpos] = unipos
    uniprot_pdb_chain_mutations_dict = seqstrmap[4][pdb_chain]
    conservation_dict = conscevomap[0]
    uniprot_sequence = seqstrmap[0]
    pdb_path = get_pdb_path(pdb_chain)
    selchain = pdb_chain.split(":")[1]
    if verbose:
        print(uniprot_sequence)
        print(selchain)
        print(pdb_path)
    if len(arguments) == 7:
        functional_mappings = []
        with open(arguments[4], "r") as new_mappings:
            lines = new_mappings.readlines()
            functional_mappings = [l.rstrip() for l in lines if l.rstrip() in allowed_mappings]
    residue_information = add_functional_mappings(residue_information, functional_mappings)
    if verbose:
        print(residue_information.mapped_annotations)
    comm_add = []
    comm_remove = []
    if len(arguments) == 5:
        comm_add, comm_remove = parse_edit_file(arguments[4])
        residue_information = filter_residue_info(comm_add, comm_remove, residue_information, "mapped_communities")
    func_add = []
    func_remove = []
    if len(arguments) == 6:
        func_add, func_remove = parse_edit_file(arguments[5])
        residue_information = filter_residue_info(func_add, func_remove, residue_information, "mapped_annotations")
    residue_information = generate_distance_to_sites(residue_information, "mapped_communities", "distance_to_communities")
    residue_information = generate_distance_to_sites(residue_information, "mapped_annotations", "distance_to_annotations")
    if verbose:
        print(residue_information.mapped_communities)
        print(residue_information.mapped_annotations)
        print(residue_information.distance_to_communities)
        print(residue_information.distance_to_annotations)
    position_names, mutation_names, mutation_names_pdb, X_predict = generate_mutation_list(uniprot_sequence, uniprot_pdb_chain_mappings_dict, uniprot_pdb_chain_mutations_dict, residue_information, conservation_dict, residue_sequence_mapping)
    if verbose:
        print(position_names[:3])
        print(mutation_names[:3])
        print(mutation_names_pdb[:3])
        print(X_predict[:3])
    X_predict = run_maestro_predictions(MAESTRO_COMMAND_LINE, MAESTRO_CONFIG_XML, pdb_path, pdb_chain, selchain, mutation_names_pdb, position_names, X_predict, verbose=verbose)
    if verbose:
        print(X_predict[:3])
    evaluate_mutation_list(uniprot_code, pdb_chain, mutation_names, X_predict)
    
    with open(f"./RESULTS/{uniprot_code}_{pdb_chain}/positions_info.js", "w") as pos_info_handle:
        pos_info_handle.write(f"var active_site_for_chain = {residue_information.active_site_for_chain};")
        pos_info_handle.write("\n")
        pos_info_handle.write(f"var all_ligand_binding_for_chain = {residue_information.all_ligand_binding_for_chain};")
        pos_info_handle.write("\n")
        pos_info_handle.write(f"var interfacing_for_chain = {residue_information.interfacing_for_chain};")
        pos_info_handle.write("\n")
        pos_info_handle.write(f"var interfacing_for_chain_all = {residue_information.interfacing_for_chain_all};")
        pos_info_handle.write("\n")
        pos_info_handle.write(f"var interfacing_per_interface = {residue_information.interacting_residues_per_interface};")
        pos_info_handle.write("\n")
        pos_info_handle.write(f"var dna_binding_for_chain = {residue_information.dna_binding_for_chain};")
        pos_info_handle.write("\n")
        pos_info_handle.write(f"var intrachain_disulphide_list = {residue_information.intrachain_disulphide_list};")
        pos_info_handle.write("\n")
        pos_info_handle.write(f"var salt_bridges = {residue_information.salt_bridges};")
        pos_info_handle.write("\n")
        pos_info_handle.write(f"var anm_characteristics = {residue_information.anm_characteristics};")
        pos_info_handle.write("\n")
        pos_info_handle.write(f"var coevolution_pos = {residue_information.mapped_communities};")
        pos_info_handle.write("\n")
        pos_info_handle.write(f"var additional_mappings = {residue_information.additional_info};")
        pos_info_handle.write("\n")
        pos_info_handle.write(f"var uniprot_residue_mapping = {uniprot_residue_mapping};")
        pos_info_handle.write("\n")
        pos_info_handle.write(f"var pdb_chain = '{pdb_chain}';")
        pos_info_handle.write("\n")
        pos_info_handle.write(f"var uniprot_code = '{uniprot_code}';")
    # generate_shaps(mutation_names, X_predict)
    shutil.copy("./template/d3.v5.min.js", f"./RESULTS/{uniprot_code}_{pdb_chain}/")
    shutil.copy("./template/d3-legend.min.js", f"./RESULTS/{uniprot_code}_{pdb_chain}/")
    shutil.copy("./template/prediction.html", f"./RESULTS/{uniprot_code}_{pdb_chain}/")
    shutil.copy("./template/prediction.js", f"./RESULTS/{uniprot_code}_{pdb_chain}/")

if __name__ == "__main__":
    arguments = sys.argv
    main_run_prediction(arguments)