import sys
import os
import numpy as np
import networkx as nx
import pickle
import json
import math
import itertools
from scipy.stats import binom
from scipy.stats import hypergeom
from scipy.special import logsumexp

ONE_AA_TYPES_CANONICAL = [
    "A", "M",
    "C", "N",
    "D", "P",
    "E", "Q",
    "F", "R",
    "G", "S",
    "H", "T",
    "I", "V",
    "K", "W",
    "L", "Y",
]

BLOSUM62_DISTRIBUTION = {
    "A":0.078,   "R":0.051,
    "N":0.041,   "D":0.052,
    "C":0.024,   "Q":0.034,
    "E":0.059,   "G":0.083,
    "H":0.025,   "I":0.062,
    "L":0.092,   "K":0.056,
    "M":0.024,   "F":0.044,
    "P":0.043,   "S":0.059,
    "T":0.055,   "W":0.014,
    "Y":0.034,   "V":0.072,
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

def calculate_henikoff_weights_lists(msa_def_name, msa, num_of_rows, num_of_cols, verbose=False):
    seq_weights = [0.] * num_of_rows
    henikoff_path = f"./CONSCOEVO_JOBS/{msa_def_name}/henikoff/"
    make_output_dir(henikoff_path)
    henikoff_save_file = henikoff_path + msa_def_name + "_hweights.pickle"
    if verbose:
        print("Weights file:")
        print(henikoff_save_file)
        print("Exists:")
        print(os.path.exists(henikoff_save_file))
    if os.path.exists(henikoff_save_file) == False:
        for col_num in range(num_of_cols):
            # column_array = get_column_from_msa(msa_dict, col_num)
            column_array = [seq[col_num] for seq in msa]
            column_non_gap = [a for a in column_array if a != "-"]
            list_of_types = list(set(column_non_gap))
            number_of_types = len(list_of_types)
            column_non_gap_str = "".join(column_non_gap)

            weighted_count_per_type = dict((c, column_non_gap_str.count(c) * number_of_types) for c in list_of_types)
            weighted_count_per_type["-"] = 0.0 #ignore gaps

            # weighted_count_list = [column_non_gap_str.count(c) * number_of_types for c in list_of_types]
            # weighted_count_list = [a for a in weighted_count_list if a > 0]
            # weighted_count_sum = sum(1/a for a in weighted_count_list)

            for row_num, pos in enumerate(column_array):
                weight_pos = weighted_count_per_type[pos]
                if weight_pos > 0:
                    seq_weights[row_num] += 1 / weight_pos
        for iw in range(len(seq_weights)):
            seq_weights[iw] /= num_of_cols
        weight_sum = sum(seq_weights)
        seq_weights = [w/weight_sum for w in seq_weights]
        pickle_writer(henikoff_save_file, seq_weights)
    else:
        seq_weights = pickle_loader(henikoff_save_file)
    seq_weights = np.array(seq_weights)
    return seq_weights

def get_conservations_lists(msa_def_name, msa, num_of_cols, seq_weights_list=[], pseudocount=0.0, distribution="blosum62", jensen_weight=0.5, bgdist_to=None, verbose=False):
    calc_type = ""
    raw_count_per_column = {}
    count_freq_per_column = {}
    count_freq_per_column = {}
    if len(list(seq_weights_list)) == 0:
        seq_weights_list = [1.] * num_of_cols
    else:
        calc_type += "henikoff"
    if pseudocount > 0:
        if calc_type:
            calc_type += "_"
        calc_type += "pseudo"
    if calc_type:
        calc_type += "_"
    calc_type += distribution

    freq_counts_path = f"./CONSCOEVO_JOBS/{msa_def_name}/freq/"
    cs_path = f"./CONSCOEVO_JOBS/{msa_def_name}/cons/"
    make_output_dir(freq_counts_path)
    make_output_dir(cs_path)
    c_rawfile_name = freq_counts_path + msa_def_name + "_raw_freq.pickle"
    cf_file_name = freq_counts_path + msa_def_name + "_" + calc_type + "_freq.pickle"
    cs_file_name = cs_path + msa_def_name + "_" + calc_type + "_cons.pickle"
    
    conservation_per_column = {}
    conservation_per_column["js"] = {}
    conservation_per_column["sh"] = {}
    conservation_per_column["js-v"] = {}
    conservation_per_column["sh-v"] = {}
    conservation_per_column["js-vm"] = {}
    conservation_per_column["sh-vm"] = {}
    conservation_per_column["gapw"] = {}

    gap_inverse_scaling_factors = []
    if os.path.exists(cs_file_name) == False:
        sum_of_weights = sum(seq_weights_list)
        for col_num in range(0, num_of_cols):
            if verbose:
                print(f"Processing column {col_num+1} of {num_of_cols}")
            # column_array = list(msa_np[:,col_num])
            # column_array = column_array_original.copy()
            column_array = [seq[col_num] for seq in msa]

            per_aa_dict_counts = {}
            bg_dist = {}
            # raw_count_per_column[col_num] = dict((pos, column_array.count(pos)) for pos in set(column_array))
            raw_count_per_column[col_num] = dict((pos, column_array.count(pos)) for pos in ONE_AA_TYPES_CANONICAL)

            gappy_array = [i for i, pos in enumerate(column_array) if pos == "-"]
            gappy_weights_list = [seq_weights_list[i] for i, pos in enumerate(column_array) if pos == "-"]
            sum_of_gap_weights = sum(gappy_weights_list)

            gap_count = len(gappy_array)
            gap_freq = len(gappy_array) / len(column_array)

            raw_count_per_column[col_num]["-"] = gap_count

            for aa_type1 in ONE_AA_TYPES_CANONICAL:
                bg_dist[aa_type1] = 1.0/20

            if seq_weights_list.shape[0]  != len(column_array):
                raise Exception("weight num - seq num differente")

            column_array_weights = [seq_weights_list[i] for i,pos in enumerate(column_array) if pos != "-"]
            column_array = [pos for pos in column_array if pos != "-"]

            for aa_type in ONE_AA_TYPES_CANONICAL:
                column_array_aai = [i for i, aa in enumerate(column_array) if aa == aa_type]
                per_aa_dict_counts[aa_type] = sum(column_array_weights[i] for i in column_array_aai)

            per_aa_dict_freqs = {}
            r_dist = {}

            jensen_shannon_divergence = 0.0
            shannon_divergence = 0.0

            if distribution == "blosum62":
                bg_dist = BLOSUM62_DISTRIBUTION
            elif distribution == "custom":
                if 0 in bgdist_to.keys():
                    bg_dist = bgdist_to[col_num]
                else:  
                    bg_dist = bgdist_to
            
            pseudocount_value = pseudocount

            jensen_shannon_divergence_marcelo = 0.0
            shannon_divergence_marcelo = 0.0
            # div_factor = 2
            if len(column_array) > 0:
                for k, v in per_aa_dict_counts.items():
                    # per_aa_dict_freqs[k] = v / (sum(seq_weights_list) + len(ONE_AA_TYPES_CANONICAL) * pseudocount)

                    per_aa_dict_freqs[k] = v / (sum_of_weights + len(ONE_AA_TYPES_CANONICAL) * pseudocount_value)
                    
                    r_dist[k] = jensen_weight * per_aa_dict_freqs[k] + (1-jensen_weight) * bg_dist[k]
                    
                    fv = per_aa_dict_freqs[k]
                    r = r_dist[k]
                    bg = bg_dist[k]
                    calc_type = "None"
                    
                    if fv > 0:
                        shannon_divergence += (fv * math.log(fv, 2))#*(-1)
                    if r != 0:
                        if fv == 0:
                            calc_type = "from bg only"
                            jensen_shannon_divergence += bg * math.log(bg/r, 2)
                        elif bg == 0.0:
                            calc_type = "from fv only"
                            jensen_shannon_divergence += fv * math.log(fv/r, 2) 
                        else:
                            calc_type = "from bg & fv"
                            jensen_shannon_divergence += fv * math.log(fv/r, 2) + bg * math.log(bg/r, 2)
                            # div_factor = 2
                    if v > 0 and verbose:
                        print(k, v, calc_type)
                # shannon_div = math.log(min(len(per_aa_dict_freqs.keys()), len(column_array)), 2)
                shannon_div = math.log(min(len(per_aa_dict_freqs.keys()), len(column_array)), 2)
                if shannon_div == 0:
                    shannon_divergence = -1.0
                else:
                    shannon_divergence /= shannon_div
                shannon_divergence = 1 - (-1 * shannon_divergence)

                # jensen_shannon_divergence /= div_factor
                jensen_shannon_divergence /= 2

                jensen_shannon_divergence_marcelo = jensen_shannon_divergence
                shannon_divergence_marcelo = shannon_divergence

                count_freq_per_column[col_num] = {"count":per_aa_dict_counts, "freq":per_aa_dict_freqs}
            gap_inverse_scaling_factors.append(1-(sum_of_gap_weights/sum_of_weights))

            conservation_per_column["js"][col_num] = jensen_shannon_divergence_marcelo
            conservation_per_column["sh"][col_num] = shannon_divergence_marcelo
        
        min_gap_scaling_factor = max(gap_inverse_scaling_factors) 
        gap_inverse_scaling_factors_rescaled = [v/min_gap_scaling_factor for v in gap_inverse_scaling_factors]
        for col_num in range(0, num_of_cols):
            conservation_per_column["js-v"][col_num] = conservation_per_column["js"][col_num] * gap_inverse_scaling_factors[col_num]
            conservation_per_column["sh-v"][col_num] = conservation_per_column["sh"][col_num] * gap_inverse_scaling_factors[col_num]
            conservation_per_column["js-vm"][col_num] = conservation_per_column["js"][col_num] * gap_inverse_scaling_factors_rescaled[col_num]
            conservation_per_column["sh-vm"][col_num] = conservation_per_column["sh"][col_num] * gap_inverse_scaling_factors_rescaled[col_num]
            conservation_per_column["gapw"][col_num] = gap_inverse_scaling_factors[col_num]
        pickle_writer(cf_file_name, count_freq_per_column)
        pickle_writer(cs_file_name, conservation_per_column)
        pickle_writer(c_rawfile_name, raw_count_per_column)
    else:
        conservation_per_column = pickle_loader(cs_file_name)
    # return conservation_per_column["js"], conservation_per_column["sh"],conservation_per_column["gapw"]
    return conservation_per_column

def calc_gap_stats_lists(msa_def_name, msa, num_of_cols, gap_type,seq_weights_list=[], verbose=False):
    gapstats_per_column = {}
    gap_counts_path = f"./CONSCOEVO_JOBS/{msa_def_name}/gap/"
    make_output_dir(gap_counts_path)
    gp_file_name = gap_counts_path + msa_def_name + "_" + gap_type + "_gap_info.pickle"
    if verbose:
        print("gap calc start!")
    if os.path.exists(gp_file_name) == False:
        sum_of_weights = sum(seq_weights_list)
        for col_num in range(0, num_of_cols):
            if verbose:
                print(f"gap calc: {col_num+1} of {num_of_cols}")
            column_array = [seq[col_num] for seq in msa]
            gappy_array = [i for i, pos in enumerate(column_array) if pos == "-"]
            gappy_weights_list = [seq_weights_list[i] for i, pos in enumerate(column_array) if pos == "-"]
            sum_of_gap_weights = sum(gappy_weights_list)
            gap_count = len(gappy_array)
            gap_freq = len(gappy_array) / len(column_array)
            gapstats_per_column[col_num] = {"gap_freq": gap_freq, "gap_weight":sum_of_gap_weights, "gap_weight_freq":sum_of_gap_weights/sum_of_weights}
        pickle_writer(gp_file_name, gapstats_per_column)
    else:
        gapstats_per_column = pickle_loader(gp_file_name)
    if verbose:
        print("gap calc done!")
    return gapstats_per_column

def cumulative_binon(calc_data):
    count_freq_ab, count_freq_a, count_freq_b, num_of_rows = calc_data
    # binon_logsum = binom.logcdf(sucessos, n_total, prob_sucesso)
    prob_b = count_freq_b / num_of_rows
    binon_logsum = binom.logcdf(count_freq_ab, num_of_rows, prob_b)
    return binon_logsum

def cumulative_hypergeom(calc_data):
    count_freq_1_2, count_freq_1, count_freq_2, num_of_rows = calc_data
    possible_num_draws = np.arange(count_freq_1_2 + 1, count_freq_2 + 1)
    hypergeom_logpmf = hypergeom._logpmf(possible_num_draws, num_of_rows, count_freq_1, count_freq_2)
    hypergeom_logsf = logsumexp(hypergeom_logpmf)
    return hypergeom_logsf

def calc_correlation(count_freq_a, count_freq_b, count_freq_ab, num_of_rows, similarity_filter, pvalue_filter, similarity_type="overlap", pvalue_type="hypergeom", weighted=False):
    ci = 0
    if weighted:
        ci = 1
    similarity_p_value = []
    similarity = None
    jaccard = count_freq_ab[ci]/(count_freq_a[ci]+count_freq_b[ci]-count_freq_ab[ci])
    overlap = count_freq_ab[ci]/min(count_freq_a[ci],count_freq_b[ci])
    if similarity_type == "jaccard":
        similarity = jaccard
    elif similarity_type == "overlap":
        similarity = overlap
    if similarity < similarity_filter:
        return []
    similarity_p_value.append(jaccard)
    similarity_p_value.append(overlap)
    # p_value_per_pair
    pvalue = None
    calc_data = (count_freq_ab[0],count_freq_a[0],count_freq_b[0])
    if weighted:
        calc_data = [count_freq_ab[0]*count_freq_ab[1],count_freq_a[0]*count_freq_a[1],count_freq_b[0]*count_freq_b[1]]
        calc_data = (math.floor(calc_data[0]),math.floor(calc_data[1]),math.floor(calc_data[2]))
    if pvalue_type == "hypergeom":
        if calc_data in hypergeom_dict:
            pvalue = hypergeom_dict[calc_data]
        else:
            data_to_calc = (count_freq_ab[0]-1,count_freq_a[0],count_freq_b[0],num_of_rows)
            pvalue = cumulative_hypergeom(data_to_calc)
            hypergeom_dict[calc_data] = pvalue
    elif pvalue_type == "binomial":
        if calc_data in binom_dict:
            pvalue = binom_dict[calc_data]
        else:
            data_to_calc = (count_freq_ab[0],count_freq_a[0],count_freq_b[0],num_of_rows)
            pvalue = cumulative_binon(data_to_calc)
            binom_dict[calc_data] = pvalue
    pvalue *= -1
    if pvalue < pvalue_filter:
        return []
    similarity_p_value.append(pvalue)
    return similarity_p_value

def get_filtered_correlations_lists(msa_def_name, freq_calc_type, msa_columns, min_freq, num_of_rows, num_of_cols, seq_weights_list, similarity_filter, pvalue_filter, from_respos=[], to_respos=[], similarity_type="overlap", pvalue_type="hypergeom", weighted_calc=True, verbose=False):
    # if min_freq >= 1:
        # raise Exception("min frequency must be between 0 and 1")

    freq_counts_path = f"./CONSCOEVO_JOBS/{msa_def_name}/freq/"
    gap_counts_path = f"./CONSCOEVO_JOBS/{msa_def_name}/gap/"
    correlations_path = f"./CONSCOEVO_JOBS/{msa_def_name}/corr/"
    make_output_dir(correlations_path)
    
    # min_count = min_freq * num_of_rows
    # max_gap_count = num_of_rows - min_count
    max_gap_freq = 1 - min_freq

    correlation_dict = {"type":freq_calc_type, "list":[]}
    correlation_dict["type"] += "_" + similarity_type
    correlation_dict["type"] += "_" + pvalue_type

    c_rawfile_name = freq_counts_path + msa_def_name + "_raw_freq.pickle"

    raw_count_per_column = pickle_loader(c_rawfile_name)
    cf_file_name = freq_counts_path + msa_def_name + "_" + freq_calc_type + "_freq.pickle"
    wf_count_per_column = pickle_loader(cf_file_name)

    gp_file_name = gap_counts_path + msa_def_name + "_" + freq_calc_type + "_gap_info.pickle"
    gap_stats_per_column = pickle_loader(gp_file_name)

    corr_file_name = correlations_path + msa_def_name + "_" + correlation_dict["type"] + "_corr.pickle"
    if os.path.exists(corr_file_name) == False:
        sum_of_weights = sum(seq_weights_list)
        seq_weights_np = np.array(seq_weights_list)
        all_rows = list(range(num_of_rows))
        all_columns = list(range(num_of_cols))
        if verbose:
            print(f"All columns: {len(all_columns)}")

        allowed_columns = [i_column for i_column in all_columns if gap_stats_per_column[i_column]["gap_weight_freq"] <= max_gap_freq]

        if verbose:
            print("Calc list of allowed res+column index pairs")
            print(f"Allowed columns: {len(allowed_columns)}")
        allowed_resids_in_columns = set()
        for i_column in allowed_columns:
            # column_array = get_column_from_msa(msa_dict, i_column)
            # column_array = msa_np[:,]
            # column_array = [seq[i_column] for seq in msa]
            column_array = msa_columns[i_column]
            res_types = set(column_array)
            if "-" in res_types:
                res_types.remove("-")
            res_types = list(res_types)

            #TODO: remake wf_count_per_column per filtered aln
            res_types = [(res_type,i_column,wf_count_per_column[i_column]["freq"][res_type]) for res_type in res_types if wf_count_per_column[i_column]["freq"][res_type] >= min_freq]
            allowed_resids_in_columns.update(res_types)
        if verbose:
            print("Done!")

        if verbose:
            print("Calc list of allowed res+column index,res+column index combinations")
        allowed_resids_combinations = sorted(list(itertools.combinations(allowed_resids_in_columns, 2)))
        if verbose:
            print("Done!")
            print("length of combinations:")
            print(len(allowed_resids_combinations))
            print("Done end!")

        dn = 0
        zx = 0
        for resid_combination in allowed_resids_combinations:
            dn += 1
            zx += 1
            if zx == 1000 or dn == len(allowed_resids_combinations):
                if verbose:
                    print(f"processing {dn} of {len(allowed_resids_combinations)}")
                    print(resid_combination)
                zx = 0
            res1 = resid_combination[0][0]
            i_column1 = resid_combination[0][1]
            res2 = resid_combination[1][0]
            i_column2 = resid_combination[1][1]

            column_array1 = np.array(msa_columns[i_column1])
            column_array2 = np.array(msa_columns[i_column2])

            has_one = column_array1 == res1
            index_one = np.where(has_one == True)[0]
            has_two = column_array2 == res2
            index_two = np.where(has_two == True)[0]

            has_both = np.logical_and(has_one, has_two)
            index_both = np.where(has_both == True)[0]

            total_num = has_one.size
            num_of_both = index_both.size

            sum_of_weights1 = np.sum(seq_weights_np[index_one])
            sum_of_weights2 = np.sum(seq_weights_np[index_two])
            sum_of_weightsB = np.sum(seq_weights_np[index_both])

            count_freq_res_a = raw_count_per_column[i_column1][res1]
            count_freq_res_a = [count_freq_res_a, sum_of_weights1]
            count_freq_res_b = raw_count_per_column[i_column2][res2]
            count_freq_res_b = [count_freq_res_b, sum_of_weights2]
            count_freq_res_ab = num_of_both
            count_freq_res_ab = [count_freq_res_ab, sum_of_weightsB]
            putative_correlation = calc_correlation(count_freq_res_a, count_freq_res_b, count_freq_res_ab, num_of_rows, similarity_filter, pvalue_filter, similarity_type=similarity_type, pvalue_type=pvalue_type, weighted=weighted_calc)
            if putative_correlation:
                parsed_correlation = JSONObject()
                parsed_correlation.res_from = res1
                parsed_correlation.pos_from = i_column1+1
                parsed_correlation.respos_from = f"{res1}{i_column1+1}"
                parsed_correlation.respos_from_fc = count_freq_res_a
                parsed_correlation.res_to = res2
                parsed_correlation.pos_to = i_column2+1
                parsed_correlation.respos_to = f"{res2}{i_column2+1}"
                parsed_correlation.respos_to_fc = count_freq_res_b
                parsed_correlation.respos_both_fc = count_freq_res_ab
                parsed_correlation.jaccard = putative_correlation[0]
                parsed_correlation.overlap = putative_correlation[1]
                parsed_correlation.p_value = putative_correlation[2]
                correlation_dict["list"].append(parsed_correlation)
        if verbose:
            print("Done!")
        pickle_writer(corr_file_name, correlation_dict)
    else:
        correlation_dict = pickle_loader(corr_file_name)
    return correlation_dict

def generate_communities_simple(correlation_list, ncols, verbose=False):
    correlation_list_as_edgelist =[(corr.respos_from, corr.respos_to, {'weight':corr.jaccard, "pvalue":corr.p_value}) for corr in correlation_list]
    G = nx.Graph()
    G.add_edges_from(correlation_list_as_edgelist)

    remove = [node for node,degree in dict(G.degree()).items() if degree < 2]
    G.remove_nodes_from(remove)

    CCG = nx.connected_components(G)
    communities = list(CCG)
    comm_more_2 = [comm for comm in communities if len(comm)>2]
    all_nodes_2_lists = []
    comm_more_2_lists = [list(comm) for comm in comm_more_2]
    for comm in comm_more_2_lists:
        all_nodes_2_lists.extend(comm)
    remove2 = [node for node in list(G.nodes) if str(node) not in all_nodes_2_lists]
    if remove2:
        G.remove_nodes_from(remove2)
    return comm_more_2_lists, G

def generate_communities(correlation_list, num_of_cols, max_col_perc=0.2, verbose=False):
    correlation_list_as_edgelist =[(corr.respos_from, corr.respos_to, {'weight':corr.jaccard, "pvalue":corr.p_value}) for corr in correlation_list]
    G = nx.Graph()
    G.add_edges_from(correlation_list_as_edgelist)
    node_num = len(list(G.nodes))
    node_num_final = node_num
    if node_num <= 30:
        CCG = nx.connected_components(G)
        communities = list(CCG)
        comm_lists = [list(comm) for comm in communities]
        if verbose:
            print(f"N cols:{num_of_cols}, N nodes: {node_num}, N nodes f: {node_num_final}")
        return comm_lists, G
    # print(f"N cols:{num_of_cols}, N nodes: {node_num}")
    #Two possibilities: very low number of correlations -> reduce jaccard
    #very high number of correlations -> raise jaccard

    comm_more_2_lists = list(G.nodes)
    remove = [node for node,degree in dict(G.degree()).items() if degree < 2]
    G.remove_nodes_from(remove)

    node_num_f1 = len(list(G.nodes))
    node_num_final = node_num_f1
    if node_num_f1 > 0:
        CCG = nx.connected_components(G)
        communities = list(CCG)
        comm_more_2 = [comm for comm in communities if len(comm)>2]
        all_nodes_2_lists = []
        comm_more_2_lists = [list(comm) for comm in comm_more_2]
        for comm in comm_more_2_lists:
            all_nodes_2_lists.extend(comm)
        remove2 = [node for node in list(G.nodes) if str(node) not in all_nodes_2_lists]
        if remove2:
            G.remove_nodes_from(remove2)
        node_num_final = len(list(G.nodes))
        if node_num_final >= num_of_cols*max_col_perc:
            if verbose:
                print("filtered by jaccard 0.8")
            G.remove_edges_from([(n1, n2) for n1, n2, w in G.edges(data="weight") if w < 0.8])
            node_num_final = len(list(G.nodes))
            if node_num_final >= num_of_cols*max_col_perc:
                if verbose:
                    print("filtered by jaccard 0.9")
                G.remove_edges_from([(n1, n2) for n1, n2, w in G.edges(data="weight") if w < 0.9])
                node_num_final = len(list(G.nodes))
            if node_num_final > 0:
                CCG = nx.connected_components(G)
                communities = list(CCG)
                comm_more_2_lists = [list(comm) for comm in communities]
            if node_num_final >= num_of_cols*max_col_perc:
                if verbose:
                    print("filtered by degree 3")
                remove = [node for node,degree in dict(G.degree()).items() if degree < 3]
                if node_num_final - len(remove) >= 30:
                    if verbose:
                        print("filtering by degree 3 yielded 30+ nodes")
                    G.remove_nodes_from(remove)
                    node_num_final = len(list(G.nodes))
                    CCG = nx.connected_components(G)
                    communities = list(CCG)
                    comm_more_2_lists = [list(comm) for comm in communities]
    if verbose:
        print(f"N cols:{num_of_cols}, N nodes: {node_num}, N nodes f: {node_num_final}")
    return comm_more_2_lists, G

def get_communities(msa_def_name, correlation_dictionary, num_of_cols, comm_type="filter", verbose=False):
    communities_function = generate_communities
    if comm_type == "simple":
        communities_function = generate_communities_simple
    communities_data, G = communities_function(correlation_dictionary['list'], num_of_cols, verbose=verbose)
    comm_parsed = (communities_data, G)
    communities_path = f"./CONSCOEVO_JOBS/{msa_def_name}/comm/"
    make_output_dir(communities_path)
    pickle_writer(f"{communities_path}{msa_def_name}_comms_{comm_type}.pickle", comm_parsed)
    return comm_parsed

def parse_comms(header_respos_list, communities_data, any_pos="verystrict"):
    all_comm_pos = []
    if any_pos == "verystrict":
        for comm in communities_data:
            positions = [int(respos[1:]) for respos in comm]
            if set(comm).issubset(header_respos_list):
                all_comm_pos.extend(positions)
    elif any_pos == "strict":
        for comm in communities_data:
            in_comm = [respos for respos in comm if respos in header_respos_list]
            positions = [int(respos[1:]) for respos in in_comm]
            all_comm_pos.extend(positions)
    return all_comm_pos

def main_run_cons_coevo(arguments, verbose=True):
    #MSAfile
    if len(arguments) != 3:
        return([0, "Missing arguments", None])
    msafile = arguments[1]
    # header_list, sequence_list = readFastaMSAlists(aln_type_alnfile, clean=True)
    # header_name = uniprot_to_header_dict[uniprot_code]
    # index_header = header_list.index(header_name)
    header_list, sequence_list = readFastaMSA(msafile, mode="list")
    header_name = header_list[0]
    if verbose:
        print("header_name")
        print(header_name)
    header_sequence = sequence_list[0]
    header_respos_list = [pos+str(i+1) for i,pos in enumerate(header_sequence)]

    num_of_rows = len(header_list)
    num_of_cols = len(sequence_list[0])
    if verbose:
        print("num_of_rows")
        print(num_of_rows)
        print("num_of_cols")
        print(num_of_cols)
    henikoff_weights_list = calculate_henikoff_weights_lists(header_name, sequence_list, num_of_rows, num_of_cols, verbose=verbose)

    if verbose:
        print("henikoff_weights_list")
        print(henikoff_weights_list)

    freq_calc_type = "henikoff_blosum62"
    conservation_dictionary = get_conservations_lists(header_name, sequence_list, num_of_cols, seq_weights_list=henikoff_weights_list, pseudocount=0.0, distribution="blosum62", jensen_weight=0.5, verbose=verbose)

    if verbose:
        print('conservation_dictionary["sh-vm"]')
        print(conservation_dictionary["sh-vm"])

    gapstats_per_column = calc_gap_stats_lists(header_name, sequence_list, num_of_cols, freq_calc_type, henikoff_weights_list, verbose=verbose)
    msa_columns = []
    for icol in range(num_of_cols):
        msa_columns.append([seq[icol] for seq in sequence_list])
    if verbose:
        print('gapstats_per_column')
        print(gapstats_per_column)
    min_freq = 0.2
    similarity_filter = 0.7
    pvalue_filter = 20
    similarity_type = "jaccard"
    global binom_dict, hypergeom_dict
    binom_dict = {}
    hypergeom_dict = {}
    correlation_dictionary = get_filtered_correlations_lists(header_name, freq_calc_type, msa_columns, min_freq, num_of_rows, num_of_cols, henikoff_weights_list, similarity_filter, pvalue_filter, from_respos=[], to_respos=[], similarity_type=similarity_type, pvalue_type="hypergeom", weighted_calc=True, verbose=verbose)
    if verbose:
        print(len(correlation_dictionary['list']))
    communities_data = get_communities(header_name, correlation_dictionary, num_of_cols, comm_type="filter", verbose=verbose)
    communities_data = communities_data[0]
    if verbose:
        print("communities_data")
        print(communities_data)
    all_comm_pos = parse_comms(header_respos_list, communities_data)
    back_comm_pos = parse_comms(header_respos_list, communities_data, any_pos="strict")
    if len(all_comm_pos) == 0:
        all_comm_pos = back_comm_pos.copy()
    if verbose:
        print("all_comm_pos")
        print(all_comm_pos)
    #sent to results and save export
    results = [conservation_dictionary["sh-vm"], all_comm_pos, back_comm_pos]
    return(1, "Cons. and Coevo calculated", results)

if __name__ == "__main__":
    arguments = sys.argv
    status, status_string, results = main_run_cons_coevo(arguments)
    if results[0] == 0:
        raise Exception(results[1])
    pickle_writer(f"{arguments[2]}_conscoevo.pickle", results)