#!/usr/bin/env python

####################################
### Structural Features          ###
### Written by: Nicole Zatorski  ###
### Last date modified: 22/7/22  ###
####################################

# import statements
from itertools import filterfalse
import sys
import glob
import sqlite3
import math
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import fdrcorrection
from statsmodels.stats.weightstats import DescrStatsW
from numpy import std
import scipy.stats

# functions
def check_if_thresholds_met(sample_info,prob,evalue,pvalue,coverage,percent_identity,len_template):
    '''
    For Scope and Interproscan classification, checks if the class aplies based on the thresholds
    Inputs:
        sample_info (list): values for the sample to compare against the threshold
        prob (int): between 0-100
        evalue (float): e value
        pvalue (float): p value
        coverage (float): sequence coverage
        percent_identity (int): between 0-100
        len_template (int): length in amino acid of sequence to compare
    Outputs:
        True or False (bool)
    '''
    try:
        sample_prob = sample_info[10]
        if sample_prob < prob:
            return False
        sample_evalue= sample_info[11]
        if sample_evalue > evalue:
            return False
        sample_pvalue= sample_info[12]
        if sample_pvalue > pvalue:
            return False
        sample_coverage= sample_info[13]
        if sample_coverage < coverage:
            return False
        sample_percent_identity= sample_info[16]
        if sample_percent_identity < percent_identity:
            return False
        sample_len_template= sample_info[18]-sample_info[17]
        if sample_len_template < len_template:
            return False
        return True
    except:
        return False

def add_to_freq_dict(freq_dict, key_to_add):
    '''
    Adds one to the frequency count of a dictionary key
    Inputs:
        freq_dict (dict): dictionary to add to
        key_to_add (str): name of key
    Outputs:
        freq_dict (dict): dictionary that was added to
    '''
    if key_to_add not in freq_dict:
        freq_dict[key_to_add] = 1
    else:
        freq_dict[key_to_add] = freq_dict[key_to_add] + 1
    return freq_dict

def get_substructs_from_oneid(id, path_to_db,prob,evalue,pvalue,coverage,percent_identity,len_template, out_dict):
    '''
    Gets scope and interproscan classes from a protien id based on the input cutoff criteria
    Inputs:
        id (str): protien name
        path_to_db (str): path to the database contianing the scope and interproscan information 
        prob (int): between 0-100 cutoff value
        evalue (float): e value cutoff value
        pvalue (float): p value cutoff value
        coverage (float): sequence coverage cutoff value
        percent_identity (int): between 0-100 cutoff value
        len_template (int): length in amino acid of sequence to compare
        out_dict (dict): collecting variable for sub structures of all proteins 
    Outputs:
        out_dict (dict): collecting variable for sub structures of all proteins 
    '''
    con = sqlite3.connect(path_to_db)
    cursor = con.cursor()
    g_uniprot_info = cursor.execute('SELECT * from uniprot_info WHERE gname like ?', ('%' + id + '%',))
    r_uid = g_uniprot_info.fetchall()
    if r_uid != []:
        uid = r_uid[0][0]
    else:
        uid = id
    domain_info = cursor.execute("SELECT * from domain WHERE uid like ?", ('%' + uid + '%',))
    r_domain = domain_info.fetchall()
    for list_index in range(len(r_domain)):
        ips_list = r_domain[list_index][5].split(';')
        descript_list = r_domain[list_index][7].split(';')
        for i in range(len(ips_list)):
            key = (ips_list[i],descript_list[i])
            out_dict['domain'] = add_to_freq_dict(out_dict['domain'], key)
    fold_info = cursor.execute("SELECT * from fold WHERE uid like ?", ('%' + uid + '%',))
    r_fold = fold_info.fetchall()
    for i in r_fold:
        if check_if_thresholds_met(i,prob,evalue,pvalue,coverage,percent_identity,len_template):
            out_dict['fold'] = add_to_freq_dict(out_dict['fold'], (i[4], i[5]))
            out_dict['superfamily'] = add_to_freq_dict(out_dict['superfamily'], (i[6], i[7]))
            out_dict['family'] = add_to_freq_dict(out_dict['family'], (i[8], i[9]))
    con.close()
    return out_dict


def init_track_dict(header_file):
    '''
    Creates a dictionary with headers as keys and 0 as values
    Inputs:
        header_file (list): list of keys for the dictionary 
    Outputs:
        track_dict (dict): dictionary with headers as keys and 0 as values
    '''
    track_dict = {}
    for key in header_file:
        track_dict[key] = 0
    return track_dict

def yn_update(out, normout, track_dict, struct_type):
    '''
    Indicates if a feature is present in any amount
    Inputs:
        out (str): string representing the collected output
        normout (str): string representing the normalized collected output
        track_dict (dict): dictionary with feature frequency information
        struct_type (str): corresponds to a header in track_dict, the feature being investigated
    Outputs:
        out (str): string representing the collected output
        normout (str): string representing the normalized collected output
    '''
    if track_dict['Number of ' + struct_type] >= 1:
        out = out + str(1)+ ','
        normout = normout + str(1)+ ','
    else:
        out = out + str(0)+ ','
        normout = normout + str(0)+ ','
    return out, normout

def check_if_found(database_dir, found_dict, gnuid, weight, unfound, found, tot_weight):
    '''
    Indicates if a protien or gene name is found in the database
    Inputs:
        database_dir (str): directory for the structural features database
        found_dict (dict): collecting variable of found gene names and protein ids with their associated weights (according to their expression levels)
        gnuid (str): id of gene or protien
        weight (int): expression level of the gene or protein
        unfound (lis): list of unfound genes or proteins 
        found (int): number of found genes or proteins
        tot_weight (int): sum of all expresion magnitudes for entire sample
    Outputs:
        tot_weight (int): sum of all expresion magnitudes for entire sample
        found_dict (dict): collecting variable of found gene names and protein ids with their associated weights (according to their expression levels)
        unfound (lis): list of unfound genes or proteins 
        found (int): number of found genes or proteins
    '''
    try:
        with open(database_dir+ gnuid + '.txt') as prerun_file:
            if gnuid not in found_dict:
                found +=1
                found_dict[gnuid] = float(weight)
                tot_weight+= float(weight)
            else:
                found_dict[gnuid] = found_dict[gnuid] + float(weight)
    except:
        unfound.append(gnuid)
    return tot_weight, found_dict, unfound, found


def num_find_sig_for_one_file(sample_file, database_dir, header, unfound, found, use_weight):
    '''
    Indicates if a protien or gene name is found in the database
    Inputs:
        sample_file (str): name of file containing list of line separaged gene or protien names with or without weights
        database_dir (str): directory for the structural features database
        unfound (lis): list of unfound genes or proteins 
        found (int): number of found genes or proteins
        use_weight (bool): should feature counts be weighted by expression levels?
    Outputs:
        sd_out (dict): standard deviations of continous feature variables
        track_dict (dict): features and their frequencies
        out_dict (dict): scope and interproscan feature frequencies
        unfound (lis): list of unfound genes or proteins 
        found (int): number of found genes or proteins
    '''
    track_dict = init_track_dict(header)
    prob= 50
    evalue= 1e-5
    pvalue= 1e-5
    coverage= 0.3
    percent_identity= 30
    len_template= 30
    out_dict = {
        'domain':{},
        'fold':{},
        'superfamily':{},
        'family':{},}
    path_to_db= 'databases/structure_database.db'
    tot_weight = 0
    found_dict = {}
    sd_dict = {'Length of protein': [], 
                'Negative region lengths' : [],
                'Positive region lengths' : [],
                'Total length of anchor regions' : [],
                'Total length of coil regions' : [],
                'Total length of conserved regions' : [],
                'Total length of disordered regions' : [],
                'Total length of globular regions' : [],
                'Total length of helix regions' : [],
                'Total length of loop regions' : [],
                'Total length of nonconserved regions' : [],
                'Total length of sheet regions' : [],
                'Total length tmh regions' : []}
    sd_index = {2: 'Length of protein', 4: 'Negative region lengths', 220: 'Positive region lengths', 222: 'Total length of anchor regions', 223: 'Total length of coil regions', 224: 'Total length of conserved regions', 225: 'Total length of disordered regions', 226: 'Total length of globular regions', 227: 'Total length of helix regions', 228: 'Total length of loop regions', 229: 'Total length of nonconserved regions', 230: 'Total length of sheet regions', 231: 'Total length tmh regions'}
    
    with open(sample_file) as sample:
        for line in sample:
            split_line = line[:-1].split(',')
            gnuid = split_line[0]
            if use_weight:
                if len(split_line) == 2:
                    weight = split_line[1]
                else:
                    weight = 1
            else:
                weight = 1
            tot_weight, found_dict, unfound, found = check_if_found(database_dir, found_dict, gnuid, weight, unfound, found, tot_weight)
            # if found >=250:
            #     break
    
    weight_list = []
    for gnuid in found_dict:
        weight_list.append(found_dict[gnuid])
        out_dict = get_substructs_from_oneid(gnuid, path_to_db,prob,evalue,pvalue,coverage,percent_identity,len_template,out_dict)
        with open(database_dir+ gnuid + '.txt') as prerun_file:
            for line in prerun_file:
                split_line = line.split(',')
                for index in range(1, len(split_line)):
                    if index in sd_index:
                        sd_dict[sd_index[index]] = sd_dict[sd_index[index]] + [float(split_line[index])]
                    track_dict[header[index]] = track_dict[header[index]] + (float(split_line[index])*found_dict[gnuid]/tot_weight)
    sd_out = {}
    for sd in sd_dict:
        if use_weight:
            sd_out[sd]=DescrStatsW(sd_dict[sd], weights=weight_list, ddof=1).std
        else:
            sd_out[sd] = std(sd_dict[sd])
    return sd_out,track_dict, out_dict,found, unfound

def write_output(out_name, out_str):
    '''
    Write output file from string
    Inputs:
        out_name (str): file to write to
        out_str (str): content to write to file
    Outputs:
        None
    '''
    f = open(out_name, 'a+')
    f.write(out_str)
    f.close()

def make_background_dict(path_to_background, background_file, background_dict):
    '''
    Make background dictionary from background frequency file
    Inputs:
        path_to_background (str): path to the file with the background frequency
        background_file (str): name of the file iwth the background frequency
        background_dict (dict): collecting variable for background frequencies
    Outputs:
        background_dict (dict): collecting variable for background frequencies
    '''
    with open(path_to_background+background_file) as fo:
        for line in fo:
            split_line = line[:-1].split(',')
            if len(split_line) == 2:
                background_dict[split_line[0]] = float(split_line[1])
    return background_dict

def make_average_background_dict(path_to_background, background_file, background_dict):
    '''
    Make background dictionary from background averages file
    Inputs:
        path_to_background (str): path to the file with the background averages
        background_file (str): name of the file iwth the background averages
        background_dict (dict): collecting variable for background averages
    Outputs:
        background_dict (dict): collecting variable for background averages
    '''
    with open(path_to_background+background_file) as fo:
        for line in fo:
            split_line = line[:-1].split(',')
            background_dict[split_line[0]] = (float(split_line[1]), float(split_line[2]))
    return background_dict

def compare_frequency_to_background(len_out_dict_sub_dict, feature_name, sample_sf_frequency, proteins_found_in_sample, background_dict, proteins_in_proteome):
    '''
    Statistical tests for frequency variables
    Inputs:
        len_out_dict_sub_dict (int): number of features found (used for p value correction)
        feature_name (str): feature being statistically compared between sample and background
        sample_sf_frequency (int): number of times the feature appears in the experimental signature
        proteins_found_in_sample (int): number of proteins found in the sample
        background_dict (dict): collecting variable for background averages
    Outputs:
        sample_sf_frequency (int): number of times the feature appears in the experimental signature
        background_sf_frequency (int): number of times the feature appears in the background signature
        p_value (float): from chi squared test between the sample and background
        corrected_p (float): bonferroni corrected p value
        fc (float): logfold change of sample/background
    '''
    sample_wo = proteins_found_in_sample - sample_sf_frequency
    if feature_name in background_dict:
        background_sf_frequency = background_dict[feature_name]
    else:
        background_sf_frequency = 0
    background_wo = proteins_in_proteome - background_sf_frequency
    oddsratio,p_value = fisher_exact([[abs(sample_sf_frequency),abs(sample_wo)],[abs(background_sf_frequency),abs(background_wo)]], alternative='greater')
    corrected_p = 0.05/len_out_dict_sub_dict
    if proteins_found_in_sample != 0:
        sample_freq_percent = sample_sf_frequency/proteins_found_in_sample 
    else:
        sample_freq_percent = 0
    background_freq_percent = background_sf_frequency/proteins_in_proteome
    if background_freq_percent > 0 and sample_freq_percent > 0:
        fc = math.log(sample_freq_percent/background_freq_percent, 10)
    else:
        if background_freq_percent == 0 or sample_freq_percent == 0:
            fc = 'Divide by zero error: ' + feature_name + ' not in background'
        else:
            fc = 'Error'
    return sample_sf_frequency, background_sf_frequency, p_value, corrected_p, fc


def compare_quant_to_background(len_out_dict_sub_dict, feature_name, sample_sf_frequency, proteins_found_in_sample, background_dict, proteins_in_proteome):
    '''
    Statistical tests for frequency variables
    Inputs:
        len_out_dict_sub_dict (int): number of features found (used for p value correction)
        feature_name (str): feature being statistically compared between sample and background
        sample_sf_frequency (tuple): value of feature in the experimental signature, standard deviation of feature values
        proteins_found_in_sample (int): number of proteins found in the sample
        background_dict (dict): collecting variable for background averages
    Outputs:
        sample_sf_frequency (tuple): value of feature in the experimental signature, standard deviation of feature values
        background_sf_frequency (tuple): value of feature in the background signature, standard deviation of feature values
        p_value (float): from t-test between the sample and background
        corrected_p (float): bonferroni corrected p value
    '''
    if feature_name in background_dict:
        background_sf_frequency = background_dict[feature_name]
    else:
        background_sf_frequency = (0,0)
    t_value = (background_sf_frequency[0]-sample_sf_frequency[0])/((((background_sf_frequency[1]**2)/proteins_in_proteome) +((sample_sf_frequency[1]**2)/proteins_found_in_sample))**0.5)
    df = proteins_found_in_sample-2
    p_value = scipy.stats.t.sf(abs(t_value), df=df)
    corrected_p = 0.05/ len_out_dict_sub_dict
    return sample_sf_frequency, background_sf_frequency, p_value, corrected_p




def find_alpha_features_for_one(sample_file, header, use_weight):
    sd_dict = {
        'Length of S Regions':[], 
        'Length of E Regions':[],
        'Length of T Regions':[],
        'Length of B Regions':[],
        'Length of G Regions':[],
        'Length of H Regions':[],
        'Length of Aggregation Prone Regions':[],
        'Minimum Distance to Center of Mass':[],
        'Maximum Distance to Center of Mass':[],
        'Average Distance to Center of Mass':[]
    }
    sd_index = {
        0: 'Length of S Regions', 
        22: 'Length of E Regions',
        44: 'Length of T Regions',
        66: 'Length of B Regions',
        88: 'Length of G Regions',
        110: 'Length of H Regions',
        132: 'Length of Aggregation Prone Regions',
        154: 'Minimum Distance to Center of Mass',
        155: 'Maximum Distance to Center of Mass',
        156: 'Average Distance to Center of Mass'
    }
    
    sample_id_dict = {}
    multiples = {}
    out_dict = {}
    for index in range(len(header)):
        out_dict[header[index]] = 0
    with open(sample_file) as fo:
        counter = set([])
        for line in fo:
            split_line = line[:-1].split(',')
            gnuid = split_line[0]
            try: 
                f = open('databases/precounted_alpha_fold/' + gnuid + '.txt', 'r')
                f.close()
                counter.add(gnuid)
                if use_weight:
                    if len(split_line) == 2:
                        weight = split_line[1]
                    else:
                        weight = 1
                else:
                    weight = 1
                if gnuid not in sample_id_dict:
                    sample_id_dict[gnuid] = float(weight)
                elif gnuid not in multiples:
                    multiples[gnuid] = [sample_id_dict[gnuid], float(weight)]
                else:
                    multiples[gnuid] = multiples[gnuid] + [float(weight)]
            except:
                pass
    for elt in multiples:
        sample_id_dict[elt] = sum(multiples[elt])/len(multiples[elt])
    tot_weight = sum(sample_id_dict.values())
    print(len(counter))
    weight_list = []
    for gnuid in sample_id_dict:
        weight_list.append(sample_id_dict[gnuid])
        f = open('databases/precounted_alpha_fold/' + gnuid + '.txt', 'r')
        vals = f.read()[:-1].split(',')[1:]
        for index in range(len(vals)):
            head = header[index]
            out_dict[head] = out_dict[head] + (float(vals[index])*sample_id_dict[gnuid]/tot_weight)
            if index in sd_index:
                    sd_dict[sd_index[index]] = sd_dict[sd_index[index]] + [float(vals[index])]
    sd_out = {}
    for sd in sd_dict:
        if use_weight:
            sd_out[sd]=DescrStatsW(sd_dict[sd], weights=weight_list, ddof=1).std
        else:
            sd_out[sd] = std(sd_dict[sd])
    return sd_out, out_dict

def run_for_all_files_in_folder(input_dir, folder_out, use_weight=False, background_folder_name='human_background'):
    '''
    Generates structural features for all files in a directory
    Inputs:
        input_dir (str): directory containing all input files
        folder_out (str): directory containing all output files
        use_weight (bool): default False, weigh structural feature output by corresponding input expression levels
        background_folder_name (str): default 'human_backgrounds', file name of background to use
    Outputs:
        (float) percentage of gene names and proteins found in the structural features database
    '''
    database_dir = './databases/precounted_human_genome/'
    header = ['Name','Crowd predictions', 'Length of protein', 'NHTM Best from query.phdPred', 'Negative region lengths', 'Number amino acid in anchor region A', 'Number amino acid in anchor region C', 'Number amino acid in anchor region D', 'Number amino acid in anchor region E', 'Number amino acid in anchor region F', 'Number amino acid in anchor region G', 'Number amino acid in anchor region H', 'Number amino acid in anchor region I', 'Number amino acid in anchor region K', 'Number amino acid in anchor region L', 'Number amino acid in anchor region M', 'Number amino acid in anchor region N', 'Number amino acid in anchor region P', 'Number amino acid in anchor region Q', 'Number amino acid in anchor region R', 'Number amino acid in anchor region S', 'Number amino acid in anchor region T', 'Number amino acid in anchor region V', 'Number amino acid in anchor region W', 'Number amino acid in anchor region Y', 'Number amino acid in coil region A', 'Number amino acid in coil region C', 'Number amino acid in coil region D', 'Number amino acid in coil region E', 'Number amino acid in coil region F', 'Number amino acid in coil region G', 'Number amino acid in coil region H', 'Number amino acid in coil region I', 'Number amino acid in coil region K', 'Number amino acid in coil region L', 'Number amino acid in coil region M', 'Number amino acid in coil region N', 'Number amino acid in coil region P', 'Number amino acid in coil region Q', 'Number amino acid in coil region R', 'Number amino acid in coil region S', 'Number amino acid in coil region T', 'Number amino acid in coil region V', 'Number amino acid in coil region W', 'Number amino acid in coil region Y', 'Number amino acid in conserved region A', 'Number amino acid in conserved region C', 'Number amino acid in conserved region D', 'Number amino acid in conserved region E', 'Number amino acid in conserved region F', 'Number amino acid in conserved region G', 'Number amino acid in conserved region H', 'Number amino acid in conserved region I', 'Number amino acid in conserved region K', 'Number amino acid in conserved region L', 'Number amino acid in conserved region M', 'Number amino acid in conserved region N', 'Number amino acid in conserved region P', 'Number amino acid in conserved region Q', 'Number amino acid in conserved region R', 'Number amino acid in conserved region S', 'Number amino acid in conserved region T', 'Number amino acid in conserved region V', 'Number amino acid in conserved region W', 'Number amino acid in conserved region Y', 'Number amino acid in disordered region A', 'Number amino acid in disordered region C', 'Number amino acid in disordered region D', 'Number amino acid in disordered region E', 'Number amino acid in disordered region F', 'Number amino acid in disordered region G', 'Number amino acid in disordered region H', 'Number amino acid in disordered region I', 'Number amino acid in disordered region K', 'Number amino acid in disordered region L', 'Number amino acid in disordered region M', 'Number amino acid in disordered region N', 'Number amino acid in disordered region P', 'Number amino acid in disordered region Q', 'Number amino acid in disordered region R', 'Number amino acid in disordered region S', 'Number amino acid in disordered region T', 'Number amino acid in disordered region V', 'Number amino acid in disordered region W', 'Number amino acid in disordered region Y', 'Number amino acid in globular region A', 'Number amino acid in globular region C', 'Number amino acid in globular region D', 'Number amino acid in globular region E', 'Number amino acid in globular region F', 'Number amino acid in globular region G', 'Number amino acid in globular region H', 'Number amino acid in globular region I', 'Number amino acid in globular region K', 'Number amino acid in globular region L', 'Number amino acid in globular region M', 'Number amino acid in globular region N', 'Number amino acid in globular region P', 'Number amino acid in globular region Q', 'Number amino acid in globular region R', 'Number amino acid in globular region S', 'Number amino acid in globular region T', 'Number amino acid in globular region V', 'Number amino acid in globular region W', 'Number amino acid in globular region Y', 'Number amino acid in helix region A', 'Number amino acid in helix region C', 'Number amino acid in helix region D', 'Number amino acid in helix region E', 'Number amino acid in helix region F', 'Number amino acid in helix region G', 'Number amino acid in helix region H', 'Number amino acid in helix region I', 'Number amino acid in helix region K', 'Number amino acid in helix region L', 'Number amino acid in helix region M', 'Number amino acid in helix region N', 'Number amino acid in helix region P', 'Number amino acid in helix region Q', 'Number amino acid in helix region R', 'Number amino acid in helix region S', 'Number amino acid in helix region T', 'Number amino acid in helix region V', 'Number amino acid in helix region W', 'Number amino acid in helix region Y', 'Number amino acid in loop region A', 'Number amino acid in loop region C', 'Number amino acid in loop region D', 'Number amino acid in loop region E', 'Number amino acid in loop region F', 'Number amino acid in loop region G', 'Number amino acid in loop region H', 'Number amino acid in loop region I', 'Number amino acid in loop region K', 'Number amino acid in loop region L', 'Number amino acid in loop region M', 'Number amino acid in loop region N', 'Number amino acid in loop region P', 'Number amino acid in loop region Q', 'Number amino acid in loop region R', 'Number amino acid in loop region S', 'Number amino acid in loop region T', 'Number amino acid in loop region V', 'Number amino acid in loop region W', 'Number amino acid in loop region Y', 'Number amino acid in nonconserved region A', 'Number amino acid in nonconserved region C', 'Number amino acid in nonconserved region D', 'Number amino acid in nonconserved region E', 'Number amino acid in nonconserved region F', 'Number amino acid in nonconserved region G', 'Number amino acid in nonconserved region H', 'Number amino acid in nonconserved region I', 'Number amino acid in nonconserved region K', 'Number amino acid in nonconserved region L', 'Number amino acid in nonconserved region M', 'Number amino acid in nonconserved region N', 'Number amino acid in nonconserved region P', 'Number amino acid in nonconserved region Q', 'Number amino acid in nonconserved region R', 'Number amino acid in nonconserved region S', 'Number amino acid in nonconserved region T', 'Number amino acid in nonconserved region V', 'Number amino acid in nonconserved region W', 'Number amino acid in nonconserved region Y', 'Number amino acid in protein A', 'Number amino acid in protein C', 'Number amino acid in protein D', 'Number amino acid in protein E', 'Number amino acid in protein F', 'Number amino acid in protein G', 'Number amino acid in protein H', 'Number amino acid in protein I', 'Number amino acid in protein K', 'Number amino acid in protein L', 'Number amino acid in protein M', 'Number amino acid in protein N', 'Number amino acid in protein P', 'Number amino acid in protein Q', 'Number amino acid in protein R', 'Number amino acid in protein S', 'Number amino acid in protein T', 'Number amino acid in protein V', 'Number amino acid in protein W', 'Number amino acid in protein Y', 'Number amino acid in sheet region A', 'Number amino acid in sheet region C', 'Number amino acid in sheet region D', 'Number amino acid in sheet region E', 'Number amino acid in sheet region F', 'Number amino acid in sheet region G', 'Number amino acid in sheet region H', 'Number amino acid in sheet region I', 'Number amino acid in sheet region K', 'Number amino acid in sheet region L', 'Number amino acid in sheet region M', 'Number amino acid in sheet region N', 'Number amino acid in sheet region P', 'Number amino acid in sheet region Q', 'Number amino acid in sheet region R', 'Number amino acid in sheet region S', 'Number amino acid in sheet region T', 'Number amino acid in sheet region V', 'Number amino acid in sheet region W', 'Number amino acid in sheet region Y', 'Number of anchor regions', 'Number of coils', 'Number of conserved regions', 'Number of disordered regions', 'Number of globular regions', 'Number of helix', 'Number of loops', 'Number of negative regions', 'Number of negative regions with length >=30', 'Number of nonconserved regions', 'Number of positive regions', 'Number of positive regions with length >=30', 'Number of predictions', 'Number of sheets', 'Number of transmembrane helices', 'Positive region lengths', 'Stretch', 'Total length of anchor regions', 'Total length of coil regions', 'Total length of conserved regions', 'Total length of disordered regions', 'Total length of globular regions', 'Total length of helix regions', 'Total length of loop regions', 'Total length of nonconserved regions', 'Total length of sheet regions', 'Total length tmh regions', 'Y/n anchor regions', 'Y/n disordered regions', 'Y/n globular regions', 'Y/n tmh regions']
    headeralpha = ['Length of S Regions', 'Number of S Regions', 'Number of Amino Acids A in S Regions', 'Number of Amino Acids R in S Regions', 'Number of Amino Acids N in S Regions', 'Number of Amino Acids D in S Regions', 'Number of Amino Acids C in S Regions', 'Number of Amino Acids E in S Regions', 'Number of Amino Acids Q in S Regions', 'Number of Amino Acids G in S Regions', 'Number of Amino Acids H in S Regions', 'Number of Amino Acids I in S Regions', 'Number of Amino Acids L in S Regions', 'Number of Amino Acids K in S Regions', 'Number of Amino Acids M in S Regions', 'Number of Amino Acids F in S Regions', 'Number of Amino Acids P in S Regions', 'Number of Amino Acids S in S Regions', 'Number of Amino Acids T in S Regions', 'Number of Amino Acids W in S Regions', 'Number of Amino Acids Y in S Regions', 'Number of Amino Acids V in S Regions', 'Length of E Regions', 'Number of E Regions', 'Number of Amino Acids A in E Regions', 'Number of Amino Acids R in E Regions', 'Number of Amino Acids N in E Regions', 'Number of Amino Acids D in E Regions', 'Number of Amino Acids C in E Regions', 'Number of Amino Acids E in E Regions', 'Number of Amino Acids Q in E Regions', 'Number of Amino Acids G in E Regions', 'Number of Amino Acids H in E Regions', 'Number of Amino Acids I in E Regions', 'Number of Amino Acids L in E Regions', 'Number of Amino Acids K in E Regions', 'Number of Amino Acids M in E Regions', 'Number of Amino Acids F in E Regions', 'Number of Amino Acids P in E Regions', 'Number of Amino Acids S in E Regions', 'Number of Amino Acids T in E Regions', 'Number of Amino Acids W in E Regions', 'Number of Amino Acids Y in E Regions', 'Number of Amino Acids V in E Regions', 'Length of T Regions', 'Number of T Regions', 'Number of Amino Acids A in T Regions', 'Number of Amino Acids R in T Regions', 'Number of Amino Acids N in T Regions', 'Number of Amino Acids D in T Regions', 'Number of Amino Acids C in T Regions', 'Number of Amino Acids E in T Regions', 'Number of Amino Acids Q in T Regions', 'Number of Amino Acids G in T Regions', 'Number of Amino Acids H in T Regions', 'Number of Amino Acids I in T Regions', 'Number of Amino Acids L in T Regions', 'Number of Amino Acids K in T Regions', 'Number of Amino Acids M in T Regions', 'Number of Amino Acids F in T Regions', 'Number of Amino Acids P in T Regions', 'Number of Amino Acids S in T Regions', 'Number of Amino Acids T in T Regions', 'Number of Amino Acids W in T Regions', 'Number of Amino Acids Y in T Regions', 'Number of Amino Acids V in T Regions', 'Length of B Regions', 'Number of B Regions', 'Number of Amino Acids A in B Regions', 'Number of Amino Acids R in B Regions', 'Number of Amino Acids N in B Regions', 'Number of Amino Acids D in B Regions', 'Number of Amino Acids C in B Regions', 'Number of Amino Acids E in B Regions', 'Number of Amino Acids Q in B Regions', 'Number of Amino Acids G in B Regions', 'Number of Amino Acids H in B Regions', 'Number of Amino Acids I in B Regions', 'Number of Amino Acids L in B Regions', 'Number of Amino Acids K in B Regions', 'Number of Amino Acids M in B Regions', 'Number of Amino Acids F in B Regions', 'Number of Amino Acids P in B Regions', 'Number of Amino Acids S in B Regions', 'Number of Amino Acids T in B Regions', 'Number of Amino Acids W in B Regions', 'Number of Amino Acids Y in B Regions', 'Number of Amino Acids V in B Regions', 'Length of G Regions', 'Number of G Regions', 'Number of Amino Acids A in G Regions', 'Number of Amino Acids R in G Regions', 'Number of Amino Acids N in G Regions', 'Number of Amino Acids D in G Regions', 'Number of Amino Acids C in G Regions', 'Number of Amino Acids E in G Regions', 'Number of Amino Acids Q in G Regions', 'Number of Amino Acids G in G Regions', 'Number of Amino Acids H in G Regions', 'Number of Amino Acids I in G Regions', 'Number of Amino Acids L in G Regions', 'Number of Amino Acids K in G Regions', 'Number of Amino Acids M in G Regions', 'Number of Amino Acids F in G Regions', 'Number of Amino Acids P in G Regions', 'Number of Amino Acids S in G Regions', 'Number of Amino Acids T in G Regions', 'Number of Amino Acids W in G Regions', 'Number of Amino Acids Y in G Regions', 'Number of Amino Acids V in G Regions', 'Length of H Regions', 'Number of H Regions', 'Number of Amino Acids A in H Regions', 'Number of Amino Acids R in H Regions', 'Number of Amino Acids N in H Regions', 'Number of Amino Acids D in H Regions', 'Number of Amino Acids C in H Regions', 'Number of Amino Acids E in H Regions', 'Number of Amino Acids Q in H Regions', 'Number of Amino Acids G in H Regions', 'Number of Amino Acids H in H Regions', 'Number of Amino Acids I in H Regions', 'Number of Amino Acids L in H Regions', 'Number of Amino Acids K in H Regions', 'Number of Amino Acids M in H Regions', 'Number of Amino Acids F in H Regions', 'Number of Amino Acids P in H Regions', 'Number of Amino Acids S in H Regions', 'Number of Amino Acids T in H Regions', 'Number of Amino Acids W in H Regions', 'Number of Amino Acids Y in H Regions', 'Number of Amino Acids V in H Regions', 'Length of Aggregation Prone Regions', 'Number of Aggregation Prone Regions', 'Number of Amino Acids A in Aggregation Prone Regions', 'Number of Amino Acids R in Aggregation Prone Regions', 'Number of Amino Acids N in Aggregation Prone Regions', 'Number of Amino Acids D in Aggregation Prone Regions', 'Number of Amino Acids C in Aggregation Prone Regions', 'Number of Amino Acids E in Aggregation Prone Regions', 'Number of Amino Acids Q in Aggregation Prone Regions', 'Number of Amino Acids G in Aggregation Prone Regions', 'Number of Amino Acids H in Aggregation Prone Regions', 'Number of Amino Acids I in Aggregation Prone Regions', 'Number of Amino Acids L in Aggregation Prone Regions', 'Number of Amino Acids K in Aggregation Prone Regions', 'Number of Amino Acids M in Aggregation Prone Regions', 'Number of Amino Acids F in Aggregation Prone Regions', 'Number of Amino Acids P in Aggregation Prone Regions', 'Number of Amino Acids S in Aggregation Prone Regions', 'Number of Amino Acids T in Aggregation Prone Regions', 'Number of Amino Acids W in Aggregation Prone Regions', 'Number of Amino Acids Y in Aggregation Prone Regions', 'Number of Amino Acids V in Aggregation Prone Regions', 'Minimum Distance to Center of Mass', 'Maximum Distance to Center of Mass', 'Average Distance to Center of Mass', 'Number of Contacts']
    all_headers = header[1:] + headeralpha
    unfound = []
    current_num_files = 0
    input_files = glob.glob(input_dir + '*')
    tot_num_files = len(input_files)

    for sample_file in input_files:

        path_to_background = './databases/'+background_folder_name+'/'
        with open(path_to_background+'number_proteins_found.csv') as fo:
            for line in fo:
                background_found = int(line)

        background_dict = make_background_dict(path_to_background, 'ipr.domain.csv', {})
        background_dict = make_background_dict(path_to_background, 'scop.family.csv', background_dict)
        background_dict = make_background_dict(path_to_background, 'scop.fold.csv', background_dict)
        background_dict = make_background_dict(path_to_background, 'scop.superfam.csv', background_dict)
        frequency_background = make_background_dict(path_to_background, 'frequency_background.csv', {})
        average_background = make_average_background_dict(path_to_background, 'average_background.csv', {})
        try:
            print(str(current_num_files/tot_num_files)+'% done')
            current_num_files +=1
            sd_out, track_dict, out_dict, found, unfound = num_find_sig_for_one_file(sample_file, database_dir, header, unfound, 0, use_weight)
            sd_alpha, out_dict_alpha = find_alpha_features_for_one(sample_file, headeralpha, use_weight)
            sd_out.update(sd_alpha)
            track_dict.update(out_dict_alpha)
            all_p_vals1 = []
            temp_lines1 = []
            all_p_vals2 = []
            temp_lines2 = []
            write_output( folder_out + 'frequency_' +sample_file.split('/')[-1],'structure_id,structure,structure_type,counts_observed,background_counts,pvalue,bonforroni_cutoff,log_fold_change,fdr\n')
            write_output( folder_out + 'average_' +sample_file.split('/')[-1],'label,average_observed,observed_standard_deviation,background_average,background_standard_deviation,pvalue,bonforroni_cutoff,fdr\n')
            for index in range(len(all_headers)):
                if all_headers[index].split(' ')[0] == "Number" or all_headers[index].split(' ')[0] == 'Y/n': 
                    sample_sf_frequency, background_sf_frequency, p_value, corrected_p, fc = compare_frequency_to_background(367, all_headers[index], track_dict[all_headers[index]], found, frequency_background, background_found)
                    temp_lines1.append('N/A,' + all_headers[index]+ ',N/A,' + str(sample_sf_frequency) + ',' +str(background_sf_frequency)+ ',' +str(p_value) + ',' +str(corrected_p)+ ',' +str(fc)+ ',')
                    all_p_vals1.append(p_value)
                elif all_headers[index] != 'Crowd predictions' and all_headers[index] != 'Stretch' and all_headers[index] != 'NHTM Best from query.phdPred':
                    sample_sf_frequency, background_sf_frequency, p_value, corrected_p = compare_quant_to_background(23, all_headers[index], (track_dict[all_headers[index]],sd_out[all_headers[index]]), found, average_background, background_found)
                    temp_lines2.append(all_headers[index] + ',' + str(track_dict[all_headers[index]])+ ',' +str(sd_out[all_headers[index]])+ ',' +str(background_sf_frequency[0])+ ',' +str(background_sf_frequency[1])+ ',' +str(p_value) + ',' +str(corrected_p) + ',')
                    all_p_vals2.append(p_value)
            rejected, fdr_list1 =fdrcorrection(all_p_vals1)
            rejected, fdr_list2 =fdrcorrection(all_p_vals2)
            for index in range(len(fdr_list1)):
                write_output(folder_out + 'frequency_' +sample_file.split('/')[-1], temp_lines1[index] + str(fdr_list1[index]) + '\n')
            for index in range(len(temp_lines2)):
                write_output(folder_out + 'average_' +sample_file.split('/')[-1], temp_lines2[index] + str(fdr_list2[index]) +'\n')
            for sub_dict in out_dict:
                all_p_vals = []
                temp_lines = []
                for elt in out_dict[sub_dict]:
                    if elt[0] != 'NULL':
                        key_temp = str(elt)[1:-1].replace(',', '-')
                        key_temp = key_temp.replace('-', ',',1)
                        key_temp = key_temp.replace("'", '')
                        sample_sf_frequency, background_sf_frequency, p_value, corrected_p, fc = compare_frequency_to_background(len(out_dict[sub_dict]), elt[0], out_dict[sub_dict][elt], found, background_dict, background_found)
                        temp_lines.append(key_temp + ',' + sub_dict+ ',' + str(sample_sf_frequency) + ',' +str(background_sf_frequency)+ ',' +str(p_value) + ',' +str(corrected_p)+ ',' +str(fc)+ ',')
                        all_p_vals.append(p_value)
                rejected, fdr_list =fdrcorrection(all_p_vals)
                for index in range(len(fdr_list)):   
                    write_output( folder_out + 'frequency_' +sample_file.split('/')[-1],temp_lines[index] + str(fdr_list[index]) + '\n')
        except:
            write_output(folder_out+'errors.csv', sample_file + '\n')
    return(found/tot_num_files)

if __name__ == '__main__':
    run_for_all_files_in_folder(*sys.argv[1:])