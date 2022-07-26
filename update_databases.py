#!/usr/bin/env python

####################################
### Update Structural Features   ###
### Written by: Nicole Zatorski  ###
### Last date modified: 22/7/22  ###
####################################

# import statements
import pandas as pd
import glob

# functions
def make_file_dict(directory, file_ending):
    '''
    Loads data from file into dictionary format
    Inputs:
        directory (str): name of directory where the files to load can be found
        file_ending (str): the type of files to be loaded
    Outputs:
        out (dict): dictionary of file contents
    '''
    file_dict = {}
    files = glob.glob(directory+'*')
    for f in files:
        fn = f.split('/')[-1]
        names = fn[:-len(file_ending)].split('--')
        for name in names:
            if name not in file_dict.keys():
                file_dict[name] = f
            else:
                if len(file_dict[name])> len(fn):
                    file_dict[name] = f
    return file_dict

def str_to_dict(to_transform, added_key_label):
    '''
    Takes a string representation of a dict from a file and makes it into an actual dict
    Inputs:
        to_transform (str): string to make into a dictionary
        add_key_label (str): additional information to add to the dictionary key label
    Outputs:
        out (dict): dictionary of input contents
    '''
    new_num = ''
    new_key = ''
    out = {}
    for elt in to_transform[1:-1]:
        if elt in 'ARNDCEQGHILKMFPSTWYV':
            new_key = elt
        elif elt in "0123456789":
            new_num = new_num + elt
        elif elt == ',':
            out[added_key_label + ' ' + new_key] = int(new_num)
            new_key = ''
            new_num = ''
    return(out)

def num_predict_pro_parse(file_path):
    '''
    Parses predict protein data for inclusion in the structural features database
    Inputs:
        file_path (str): path to the file with the contents of predict protein for a gene or protein of interest
    Outputs:
        collect (dict): dictionary of predict protein features for the input gene or protein
    '''
    line_num = -1
    collect = {}
    check_dict = {}
    with open(file_path) as fo:
        for line in fo:
            line_num += 1
            split_line = line[:-1].split(',,')
            if line_num == 0:
                collect['Number of transmembrane helices']= int(split_line[0])
                if split_line[1] == '-':
                    collect['NHTM Best from query.phdPred'] = 0
                else:
                    collect['NHTM Best from query.phdPred'] = int(split_line[1])
            elif line_num == 2:
                key = ['Stretch', 'Crowd predictions', 'Number of predictions', 'Number of positive regions', 'Number of positive regions with length >=30','Positive region lengths', 'Number of negative regions', 'Number of negative regions with length >=30','Negative region lengths']
                key_index = [0,1,2,3,5,6,7,9,10]
                for i in range(len(key)):
                    ki = key_index[i]
                    key_name = key[i]
                    if ki == 6 or ki == 10:
                        charged_list = split_line[ki][1:-1].split(', ')
                        charged_len = 0
                        for elt in charged_list:
                            charged_len += int(elt)
                        collect[key_name] = charged_len
                    else:
                        collect[key_name] = int(split_line[ki])
            elif line_num == 4:
                collect['Number of coils'] = split_line[0]
                temp_dict = str_to_dict(split_line[2], 'Number amino acid in coil region')
                collect['Total length of coil regions'] = max(int(split_line[1]), sum(temp_dict.values()))

                collect = {**collect, **temp_dict}
            elif line_num == 6:
                key = ['Number of helix','Total length of helix regions', 'Number amino acid in helix region', 'Number of sheets', 'Total length of sheet regions', 'Number amino acid in sheet region', 'Number of loops' , 'Total length of loop regions', 'Number amino acid in loop region']
                key_index = [0,1,3,4,5,7,8,9,11]
                key.reverse()
                key_index.reverse()
                for i in range(len(key)):
                    ki = key_index[i]
                    key_name = key[i]
                    if ki == 3 or ki == 7 or ki == 11:
                        temp_dict = str_to_dict(split_line[ki], key_name)
                        check_key = key_name[len('Number amino acid in '):-len(' region')]
                        check_dict[check_key] = sum(temp_dict.values())
                        collect = {**collect, **temp_dict}
                    elif ki ==1 or ki==5 or ki==9:
                        check_key = key_name[len('Total length of '): -len(' regions')]
                        collect[key_name] = max(int(split_line[ki]), check_dict[check_key])
                    else:
                        collect[key_name] = int(split_line[ki])
            elif line_num == 8:
                key = ['Number of conserved regions', 'Total length of conserved regions', 'Number amino acid in conserved region', 'Number of nonconserved regions', 'Total length of nonconserved regions', 'Number amino acid in nonconserved region']
                for i in reversed(range(len(key))):
                    key_name = key[i]
                    if i == 2 or i == 5:
                        temp_dict = str_to_dict(split_line[i], key_name)
                        check_key = key_name[len('Number amino acid in '):-len(' region')]
                        check_dict[check_key] = sum(temp_dict.values())
                        collect = {**collect, **temp_dict}
                    elif i ==1 or i==4:
                        check_key = key_name[len('Total length of '): -len(' regions')]
                        collect[key_name] = max(int(split_line[ki]), check_dict[check_key])
                    else:
                        collect[key_name] = int(split_line[i])
    return collect

def get_glob(glob_file_name):
    '''
    Parses glob data for inclusion in the structural features database
    Inputs:
        glob_file_name (str): path to the file with the contents of glob for a gene or protein of interest
    Outputs:
        number_globs (int): number of globular regions in the protein
        glob_len_list (list): list of the length of each globular region
        glob_comp_dict (dict): dictionary of amino acid frequencies in globular regions for the input gene or protein
    '''
    number_globs = 0
    glob_len_list = []
    seen_globs = False
    glob_flag = None
    glob_comp_dict = {'R':0, 'H':0, 'K':0, 'D':0, 'E':0, 'S':0, 'T':0, 'N':0, 'Q':0, 'C':0, 'U':0, 'G':0, 'P':0, 'A':0, 'V':0, 'I':0, 'L':0, 'M':0, 'F':0, 'Y':0, 'W':0}
    with open(glob_file_name) as fileobject:
        for line in fileobject:
            if 'globular domain ' in line:
                number_globs += 1
                index_glob, len_glob = line[25:].split('.')
                start, stop = len_glob.split('-')
                len_segment = int(stop)-int(start) +1
                glob_len_list.append(len_segment)
                seen_globs = True
            else:
                if seen_globs:
                    if '#' in line:
                        break
                    for char in line:
                        if char == char.upper():
                            try:
                                glob_comp_dict[char] = glob_comp_dict[char] + 1
                            except:
                                pass
                                # glob_comp_dict['X'] = glob_comp_dict['X'] + 1
    return number_globs, glob_len_list, glob_comp_dict

def get_disorder(iupred_file_name ,disorder_cutoff):
    '''
    Parses iupred data for inclusion in the structural features database
    Inputs:
        iupred_file_name (str): path to the file with the contents of iupred for a gene or protein of interest
        disorder_cutoff (float): precent prediction of being an amino acid in a disordered region
    Outputs:
        disorder_list (list): list of the length of each disordered region
        tot_num_disorder (int): number of disordered regions in the protein
        disorder_aa_comp (dict): dictionary of amino acid frequencies in disordered regions for the input gene or protein
        anchor_list (list): list of the length of each anchor region
        tot_num_anchor (int): number of anchor regions in the protein
        anchor_aa_comp (dict): dictionary of amino acid frequencies in anchor regions for the input gene or protein
    '''
    disorder_flag = None
    len_disorder = 0
    disorder_list =[]
    tot_num_disorder = 0
    disorder_aa_comp = {'R':0, 'H':0, 'K':0, 'D':0, 'E':0, 'S':0, 'T':0, 'N':0, 'Q':0, 'C':0, 'U':0, 'G':0, 'P':0, 'A':0, 'V':0, 'I':0, 'L':0, 'M':0, 'F':0, 'Y':0, 'W':0}
    
    anchor_flag = None
    len_anchor = 0
    anchor_list = []
    tot_num_anchor = 0
    anchor_aa_comp = {'R':0, 'H':0, 'K':0, 'D':0, 'E':0, 'S':0, 'T':0, 'N':0, 'Q':0, 'C':0, 'U':0, 'G':0, 'P':0, 'A':0, 'V':0, 'I':0, 'L':0, 'M':0, 'F':0, 'Y':0, 'W':0}
    
    with open(iupred_file_name) as fileobject:
        for line in fileobject:
            if line[0] == '#':
                continue
            row = line.split('\t')
            aa = row[1]
            if float(row[2]) >= disorder_cutoff:
                disorder_flag = True
                len_disorder += 1
                try:
                    disorder_aa_comp[aa] = disorder_aa_comp[aa] + 1
                except:
                    pass
            else:
                if disorder_flag == True:
                    disorder_list.append(len_disorder)
                    tot_num_disorder += 1
                    len_disorder = 0
                disorder_flag = False
            if float(row[3][:-1]) >= disorder_cutoff:
                anchor_flag = True
                len_anchor += 1
                try:
                    anchor_aa_comp[aa] = anchor_aa_comp[aa] +1
                except:
                    pass
            else:
                if anchor_flag == True:
                    anchor_list.append(len_anchor)
                    tot_num_anchor +=1
                    len_anchor = 0
                anchor_flag = False
    return disorder_list,tot_num_disorder, disorder_aa_comp, anchor_list, tot_num_anchor, anchor_aa_comp         

def get_stats_from_hhpred(file_name):
    '''
    Parses hhpred data for inclusion in the structural features database
    Inputs:
        file_name (str): path to the file with the contents of hhpred for a gene or protein of interest
    Outputs:
        total_length (int): length of the protein in amino acids
        totAA_dict (dict): dictionary of amino acid frequencies in the input gene or protein
    '''
    total_length = 0
    totAA_dict = {'R':0, 'H':0, 'K':0, 'D':0, 'E':0, 'S':0, 'T':0, 'N':0, 'Q':0, 'C':0, 'U':0, 'G':0, 'P':0, 'A':0, 'V':0, 'I':0, 'L':0, 'M':0, 'F':0, 'Y':0, 'W':0}
    with open(file_name) as fileobject:
        for line in fileobject:
            if line[0] != '>':
                for elt in line:
                    if elt != '\n':
                        total_length += 1
                        if elt in totAA_dict:
                            totAA_dict[elt] = totAA_dict[elt] + 1   
    return total_length, totAA_dict

def get_tm(file_path):
    '''
    Parses transmembrane data for inclusion in the structural features database
    Inputs:
        file_path (str): path to the file with the contents of tm for a gene or protein of interest
    Outputs:
        tmyes (1 or 0): 1 if any tms present
        inside (list): list of lengths of helical regions inside membrane
        outside (list): list of lengths of helical regions outside membrane
        tmhelix (list): list of lengths of helical regions within membrane
        len(inside) (int): number of helicies inside the membrane
        len(outside) (int): number of helicies outside the membrane
        len(tmhelix) (int): number of helicies within the membrane
        sum(inside) (int): length of all helicies inside the membrane
        sum(outside) (int): length of all helicies outside the membrane
        sum(tmhelix) (int): length of all helicies within the membrane
    '''
    inside = []
    outside = []
    tmhelix = []
    with open(file_path) as fileobject:
        for line in fileobject:
            if line[0] != '#':
                parsed_line = line[:-1].split('\t')
                state = parsed_line[2]
                val_str = parsed_line[3].split(' ')
                first = True
                len_state = 1
                for element in val_str:
                    try:
                        pos = int(element)
                        if first:
                            first = False
                            len_state -= pos
                        else:
                            len_state += pos
                    except:
                        pass
                if state == 'inside':
                    inside.append(len_state)
                elif state == 'outside':
                    outside.append(len_state)
                elif state == 'TMhelix':
                    tmhelix.append(len_state)
    if len(tmhelix) > 0:
        tmyes = 1
    else:
        tmyes = 0
    return (tmyes, inside, outside, tmhelix, len(inside), len(outside), len(tmhelix), sum(inside), sum(outside), sum(tmhelix))

def yn_binary(num_to_check):
    '''
    Is a feature present?
    Inputs:
        num_to_check (int): feature frequency to check
    Outputs:
        1 or 0: 1 if the feature is present
    '''
    if num_to_check > 0:
        return 1
    return 0

def write_outputs(input_file):
    '''
    writes the files for the structural features database
    Inputs:
        input_file (str): path to file with all gene names or protiens to be included in the structural features database
    Outputs:
        error_report (list): list of lines containing information about which genes or protein data was not found when creating the database
    '''
    # # make look up dicts
    predict_pro_num_dir = './databases/numerical_predict_protein_values/'
    num_predict_pro_file_dict = make_file_dict(predict_pro_num_dir, '.txt')
    iupred_glob_dir = './databases/iupred2a_glob/'
    iupred_long_dir = './databases/iupred2a_long/'
    glob_file_dict = make_file_dict(iupred_glob_dir, '.fas.txt')
    iupred_file_dict = make_file_dict(iupred_long_dir, '.fas.txt') 
    hhpred_file_dir = './databases/human_proteome_hhpred/'
    hhpred_file_dict = make_file_dict(hhpred_file_dir, '.fas')
    tm_file_dir = './databases/tm_output/'
    tm_file_dict = make_file_dict(tm_file_dir, '.fas.txt')
    pred_pro_keys = ['Number of transmembrane helices', 'NHTM Best from query.phdPred', 'Stretch', 'Crowd predictions', 'Number of predictions', 'Number of positive regions', 'Number of positive regions with length >=30', 'Positive region lengths', 'Number of negative regions', 'Number of negative regions with length >=30', 'Negative region lengths', 'Number of coils', 'Total length of coil regions', 'Number amino acid in coil region A', 'Number amino acid in coil region R', 'Number amino acid in coil region N', 'Number amino acid in coil region D', 'Number amino acid in coil region C', 'Number amino acid in coil region E', 'Number amino acid in coil region Q', 'Number amino acid in coil region G', 'Number amino acid in coil region H', 'Number amino acid in coil region I', 'Number amino acid in coil region L', 'Number amino acid in coil region K', 'Number amino acid in coil region M', 'Number amino acid in coil region F', 'Number amino acid in coil region P', 'Number amino acid in coil region S', 'Number amino acid in coil region T', 'Number amino acid in coil region W', 'Number amino acid in coil region Y', 'Number amino acid in coil region V', 'Number amino acid in loop region A', 'Number amino acid in loop region R', 'Number amino acid in loop region N', 'Number amino acid in loop region D', 'Number amino acid in loop region C', 'Number amino acid in loop region E', 'Number amino acid in loop region Q', 'Number amino acid in loop region G', 'Number amino acid in loop region H', 'Number amino acid in loop region I', 'Number amino acid in loop region L', 'Number amino acid in loop region K', 'Number amino acid in loop region M', 'Number amino acid in loop region F', 'Number amino acid in loop region P', 'Number amino acid in loop region S', 'Number amino acid in loop region T', 'Number amino acid in loop region W', 'Number amino acid in loop region Y', 'Number amino acid in loop region V', 'Total length of loop regions', 'Number of loops', 'Number amino acid in sheet region A', 'Number amino acid in sheet region R', 'Number amino acid in sheet region N', 'Number amino acid in sheet region D', 'Number amino acid in sheet region C', 'Number amino acid in sheet region E', 'Number amino acid in sheet region Q', 'Number amino acid in sheet region G', 'Number amino acid in sheet region H', 'Number amino acid in sheet region I', 'Number amino acid in sheet region L', 'Number amino acid in sheet region K', 'Number amino acid in sheet region M', 'Number amino acid in sheet region F', 'Number amino acid in sheet region P', 'Number amino acid in sheet region S', 'Number amino acid in sheet region T', 'Number amino acid in sheet region W', 'Number amino acid in sheet region Y', 'Number amino acid in sheet region V', 'Total length of sheet regions', 'Number of sheets', 'Number amino acid in helix region A', 'Number amino acid in helix region R', 'Number amino acid in helix region N', 'Number amino acid in helix region D', 'Number amino acid in helix region C', 'Number amino acid in helix region E', 'Number amino acid in helix region Q', 'Number amino acid in helix region G', 'Number amino acid in helix region H', 'Number amino acid in helix region I', 'Number amino acid in helix region L', 'Number amino acid in helix region K', 'Number amino acid in helix region M', 'Number amino acid in helix region F', 'Number amino acid in helix region P', 'Number amino acid in helix region S', 'Number amino acid in helix region T', 'Number amino acid in helix region W', 'Number amino acid in helix region Y', 'Number amino acid in helix region V', 'Total length of helix regions', 'Number of helix', 'Number amino acid in nonconserved region A', 'Number amino acid in nonconserved region R', 'Number amino acid in nonconserved region N', 'Number amino acid in nonconserved region D', 'Number amino acid in nonconserved region C', 'Number amino acid in nonconserved region E', 'Number amino acid in nonconserved region Q', 'Number amino acid in nonconserved region G', 'Number amino acid in nonconserved region H', 'Number amino acid in nonconserved region I', 'Number amino acid in nonconserved region L', 'Number amino acid in nonconserved region K', 'Number amino acid in nonconserved region M', 'Number amino acid in nonconserved region F', 'Number amino acid in nonconserved region P', 'Number amino acid in nonconserved region S', 'Number amino acid in nonconserved region T', 'Number amino acid in nonconserved region W', 'Number amino acid in nonconserved region Y', 'Number amino acid in nonconserved region V', 'Total length of nonconserved regions', 'Number of nonconserved regions', 'Number amino acid in conserved region A', 'Number amino acid in conserved region R', 'Number amino acid in conserved region N', 'Number amino acid in conserved region D', 'Number amino acid in conserved region C', 'Number amino acid in conserved region E', 'Number amino acid in conserved region Q', 'Number amino acid in conserved region G', 'Number amino acid in conserved region H', 'Number amino acid in conserved region I', 'Number amino acid in conserved region L', 'Number amino acid in conserved region K', 'Number amino acid in conserved region M', 'Number amino acid in conserved region F', 'Number amino acid in conserved region P', 'Number amino acid in conserved region S', 'Number amino acid in conserved region T', 'Number amino acid in conserved region W', 'Number amino acid in conserved region Y', 'Number amino acid in conserved region V', 'Total length of conserved regions', 'Number of conserved regions']
    output_keys = ['Crowd predictions', 'Length of protein', 'NHTM Best from query.phdPred', 'Negative region lengths', 'Number amino acid in anchor region A', 'Number amino acid in anchor region C', 'Number amino acid in anchor region D', 'Number amino acid in anchor region E', 'Number amino acid in anchor region F', 'Number amino acid in anchor region G', 'Number amino acid in anchor region H', 'Number amino acid in anchor region I', 'Number amino acid in anchor region K', 'Number amino acid in anchor region L', 'Number amino acid in anchor region M', 'Number amino acid in anchor region N', 'Number amino acid in anchor region P', 'Number amino acid in anchor region Q', 'Number amino acid in anchor region R', 'Number amino acid in anchor region S', 'Number amino acid in anchor region T', 'Number amino acid in anchor region V', 'Number amino acid in anchor region W', 'Number amino acid in anchor region Y', 'Number amino acid in coil region A', 'Number amino acid in coil region C', 'Number amino acid in coil region D', 'Number amino acid in coil region E', 'Number amino acid in coil region F', 'Number amino acid in coil region G', 'Number amino acid in coil region H', 'Number amino acid in coil region I', 'Number amino acid in coil region K', 'Number amino acid in coil region L', 'Number amino acid in coil region M', 'Number amino acid in coil region N', 'Number amino acid in coil region P', 'Number amino acid in coil region Q', 'Number amino acid in coil region R', 'Number amino acid in coil region S', 'Number amino acid in coil region T', 'Number amino acid in coil region V', 'Number amino acid in coil region W', 'Number amino acid in coil region Y', 'Number amino acid in conserved region A', 'Number amino acid in conserved region C', 'Number amino acid in conserved region D', 'Number amino acid in conserved region E', 'Number amino acid in conserved region F', 'Number amino acid in conserved region G', 'Number amino acid in conserved region H', 'Number amino acid in conserved region I', 'Number amino acid in conserved region K', 'Number amino acid in conserved region L', 'Number amino acid in conserved region M', 'Number amino acid in conserved region N', 'Number amino acid in conserved region P', 'Number amino acid in conserved region Q', 'Number amino acid in conserved region R', 'Number amino acid in conserved region S', 'Number amino acid in conserved region T', 'Number amino acid in conserved region V', 'Number amino acid in conserved region W', 'Number amino acid in conserved region Y', 'Number amino acid in disordered region A', 'Number amino acid in disordered region C', 'Number amino acid in disordered region D', 'Number amino acid in disordered region E', 'Number amino acid in disordered region F', 'Number amino acid in disordered region G', 'Number amino acid in disordered region H', 'Number amino acid in disordered region I', 'Number amino acid in disordered region K', 'Number amino acid in disordered region L', 'Number amino acid in disordered region M', 'Number amino acid in disordered region N', 'Number amino acid in disordered region P', 'Number amino acid in disordered region Q', 'Number amino acid in disordered region R', 'Number amino acid in disordered region S', 'Number amino acid in disordered region T', 'Number amino acid in disordered region V', 'Number amino acid in disordered region W', 'Number amino acid in disordered region Y', 'Number amino acid in globular region A', 'Number amino acid in globular region C', 'Number amino acid in globular region D', 'Number amino acid in globular region E', 'Number amino acid in globular region F', 'Number amino acid in globular region G', 'Number amino acid in globular region H', 'Number amino acid in globular region I', 'Number amino acid in globular region K', 'Number amino acid in globular region L', 'Number amino acid in globular region M', 'Number amino acid in globular region N', 'Number amino acid in globular region P', 'Number amino acid in globular region Q', 'Number amino acid in globular region R', 'Number amino acid in globular region S', 'Number amino acid in globular region T', 'Number amino acid in globular region V', 'Number amino acid in globular region W', 'Number amino acid in globular region Y', 'Number amino acid in helix region A', 'Number amino acid in helix region C', 'Number amino acid in helix region D', 'Number amino acid in helix region E', 'Number amino acid in helix region F', 'Number amino acid in helix region G', 'Number amino acid in helix region H', 'Number amino acid in helix region I', 'Number amino acid in helix region K', 'Number amino acid in helix region L', 'Number amino acid in helix region M', 'Number amino acid in helix region N', 'Number amino acid in helix region P', 'Number amino acid in helix region Q', 'Number amino acid in helix region R', 'Number amino acid in helix region S', 'Number amino acid in helix region T', 'Number amino acid in helix region V', 'Number amino acid in helix region W', 'Number amino acid in helix region Y', 'Number amino acid in loop region A', 'Number amino acid in loop region C', 'Number amino acid in loop region D', 'Number amino acid in loop region E', 'Number amino acid in loop region F', 'Number amino acid in loop region G', 'Number amino acid in loop region H', 'Number amino acid in loop region I', 'Number amino acid in loop region K', 'Number amino acid in loop region L', 'Number amino acid in loop region M', 'Number amino acid in loop region N', 'Number amino acid in loop region P', 'Number amino acid in loop region Q', 'Number amino acid in loop region R', 'Number amino acid in loop region S', 'Number amino acid in loop region T', 'Number amino acid in loop region V', 'Number amino acid in loop region W', 'Number amino acid in loop region Y', 'Number amino acid in nonconserved region A', 'Number amino acid in nonconserved region C', 'Number amino acid in nonconserved region D', 'Number amino acid in nonconserved region E', 'Number amino acid in nonconserved region F', 'Number amino acid in nonconserved region G', 'Number amino acid in nonconserved region H', 'Number amino acid in nonconserved region I', 'Number amino acid in nonconserved region K', 'Number amino acid in nonconserved region L', 'Number amino acid in nonconserved region M', 'Number amino acid in nonconserved region N', 'Number amino acid in nonconserved region P', 'Number amino acid in nonconserved region Q', 'Number amino acid in nonconserved region R', 'Number amino acid in nonconserved region S', 'Number amino acid in nonconserved region T', 'Number amino acid in nonconserved region V', 'Number amino acid in nonconserved region W', 'Number amino acid in nonconserved region Y', 'Number amino acid in protein A', 'Number amino acid in protein C', 'Number amino acid in protein D', 'Number amino acid in protein E', 'Number amino acid in protein F', 'Number amino acid in protein G', 'Number amino acid in protein H', 'Number amino acid in protein I', 'Number amino acid in protein K', 'Number amino acid in protein L', 'Number amino acid in protein M', 'Number amino acid in protein N', 'Number amino acid in protein P', 'Number amino acid in protein Q', 'Number amino acid in protein R', 'Number amino acid in protein S', 'Number amino acid in protein T', 'Number amino acid in protein V', 'Number amino acid in protein W', 'Number amino acid in protein Y', 'Number amino acid in sheet region A', 'Number amino acid in sheet region C', 'Number amino acid in sheet region D', 'Number amino acid in sheet region E', 'Number amino acid in sheet region F', 'Number amino acid in sheet region G', 'Number amino acid in sheet region H', 'Number amino acid in sheet region I', 'Number amino acid in sheet region K', 'Number amino acid in sheet region L', 'Number amino acid in sheet region M', 'Number amino acid in sheet region N', 'Number amino acid in sheet region P', 'Number amino acid in sheet region Q', 'Number amino acid in sheet region R', 'Number amino acid in sheet region S', 'Number amino acid in sheet region T', 'Number amino acid in sheet region V', 'Number amino acid in sheet region W', 'Number amino acid in sheet region Y', 'Number of anchor regions', 'Number of coils', 'Number of conserved regions', 'Number of disordered regions', 'Number of globular regions', 'Number of helix', 'Number of loops', 'Number of negative regions', 'Number of negative regions with length >=30', 'Number of nonconserved regions', 'Number of positive regions', 'Number of positive regions with length >=30', 'Number of predictions', 'Number of sheets', 'Number of transmembrane helices', 'Positive region lengths', 'Stretch', 'Total length of anchor regions', 'Total length of coil regions', 'Total length of conserved regions', 'Total length of disordered regions', 'Total length of globular regions', 'Total length of helix regions', 'Total length of loop regions', 'Total length of nonconserved regions', 'Total length of sheet regions', 'Total length tmh regions', 'Y/n anchor regions', 'Y/n disordered regions', 'Y/n globular regions', 'Y/n tmh regions']
    aa_list = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

    error_report = []
    current_line = 0
    tot_lines = 116681
    # Look through list of gns
    with open(input_file) as fileobject:
        for line in fileobject:
            print(str(current_line/tot_lines) + '% done')
            current_line += 1
            name = line[:-1]
            out_dict = {}
            try:
                length, totAA_dict = get_stats_from_hhpred(hhpred_file_dict[name])
                length = max(length, sum(totAA_dict.values()))
                out_dict['Length of protein'] = length
                for elt in aa_list:
                    out_dict['Number amino acid in protein ' + elt] = totAA_dict[elt]
            except:
                length = 0
                out_dict['Length of protein'] = 0
                for elt in aa_list:
                    out_dict['Number amino acid in protein ' + elt] = 0
                error_report.append([name + 'hhpred_stats'])
            
            try:
                glob_file = glob_file_dict[name]
                number_globs, glob_len_list, glob_comp_dict = get_glob(glob_file)
                globsyes = yn_binary(number_globs)
                tot_glob_len = max(sum(glob_len_list), sum(glob_comp_dict.values()))
                if tot_glob_len <= length:
                    out_dict['Y/n globular regions'] = globsyes
                    out_dict['Number of globular regions'] = number_globs
                    out_dict['Total length of globular regions'] = tot_glob_len
                    for element in aa_list:
                        out_dict['Number amino acid in globular region ' + element] = glob_comp_dict[element]
                else:
                    for element in aa_list:
                        out_dict['Number amino acid in globular region ' + element] = 0
                    out_dict['Y/n globular regions'] = 0
                    out_dict['Number of globular regions'] = 0
                    out_dict['Total length of globular regions'] = 0
                    error_report.append([name + ',glob'])
            except:
                for element in aa_list:
                    out_dict['Number amino acid in globular region ' + element] = 0
                out_dict['Y/n globular regions'] = 0
                out_dict['Number of globular regions'] = 0
                out_dict['Total length of globular regions'] = 0
                error_report.append([name + ',glob'])

            try:
                disorder_file = iupred_file_dict[name]
                disorder_list,tot_num_disorder, disorder_aa_comp, anchor_list, tot_num_anchor, anchor_aa_comp = get_disorder(disorder_file, 0.5)
                disorderyes = yn_binary(tot_num_disorder)
                tot_disorder_len = max(sum(disorder_list), sum(disorder_aa_comp.values()))
                if tot_disorder_len <= length:
                    out_dict['Y/n disordered regions'] = disorderyes
                    out_dict['Number of disordered regions'] = tot_num_disorder
                    out_dict['Total length of disordered regions'] = tot_disorder_len
                    for element in aa_list:
                        out_dict['Number amino acid in disordered region ' + element] = disorder_aa_comp[element]
                else:
                    out_dict['Y/n disordered regions'] = 0
                    out_dict['Number of disordered regions'] = 0
                    out_dict['Total length of disordered regions'] = 0
                    for element in aa_list:
                        out_dict['Number amino acid in disordered region ' + element] = 0
                    error_report.append([name + ',disordered'])
                
                anchoryes = yn_binary(tot_num_anchor)
                tot_anchor_len = max(sum(anchor_list), sum(anchor_aa_comp.values()))
                if tot_anchor_len <= length:
                    out_dict['Y/n anchor regions'] = anchoryes
                    out_dict['Number of anchor regions'] = tot_num_anchor
                    out_dict['Total length of anchor regions'] = tot_anchor_len
                    for element in aa_list:
                        out_dict['Number amino acid in anchor region ' + element] = anchor_aa_comp[element]
                else:
                    out_dict['Y/n anchor regions'] = 0
                    out_dict['Number of anchor regions'] = 0
                    out_dict['Total length of anchor regions'] = 0
                    for element in aa_list:
                        out_dict['Number amino acid in anchor region ' + element] = 0
                    error_report.append([name + ',anchor'])
            except:
                out_dict['Y/n disordered regions'] = 0
                out_dict['Number of disordered regions'] = 0
                out_dict['Total length of disordered regions'] = 0
                for element in aa_list:
                    out_dict['Number amino acid in disordered region ' + element] = 0
                out_dict['Y/n anchor regions'] = 0
                out_dict['Number of anchor regions'] = 0
                out_dict['Total length of anchor regions'] = 0
                for element in aa_list:
                    out_dict['Number amino acid in anchor region ' + element] = 0
                error_report.append([name + ',iupred_disorder'])
            
            if name in num_predict_pro_file_dict:
                output_predpro = num_predict_pro_parse(num_predict_pro_file_dict[name])
                checks = max([
                    output_predpro['Positive region lengths'], 
                    output_predpro['Negative region lengths'],
                    output_predpro['Total length of coil regions'],
                    output_predpro['Total length of helix regions'],
                    output_predpro['Total length of sheet regions'],
                    output_predpro['Total length of loop regions'],
                    output_predpro['Total length of conserved regions'],
                    output_predpro['Total length of nonconserved regions']
                ])
                if checks <= length:
                    out_dict = {**output_predpro, **out_dict}
                else:
                    error_report.append([name + ',predictpro'])
                    for k in pred_pro_keys:
                        out_dict[k] = 0
            else:
                for k in pred_pro_keys:
                    out_dict[k] = 0
                error_report.append([name + ',predictpro'])
            if name in tm_file_dict:
                tmyes, inside, outside, tmhelix, num_inside, num_outside, num_tmhelix, tot_len_inside, tot_len_outside, tot_len_tmhelix = get_tm(tm_file_dict[name])
                if tot_len_tmhelix <=length:
                    out_dict['Total length tmh regions'] = tot_len_tmhelix
                    out_dict['Y/n tmh regions'] = tmyes
                    out_dict['Number of transmembrane helices'] = num_tmhelix

                else:
                    out_dict['Total length tmh regions'] = 0
                    out_dict['Y/n tmh regions'] = 0
                    out_dict['Number of transmembrane helices'] = 0
                    error_report.append([name + ',tmh'])
            else:
                error_report.append([name + ',tmh'])
                out_dict['Total length tmh regions'] = 0
                out_dict['Y/n tmh regions'] = 0
                if 'Number of transmembrane helices' not in out_dict:
                    out_dict['Number of transmembrane helices'] = 0
            out_line = name + ','
            for key in output_keys:
                out_line = out_line + str(out_dict[key]) + ','
            f = open('databases/precounted_human_genome/'+name+'.txt', 'w')
            f.write(out_line[:-1])
            f.close() 
        return error_report

if __name__ == '__main__':
    input_file ='databases/all_ids.txt'
    error_report = write_outputs(input_file)
    pd.DataFrame(error_report).to_csv('unfound.csv', index = None, header= None)