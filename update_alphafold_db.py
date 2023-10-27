import glob

def write_output(out_name, out_str):
    '''
    Write output file from string
    Inputs:
        out_name (str): file to write to
        out_str (str): content to write to file
    Outputs:
        None
    '''
    f = open(out_name, 'w+')
    f.write(out_str)
    f.close()

def get_afdssp(gnuid):
    '''
    gets DSSP secondary structures from gene or uniprot ID
    Inputs:
        gnuid (str): gene or uniprot ID
    Outputs:
        all_feats (dict of dicts): DSSP feature dictionary
    '''
    path_to_db = 'databases/af_dssp/'
    line_index = -1
    all_feats = {
        'S':{'length':0, 'number':0},
        'E':{'length':0, 'number':0},
        'T':{'length':0, 'number':0},
        'B':{'length':0, 'number':0},
        'G':{'length':0, 'number':0},
        'H':{'length':0, 'number':0}
    }
    previous_feature = '-'
    with open(path_to_db + gnuid + '-F1.csv') as fo:
        for line in fo:
            line_index += 1
            if line_index != 0:
                split_line = line[:-1].split(',')
                if split_line[2] in all_feats:
                    feat = split_line[2]
                    aa = split_line[1]
                    if aa in all_feats[feat]:
                        all_feats[feat][aa] = all_feats[feat][aa] + 1
                    else:
                        all_feats[feat][aa] = 1
                    all_feats[feat]['length'] = all_feats[feat]['length'] + 1
                    if previous_feature != feat:
                        all_feats[feat]['number'] = all_feats[feat]['number'] + 1
                        previous_feature = feat
    return all_feats

def get_a3d(gnuid):
    '''
    gets a3d features from gene or uniprot ID
    Inputs:
        gnuid (str): gene or uniprot ID
    Outputs:
        out_dict (dict): feature dictionary
    '''
    path_to_db = 'databases/A3D_scores/AF-'
    line_index = -1
    out_dict = {'length':0, 'number':0}
    previous_feature = -1
    with open(path_to_db + gnuid + '-F1_A3D.csv') as fo:
        for line in fo:
            line_index += 1
            if line_index != 0:
                split_line = line[:-1].split(',')
                if float(split_line[3]) > 0:
                    aa = split_line[1][-1]
                    if aa in out_dict:
                        out_dict[aa] = out_dict[aa] + 1
                    else:
                        out_dict[aa] = 1
                    if previous_feature <= 0 :
                        out_dict['number'] = out_dict['number'] + 1
                    out_dict['length'] = out_dict['length'] + 1
                previous_feature = float(split_line[3])
    return out_dict

def get_com(gnuid):
    '''
    gets center of mass features from gene or uniprot ID
    Inputs:
        gnuid (str): gene or uniprot ID
    Outputs:
        min (float): minimum distance of a residue to the center of mass
        max (float): maximum distance of a residue to the center of mass
        average/line_index (float): average distance of a residue to the center of mass
    '''
    path_to_db = 'databases/centerofmass/'
    line_index = -1
    min = 100000000000000
    max = 0
    average = 0
    previous_feature = -1
    with open(path_to_db + gnuid + '_dcom.csv') as fo:
        for line in fo:
            line_index += 1
            if line_index != 0:
                split_line = line[:-1].split(',')
                val = float(split_line[1])
                if val < min:
                    min= val
                if val > max:
                    max = val
                average += val
    return min,max,average/line_index

def get_num_contacts(gnuid):
    '''
    gets nunber of amino acid contacts from gene or uniprot ID
    Inputs:
        gnuid (str): gene or uniprot ID
    Outputs:
        num_contacts (int): nunber of amino acid contacts in the protein
    '''
    path_to_db = 'databases/contacts/'
    line_index = -1
    num_contacts = 0
    tot = 0
    with open(path_to_db + gnuid + '_map.csv') as fo:
        for line in fo:
            line_index += 1
            if line_index != 0:
                split_line = line[:-1].split(',')
                tot += len(split_line)
                for pos in split_line[line_index:]:
                    if pos == 'True':
                        num_contacts += 1
    return num_contacts


aa_list = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
files_run = glob.glob('databases/precounted_alpha_fold/*')
with open('databases/uniprot-gn-map.txt') as fo:
    for line in fo:
        split_line = line[:-1].split('__')
        all_row_ids = split_line[0].split('--') + split_line[1].split('--')
        found = False
        ids_to_fix = []
        for id in all_row_ids:
            try:
                fo2 = open('databases/precounted_alpha_fold/'+id + '.txt', 'r')
                out = fo2.read()
                fo2.close()
                found = True
                break
            except:
                pass
        for round2id in all_row_ids:
            if found:
                try:
                    fo3 = open('databases/precounted_alpha_fold/'+round2id + '.txt', 'r')
                    fo3.close()
                except:
                    new_out = round2id + ',' + out.split(',',1)[-1]
                    write_output('databases/precounted_alpha_fold/'+round2id + '.txt', new_out)
