#!/usr/bin/env python

###################################################
### Structural Features Background Generation   ###
### Written by: Nicole Zatorski                 ###
### Last date modified: 22/7/22                 ###
###################################################

# import statements
import sys
import os
import glob

# functions
def write_output(out_name, out_str):
    '''
    Appends the out_str to the end of the file described by out_name
    Inputs:
        out_name: (str) path and name of file to write to
        out_str: (str) information to append to the end of the file
    Outputs:
        None
    '''
    f = open(out_name, 'a+')
    f.write(out_str)
    f.close()

def make_background(input_folder, file_id, output_dir_name):
    '''
    Converts the output from structural features into background that structural features can use
    Inputs:
        input_folder: (str) path where to find the file the background will be generated from
        file_id: (str) name of the file the bacground will be generated from
        output_dir_name: (str) path and name of file to write to, doesn't need to already exist but can
    Outputs:
        1 (int) when complete
    '''
    if not os.path.exists('./databases/'+output_dir_name + '/'):
        os.makedirs('./databases/'+output_dir_name + '/')
    output_dir = './databases/'+output_dir_name + '/'
    out_average = ''
    index = 0
    with open(input_folder + '/average_'+file_id+'.csv') as fo:
        for line in fo:
            if index != 0:
                split_line = line[:-1].split(',')
                out_average = out_average + split_line[0] + ',' + split_line[1] + ',' + split_line[2] + '\n'
            index +=1
    frequency = ''
    domain = ''
    fam = ''
    superfam =''
    fold = ''
    with open(input_folder + '/frequency_'+file_id+'.csv') as fo:
        for line in fo:
            split_line = line[:-1].split(',')
            if split_line[2] == 'N/A':
                frequency = frequency + split_line[1] + ',' + split_line[3] + '\n'
            elif split_line[2] == 'domain':
                domain = domain + split_line[0] + ',' + split_line[3] + '\n'
            elif split_line[2] == 'fold':
                fold = fold + split_line[0] + ',' + split_line[3] + '\n'
            elif split_line[2] == 'family':
                fam = fam + split_line[0] + ',' + split_line[3] + '\n'
            elif split_line[2] == 'superfamily':
                superfam = superfam + split_line[0] + ',' + split_line[3] + '\n'
    write_output(output_dir +'average_background.csv' ,out_average)
    write_output(output_dir +'frequency_background.csv' ,frequency)
    write_output(output_dir +'scop.fold.csv' ,fold)
    write_output(output_dir +'ipr.domain.csv' ,domain)
    write_output(output_dir +'scop.family.csv' ,fam)
    write_output(output_dir +'scop.superfam.csv' ,superfam)

    return 1

if __name__ == '__main__':
    make_background('prot_backout_cmap', 'prot_over', 'prot_over')
    make_background('prot_backout_cmap', 'prot_under', 'prot_under')
    # make_background('norm_breast', 'gtex_breast', 'gtex_breast250')
    # make_background(*sys.argv[1:])
    # make_background('nonbreast_gtex_out/', 'all_ids', 'updated_human_background')
    # files = glob.glob('rna_backout/*normal*')

    # # files = glob.glob('proteomics_backout/*')
    # for f in files:
    #     # id = f.split('_')[-1].split('.')[0]
    #     # make_background('proteomics_backout', id, id)


    #     id = f.split('_')[-1].split('-')[0]
    #     make_background('rna_backout', id+'-normal-tissue', id)