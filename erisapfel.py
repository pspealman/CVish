# -*- coding: utf-8 -*-
"""
erisapfel - "apple of eris, goddess of discord"
# Purpose: Uses split and discordant reads to identify genomic breakpoints
# further processing is used to extract context and sequence.

#change log
(Initial Release - Public Beta: Rider Theory)

# 12.13.2021
Version 1.0 (Write Reflection)
    Change log
    _x_ 'gff' error
    _x_ updated requirement list
    _x_ corrected RPKM

Version 1.1 (Snatch Hen) 01.19.2022
    _x_ mkdir error
    _x_ gzip file handling in subfolders
    _x_ -filter argument on depth command now optional
    _x_ corrected module integration
    _x_ removed soft reads 
    _x_ Mobile genetic element (MGE) mapping
    _x_ Added 'separate_by_strand' option
        This seperates reads into new sam/bam files based on strand identificaiton.
        Visualization option
        
Version 1.2 (Area Agenda) 02.23.2022
    _x_ ML model support for MGE prediction
        _x_ model selection
    _x_ vcf format output (v1.2.2 - 03.12.20)
        _x_ with bgzip and tbi indexing
    
    
Future versions:
    ___ Added 'find_breakpoints'
    ___ Add CLI score filter so lines aren't added to the gff and filtered from the vcf
    ___ Check contig reporting in gff, vcf

@author: Pieter Spealman ps163@nyu.edu
"""

"""Requirements:
bwa/intel   0.7.17
samtools    1.14
bedtools    2.29.2
velvet      1.2.10
blast+      2.11.0
samblaster  0.1.26
mafft       7.475
emboss      6.6.0
"""
#
import os
import argparse
import subprocess
import pickle
import pandas as pd
import numpy as np
import json
import re
import time
from datetime import date
import scipy.stats as stats

import warnings
warnings.simplefilter("ignore")

#01.08.2021 - loading modules 
#exec(open('/usr/share/Modules/init/python.py').read())

#function handles help commands
def help_dialog():
    monolog=('Manual for erisapfel\n'+
    '#=============================================================#\n'+
    'Created by Pieter Spealman\n ps163@nyu.edu \n'+
    'Release version: 0.9 \n'+
    'Release date: 10.05.19 \n'+
    '\tThis is a multistep package that identifies CNV breakpoints using split\n'+
    'and discordant sequencing reads.\n'+
    '#=============================================================#\n'+
    'Usage:\n\tCNV identification with ancestor.\n'+
    'This method should be used for experiments involving evolved, adapted, or\n'+
    'normal/tumor samples.\n'+
    'For more information refer to docs/README.txt\n'+
    'For demonstration use:\n\t python erisapfel.py -demo\n'+
    'To run a install test using defaults, use:\n\t python erisapfel.py -test\n'+
    '')
    print(monolog)

def demo():
    monolog = ('Usage:\n\tMethod 2. CNV identification with ancestor.\n'+
    'This method should be used for experiments involving evolved, adapted, or\n'+
    'normal/tumor samples.\n'+
    'Example:\n'+
    'Step 1. Generate both split and discordant read files for ancestor strain:\n\t'+
    'demo: python erisapfel.py -make -fa demo/input/demo.fna -gff demo/input/demo.gff -fastq_1 demo/input/n01_ancestor.fastq.gz -fastq_2 demo/input/n02_ancestor.fastq.gz -o demo/output/demo_anc \n\n'+
    'Step 2. Generate both split and discordant read files for evolved strain:\n\t'+
    'demo: python erisapfel.py -make -fa demo/input/demo.fna -gff demo/input/demo.gff -fastq_1 demo/input/n01_evolved.fastq.gz -fastq_2 demo/input/n02_evolved.fastq.gz -o demo/output/demo_evo \n\n'+
    'Step 3. Identify Breakpoints in Ancestor using -breakpoint\n\t'+
    'demo: python erisapfel.py -breakpoint -i demo/output/demo_anc -break_tab demo/output/demo_ancestor.tab -break_bed demo/output/demo_ancestor.bed \n\n'+
    'Step 4. Build Filter, including Ancestor breakpoints, using -filter\n\t'+ 
    'demo: python erisapfel.py -filter -gff demo/input/demo.gff -ancestor demo/output/demo_ancestor.bed -o demo/output/filter.p\n\n'+
    'Step 5. Identify Breakpoints in Adapted population, and applying filtering, using -breakpoint\n\t'+
    'demo: python erisapfel.py -breakpoint -filter_object demo/output/filter.p -i demo/output/demo_evo -break_tab demo/output/demo_evolved.tab -break_bed demo/output/demo_evolved.bed\n\n'+
    'Step 6. Derive the breakpoint contigs from split and / or discordant reads, using -localseq\n\t'+ 
    'demo: python erisapfel.py -localseq -break_tab demo/output/demo_evolved.tab -break_bed demo/output/demo_evolved.bed -fastq_1 demo/input/n01_evolved.fastq.gz -fastq_2 demo/input/n02_evolved.fastq.gz -o demo/input/localseq_evo \n\n'+
    'Step 7. Realign contigs back into reference genome to identify breakpoint origins\n\t'+
    'demo: python erisapfel.py -realign -fa demo/input/demo.fna -contigs demo/input/localseq_evo/ -min_coverage 1.5 -max_eval 0.05 -o demo/input/localseq_evo/realign_evo\n\n'+
    ''+
    '    python erisapfel.py -depth -mpileup -fa /scratch/work/cgsb/genomes/Public/Fungi/Saccharomyces_cerevisiae/NCBI/R64/GCF_000146045.2_R64_genomic.fna -bam_files grace/input/bc_c01/01 -o try_depth/01'
    '    python erisapfel.py -peaks -read_type RD -i try_depth/01 -o try_depth/01'

    'python erisapfel.py -lof -read_type discordant -train g250_39_v2 g250_40_v2 g250_41_v2 g250_42_v2 g250_49_v2 -i try_more -o try_more'+
    'python erisapfel.py -view -i try_more -o try_more_foo'+

"""

"""
    'python erisapfel.py -localseq -break_tab demo/output/demo_evolved.tab -break_bed demo/output/demo_evolved.bed -fastq_1 demo/input/n01_evolved.fastq.gz -fastq_2 demo/input/n02_evolved.fastq.gz -o demo/input/localseq_evo'+

    '')
    print(monolog)

#chromo_list = ['NC_001135.5']
#chromo_list = ['NC_001135.5', 'NC_001143.9']
#chromo_list = ['NC_001133.9','NC_001134.8','NC_001141.2','NC_001135.5', 'NC_001143.9']
#
#chromo_list = ['NC_001133.9','NC_001134.8','NC_001135.5','NC_001136.10','NC_001137.3','NC_001138.5','NC_001139.9','NC_001140.6', 
#'NC_001141.2','NC_001142.9','NC_001143.9','NC_001144.5','NC_001145.3','NC_001146.8','NC_001147.6','NC_001148.4']

def handle_outfile(p_output):
            
    if '/' in p_output:
        output_dir = str(p_output.rsplit('/',1)[0]+'/')
        output_file = p_output.rsplit('/',1)[1]
        
        if len(output_dir.strip())>=1:
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
    else:
        output_dir = ''
        output_file = args.output_file
        
    return(output_dir, output_file)  

def io_initialize():
    #Create default directories
    def make_subdirectory(is_path,subdir):
        new_dir = ('{}/{}').format(is_path, subdir)
        if not os.path.exists(new_dir):
            os.makedirs(new_dir)
    
    #Define default directories
    s_path = ('{}{}').format(output_dir, output_file)
    bam_dir = ('{}/bam/').format(s_path)
    final_output_dir = ('{}/output/').format(s_path)
    temp_dir = ('{}/temp/').format(s_path)
    pickles_dir = ('{}/pickles/').format(s_path)
    
    make_subdirectory(s_path,'bam')
    make_subdirectory(s_path,'output')
    make_subdirectory(s_path,'temp')
    make_subdirectory(s_path,'pickles')
    
    return(s_path, bam_dir, final_output_dir, temp_dir, pickles_dir)

def io_make(runmode='resource'):
    if runmode == 'resource':
        outline = ('Files and parameters used are stored in {}/resource.tab').format(final_output_dir)
        print(outline)
        
        resource_file_name = ('{}/resource.tab').format(final_output_dir)
        resource_file = open(resource_file_name, 'w')
        
        resource_file.close()
        
    if runmode == 'ml':
        ml_file_name = ('{}/ml_hypo.tab').format(temp_dir)
        ml_file = open(ml_file_name, 'w')

        for cat, val in ml_hypo_dict.items():                
            outline = ('{}\t{}\n').format(cat,val)
            ml_file.write(outline)
        ml_file.close()
        
        ml_file_name = ('{}/ml_loci.tab').format(temp_dir)
        ml_file = open(ml_file_name, 'w')

        for cat, val in ml_loci_dict.items():                
            outline = ('{}\t{}\n').format(cat,val)
            ml_file.write(outline)
        ml_file.close()
                     
def test():
    #TODO make accurate for latest version
    monolog = ('=== Currently Testing Erisapfel.py ===')
    print(monolog)
    
    monolog = ('\tTesting Step 0. make bwa for demo.fna\n')
    print(monolog)    
    bashCommand = ('bwa index demo/input/demo.fna')
    print(bashCommand)       
    subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)
    
    monolog = ('\tTesting Step 1. -make command for ancestor\n')
    print(monolog)    
    bashCommand = ('python erisapfel.py -make -fa demo/input/demo.fna -gff demo/input/demo.gff -fastq_1 demo/input/n01_ancestor.fastq.gz -fastq_2 demo/input/n02_ancestor.fastq.gz -o demo/output/demo_anc')
    print(bashCommand)       
    subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)
    
    monolog = ('\tTesting Step 2. -make command for evolved\n')
    print(monolog)    
    bashCommand = ('python erisapfel.py -make -fa demo/input/demo.fna -gff demo/input/demo.gff -fastq_1 demo/input/n01_evolved.fastq.gz -fastq_2 demo/input/n02_evolved.fastq.gz -o demo/output/demo_evo')
    print(bashCommand)       
    subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)
    
    monolog = ('\tTesting Step 3. -breakpoint command for ancestor\n')
    print(monolog)    
    bashCommand = ('python erisapfel.py -breakpoint -i demo/output/demo_anc -break_tab demo/output/demo_ancestor.tab -break_bed demo/output/demo_ancestor.bed')
    print(bashCommand)       
    subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)

    monolog = ('\tTesting Step 4. -filter command using ancestor\n')
    print(monolog)    
    bashCommand = ('python erisapfel.py -filter -gff demo/input/demo_filter.gff -ancestor demo/output/demo_ancestor.bed -o demo/output/filter.p')
    print(bashCommand)       
    subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)
    
    monolog = ('\tTesting Step 5. -breakpoint command for evolved\n')
    print(monolog)    
    bashCommand = ('python erisapfel.py -breakpoint -filter_object demo/output/filter.p -i demo/output/demo_evo -break_tab demo/output/demo_evolved.tab -break_bed demo/output/demo_evolved.bed')
    print(bashCommand)       
    subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)
    
    monolog = ('\tTesting Step 6. -localseq command\n')
    print(monolog)    
    bashCommand = ('python erisapfel.py -localseq -split_only -break_tab demo/output/demo_evolved.tab -break_bed demo/output/demo_evolved.bed -fastq_1 demo/input/n01_evolved.fastq.gz -fastq_2 demo/input/n02_evolved.fastq.gz -o demo/input/localseq_evo')
    print(bashCommand)       
    subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)
    
    monolog = ('\tTesting Step 7. -realign command\n')
    print(monolog)    
    bashCommand = ('python erisapfel.py -realign -fa demo/input/demo.fna -contigs demo/input/localseq_evo/ -o demo/input/localseq_evo/realign_evo')
    print(bashCommand)       
    (subprocess.run([bashCommand],stderr= subprocess.STDOUT,shell=True))
    
parser = argparse.ArgumentParser(add_help=True)
#handle help_dialog
parser.add_argument('-man',"--manual", action='store_true')
parser.add_argument('-demo',"--demo",action='store_true')
parser.add_argument('-test',"--test",action='store_true')
parser.add_argument('-view',"--view_resource")

parser.add_argument('-run_name', '--run_name')

parser.add_argument('-make',"--make", action='store_true')
parser.add_argument('-ends',"--paired_end")
parser.add_argument('-fa',"--fa_file")
parser.add_argument('-bam', "--bam_file")

parser.add_argument('-filter', '--make_filter', action='store_true')
parser.add_argument('-filter_tab', '--filter_tab', nargs='+')
parser.add_argument('-filter_gff','--filter_gff', nargs='+')
parser.add_argument('-filter_bed', '--filter_bed', nargs='+')

parser.add_argument('-map', '--map_reads', action='store_true')
parser.add_argument('-contig',"--contig", action='store_true')
parser.add_argument('-z_split',"--split_zscore")
parser.add_argument('-z_disco',"--disco_zscore")

parser.add_argument('-qval', '--mapq_val')

parser.add_argument('-fo', "--filter_object")
parser.add_argument('-flank','--filter_flanking')

parser.add_argument('-break_tab', "--break_tab")
parser.add_argument('-break_bed', "--break_bed")
parser.add_argument('-split_only', '--split_only', action='store_true')

parser.add_argument('-chromo_size', "--chromo_size")

parser.add_argument('-localseq', '--build_sequence', action='store_true')
parser.add_argument('-hashsize', '--hashsize')
parser.add_argument('-min_contig','--min_contig')
parser.add_argument('-max_divergence','--max_divergence')
parser.add_argument('-cov_cutoff','--cov_cutoff')

parser.add_argument('-fastq_1', "--fastq_1")
parser.add_argument('-fastq_2', "--fastq_2")

parser.add_argument('-i', "--input_file")
parser.add_argument('-o', "--output_file")
parser.add_argument('-bed', "--bed_file")
parser.add_argument('-split_score', "--split_score")
parser.add_argument('-disco_score', "--disco_score")

parser.add_argument('-sw', "--split_weight")
parser.add_argument('-dw', "--disco_weight")

parser.add_argument('-realign',"--realign", action='store_true')
parser.add_argument('-realign_contigs',"--realign_contigs", action='store_true')
parser.add_argument('-contigs', "--contigs_file")
parser.add_argument('-min_coverage', "--min_coverage")
parser.add_argument('-max_eval', "--max_eval")
parser.add_argument('-overlap_mask','--overlap_mask')
parser.add_argument('-min_score','--min_score')

parser.add_argument('-depth', '--depth_analysis', action='store_true')

parser.add_argument('-read_type','--read_type_list', nargs='+')

parser.add_argument('-peaks', '--peaks', action='store_true')
parser.add_argument('-min_length', '--CNV_min_length')
parser.add_argument('-gap', '--max_gap')
parser.add_argument('-purity', '--min_purity')
parser.add_argument('-windowy', '--window')

parser.add_argument('-breakpoints','--find_breakpoints', action='store_true')

parser.add_argument('-strand', '--separate_by_strand')

#get feature values
#python erisapfel.py -mge /scratch/ps163/Project_Carolino/metadata/te_a_NCBI_Gap1.gff -run_name ${sample_name}
#python erisapfel.py -kmge /scratch/ps163/Project_Carolino/mugio/DGY1657_gap1_mge_mugio_predicted_TEa.gff -mge /scratch/ps163/Project_Carolino/metadata/te_a_NCBI_Gap1.gff -run_name ${sample_name}
parser.add_argument('-mge', "--mge_file")
parser.add_argument('-kmge', "--known_mge")

#model functions
#python erisapfel.py -ml_select -i temp/mge_value_df.tab --metaparameters fast -o
parser.add_argument('-ml_select', '--model_selection', action='store_true')
parser.add_argument('-ml_meta', '--metaparameters')

#python erisapfel.py -ml_train -i temp/mge_value_df.tab --preprocessor MinMaxScaler --sampler AllKNN --classifier LogisticRegression_lbgfs --model_object DGY_model.p
parser.add_argument('-ml_train', '--model_train', action='store_true')
parser.add_argument('-ml_p', '--preprocessor')
parser.add_argument('-ml_s', '--sampler')
parser.add_argument('-ml_c', '--classifier')

#python erisapfel.py -ml_predict --model_object DGY_model.p -run_name ${sample_name}
parser.add_argument('-ml_predict', '--model_predict', action='store_true')
parser.add_argument('-model', '--model_object')
parser.add_argument('-all', '--return_all', action='store_true')



args = parser.parse_args()

if args.manual:
    help_dialog()
    
if args.demo:
    demo()
    
if args.test:
    test()
    
if args.view_resource:
    resource_dict = pickle.load(open(args.view_resource, 'rb'))
    
    for key, value in resource_dict.items():
        print('key: ', key)
        print('value: ', value)
        
        
''' '''        
filter_region_dict = {}
region_name_dict = {}
hypothesis_score_dict = {}
region_inclusion_dict = {}
contig_seq_dict = {}

ml_hypo_dict = {}
ml_loci_dict = {}
ml_san_check = {}

processed_reads = 0
used_reads = 0
filtered_reads = 0

if args.filter_object and not args.make_filter:
    region_filter_dict = pickle.load(open(args.filter_object, 'rb'))

if args.filter_flanking:
    flank_nt = int(args.filter_flanking)
else:
    #This is presuming that a read length is 75
    flank_nt = 75
    
if args.output_file:
    output_dir, output_file = handle_outfile(args.output_file)
    
if args.run_name:
    make_name = ("results/{}").format(args.run_name)
    output_dir, output_file = handle_outfile(make_name)
    #TODO - does this break everything?    
    s_path, bam_dir, final_output_dir, temp_dir, pickles_dir = io_initialize()

# Define scores and weights
if args.split_score:
    min_score = int(args.split_score)
else:
    min_score = int(9)   
    
if args.disco_score:
    min_score = int(args.disco_score)
else:
    min_score = int(27)  

if args.split_weight:
    sw = float(args.split_weight)
else:
    sw = float(3)

if args.disco_weight:
    dw = float(args.disco_weight)
else:
    dw = float(1)
        
if not args.min_coverage:
    min_coverage = 2
else:
    min_coverage = float(args.min_coverage)
    
if not args.max_eval:
    max_eval = 0.01
else:
    max_eval = float(args.max_eval)
    
if not args.overlap_mask:
    overlap_mask = 100
else:
    overlap_mask = int(args.overlap_mask)

compliment_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}

if args.mge_file:
    run_mge = True
else:
    run_mge = False
    
if args.mapq_val:
    mapq_val = int(args.mapq_val)
else:
    mapq_val = 1
    
#set minimum reporting score for _realign.gff output 
if args.min_score:
    c_min_score = int(args.min_score)
else:
    c_min_score = 100

def pickle_loader(file_name, runmode):
    if runmode == 'dict':
        file = open(file_name,'rb')
        object_file = pickle.load(file)
        file.close()
    
    if runmode == 'df':
        print(file_name)
        object_file = pd.read_pickle(file_name)
    return(object_file)

def io_append(resource_dict, runmode='resource'):
    if runmode == 'resource':
        pre_existing_dict = {}
        
        resource_file_name = ('{}/resource.tab').format(final_output_dir)    
        resource_file = open(resource_file_name)
        
        for line in resource_file:
            cat, val = line.split('\t')
            pre_existing_dict[cat]=val
        resource_file.close()
    
        resource_file_name = ('{}/resource.tab').format(final_output_dir)    
        resource_file = open(resource_file_name, 'a')
        
        for cat, val in resource_dict.items():
            if cat not in pre_existing_dict:
                outline = ('{}\t{}\n').format(cat,val)
                resource_file.write(outline)
                
        resource_file.close()
        
        resource_pickle_name = ('{}/resource.p').format(final_output_dir)         
        with open(resource_pickle_name, 'wb') as file:
            pickle.dump(resource_dict, file)

    if runmode == 'ml_hypo':
        ml_pickle_name = ('{}/ml_hypo.p').format(temp_dir)         
        with open(ml_pickle_name, 'wb') as file:
            pickle.dump(ml_hypo_dict, file)
            
        ml_pickle_name = ('{}/ml_loci.p').format(temp_dir)         
        with open(ml_pickle_name, 'wb') as file:
            pickle.dump(ml_loci_dict, file)
        
def io_load(runmode='resource'):
    if runmode == 'resource':
        resource_file_name = ('{}/resource.p').format(final_output_dir)
        
        with open(resource_file_name, 'rb') as resource_file:
            resource_dict = pickle.load(resource_file)
            
        return(resource_dict) 
        
    if runmode == 'ml':
        ml_pickle_name = ('{}/ml_hypo.p').format(temp_dir) 
        with open(ml_pickle_name, 'rb') as ml_file:
            ml_hypo_dict = pickle.load(ml_file)

        ml_pickle_name = ('{}/ml_loci.p').format(temp_dir) 
        with open(ml_pickle_name, 'rb') as ml_file:
            ml_loci_dict = pickle.load(ml_file)
            
        return(ml_hypo_dict, ml_loci_dict) 
        
    if runmode == 'fq':
        pickle_in = ("{}/uniuid_to_seq_dict.p").format(pickles_dir)
        with open(pickle_in, 'rb') as fq_file:
            uniuid_to_seq_dict = pickle.load(fq_file)

        pickle_in = ("{}/uniuid_to_phred_dict.p").format(pickles_dir)
        with open(pickle_in, 'rb') as fq_file:
            uniuid_to_phred_dict = pickle.load(fq_file)

        return(uniuid_to_seq_dict, uniuid_to_phred_dict) 
    
def reverse_compliment(seq):
    seq = seq.upper()
    seq = seq[::-1]
    rc_seq = ''
    for nt in seq:
        if nt in compliment_dict:
            c_nt = compliment_dict[nt]
        else:
            c_nt = nt
            
        rc_seq += c_nt
        
    return(rc_seq)

def convert_sort(outfile_name):
    bash_command = ('samtools view -Sb {output_file}.sam > {output_file}_u.bam').format(output_file=outfile_name)
    subprocess.run([bash_command],stderr=subprocess.STDOUT,shell=True)

    bash_command = ('samtools sort -T tmp_sort -o {output_file}.bam {output_file}_u.bam').format(output_file=outfile_name)
    subprocess.run([bash_command],stderr=subprocess.STDOUT,shell=True)
        
    bash_command = ('samtools index {output_file}.bam').format(output_file=outfile_name)
    subprocess.run([bash_command],stderr=subprocess.STDOUT,shell=True)
        
    bash_command = ('rm {output_file}.sam {output_file}_u.bam').format(output_file=outfile_name)
    subprocess.run([bash_command],stderr=subprocess.STDOUT,shell=True)

    #bash_command = ('bedtools genomecov -ibam {output_file}.bam -d > {output_file}.depth').format(output_file=outfile_name)
    #subprocess.run([bash_command],stderr=subprocess.STDOUT,shell=True)
        
    return()

def parse_cigar(cigar, run_mode):
    """This function calculates the offset for the read based on the match
        For the 'first_match' run_mode it always takes the furthest left M, as this is the assigned start in sam file format
        if additional nucleotide types preceed the match these need to be stripped out. 
        For the 'last_match' run_mode it always takes the fullest length from the alignement.
        Such that it includes each M and each deletion 'D', it does not include insertions, hard or soft clips. 
    """
    match = 0
    if run_mode == 'first_match':        
        if 'M' in cigar:
            match_str = cigar.split('M')[0]
            for each_char in match_str:
                if each_char.isalpha():
                    match = int(match_str.rsplit(each_char,1)[1])
        return(int(match))
                    
    if run_mode == 'last_match':      
        if 'M' in cigar:
            temp_dict = {}
            val_str = ''

            for each_char in cigar:
                if each_char.isalpha():
                    if each_char in temp_dict:
                        pre_val = int(temp_dict[each_char])
                        temp_dict[each_char] = int(val_str) + pre_val
                        val_str = ''
                    else:
                        temp_dict[each_char] = int(val_str)
                        val_str = ''
                else:
                    val_str += str(each_char)
            
            match = temp_dict['M']
            if 'D' in temp_dict:
                match += temp_dict['D']
            
        return(match)
            
    if run_mode == 'non-matching':
        #count the matching nt                   
        for index in range(1,cigar.count('M')+1):
            temp = cigar.rsplit('M',index)[0]
            temp_list = re.split('M|S|D|I|H', temp)
            match += int(temp_list[-1])
            
        #count all the nts
        remsize = re.split('M|S|D|I|H', cigar)
        tote_size=0
        for step in remsize:
            if len(step)>0:
                tote_size+=int(step)
                
        return(tote_size-match)
        
    if run_mode == 'matching':                   
        for index in range(1,cigar.count('M')+1):
            temp = cigar.rsplit('M',index)[0]
            temp_list = re.split('M|S|D|I|H', temp)
            match += int(temp_list[-1])
                            
        return(match)
        
def make_fasta(temp_split_seq_file_name, hypo, qname_list):
    temp_split_seq_file = open(temp_split_seq_file_name, 'w')
        
    for qname in qname_list:
        if qname.count('_') == 1:
            qname = qname.split('_')[1]
        else:
            #split_NB501157:106:HKFYTBGX2:4:11504:6476:10198_2.2
            preid, midid, postid = qname.split('_')
            postid = postid.split('.')[1]
            qname = ('{}.{}').format(midid, postid)
            
        if qname in uniuid_to_seq_dict:
            seq = uniuid_to_seq_dict[qname]
            outline = ('>{}\n{}\n').format(qname,seq)
            temp_split_seq_file.write(outline) 
            
        else:    
            print(uniuid_to_seq_dict.keys())
            print(qname)
            1/0
            
    temp_split_seq_file.close()
    
def degzip(full_name, file_name, pre_dir=''):
    # a lot of this is to ensure that the original fastq is left undisturbed
    # it would be more elegant to just make a temp directory and copy it there
    
    if (full_name == file_name + '.gz') or (full_name == file_name) or (full_name == pre_dir+file_name):
        if '.gz' not in file_name:
            bashCommand = ('gunzip -c {} > {}').format(full_name, file_name)
            subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)
            print(bashCommand)
            
        if '.gz' in file_name:
            bashCommand = ('gunzip -c {} > {}').format(full_name, file_name.split('.gz')[0])
            subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)
            print(bashCommand)
        
    else:
        bashCommand = ('cp {} {}.gz').format(full_name, file_name)
        subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)
        print(bashCommand)
        
        bashCommand = ('gunzip -f {}').format(file_name)
        subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)
        print(bashCommand)
        
def parse_fastq(fastq_full_name, series_number, qname_lookup, uniuid_to_seq_dict, uniuid_to_phred_dict):
    
    outline = ('Parsing {}, using {}').format(fastq_full_name, series_number)
    print(outline)
    
    if '/' in fastq_full_name:
        pre_dir, fastq_file_name = handle_outfile(fastq_full_name)
        
        if fastq_full_name[-3:]=='.gz':
            #
            degzip(fastq_full_name, fastq_file_name, pre_dir)
            fastq_file_name = fastq_file_name[:-3]
        else:
            bashCommand = ('cp -f {} {}').format(fastq_full_name, fastq_file_name)
            subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)

    else:
        if fastq_full_name[-3:]=='.gz':
            fastq_file_name = fastq_full_name[:-3]
            degzip(fastq_full_name, fastq_file_name)
        else:
            fastq_file_name = fastq_full_name
    
    fq_file = open(fastq_file_name)
    line_ct = 10
    
    for line in fq_file:
        line = line.strip()
        line_ct += 1
        
        if line[0]=='@':
            inline_qname = line.split('@')[1].split(' ')[0]
            query_qname = line.split('@')[1].split(' ')[0] + series_number
            if (inline_qname in qname_lookup) and (query_qname not in uniuid_to_seq_dict):
                line_ct = 0
            else:
                line_ct = 10
                
        if line_ct==1:            
            if query_qname in uniuid_to_seq_dict:
                print('seq collision', query_qname, uniuid_to_seq_dict[query_qname])
                1/0
            else:
                uniuid_to_seq_dict[query_qname] = line

                    
        if line_ct==3:
            if query_qname in uniuid_to_phred_dict:
                print('phred collision', query_qname, uniuid_to_phred_dict[query_qname])
                1/0
            else:
                uniuid_to_phred_dict[query_qname] = line

    fq_file.close()
    return(uniuid_to_seq_dict, uniuid_to_phred_dict)
    
def build_fastq(infile_name, uniuid_to_seq_dict, uniuid_to_phred_dict):
    infile = open(infile_name)
    outfile = open(infile_name.split('.fa')[0] + '.fastq', 'w')
    
    for line in infile:
        if line[0] == '>':
            uid = line.split('>')[1].strip()
            
            if uid not in uniuid_to_seq_dict or uid not in uniuid_to_phred_dict:
                if uid not in uniuid_to_seq_dict:
                    print(uid)
                    print((uniuid_to_seq_dict))
                    print(len(uniuid_to_seq_dict))
                    1/0
                else:
                    print(uid)
                    print(uniuid_to_phred_dict)
                    print(len(uniuid_to_phred_dict))
                    1/0
            
            else:
                if uid in uniuid_to_seq_dict:
                    seq = uniuid_to_seq_dict[uid]
                
                if uid in uniuid_to_phred_dict:
                    phred = uniuid_to_phred_dict[uid]
                
                outline = ('@{}\n{}\n+\n{}\n').format(uid, seq, phred)
                outfile.write(outline)
            
    infile.close()
    outfile.close()

def prep_qname(qname, hypo, outline):
     
    #map each hypo to a read or set of reads
    if hypo in complete_qname_dict:
        complete_qname_dict[hypo]+=[qname]

    else:
        complete_qname_dict[hypo]=[qname]
    
    #map each read to a hypo
    if qname not in qname_lookup_dict:
        qname_lookup_dict[qname] = hypo
    
    #map each hypo to a location
    if hypo not in locus_lookup_dict:
        locus_lookup_dict[hypo] = outline

def parse_brks(break_file_name):
    print('\tParsing breakpoints')
    brks = open(break_file_name)
    
    for line in brks:
        if line[0]!='#':
            line = line.strip()
            
            hypo = line.split('\t')[0]
            chromo = line.split('\t')[1]
            start = int(line.split('\t')[2])
            stop = int(line.split('\t')[3])
            outline=('{},{},{}').format(chromo, start, stop)
            qname = line.split('\t')[8]
            
            if qname.count(',')>0:
                qname_list = qname.split(',')
                for qname in qname_list:
                    prep_qname(qname, hypo, outline)
            else:
                prep_qname(qname, hypo, outline)    
                
    brks.close()
        
def unpackbits(x,num_bits=12):
    xshape = list(x.shape)
    x = x.reshape([-1,1])
    to_and = 2**np.arange(num_bits).reshape([1,num_bits])
    upb = (x & to_and).astype(bool).astype(int).reshape(xshape + [num_bits])

    #0  (rp)    read_paired
    #1  (rmp)    read_mapped_in_proper_pair
    #2  (ru)    read_unmapped
    #3  (mu)    mate_unmapped
    #4  (rrs)    read_reverse_strand
    #5  (mrs)    mate_reverse_strand
    #6  (fip)    first_in_pair
    #7  (sip)    second_in_pair
    #8  (npa)    not_primary_alignment
    #9  (rfp)    read_fails_platform
    #10 (pcr)    read_is_PCR_or_optical_duplicate
    #11 (sa)    supplementary_alignment
    
    """ DISCORDANT definition (from samblaster)
        Both side of the read pair are mapped (neither FLAG 0x4 or 0x8 is set).
        The properly paired FLAG (0x2) is not set.
        Note: We implemented an additional criteria to distinguish between strand re-orientations and distance issues
        Strand Discordant reads must be both on the same strand.
    """
        
    """ SPLIT READS
        Identify reads that have between two and --maxSplitCount [2] primary and supplemental alignments.
        Sort these alignments by their strand-normalized position along the read.
        Two alignments are output as splitters if they are adjacent on the read, and meet these criteria:
            each covers at least --minNonOverlap [20] base pairs of the read that the other does not.
            the two alignments map to different reference sequences and/or strands. 
            the two alignments map to the same sequence and strand, and represent a SV that is at least --minIndelSize [50] in length, 
            and have at most --maxUnmappedBases [50] of un-aligned base pairs between them.
        Split read alignments that are part of a duplicate read will be output unless the -e option is used.
    """
    return(upb) 
    
def run_bedgraph(bam_file, each_sample):
    print('\tRunning bedtools genomecov on '+ str(infile_name) + '...')

    bashCommand = ('bedtools genomecov -bg -ibam {} > {}{}.bedgraph').format(bam_file, final_output_dir, each_sample)
    print(bashCommand)
    subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)

def run_mpileup(infile_name, each_sample, runmode):
    print('\tRunning mpileup genomecov on '+ str(infile_name) + '...')
     
    if runmode == 'RD':
        print('\tRunning samtools to determine sample depth...')
        if filter_by_bed:
            bashCommand = ('samtools mpileup -B -f {} {} -l {} -a -o {}/{}_{}.mpileup').format(fa_file, infile_name, filter_bed, pickles_dir, each_sample, runmode)        
        else:
            bashCommand = ('samtools mpileup -B -f {} {} -a -o {}/{}_{}.mpileup').format(fa_file, infile_name, pickles_dir, each_sample, runmode)                        
        print(bashCommand)
        subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)

    if runmode == 'discordant':
        print('\tRunning samtools to determine discordant read depth...')
        if filter_by_bed:
            bashCommand = ('samtools mpileup -B -f {} {} -l {} -a -A -o {}/{}_{}.mpileup').format(fa_file, infile_name, filter_bed, pickles_dir, each_sample, runmode)            
        else:
            bashCommand = ('samtools mpileup -B -f {} {} -a -A -o {}/{}_{}.mpileup').format(fa_file, infile_name, pickles_dir, each_sample, runmode)
        print(bashCommand)
        subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)
        
    if runmode == 'split':
        print('\tRunning samtools to determine split read depth...')
        if filter_by_bed:
            bashCommand = ('samtools mpileup -B -f {} {} -l {} -a -A -o {}/{}_{}.mpileup').format(fa_file, infile_name, filter_bed, pickles_dir, each_sample, runmode)     
        else:
            bashCommand = ('samtools mpileup -B -f {} {} -a -A -o {}/{}_{}.mpileup').format(fa_file, infile_name, pickles_dir, each_sample, runmode)     
        print(bashCommand)            
        subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)
        

def populate_filter_dict(chromo, start, stop, region):    
    if chromo not in filter_region_dict:
        filter_region_dict[chromo]=set()

    if chromo in filter_region_dict:
        for index in range(start,stop+1):
            filter_region_dict[chromo].add(index)  

    if region not in region_name_dict:
        region_name_dict[region]=0 
        
    region_name_dict[region]+=1
    

def parse_filter_gff(gff_file):
    parse_filter_object = open(gff_file)
    for line in parse_filter_object:
        if line[0]!='#':
            region = 'unannotated'
            if 'ID=' in line:
                region = line.split('ID=')[1].split(';')[0].strip()
            else:
                if line.count('\t') > 7:
                    region = line.split('\t')[8].split(';')[0].strip()
                
            chromo = line.split('\t')[0]
            start = int(line.split('\t')[3])
            stop = int(line.split('\t')[4])
            populate_filter_dict(chromo, start, stop, region)
            
    parse_filter_object.close()
    
def parse_tab(tab_file, anc_QC_list):
    anc_tab = open(tab_file)
    
    for line in anc_tab:
        d_ct = 0
        s_ct = 0
        chromo = line.split('\t')[1]
        start = int(line.split('\t')[2])
        stop = int(line.split('\t')[3])
        breakpoint_name = line.split('\t')[0]
        
        uid_str = line.split('\t')[8].split(',')
        s_ct = uid_str.count('split_')
        d_ct = uid_str.count('discordant_')

        score = (sw*s_ct) + (dw*d_ct)
        
        if (s_ct >= 1) and (score >= min_score):
            anc_QC_list.append(breakpoint_name)
            populate_filter_dict(chromo, start, stop, breakpoint_name)
        
        if (s_ct == 0) and (score >= 3*min_score):
            anc_QC_list.append(breakpoint_name)
            populate_filter_dict(chromo, start, stop, breakpoint_name)
            
    anc_tab.close()
    
    return(anc_QC_list)
            
def parse_bed(bed_file, source):
    in_bed = open(bed_file)
    
    for line in in_bed:
        chromo = line.split('\t')[0]
        start = int(line.split('\t')[1])
        stop = int(line.split('\t')[2])
        breakpoint_name = line.split('\t')[3] + '_' + source
        
        #if breakpoint_name not in anc_QC_list:
        populate_filter_dict(chromo, start, stop, breakpoint_name)
    
    in_bed.close()

def check_coverage(source, prefilter_object_name, contig_seq_dict):      
    process = False

    temp_filter_object = open(prefilter_object_name)
    
    print(prefilter_object_name)
    seq_set = set()
    
    for line in temp_filter_object:
        if line[0] == '>' and not process:
            process=True
            #node_name = source
                    
        if line[0]!= '>' and process:
            seq_set.add(line.strip())
            process=False
            
    temp_filter_object.close()

    if len(seq_set)>0:
        seq_str = ''
        for seq in seq_set:
            seq_str+=seq
            
        contig_seq_dict[source]=seq_str
            
    return(contig_seq_dict)    

def long_walk(nt, long_range, ct_dict):
    miss_ct = 0
    depth = []
    
    gap = 10

    for s_step in long_range:
        if s_step in ct_dict:
            depth.append(float(ct_dict[s_step]))
            nt=s_step
            miss_ct = 0
            
        else:
            miss_ct += 1

        if miss_ct >= 3*gap:
            return(nt,depth)
            
    return(nt,depth)

def local_walk(short_range, ct_dict):
    hit_ct = 0

    for s_step in short_range:
        if s_step in ct_dict:
            hit_ct += 1
            
        if hit_ct >= (purity*min_length):
            return(True)
    else:
        return(False)
    
def build_mge_lookup_list(chromo):
    global mge_nt_lookup_dict
    
    mge_nt_list = []
    
    sub_mge_lookup = mge_nt_lookup_dict[chromo]
    
    for nt in sub_mge_lookup:
        mge_nt_list.append(nt)

    return(mge_nt_list)  

def build_mge_region_lookup_list(chromo):
    global mge_nt_lookup_dict
    
    mge_nt_list = []
    
    sub_mge_lookup = mge_nt_lookup_dict[chromo]
    
    for nt in sub_mge_lookup:
        mge_nt_list.append(nt)

    return(mge_nt_list)
    
def build_mge_trace(subz_df, suby_df, each_type, c_median, g_median, mge_g_std, chromo, mge_gff):
    global mge_nt_lookup_dict
    
    sub_mge_lookup = mge_nt_lookup_dict[chromo]
    
    subz_dict = subz_df.to_dict(orient='index')
    subz_ct_dict = {}
    for index, deets in subz_dict.items(): 
        subz_ct_dict[int(deets['nuc'])]=float(deets['ct'])
        
    suby_dict = suby_df.to_dict(orient='index')
    suby_ct_dict = {}
    for index, deets in suby_dict.items():
        suby_ct_dict[int(deets['nuc'])]=float(deets['ct'])
    
    seed_set = set()
    
    if len(subz_ct_dict) > 1:
        min_nuc = min(subz_ct_dict)
        max_nuc = max(subz_ct_dict)                            

        for nuc, ct in subz_ct_dict.items():   
            if nuc in sub_mge_lookup:
                if (nuc not in seed_set):
                    short_range = range(nuc-int(min_length/2),nuc+int(min_length/2))
                    if (local_walk(short_range, subz_ct_dict)) or (each_type != 'RD'):                        
                        left_long_range = range(min_nuc-int(min_length/2), nuc)[::-1]
                        right_long_range = range(nuc, max_nuc+int(min_length/2))
                        
                        start, depth = long_walk(nuc, left_long_range, suby_ct_dict)          
                        stop, depth_2 = long_walk(nuc, right_long_range, suby_ct_dict)
            
                        depth += depth_2
                        
                        depth_median = np.mean(depth)
                        depth_std = np.std(depth)
    
                        for nt in range(start, stop+1):
                            seed_set.add(nt)
                                                        
                        if (start >= 0) and (stop >= 0) and (depth_median >= (g_median + mge_g_std)):
                            rel_c = depth_median/float(c_median)
                            rel_g = depth_median/float(g_median)
                                                        
                            outline = ('{chromo}\t{each_type}_depth\tMGE_CNV'
                                       '\t{start}\t{stop}\t{sum_d}\t.\t.'
                                       '\tID={chromo}:{start}-{stop};rel_chromosome_RD={rel_c};rel_genome_RD={rel_g};sample={sample}\n'
                                       ).format(
                                           chromo = every_chromo, 
                                           each_type = each_type, 
                                           start = start, stop=stop, 
                                           sum_d = np.mean(depth),
                                           rel_c = round(rel_c,2), 
                                           rel_g = round(rel_g,2), 
                                           sample = each_sample)
                            
                            mge_gff.write(outline)
                                                        
                            locus = ('{}:{}-{}').format(every_chromo, start, stop)
                            
                            depth_dict[locus] = [rel_c, rel_g, depth_std]
        
def build_trace(subz_df, suby_df, each_type, c_median, g_median):
    
    subz_dict = subz_df.to_dict(orient='index')
    subz_ct_dict = {}
    for index, deets in subz_dict.items(): 
        subz_ct_dict[int(deets['nuc'])]=float(deets['ct'])
        
    suby_dict = suby_df.to_dict(orient='index')
    suby_ct_dict = {}
    for index, deets in suby_dict.items():
        suby_ct_dict[int(deets['nuc'])]=float(deets['ct'])
    
    seed_set = set()
    
    if len(subz_ct_dict) > 1:
        min_nuc = min(subz_ct_dict)
        max_nuc = max(subz_ct_dict)                            

        for nuc, ct in subz_ct_dict.items():   
            
            if (nuc not in seed_set):
                short_range = range(nuc-int(min_length/2),nuc+int(min_length/2))
                
                if (local_walk(short_range, subz_ct_dict)) or (each_type != 'RD'):                    
                    left_long_range = range(min_nuc-int(min_length/2), nuc)[::-1]
                    right_long_range = range(nuc, max_nuc+int(min_length/2))
                    
                    start, depth = long_walk(nuc, left_long_range, suby_ct_dict)          
                    stop, depth_2 = long_walk(nuc, right_long_range, suby_ct_dict)
        
                    depth += depth_2
                    
                    depth_median = np.mean(depth)
                    depth_std = np.std(depth)

                    for nt in range(start, stop+1):
                        seed_set.add(nt)
                                                
                    if (start >= 0) and (stop >= 0):
                        rel_c = depth_median/float(c_median)
                        rel_g = depth_median/float(g_median)
                        
                        outline = ('{chromo}\t{each_type}_depth\tCNV'
                                   '\t{start}\t{stop}\t{sum_d}\t.\t.'
                                   '\tID={chromo}:{start}-{stop};rel_chromosome_RD={rel_c};rel_genome_RD={rel_g};sample={sample}\n'
                                   ).format(
                                       chromo = every_chromo, 
                                       each_type = each_type, 
                                       start = start, stop=stop, 
                                       sum_d = np.mean(depth),
                                       rel_c = round(rel_c,2), 
                                       rel_g = round(rel_g,2), 
                                       sample = each_sample)
                        
                        depth_gff_outfile.write(outline)
                        
                        locus = ('{}:{}-{}').format(every_chromo, start, stop)
                        
                        depth_dict[locus] = [rel_c, rel_g, depth_std]
                            
def get_seq(seq_file_name):
    seq_file = open(seq_file_name)
    
    #seq_file = open('C:/Gresham/tiny_projects/tiny_erisapfel/temp/test.txt')
    seq_dict = {}
    seq = False
        
    for line in seq_file:
        if line[0] == '>':
            if seq:
                if seq not in seq_dict:
                    seq_dict[seq] = 0
                
                seq_dict[seq] += 1
            
            seq = ''
        
        if line[0] != '>':
            seq += line.strip()
            
    seq_file.close()
    
    if seq not in seq_dict:
        seq_dict[seq] = 0
    
    seq_dict[seq] += 1
    
    if len(seq_dict) > 1:
        print(seq_file_name)
        print(seq_dict)
        1/0
                
    max_seq = ''
    max_score = 0
    for seq in seq_dict:
        seq_score = seq_dict[seq]
        
        if seq_score == max_score:
            if len(seq) > len(max_seq):
                max_seq = seq
                        
        if seq_score > max_score:
            max_seq = seq
            max_score = seq_score
    
    #print('max_seq', max_seq)
    
    return(max_seq)


def return_most_unique_read(read_nt_dict):
    reverse_read_nt_dict = {}
    for nt in read_nt_dict:
        numb_set = read_nt_dict[nt]
        
        for numb in numb_set:
            if numb not in reverse_read_nt_dict:
                reverse_read_nt_dict[numb] = []
            reverse_read_nt_dict[numb].append(len(numb_set))
    
    max_read_unique = False
    for numb in reverse_read_nt_dict:
        max_unique = 0
        unique_map = 0
        
        for nt in reverse_read_nt_dict[numb]:
            if nt == 1:
                unique_map += 1
        
        if (unique_map >= 15) and (unique_map > max_unique):
            max_unique = unique_map
            max_read_unique = numb
        
    return(max_read_unique)

def build_otherside_dict(query_deets_dict, split_assignment, anchor, otherside_type, otherside_dict):
    for next_anchor in split_assignment[otherside_type]:
        if next_anchor != anchor:
            chromo = query_deets_dict[next_anchor]['chromo']
            next_bit_score = query_deets_dict[next_anchor]['bit_score']
            next_hf = min(query_deets_dict[next_anchor]['hit_from'], query_deets_dict[next_anchor]['hit_to'])
            next_ht = max(query_deets_dict[next_anchor]['hit_from'], query_deets_dict[next_anchor]['hit_to'])

            if chromo not in otherside_dict:
                otherside_dict[chromo] = {}

            for next_nt in range(next_hf, next_ht+1):
                if next_nt not in otherside_dict[chromo]:
                    otherside_dict[chromo][next_nt] = 0
                    
                otherside_dict[chromo][next_nt]+=next_bit_score
                
    return(otherside_dict)

def parse_json(json_file_name, seq, hypothesis_dict, anchor_contig_dict):
    '''
    

    Parameters
    ----------
    json_file_name : str
        json_file_name points to the json file output from the blastn results
    contig_seq_dict : dict
        DESCRIPTION.
    gff_list_dict : TYPE
        DESCRIPTION.
    gff_rank_dict : TYPE
        DESCRIPTION.
    blast_read_aligned_sections_dict : TYPE
        DESCRIPTION.
    blast_genome_aligned_regions_dict : TYPE
        DESCRIPTION.
    gff_uid_dict : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    
    active = False
    query_deets_dict = {}
    try:
        data = json.load(open(json_file_name))
    except:
        print('json file ' + json_file_name + ' not found')
    
    for report_index in range(len(data["BlastOutput2"])):
        data_dict = (data["BlastOutput2"][report_index])
        for each_report in data_dict.items():
            for a_key, a_value in enumerate(each_report):
                if type(a_value)==dict:                   
                    for b_key, b_value in a_value.items():
                        if type(b_value)==dict:
                            for c_key, c_value in b_value.items():
                                if ('bl2seq' in c_key) and (type(c_value)==list):
                                    hit_dict = c_value[0]
                                    
                                    for d_key, d_value in hit_dict.items():                                        
                                        q_title = str(hit_dict['query_title'])
                                          
                                        if ('hits' in d_key) and (type(d_value)==list) and (len(d_value)>0):
                                            for each_hits in d_value:
                                                
                                                for e_key, e_value in each_hits.items():
                                                
                                                    base = q_title + '.'+str(each_hits['num']) 
                                                    chromo = each_hits['description']
                                                    chromo = chromo[0]
                                                    #chromo_dict[base] = str(chromo['id'])
                                                        
                                                    if (e_key == 'hsps') and (type(e_value)==list):
                                                        for e_index in range(len(e_value)):
                                                            each_hsps = e_value[e_index]
    
                                                            numb = str(base)+'.'+str(each_hsps['num'])
                                                            
                                                            if len(numb)>1:
                                                                active = True
                                                                hit_from = int(each_hsps["hit_from"])                                                                
                                                                hit_to = int(each_hsps["hit_to"])                                                                
                                                                query_from = int(each_hsps["query_from"])                                                                
                                                                query_to = int(each_hsps["query_to"])                                                            
                                                                bit_score = float(each_hsps["bit_score"])
                                                                evalue_score = float(each_hsps["evalue"])
                                                                query_strand = str(each_hsps["query_strand"])                                                                
                                                                hit_strand = str(each_hsps["hit_strand"])                                                                
                                                                qseq = str(each_hsps["qseq"])                                                            
                                                                hseq = str(each_hsps["hseq"])
                                                            
                                                            if evalue_score > max_eval:
                                                                active = False
                                                            
                                                            if active:
                                                                active = False
                                                                query_deets_dict[numb] = {
                                                                    'q_id':base,
                                                                    'chromo':chromo['id'],
                                                                    'hit_from':hit_from,
                                                                    'hit_to':hit_to,
                                                                    'query_from':query_from,
                                                                    'query_to':query_to,
                                                                    'bit_score':bit_score,
                                                                    'query_strand':query_strand,
                                                                    'hit_strand':hit_strand,
                                                                    'qseq':qseq, 
                                                                    'hseq':hseq,
                                                                    'q_title':q_title,
                                                                    'contig':seq}
                                                                
                                                                numb = 0

    '''
    for each nt of the read (query), determine if it is part of the anchor or not
    '''
    read_nt_dict = {}
    
    
    for numb in query_deets_dict:
        qf = query_deets_dict[numb]['query_from']
        qt = query_deets_dict[numb]['query_to']
        
        for nt in range(qf, qt+1):
            if nt not in read_nt_dict:
                read_nt_dict[nt] = set()
                
            read_nt_dict[nt].add(numb)

    split_assignment = {}
    
    for numb in query_deets_dict:
        unique_map = 0
        multi_map = 0
        qf = query_deets_dict[numb]['query_from']
        qt = query_deets_dict[numb]['query_to']
        
        for nt in range(qf, qt+1):
            if len(read_nt_dict[nt]) > 1:
                multi_map += 1
            else:
                unique_map += 1
                
        if (unique_map/(unique_map+multi_map)) > 0.8:
            if 'anchor' not in split_assignment:
                split_assignment['anchor'] = set()
                
            split_assignment['anchor'].add(numb)
            
        else:
            if 'multimap' not in split_assignment:
                split_assignment['multimap'] = set()
                
            split_assignment['multimap'].add(numb)
            
    
    if 'anchor' not in split_assignment:
        most_unique_read = return_most_unique_read(read_nt_dict)
        
        if most_unique_read:
            split_assignment['anchor'] = set()
            split_assignment['anchor'].add(most_unique_read)

    if 'anchor' in split_assignment:        
        #anchors first
        for anchor in split_assignment['anchor']:
            
            #build anchor side:
            chromo = query_deets_dict[anchor]['chromo']
            hf = min(query_deets_dict[anchor]['hit_from'], query_deets_dict[anchor]['hit_to'])
            ht = max(query_deets_dict[anchor]['hit_from'], query_deets_dict[anchor]['hit_to'])
            
            if chromo not in hypothesis_dict:
                hypothesis_dict[chromo] = {}
                
            
            
            #build otherside
            for nt in range(hf, ht+1):
                if nt not in hypothesis_dict[chromo]:
                    otherside_dict = {}
                else:
                    otherside_dict = hypothesis_dict[chromo][nt]
                    
                for otherside_type in split_assignment:
                    #print(otherside_type)
                    otherside_dict = build_otherside_dict(query_deets_dict, split_assignment, anchor, otherside_type, otherside_dict)
                    
                hypothesis_dict[chromo][nt] = otherside_dict
                    
            if chromo not in anchor_contig_dict:
                anchor_contig_dict[chromo] = {}
                
                for nt in range(hf, ht+1):
                    if nt not in anchor_contig_dict[chromo]:
                        anchor_contig_dict[chromo][nt] = set()
                        
                    anchor_contig_dict[chromo][nt].add(seq)
            
    return(hypothesis_dict, anchor_contig_dict)

def make_anchor_regions(nt_list, gap):
    nt_list.sort()

    region_dict = {}

    new_region = True
    
    for nt in nt_list:
        if new_region:
            region_number = len(region_dict)
            region_dict[len(region_dict)] = {'start':nt, 'stop':0}
            new_region = False
        
        if nt + gap not in nt_list:
            region_dict[region_number]['stop'] = nt
            new_region = True
    
    return(region_dict)

def make_otherside_regions(chromo, region_dict, chromo_hypothesis_dict, chromo_anchor_contig_dict, otherside_dict, gap):

    uid_ct = 0

    for region_number in region_dict:
        start = region_dict[region_number]['start']
        stop = region_dict[region_number]['stop']
        
        contig_set = set()
        for nt in range(start, stop+1):
            if nt in chromo_anchor_contig_dict:
                for contig in chromo_anchor_contig_dict[nt]:
                    contig_set.add(contig)
        
        for nt in range(start, stop+1):
            if nt in chromo_hypothesis_dict:                
                for otherside_chromo in chromo_hypothesis_dict[nt]:
                    next_nt_list = []
                    for next_nt in chromo_hypothesis_dict[nt][otherside_chromo]:
                        next_nt_list.append(next_nt)
                        
                    next_region_dict = make_anchor_regions(next_nt_list, gap)
                    
                    for next_region in next_region_dict:                       
                        next_bs_list = []
                        next_start = next_region_dict[next_region]['start']
                        next_stop = next_region_dict[next_region]['stop']
                        
                        for next_nt in range(next_start, next_stop+1):
                            #because of gaps not all nt may be present for scoring, these will be counted as a 0
                            if next_nt in chromo_hypothesis_dict[nt][otherside_chromo]:
                                next_bs_list.append(chromo_hypothesis_dict[nt][otherside_chromo][next_nt])
                            else:
                                next_bs_list.append(0)
                        
                        score = round(np.median(next_bs_list))
                        
                        
                        if chromo not in otherside_dict:
                            otherside_dict[chromo] = {}
                        
                        otherside_dict[chromo][uid_ct] = {"anchor_chromo":chromo, 
                                                               "anchor_start":start,
                                                               "anchor_stop":stop,
                                                               "other_chromo":otherside_chromo,
                                                               "other_start":next_start,
                                                               "other_stop":next_stop,
                                                               "score":score,
                                                               "contig":contig_set}
                        
                        uid_ct += 1
                                
    return(otherside_dict) 
    
        
'''
The key premise here is that there are two types of reads 
1. will be a read with unique matches on both sides (2 anchors), this is a trivial case
2. a read with only one unique match (1 anchors) one variable.

In the second case it can be informative to ask two questions:
    a. what are all the variable loci the anchor maps to?
    b. what is the relative read depth of those loci?
    
When a read is aligned the anchor(s) nucelotides (both read and reference) should be identified as should the read depth weighted variable.
The anchor and variables should then be read out. 
''' 

def fourway_test(uid_start, next_start, uid_stop, next_stop, gap):
    # two edges
    if ((abs(uid_start - next_start) <= gap) or 
        (abs(uid_stop - next_stop) <= gap) or
        (abs(uid_start - next_stop) <= gap) or
        (abs(uid_stop - next_start) <= gap)):
        return(True)
    
    # complete overlap
    if ((uid_start <= next_start) and 
        (uid_stop >= next_stop)):
        return(True)
    
    if ((uid_start >= next_start) and 
        (uid_stop <= next_stop)):
        return(True)
    
    # body overlap
    if ((uid_start <= next_start) and 
        (uid_stop >= next_start)):
        return(True)

    if ((next_start <= uid_start) and 
        (next_stop >= uid_start)):
        return(True)
    
    if ((uid_stop <= next_stop) and 
        (uid_stop >= next_start)):
        return(True)

    if ((next_start <= uid_start) and 
        (next_stop >= uid_start)):
        return(True)

    return(False)

def collapse_other(chromoside_dict):
    anchor_ct = 0
    
    to_add_dict = {}
    
    for uid in chromoside_dict:
        anchor_ct += 1
        
        for next_uid in chromoside_dict:                  
            if (uid != next_uid):
                uid_chromo = chromoside_dict[uid]['anchor_chromo']
                next_chromo = chromoside_dict[next_uid]['anchor_chromo']
                uid_o_chromo = chromoside_dict[uid]['other_chromo']
                next_o_chromo = chromoside_dict[next_uid]['other_chromo']
                                
                if (uid_chromo == next_chromo) and (uid_o_chromo == next_o_chromo):
                    uid_start=chromoside_dict[uid]['anchor_start']
                    uid_stop=chromoside_dict[uid]['anchor_stop']
                    next_start=chromoside_dict[next_uid]['anchor_start']
                    next_stop=chromoside_dict[next_uid]['anchor_stop']
                    
                    uid_o_start=chromoside_dict[uid]['other_start']
                    uid_o_stop=chromoside_dict[uid]['other_stop']
                    next_o_start=chromoside_dict[next_uid]['other_start']
                    next_o_stop=chromoside_dict[next_uid]['other_stop'] 
                    
                    pass_both = False
                    
                    if fourway_test(uid_start, next_start, uid_stop, next_stop, gap):    
                        if fourway_test(uid_o_start, next_o_start, uid_o_stop, next_o_stop, gap):
                            pass_both = True
                        
                    if fourway_test(uid_start, next_o_start, uid_stop, next_o_start, gap): 
                        if fourway_test(uid_o_start, next_start, uid_o_stop, next_stop, gap):
                            pass_both = True
                        
                    if pass_both:
                        new_score = max(chromoside_dict[uid]['score'], chromoside_dict[next_uid]['score'])
                        
                        print('pass 4')
                        
                        print(uid, next_uid)
                        
                        print(chromoside_dict[uid])
                        print(chromoside_dict[next_uid])
                                                    
                        new_contig_set = chromoside_dict[uid]['contig']
                        for contig in chromoside_dict[next_uid]['contig']:
                            new_contig_set.add(contig)
                                                        
                        to_add_dict = {"anchor_chromo":uid_chromo, 
                                        "anchor_start":min(uid_start, next_start),
                                        "anchor_stop":max(uid_stop, next_stop),
                                        "other_chromo":uid_o_chromo,
                                        "other_start":min(uid_o_start, next_o_start),
                                        "other_stop":max(uid_o_stop, next_o_stop),
                                        "score":new_score,
                                        "contig":new_contig_set}
                        
                        print('to_add_dict', to_add_dict)
                        
                        return(True, uid, next_uid, to_add_dict)
                        
    return(False, 0, 0, to_add_dict)

def collapse_hypotheses(otherside_dict, gap):

    original_size = 0
    chromo_set = set()
    
    for chromo in otherside_dict:
        chromo_set.add(chromo)
        
        original_size += len(otherside_dict[chromo])
    
    new_uid = original_size+1
    reduced = True
    
    for chromo in chromo_set:
        og_size = len(otherside_dict[chromo])
        
        ts = time.ctime()
        print('reduced ',ts)
        reduced = True
        
        print(chromo)
        
        while reduced:
            print('starting again')
            chromoside_dict = otherside_dict[chromo]
            print(len(otherside_dict))
            
            size = len(chromoside_dict)
            if size % 1000 == 0:
                ts = time.ctime()
                print('reduced ',ts)
                outline = ('{} size, {}, {}, {}').format(chromo, size, og_size, round(size/og_size,2))
                print(outline)
            
            reduced, uid, next_uid, to_add_dict = collapse_other(chromoside_dict)
            print('now outside')
            
            if reduced:
                if new_uid in otherside_dict[chromo]:
                    print('uid collision error')
                    1/0 
                    
                otherside_dict[chromo][new_uid] = to_add_dict
            
                del otherside_dict[chromo][uid]
                    
                del otherside_dict[chromo][next_uid]
                    
                new_uid+=1
                
            print('is it reduced')
            print(len(otherside_dict))
            print()
            
    return(otherside_dict)
                              
def calculate_global_disco_stats(each_sample):
    global_disco_stats_dict ={
        'rd_df':0, 'rd_global_median': 0, 'rd_global_std': 0,
        'disco_df':0, 'disco_global_median': 0, 'disco_global_std': 0,
        'split_df':0, 'split_global_median': 0, 'split_global_std': 0}
    
    pickle_name = ('{}/mpileup_{}_RD.p').format(pickles_dir, each_sample)
    rd_df = pickle_loader(pickle_name, 'df')
    
    fzero = rd_df.replace(0, np.NaN)
    rd_global_median = fzero["ct"].median()
    rd_global_std = fzero["ct"].std()
    
    global_disco_stats_dict['rd_df'] = rd_df
    global_disco_stats_dict['rd_global_median'] = rd_global_median
    global_disco_stats_dict['rd_global_std'] = rd_global_std
    
    pickle_name = ('{}/mpileup_{}_discordant.p').format(pickles_dir, each_sample)
    disco_df = pickle_loader(pickle_name, 'df')
    
    fzero = disco_df.replace(0, np.NaN)
    disco_global_median = fzero["ct"].median()
    disco_global_std = fzero["ct"].std()
    
    global_disco_stats_dict['disco_df'] = disco_df
    global_disco_stats_dict['disco_global_median'] = disco_global_median
    global_disco_stats_dict['disco_global_std'] = disco_global_std
    
    pickle_name = ('{}/mpileup_{}_split.p').format(pickles_dir, each_sample)
    split_df = pickle_loader(pickle_name, 'df')
    
    fzero = split_df.replace(0, np.NaN)
    split_global_median = fzero["ct"].median()
    split_global_std = fzero["ct"].std()
           
    global_disco_stats_dict['split_df'] = split_df
    global_disco_stats_dict['split_global_median'] = split_global_median
    global_disco_stats_dict['split_global_std'] = split_global_std
    
    return(global_disco_stats_dict)

def collect_stats(chromo, start, stop, mpileup_df, runmode, score_store):
    anchor_store = ('{}:{}-{}').format(
        chromo, start, stop)
    
    if runmode == 'rd':
        type_df = mpileup_df['rd_df']
        global_median = mpileup_df['rd_global_median']
        global_std = mpileup_df['rd_global_std']
        
        if anchor_store in score_store[runmode]:
            sub_df = score_store[runmode][anchor_store] 
        else:
            sub_df = type_df.loc[(type_df["chromo"] == chromo) & 
               (type_df["nuc"] >= start) & 
               (type_df["nuc"] <= stop)]
            
            score_store[runmode][anchor_store] = sub_df
    
    if runmode == 'split':
        type_df = mpileup_df['split_df']
        global_median = mpileup_df['split_global_median']
        global_std = mpileup_df['split_global_std']
        
        if anchor_store in score_store[runmode]:
            sub_df = score_store[runmode][anchor_store] 
        else:
            sub_df = type_df.loc[(type_df["chromo"] == chromo) & 
               (type_df["nuc"] >= start) & 
               (type_df["nuc"] <= stop)]
            
            score_store[runmode][anchor_store] = sub_df
        
    if runmode == 'disco':
        type_df = mpileup_df['disco_df']
        global_median = mpileup_df['disco_global_median']
        global_std = mpileup_df['disco_global_std']  
            
        if anchor_store in score_store[runmode]:
            sub_df = score_store[runmode][anchor_store] 
        else:
            sub_df = type_df.loc[(type_df["chromo"] == chromo) & 
               (type_df["nuc"] >= start) & 
               (type_df["nuc"] <= stop)]
            
            score_store[runmode][anchor_store] = sub_df
                            
    sub_median = sub_df["ct"].median()
    sub_sum = sub_df["ct"].sum()
    
    if sub_sum > 3:

        if global_std == 0:
            global_std = 1
            
        difference_ratio = (sub_median-global_median)/global_std
    
        return(difference_ratio, score_store)
    
    return(0, score_store)

def disco_depth_stats(uid_deets, global_disco_stats_dict, score, score_store):        
    
    anchor_chromo = uid_deets['anchor_chromo']
    anchor_start = uid_deets['anchor_start']
    anchor_stop = uid_deets['anchor_stop']
    
    other_chromo = uid_deets['other_chromo']
    other_start = uid_deets['other_start']
    other_stop = uid_deets['other_stop']
    
    rd_anchor_term, score_store = collect_stats(
        anchor_chromo, anchor_start, 
        anchor_stop, global_disco_stats_dict, "rd", score_store)
                
    rd_other_term, score_store = collect_stats(
        other_chromo, other_start, 
        other_stop, global_disco_stats_dict, "rd", score_store)
    
    split_anchor_term, score_store = collect_stats(
        anchor_chromo, anchor_start, 
        anchor_stop, global_disco_stats_dict, "split", score_store)
            
    split_other_term, score_store = collect_stats(
        other_chromo, other_start, 
        other_stop, global_disco_stats_dict, "split", score_store)
            
    disco_anchor_term, score_store = collect_stats(
        anchor_chromo, anchor_start, 
        anchor_stop, global_disco_stats_dict, "disco", score_store)
            
    disco_other_term, score_store = collect_stats(
        other_chromo, other_start, 
        other_stop, global_disco_stats_dict, "disco", score_store)
    
    relative_score = (rd_anchor_term * score * rd_other_term + 
        split_anchor_term * score* split_other_term +
        disco_anchor_term * score* disco_other_term)
    
    return(relative_score, score_store)
    
def summarize_hypotheses(hypothesis_dict, anchor_contig_dict, gap, resource_dict):
    '''
    Summarize hypotheses 
    '''    
    fa_file = resource_dict['genome_fa']
    bam_file = resource_dict['bam_file']
    
    chromo_size_p = resource_dict['chromo_size']
    with open(chromo_size_p, 'rb') as fp:
        chromo_size_dict = pickle.load(fp)
        
    chromo_list = list(chromo_size_dict.keys())
    chromo_list.sort()
    
    chromo_size_str = ''
    
    for chromo in chromo_list:
        chromo_size = chromo_size_dict[chromo]
        
        chromo_size = chromo_size.strip()
        
        outline = ('##contig=<ID={chromo},length={chromo_size}>\n').format(
            chromo = chromo, chromo_size = chromo_size)
        
        chromo_size_str+=outline
    
    
    #vcf format requires a static header
    today = date.today().strftime("%d_%m_%Y")
    vcf_header = ('##fileformat=VCFv4.2\n'
                  '##fileDate={today}\n'
                  '##source=erisapfelv1.4\n'
                  '##reference={fa_file}\n'
                  '{contig}'
                  '##ALT=<ID=SV,Description="Structural Variant">\n'
                  '##INFO=<ID=AS,Number=1,Type=Float,Description="Absolute Score">\n'
                  '##INFO=<ID=RF,Number=1,Type=Float,Description="Relative Score">\n'
                  '##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">\n'
                  '##INFO=<ID=OTHERSIDE,Number=1,Type=String,Description="Coordinate of otherside of breakpoint">\n'
                  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                  '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\n'
                  '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n'
                  '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{bam_file}\n').format(
                      today=today, fa_file=fa_file, bam_file=bam_file, contig = chromo_size_str)
                      
    gff_header = ('###fileDate={today}\n'
                  '##source=erisapfelv1.4\n'
                  '##reference={fa_file}\n').format(
                      today=today, fa_file=fa_file)
    
    #while the vcf records are still handled by the vcf_set
    vcf_set = set()
    #while the vcf records are still handled by the vcf_set
    gff_set = set()                  
    otherside_dict = {}
        
    for chromo in hypothesis_dict:
                            
        nt_list = []
        for nt in hypothesis_dict[chromo]:
            nt_list.append(nt)
            
        region_dict = make_anchor_regions(nt_list, gap)
        
        otherside_dict = make_otherside_regions(chromo, region_dict, hypothesis_dict[chromo], anchor_contig_dict[chromo], otherside_dict, gap)
        
        original_size = 0
        for chromo in otherside_dict:
            original_size += len(otherside_dict[chromo])   
    
    otherside_dict = collapse_hypotheses(otherside_dict, gap)
    
    changed_size = 0
    for chromo in otherside_dict:
        changed_size += len(otherside_dict[chromo])
    
    resource_dict = io_load()
        
    each_sample = resource_dict['run_name']
                
    global_disco_stats_dict = calculate_global_disco_stats(each_sample)
    
    score_store = {'rd':{}, 'split':{}, 'disco':{}}
    top_performer_dict = {}
    
    for chromo in otherside_dict:
        chromo_size = len(otherside_dict[chromo])
        pissct = 0
        for uid in otherside_dict[chromo]:
            # print('otherside_dict[chromo][uid]')
            # print(otherside_dict[chromo][uid])
            
            outline = ('chromo: {}, uid: {}, ct: {} out of {}').format(
                chromo, uid, pissct, chromo_size)
            pissct+=1
                        
            anchor_chromo = otherside_dict[chromo][uid]['anchor_chromo']
            anchor_start = otherside_dict[chromo][uid]['anchor_start']
            anchor_stop = otherside_dict[chromo][uid]['anchor_stop']
            other_chromo = otherside_dict[chromo][uid]['other_chromo']
            other_start = otherside_dict[chromo][uid]['other_start']
            other_stop = otherside_dict[chromo][uid]['other_stop']
            score = otherside_dict[chromo][uid]['score']
            contig_set = otherside_dict[chromo][uid]['contig']
            
            relative_score, score_store = disco_depth_stats(otherside_dict[chromo][uid],
                        global_disco_stats_dict, score, score_store)
            
            if relative_score > 0:
                if chromo not in top_performer_dict:
                    top_performer_dict[chromo] = {}
                    
                anchor_site = ('{anchor_start}-{anchor_stop}').format(
                    anchor_start = anchor_start,
                    anchor_stop = anchor_stop)
            
                if anchor_site in top_performer_dict[chromo]:
                    if relative_score > top_performer_dict[chromo][anchor_site]['relative_score']:
                        print('Replacing: ', top_performer_dict[chromo][anchor_site])
                        top_performer_dict[chromo][anchor_site]['relative_score'] = relative_score
                        top_performer_dict[chromo][anchor_site]['uid'] = uid
                        print('with: ', top_performer_dict[chromo][anchor_site])
                        
                if anchor_site not in top_performer_dict[chromo]:
                    top_performer_dict[chromo][anchor_site] = {'relative_score': relative_score, 'uid': uid}

    top_performer_set = set()
    for chromo in top_performer_dict:
        for anchor_site in top_performer_dict[chromo]:
            uid = top_performer_dict[chromo][anchor_site]['uid']
            top_performer_set.add(uid)
            print('top_performer_set')
            print(uid)
    
    for uid in top_performer_set:
        for chromo in otherside_dict:
            if uid in otherside_dict[chromo]:
                anchor_chromo = otherside_dict[chromo][uid]['anchor_chromo']
                anchor_start = otherside_dict[chromo][uid]['anchor_start']
                anchor_stop = otherside_dict[chromo][uid]['anchor_stop']
                other_chromo = otherside_dict[chromo][uid]['other_chromo']
                other_start = otherside_dict[chromo][uid]['other_start']
                other_stop = otherside_dict[chromo][uid]['other_stop']
                score = otherside_dict[chromo][uid]['score']
                contig_set = otherside_dict[chromo][uid]['contig']
                                
                otherside = ('{other_chromo}:{other_start}-{other_stop}').format(
                    other_chromo = other_chromo,
                    other_start = other_start,
                    other_stop = other_stop)
                
                if len(contig_set) > 0:
                    contig = ''
                    for each in contig_set:
                        outline = ('{};').format(each)
                        contig += (outline)
                    contig = contig[:-1]
                else:
                    contig = '<SV>'

                process = False
                if args.return_all:
                    process = True
                else:
                    if score >= c_min_score:
                        process = True
                        
                if process:
                
                    gff_line = ('{chromo}\terisapfel\t{uid}_anchor_split'
                                '\t{start}\t{stop}\t.\t.\t{score}'
                                '\tnode_uid={uid};otherside={otherside}_breeze;contig={contig}\n').format(
                                    chromo = anchor_chromo,
                                    uid = uid,
                                    start = anchor_start,
                                    stop = anchor_stop,
                                    score = score,
                                    otherside = otherside,
                                    contig = contig)
                                    
                    gff_set.add(gff_line)
                    
                    vcf_line = ('{chromo}\t{pos}\t{uid}'
                                '\tN\t{alt}\t60\tPASS\tSVMETHOD=erisapfel_1.2;OTHERSIDE={info}'
                                '\tGT:GQ:DP\t.:.:{dp}\n').format(
                                    chromo = anchor_chromo,
                                    pos = anchor_start,
                                    uid = uid,
                                    info = otherside,
                                    alt = contig,
                                    dp = int(score))
                                    
                    vcf_set.add(vcf_line)
                    
                    anchorside = ('{chromo}:{start}-{stop}').format(
                        chromo = anchor_chromo,
                        start = anchor_start,
                        stop = anchor_stop)
            
                    rev_gff_line = ('{chromo}\terisapfel\t{uid}_breeze_split'
                                '\t{start}\t{stop}\t.\t.\t{score}'
                                '\tnode_uid={uid}_anchor;otherside={anchorside};contig={contig}\n').format(
                                    chromo = other_chromo,
                                    uid = uid,
                                    start = other_start,
                                    stop = other_stop,
                                    score = score,
                                    anchorside = anchorside,
                                    contig = contig)
                    
                    gff_set.add(rev_gff_line)
                    
                    vcf_line = ('{chromo}\t{pos}\t{uid}'
                                '\tN\t{alt}\t60\tPASS\tSVMETHOD=erisapfel_1.2;OTHERSIDE={info}'
                                '\tGT:GQ:DP\t.:.:{dp}\n').format(
                                    chromo = other_chromo,
                                    pos = other_start,
                                    uid = uid,
                                    info = anchorside,
                                    alt = contig,
                                    dp = int(score))
                                    
                    vcf_set.add(vcf_line)
            
    return(gff_set, gff_header, vcf_set, vcf_header)

def parse_msa_file(filename):
    '''
    Parse multiple sequence alignment file to dictionaries

    Parameters
    ----------
    filename : string
        multiple sequence alignment file name

    Returns
    -------
    msa_dict : dictionary containing hashed name (key) to sequence
    hash_to_name : dictionary containg hashed name (key) to Read ID 

    '''
    infile = open(filename)
    hash_to_name = {}
    msa_dict = {}
    seq = ''
    
    for line in infile:
        if line[0] == '>':
            name = (line.split('>')[1].split('_')[0].strip())
            hname = str(hash(name))
            #score = (line.split('_')[1].strip())
            
            #sname = ('{}_{}').format(name, score)
            
            if seq != '':
                msa_dict[hname] = seq
    
            hash_to_name[hname] = name
            
            seq = ''
        else:
            seq += line.strip()
            
    #picks up any trailers
    msa_dict[hname] = seq
    infile.close()
    
    return(msa_dict, hash_to_name)
    
def parse_msa_kmer(seq_1, seq_2):
    '''
    

    Parameters
    ----------
    seq_1 : string
        query string aligned
    seq_2 : string
        source string aligned
        
    Returns
    -------
    mismatch_ct: int representing mismatch between two sequences
    match_ct : int representing match 

    '''
    match_ct = 0
    mismatch_ct = 0
    phase = 0
    for each in range(len(seq_1)):
        
        if phase == 0:
            if ((seq_1[each] != '-') and (seq_2[each] != '-')):
                phase = 1
            
        if phase == 1:
            if (len(seq_1[each:]) == seq_1[each:].count('-')) or (len(seq_2[each:]) == seq_2[each:].count('-')):
                phase = 2
                
            else: 
                if (seq_1[each] != seq_2[each]):
                    mismatch_ct += 1 
                else:
                    match_ct += 1
                         
    return(mismatch_ct, match_ct)
    
def build_clusters(distance_dict, msa_dict, hash_to_name, contig_seq_dict):
    cluster = 0
    name_to_cluster_dict = {}
    cluster_to_name_dict = {}
    
    for compare, deets in distance_dict.items():
        mismatch, match = deets
        if ((mismatch <= 1) or ((mismatch/max(1,float(match))) <= 0.04)) and (match >= 25):
            name1, name2 = compare.split('_')
            
            if (name1 in name_to_cluster_dict) and (name2 in name_to_cluster_dict):
                name1_cluster_set = name_to_cluster_dict[name1]
                for eachc in name1_cluster_set:
                    name1_list = cluster_to_name_dict[eachc]
                
                    if name2 not in name1_list:
                        name_to_cluster_dict[name2].union(name1_cluster_set)
                        cluster_to_name_dict[eachc].append(name2)
                    
                name2_cluster_set = name_to_cluster_dict[name2]
                for eachc in name2_cluster_set:
                    name2_list = cluster_to_name_dict[eachc]
                
                    if name1 not in name2_list:
                        name_to_cluster_dict[name1].union(name2_cluster_set)
                        cluster_to_name_dict[eachc].append(name1)
            
            if (name1 in name_to_cluster_dict) and (name2 not in name_to_cluster_dict):
                temp_cluster_set = name_to_cluster_dict[name1]
                name_to_cluster_dict[name2] = temp_cluster_set
                for eachc in temp_cluster_set:
                    cluster_to_name_dict[eachc].append(name2)
                
            if (name1 not in name_to_cluster_dict) and (name2 in name_to_cluster_dict):
                temp_cluster_set = name_to_cluster_dict[name2]
                name_to_cluster_dict[name1] = temp_cluster_set
                for eachc in temp_cluster_set:
                    cluster_to_name_dict[eachc].append(name1)
            
            if (name1 not in name_to_cluster_dict) and (name2 not in name_to_cluster_dict) :
                new_cluster = set()
                new_cluster.add(cluster)
                name_to_cluster_dict[name1] = new_cluster
                name_to_cluster_dict[name2] = new_cluster
                cluster_to_name_dict[cluster] = [name1, name2]
                
                cluster+=1
     
    prefix_list = []
    for cluster, hnames in cluster_to_name_dict.items():
        node_name = ('{hypo}_{cluster}').format(cluster=cluster, hypo=hypo)
        cluster_seq_name = ('{temp}/{node_name}_cluster_seq.fa').format(node_name=node_name, temp=temp_dir)
        cluster_seq = open(cluster_seq_name, 'w')

        if len(hnames) > 1:
            for hname in hnames:
                name = hash_to_name[hname]
                seq = msa_dict[hname].replace('-','')
                
                outline = ('>{}\n{}\n').format(name, seq)
                cluster_seq.write(outline)
            cluster_seq.close()
            
            cluster_align_name = ('{temp}/{node_name}_cluster_aligned.msf').format(node_name=node_name, temp=temp_dir)
            bashCommand = ('mafft --auto --adjustdirection --globalpair --quiet {} > {}').format(cluster_seq_name, cluster_align_name)
            subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)
            print(bashCommand)
            
            bashCommand = ('cons -name {}_{} -plurality 1 {} -outseq {}').format(node_name, len(hnames), cluster_align_name, cluster_seq_name)
            subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)
            print(bashCommand)
            
            contig_seq_dict = check_coverage(node_name, cluster_seq_name, contig_seq_dict)
            
            prefix = ('{temp}/{node_name}').format(temp=temp_dir, node_name=node_name)
            prefix_list.append(prefix)

    return(prefix_list, contig_seq_dict)

def parse_mge(mge_file_name):
    mge_nt_lookup_dict = {}
    mge_region_lookup_dict = {}
    mge_name_to_region_dict = {}
    mge_region_to_name_dict = {}
    
    
    mge_file = open(mge_file_name)
    
    for line in mge_file:
        line = line.strip()
        
        if line[0] != '#':
            chromo = line.split('\t')[0]
            start = int(line.split('\t')[3])
            stop = int(line.split('\t')[4])
            name = line.split('\t')[8]
                        
            mge_name_to_region_dict[name] = (start, stop)
            
            if chromo not in mge_region_to_name_dict:
                mge_region_to_name_dict[chromo] = {}
            
            mge_region_to_name_dict[chromo][(start,stop)] = name
                        
            if chromo not in mge_nt_lookup_dict:
                mge_nt_lookup_dict[chromo] = {}
                
            for nt in range((start - flank_nt), (stop + flank_nt + 1)):
                if nt not in mge_nt_lookup_dict[chromo]:
                    mge_nt_lookup_dict[chromo][nt] = set()
                
                mge_nt_lookup_dict[chromo][nt].add(name)
                
            if chromo not in mge_region_lookup_dict:
                mge_region_lookup_dict[chromo] = set()
                
            mge_region_lookup_dict[chromo].add((start, stop))
                
    mge_file.close()
    
    return(mge_nt_lookup_dict, mge_region_lookup_dict, mge_name_to_region_dict, mge_region_to_name_dict)

if args.mge_file:
    run_mge = True
    mge_nt_lookup_dict, mge_region_lookup_dict, mge_name_to_region_dict, mge_region_to_name_dict = parse_mge(args.mge_file)
               

""" Step One """
if args.make:
    """ This command uses bwa mem to align the paired end (required) fastq files to the reference genome
    These are passed to samblaster which seperates out discordant and split reads
    and then samtools for conversion to bam and bam index generation.
    
    Resource file is created to record and maintain run associated 
    information such as file names and parameters. This is located in 
    /idx/resources.tab

    """
    
    chromo_size_dict = {}
    resource_dict = {}
    
    io_make()
    
    resource_dict['run_name']=args.run_name
    resource_dict['results_dir']=args.run_name
    resource_dict['genome_fa']=args.fa_file
    
    if not args.paired_end:
        end_is = 'pe'
        resource_dict['end_is']=end_is
        resource_dict['fastq_1']=args.fastq_1
        resource_dict['fastq_2']=args.fastq_2
    else:
        end_is = args.paired_end    
        resource_dict['end_is']=end_is
        resource_dict['fastq_1']=args.fastq_1
    
    command_file_name=('{}_command_file.sh').format(s_path)
    command_file = open(command_file_name,'w')
    command_file.write('module load bwa/intel/0.7.17\n'
                       'module load samtools/intel/1.14\n'
                       'module load bedtools/intel/2.29.2\n'
                       'module load velvet/1.2.10\n'
                       'module load blast+/2.11.0\n'
                       'module load samblaster/0.1.26\n'
                       'module load mafft/intel/7.475\n'
                       'module load emboss/intel/6.6.0\n\n\n')
    
    ref_fa=('genome_fa={}\n').format(args.fa_file)
    #ref_gff=('\tgenome_gff={}\n').format(args.gff_file)
         
    if (end_is == 'pe') or (end_is == 'paired_end'):
        infastq_1=('fastq_1={}\n').format(args.fastq_1)
        infastq_2=('fastq_2={}\n').format(args.fastq_2)
    else:
        infastq_1=('fastq_1={}\n').format(args.fastq_1)
        infastq_2=('fastq_2={}\n').format("")
        
    output_filename = ('output_file={bam_dir}/{output_file}\n'
                     'mkdir -p {bam_dir}\n'
                     'echo bam_dir $PWD + {bam_dir}\n').format(
        bam_dir=bam_dir, 
        output_file=output_file)
                         
    #command_file.write(ref_fa + ref_gff + infastq_1 + infastq_2 + output_filename)
    command_file.write(ref_fa + infastq_1 + infastq_2 + output_filename)
    
    #using samblaster to extract discordant (-d) and split (-s) reads, repeat reads are excluded (-e) 
    if (end_is == 'int') or (end_is == 'interleaved'):   
        outline = ('bwa mem -M -t 16 -p ${genome_fa} ${fastq_1} | samblaster -M -e -d ${output_file}_discordant.sam -s ${output_file}_split.sam | samtools sort -@16 -o ${output_file}.bam\n')
        command_file.write(outline)
    
    if (end_is == 'pe') or (end_is == 'paired_end'):   
        outline = ('bwa mem -M -t 16 ${genome_fa} ${fastq_1} ${fastq_2} | samblaster -M -e -d ${output_file}_discordant.sam -s ${output_file}_split.sam | samtools sort -@16 -o ${output_file}.bam\n')
        command_file.write(outline)

    if (end_is == 'se') or (end_is == 'single_end'):   
        outline = ('bwa mem -M -t 16 ${genome_fa} ${fastq_1} | samblaster -M -e -s ${output_file}_split.sam | samtools sort -@16 -o ${output_file}.bam\n')
        command_file.write(outline)   

    outline = ('bedtools genomecov -ibam ${output_file}.bam -bg > ${output_file}.bedGraph\n'+
    'samtools index ${output_file}.bam\n'+
    'samtools view -h ${output_file}.bam -o ${output_file}.sam\n'+
    'samtools view -bT ${genome_fa} ${output_file}_discordant.sam > ${output_file}_discordant_u.bam\n'+
    'samtools sort -T tmp_sort -o ${output_file}_discordant.bam ${output_file}_discordant_u.bam\n'+
    'samtools index ${output_file}_discordant.bam\n'+
    'samtools view -bT ${genome_fa} ${output_file}_split.sam > ${output_file}_split_u.bam\n'+
    'samtools sort -T tmp_sort -o ${output_file}_split.bam ${output_file}_split_u.bam\n'+
    'samtools index ${output_file}_split.bam\n'+
    'rm ${output_file}_discordant_u.bam\n'+
    'rm ${output_file}_split_u.bam\n\n')
    command_file.write(outline)
    command_file.close()
    
    bash_command = ('bash {}').format(command_file_name)
    subprocess.run([bash_command],stderr=subprocess.STDOUT,shell=True)
    
    header_file_name = ('{}{}_split.sam').format(bam_dir, output_file)
    read_file = open(header_file_name)
    ct = 0
    for line in read_file:
        if line.split('\t')[0] == '@SQ':
            chromo_name = line.split('\t')[1].split(':')[1]
            chromo_length = line.split('\t')[2].split(':')[1]
            chromo_size_dict[chromo_name] = chromo_length

    pickle_name = ('{}/{}_chromo_size.p').format(pickles_dir, output_file)
    pickle.dump(chromo_size_dict, open(pickle_name, 'wb'))
    
    resource_dict['bam_file'] = (bam_dir+output_file+str('.bam'))
    resource_dict['discordant_file'] = (bam_dir+output_file+str('_discordant.bam'))
    resource_dict['split_file'] = (bam_dir+output_file+str('_split.bam'))    
    #resource_dict['soft_file'] = (bam_dir+output_file+str('_soft.bam'))
    
    resource_dict['sam_file'] = (bam_dir+output_file+str('.sam'))
    resource_dict['discordant_sam_file'] = (bam_dir+output_file+str('_discordant.sam'))
    resource_dict['split_sam_file'] = (bam_dir+output_file+str('_split.sam'))
    #resource_dict['soft_sam_file']=(bam_dir+output_file+str('_soft.sam'))
    
    resource_dict['chromo_size']=(pickle_name)
    
    #handle read_type
    if not args.read_type_list:
        read_type_list = ['RD','discordant','split']
    else:
        read_type_list = args.read_type_list
        
    resource_dict['read_types']=read_type_list

    io_append(resource_dict)
    
""" Step Two"""
if args.depth_analysis:    
    """
    Depth_analysis
        Function: For each read type ('RD','discordant','split')
        1. Open corresponding bam file
        2. Pulls read counts from bam files into chromosomes.
        3. Depth statistics, medians and std are calculated. 
        * NB - these stats 
    """    
    resource_dict = io_load()
    print(resource_dict)
        
    fa_file = resource_dict['genome_fa']
    each_sample = resource_dict['run_name']
    bam_file = resource_dict['bam_file']
    sam_file = resource_dict['sam_file']
    # soft_file = resource_dict['soft_file']
    # soft_sam_file = resource_dict['soft_sam_file']
    read_type_list = resource_dict['read_types']
    
    if not args.filter_bed:
        filter_by_bed = False
        
    if args.filter_bed:
        if len(args.filter_bed) == 0:
            filter_by_bed = False
        else:
            filter_by_bed = True
            if len(args.filter_bed):
                filter_bed = args.filter_bed[0]
                resource_dict['depth_region_filter']=filter_bed
            else:
                print('Please only specify a single bed file for filtering')
                1/0
                                        
    for each_type in read_type_list:
        if each_type == 'RD':
            infile_name = resource_dict['bam_file']
        if each_type == 'split':
            infile_name = resource_dict['split_file']
        if each_type == 'discordant':
            infile_name = resource_dict['discordant_file']
        
        outline = ('\nProcessing sample {}...').format(each_sample)
        print(outline)
        
        run_bedgraph(bam_file, each_sample)
        
        print('\tParsing mpileup into chromosomes...')
        run_mpileup(infile_name, each_sample, each_type)

        populated=False
        mpileup_df = []
        raw_mpileup = str(pickles_dir + '/' + each_sample + '_' + each_type + '.mpileup')
        temp_df = pd.read_csv(raw_mpileup, sep='\t', header=None, names = ['chromo', 'nuc', 'base', 'ct', 'val1', 'val2'])

        outline = ('\tProcessing {}...').format(each_type)
        print(outline)
                
        if populated:
            mpileup_df['chromo']=temp_df['chromo']
            
        if not populated:
            populated = True
            mpileup_df = temp_df[['ct']].copy()
            mpileup_df['nuc']=temp_df['nuc']
            mpileup_df['chromo']=temp_df['chromo']
                                   
        pickle_name = ('{}/mpileup_{}_{}.p').format(pickles_dir, each_sample, each_type)
        mpileup_df.to_pickle(pickle_name)
        print(mpileup_df.head())

        fg_median = mpileup_df["ct"].mean()
        fg_std = mpileup_df["ct"].std()
        
        zscore_read_type = ('{}_zscore_median').format(each_type)
        resource_dict[zscore_read_type] = fg_median
        zscore_read_type = ('{}_zscore_std').format(each_type)
        resource_dict[zscore_read_type] = fg_std
            
    io_append(resource_dict)
        
""" Step Three """
if args.peaks:
    """
    Peak detector
    1. load 'mpileup' files into dataframe pickles.    additional feature, filter regions of chromosome with large deletions
        
    """
    resource_dict = io_load()
    read_type_list = resource_dict['read_types']
    
    each_sample = resource_dict['run_name']
    
    chromo_size_p = pickle_loader(resource_dict['chromo_size'], 'dict')
    
    chromo_list = []
    for chromo in chromo_size_p:
        chromo_list.append(chromo)     
        
    if args.CNV_min_length:
        min_length = int(args.CNV_min_length)
    else:
        min_length = 300
        
    if args.max_gap:
        gap = int(args.max_gap)
    else:
        gap = 300

    if args.min_purity:
        purity = (float(args.min_purity))
    else:
        purity = 0.9
                        
    file_name_lookup = []
    depth_dict = {}
    # chromo_depth_dict = {}
    
    for each_type in read_type_list:
        print('Processing ', each_type)
        
        depth_gff_file_name = ('{}/{}_{}_depth.gff').format(final_output_dir, each_sample, each_type)
        depth_gff_outfile = open(depth_gff_file_name,'w')
        
        chromo_gff_file_name = ('{}/{}_{}_chromosome.gff').format(final_output_dir, each_sample, each_type)
        chromo_gff_outfile = open(chromo_gff_file_name,'w')
    
        pickle_name = ('{}/mpileup_{}_{}.p').format(pickles_dir, each_sample, each_type)
        mpileup_df = pickle_loader(pickle_name, 'df')
                
        zscore_read_type = ('{}_zscore_median').format(each_type)
        g_median = resource_dict[zscore_read_type]
        zscore_read_type = ('{}_zscore_std').format(each_type)
        g_std = resource_dict[zscore_read_type]
                        
        outline = ('Global Median:\t{}\nGlobal Std:\t{}\n').format(g_median,g_std)
        print(outline)
         
        for every_chromo in chromo_list:
            m_chromo_df = mpileup_df.loc[mpileup_df["chromo"] == every_chromo]
            
            c_median = m_chromo_df["ct"].mean()
            c_std = m_chromo_df["ct"].std()
            chromo_annu = c_median+c_std
            
            start, stop = 1, len(m_chromo_df)
            outline = ('{chromo}\t{each_type}_depth\tChromosome\t{start}\t{stop}\t{c_median}\t.\t.\tID={chromo}:{start}-{stop};Median_RD={c_median};Std_RD={c_std};rel_genome_RD={rel_g};sample={sample}\n').format(chromo=every_chromo, each_type=each_type, start=start, stop=stop, c_median=c_median, c_std=c_std, rel_g=round(c_median/float(g_median),2), sample=each_sample)
            chromo_gff_outfile.write(outline)
        
        for every_chromo in chromo_list:
            
            m_chromo_df = mpileup_df.loc[mpileup_df["chromo"] == every_chromo]
    
            if each_type == 'RD':
                fzero = m_chromo_df.replace(0, np.NaN)
                c_median = fzero["ct"].mean()
                c_std = fzero["ct"].std()
                
            else:
                c_median = m_chromo_df["ct"].mean()
                c_std = m_chromo_df["ct"].std()              
            
            #this selects nts with more than some zscore of depth
            subz_df = m_chromo_df.loc[(m_chromo_df["ct"] >= (c_median + 3*c_std))]
            suby_df = m_chromo_df.loc[m_chromo_df["ct"] >= (c_median + c_std)]
            
            build_trace(subz_df, suby_df, each_type, c_median, g_median)
            
            # this detects deletions in the DNA depth, independant of discordance
            if each_type == 'RD':
                subz_df = m_chromo_df.loc[(m_chromo_df["ct"]) == 0]
                suby_df = m_chromo_df.loc[(m_chromo_df["ct"]) == 0]
                
                build_trace(subz_df, suby_df, each_type, c_median, g_median)
                                
        depth_gff_outfile.close()
        chromo_gff_outfile.close()
        
        depth_pickle_name = ('{}/peaks_{}_{}.p').format(pickles_dir, each_sample, each_type)         
        with open(depth_pickle_name, 'wb') as file:
            pickle.dump(depth_dict, file)
            
def calculate_over_frame(start, stop, chromo, df):
    mge_nt_list = []
    
    for each in range(start,stop +1):
        mge_nt_list.append(each)
    
    chromo_df = df.loc[df["chromo"] == chromo]
    mge_df = chromo_df.query('nuc in @mge_nt_list')
        
    return(mge_df) 

def nan_inf_cleaner(val):
    if np.isnan(val):
        val = 0
    
    if np.isinf(val):
        val = 0
        
    return(val)

def weighted_outlier(df, region):
    mean, std = df.mean(), df.std()
    
    top = df[df >= mean+std]
    bottom = df[df <= mean-std]
        
    if bottom.mean(skipna=True) and bottom.std(skipna=True):
    #if (sum(bottom) > 1):
        diff = ( (top.mean(skipna=True) - top.std(skipna=True)) - (bottom.mean(skipna=True) + bottom.std(skipna=True)) )
        rel_outlier = (len(top)+len(bottom))/max(len(df),1)
        wo = nan_inf_cleaner( rel_outlier * diff)
        
        
    else:
        diff = ((top.mean(skipna=True) - top.std(skipna=True)))
        rel_outlier = (len(top))/max(len(df),1)
        wo = nan_inf_cleaner(rel_outlier * diff)    
   
    return(wo)


def tea_stats(df_num, df_den, region):
    num_cov = nan_inf_cleaner(stats.variation(np.array(df_num)))
    
    num_sum = df_num.sum()
    
    num_wo = weighted_outlier(df_num, region)
    
    num_diff = max(df_num)-min(df_num)
    num_entropy = nan_inf_cleaner(stats.entropy(np.array(df_num)))
    
    den_cov = nan_inf_cleaner(stats.variation(np.array(df_den)))
    
    den_sum = df_den.sum()
    
    den_wo = weighted_outlier(df_den, region)
    
    den_diff = max(df_den)-min(df_den)
    den_entropy = nan_inf_cleaner(stats.entropy(np.array(df_den)))
    
    cov = num_cov/max(den_cov,1)
    rsum = num_sum/max(den_sum,1)
    wo = num_wo/max(den_wo,1)
    diff = num_diff/max(den_diff,1)
    entropy = num_entropy/max(den_entropy,1)
        
    return(cov, rsum, wo, diff, entropy)

def tea_features(region, chromo, rd_mpileup_df, runmode):
    
    temp_dict = {}
    
    if runmode == 'rd':
        adjust = 300
        length = max(region) - min(region) + 2*adjust
        
        left_outside_start = min(region)-int(adjust)
        left_outside_stop = min(region)
        
        inside_start = min(region)  
        inside_stop = max(region)
        
        right_outside_start = max(region)  
        right_outside_stop = max(region)+int(adjust)
        
        left_outside = calculate_over_frame(left_outside_start, left_outside_stop, chromo, rd_mpileup_df)
        inside = calculate_over_frame(inside_start, inside_stop, chromo, rd_mpileup_df)
        
        left_cov, left_sum, left_wo, left_diff, left_entropy = tea_stats(inside['ct'], left_outside['ct'], region)
    
        right_outside = calculate_over_frame(right_outside_start, right_outside_stop, chromo, rd_mpileup_df)
        #right_inside = calculate_over_frame(right_inside_start, right_inside_stop, chromo, rd_mpileup_df)
        
        right_cov, right_sum, right_wo, right_diff, right_entropy = tea_stats(inside['ct'], right_outside['ct'], region)
        
        df_cov = (left_cov + right_cov)/length
        df_sum = (left_sum + right_sum)/length
        df_wo = (left_wo + right_wo)/length
        df_diff = (left_diff + right_diff)/length
        df_entropy = (left_entropy + right_entropy)/length
    
    if runmode == 'disco':
        adjust = 300
        
        left_outside_start = min(region)-int(adjust)
        left_outside_stop = min(region)-75
        left_inside_start = min(region)-75  
        left_inside_stop = min(region)+int(adjust)
        
        right_inside_start = max(region)-int(adjust)
        right_inside_stop = max(region)+75
        right_outside_start = max(region)+75  
        right_outside_stop = max(region)+int(adjust)
    
        left_outside = calculate_over_frame(left_outside_start, left_outside_stop, chromo, rd_mpileup_df)
        left_inside = calculate_over_frame(left_inside_start, left_inside_stop, chromo, rd_mpileup_df)
        
        left_cov, left_sum, left_wo, left_diff, left_entropy = tea_stats(left_inside['ct'], left_outside['ct'],region)
    
        right_outside = calculate_over_frame(right_outside_start, right_outside_stop, chromo, rd_mpileup_df)
        right_inside = calculate_over_frame(right_inside_start, right_inside_stop, chromo, rd_mpileup_df)
        
        right_cov, right_sum, right_wo, right_diff, right_entropy = tea_stats(right_inside['ct'], right_outside['ct'], region)
        
        df_cov = left_cov + right_cov
        df_sum = left_sum + right_sum
        df_wo = left_wo + right_wo
        df_diff = left_diff + right_diff
        df_entropy = left_entropy + right_entropy
    
    temp_dict = {'df_cov': df_cov, 'df_sum': df_sum,
                 'df_wo': df_wo, 'df_diff': df_diff,
                 'df_entropy': df_entropy}
    
    return(temp_dict)
            
def mge_trace_three_frames(mpileup_df, rd_mpileup_df, split_mpileup_df):
    global mge_region_lookup_dict
    
    trace_dict = {}
    
    for chromo in chromo_list:
        if chromo in mge_region_lookup_dict:
            sub_mg_set = mge_region_lookup_dict[chromo]
            
            if chromo not in trace_dict:
                trace_dict[chromo] = {}
                            
            for region in sub_mg_set:            
                if region in trace_dict[chromo]:
                    print('Region error: ', region)
                    1/0
                    
                trace_dict[chromo][region] = {}
                trace_dict[chromo][region] = {'rd': {}, 'disco': {}, 'split':{}}
                
                trace_dict[chromo][region]['rd'] = tea_features(region, chromo, rd_mpileup_df, 'rd')
                trace_dict[chromo][region]['disco'] = tea_features(region, chromo, mpileup_df, 'disco')
                trace_dict[chromo][region]['split'] = tea_features(region, chromo, split_mpileup_df, 'disco')
                
    return(trace_dict)

def populate_confusion_matrix(val, tval, rdict):
    if val == tval:
        if tval:
            rdict['TP']+=1
        else:
            rdict['TN']+=1
    else:
        if tval:
            rdict['FN']+=1
        else:
            rdict['FP']+=1    
    return(rdict)

def print_confusion_matrix(name, rdict, read_type, outfile):
    outline = ('{name}\t{read_type}\tTP:{tp}\tFP:{fp}\tFN:{fn}\tTN:{tn}\n').format(
        name = name,
        read_type = read_type,
        tp = rdict['TP'],
        fp = rdict['FP'],
        fn = rdict['FN'],
        tn = rdict['TN'])
            
    outfile.write(outline)
    #return(outline)
    
def calc_value_weight(result_val, result_dict):
    denom = (result_dict['FP']+result_dict['TP'])
    #denom = (result_dict['FP']+result_dict['TP']+result_dict['FN'])
    s_val = int(result_val)*(result_dict['TP']/max(denom,1))
    return (s_val)
    
def norm_value(val_list, total, runmode):
    
    if runmode == 'count':
        total = max(sum(val_list),1)
        val_list = [x/total for x in val_list]
    
        val_median = np.median(val_list)
        val_std = np.std(val_list)
        
        return(val_median, val_std, total)
    
    if runmode == 'distribution':
        val_median = np.median(val_list)
        val_std = np.std(val_list)
        
        return(val_median, val_std)
    
    if runmode == 'norm':
        val = val_list / max(total, 1)
        
        return(val)
    
def mge_to_model_parser(mge_file_name):
    value_df_object = {}

    mge_file = open(mge_file_name)
    
    for line in mge_file:
        if line[0] != '#':
            line = line.strip()
            sample_name, sample_path = line.split('\t')
            
            if sample_name in value_df_object:
                if sample_path != value_df_object[sample_name]:
                    print('Error: ensure that sample_name (first column) is unique')
                    exit()
                    
            value_df_object[sample_name] = sample_path
    
    mge_file.close()
    
    return(value_df_object) 
    
    
def model_pipe(sample, X_train, X_test, y_train, y_test, 
               preprocess, sampler, classifier, performance_dict, log_file):
    
    pipeline = [make_pipeline(preprocess[0], sampler[0], classifier[0])]
    
    for model in pipeline:
        perf_name = ('{}\t{}\t{}\t').format(preprocess[1], sampler[1], classifier[1])
        print(perf_name)
        
        try:
            toc = time.perf_counter()
            
            model.fit(X_train, y_train)
            y_pred = model.predict(X_test)
            
            tic = time.perf_counter()
            
            gmean = geometric_mean_score(y_test, y_pred)
            fval = f1_score(y_test, y_pred)
            tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
            time_elapsed = tic-toc
            
            outline = ('{strain}\t{perf_name}\t'
                       '{gmean}\t{fval}\t{tp}\t{fp}\t{fn}\t{tn}\t{time_elapsed}\n').format(
                           strain=sample, perf_name=perf_name,
                           gmean=round(gmean,3), fval = round(fval,3),
                           tp = tp, fp = fp, fn = fn, tn = tn,
                           time_elapsed = time_elapsed
                           )
            #print(outline)
            log_file.write(outline)
            
            if perf_name not in performance_dict:
                performance_dict[perf_name] = {'gmean':[],'fval':[], 'time_elapsed':[]}
                
            performance_dict[perf_name]['gmean'].append(gmean)
            performance_dict[perf_name]['fval'].append(fval)
            performance_dict[perf_name]['time_elapsed'].append(time_elapsed)
            
        except:
            outline = ('{}Exception occurred applying model\n').format(perf_name)
            #print(outline)
            log_file.write(outline)
            
    return(performance_dict)
    
def calc_final_output_performance_values(perfo, performance_dict):
    gmean = round(np.median(performance_dict[perfo]['gmean']),3)
    gstd = round(np.std(performance_dict[perfo]['gmean']),3)
    fval = round(np.median(performance_dict[perfo]['fval']),3)
    fstd = round(np.std(performance_dict[perfo]['fval']),3)
    time_elapsed = round(np.median(performance_dict[perfo]['time_elapsed']),3)
    
    score=(gmean-gstd)+(fval-fstd)
    
    return(gmean, gstd, fval, fstd, time_elapsed, score)

def meta_selection(sample, X, y, performance_dict, log_file):
    test_number = 100
    if args.metaparameters == 'fast':
        test_number = 10
    
    outline = ("Running {sample} for {test_number} runs ...").format(
        sample=sample, test_number = test_number)
    
    print(outline)
    
    for test in range(test_number + 1):
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20, stratify=y)
        
        wval = len(y_train[y_train["onehot"]>0])/len(y_train)
                                     
        weights = {0:wval, 1:1.0}
        
        if args.metaparameters == 'fast':
            for preprocess in [('passthrough','none'),
                               (MinMaxScaler(),'MinMaxScaler'),
                               ]:
            
                for classifier in [(LogisticRegression(class_weight=weights, max_iter=1000), 'LogisticRegression_lbgfs'),
                                   (LogisticRegression(solver='liblinear', class_weight=weights, max_iter=1000), 'LogisticRegression_liblinear'),
                                   ]:
                    #
                    for sampler in [('passthrough','none'),
                                    (FunctionSampler(), 'FunctionSampler'),
                                    (AllKNN(), 'AllKNN')
                                    ]:
                        #
                        performance_dict = model_pipe(sample, X_train, X_test, y_train, y_test, 
                                       preprocess, sampler, classifier, performance_dict, log_file)
        
        if args.metaparameters == 'robust':
            for preprocess in [('passthrough','none'),
                               (MinMaxScaler(),'MinMaxScaler'),
                               (StandardScaler(), 'StandardScaler'),
                               (QuantileTransformer(),'QuantileTransformer'),
                               ]:
            
                for classifier in [(KNeighborsClassifier(), 'KNeighborsClassifier'),
                                   (LogisticRegression(class_weight=weights, max_iter=1000), 'LogisticRegression_lbgfs'),
                                   (ComplementNB(), 'ComplementNB'),
                                   ]:
                    #
                    for sampler in [('passthrough','none'),
                                    (FunctionSampler(), 'FunctionSampler'),
                                    (ADASYN(), 'ADASYN'),
                                    (AllKNN(), 'AllKNN')
                                    ]:
                        #
                        performance_dict = model_pipe(sample, X_train, X_test, y_train, y_test, 
                                       preprocess, sampler, classifier, performance_dict, log_file)
                        
        if args.metaparameters == 'complete':
            for preprocess in [('passthrough','none'),
                               (MinMaxScaler(),'MinMaxScaler'),
                               (StandardScaler(), 'StandardScaler'),
                               (Normalizer(), 'Normalizer'),
                               (PowerTransformer(),'PowerTransformer'),
                               (QuantileTransformer(),'QuantileTransformer'),
                               ]:
            
                for classifier in [(KNeighborsClassifier(), 'KNeighborsClassifier'),
                                   (GaussianNB(), 'GaussianNB'),
                                   (LogisticRegression(class_weight=weights, max_iter=1000), 'LogisticRegression_lbgfs'),
                                   (LogisticRegression(solver='liblinear', class_weight=weights, max_iter=1000), 'LogisticRegression_liblinear'),
                                   (ComplementNB(), 'ComplementNB'),
                                   (RandomForestClassifier(), 'RandomForestClassifier'),
                                   ]:
                    #
                    for sampler in [('passthrough','none'),
                                    (FunctionSampler(), 'FunctionSampler'),
                                    (RandomOverSampler(), 'RandomOverSampler'),
                                    (ADASYN(), 'ADASYN'),
                                    (SMOTE(), 'SMOTE'), 
                                    (RandomUnderSampler(), 'RandomUnderSampler'),
                                    (NearMiss(version=1), 'NearMiss'),
                                    (AllKNN(), 'AllKNN')
                                    ]:
                        #
                        performance_dict = model_pipe(sample, X_train, X_test, y_train, y_test, 
                                       preprocess, sampler, classifier, performance_dict, log_file)
                    
    return(performance_dict)
        
def model_selection(value_df_object):
        
    #TODO make into standard temp_dir 
    log_file_name = ('{}{}_model_selection_run.log').format(output_dir, output_file)
    log_file = open(log_file_name, 'w')
    
    header = ('sample\tmodel\tgmean\tfval\tTP\tFP\tFN\tTN\tTime_elapsed\n')
    log_file.write(header)
    
    #TODO remove
    performance_dict = {}
    #
    for sample in value_df_object:
        
        inpickle = value_df_object[sample]
        print(sample, inpickle)
        value_df = pd.read_pickle(inpickle) 
                     
        X = value_df.drop(['onehot'], axis=1)
        y = value_df[["onehot"]]
        
        performance_dict = meta_selection(sample, X, y, performance_dict, log_file)
        
    log_file.close()

    results_name = ('{}{}_model_selection_best.log').format(output_dir, output_file)
    results_file = open(results_name, 'w')

    rank_dict = {}
    score_list = []
    
    for perfo in performance_dict:
        _gmean, _gstd, _fval, _fstd, _time_elapsed, score = calc_final_output_performance_values(perfo, performance_dict)
        
        score_list.append(score)
        
    percentile_array = np.percentile(np.array(score_list), range(0,101,10))
    
    for perfo in performance_dict:
        gmean, gstd, fval, fstd, time_elapsed, score = calc_final_output_performance_values(perfo, performance_dict)
        
        outline = ('{perfo}\tscore: {score}\tgmean: {gmean}\t'
                   'gmean_std: {gstd}\tf1: {fval}\tf1_std: {fstd}\ttime_elapsed:{time_elapsed}\n').format(
                       perfo = perfo, score = score, gmean = gmean, gstd = gstd,
                       fval = fval, fstd = fstd, 
                       time_elapsed = time_elapsed)

        pct = calc_percentile(round(score, 5), percentile_array)
        
        if pct > 90:
            print(outline)
            
        results_file.write(outline)
            
        if score not in rank_dict:
            rank_dict[score] = set()
            
        rank_dict[score].add(perfo)
        
    results_file.close()
        
    rank_list = list(rank_dict.keys())
    rank_list.sort(reverse=True)
    
    if len(rank_dict[rank_list[0]]) > 1:
        new_rank_dict = {}
        for perfo in rank_dict[rank_list[0]]:
            gmean, gstd, fval, fstd, time_elapsed, _score = calc_final_output_performance_values(perfo, performance_dict)
            score=(gmean-gstd)+(fval-fstd)-time_elapsed
            
            if score not in new_rank_dict:
                new_rank_dict[score] = set()
                
            new_rank_dict[score].add(perfo)
            
        new_rank_list = list(new_rank_dict.keys())
        new_rank_list.sort(reverse=True)
        
        for perfo in new_rank_dict[new_rank_list[0]]:
            gmean, gstd, fval, fstd, time_elapsed, _score = calc_final_output_performance_values(perfo, performance_dict)
            score=(gmean-gstd)+(fval-fstd)-time_elapsed
        
            outline = ('Best performer: {perfo}\tscore: {score}\tgmean: {gmean}\t'
                       'gmean_std: {gstd}\tf1: {fval}\tf1_std: {fstd}\ttime_elapsed: {time_elapsed}\n').format(
                           perfo = perfo, score = score, gmean = gmean, gstd = gstd,
                           fval = fval, fstd = fstd,
                           time_elapsed = time_elapsed)
            
            print(outline)
        
    if len(rank_dict[rank_list[0]]) == 1:
        for perfo in rank_dict[rank_list[0]]:
            gmean = round(np.median(performance_dict[perfo]['gmean']),3)
            gstd = round(np.std(performance_dict[perfo]['gmean']),3)
            fval = round(np.median(performance_dict[perfo]['fval']),3)
            fstd = round(np.std(performance_dict[perfo]['fval']),3)
            time_elapsed = round(np.std(performance_dict[perfo]['time_elapsed']),3)
            
            score=(gmean-gstd)+(fval-fstd)
            
            outline = ('Best performer: {perfo}\tscore: {score}\tgmean: {gmean}\t'
                       'gmean_std: {gstd}\tf1: {fval}\tf1_std: {fstd}\ttime_elapsed: {time_elapsed}\n').format(
                           perfo = perfo, score = score, gmean = gmean, gstd = gstd,
                           fval = fval, fstd = fstd,
                           time_elapsed = time_elapsed)
            
            print(outline)
            
def model_train(value_df_object):
    
    preprocessor_name = args.preprocessor
    sampler_name = args.sampler
    classifier_name = args.classifier
    
    master_value_df = pd.DataFrame()
    
    for sample in value_df_object:    
        inpickle = value_df_object[sample]
        with open(inpickle, 'rb') as fp:
            value_df = pickle.load(fp) 
            
        outline = ('Loading {}...\n Addding {} events...').format(sample, len(value_df))
        print(outline)
        
        master_value_df = pd.concat([master_value_df, value_df], ignore_index=True)
                     
    print('Total events evaluated: ', str(len(master_value_df)))
    #master_value_df
    X = master_value_df.drop(['onehot'], axis=1)
    y = master_value_df[["onehot"]]
        
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20, stratify=y)
    
    wval = len(y_train[y_train["onehot"]>0])/len(y_train)
                                 
    weights = {0:wval, 1:1.0}
    
    preprocess_lookup = {'none': 'passthrough',
                       'MinMaxScaler': MinMaxScaler(),
                       'StandardScaler': StandardScaler(),
                       'Normalizer': Normalizer(),
                       'PowerTransformer': PowerTransformer(),
                       'QuantileTransformer': QuantileTransformer(),
                       }
    
    classifier_lookup = {'KNeighborsClassifier': KNeighborsClassifier(),
                         'GaussianNB': GaussianNB(),
                         'LogisticRegression_lbgfs': LogisticRegression(class_weight=weights, max_iter=1000),
                         'LogisticRegression_liblinear': LogisticRegression(solver='liblinear', class_weight=weights, max_iter=1000),
                         'ComplementNB': ComplementNB(),
                         'RandomForestClassifier': RandomForestClassifier(),
                        }
                    #
    sampler_lookup = {'none':'passthrough',
                    'FunctionSampler': FunctionSampler(),
                    'RandomOverSampler': RandomOverSampler(),
                    'ADASYN': ADASYN(),
                    'SMOTE': SMOTE(),
                    'RandomUnderSampler': RandomUnderSampler(),
                    'NearMiss': NearMiss(version=1),
                    'AllKNN': AllKNN(),
                    }
    
    if preprocessor_name in preprocess_lookup:
        preprocess = preprocess_lookup[preprocessor_name]
    else:
        print('Error identifying preprocesor.')
        print('Please specify: ', preprocess_lookup.keys())
        
    if sampler_name in sampler_lookup:
        sampler = sampler_lookup[sampler_name]
    else:
        print('Error identifying sampler.')
        print('Please specify: ', sampler_lookup.keys())
        
    if classifier_name in classifier_lookup:
        classifier = classifier_lookup[classifier_name]
    else:
        print('Error identifying classifier.')
        print('Please specify: ', classifier_lookup.keys())
        
    pipeline = [make_pipeline(preprocess, sampler, classifier)]
    
    for model in pipeline:
        perf_name = ('{}\t{}\t{}\t').format(preprocessor_name, sampler_name, classifier_name)
        print(perf_name)
        
        try:
            toc = time.perf_counter()
            
            model.fit(X_train, y_train)
            y_pred = model.predict(X_test)
            
            tic = time.perf_counter()
            
            gmean = geometric_mean_score(y_test, y_pred)
            fval = f1_score(y_test, y_pred)
            tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
            time_elapsed = tic-toc
            
            outline = ('Top_model:\t{perf_name}\t'
                       'gmean: {gmean}\tfval: {fval}\tTP: {tp}\tFP: {fp}\tFN: {fn}\tTN: {tn}\ttime_elapsed: {time_elapsed}\n').format(
                           perf_name=perf_name,
                           gmean=round(gmean,3), fval = round(fval,3),
                           tp = tp, fp = fp, fn = fn, tn = tn,
                           time_elapsed = time_elapsed
                           )
            print(outline)
            
            pickle_out = ("{}").format(args.model_object)
            print('Model saved to: ', pickle_out)
            pickle.dump(model, open(pickle_out, 'wb'))
            
            return(True)
            
        except:
            outline = ('{} Exception occurred applying model\n').format(perf_name)
            1/0
            
    return(False)


def model_predict(value_df_object):    
    model_name = args.model_object
    with open(model_name, 'rb') as fp:
        model = pickle.load(fp)
     
    master_value_df = pd.DataFrame()
    for sample in value_df_object:    
        inpickle = value_df_object[sample]
        
        with open(inpickle, 'rb') as fp:
            master_value_df = pickle.load(fp)
          
        try:
            outfile_name = ('{}/{}_predicted_mge.gff').format(final_output_dir, sample) 
            outfile = open(outfile_name, 'w')
                
            outline = ('Loading {}...\n with {} events...').format(sample, len(master_value_df))
            print(outline)
            
            print('Total events evaluated: ', str(len(master_value_df)))
            
            X = master_value_df.drop(['onehot'], axis=1)
            y_pred = model.predict(X)
            
            y_pred_prob = model.predict_proba(X)[:,1]
            
            master_value_df['predict'] = y_pred
            master_value_df['confidence'] = y_pred_prob
            
            csv_name = ('{}/{}_predicted_mge.csv').format(final_output_dir, sample)
            master_value_df.to_csv(csv_name)
            
            for index, row in master_value_df.iterrows():
                chromo, regions = index.split(':')
                start, stop = regions.split('-')
                name = index
                
                process = False
                
                if args.return_all:
                    process = True
                else:
                    if (row['predict'] == 1) or (row['onehot'] == 1):
                        process = True
                        
                if process:                
                    outline = ('{chromo}\tmge\tmge_prediction\t{start}\t{stop}\t{predict}'
                               '\t.\t.\tID={name}; score={score}; predict={predict};'
                               'confidence={confidence}; exp={onehot}\n').format(
                                   chromo = chromo,
                                   start = start, stop = stop,
                                   predict = row['predict'], name = name, 
                                   score = row['total_score'], confidence = row['confidence'],
                                   onehot = row['onehot'])
                    
                    print(outline)
                    outfile.write(outline)
                
            outfile.close()
            
        except:
            outline = ('Error occurred applying model to sample {}\n').format(sample)
            print(outline)
            1/0
                    
def calculate_score(value_df):
        
    #calc score
    temp_df = value_df.drop(['length', 'onehot'], axis=1)
    value_df['total_score'] = (temp_df.sum(axis=1))
            
    for col_name in value_df.columns:        
        if col_name not in ['length', 'onehot']:
            
            new_evl_name = col_name + '_eval'
            new_zsc_name = col_name + '_zscore'
            new_pval_name = col_name + '_zscore_pval'
            g_median = value_df[col_name].mean()
            g_std = value_df[col_name].std()
            
            eval_case = (value_df[col_name] > (g_median + g_std))
            zscore = (value_df[col_name] - g_median) / max(g_std, 1)
            
            value_df[new_evl_name] = eval_case
            value_df[new_zsc_name] = zscore
            value_df[new_pval_name] = stats.norm.sf(abs(zscore))
    
    return(value_df)

    
def value_log_file(sample, trace_dict, defined_dict, region_to_feature_dict):
    ##TODO change values    
    pickle_out = ("{}/{}_trace_dict.p").format(pickles_dir, sample)
    print(pickle_out)
    pickle.dump(trace_dict, open(pickle_out, 'wb'))
        
    pickle_out = ("{}/{}_defined_dict.p").format(pickles_dir, sample)
    print(pickle_out)
    pickle.dump(defined_dict, open(pickle_out, 'wb'))
    
    pickle_out = ("{}/{}_region_to_feature_dict.p").format(pickles_dir, sample)
    print(pickle_out)
    pickle.dump(region_to_feature_dict, open(pickle_out, 'wb')) 
    
    value_dict = {}
    
    uid_list = []
    abs_list = []
    for chromo in trace_dict:
        for region in trace_dict[chromo]:
            uid = ('{}:{}-{}').format(chromo, min(region), max(region))
            uid_list.append(uid)
            
            if region in region_to_feature_dict[chromo]:
                value_dict[uid] = region_to_feature_dict[chromo][region]
                
                for read_type in ['rd', 'disco', 'split']:
                    for feature in trace_dict[chromo][region][read_type].keys():
                        feature_name = ('{}_{}').format(read_type, feature)
                        value_dict[uid][feature_name] = trace_dict[chromo][region][read_type][feature]
                            
                if args.known_mge:
                    if region in defined_dict[chromo]:
                        value_dict[uid]['onehot'] = int(defined_dict[chromo][region])
                    else:
                        print('missing target value: ', chromo, region)
                else:
                    target_value = 'unknown'
                
            else:
                abs_list.append([chromo, region])
                
    value_df = pd.DataFrame.from_dict(value_dict, orient='index')
    
    value_df = value_df.fillna(0)
    
    value_df = calculate_score(value_df)
    outfile_name = ('{}/{}_feature_values.tab').format(final_output_dir, sample)
    value_df.to_csv(outfile_name)
    
    pickle_out = ('{}/{}_feature_values.p').format(final_output_dir, sample)
    print(pickle_out)
    pickle.dump(value_df, open(pickle_out, 'wb')) 
    
    return(value_df)

def calc_percentile(score, percentile_array):
    index = 0
    for pct in percentile_array:
        
        if score < pct:
            return(index)
        
        index+=10
    return(index)

def anchor_properties(anchor_dict):
    temp_dict = {}
    
    return_line = 'anchors;'
    
    for num in anchor_dict:
        (locus, ct) = anchor_dict[num]
        if ct not in temp_dict:
            temp_dict[ct] = set()
        
        temp_dict[ct].add(locus)
        
    total = sum(temp_dict.keys())
    
    slist = sorted(temp_dict.keys())
    
    for ct in slist:
        for locus in temp_dict[ct]:
            outline = ('name={}, pct={}, ct={};').format(locus, 
            round(100*ct/total), 
            ct)
            return_line += outline
            
    return(return_line[:-1], temp_dict)               
            

def make_mge_gff(scored_uid_df, feature_dict):
    scored_uid_dict = scored_uid_df.to_dict(orient='index')
    
    global mge_region_to_name_dict

    outfile_name = ('{}/{}_candidate_TEa.gff').format(final_output_dir, each_sample)
    outfile = open(outfile_name, 'w')
    
    # pre_percentile_array = []
    
    # for uid in scored_uid_dict:
    #     for region in scored_uid_dict[chromo]:
    #         pre_percentile_array.append(scored_uid_dict[chromo][region])
            
    percentile_array = np.percentile(np.array(scored_uid_df['total_score_zscore']), range(0,101,10))
           
    for uid in scored_uid_dict:
        chromo, regions = uid.split(':')
        start = int(regions.split('-')[0])
        stop = int(regions.split('-')[1])
        
        region = (start, stop)

        mge_name = mge_region_to_name_dict[chromo][region]
        
        if mge_name in feature_dict:
            anchor_line, anchor_dict = anchor_properties(feature_dict[mge_name]['anchors'])
        else:
            anchor_line, anchor_dict = 'no_significant_anchors', {}
            
        score = scored_uid_dict[uid]['total_score_zscore']
        pval = scored_uid_dict[uid]['total_score_zscore_pval']
        
        pct = calc_percentile(score, percentile_array)
        
        #chrXVI	ps_sgd	ltr	62389	62720	.	.	.	ID=S000007192_YPLWdelta5
        outline = ('{chromo}\tmge\tcandidate_TEa\t{start}\t{stop}\t{score}'
                   '\t.\t.\t{name}; percentile={pct}; zscore_pval={pval}; {anchor_line}\n').format(
                       chromo = chromo,
                       start = start, stop = stop,
                       pct = pct, name = mge_name, score = score,
                       pval = pval, anchor_line = anchor_line)
        
        outfile.write(outline)
        
        for ct in anchor_dict:
            for locus in anchor_dict[ct]:
                achromo, nts = locus.split(':')
                astart, astop = nts.split('-')
                
                outline = ('{chromo}\tmge\tcandidate_anchor\t{start}\t{stop}\t{ct}'
                           '\t.\t.\t{name}_{locus}; score={score}; zscore_pval={pval}; percentile={pct}; {anchor_line}\n').format(
                               chromo = achromo,
                               start = astart, stop = astop,
                               ct = ct, pct = pct, name = mge_name, 
                               locus = locus, score = score,
                               pval = pval, anchor_line = anchor_line)
                
                outfile.write(outline)
                
    outfile.close()


def mge_trace_builder(defined_dict, region_to_feature_dict, feature_dict):
    
    pickle_name = ('{}/mpileup_{}_RD.p').format(pickles_dir, each_sample)
    rd_mpileup_df = pickle_loader(pickle_name, 'df')
    
    pickle_name = ('{}/mpileup_{}_discordant.p').format(pickles_dir, each_sample)
    discordant_mpileup_df = pickle_loader(pickle_name, 'df')  
    
    pickle_name = ('{}/mpileup_{}_split.p').format(pickles_dir, each_sample)
    split_mpileup_df = pickle_loader(pickle_name, 'df')  
              
    trace_dict = mge_trace_three_frames(discordant_mpileup_df, rd_mpileup_df, split_mpileup_df)
    
    scored_uid_df = value_log_file(each_sample, trace_dict, defined_dict, region_to_feature_dict)
                
    make_mge_gff(scored_uid_df, feature_dict)

def count_uids_in_region(start, stop, nt_dict):
    uid_set = set()
    
    for nt in range(start, stop + 1):
        if nt in nt_dict:
            for uid in nt_dict[nt]:
                uid_set.add(uid)
                
    return(len(uid_set))

def mge_combine_regions(nt_deets):
    start, stop = min(nt_deets), max(nt_deets) 
    
    #region_set collects paired start, stops (tuple)
    region_set = set()
    #covered_set records which nts have already been assigned
    covered_set = set()
    
    # region_start should always be the min(nt_deets) not in covered_set
    for start_nt in range(start, stop + 1):        
        if (start_nt not in covered_set) and (start_nt in nt_deets):

            region_start = start_nt
            region_stop = start_nt
            
            gap_ct = 0
            
            for scan_nt in range(region_start, stop + 1):
                
                if gap_ct < gap:
                    if (scan_nt not in covered_set):
                        covered_set.add(scan_nt)
                        
                    if scan_nt in nt_deets:
                        region_stop = scan_nt
                        gap_ct = 0
                
                gap_ct += 1
                
            region_set.add((region_start, region_stop))
    
    return(region_set)

def build_edge_read_set(discordant_sam_file, split_sam_file):
    global mge_nt_lookup_dict, mge_name_to_region_dict
    
    feature_dict = {}
    selected_uid_dict = {}
    
    for runmode in ['discordant', 'split']:        
        if runmode == 'discordant':
            print('Parsing discordant reads...')
            read_file_name = discordant_sam_file        
            
        if runmode == 'split':
            print('Parsing split reads...')
            read_file_name = split_sam_file
            
        read_file = open(read_file_name)

        for line in read_file:
            line = line.strip()
            
            if line[0]!='@':                
                uid = line.split('\t')[0]
                chromo = line.split('\t')[2]
                
                if chromo in mge_nt_lookup_dict:
                    start = int(line.split('\t')[3])
                    mapq = int(line.split('\t')[4])
                    
                    if mapq >= mapq_val:
                        cigar = line.split('\t')[5]
                        
                        try:
                            stop = start + parse_cigar(cigar, 'last_match')
                        except:
                            print('parse_cigar error')
                            print(start, cigar, parse_cigar(cigar, 'last_match'))
                            1/0
    
                        for nt in range(start, stop + 1):
                            if nt in mge_nt_lookup_dict[chromo]: 
                                mge_name_set = mge_nt_lookup_dict[chromo][nt]
                                                                
                                for mge_name in mge_name_set:                                    
                                    if mge_name not in feature_dict:
                                        feature_dict[mge_name] = {
                                            'chromo' : chromo,
                                            'mge_region': set(), 
                                            'anchor_loci': {},
                                            'uid': set(),
                                            'mge_mapq':[],
                                            'anchors': {},
                                            'anchor_mapq':[]}
                                        
                                    feature_dict[mge_name]['mge_region'] = mge_name_to_region_dict[mge_name]
                                    feature_dict[mge_name]['uid'].add(uid)
                                    feature_dict[mge_name]['mge_mapq'].append(mapq)
                                    
                                    
                                    if uid not in selected_uid_dict:
                                        selected_uid_dict[uid] = set()
                                    
                                    selected_uid_dict[uid].add(mge_name)
                                                                                    
        read_file.close()
        
        print('Reparsing ', read_file_name)
        
        read_file = open(read_file_name)

        for line in read_file:
            line = line.strip()
            
            if line[0]!='@':                
                uid = line.split('\t')[0]
                
                if uid in selected_uid_dict:
                    chromo = line.split('\t')[2]
                    start = int(line.split('\t')[3])
                    mapq = int(line.split('\t')[4])
                    cigar = line.split('\t')[5]
                                        
                    #filter by qscore should suppress non-unique mapping
                    if mapq >= mapq_val:
                        try:
                            stop = start + parse_cigar(cigar, 'last_match')
                        except:
                            print('parse_cigar error')
                            print(start, cigar, parse_cigar(cigar, 'last_match'))
                            1/0
                            
                        for mge_name in selected_uid_dict[uid]:
                            mge_chromo = feature_dict[mge_name]['chromo']
                            mge_region = feature_dict[mge_name]['mge_region']
                                                    
                            #edge case to make sure we aren't counting those reads that self map.
                            process = True
                            if chromo == mge_chromo:
                                for nt in range((start - gap), (stop + gap +1)):
                                    if nt in range(min(mge_region), max(mge_region)):
                                        process = False                           
                            
                            if process:
                                feature_dict[mge_name]['anchor_mapq'].append(mapq)
                                
                                if chromo not in feature_dict[mge_name]['anchor_loci']:
                                    feature_dict[mge_name]['anchor_loci'][chromo] = {}
                                    
                                for nt in range(start, stop + 1):
                                    if nt not in feature_dict[mge_name]['anchor_loci'][chromo]:
                                        feature_dict[mge_name]['anchor_loci'][chromo][nt] = set()
                                        
                                    feature_dict[mge_name]['anchor_loci'][chromo][nt].add(uid)
                                        
    print('Associating reads ... ')
                            
    for mge_name in feature_dict:
        for chromo in feature_dict[mge_name]['anchor_loci']:
            nt_deets = feature_dict[mge_name]['anchor_loci'][chromo]
            
            if min(nt_deets) + gap >= max(nt_deets):
                num = len(feature_dict[mge_name]['anchors'])
                locus = ('{}:{}-{}').format(chromo, min(nt_deets), max(nt_deets))
                ct = count_uids_in_region(min(nt_deets), 
                                          max(nt_deets), 
                                          feature_dict[mge_name]['anchor_loci'][chromo])
                feature_dict[mge_name]['anchors'][num] = (locus, ct)                
            else:
                region_set = mge_combine_regions(nt_deets)
                
                for region in region_set:
                    num = len(feature_dict[mge_name]['anchors'])
                    locus = ('{}:{}-{}').format(chromo, min(region), max(region))
                    ct = count_uids_in_region(min(nt_deets), 
                                          max(nt_deets), 
                                          feature_dict[mge_name]['anchor_loci'][chromo])
                
                feature_dict[mge_name]['anchors'][num] = (locus, ct)
                
    region_to_feature_dict = {}
    
    for mge_name in feature_dict:
        chromo = feature_dict[mge_name]['chromo']
        mge_region = feature_dict[mge_name]['mge_region']

        median_mge_mapq = np.median(feature_dict[mge_name]['mge_mapq'])        
        mean_mge_mapq = np.mean(feature_dict[mge_name]['mge_mapq'])
        std_mge_mapq = np.std(feature_dict[mge_name]['mge_mapq'])

        median_anchor_mapq = np.median(feature_dict[mge_name]['anchor_mapq'])
        mean_anchor_mapq = np.mean(feature_dict[mge_name]['anchor_mapq'])
        std_anchor_mapq = np.std(feature_dict[mge_name]['anchor_mapq'])
        
        if chromo not in region_to_feature_dict:
            region_to_feature_dict[chromo] = {}
            
        if mge_region not in region_to_feature_dict[chromo]:
            region_to_feature_dict[chromo][mge_region] = {'rel_anchor':0,
                                                          'ttl_anchor':0,
                                                          'mge':0}
        read_ct = 0
        max_read_ct = 0
        for num in feature_dict[mge_name]['anchors']:
            num_reads = feature_dict[mge_name]['anchors'][num][1]
            read_ct += num_reads
            
            if num_reads > max_read_ct:
                max_read_ct = num_reads
                
        region = feature_dict[mge_name]['mge_region']
        length = max(max(region) - min(region),1)
        
        region_to_feature_dict[chromo][mge_region]['median_mge_mapq'] = median_mge_mapq
        region_to_feature_dict[chromo][mge_region]['mean_mge_mapq'] = mean_mge_mapq
        region_to_feature_dict[chromo][mge_region]['std_mge_mapq'] = std_mge_mapq
        region_to_feature_dict[chromo][mge_region]['median_anchor_mapq'] = median_anchor_mapq
        region_to_feature_dict[chromo][mge_region]['mean_anchor_mapq'] = mean_anchor_mapq
        region_to_feature_dict[chromo][mge_region]['std_anchor_mapq'] = std_anchor_mapq
        region_to_feature_dict[chromo][mge_region]['length'] = length
        region_to_feature_dict[chromo][mge_region]['rel_anchor'] = (max_read_ct**2) / (max(read_ct,1) * length)
        region_to_feature_dict[chromo][mge_region]['rel_ttl_anchor'] = read_ct/length
        region_to_feature_dict[chromo][mge_region]['rel_mge'] = len(feature_dict[mge_name]['uid'])/length
        region_to_feature_dict[chromo][mge_region]['anchor'] = (max_read_ct**2) / (max(read_ct,1))
        region_to_feature_dict[chromo][mge_region]['ttl_anchor'] = read_ct
        region_to_feature_dict[chromo][mge_region]['mge'] = len(feature_dict[mge_name]['uid'])


            
    return(feature_dict, region_to_feature_dict)
                        
""" Step Three """
if args.mge_file:
    """
    #    
    """
    resource_dict = io_load()
    read_type_list = resource_dict['read_types']
    discordant_sam_file = resource_dict['discordant_sam_file']
    split_sam_file = resource_dict['split_sam_file']
        
    each_sample = resource_dict['run_name']
    
    chromo_size_p = pickle_loader(resource_dict['chromo_size'], 'dict')
    
    chromo_list = []
    for chromo in chromo_size_p:
        chromo_list.append(chromo)     
        
    if args.CNV_min_length:
        min_length = int(args.CNV_min_length)
    else:
        min_length = 300
        
    if args.max_gap:
        gap = int(args.max_gap)
    else:
        gap = 300
                        
    file_name_lookup = []
    depth_dict = {}
    
    defined_dict = {}
    
    if args.known_mge:
        infile = open(args.known_mge)
        
        for line in infile:
            #{chromo}\tmugio\tcandidate_TE\t'
            #               '{start}\t{stop}\t{tea}\t.\t.\t'
            #               'ID={name}; score={score}; TEa={tea}; is uid={is_uid}
            
            line = line.strip()
            chromo, _a, _b, start, stop, tea, _sign, _dot, deets = line.split('\t')
            
            tea = bool(int(tea))
            
            region = (int(start), int(stop))
            
            if chromo not in defined_dict:
                defined_dict[chromo] = {}
                
            if region not in defined_dict[chromo]:
                defined_dict[chromo][region] = tea
                
            else:
                if defined_dict[chromo][region] == False:
                    defined_dict[chromo][region] = tea
    
    feature_dict, region_to_feature_dict = build_edge_read_set(discordant_sam_file, split_sam_file)
            
    mge_trace_builder(defined_dict, region_to_feature_dict, feature_dict)

        
        
""" Step Four"""
if args.map_reads:
    #load resources:
    resource_dict = io_load()
    
    chromo_size_pickle = resource_dict['chromo_size']
    fastq_1_filename = resource_dict['fastq_1']
    fastq_2_filename = resource_dict['fastq_2']
    
    bam_file = resource_dict['bam_file']
    sam_file = resource_dict['sam_file']
    
    discordant_sam_file = resource_dict['discordant_sam_file']
    split_sam_file = resource_dict['split_sam_file']
    #soft_sam_file = resource_dict['soft_sam_file']
    
    read_type_list = resource_dict['read_types']
    read_type_list.remove('RD')
        
    with open(chromo_size_pickle, 'rb') as fp:
        chromo_size_dict = pickle.load(fp)
    
    score = (resource_dict['split_zscore_median'] + resource_dict['split_zscore_std']) * sw
    outline = ('Using split z-score: {score}\nSplit median: {median}\t'
               'Split standard deviation: {std}\nSplit weight {weight}').format(
        score = score, 
        median = resource_dict['split_zscore_median'],
        std = resource_dict['split_zscore_std'],
        weight = sw)
                   
    print(outline)
                
    if 'discordant_zscore_median' in resource_dict:
        disco_score = (resource_dict['discordant_zscore_median'] + resource_dict['discordant_zscore_std']) * dw
        
        outline = ('Using discordant z-score: {score}\nDisco median: {median}\t'
                   'Disco standard deviation: {std}\nDisco weight {weight}').format(
            score = disco_score, 
            median = resource_dict['discordant_zscore_median'],
            std = resource_dict['discordant_zscore_std'],
            weight = dw)
                       
        print(outline)
    
    #initialize variables
    fastq_1_list = []
    fastq_2_list = []
    uid_to_seq_dict = {}
    uid_deets_dict = {}
    uid_source_dict = {}
    chromo_walk_dict = {}
    uniuid_to_seq_dict = {}
    uniuid_to_phred_dict = {}
    ct = 0

    alignment_map = {}
    hypothesis_map = {}
    chromo_map = {}
    
    if not args.min_coverage:
        min_coverage = 3
    else:
        min_coverage = float(args.min_coverage)
                            
    def correct_size(chromo, start, stop, chromo_size_dict):
        #This function adds a correction factor (flanking) when possible
        correction = int(flank_nt)
        if chromo in chromo_size_dict:
            if (start - correction) > 0:
                new_start = (start - correction)
            else:
                new_start = 1

            if (stop + correction) < int(chromo_size_dict[chromo]):
                new_stop = stop + correction
            else:
                new_stop = int(chromo_size_dict[chromo])-1
                
            if new_stop <= new_start:
                print(line)
                1/0
        return(new_start, new_stop)
        
    def parse_deets(deets):
        """ This function generates a stop codon by scraping start and cigar 
        values and passing them to parse_cigar, and then correcting the size to 
        enable flanking.
        """
        chromo, start, qscore, cigar = deets.split('\t')
        start = int(start)
        start, stop = correct_size(chromo, start, (start + parse_cigar(cigar, 'first_match')), chromo_size_dict)
        return(chromo, start, stop, qscore, cigar)
            
    def find_hypo(locus_to_hypo_lookup, chromo, start, stop):
        for each_nt in range(start, stop+1):
            region_name = str(chromo + '_' + str(each_nt))

            if region_name in locus_to_hypo_lookup:
                hypo_id = locus_to_hypo_lookup[region_name]
                return(hypo_id)
                
        #if that doesn't work:
        print('error attempting to map to unassigned region')
        print(locus_to_hypo_lookup)
        print(chromo, start, stop)
        
        1/0
            
    def filter_regions(chromo, start, stop, region_filter_dict, ancestor_filter):
        process_brks = True
        if chromo in region_filter_dict:
            region_set = region_filter_dict[chromo]
            for nt in range((start-flank_nt), (stop+flank_nt+1)):
                if nt in region_set:
                    process_brks=False
                    ancestor_filter+=1
                    return(process_brks, ancestor_filter)
        return(process_brks, ancestor_filter)
    
    def break_walk(chromo_walk_dict, chromo, start, stop):
        #This function updates the chromo walk dict with the start and stop of every breakpoint associated read.
        if chromo in chromo_walk_dict:
            update = False
            
            aligned_set = chromo_walk_dict[chromo]
            
            for nt in range(start, stop+1):
                if nt not in aligned_set:
                    aligned_set.add(nt)
                    update = True
                    
            if update:
                chromo_walk_dict[chromo] = aligned_set
        
        if chromo not in chromo_walk_dict:
            aligned_set = set()
            for nt in range(start, stop+1):
                aligned_set.add(nt)
                    
            chromo_walk_dict[chromo] = aligned_set
    
        return(chromo_walk_dict)
        
    def load_peaks(peaks_file_name):
        peaks_file = open(peaks_file_name)
        
        peaks_dict = {}
        for line in peaks_file:
            #'{chromo}\tRead_Depth\tCNV\t{start}\t{stop}\t{rel_c}\t.\t.\tID={chromo}:{start}-{stop};rel_chromosome_RD={rel_c};rel_genome_RD={rel_g};sample={sample}\n'
            chromo = line.split('\t')[0]
            start = int(line.split('\t')[3])
            stop = int(line.split('\t')[4])

            #section handles flanking and chromosome ends
            chromo_length = int(chromo_size_dict[chromo])
            if start-flank_nt > 0:
                start -= flank_nt
            else:
                start = 0
            
            if stop + flank_nt < chromo_length:
                stop += flank_nt
            else:
                stop = chromo_length
            
            #section handles loading range dictionary
            if chromo in peaks_dict:
                loci_set = peaks_dict[chromo]
                
                for nt in range(start,stop+1):
                    loci_set.add(nt)
                    
                peaks_dict[chromo]=loci_set
                
            else:
                loci_set = set()
                
                for nt in range(start,stop+1):
                    loci_set.add(nt)
                    
                peaks_dict[chromo]=loci_set
            
        return(peaks_dict)
        
    def map_alignments(uni_uid, uid_align, alignment_map, hypothesis_map, ancestor_filter):
        uid = uni_uid.split('.')[0]
        process_brks = True
        chromo, start, stop, qscore, cigar, instance = uid_align
        uni_uid_instance = ('{}~{}').format(uni_uid, instance)
                                        
        if args.filter_object:
            process_brks, ancestor_filter = filter_regions(chromo, start, stop, region_filter_dict, ancestor_filter)

        if process_brks:
            if uid not in hypothesis_map:
                hypothesis_map[uid] = [uni_uid_instance]
            else:
                hypothesis_map[uid] += [uni_uid_instance]
            
            if uni_uid_instance not in alignment_map:
                alignment_map[uni_uid_instance] = [chromo, start, stop, qscore, cigar]
            else:
                p_chromo = uid_align[0]
                p_start = uid_align[1]
                p_stop = uid_align[2]
                if (p_chromo != chromo) or (p_start != start) or (p_stop != stop):
                    print('uni_uid_instance collision', uni_uid_instance, uid_align, alignment_map[uni_uid_instance])
        
        return(alignment_map, hypothesis_map, ancestor_filter)
    
    def refine_alignment(uni_uid, alignments, alignment_map, runmode, locus_to_hypo_lookup, refined_map):
        q_uid = ('{}_{}').format(runmode, uni_uid)
        
        if runmode == 'discordant':
            if len(alignments) == 2:
                hypo_id_list = []

                for each_instance in alignments:
                    chromo, start, stop, qscore, cigar = alignment_map[each_instance]
                    hypo_id_list.append(find_hypo(locus_to_hypo_lookup, chromo, start, stop))
                                
                hypo_id_list.sort()
                breakpoint_hypothesis = ('{}_{}').format(hypo_id_list[0],hypo_id_list[1])
                                
                if breakpoint_hypothesis in refined_map:
                    refined_map[breakpoint_hypothesis] += [q_uid]
                else:
                    refined_map[breakpoint_hypothesis] = [q_uid]
                    
        if runmode == 'split':
            if len(alignments) > 1:
                temp_uni_uid_dict = {}
    
                for each_instance in alignments:
                    chromo, start, stop, qscore, cigar = alignment_map[each_instance]
                    hypo_id = find_hypo(locus_to_hypo_lookup, chromo, start, stop)
                    
                    uni_uid, instance = each_instance.split('~')

                    if uni_uid not in temp_uni_uid_dict:
                        temp_uni_uid_dict[uni_uid] = [hypo_id]
                    else:
                        temp_uni_uid_dict[uni_uid] += [hypo_id]
                
                for uni_uid, hypo_id_list in temp_uni_uid_dict.items():
                    for x in range(len(hypo_id_list)-1):
                        for y in range(x+1,len(hypo_id_list)):
                            id_list = [hypo_id_list[x], hypo_id_list[y]]
                            id_list.sort()
                            breakpoint_hypothesis = ('{}_{}').format(hypo_id_list[0],hypo_id_list[1])
                
                            if breakpoint_hypothesis in refined_map:
                                refined_map[breakpoint_hypothesis] += [q_uid]
                            else:
                                refined_map[breakpoint_hypothesis] = [q_uid]
                                   
        return(refined_map)
                
    for runmode in read_type_list:        
        if runmode == 'discordant':
            print('Parsing discordant reads...')
            read_file_name = discordant_sam_file        
            peaks_file_name = ('{}/{}_discordant_depth.gff').format(final_output_dir, output_file)
        
        if runmode == 'split':
            print('Parsing split reads...')
            read_file_name = split_sam_file
            peaks_file_name = ('{}/{}_split_depth.gff').format(final_output_dir, output_file)
        
        #filter_by_peaks:
        peaks_dict = load_peaks(peaks_file_name)
            
        ct = 0
        uni_uid_set = set()
        process_list = []
        read_file = open(read_file_name)

        for line in read_file:
            process = False
            ct+=1
            
            if line[0]!='@':                
                uid = line.split('\t')[0]
                flag = line.split('\t')[1]
                chromo = line.split('\t')[2]
                start = int(line.split('\t')[3])
                qscore = int(line.split('\t')[4])
                cigar = line.split('\t')[5]
                
                try:
                    stop = start + parse_cigar(cigar, 'last_match')
                except:
                    print('parse_cigar error')
                    print(start, cigar, parse_cigar(cigar, 'last_match'))
                    1/0
                                    
                if qscore >= mapq_val:
                    process = False
                    if chromo in peaks_dict:
                        loci_set = peaks_dict[chromo]
                        
                        for nt in range(start,stop):
                            if nt in loci_set:
                                process = True
                    
                if process:
                    """ here we make the origin of the fastq file explicit, first asking is the reads are paired,
                    then asking if they are from the first ('forward') or second ('reverse') fastq files.
                    NOTE: I haven't tested this with interlaced fastq, and I'm not sure how it would work.
                    """                    
                    if np.unpackbits(np.array([flag], dtype=np.uint8))[-1]==1:
                        if np.unpackbits(np.array([flag], dtype=np.uint8))[-8]==1:
                            fastq_2_list.append(uid)
                            side = '.2'
                        else:
                            fastq_1_list.append(uid)
                            side = '.1'
                    else:
                        side = '.1'
                        
                    uni_uid = str(uid)+side
                    
                    if uni_uid in uid_deets_dict:
                        instance = len(uid_deets_dict[uni_uid])
                        uid_deets_dict[uni_uid]+=[[chromo, start, stop, qscore, cigar, instance]]
                    else:
                        uid_deets_dict[uni_uid]=[[chromo, start, stop, qscore, cigar, 0]]
                        
                    if runmode not in uid_source_dict:
                        uni_uid_set.add(uni_uid)
                        uid_source_dict[runmode] = uni_uid_set
                    else:
                        uid_source_dict[runmode].add(uni_uid)
                        
        read_file.close()
        
        
    for runmode in read_type_list:
        pull_from_rd_1 = set()
        pull_from_rd_2 = set()
        
        if runmode == 'split' or runmode == 'discordant':            
            uni_uid_set = uid_source_dict[runmode]
            
            for uni_uid in uni_uid_set:
                uid = uni_uid.split('.')[0]
                               
                if uid.count('_') == 1:
                    uid = uid.split('_')[0]
                    
                if uid.count('_') > 1:
                    uid = uid.split('_')[0] 
                    print('uni_uid ... ')
                    print(uid)
                    1/0
                
                if (uid+'.1' not in uid_deets_dict):
                    pull_from_rd_1.add(uid)
                    fastq_1_list.append(uid)
                    
                    
                if (uid+'.2' not in uid_deets_dict):
                    pull_from_rd_2.add(uid)
                    fastq_2_list.append(uid)
                
                    
        rd_file = open(sam_file)
        
        for line in rd_file:
            process_1 = False
            process_2 = False
            
            if line[0]!='@':                
                uid = line.split('\t')[0]
                flag = line.split('\t')[1]
                                    
                if uid in pull_from_rd_1:
                    if np.unpackbits(np.array([flag], dtype=np.uint8))[-8]==0:
                        process_1 = True
                        uni_uid = uid+'.1'
                
                if uid in pull_from_rd_2:
                    if np.unpackbits(np.array([flag], dtype=np.uint8))[-8]==1:
                        process_2 = True
                        uni_uid = uid+'.2'
                        
                if (process_1 or process_2) and (uni_uid not in uid_deets_dict):
                    chromo = line.split('\t')[2]
                    start = int(line.split('\t')[3])
                    cigar = line.split('\t')[5]

                    try:
                        stop = start + parse_cigar(cigar, 'last_match')
                    except:
                        print('parse_cigar error')
                        print(start, cigar, parse_cigar(cigar, 'last_match'))
                        1/0
                        
                    if uni_uid in uid_deets_dict:
                        instance = len(uid_deets_dict[uni_uid])
                        uid_deets_dict[uni_uid]+=[[chromo, start, stop, qscore, cigar, instance]]
                    else:
                        uid_deets_dict[uni_uid]=[[chromo, start, stop, qscore, cigar, 0]]
                        
                    if runmode not in uid_source_dict:
                        uni_uid_set.add(uni_uid)
                        uid_source_dict[runmode] = uni_uid_set
                    else:
                        uid_source_dict[runmode].add(uni_uid)
                            
        rd_file.close() 
        
    if len(fastq_1_list) > 0:
        fastq_1_file = open(fastq_1_filename)
        uniuid_to_seq_dict, uniuid_to_phred_dict = parse_fastq(fastq_1_filename, '.1', fastq_1_list, uniuid_to_seq_dict, uniuid_to_phred_dict)

    if len(fastq_2_list) > 0:
        fastq_2_file = open(fastq_2_filename)
        uniuid_to_seq_dict, uniuid_to_phred_dict = parse_fastq(fastq_2_filename, '.2', fastq_2_list, uniuid_to_seq_dict, uniuid_to_phred_dict)
                      
    pickle_out = ("{}/uniuid_to_seq_dict.p").format(pickles_dir)
    print(pickle_out)
    pickle.dump(uniuid_to_seq_dict, open(pickle_out, 'wb'))
    
    pickle_out = ("{}/uniuid_to_phred_dict.p").format(pickles_dir)
    pickle.dump(uniuid_to_phred_dict, open(pickle_out, 'wb'))
    
    pickle_out = ("{}/uid_source_dict.p").format(pickles_dir)
    pickle.dump(uid_source_dict, open(pickle_out, 'wb'))
    
    print('Mapping regions...')
    """ This component converts the two sides of the split reads into a pair of addresses
    (ie: chromo_1, start_1, stop_1 and chromo_2, start_2, stop_2) if these are not filtered 
    by the filter regions command, they are stored in the hypothesis_map dictionary as a pair
    """
    ancestor_filter = 0
    
    for runmode in read_type_list:
        uni_uid_set = uid_source_dict[runmode]
        
        if runmode == 'discordant':
            print('Running discordant read mode on ' + str(len(uni_uid_set)) + ' reads ...')
                                
        if runmode == 'split':
            print('Running split read mode on ' + str(len(uni_uid_set)) + ' reads ...')
                        
        for uni_uid in uni_uid_set:
            deets = uid_deets_dict[uni_uid]

            if len(deets) > 1:
                for uid_align in deets:
                    alignment_map, hypothesis_map, ancestor_filter = map_alignments(uni_uid, uid_align, alignment_map, hypothesis_map, ancestor_filter)
            else:
                uid_align = deets[0]
                alignment_map, hypothesis_map, ancestor_filter = map_alignments(uni_uid, uid_align, alignment_map, hypothesis_map, ancestor_filter)
    
    outline = ('\t{} reads mapped to {} non-unique locations with {} reads filtered.').format(len(hypothesis_map), len(alignment_map), ancestor_filter)           
    print(outline)
       
    print('Defining regions ...')
    
    """ Here make a chromosome specific set of nucleotides that have an alignment on them."""
    for uni_uid_instance, deets in alignment_map.items():
        chromo, start, stop, qscore, cigar = deets
        chromo_walk_dict = break_walk(chromo_walk_dict, chromo, int(start), int(stop))      
            
    """ Here we take those discrete nucleotides and connect them into regions. 
    the hypo_to_locus_lookup allows the look up of the location each hypothesis refers to, 
    while the locus_to_hypo_lookup performs the inverse. """
    
    hypo_to_locus_lookup = {}
    hypo_ct = 0
    locus_to_hypo_lookup = {} 
    
    for each_chromo, aligned_set in chromo_walk_dict.items():               
        if len(aligned_set) > 0:
            min_chromo = min(aligned_set)
            max_chromo = max(aligned_set)
            print('\tProcessing '+ each_chromo + '...')

            start, stop, open_range = 0, 0, False
            
            for nt in range(min_chromo, max_chromo+1):
                if nt in aligned_set:                       
                    if not open_range:
                        open_range = True
                        start = nt
                        hypo_ct += 1
                            
                if open_range and ((nt+1) not in aligned_set):
                    open_range = False
                    stop = nt
                    hypo_to_locus_lookup[hypo_ct] = (each_chromo, start, stop)
                        
                    for each_nt in range(start, stop):
                        range_name = str(each_chromo + '_' + str(each_nt))
                        locus_to_hypo_lookup[range_name] = hypo_ct
        else:
            print('error in aligned set')
            print(each_chromo, aligned_set)
            1/0

    print('Refining hypotheses ...')                 
    refined_map = {}
    uids_that_were_filtered = []
    
    for runmode in read_type_list:
        uni_uid_set = uid_source_dict[runmode]
        
        for uni_uid in uni_uid_set:
            uid = uni_uid.split('.')[0]
            if uid in hypothesis_map:
                alignments = hypothesis_map[uid]
                refined_map = refine_alignment(uni_uid, hypothesis_map[uid], alignment_map, runmode, locus_to_hypo_lookup, refined_map)
            else:
                uids_that_were_filtered.append(uid)
                            
    #define outputs
    break_tab = ('{}/break_bt.tab').format(temp_dir)
    break_bed = ('{}/break_bd.bed').format(temp_dir)
    dbreak_tab = ('{}/dbreak_dbt.tab').format(temp_dir)
    dbreak_bed = ('{}/dbreak_dbd.bed').format(temp_dir)
    #
    break_tab_file = open(break_tab,'w')
    break_bed_file = open(break_bed,'w')
    dbreak_tab_file = open(dbreak_tab, 'w')
    dbreak_bed_file = open(dbreak_bed, 'w')
    
    refined_ct = 0
    dbreak_ct = 0
    
    for breakpoint_hypothesis, uid_list in refined_map.items():
        process_brks_1, process_brks_2 = True, True
        
        hypo_side_left = int(breakpoint_hypothesis.split('_')[0])
        hypo_side_right = int(breakpoint_hypothesis.split('_')[1])
        
        chromo_l, start_l, stop_l = hypo_to_locus_lookup[hypo_side_left]
        chromo_r, start_r, stop_r = hypo_to_locus_lookup[hypo_side_right]
        
        if args.filter_object:
            process_brks_1, ancestor_filter = filter_regions(chromo_l, start_l, stop_l, region_filter_dict, ancestor_filter)
            process_brks_2, ancestor_filter = filter_regions(chromo_r, start_r, stop_r, region_filter_dict, ancestor_filter)
                
        if process_brks_1 and process_brks_2:
            
            
            if len(uid_list) > 1:
                split_ct = 0
                disco_ct = 0
                #soft_ct = 0
                processed_list = []
                uid_str = ''
                
                temp_uid_dict = {}

                for each in uid_list:
                    if each.count('_') == 1:
                        runmode, uid = each.split('_')
                        
                    else:
                        runmode, uid = each.split('_', 1)
                        
                    if uid not in temp_uid_dict:
                        temp_uid_dict[uid] = set()
                        
                    temp_uid_dict[uid].add(runmode)
                                            
                for uid, runmodes in temp_uid_dict.items():
                    if 'split' in runmodes:
                        processed_list.append(uid)
                        split_ct += 1
                        uid_str += ('split_{},').format(uid)
                        
                    if 'discordant' in runmodes:
                        disco_ct += 1
                        uid_str += ('discordant_{},').format(uid)
                                        
                uid_str = uid_str.rsplit(',',1)[0]
                
                if len(temp_uid_dict) >= 2:                    
                    s_weight = sw * split_ct
                    d_weight = dw * disco_ct
                    
                    if args.split_only:
                        combined_weight = s_weight
                    else:
                        combined_weight = (d_weight + s_weight)
                                            
                    """generate two file sets - one set contains those hypothetical breaks with split read support
                    while the other set only has significant values of discordant reads available"""
        
                    if ((s_weight > 0) and (combined_weight > score)):
                        refined_ct += 1
                        if hypo_side_left != hypo_side_right:                
                            outline = ('{}_left\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n').format(breakpoint_hypothesis, chromo_l, start_l, stop_l, chromo_r, start_r, stop_r, combined_weight, uid_str)
                            break_tab_file.write(outline)
                            outline = ('{}\t{}\t{}\t{}_left\t{}\t.\n').format(chromo_l, start_l, stop_l, breakpoint_hypothesis, combined_weight)
                            break_bed_file.write(outline)
                            
                            outline = ('{}_right\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n').format(breakpoint_hypothesis, chromo_r, start_r, stop_r, chromo_l, start_l, stop_l, combined_weight, uid_str)
                            break_tab_file.write(outline)
                            outline = ('{}\t{}\t{}\t{}_right\t{}\t.\n').format(chromo_r, start_r, stop_r, breakpoint_hypothesis, combined_weight)
                            break_bed_file.write(outline)
                            
                        else:
                            outline = ('{}_both\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n').format(breakpoint_hypothesis, chromo_l, start_l, stop_l, chromo_l, start_l, stop_l, combined_weight, uid_str)
                            break_tab_file.write(outline)
                            outline = ('{}\t{}\t{}\t{}_both\t{}\t.\n').format(chromo_l, start_l, stop_l, breakpoint_hypothesis, combined_weight)
                            break_bed_file.write(outline)
                            
                    if ((split_ct == 0) and d_weight > 3*disco_score):
                        dbreak_ct += 1
                        if hypo_side_left != hypo_side_right:                
                            outline = ('{}_left\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n').format(breakpoint_hypothesis, chromo_l, start_l, stop_l, chromo_r, start_r, stop_r, combined_weight, uid_str)
                            dbreak_tab_file.write(outline)
                            outline = ('{}\t{}\t{}\t{}_left\t{}\t.\n').format(chromo_l, start_l, stop_l, breakpoint_hypothesis, combined_weight)
                            dbreak_bed_file.write(outline)
                            
                            outline = ('{}_right\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n').format(breakpoint_hypothesis, chromo_r, start_r, stop_r, chromo_l, start_l, stop_l, combined_weight, uid_str)
                            dbreak_tab_file.write(outline)
                            outline = ('{}\t{}\t{}\t{}_right\t{}\t.\n').format(chromo_r, start_r, stop_r, breakpoint_hypothesis, combined_weight)
                            dbreak_bed_file.write(outline)
                            
                        else:
                            outline = ('{}_both\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n').format(breakpoint_hypothesis, chromo_l, start_l, stop_l, chromo_l, start_l, stop_l, combined_weight, uid_str)
                            dbreak_tab_file.write(outline)
                            outline = ('{}\t{}\t{}\t{}_both\t{}\t.\n').format(chromo_l, start_l, stop_l, breakpoint_hypothesis, combined_weight)
                            dbreak_bed_file.write(outline)
            
    break_tab_file.close()
    break_bed_file.close()
    dbreak_tab_file.close()
    dbreak_bed_file.close()    

    outline = ('\t{} significant split hypotheses have been resolved to {} split read supported breakpoints.').format(len(hypothesis_map), refined_ct)
    print(outline)
    if dbreak_ct > 0:
        outline = ('\tIncluding {} discordant read supported breakpoints.').format(dbreak_ct)
        print(outline)
            
    pickle_out = ("{}/refined_map.p").format(pickles_dir)
    pickle.dump(refined_map, open(pickle_out, 'wb'))
            
""" Step Five """
if args.make_filter:
    
    if args.filter_gff:
        gff_file = set(args.filter_gff)

        for each_gff_file in args.filter_gff:
            parse_filter_gff(each_gff_file)
                            
    if args.filter_bed:        
        for each in args.filter_bed:
            if '/' in each:
                source_name = each.rsplit('/',1)[1]
            else:
                source_name = each
                
            parse_bed(each, source_name)
    
    print('Number of filtered elements:')
    print(len(region_name_dict))

    is_list = []
    for each, ct in region_name_dict.items():
        if each not in is_list:
            is_list.append(each)
    
    is_list.sort()
    
    for each in is_list:
        ct = region_name_dict[each]
        outline = ('\t{}:\t{}').format(each,ct)
        print(outline)
                
    pickle_out = ("{}").format(args.filter_object)
    pickle.dump(filter_region_dict, open(pickle_out, 'wb'))
        
""" Step Six """
if args.build_sequence:
    if args.hashsize:
        hashsize = args.hashsize
    else:
        hashsize = 15

    if args.min_contig:
        min_contig = args.min_contig
    else:
        min_contig = 25 
        
    if args.max_divergence:
        max_divergence = args.max_divergence
    else:
        max_divergence = 0.05
        
    if args.cov_cutoff:
        cov_cutoff = args.cov_cutoff
    else:
        cov_cutoff = 1
        
    if args.max_gap:
        gap = int(args.max_gap)
    else:
        gap = 5
        
    complete_qname_dict = {}
    qname_lookup_dict = {}

    locus_lookup_dict={}
    disco_outlines = []
    contig_seq_dict = {}
    gff_list_dict = {}
    gff_rank_dict = {}
    blast_read_aligned_sections_dict = {}
    blast_genome_aligned_regions_dict = {}
    gff_uid_dict ={} 
    best_uid_dict = {}
    
    resource_dict = io_load()
    fa_file = resource_dict['genome_fa']
    fastq_1 = resource_dict['fastq_1']
    fastq_2 = resource_dict['fastq_2']
    bam_file = resource_dict['bam_file']
    
    uniuid_to_seq_dict, uniuid_to_phred_dict = io_load('fq')
        
    break_tab = ('{}/break_bt.tab').format(temp_dir)   
    dbreak_tab = ('{}/dbreak_dbt.tab').format(temp_dir)
                               
    """parse_brks:
    This recovers the loci from the hypothesized breakpoints identified in earlier
    steps.
    """    
    print('Starting breakpoint testing...')
    parse_brks(break_tab)

    # hypothesis_dict stores anchors and othersides by nt
    hypothesis_dict = {}
    # each anchor has a single contig, these contigs are collected by nt
    anchor_contig_dict = {}
    
    cluster = 0
    qname_ct = 0
    for hypo, qname_list in complete_qname_dict.items():
        qname_ct +=1
        if qname_ct % 1 == 0:
            outline = ('Breakpoint Hypotheses Tested: {} \t Number remaining {} \t Percent Complete {}\n').format(qname_ct, len(complete_qname_dict)-qname_ct, round(qname_ct/float(len(complete_qname_dict)),2))
            print(outline)
        
        temp_seq_file_name = ('{}/{}_temp_seq.fa').format(temp_dir, hypo)
        make_fasta(temp_seq_file_name, hypo, qname_list)
        
        msa_dict, hash_to_name = parse_msa_file(temp_seq_file_name)
        
        distance_dict = {}
        # distance dict represents distance (mismatch) between two sequences
        for name1, seq1 in msa_dict.items():
            for name2, seq2 in msa_dict.items():
                check_compare = ('{}_{}').format(min(name1, name2), max(name1, name2))
                
                if (name1 != name2) and (check_compare not in distance_dict):
                    if name1 == min(name1, name2):
                        seqA = seq1
                        seqB = seq2
                    else:
                        seqA = seq2
                        seqB = seq1
                        
                    distance_dict[check_compare] = parse_msa_kmer(seqA, seqB)
                    
        prefix_list, contig_seq_dict = build_clusters(distance_dict, msa_dict, hash_to_name, contig_seq_dict)  

        for prefix in prefix_list:
            json_short_file_name = ('{}_short_blastn.json').format(prefix)
            json_long_file_name = ('{}_long_blastn.json').format(prefix)
            cluster_seq_name = ('{}_cluster_seq.fa').format(prefix)
            
            seq = get_seq(cluster_seq_name)
            
            bashCommand = ('blastn -task "blastn-short" -query {} -subject {fa_file} -outfmt 15 -out {}').format(cluster_seq_name, json_short_file_name, fa_file=fa_file)
            subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)
            hypothesis_dict, anchor_contig_dict = parse_json(json_short_file_name, seq, hypothesis_dict, anchor_contig_dict)
            #gff_list_dict, gff_rank_dict, blast_read_aligned_sections_dict, blast_genome_aligned_regions_dict, gff_uid_dict = parse_json(json_short_file_name, contig_seq_dict, gff_list_dict, gff_rank_dict, blast_read_aligned_sections_dict, blast_genome_aligned_regions_dict, gff_uid_dict)
       
            bashCommand = ('blastn -query {} -subject {fa_file} -outfmt 15 -out {}').format(cluster_seq_name, json_long_file_name, fa_file=fa_file)
            subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)
            hypothesis_dict, anchor_contig_dict = parse_json(json_long_file_name, seq, hypothesis_dict, anchor_contig_dict)
            #gff_list_dict, gff_rank_dict, blast_read_aligned_sections_dict, blast_genome_aligned_regions_dict, gff_uid_dict = parse_json(json_long_file_name, contig_seq_dict, gff_list_dict, gff_rank_dict, blast_read_aligned_sections_dict, blast_genome_aligned_regions_dict, gff_uid_dict)

    gff_set, gff_header, vcf_set, vcf_header = summarize_hypotheses(hypothesis_dict, anchor_contig_dict, gap, resource_dict)
    

    resource_pickle_name = ('{}/hypothesis_dict.p').format(pickles_dir)
    with open(resource_pickle_name, 'wb') as file:
         pickle.dump(hypothesis_dict, file)
            
    gff_file_name = ('{}/{}_SV_CNV.gff').format(final_output_dir, output_file)
    gff_file = open(gff_file_name,'w')    
    
    gff_file.write(gff_header)
    
    for each_gff in gff_set:
        gff_file.write(each_gff)
    
    gff_file.close()
    
    vcf_file_name = ('{}/{}_SV_CNV.vcf').format(final_output_dir, output_file)
    vcf_file = open(vcf_file_name,'w')    
    
    vcf_file.write(vcf_header)
    
    for each_vcf in vcf_set:
        vcf_file.write(each_vcf)
    
    vcf_file.close()
    
    bashCommand = ('bcftools sort {vcf_file} > {vcf_file}.sorted').format(
        vcf_file = vcf_file_name)
    print(bashCommand)       
    subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)
    
    bashCommand = ('bgzip {vcf_file}.sorted -cf > {vcf_file}.gz').format(
        vcf_file = vcf_file_name)
    print(bashCommand)       
    subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)
    
    bashCommand = ('tabix -p vcf {vcf_file}.gz').format(
        vcf_file = vcf_file_name)
    print(bashCommand)       
    subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)
    
    outline = ('\t{} nodes realigned to {} regions.\n\n\t'
               'For a complete representation refer to:\n\t\t{} or {}').format(
                   len(complete_qname_dict), len(gff_set), gff_file_name, vcf_file_name)
    print(outline)
    
if args.separate_by_strand:
    infile_name = args.input_file
    outfile = args.output_file
    
    if infile_name[-4:] == '.bam':
        infile_name = infile_name.split('.bam')[0]
        bashCommand = ('samtools view -h {infile}.bam -o {infile}.sam').format(infile=infile_name)
        print(bashCommand)       
        subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)
        
        infile_name = ('{infile}.sam').format(infile=infile_name)
        
    if outfile[-4:] == '.sam':
        outfile = outfile.split('.sam'[0])
        
    infile = open(infile_name)
        
    plus_file_name = ('{outfile}_plus').format(outfile=outfile)
    plus_file = open(plus_file_name+'.sam', 'w')
    
    minus_file_name = ('{outfile}_minus').format(outfile=outfile)
    minus_file = open(minus_file_name+'.sam', 'w')
    
    for line in infile:
        if line[0] == '@':
            plus_file.write(line)
            minus_file.write(line)
            
        if line[0] != '@':            
            strand = unpackbits(np.array([int(line.split('\t')[1])]))[0][4]
            
            if strand == 0:
                plus_file.write(line)
                
            if strand == 1:
                minus_file.write(line)
    
    infile.close()        
    plus_file.close()
    minus_file.close()
    
    convert_sort(plus_file_name)
    convert_sort(minus_file_name)
    
if args.find_breakpoints:
    #output_dir, output_file
    run_name = args.run_name
    
    bash_file_name = output_file.rsplit('.', 1)[0] + '.sh'
    
    bash_file = open(bash_file_name, 'w')
            
    outline = ('python erisapfel.py -make -fa {ref_fa}'
               ' -fastq_1 {fastq_1} -fastq_2 {fastq_2}'
               ' -run_name {run_name}\n').format(
                   ref_fa = args.fa_file, 
                   fastq_1 = args.fastq_1,
                   fastq_2 = args.fastq_2,
                   run_name = run_name)
                  
    bash_file.write(outline)
    
    outline = ('python erisapfel.py -depth'
               ' -run_name {run_name}\n').format(
                   run_name = run_name)
    bash_file.write(outline)
    
    outline = ('python erisapfel.py -peaks'
               ' -run_name {run_name}\n').format(
                   run_name = run_name)
    bash_file.write(outline)
                   
    
    if args.make_filter:
        filter_bed = args.filter_bed
        
        outline = ('python erisapfel.py -filter -filter_bed {filter_bed} -o {run_name}_filter.p\n').format(
            filter_bed = filter_bed, 
            run_name = run_name)
        
        bash_file.write(outline)
        
        outline = ('python erisapfel.py -map -run_name {run_name} --filter_object {run_name}_filter.p\n').format(
            run_name = run_name)
        
        bash_file.write(outline)
    
    outline = ('python erisapfel.py -map'
               ' -run_name {run_name}\n').format(
                   run_name = run_name)
    bash_file.write(outline)
    
    
    outline = ('python erisapfel.py -localseq'
               ' -run_name {run_name}\n').format(
                   run_name = run_name)
    bash_file.write(outline)
                                    
    monolog = ('\tRunning breakpoint identification using {bfn} ...\n').format(bfn = bash_file_name)
    print(monolog)    
    bashCommand = ('bash {bfn}').format(bfn = bash_file_name) 
    print(bashCommand)       
    subprocess.run([bashCommand], stderr=subprocess.STDOUT, shell=True)
    
if args.model_selection or args.model_train or args.model_predict:
    #TODO future version    
    from sklearn.preprocessing import MinMaxScaler
    from sklearn.preprocessing import StandardScaler
    from sklearn.preprocessing import Normalizer
    from sklearn.preprocessing import PowerTransformer
    from sklearn.preprocessing import QuantileTransformer
    
    from sklearn.model_selection import train_test_split

    from sklearn.metrics import confusion_matrix
    from sklearn.metrics import f1_score
    from imblearn.metrics import geometric_mean_score
    
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.naive_bayes import GaussianNB
    from sklearn.naive_bayes import ComplementNB
    from sklearn.linear_model import LogisticRegression
    from sklearn.ensemble import RandomForestClassifier
    
    from imblearn import FunctionSampler
    from imblearn.under_sampling import RandomUnderSampler, NearMiss, AllKNN
    from imblearn.over_sampling import ADASYN, RandomOverSampler, SMOTE
    from imblearn.pipeline import make_pipeline

if args.model_selection:    
    mge_file_name = (args.input_file)
    
    value_df_object = mge_to_model_parser(mge_file_name)
                    
    model_selection(value_df_object)
    
if args.model_train:    
    mge_file_name = (args.input_file)
    
    value_df_object = mge_to_model_parser(mge_file_name)
                    
    topmodel = model_train(value_df_object) 
    
if args.model_predict:
    resource_dict = io_load()        
    sample = resource_dict['run_name']
    
    feature_pickle = ('{}/{}_feature_values.p').format(final_output_dir, sample)
    value_df_object = {}
    value_df_object[sample] = feature_pickle
                                
    topmodel = model_predict(value_df_object) 
    

                            
