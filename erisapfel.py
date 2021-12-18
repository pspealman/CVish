# -*- coding: utf-8 -*-
"""
erisapfel - "apple of eris, goddess of discord"
# Purpose: Uses split, discordant, and soft reads to identify genomic breakpoints
# further processing is used to extract context and sequence.

#change log
(Initial Release - Public Beta: Rider Theory)

# 12.13.2021
Version 1.0 (Write Reflection)
    Change log
    _x_ 'gff' error
    _x_ updated requirement list
    _x_ corrected RPKM

Version 1.0.1 (...)
    _x_ mkdir error


@author: Pieter Spealman 
"""

"""Requirements :
bwa/intel   0.7.17
samtools    1.14
bedtools    2.29.2
velvet      1.2.10
blast+      2.11.0
samblaster  0.1.26
mafft       7.475
emboss      6.6.0
"""
#reduced_file
import os
import argparse
import subprocess
import pickle
import pandas as pd
import numpy as np
import json
import re   
import scipy.stats as stats

import warnings
warnings.simplefilter("ignore")

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
    'demo: python erisapfel.py -breakpoint -load_filter demo/output/filter.p -i demo/output/demo_evo -break_tab demo/output/demo_evolved.tab -break_bed demo/output/demo_evolved.bed\n\n'+
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

#chromo_list = ['NC_001133.9','NC_001134.8','NC_001135.5','NC_001136.10','NC_001137.3','NC_001138.5','NC_001139.9','NC_001140.6', 
#'NC_001141.2','NC_001142.9','NC_001143.9','NC_001144.5','NC_001145.3','NC_001146.8','NC_001147.6','NC_001148.4']

def output_handler(output):
    if len(output.strip()) > 1:
        print(output)

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
    
def test():
    monolog = ('=== Currently Testing Erisapfel.py ===')
    print(monolog)
    
    monolog = ('\tTesting Step 0. make bwa for demo.fna\n')
    print(monolog)    
    bashCommand = ('bwa index demo/input/demo.fna')
    print(bashCommand)       
    output_handler(subprocess.check_output([bashCommand],stderr=subprocess.STDOUT,shell=True))
    

    monolog = ('\tTesting Step 1. -make command for ancestor\n')
    print(monolog)    
    bashCommand = ('python erisapfel.py -make -fa demo/input/demo.fna -gff demo/input/demo.gff -fastq_1 demo/input/n01_ancestor.fastq.gz -fastq_2 demo/input/n02_ancestor.fastq.gz -o demo/output/demo_anc')
    print(bashCommand)       
    output_handler(subprocess.check_output([bashCommand],stderr=subprocess.STDOUT,shell=True))
    
    monolog = ('\tTesting Step 2. -make command for evolved\n')
    print(monolog)    
    bashCommand = ('python erisapfel.py -make -fa demo/input/demo.fna -gff demo/input/demo.gff -fastq_1 demo/input/n01_evolved.fastq.gz -fastq_2 demo/input/n02_evolved.fastq.gz -o demo/output/demo_evo')
    print(bashCommand)       
    output_handler(subprocess.check_output([bashCommand],stderr=subprocess.STDOUT,shell=True))
    
    monolog = ('\tTesting Step 3. -breakpoint command for ancestor\n')
    print(monolog)    
    bashCommand = ('python erisapfel.py -breakpoint -i demo/output/demo_anc -break_tab demo/output/demo_ancestor.tab -break_bed demo/output/demo_ancestor.bed')
    print(bashCommand)       
    output_handler(subprocess.check_output([bashCommand],stderr=subprocess.STDOUT,shell=True))

    monolog = ('\tTesting Step 4. -filter command using ancestor\n')
    print(monolog)    
    bashCommand = ('python erisapfel.py -filter -gff demo/input/demo_filter.gff -ancestor demo/output/demo_ancestor.bed -o demo/output/filter.p')
    print(bashCommand)       
    output_handler(subprocess.check_output([bashCommand],stderr=subprocess.STDOUT,shell=True))
    
    monolog = ('\tTesting Step 5. -breakpoint command for evolved\n')
    print(monolog)    
    bashCommand = ('python erisapfel.py -breakpoint -load_filter demo/output/filter.p -i demo/output/demo_evo -break_tab demo/output/demo_evolved.tab -break_bed demo/output/demo_evolved.bed')
    print(bashCommand)       
    output_handler(subprocess.check_output([bashCommand],stderr=subprocess.STDOUT,shell=True))
    
    monolog = ('\tTesting Step 6. -localseq command\n')
    print(monolog)    
    bashCommand = ('python erisapfel.py -localseq -split_only -break_tab demo/output/demo_evolved.tab -break_bed demo/output/demo_evolved.bed -fastq_1 demo/input/n01_evolved.fastq.gz -fastq_2 demo/input/n02_evolved.fastq.gz -o demo/input/localseq_evo')
    print(bashCommand)       
    output_handler(subprocess.check_output([bashCommand],stderr=subprocess.STDOUT,shell=True))
    
    monolog = ('\tTesting Step 7. -realign command\n')
    print(monolog)    
    bashCommand = ('python erisapfel.py -realign -fa demo/input/demo.fna -contigs demo/input/localseq_evo/ -o demo/input/localseq_evo/realign_evo')
    print(bashCommand)       
    output_handler(subprocess.check_output([bashCommand],stderr=subprocess.STDOUT,shell=True))
    
parser = argparse.ArgumentParser()
#handle help_dialog
parser.add_argument('-man',"--manual", action='store_true')
parser.add_argument('-demo',"--demo",action='store_true')
parser.add_argument('-test',"--test",action='store_true')
parser.add_argument('-run_name', '--run_name')

parser.add_argument('-make',"--make", action='store_true')
parser.add_argument('-ends',"--paired_end")
parser.add_argument('-fa',"--fa_file")
parser.add_argument('-bam', "--bam_file")

parser.add_argument('-filter', '--make_filter', action='store_true')
parser.add_argument('-filter_gff','--filter_gff', nargs='+')
parser.add_argument('-filter_bed', '--filter_bed', nargs='+')
parser.add_argument('-filter_tab', '--filter_tab', nargs='+')

parser.add_argument('-map', '--map_reads', action='store_true')
parser.add_argument('-contig',"--contig", action='store_true')
parser.add_argument('-z_split',"--split_zscore")
parser.add_argument('-z_disco',"--disco_zscore")

parser.add_argument('-qval', '--mapq_val')

parser.add_argument('-load_filter', "--load_filter")
parser.add_argument('-flank','--filter_flanking')

parser.add_argument('-break_tab', "--break_tab")
parser.add_argument('-break_bed', "--break_bed")
parser.add_argument('-split_only', '--split_only', action='store_true')

parser.add_argument('-chromo_size', "--chromo_size")

parser.add_argument('-localseq', '--build_sequence', action='store_true')
parser.add_argument('-soft', "--soft_file")
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

parser.add_argument('-depth', '--depth_analysis', action='store_true')

parser.add_argument('-read_type','--read_type_list', nargs='+')

parser.add_argument('-peaks', '--peaks', action='store_true')
parser.add_argument('-min_length', '--CNV_min_length')
parser.add_argument('-gap', '--max_gap')
parser.add_argument('-purity', '--min_purity')
parser.add_argument('-windowy', '--window')

parser.add_argument('-breaks','--breakpoints', action='store_true')

args = parser.parse_args()

if args.manual:
    help_dialog()
    
if args.demo:
    demo()
    
if args.test:
    test()

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

if args.load_filter:
    region_filter_dict = pickle.load(open(args.load_filter, 'rb'))

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
    
#Define default directories
s_path = ('{}{}').format(output_dir, output_file)
bam_dir = ('{}/bam/').format(s_path)
idx_dir = ('{}/idx/').format(s_path)
# TODO rename cluster 'temp'
temp_dir = ('{}/cluster/').format(s_path)
pickles_dir = ('{}/pickles/').format(s_path)
    
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

def pickle_loader(file_name, runmode):
    if runmode == 'dict':
        file = open(file_name,'rb')
        object_file = pickle.load(file)
        file.close()
    
    if runmode == 'df':
        object_file = pd.read_pickle(file_name)
    return(object_file)

def io_make(runmode='resource'):
    if runmode == 'resource':
        resource_file_name = ('{}/resource.tab').format(idx_dir)
        resource_file = open(resource_file_name, 'w')
        
        resource_file.close()
        
    if runmode == 'ml':
        ml_file_name = ('{}/ml_hypo.tab').format(idx_dir)
        ml_file = open(ml_file_name, 'w')

        for cat, val in ml_hypo_dict.items():                
            outline = ('{}\t{}\n').format(cat,val)
            ml_file.write(outline)
        ml_file.close()
        
        ml_file_name = ('{}/ml_loci.tab').format(idx_dir)
        ml_file = open(ml_file_name, 'w')

        for cat, val in ml_loci_dict.items():                
            outline = ('{}\t{}\n').format(cat,val)
            ml_file.write(outline)
        ml_file.close()
        
def io_append(resource_dict, runmode='resource'):
    if runmode == 'resource':
        pre_existing_dict = {}
        
        resource_file_name = ('{}/resource.tab').format(idx_dir)    
        resource_file = open(resource_file_name)
        
        for line in resource_file:
            cat, val = line.split('\t')
            pre_existing_dict[cat]=val
        resource_file.close()
    
        resource_file_name = ('{}/resource.tab').format(idx_dir)    
        resource_file = open(resource_file_name, 'a')
        
        for cat, val in resource_dict.items():
            if cat not in pre_existing_dict:
                outline = ('{}\t{}\n').format(cat,val)
                resource_file.write(outline)
                
        resource_file.close()
        
        resource_pickle_name = ('{}/resource.p').format(idx_dir)         
        with open(resource_pickle_name, 'wb') as file:
            pickle.dump(resource_dict, file)

    if runmode == 'ml_hypo':
        ml_pickle_name = ('{}/ml_hypo.p').format(idx_dir)         
        with open(ml_pickle_name, 'wb') as file:
            pickle.dump(ml_hypo_dict, file)
            
        ml_pickle_name = ('{}/ml_loci.p').format(idx_dir)         
        with open(ml_pickle_name, 'wb') as file:
            pickle.dump(ml_loci_dict, file)
        
def io_load(runmode='resource'):
    if runmode == 'resource':
        resource_file_name = ('{}/resource.p').format(idx_dir)
        
        with open(resource_file_name, 'rb') as resource_file:
            resource_dict = pickle.load(resource_file)
            
        return(resource_dict) 
        
    # TODO change idx to pickle
    if runmode == 'ml':
        ml_pickle_name = ('{}/ml_hypo.p').format(idx_dir) 
        with open(ml_pickle_name, 'rb') as ml_file:
            ml_hypo_dict = pickle.load(ml_file)

        ml_pickle_name = ('{}/ml_loci.p').format(idx_dir) 
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
        qname = qname.split('_')[1]
        #print('pf_q', qname)
        if qname in uniuid_to_seq_dict:
            seq = uniuid_to_seq_dict[qname]
            outline = ('>{}\n{}\n').format(qname,seq)
            temp_split_seq_file.write(outline) 
            #print(outline)
        else:
            1/0
            
    temp_split_seq_file.close()
    
def degzip(full_name, file_name):
    if (full_name == file_name + '.gz') or (full_name == file_name):
        bashCommand = ('gunzip -c {} > {}').format(full_name, file_name)
        output_handler(subprocess.check_output([bashCommand],stderr=subprocess.STDOUT,shell=True))
        print(bashCommand)
        
    else:
        bashCommand = ('cp {} {}.gz').format(full_name, file_name)
        output_handler(subprocess.check_output([bashCommand],stderr=subprocess.STDOUT,shell=True))
        print(bashCommand)
        
        bashCommand = ('gunzip -f {}').format(file_name)
        output_handler(subprocess.check_output([bashCommand],stderr=subprocess.STDOUT,shell=True))
        print(bashCommand)
        
def parse_fastq(fastq_full_name, series_number, qname_lookup, uniuid_to_seq_dict, uniuid_to_phred_dict):
    
    outline = ('Parsing {}, using {}').format(fastq_full_name, series_number)
    print(outline)
    
    if '/' in fastq_full_name:
        nope_dir, fastq_file_name = handle_outfile(fastq_full_name)
        
        if fastq_full_name[-3:]=='.gz':
            #
            degzip(fastq_full_name, fastq_file_name)
            fastq_file_name = fastq_file_name[:-3]
        else:
            bashCommand = ('cp -f {} {}').format(fastq_full_name, fastq_file_name)
            output_handler(subprocess.check_output([bashCommand],stderr=subprocess.STDOUT,shell=True))

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
                    #print('seq', seq)
                
                if uid in uniuid_to_phred_dict:
                    phred = uniuid_to_phred_dict[uid]
                    #print('phred', phred)
                
                outline = ('@{}\n{}\n+\n{}\n').format(uid, seq, phred)
                outfile.write(outline)
                #print(outline)
            
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

def parse_brks(break_file_name, each_type):
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
    
def unpackbits(x,runmode='upb'):
    #TODO replace the unpack bits functions
    #    
    #0  1    0x1   (rp)    read_paired
    #1  2    0x2   (rmp)   read_mapped_in_proper_pair
    #2  4    0x4   (ru)    read_unmapped
    #3  8    0x8   (mu)    mate_unmapped
    #4  16   0x10  (rrs)   read_reverse_strand
    #5  32   0x20  (mrs)   mate_reverse_strand
    #6  64   0x40  (fip)   first_in_pair
    #7  128  0x80  (sip)   second_in_pair
    #8  256  0x100 (npa)   not_primary_alignment
    #9  512  0x200 (rfp)   read_fails_platform
    #10 1024 0x400 (pcr)   read_is_PCR_or_optical_duplicate
    #11 2058 0x800 (sa)    supplementary_alignment
     
    upb = []
    
    for exp in range(11,-1,-1):
        bit = 2**exp
        
        if (x - bit) >= 0:
            upb.append(1)
            x-=(bit)
        else:
            upb.append(0)

    if runmode == 'upb':
        return(upb[:8])
    
    # psb -> pseudo-samblaster
    if runmode == 'psb':
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
    
        psb = upb[::-1]
        return(psb)
        # if (psb[1] == 0) and (psb[2] == 0) and (psb[3] == 0):
        #     return('evaluate')
        
        # if 
              

""" Step One """
if args.make:
    """ This command uses bwa mem to align the fastq files to the reference genome
    These are passed to samblaster which seperates out discordant and split reads
    and then samtools for conversion and index generation.
    """
    
    chromo_size_dict = {}
    resource_dict = {}
    
    #Create default directories
    
    def make_subdirectory(is_path,subdir):
        new_dir = ('{}/{}').format(is_path, subdir)
        if not os.path.exists(new_dir):
            os.makedirs(new_dir)
    
    make_subdirectory(s_path,'bam')
    make_subdirectory(s_path,'idx')
    make_subdirectory(s_path,'cluster')
    make_subdirectory(s_path,'pickles')
    io_make()
    
    resource_dict['run_name']=args.run_name
    resource_dict['results_dir']=args.run_name
    resource_dict['genome_fa']=args.fa_file
    #resource_dict['genome_gff']=args.gff_file
           
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
    command_file.write('#!/bin/bash\n')
    
    ref_fa=('\tgenome_fa={}\n').format(args.fa_file)
    #ref_gff=('\tgenome_gff={}\n').format(args.gff_file)
         
    if (end_is == 'pe') or (end_is == 'paired_end'):
        infastq_1=('fastq_1={}\n').format(args.fastq_1)
        infastq_2=('fastq_2={}\n').format(args.fastq_2)
    else:
        infastq_1=('fastq_1={}\n').format(args.fastq_1)
        infastq_2=('fastq_2={}\n').format("")
        
    output_filename=('output_file={bam_dir}/{output_file}\n'
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
    output = subprocess.check_output([bash_command],stderr=subprocess.STDOUT,shell=True)
    print(output)
    
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
    resource_dict['soft_file'] = (bam_dir+output_file+str('_soft.bam'))
    
    resource_dict['sam_file'] = (bam_dir+output_file+str('.sam'))
    resource_dict['discordant_sam_file'] = (bam_dir+output_file+str('_discordant.sam'))
    resource_dict['split_sam_file'] = (bam_dir+output_file+str('_split.sam'))
    resource_dict['soft_sam_file']=(bam_dir+output_file+str('_soft.sam'))
    
    resource_dict['chromo_size']=(pickle_name)

    io_append(resource_dict)
    
""" Step Two"""
if args.depth_analysis:    
    """
    02.22.19 - depth_analysis
        Function: pulls read depths from bam files into chromosomes.
    """    
    resource_dict = io_load()
    print(resource_dict)
    if 'genome_fa' not in resource_dict:
        1/0
        
    fa_file = resource_dict['genome_fa']
    each_sample = resource_dict['run_name']
    bam_file = resource_dict['bam_file']
    sam_file = resource_dict['sam_file']
    soft_file = resource_dict['soft_file']
    soft_sam_file = resource_dict['soft_sam_file']
    

    #handle read_type
    if not args.read_type_list:
        read_type_list = ['RD','discordant','split','soft']
    else:
        read_type_list = args.read_type_list
        
    resource_dict['read_types']=read_type_list

    #handle run_name    
    #if args.run_name:
    #    path = bam_dir
    #    each_sample = output_file
           
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
    
#    total_reads_dict = {}
#    read_map_dict = {}
    
    def run_soft(bam_file, soft_sam_file, soft_file):
        #handling soft_reads
        #bashCommand = ('samtools view -h {} > {}/temp_global.sam').format(bam_file, bam_dir)
        #print(bashCommand)
        #output_handler(subprocess.check_output([bashCommand],stderr=subprocess.STDOUT,shell=True))

        #softfile_name = ('{}{}.sam').format(bam_dir, each_sample)
        is_softfile = open(soft_sam_file,'w')
        
        #sam_file_name = ('{}/temp_global.sam').format(bam_dir)
        infile = open(sam_file)
        for line in infile:
            if line[0] == '@':
                is_softfile.write(line)
            if line[0] != '@':
                cigar = line.split('\t')[5]
                if '*' not in cigar: 
                    size = parse_cigar(cigar, 'non-matching')
                    if size >= 10:
                        is_softfile.write(line)
        infile.close()
        is_softfile.close()
                    
        bashCommand = ('samtools view -Sb {} -o {}/temp_unsorted.bam').format(soft_sam_file, bam_dir)            
        print(bashCommand)
        output_handler(subprocess.check_output([bashCommand],stderr=subprocess.STDOUT,shell=True))
        bashCommand = ('samtools sort {}/temp_unsorted.bam -o {}').format(bam_dir, soft_file)
        print(bashCommand)
        output_handler(subprocess.check_output([bashCommand],stderr=subprocess.STDOUT,shell=True))
        bashCommand = ('samtools index {}').format(soft_file)
        print(bashCommand)
        output_handler(subprocess.check_output([bashCommand],stderr=subprocess.STDOUT,shell=True))

    def run_bedgraph(bam_file, each_sample):
        print('\tRunning bedtools genomecov on '+ str(infile_name) + '...')
        bashCommand = ('bedtools genomecov -bg -ibam {} > {}{}.bedgraph').format(bam_file, idx_dir, each_sample)
        print(bashCommand)
        output_handler(subprocess.check_output([bashCommand],stderr=subprocess.STDOUT,shell=True))

    def run_mpileup(infile_name, each_sample, runmode):            
        if runmode == 'RD':
            print('\tRunning samtools to determine sample depth...')
            if filter_by_bed:
                bashCommand = ('samtools mpileup -B -C 50 -f {} {} -l {} -a -o {}/{}_{}.mpileup').format(fa_file, infile_name, filter_bed, pickles_dir, each_sample, runmode)        
            else:
                bashCommand = ('samtools mpileup -B -C 50 -f {} {} -a -o {}/{}_{}.mpileup').format(fa_file, infile_name, pickles_dir, each_sample, runmode)                        
            print(bashCommand)
            output_handler(subprocess.check_output([bashCommand],stderr=subprocess.STDOUT,shell=True))

        if runmode == 'discordant':
            print('\tRunning samtools to determine discordant read depth...')
            if filter_by_bed:
                bashCommand = ('samtools mpileup -B -C 25 -f {} {} -l {} -a -A -o {}/{}_{}.mpileup').format(fa_file, infile_name, filter_bed, pickles_dir, each_sample, runmode)            
            else:
                bashCommand = ('samtools mpileup -B -C 25 -f {} {} -a -A -o {}/{}_{}.mpileup').format(fa_file, infile_name, pickles_dir, each_sample, runmode)
            print(bashCommand)
            output_handler(subprocess.check_output([bashCommand],stderr=subprocess.STDOUT,shell=True))
            
        if runmode == 'split':
            print('\tRunning samtools to determine split read depth...')
            if filter_by_bed:
                bashCommand = ('samtools mpileup -B -C 10 -f {} {} -l {} -a -A -o {}/{}_{}.mpileup').format(fa_file, infile_name, filter_bed, pickles_dir, each_sample, runmode)     
            else:
                bashCommand = ('samtools mpileup -B -C 10 -f {} {} -a -A -o {}/{}_{}.mpileup').format(fa_file, infile_name, pickles_dir, each_sample, runmode)     
            print(bashCommand)            
            output_handler(subprocess.check_output([bashCommand],stderr=subprocess.STDOUT,shell=True))
            
        if runmode == 'soft':
            print('\tRunning samtools to determine soft clip read depth...')
            if filter_by_bed:
                bashCommand = ('samtools mpileup -B -C 25 -f {} {} -l {} -a -A -o {}/{}_{}.mpileup').format(fa_file, infile_name, filter_bed, pickles_dir, each_sample, runmode)            
            else:
                bashCommand = ('samtools mpileup -B -C 25 -f {} {} -a -A -o {}/{}_{}.mpileup').format(fa_file, infile_name, pickles_dir, each_sample, runmode)
            print(bashCommand)
            output_handler(subprocess.check_output([bashCommand],stderr=subprocess.STDOUT,shell=True))
            
    def process_bam(infile_name, read_type):        
        print('\tRunning samtools to determine total mapped reads...')
        bashCommand = ('samtools view -F 0x4 {} | cut -f 1 | sort | uniq | wc -l > temp.ct').format(infile_name)
        print(bashCommand) 
        output_handler(subprocess.check_output([bashCommand],stderr=subprocess.STDOUT,shell=True))
        
        temp_ct_file = open('temp.ct')
        for line in temp_ct_file:
            line = line.strip()
            total_ct = int(line)
        temp_ct_file.close()
        
        if read_type == 'RD':
            resource_dict['total_RD_aligned_reads']=total_ct
        if read_type == 'split':
            resource_dict['total_split_aligned_reads']=total_ct
        if read_type == 'discordant':
            resource_dict['total_discordant_aligned_reads']=total_ct
        if read_type == 'soft':
            resource_dict['total_soft_aligned_reads']=total_ct
                
        return(total_ct)
                                
    for each_type in read_type_list:
        if each_type == 'RD':
            infile_name = resource_dict['bam_file']
        if each_type == 'split':
            infile_name = resource_dict['split_file']
        if each_type == 'discordant':
            infile_name = resource_dict['discordant_file']
        if each_type == 'soft':
            infile_name = resource_dict['soft_file']
        
        outline = ('\nProcessing sample {}...').format(each_sample)
        print(outline)
        
        if each_type == 'soft':
            run_soft(bam_file, soft_sam_file, soft_file)
        
        total_reads = process_bam(infile_name, each_type)
        
        run_bedgraph(bam_file, each_sample)
        
        print('\tParsing mpileup into chromosomes...')

        populated=False
        mpileup_df = []
        
        run_mpileup(infile_name, each_sample, each_type)

        raw_mpileup = str(pickles_dir + '/' + each_sample + '_' + each_type + '.mpileup')
        temp_df = pd.read_csv(raw_mpileup, sep='\t', header=None, names = ['chromo', 'nuc', 'base', 'ct', 'val1', 'val2'])

        outline = ('\tProcessing {}...').format(each_type)
        print(outline)
        
        if total_reads > 0:
            rpkm = lambda x: float(150*int(x)*10**9)/float(total_reads)            
        else:
            1/0
            rpkm = lambda x: float(x)
        
        if populated:
            mpileup_df['rpkm'] = temp_df['ct'].map(rpkm)
            mpileup_df['chromo']=temp_df['chromo']
            
        else:
            populated = True
            mpileup_df = temp_df[['ct']].copy()
            mpileup_df['nuc']=temp_df['nuc']
            mpileup_df['rpkm'] = temp_df['ct'].map(rpkm)
            mpileup_df['chromo']=temp_df['chromo']
                       
        pickle_name = ('{}/mpileup_{}_{}.p').format(pickles_dir, each_sample, each_type)
        mpileup_df.to_pickle(pickle_name)
        print(mpileup_df.head())

        fzero = mpileup_df.replace(0, np.NaN)
        fg_median = fzero["ct"].median()
        fg_std = fzero["ct"].std()
        
        zscore_read_type = ('{}_zscore_median').format(each_type)
        resource_dict[zscore_read_type] = fg_median
        zscore_read_type = ('{}_zscore_std').format(each_type)
        resource_dict[zscore_read_type] = fg_std
            
    io_append(resource_dict)
        
""" Step Three """
if args.peaks:
    """
    02.28.19 - peak detector
        additional feature, filter regions of chromosome with large deletions
        
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
        
    #handle run_name    
#    if args.run_name:
#        path = bam_dir
#        each_sample = output_file
                
    file_name_lookup = []
    depth_dict = {}
    chromo_depth_dict = {}
    
    def long_walk(seed, long_range, rpkm_dict):
        nt = seed
        miss_ct = 0
        depth = []

        for s_step in long_range:
            if s_step in rpkm_dict:
                depth.append(float(rpkm_dict[s_step]))
                nt=s_step
                miss_ct = 0
                
            else:
                miss_ct += 1

            if miss_ct >= 3*gap:
                return(nt,depth)
                
        return(nt,depth)
    
    def local_walk(short_range, rpkm_dict):
        hit_ct = 0

        for s_step in short_range:
            if s_step in rpkm_dict:
                hit_ct += 1
                
            if hit_ct >= (purity*min_length):
                return(True)
        else:
            return(False)
            
    def build_trace(subz_df, suby_df, each_type, c_median, g_median):
        print('running build_trace')
        
        subz_dict = subz_df.to_dict(orient='index')
        subz_rpkm_dict = {}
        for index, deets in subz_dict.items():
            nuc = int(deets['nuc'])
            rpkm = float(deets['ct'])
            subz_rpkm_dict[nuc]=rpkm
            
        suby_dict = suby_df.to_dict(orient='index')
        suby_rpkm_dict = {}
        for index, deets in suby_dict.items():
            nuc = int(deets['nuc'])
            rpkm = float(deets['ct'])
            suby_rpkm_dict[nuc]=rpkm
        
        seed_set = set()
        
        if len(subz_rpkm_dict) > 1:
            min_nuc = min(subz_rpkm_dict)
            max_nuc = max(subz_rpkm_dict)                            
    
            for nuc, rpkm in subz_rpkm_dict.items():            
                if (nuc not in seed_set):
                    short_range = range(nuc-int(min_length/2),nuc+int(min_length/2))
                    #print(short_range)
                    
                    if(local_walk(short_range, subz_rpkm_dict)) or (each_type != 'RD'):
                        left_long_range = range(min_nuc-int(min_length/2),nuc)[::-1]
                        right_long_range = range(nuc,max_nuc+int(min_length/2))
                        
                        start, depth = long_walk(nuc, left_long_range, suby_rpkm_dict)          
                        stop, depth_2 = long_walk(nuc, right_long_range, suby_rpkm_dict)
            
                        depth += depth_2
                        
                        depth_median = np.median(depth)
                        depth_std = np.std(depth)
    
                        for nt in range(start, stop+1):
                            seed_set.add(nt)
                            
                        # if (c_median <= 0):
                        #     print('c_median', c_median)
                        #     1/0
                            
                        if (start >= 0) and (stop >= 0):
                            rel_c = depth_median/float(c_median)
                            rel_g = depth_median/float(g_median)
                            print('rel_c', depth_median, c_median)
                            
                            outline = ('{chromo}\t{each_type}_depth\tCNV'
                                       '\t{start}\t{stop}\t{rel_c}\t.\t.'
                                       '\tID={chromo}:{start}-{stop};rel_chromosome_RD={rel_c};rel_genome_RD={rel_g};sample={sample}\n'
                                       ).format(
                                           chromo=every_chromo, 
                                           each_type=each_type, 
                                           start=start, stop=stop, 
                                           rel_c=round(rel_c,2), 
                                           rel_g = round(rel_g,2), 
                                           sample = each_sample)
                            depth_gff_outfile.write(outline)
                            print(outline)
                            
                            locus = ('{}:{}-{}').format(every_chromo, start, stop)
                            
                            depth_dict[locus] = [rel_c, rel_g, depth_std]

    for each_type in read_type_list:
        print(each_type)
        depth_gff_file_name = ('{}/{}_{}_depth.gff').format(idx_dir, each_sample, each_type)
        depth_gff_outfile = open(depth_gff_file_name,'w')
        
        chromo_gff_file_name = ('{}/{}_{}_chromosome.gff').format(idx_dir, each_sample, each_type)
        chromo_gff_outfile = open(chromo_gff_file_name,'w')
    
        pickle_name = ('{}/mpileup_{}_{}.p').format(pickles_dir, each_sample, each_type)
        mpileup_df = pickle_loader(pickle_name, 'df')
        #print(mpileup_df.head())
        #1/0
                
        zscore_read_type = ('{}_zscore_median').format(each_type)
        g_median = resource_dict[zscore_read_type]
        zscore_read_type = ('{}_zscore_std').format(each_type)
        g_std = resource_dict[zscore_read_type]
                        
        outline = ('Global Median:\t{}\nGlobal Std:\t{}\n').format(g_median,g_std)
        #print(outline)
        #1/0
        
        for every_chromo in chromo_list:
            print('chromo', every_chromo)
            m_chromo_df = mpileup_df.loc[mpileup_df["chromo"] == every_chromo]
            print(m_chromo_df.head())
            #1/0
            fzero = m_chromo_df.replace(0, np.NaN)
            c_median = fzero["ct"].median()
            c_std = fzero["ct"].std()
            chromo_annu = c_median+c_std
            
            print(g_median, g_std, c_median, c_std, chromo_annu)
            
            start, stop = 1, len(m_chromo_df)
            outline = ('{chromo}\t{each_type}_depth\tChromosome\t{start}\t{stop}\t{c_median}\t.\t.\tID={chromo}:{start}-{stop};Median_RD={c_median};Std_RD={c_std};rel_genome_RD={rel_g};sample={sample}\n').format(chromo=every_chromo, each_type=each_type, start=start, stop=stop, c_median=c_median, c_std=c_std, rel_g=round(c_median/float(g_median),2), sample=each_sample)
            chromo_gff_outfile.write(outline)
            print(outline)
            
            if each_type == 'RD':
                chromo_depth_dict[every_chromo]=[(g_median/float(c_median)), (g_std/float(c_std))]
        
        for every_chromo in chromo_list:
            
            m_chromo_df = mpileup_df.loc[mpileup_df["chromo"] == every_chromo]
  
            if each_type == 'RD':
                fzero = m_chromo_df.replace(0, np.NaN)
                c_median = fzero["ct"].median()
                c_std = fzero["ct"].std()
            else:
                c_median = max(m_chromo_df["ct"].median(),1)*chromo_depth_dict[every_chromo][0]
                c_std = m_chromo_df["ct"].std()*chromo_depth_dict[every_chromo][1]
                
            subz_df = m_chromo_df.loc[(mpileup_df["ct"] >= (c_median + 3*c_std))]
            suby_df = m_chromo_df.loc[mpileup_df["ct"] >= (c_median + c_std)]
            
            build_trace(subz_df, suby_df, each_type, c_median, g_median)
            
            if each_type == 'RD':
                subz_df = m_chromo_df.loc[(mpileup_df["ct"] <= (c_median - 3*c_std))]
                subz_df = m_chromo_df.loc[(mpileup_df["ct"]) == 0]
                suby_df = m_chromo_df.loc[mpileup_df["ct"] <= (c_median - c_std)]
                suby_df = m_chromo_df.loc[(mpileup_df["ct"]) == 0]
                
                build_trace(subz_df, suby_df, each_type, c_median, g_median)
    
        depth_gff_outfile.close()
        chromo_gff_outfile.close()
        
        depth_pickle_name = ('{}/peaks_{}_{}.p').format(pickles_dir, each_sample, each_type)         
        with open(depth_pickle_name, 'wb') as file:
            pickle.dump(depth_dict, file)
        
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
    soft_sam_file = resource_dict['soft_sam_file']
    
    read_type_list = resource_dict['read_types']
    read_type_list.remove('RD')
        
    with open(chromo_size_pickle, 'rb') as fp:
        chromo_size_dict = pickle.load(fp)
    
    score = (resource_dict['split_zscore_median'] + resource_dict['split_zscore_std']) * sw
    print('Using split z-score: ' + str(score))
                
    if 'discordant_zscore_median' in resource_dict:
        disco_score = (resource_dict['discordant_zscore_median'] + resource_dict['discordant_zscore_std']) * dw
        print('Using disco z-score: ' + str(disco_score))
    
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

    filter_by_peaks = True

    alignment_map = {}
    hypothesis_map = {}
    chromo_map = {}

    if not args.min_coverage:
        min_coverage = 3
    else:
        min_coverage = float(args.min_coverage)
        
    if args.mapq_val:
        mapq_val = int(args.mapq_val)
    else:
        mapq_val = 30
                    
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
        #return(-1)
            
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
            
            #section handle loding range dictionary
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
                                        
        if args.load_filter:
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
        #print(runmode)
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
                    
        if runmode == 'soft':
            if len(alignments) == 1:
                chromo, start, stop, qscore, cigar = alignment_map[alignments[0]]
                hypo_id = (find_hypo(locus_to_hypo_lookup, chromo, start, stop))
                #
                breakpoint_hypothesis = ('{hid}_{hid}').format(hid=hypo_id)
                if breakpoint_hypothesis in refined_map:
                    refined_map[breakpoint_hypothesis] += [q_uid]
                else:
                    refined_map[breakpoint_hypothesis] = [q_uid]
                
            else:
                for each_instance in alignments:
                    chromo, start, stop, qscore, cigar = alignment_map[each_instance]
                    hypo_id = find_hypo(locus_to_hypo_lookup, chromo, start, stop)
                    #
                    breakpoint_hypothesis = ('{hid}_{hid}').format(hid=hypo_id)
                    if breakpoint_hypothesis in refined_map:
                        refined_map[breakpoint_hypothesis] += [q_uid]
                    else:
                        refined_map[breakpoint_hypothesis] = [q_uid]
               
        return(refined_map)
                
    for runmode in read_type_list:        
        if runmode == 'discordant':
            print('Parsing discordant reads...')
            read_file_name = discordant_sam_file        
            peaks_file_name = ('{}/{}_discordant_depth.gff').format(idx_dir, output_file)
        
        if runmode == 'split':
            print('Parsing split reads...')
            read_file_name = split_sam_file
            peaks_file_name = ('{}/{}_split_depth.gff').format(idx_dir, output_file)
            
        if runmode == 'soft':
            print('Parsing soft reads...')
            read_file_name = soft_sam_file
            peaks_file_name = ('{}/{}_soft_depth.gff').format(idx_dir, output_file)
        
        if filter_by_peaks:
            peaks_dict = load_peaks(peaks_file_name)
            
        ct = 0
        uni_uid_set = set()
        process_list = []
        read_file = open(read_file_name)

        for line in read_file:
            process = False
            ct+=1
            #if (ct % 1000 == 0):
                #print('\t'+str(ct)+' reads processed.')
            
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

                if runmode == 'split':
                    uid = line.split('\t')[0].split('_')[0]

                if qscore >= mapq_val:
                    process = True
                    
                    if filter_by_peaks:
                        process = False
                        if chromo in peaks_dict:
                            loci_set = peaks_dict[chromo]
                            
                            for nt in range(start,stop):
                                if nt in loci_set:
                                    process = True
                # if not process:
                #     print('nothing found')
                    
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
            #print(uid_source_dict)
            
            uni_uid_set = uid_source_dict[runmode]
            
            for uni_uid in uni_uid_set:
                uid = uni_uid.split('.')[0]
                
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
        
    #print(fastq_1_list)
    if len(fastq_1_list) > 0:
        fastq_1_file = open(fastq_1_filename)
        uniuid_to_seq_dict, uniuid_to_phred_dict = parse_fastq(fastq_1_filename, '.1', fastq_1_list, uniuid_to_seq_dict, uniuid_to_phred_dict)

    #print(fastq_2_list)
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
                        
        if runmode == 'soft':
            print('Running soft read mode on ' + str(len(uni_uid_set)) + ' reads ...')
        
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
                            
            #q_uid = ('{}_{}').format(runmode, uid)
    print(len(uids_that_were_filtered))
    print(uids_that_were_filtered[:6])
    #define outputs
    break_tab = ('{}/break_bt.tab').format(idx_dir)
    break_bed = ('{}/break_bd.bed').format(idx_dir)
    dbreak_tab = ('{}/dbreak_dbt.tab').format(idx_dir)
    dbreak_bed = ('{}/dbreak_dbd.bed').format(idx_dir)
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
        
        if args.load_filter:
            process_brks_1, ancestor_filter = filter_regions(chromo_l, start_l, stop_l, region_filter_dict, ancestor_filter)
            process_brks_2, ancestor_filter = filter_regions(chromo_r, start_r, stop_r, region_filter_dict, ancestor_filter)
                
        if process_brks_1 and process_brks_2:
            if len(uid_list) > 1:
                split_ct = 0
                disco_ct = 0
                soft_ct = 0
                processed_list = []
                uid_str = ''
                
                temp_uid_dict = {}

                for each in uid_list:
                    runmode, uid = each.split('_')

                    if uid not in temp_uid_dict:
                        temp_uid_dict[uid] = [runmode]
                        
                    else:
                        temp_uid_dict[uid] += [runmode]
                        
                for uid, runmodes in temp_uid_dict.items():
                    if len([i for i in runmodes if i == 'split']) > 0:
                        processed_list.append(uid)
                        split_ct += 1
                        uid_str += ('split_{},').format(uid)
                        
                    if len([i for i in runmodes if i == 'discordant']) > 0 and (uid not in processed_list):
                        disco_ct += 1
                        uid_str += ('discordant_{},').format(uid)
                        
                    if len([i for i in runmodes if i == 'soft']) > 0 and (uid not in processed_list):
                        soft_ct += 1
                        uid_str += ('soft_{},').format(uid)
                
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
    
    def populate_filter_dict(chromo, start, stop, region):
        print('pfd',chromo, start, stop, region)
        print(len(region_name_dict))
        
        if chromo not in filter_region_dict:
            filter_region_dict[chromo]=set()
    
        if chromo in filter_region_dict:
            for index in range(start,stop+1):
                filter_region_dict[chromo].add(index)  
    
        if region in region_name_dict:
            region_name_dict[region]+=1
            
        if region not in region_name_dict:
            region_name_dict[region]=1 
    
    def parse_filter_gff(gff_file):
        filter_file = open(gff_file)
        for line in filter_file:
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
                
        filter_file.close()
        
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

    if args.run_name:
        break_tab = ('{}/break_bt.tab').format(idx_dir)
        parse_tab(break_tab)
        break_bed = ('{}/break_bd.bed').format(idx_dir)
        parse_bed(break_bed)
        dbreak_tab = ('{}/dbreak_dbt.tab').format(idx_dir)
        parse_tab(dbreak_tab)
        dbreak_bed = ('{}/dbreak_dbd.bed').format(idx_dir)
        parse_bed(dbreak_bed)
    
    if args.filter_gff:
        if len(args.filter_gff) > 1:
            for each_gff_file in args.filter_gff:
                parse_filter_gff(each_gff_file)
        else:
            parse_filter_gff(args.filter_gff)
        
    if args.filter_tab:
        anc_QC_list = []
        
        for each in args.filter_tab:
            if '/' in each:
                source_name = each.rsplit('/',1)[1]
            else:
                source_name = each
            
            anc_QC_list = parse_tab(each, anc_QC_list)
                    
        print('Filtering '+str(len(anc_QC_list)))
        print(anc_QC_list)
                    
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
        
    pickle_out = ("{}{}").format(output_dir, output_file)
    pickle.dump(filter_region_dict, open(pickle_out, 'wb'))

    # pickle_out = ("{}_filter.p").format(pickles_dir)
    # pickle.dump(filter_region_dict, open(pickle_out, 'wb'))

def check_coverage(source, prefilter_file_name, contig_seq_dict):      
    process = False

    temp_filter_file = open(prefilter_file_name)
    
    for line in temp_filter_file:
        print(line)
        if line[0] == '>' and not process:
            eval_cov = float(line.rsplit('_',1)[1])
            if eval_cov >= float(min_coverage):
                    process=True
                    seq_str=''
                    node_name = source+'_'+str(int(eval_cov))
                    
        if line[0]!= '>' and process:
            seq_str += line.strip()
            
    temp_filter_file.close()

    if len(seq_str)>1:
        contig_seq_dict[node_name]=seq_str
        #addend = True
        
    print('node name '+str(node_name))
    
    return(contig_seq_dict)    
        
def gff_from_json(json_file_name, contig_seq_dict, gff_list_dict, gff_rank_dict, blast_read_aligned_sections_dict, blast_genome_aligned_regions_dict, gff_uid_dict):
    active = False
    query_deets_dict = {}
    chromo_dict={}
    try:
        data = json.load(open(json_file_name))
    except:
        print('json file ' + json_file_name + ' not found')
        return(gff_list_dict, gff_rank_dict, blast_read_aligned_sections_dict, blast_genome_aligned_regions_dict, gff_uid_dict)

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
                                                    chromo_dict[base] = str(chromo['id'])
                                                        
                                                    if (e_key == 'hsps') and (type(e_value)==list):
                                                        for e_index in range(len(e_value)):
                                                            each_hsps = e_value[e_index]
    
                                                            numb = str(base)+'.'+str(each_hsps['num'])
                                                            
                                                            if len(numb)>1:
                                                                active = True
                                                                
                                                                hit_from = int(each_hsps["hit_from"])                                                                
                                                                hit_to = int(each_hsps["hit_to"])                                                                
                                                                query_from = str(each_hsps["query_from"])                                                                
                                                                query_to = str(each_hsps["query_to"])                                                            
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
                                                                query_deets_dict[numb] = ['q_id','hit_from','hit_to','query_from','query_to','bit_score','query_strand','hit_strand','qseq', 'hseq','q_title']
                                                                query_deets_dict[numb][0] = base
                                                                query_deets_dict[numb][1] = hit_from
                                                                query_deets_dict[numb][2] = hit_to
                                                                query_deets_dict[numb][3] = query_from
                                                                query_deets_dict[numb][4] = query_to
                                                                query_deets_dict[numb][5] = bit_score
                                                                query_deets_dict[numb][6] = query_strand
                                                                query_deets_dict[numb][7] = hit_strand
                                                                query_deets_dict[numb][8] = qseq
                                                                query_deets_dict[numb][9] = hseq
                                                                query_deets_dict[numb][10] = q_title
                                                                numb = 0

    ct = 0    
    for numb, deets in query_deets_dict.items():
        q_title = query_deets_dict[numb][10]
                
        if q_title not in contig_seq_dict:
            contig_seq = ''
            print(numb, q_title, contig_seq_dict)
            1/0
        else:
            contig_seq = contig_seq_dict[q_title]

        base = query_deets_dict[numb][0]
        reduced_base = query_deets_dict[numb][0].rsplit('.', 1)[0]
        hypo_id = str(reduced_base + '_' + str(ct))

        if reduced_base in gff_uid_dict:
            gff_uid_dict[reduced_base].add(hypo_id)
        else:
            hypo_id_set = set()
            hypo_id_set.add(hypo_id)
            gff_uid_dict[reduced_base] = hypo_id_set
                
        chromo = chromo_dict[base]
        
        hit_from = query_deets_dict[numb][1]
        hit_to = query_deets_dict[numb][2]
        
        q_from = int(query_deets_dict[numb][3])
        q_to = int(query_deets_dict[numb][4])
        
        if hit_from < hit_to:
            start = hit_from
            stop = hit_to
        else:
            start = hit_to
            stop = hit_from
                    
        if query_deets_dict[numb][7] == 'Plus':
            sign = '+'
        else:
            sign = '-'
            
        bit_score = float(query_deets_dict[numb][5])
        
        if query_deets_dict[numb][6] != query_deets_dict[numb][7]:
            orient = 'reverse'
        else:
            orient = 'forward'
            
        mod_seq = ''
        qseq = query_deets_dict[numb][8]
        hseq = query_deets_dict[numb][9]
        
        for index in range(len(qseq)):
            if qseq[index] == hseq[index]:
                mod_seq += qseq[index].upper()
            if qseq[index] != hseq[index]:
                mod_seq += qseq[index].lower()
        
        node_loci = ('{},{},{},{}').format(q_title, query_deets_dict[numb][3], query_deets_dict[numb][4], float(query_deets_dict[numb][5]))
        gff_line = ('{}\terisapfel\t{}_aligned_split\t{}\t{}\t.\t{}\t{}\tnode_uid={}; orient={}; alt_seq={}; ref_seq={}; contig={}\n').format(chromo, hypo_id, start, stop, sign, int(round(bit_score)), node_loci, orient, mod_seq, hseq, contig_seq)

        aligned_region_set = set()
        
        #set internal check
        for each in range(q_from,q_to):
            aligned_region_set.add(reduced_base + '.' + str(each))
        
        #set genome wide check
        genome_region_set = set()
        for each in range(start, stop):
            genome_region_set.add(chromo + '-' + str(each))
        
        if hypo_id in blast_read_aligned_sections_dict:
            blast_read_aligned_sections_dict[hypo_id] += [aligned_region_set]
        else:
            blast_read_aligned_sections_dict[hypo_id] = [aligned_region_set]
        
        if hypo_id in blast_genome_aligned_regions_dict:
            blast_genome_aligned_regions_dict[hypo_id] += [genome_region_set]
        else:
            blast_genome_aligned_regions_dict[hypo_id] = [genome_region_set]
            
        if hypo_id in gff_rank_dict:
            gff_rank_dict[hypo_id] += [bit_score]
        else:
            gff_rank_dict[hypo_id] = [bit_score]
        
        if hypo_id in gff_list_dict:
            gff_list_dict[hypo_id] += [gff_line]
        else:
            gff_list_dict[hypo_id] = [gff_line]
            
        ct+=1

    return(gff_list_dict, gff_rank_dict, blast_read_aligned_sections_dict, blast_genome_aligned_regions_dict, gff_uid_dict)    
 
def score_uids(hypo_id, gff_line, gff_score_dict, gff_uni_uid_dict):
    chromo, efield, inline_hypo, start, stop, dot, sign, bit_score, deets = gff_line.split('\t')
    bit_score = float(bit_score)
    outline = ('{}:{}-{},{}').format(chromo, start, stop, hypo_id)
    
    if outline not in gff_score_dict:
        gff_score_dict[outline] = bit_score
        gff_uni_uid_dict[hypo_id] = gff_line
    else:
        old_score = gff_score_dict[outline]
        if old_score < bit_score:
            gff_score_dict[outline] = bit_score
            gff_uni_uid_dict[hypo_id] = gff_line
    
    return(gff_score_dict, gff_uni_uid_dict)
    
def define_alignment(gff_list_dict, gff_rank_dict, blast_read_aligned_sections_dict, blast_genome_aligned_regions_dict, gff_uid_dict, best_uid_dict):
    
    """Remove duplicates gff calls"""
    gff_score_dict = {}
    gff_uni_uid_dict = {}
    
    for each_hypo, gff_line in gff_list_dict.items():
        print(each_hypo, gff_line)
        #TODO - not sure about this hypo_id
        hypo_id = each_hypo.split('_')[0] + '_' + each_hypo.split('_')[1] + '_' + each_hypo.split('_')[2]
        print(hypo_id)
        
        if len(gff_line) < 2:
            gff_score_dict, gff_uni_uid_dict = score_uids(hypo_id, gff_line[0], gff_score_dict, gff_uni_uid_dict)
        else:
            for each_line in gff_line:
                gff_score_dict, gff_uni_uid_dict = score_uids(hypo_id, each_line, gff_score_dict, gff_uni_uid_dict)
    
    print(len(gff_list_dict))
    print(len(gff_uni_uid_dict))
    
    gff_list = []
    for each, line in gff_uni_uid_dict.items():
        gff_list.append(line)       
    
    return(gff_list)
        
def make_range_set(start, stop):
    temp_set = set()
    for each in range((start-overlap_mask),(stop+overlap_mask)):
        temp_set.add(each)
    return(temp_set)     
    
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
    
    uniuid_to_seq_dict, uniuid_to_phred_dict = io_load('fq')
        
    break_tab = ('{}/break_bt.tab').format(idx_dir)   
    dbreak_tab = ('{}/dbreak_dbt.tab').format(idx_dir)
            
    def parse_msa_file(filename):
        
        infile = open(filename)
        hash_to_name = {}
        msa_dict = {}
        seq = ''
        
        for line in infile:
            if line[0] == '>':
                name = (line.split('>')[1].split('_')[0].strip())
                hname = str(hash(name))
                
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
        #print('pmk')
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
                output_handler(subprocess.check_output([bashCommand],stderr=subprocess.STDOUT,shell=True))
                print(bashCommand)
                
                bashCommand = ('cons -name {}_{} -plurality 1 {} -outseq {}').format(node_name, len(hnames), cluster_align_name, cluster_seq_name)
                output_handler(subprocess.check_output([bashCommand],stderr=subprocess.STDOUT,shell=True))
                print(bashCommand)
                
                contig_seq_dict = check_coverage(node_name, cluster_seq_name, contig_seq_dict)
                
                prefix = ('{temp}/{node_name}').format(temp=temp_dir, node_name=node_name)
                prefix_list.append(prefix)

        return(prefix_list, contig_seq_dict)
                
    parse_this = True

    if parse_this:     
        """parse_brks:
        This recovers the loci from the hypothesized breakpoints identified in earlier
        steps.
        """    
        parse_brks(break_tab, 'seq')
        #parse_brks(dbreak_tab, 'discordant')
        
        print('Starting breakpoint testing...')
        #spades_file_name = ('{}/spades_output/scaffolds.fasta').format(temp_dir)
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
                
                bashCommand = ('blastn -task "blastn-short" -query {} -subject {fa_file} -outfmt 15 -out {}').format(cluster_seq_name, json_short_file_name, fa_file=fa_file)
                output_handler(subprocess.check_output([bashCommand],stderr=subprocess.STDOUT,shell=True))
                gff_list_dict, gff_rank_dict, blast_read_aligned_sections_dict, blast_genome_aligned_regions_dict, gff_uid_dict = gff_from_json(json_short_file_name, contig_seq_dict, gff_list_dict, gff_rank_dict, blast_read_aligned_sections_dict, blast_genome_aligned_regions_dict, gff_uid_dict)
           
                bashCommand = ('blastn -query {} -subject {fa_file} -outfmt 15 -out {}').format(cluster_seq_name, json_long_file_name, fa_file=fa_file)
                output_handler(subprocess.check_output([bashCommand],stderr=subprocess.STDOUT,shell=True))
                gff_list_dict, gff_rank_dict, blast_read_aligned_sections_dict, blast_genome_aligned_regions_dict, gff_uid_dict = gff_from_json(json_long_file_name, contig_seq_dict, gff_list_dict, gff_rank_dict, blast_read_aligned_sections_dict, blast_genome_aligned_regions_dict, gff_uid_dict)

        resource_pickle_name = ('{}/gff_list_dict.p').format(pickles_dir)
        with open(resource_pickle_name, 'wb') as file:
            pickle.dump(gff_list_dict, file)
            
        resource_pickle_name = ('{}/gff_rank_dict.p').format(pickles_dir)         
        with open(resource_pickle_name, 'wb') as file:
            pickle.dump(gff_rank_dict, file)
            
        resource_pickle_name = ('{}/blast_read_aligned_sections_dict.p').format(pickles_dir)         
        with open(resource_pickle_name, 'wb') as file:
            pickle.dump(blast_read_aligned_sections_dict, file)
            
        resource_pickle_name = ('{}/blast_genome_aligned_regions_dict.p').format(pickles_dir)         
        with open(resource_pickle_name, 'wb') as file:
            pickle.dump(blast_genome_aligned_regions_dict, file)
            
        resource_pickle_name = ('{}/gff_uid_dict.p').format(pickles_dir)         
        with open(resource_pickle_name, 'wb') as file:
            pickle.dump(gff_uid_dict, file)
    
    pickle_in = ('{}/gff_list_dict.p').format(pickles_dir)
    gff_list_dict = pickle_loader(pickle_in, 'dict')

    pickle_in = ('{}/gff_rank_dict.p').format(pickles_dir)
    gff_rank_dict = pickle_loader(pickle_in, 'dict')

    pickle_in = ('{}/blast_read_aligned_sections_dict.p').format(pickles_dir)
    blast_read_aligned_sections_dict = pickle_loader(pickle_in, 'dict')
    
    pickle_in = ('{}/blast_genome_aligned_regions_dict.p').format(pickles_dir)
    blast_genome_aligned_regions_dict = pickle_loader(pickle_in, 'dict')
        
    pickle_in = ('{}/gff_uid_dict.p').format(pickles_dir)
    gff_uid_dict = pickle_loader(pickle_in, 'dict')
        
    print('Protential breakpoints: '+ str(len(gff_uid_dict)))
    gff_list = define_alignment(gff_list_dict, gff_rank_dict, blast_read_aligned_sections_dict, blast_genome_aligned_regions_dict, gff_uid_dict, best_uid_dict)
       
    gff_list = []
    for each, sets_line in gff_list_dict.items():
        print('each', each)
        for each_line in sets_line:
            print('each_lise', each_line)
            gff_list.append(each_line)
    
    gff_file_name = ('{}/{}_realigned.gff').format(idx_dir, output_file)
    gff_file = open(gff_file_name,'w')    
    
    for each_gff in gff_list:
        gff_file.write(each_gff)
    
    gff_file.close()
    #reduced_file.close()
            
    outline = ('\t{} nodes realigned to {} regions.\n\n\tFor a complete representation refer to:\n\t\t{}').format(len(complete_qname_dict), len(gff_list), gff_file_name)
    print(outline)
    
    """ Discordant Analysis 
    In this section we reprocess the discordant hits to identify topological changes.
    """
    #parse_brks(dbreak_tab, 'discordant')
    
""" Step Seven """
""" Let's build an ml! """
if args.breakpoints:
    def get_genome_depths(resource_dict):
        rd_fg_median = float(resource_dict['RD_zscore_median'])
        rd_fg_std = float(resource_dict['RD_zscore_std'])
        split_fg_median = float(resource_dict['split_zscore_median'])
        split_fg_std = float(resource_dict['split_zscore_std'])
        discordant_fg_median = float(resource_dict['discordant_zscore_median'])
        discordant_fg_std = float(resource_dict['discordant_zscore_std'])
        soft_fg_median = float(resource_dict['soft_zscore_median'])
        soft_fg_std = float(resource_dict['soft_zscore_std'])
            
        return([rd_fg_median, rd_fg_std, split_fg_median, split_fg_std, 
                discordant_fg_median, discordant_fg_std, soft_fg_median, soft_fg_std])
            
    def calc_depth(depth_df, chromo, start, stop):
        t_median, t_std = 0, 0
        
        chromo_df = depth_df.loc[(depth_df["chromo"] == chromo)]
        fzero = chromo_df.replace(0, np.NaN)
        fc_median = fzero["ct"].median()
        #fc_std = fzero["ct"].std()
        
        target_range = depth_df.loc[(depth_df["chromo"] == chromo) & (depth_df["nuc"] >= start) & (depth_df["nuc"] <= stop)]
        
        t_median = target_range["ct"].median()
        t_std = target_range["ct"].std()
        
        rc_median = t_median/float(fc_median)
        
        return(rc_median, t_median, t_std)
        
    def parse_sam(sam_file_name, uid_set, runmode='in_uid_set'):
        if runmode == 'in_uid_set':
            sam_file = open(sam_file_name)
            uid_dict = {}
            
            for line in sam_file:
                if line[0] != '@':
                    uid, bit_score, chromo, nuc, mapq, cigar = line.split('\t')[0:6]
                    bit_score = int(bit_score)
                    start = int(nuc)
                    stop = start + parse_cigar(cigar, 'last_match')
                    if uid in uid_set:
                        if uid in uid_dict:
                            uid_dict[uid]+=[[bit_score, chromo, start, stop]]
                        else:
                            uid_dict[uid]=[[bit_score, chromo, start, stop]]
            sam_file.close()
            
            return(uid_dict)
            
        if runmode == 'feature_load':
            sam_file = open(sam_file_name)
            uid_dict = {}
            
            for line in sam_file:
                if line[0] != '@':
                    uid, bit_score, chromo, nuc = line.split('\t')[0:4]
                    bit_score = int(bit_score)
                    nuc = int(nuc)
                    
                    if uid in uid_dict:
                        uid_dict[uid] += [[bit_score, chromo, nuc]]
                    else:
                        uid_dict[uid] = [[bit_score, chromo, nuc]]
            sam_file.close()
            
            return(uid_dict)
    
        
    def build_bit_loci_map(hypo, locus, list_hypo_uids, uid_dict, path):
        print(path)
        loci_to_reads_dict = {}
        
        locus_chromo = locus.split(':')[0]
        l_start = int(locus.split(':')[1].split('-')[0])
        l_stop = int(locus.split(':')[1].split('-')[1])
        print('lstart, lstop', l_start, l_stop)
        
        for each_uid in list_hypo_uids:
            print('each_uid', each_uid)
            if each_uid in uid_dict:
                pairs = uid_dict[each_uid]
                print('each_uid, pairs', each_uid, pairs)
                if pairs > 1:
                    for each in pairs:
                        bit_score, chromo, start, stop = each
                        print('bit_score, chromo, start, stop', bit_score, chromo, start, stop)
                        if chromo == locus_chromo:
                            if not (stop < l_start) or not (start > l_stop):
                                print('success')
                                1/0
                                if locus not in loci_to_reads_dict:
                                    loci_to_reads_dict[locus] = [bit_score]
                                else:
                                    loci_to_reads_dict[locus] += [bit_score]
                else:
                    print('single')
                    1/0
                
        return(loci_to_reads_dict)
    
    def parse_stats(hypo, iter_ml_hypo_dict, stats):
        locus, contig, list_hypo_uids = stats
        
        disco_loci_to_reads_dict = build_bit_loci_map(hypo, locus, list_hypo_uids, disco_uid_dict, 'discordant')
        split_loci_to_reads_dict = build_bit_loci_map(hypo, locus, list_hypo_uids, split_uid_dict, 'split')
        soft_loci_to_reads_dict = build_bit_loci_map(hypo, locus, list_hypo_uids, soft_uid_dict, 'soft')
        
        iter_ml_hypo_dict[hypo] += [disco_loci_to_reads_dict, split_loci_to_reads_dict, soft_loci_to_reads_dict]
        print(len(disco_loci_to_reads_dict),len(split_loci_to_reads_dict),len(soft_loci_to_reads_dict))
        if (len(disco_loci_to_reads_dict) > 1):
            print('yep disco')
            1/0
        if (len(soft_loci_to_reads_dict) > 1):
            print('yep soft')
            1/0 
    
        return(ml_hypo_dict)
        
    def parse_vcf(vcf_name):
    
        vcf_file = open(vcf_name)
        vcf_locus_dict = {}
        vcf_id_dict = {}
        
        for line in vcf_file:
            if line[0]!='#':
                chromo = line.split('\t')[0]
                start = int(line.split('\t')[1])
                 
                v_id = int(line.split('\t')[2])
                
                nstart = start - 150
                nstop = start + 150
                locus = ('{}:{}-{}').format(chromo, nstart, nstop)
                
                if chromo not in vcf_locus_dict:
                    vcf_locus_dict[chromo] = [[nstart, nstop]]
                else:
                    vcf_locus_dict[chromo] += [[nstart, nstop]]
                    
                if v_id not in vcf_id_dict:
                    vcf_id_dict[locus] = v_id
                
        vcf_file.close()
        
        return(vcf_locus_dict, vcf_id_dict)
        
    """
    
    """   

    resource_dict = io_load()
     
    fa_file = resource_dict['genome_fa']
    each_sample = resource_dict['run_name']
    read_type_list = resource_dict['read_types']
    
    fg_depths = get_genome_depths(resource_dict)

    do_this = True
    if do_this:    
        # Load refined_map
        pickle_out = ("{}/refined_map.p").format(pickles_dir)
        refined_map = pickle.load(open(pickle_out))
        uid_set = set()
            
        """ First, build the loci information
        Each line in the gff contains a single locus. We want to extract all read
        all read depth measures as well as the reference_sequence of these regions.
        Some of these reads have already been processed, we can load the previous peak data for those files that have already been measured
        This may be a chunk of code but it saves some time.
        
        For those loci that haven't been preocessed yet we need to load the relevant mpileups and and take those measures.
        """
            
        depth_pickle_name = ('{}/peaks_{}_RD.p').format(pickles_dir, each_sample)
        rd_depth_dict = pickle_loader(depth_pickle_name, 'dict')
        
        depth_pickle_name = ('{}/mpileup_{}_RD.p').format(pickles_dir, each_sample)
        rd_mpileup_df = pickle_loader(depth_pickle_name, 'df')
    
        depth_pickle_name = ('{}/peaks_{}_split.p').format(pickles_dir, each_sample)
        split_depth_dict = pickle_loader(depth_pickle_name, 'dict')
        
        depth_pickle_name = ('{}/mpileup_{}_split.p').format(pickles_dir, each_sample)
        split_mpileup_df = pickle_loader(depth_pickle_name, 'df')
        
        depth_pickle_name = ('{}/peaks_{}_discordant.p').format(pickles_dir, each_sample)
        discordant_depth_dict = pickle_loader(depth_pickle_name, 'dict')
        
        depth_pickle_name = ('{}/mpileup_{}_discordant.p').format(pickles_dir, each_sample)
        discordant_mpileup_df = pickle_loader(depth_pickle_name, 'df')
    
        depth_pickle_name = ('{}/peaks_{}_soft.p').format(pickles_dir, each_sample)
        soft_depth_dict = pickle_loader(depth_pickle_name, 'dict')
        
        depth_pickle_name = ('{}/mpileup_{}_soft.p').format(pickles_dir, each_sample)
        soft_mpileup_df = pickle_loader(depth_pickle_name, 'df')
        
        line_ct = 0
        dict_ct = 0
        df_ct = 0
        loci_to_reads_dict = {}
        
        gff_file_name = ('{}/{}_realigned.gff').format(idx_dir, output_file)
        gff_file = open(gff_file_name)
        for line in gff_file:
            print(line)
            line_ct += 1 
            if line[0] != '#':
                #print(line_ct)
                #NC_001147.6     erisapfel       hypo:267_774_right_8622470508732422200_9_aligned_split  1091195 1091286 .       +       183     node_uid=hypo:267_774_right_8622470508732422200,11,102,182.869; orient=forward; alt_seq=ATTTTCAttttttttttttAATTTCGGTCAGAAAGCCGGGTAAGGTATGACAGCGAGAGTAGAGGTAGATGTGAGAGAgtgtgtgggtgtgg; ref_seq=ATTTTCATTTTTTTTTTTTAATTTCGGTCAGAAAGCCGGGTAAGGTATGACAGCGAGAGTAGAGGTAGATGTGAGAGAGTGTGTGGGTGTGG; contig=ACATGGGTTTATTTTCATTTTTTTTTTTTAATTTCGGTCAGAAAGCCGGGTAAGGTATGACAGCGAGAGTAGAGGTAGATGTGAGAGAGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTG;type=",transposition,transposition,ODIRA_CNV,transposition,transposition,transposition,transposition,transposition,local_rearrangement_or_CNV,transposition,transposition,transposition,ODIRA_CNV,transposition,transposition,local_rearrangement_or_CNV,transposition,transposition,local_rearrangement_or_CNV,transposition,transposition,transposition,transposition,ODIRA_CNV,transposition,local_rearrangement_or_CNV,transposition,transposition,transposition,transposition,transposition,ODIRA_CNV,transposition,ODIRA_CNV,transposition,transposition,transposition,local_rearrangement_or_CNV,ODIRA_CNV,transposition,transposition"
    
                chromo, source, hypo, start, stop, dot, dash, score, deets = line.split('\t')
                
                node_uid, orient, alt_seq, ref_seq, contig = deets.split(';')
                start = int(start)
                stop = int(stop)
                locus = ('{}:{}-{}').format(chromo, start, stop)
                
                hypo_single = hypo.split('_')[0] + '_' + hypo.split('_')[1]
                
                uid_raw = refined_map[hypo_single]
                list_hypo_uids = set()
                for each in uid_raw:
                    uid_set.add(each)
                    uid_set.add(each.split('_')[0])
                    list_hypo_uids.add(each)
                    list_hypo_uids.add(each.split('_')[0])
                    
                if 'reverse' in orient:
                    loci_seq = reverse_compliment(ref_seq.split('=')[1])
                else:
                    loci_seq = ref_seq.split('=')[1]
                    
                #swapped out hypo single
                if hypo_single not in ml_hypo_dict:
                    ml_hypo_dict[hypo_single] = [[locus, contig, list_hypo_uids]]
                else:
                    ml_hypo_dict[hypo_single] += [[locus, contig, list_hypo_uids]]
                    
                if hypo_single not in ml_san_check:
                    ml_san_check[hypo_single] = [[line, locus, contig, list_hypo_uids]]
                else:
                    ml_san_check[hypo_single] += [[line, locus, contig, list_hypo_uids]]                    
                
                
                if locus not in ml_loci_dict:               
                    if locus in rd_depth_dict:
                        rd_rel_c, rd_rel_g, rd_t_std = rd_depth_dict[locus]
                        dict_ct+=1
                    else:
                        rd_rel_c, raw_g, rd_t_std = calc_depth(rd_mpileup_df, chromo, start, stop)
                        rd_rel_g = raw_g / float(fg_depths[0])
                        df_ct+=1
                    #print('rd')
        
                    if locus in split_depth_dict:
                        sp_rel_c, sp_rel_g, sp_t_std = split_depth_dict[locus]
                        dict_ct+=1
                    else:
                        sp_rel_c, raw_g, sp_t_std = calc_depth(split_mpileup_df, chromo, start, stop)
                        sp_rel_g = raw_g / float(fg_depths[2])
                        df_ct+=1
                    #print('sp')
        
                    if locus in discordant_depth_dict:
                        ds_rel_c, ds_rel_g, ds_t_std = discordant_depth_dict[locus]
                        dict_ct+=1
                    else:
                        ds_rel_c, raw_g, ds_t_std = calc_depth(discordant_mpileup_df, chromo, start, stop)
                        ds_rel_g = raw_g / float(fg_depths[4])
                        df_ct+=1
                    #print('ds')
                        
                    if locus in soft_depth_dict:
                        so_rel_c, so_rel_g, so_t_std = soft_depth_dict[locus]
                        dict_ct+=1
                    else:
                        so_rel_c, raw_g, so_t_std = calc_depth(soft_mpileup_df, chromo, start, stop)
                        so_rel_g = raw_g / float(fg_depths[6])
                        df_ct+=1
                    #print('so')
                
                    ml_loci_dict[locus]=[rd_rel_c, rd_rel_g, rd_t_std, fg_depths[1],
                                 sp_rel_c, sp_rel_g, sp_t_std, fg_depths[3],
                                 ds_rel_c, ds_rel_g, ds_t_std, fg_depths[5],
                                 so_rel_c, so_rel_g, so_t_std, fg_depths[7], loci_seq]
        gff_file.close()
        
        print('this far')
        
        resource_pickle_name = ('{}/uid_set.p').format(idx_dir)         
        with open(resource_pickle_name, 'wb') as file:
            pickle.dump(uid_set, file)
        
        io_append(ml_loci_dict, 'ml_loci')
        io_append(ml_hypo_dict, 'ml_hypo')
    
    """Second, build the hypo information
    Unlike the locus information the hypo data is specific to certian connected reads not individual locations.
    To get information about these reads we need to keep the connection between them intact.
    Each line has a hypo which should be present in the refined_map which stores the uid_list. 
    These uid's allow us to reach into the sam files and extract the bit encoded information. 
    These bit encodings can (i hope) instuct us as to the nature of the underlying structure.
    
    Beyond that each line has a locus which can be used to look up the locus stats in ml_loci_dict and the aligned contig.
    """

    if not do_this:
        ml_hypo_dict, ml_loci_dict = io_load('ml')
    
    resource_pickle_name = ('{}/uid_set.p').format(idx_dir)
    with open(resource_pickle_name, 'rb') as uid_file:
        uid_set = pickle.load(uid_file)
    
    pickle_out = ("{}/refined_map.p").format(pickles_dir)
    refined_map = pickle.load(open(pickle_out))

#    refct = 0
#    allrefset = set()
#    hitrefset = set()
#    is_ref = open('refined_map.tab','w')
#    for hypo, uid_list in refined_map.items():
#        outline = ('{}\t{}\n').format(hypo, uid_list)
#        is_ref.write(outline)
#        for each in uid_list:
#            allrefset.add(each)            
#        if len(uid_list) > 1:
#            refct += 1
#            hitrefset.add(each)
#    is_ref.close()
    
    base_filename=('{}/{}').format(bam_dir, output_file)
    for each_type in read_type_list:
        if each_type == 'discordant':
            #sam_file_name = base_filename + '_discordant.sam'
            #disco_uid_dict = parse_sam(sam_file_name, uid_set)
            # TODO shift to the dictionary type
            disco_uid_dict = parse_sam(resource_dict['discordant_sam_file'], uid_set) 
        if each_type == 'split':
            #sam_file_name = base_filename + '_split.sam'
            #split_uid_dict = parse_sam(sam_file_name, uid_set)      
            # TODO shift to the dictionary type
            split_uid_dict = parse_sam(resource_dict['split_sam_file'], uid_set) 
        if each_type == 'soft':
            soft_uid_dict = parse_sam(resource_dict['soft_sam_file'], uid_set)
            
    #print('outside disco_uid_dict', len(disco_uid_dict))
    #print('outside split_uid_dict', len(split_uid_dict))
    #print(len(soft_uid_dict))
    
    iter_ml_hypo_dict = ml_hypo_dict.copy()

    for hypo, all_stats in ml_hypo_dict.items():
        print('ml', hypo, ml_san_check[hypo])
        print('as', hypo, all_stats)
        for index in range(len(all_stats)):
            stats = all_stats[index]                    
            iter_ml_hypo_dict = parse_stats(hypo, iter_ml_hypo_dict, stats)
            print()
        
#        if len(all_stats) > 1:
#            print('all_stats', all_stats)
#            for index in range(len(all_stats)):
#                stats = all_stats[index]
#                if len(stats) < 2:
#                    print(index, all_stats)
#                    print(stats)
#                    
#                ml_hypo_dict = parse_stats(hypo, ml_hypo_dict, stats, len(all_stats))
#        else:
##            print('small', all_stats)
#            ml_hypo_dict = parse_stats(hypo, ml_hypo_dict, all_stats, len(all_stats))

    #print('outside uid_list', list_hypo_uids)
    #print('outside disco_uid_dict', disco_uid_dict)
    #print(soft_uid_dict)
        
    #print('485_485', ml_hypo_dict['485_485'])
        
    print('nope')
    1/0

    io_append(ml_hypo_dict, 'ml_hypo')
    
    """Third, build the connections 
    """   
    
    ml_hypo_dict, ml_loci_dict = io_load('ml')
    
    #print(len(ml_hypo_dict), len(ml_loci_dict))
    
    for hypo, deets in ml_hypo_dict.items():
        
        locus, contig, uid_list, dict_list = deets
        disco_loci_to_reads_dict, split_loci_to_reads_dict, soft_loci_to_reads_dict = dict_list
        #print(hypo)
        #print(locus, contig, uid_list, disco_loci_to_reads_dict, split_loci_to_reads_dict, soft_loci_to_reads_dict)
        #1/0

  
    

    
    sam_file_name = 'C:\\Gresham\\MinIon\\xxx.sam'
    vcf_name = 'C:\\Gresham\\Project_Erisapfel\\bam\\BC03_1731_sniffles.vcf'
    
    uid_dict = parse_sam(sam_file_name, '', 'feature_load')
    vcf_locus_dict, vcf_id_dict = parse_vcf(vcf_name)
    
    vcf_scores = {}
    vcf_total = 0
    out_scores = {}
    out_total = 0
    
    for uid, deets in uid_dict.items():
        if len(deets) > 1:
            for each_d in deets:
                bs, ch, nu = each_d
                
                if ch in vcf_locus_dict:
                    ranges = vcf_locus_dict[ch]
                    for each_r in ranges:
                        if (nu >= each_r[0]) and (nu <= each_r[1]):
                            vcf_total+=1
                            if bs not in vcf_scores:
                                vcf_scores[bs] = 1
                            else:
                                vcf_scores[bs] +=1  
                        else:
                            out_total+=1
                            if bs not in out_scores:
                                out_scores[bs] = 1
                            else:
                                out_scores[bs] +=1
                                
    #print(vcf_scores)
    #print(out_scores)
    
    sig_list = []
    
    for bs, ct in out_scores.items():
        if bs in vcf_scores:    
            vcf_ct = vcf_scores[bs]        
            oddsratio, pvalue = stats.fisher_exact([[vcf_ct, vcf_total], [ct, out_total]])
            if (pvalue < 0.05) and (oddsratio >= 2) and min(ct, vcf_ct) > 5:
                print(bs, oddsratio, pvalue)
                sig_list.append([bs, oddsratio, pvalue, vcf_ct, ct])
    
    
    """
    End of script
    """
