# -*- coding: utf-8 -*-
"""
CVish
# Purpose: Uses split and discordant reads to identify potetnial genomic breakpoints
# further processing is used to extract context and sequence.

#change log
(Initial Release - Public Beta: Down to Earth)
Version 1.0 2023.05.09 (arderley)
    Rename project from 'erisapfel' to 'CVish'
    _x_ Add template() function
    _x_ Add config() function
    _x_ Improved sequence clustering
    _x_ Improved breakpoint overlap reduction
    _x_ Added gff dependant powered feature depth calculator 
    _x_ update help()
    _x_ update demo()
    _x_ update test()
    _x_ update run()
    
Version 1.1 2023.05.20 (mitaclau)
    _x_ Check contig reporting in gff, vcf
    _x_ Add CLI score filter so lines aren't added to the gff and filtered from the vcf
    _x_ Add clean() function
    _x_ Add disco optimization improvement
    
Version 1.2 2023.07.07 (canorting)
    _x_ ploidy option
    _x_ runmode option
    _x_ bwa index, note in document that this is not bwa-mem2
    _x_ removed orphaned bcftools, tabix
    _x_ added cn_weights and sv_weights to final score 

Version 1.3 2024.01.29 (heimaless)
    _x_ fixed sparse discordant reads normalization (-call) 

Future versions:
    ___ Use better demo data (chromo VI, XI)
    ___ Stop double gunzipping fastq.gz
    ___ Progress tracker
    ___ Add verbose() function
    ___ Unify and import filters further upstream in analysis 
    ___ Add gff dependant filter maker
    ___ Implement Single End option

@author: Pieter Spealman ps163@nyu.edu
"""

"""Requirements:
bwa         0.7.17
samtools    1.14
bedtools    2.29.2
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
from datetime import date
import scipy.stats as stats

import warnings
warnings.simplefilter("ignore")

#function handles help commands
def help_dialog():
    monolog=('Manual for cvish\n'+
    '#=============================================================#\n'+
    'Created by Pieter Spealman\n ps163@nyu.edu \n'+
    'Release version: 1.2 \n'+
    'Release date: 06.07.2023 \n'+
    '\tThis is a multistep package that identifies CNV and  breakpoints using split\n'+
    'and discordant sequencing reads.\n'+
    '#=============================================================#\n'+
    'Usage:\n\tCopy number and structural variant identification.\n'+
    'This method should be used for experiments involving evolved, adapted, or\n'+
    'normal/tumor samples.\n'+
    'For more information refer:github.com/pspealman/CVish\n'+
    'For demonstration use:\n\t python cvish.py -demo\n'+
    'To run a test using defaults, use:\n\t python cvish.py -test\n'+
    '')
    print(monolog)

def demo():
    monolog = ('\n\nUsage:\n\tCNV identification with ancestor.\n'+
    'This method should be used for experiments involving evolved, adapted, or\n'+
    'normal/tumor samples.\n'+
    'Example:\n\n'+
    'Step 1. Make breakpoint predictions in ancestor:\n'+
    '    python cvish.py -run -fa demo/demo.fna -fastq_1 demo/n01_ancestor.fastq.gz -fastq_2 demo/n02_ancestor.fastq.gz -config demo/demo_config.tsv -run_name demo_anc\n\n'+
    'Step 2. Verify ancestor predictions. \n'+
    '    #NOTE: By default the output is located in results/<run_name>/output/ so :\n'+
    '    nano results/demo_anc/output/demo_anc_SV_CNV.gff\n\n'+
    'Step 3. Make breakpoint predictions in evolved:\n'+
    '    #NOTE: we filter the regions already identified in the ancestor using "-filter_gff results/demo_anc/output/demo_anc_SV_CNV.gff"\n'+
    '    python cvish.py -run -fa demo/demo.fna -fastq_1 demo/n01_evolved.fastq.gz -fastq_2 demo/n02_evolved.fastq.gz -config demo/demo_config.tsv -run_name demo_evo -filter_gff results/demo_anc/output/demo_anc_SV_CNV.gff\n\n'+
    'Step 4. Verify evolved predictions. \n'+
    '    nano results/demo_evo/output/demo_evo_SV_CNV.gff')
    print(monolog)
    
def convert_bool(x):
    if isinstance(x, str):
        y = x.lower()
        if y in ['t', 'true', '1', 'y', 'yes']:
            return(True)
        
        if y in ['f', 'false', '0', 'none', 'no', 'na']:
            return(False)
        
        return(x)
    else:
        return(x)
    
def eval_type(resource_dict, param, value):
    bool_list = ['ref_gff', 'filter_object', 'filter_gff', 'filter_bed',
                 'filter_chromosome', 'depth_region_filter', 'with_disco',
                 'ref_gff_feature', 'gff_feature_name',
                 'sbatch_filename', 'module_filename', 'verbose',
                 'no_post_run_clean_up', 'skip_fasta_index']
    
    int_list = ['mapq_val', 'filter_flanking', 'split_score', 'disco_score', 
                'max_gap', 'resolution_gap', 'cnv_min_length', 'expected_ploidy']
    
    float_list = ['split_weight', 'disco_weight', 'max_eval', 'min_confidence_score',
                  'min_purity']
        
    if param in resource_dict:
        print('encountered duplicate parameter name: ', resource_dict)
        return(resource_dict)
    
    if param not in resource_dict:
        if param in bool_list:
            resource_dict[param] = convert_bool(value)
            return(resource_dict)
        
        if param in int_list:
            resource_dict[param] = int(value)
            return(resource_dict)
        
        if param in float_list:
            resource_dict[param] = float(value)
            return(resource_dict)
        
        resource_dict[param]=value
        return(resource_dict)
    
def handle_outfile(p_output):
    #TODO check that this works
    if '/' in p_output:
        output_dir = str(p_output.rsplit('/',1)[0]+'/')
        output_file = p_output.rsplit('/',1)[1]
        
        if len(output_dir.strip())>=1:
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
    else:
        output_dir = ''
        output_file = p_output
        
    return(output_dir, output_file)  

def io_initialize():
    #global resource_dict
    
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

def io_make(runmode='resource', outfile_name='na', config_dict = {}):    
    if runmode == 'resource':
        outline = ('Files and parameters used are stored in {}/resource.tab').format(final_output_dir)
        print(outline)
        
        resource_file_name = ('{}/resource.tab').format(final_output_dir)
        resource_file = open(resource_file_name, 'w')
        
        resource_file.close()
        
    if runmode == 'template':
        outline = ('Configuration template and default parameters used are stored in {}').format(outfile_name)
        print(outline)
        
        resource_file_name = ('{}').format(outfile_name)
        resource_file = open(resource_file_name, 'w')
        
        header = ('#parameter_name\tvalue\n')
        resource_file.write(header)
        
        for param in config_dict:
            outline = ('{param}\t{value}\n').format(param = param,
                                                    value = config_dict[param])
            resource_file.write(outline)
            
        resource_file.close()
        
    return(resource_file_name)
                                     
def test():
    #TODO make accurate for latest version
    monolog = ('=== Currently testing cvish.py ===')
    print(monolog)
    
    monolog = ('\tTesting Step 1. -run command for ancestor\n')
    print(monolog)    
    bashCommand = ('python cvish.py -run -fa demo/demo.fna -fastq_1 demo/n01_ancestor.fastq.gz -fastq_2 demo/n02_ancestor.fastq.gz -config demo/demo_config.tsv -run_name demo_anc')
    print(bashCommand)       
    subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)
    
    monolog = ('\tTesting Step 2. -run command for evolved\n')
    print(monolog)    
    bashCommand = ('python cvish.py -run -fa demo/demo.fna -fastq_1 demo/n01_evolved.fastq.gz -fastq_2 demo/n02_evolved.fastq.gz -config demo/demo_evo_config.tsv -filter_gff results/demo_anc/output/demo_anc_SV_CNV.gff -run_name demo_evo')
    print(bashCommand)       
    subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)
    
    
'''
# PARSER BLOCK
'''
parser = argparse.ArgumentParser(add_help=True)

''' Help '''
parser.add_argument('-man',"--manual", help="Manual screen", action='store_true')
parser.add_argument('-demo',"--demo", help="Run Demo", action='store_true')
parser.add_argument('-test',"--test", help="Run Test", action='store_true')

''' Commands '''
parser.add_argument('-skip_index',"--skip_fasta_index", 
                    help = "use bwa to index reference fasta file",
                    default = False,
                    action='store_true')

parser.add_argument('-template',"--make_template_file", 
                    help = "Generate tab-separated config template")

parser.add_argument('-config',"--load_configuration_file", 
                    help = "Define path to the configuration file to use")

parser.add_argument('-cli_config',"--configuration", 
                    help = "Set / Modify configuration parameters from the command line",
                    action='store_true')

parser.add_argument('-load',"--load_sequences",
                    help = "Load fastq files, parse into splits and breaks ",
                    action='store_true')

parser.add_argument('-depth', '--depth_analysis', 
                    help = "Calculate read depth and stats",
                    action='store_true')

parser.add_argument('-peaks', '--peaks', 
                    help = "Calculate read peaks using depth stats",
                    action='store_true')

parser.add_argument('-filter', '--make_filter', 
                    help = "[optional] Region filter command,"
                    "must be used with -filter_gff and or -filter_bed",
                    default = False,
                    action='store_true')

parser.add_argument('-map', '--map_reads',
                    help = "Maps reads into contigs",
                    action='store_true')

parser.add_argument('-call', '--call_breakpoints',
                    help = "Call breakpoints using reduced contigs", 
                    action='store_true')

''' Required '''
parser.add_argument('-fa',"--fa_file", help = 'Reference genome fasta file')
parser.add_argument('-fastq_1', "--fastq_1", help = 'Single file, R1 of paired reads')
parser.add_argument('-fastq_2', "--fastq_2", help = 'Single file, R1 of paired reads')
parser.add_argument('-run_name', '--run_name', help = 'Will be run identifier and results directory')

''' Optional '''
# DEBUG options
parser.add_argument('-verbose', "--verbose", 
                    help = '[DEBUG] print optional text to screen',
                    default = False,
                    action='store_true')

parser.add_argument('-read_type','--read_type_list', 
                    help = '[DEBUG] Read types to be used in run (use for debug purposes)',
                    default=set(['RD','discordant','split']),
                    nargs='+')

parser.add_argument('-temp_dir','--temp_dir', 
                    help = '[DEBUG] Read types to be used in run (use for debug purposes)',
                    default=set(['RD','discordant','split']),
                    nargs='+')


parser.add_argument('-no_clean','--no_post_run_clean_up', 
                    help = '[DEBUG] Prevent the deletion of intermediary files',
                    default = False,
                    action='store_true')

# GFF objects
parser.add_argument('-gff', "--ref_gff_file", 
                    help = '[optional] reference genome gff for relative depth calls',
                    default = False)

parser.add_argument('-gff_feature', "--ref_gff_feature",
                    help = '[optional, --ref_gff_file] gff feature type to use for depth calls',
                    default = 'CDS')

# parser.add_argument('-gff_feature_name', '--gff_feature_name',
#                     help = '[optional, --ref_gff_file] text used to extract feature name from details column of gff',
#                     default = 'CDS:')

parser.add_argument('-gff_feature_name', '--gff_feature_name',
                    help = '[optional, --ref_gff_file] text used to extract feature name from details column of gff',
                    default = 'gene=')

# Filter objects
parser.add_argument('-fo', "--filter_object",
                    help = "[optional, --make_filter] output file name of 'make_filter' command",
                    default = False)

parser.add_argument('-filter_gff','--filter_gff',
                    default = False,
                    help = "[optional, --make_filter] gff input file name for 'make_filter' command",
                    nargs='+')

parser.add_argument('-filter_bed', '--filter_bed',
                    help = "[optional, --make_filter] bed input file name for 'make_filter' command",
                    default = False,
                    nargs='+')

parser.add_argument('-filter_chr', '--filter_chromosome', 
                    help = "[optional, --load_sequences] a space seperated list for chromosome to include during 'load_sequences'",
                    default = False, 
                    nargs='+')

parser.add_argument('-depth_filter', '--depth_filter_bed',
                    help = "[optional] bed file with regions to include during 'depth_analysis'",
                    default = False)

# Identification weights
parser.add_argument('-qval', '--mapq_val',
                    help = 'aligned minimum MAPQ quality score',
                    default = 50)

parser.add_argument('-flank','--filter_flanking',
                    help = 'flanking size (bp) for peak calls',
                    default = 75)

parser.add_argument('-with_disco', '--with_disco',
                    help = 'include discordant reads in calculations',
                    default=True, action='store_false')

parser.add_argument('-split_score', "--split_score",
                    help = 'Minimum score needed by split reads to trigger sensitivity calls',
                    default = 9)

parser.add_argument('-disco_score', "--disco_score", 
                    help = 'Minimum score needed by discordant reads to trigger sensitivity calls',
                    default = 27)

parser.add_argument('-sw', "--split_weight", 
                    help = 'Weighted value for a split read',
                    default = 3, type=float)

parser.add_argument('-dw', "--disco_weight", 
                    help = 'Weighted value for a discordant read',
                    default = 1, type=float)

parser.add_argument('-max_eval', "--max_eval", 
                    help="Max Eval for blast alignment for recombination",
                    default = 0.01, type=float)

parser.add_argument('-min_score','--min_confidence_score',
                    help = 'Called breakpoints with values less than this are not output',
                    default=1, 
                    type=float)

parser.add_argument('-min_length','--cnv_min_length',
                    help = 'Minimum length of region to evaluate for variants',
                    default=300, 
                    type=int)

parser.add_argument('-purity', '--min_purity', 
                    help = 'Minimum percentage of sensitivity hits versus misses over a candidate peak region',
                    default = 0.9, type=float)

parser.add_argument('-gap', '--max_gap',
                    help = 'Maximum distance for peak sensitive combinations',
                    default = 300, type=int)

parser.add_argument('-rgap', '--resolution_gap',
                    help="Max distance between sequence fragments to be combined", 
                    default = 5, type=int)

parser.add_argument('-high', '--high_sensitivity_mode',
                    help = "[optional] run with high sensitivity (and high FDR) parameters'",
                    default = False)

parser.add_argument('-low', '--low_sensitivity_mode',
                    help = "[optional] run with low sensitivity (and low FDR) parameters'",
                    default = False)

parser.add_argument('-ploidy', "--expected_ploidy", 
                    help="[optional] This effects the '_feature_copy_number.tsv' output by dividing the relative copy number by the ploidy number",
                    default = 1, type=int)


''' One Line Run '''
parser.add_argument('-run',"--run", help="Single line command", action='store_true')
parser.add_argument('-sbatch', '--sbatch_filename', help="[optional, --run] Run as SBATCH", default = False)
parser.add_argument('-module', '--module_filename', help="[optional, --run] Load Modules", default = False)
parser.add_argument('-no_run',"--no_run", 
                    help="[optional, --run] Generates shell script file but does not execute it",
                    default = False,
                    action='store_true')

args = parser.parse_args()

if args.manual:
    help_dialog()
    
if args.demo:
    demo()
    
if args.test:
    test()
                    
''' globals '''        
filter_region_dict = {}
region_name_dict = {}
contig_seq_dict = {}

compliment_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
global_hypothesis_dict = {}
global_contig_dict = {}

if args.run_name:
    make_name = ("results/{}").format(args.run_name)
    output_dir, output_file = handle_outfile(make_name)   
    s_path, bam_dir, final_output_dir, temp_dir, pickles_dir = io_initialize()
    
''' 
Function definition block
'''
def cleanup_chromo_set(is_str):
    #Cleaning up the string and set jumble
    is_str = is_str.strip().replace("',","").replace("{","").replace("}","").replace("'","")
    
    #handling different seperators
    if (', ' in is_str) or (',' in is_str) or (' ' in is_str):
        if ', ' in is_str:
            is_list = is_str.split(', ')
        
        if ',' in is_str:
            is_list = is_str.split(',')
        
        if ' ' in is_str:
            is_list = is_str.split(' ')
            
        return(set(is_list))
    
    return(set([is_str]))

def load_configuration_file(config_file_name):
    config_file = open(config_file_name)
                
    resource_dict = {}
    for line in config_file:
        if line[0] != '#':
            line = line.strip()
            param = line.split('\t')[0]
            value = line.split('\t')[1]
            
            resource_dict = eval_type(resource_dict, param, value)
    
    resource_file_name = io_make()
    print('Configuration file loaded ', resource_file_name)
    
    if not args.run_name:
        print('Each run requires a unique run name. Use -run_name <name>\n')
        
    resource_dict['run_name'] = args.run_name
                    
    if 'results_dir' not in resource_dict:
        resource_dict['results_dir']=resource_dict['run_name']
            
    if args.fa_file:
        resource_dict['genome_fa']=args.fa_file
            
    if args.fastq_1:
        resource_dict['fastq_1']=args.fastq_1
    
    if args.fastq_2:
        resource_dict['fastq_2']=args.fastq_2
        
    if args.filter_chromosome:
        resource_dict['filter_chromosome'] = args.filter_chromosome
            
    if resource_dict['filter_chromosome']:
        if not isinstance(resource_dict['filter_chromosome'], set):
            if not isinstance(resource_dict['filter_chromosome'], list):
                clean_fc = cleanup_chromo_set(resource_dict['filter_chromosome'])                
                resource_dict['filter_chromosome'] = clean_fc
                
            if isinstance(resource_dict['filter_chromosome'], list):
                clean_fc = resource_dict['filter_chromosome']
                resource_dict['filter_chromosome'] = set(clean_fc)             
                    
    if args.with_disco:
        resource_dict['with_disco'] = args.with_disco
        
    if args.no_post_run_clean_up:
        resource_dict['no_post_run_clean_up'] = args.no_post_run_clean_up
        
    if 'no_post_run_clean_up' not in resource_dict:
        resource_dict['no_post_run_clean_up'] = args.filter_object
        
    if resource_dict['high_sensitivity_mode']:
        resource_dict['mapq_val'] = 10
        resource_dict['filter_flanking'] = 100
        resource_dict['split_score'] = 3
        resource_dict['disco_score'] = 9
        resource_dict['max_eval'] = 0.05
        resource_dict['min_confidence_score'] = 0.5
        resource_dict['cnv_min_length'] = 100
        resource_dict['min_purity'] = 0.8

    if resource_dict['low_sensitivity_mode']:
        resource_dict['mapq_val'] = 50
        resource_dict['filter_flanking'] = 50
        resource_dict['split_score'] = 18
        resource_dict['disco_score'] = 42
        resource_dict['max_eval'] = 0.001
        resource_dict['min_confidence_score'] = 1.5
        resource_dict['cnv_min_length'] = 400
        resource_dict['min_purity'] = 0.9
 
    print('### Configuration Parameters ###')
    for param in resource_dict:
        print(param, resource_dict[param], type(resource_dict[param]))
    print('###')
  
    #Add for DEBUG
    resource_dict['read_types'] = set(['RD', 'split', 'discordant'])
        
    io_append(resource_dict)
    
    return(resource_dict)

def generate_config_file():
    #initialize resources object and file:
    resource_file_name = io_make()
    #Required parameters:
    resource_dict = {'with_disco':args.with_disco,
                     'depth_region_filter':args.depth_filter_bed,
                     'ref_gff':args.ref_gff_file,
                     'ref_gff_feature':args.ref_gff_feature,
                     'gff_feature_name':args.gff_feature_name,
                     'module_filename':args.module_filename,
                     'max_gap':args.max_gap,
                     'min_purity':args.min_purity,
                     'filter_gff':args.filter_gff,
                     'filter_bed':args.filter_bed,
                     'filter_object':args.filter_object,
                     'resolution_gap':args.resolution_gap,
                     'max_eval':args.max_eval,
                     'mapq_val':args.mapq_val,
                     'filter_flanking':args.filter_flanking,
                     'split_score':args.split_score,
                     'disco_score':args.disco_score,
                     'split_weight':args.split_weight,
                     'disco_weight':args.disco_weight,
                     'sbatch_filename':args.sbatch_filename,
                     'min_confidence_score':args.min_confidence_score,
                     'cnv_min_length':args.cnv_min_length,
                     'no_post_run_clean_up':args.no_post_run_clean_up,
                     'expected_ploidy':args.expected_ploidy,
                     'high_sensitivity_mode':args.high_sensitivity_mode,
                     'low_sensitivity_mode':args.low_sensitivity_mode,
                     'skip_fasta_index':args.skip_fasta_index
                     }
    
    resource_dict['run_name']=args.run_name
    resource_dict['results_dir']=args.run_name
    resource_dict['genome_fa']=args.fa_file
    resource_dict['fastq_1']=args.fastq_1
    resource_dict['fastq_2']=args.fastq_2
    
    if args.filter_chromosome:
        if not isinstance(args.filter_chromosome,set):
            resource_dict['filter_chromosome']=set(args.filter_chromosome)
        resource_dict['filter_chromosome']=args.filter_chromosome
    else:
        resource_dict['filter_chromosome']=args.filter_chromosome

    #handle read_type
    if args.read_type_list:
        resource_dict['read_types'] = set(args.read_type_list)
    else:
        resource_dict['read_types'] = set(['RD','discordant','split'])
            
    io_append(resource_dict)
    
    return(resource_dict, resource_file_name)

def generate_run_file():
    ### creates a shell script enabling the whole pipeline to be run
    '''
    #python cvish.py -run -fa ${ref_fa} -fastq_1 ${new_fastq_dir}/${sample_name}_n01.fastq -fastq_2 ${new_fastq_dir}/${sample_name}_n02.fastq -run_name ${sample_name} -config ${sample_name}_config.tsv -run_name ${sample_name}
    '''
    if args.load_configuration_file:
        resource_dict = load_configuration_file(args.load_configuration_file)
        config_file_name = ('results/{}/output//resource.tab').format(resource_dict['results_dir'])
    else:
        resource_dict, config_file_name = generate_config_file()

    outfile_name = ('run_cvish_{}.sh').format(args.run_name)
    outfile = open(outfile_name,'w')
    print('Generating shell file ... ', outfile_name)
    
    if 'sbatch_filename' in resource_dict:
        sbatch_file_name = resource_dict['sbatch_filename']
        if sbatch_file_name:
            sbatch_file = open(sbatch_file_name)
            outfile.write(sbatch_file.read())
     
    # module file is useful for executing in a HPC environment where installed
    # components can be selectively loaded. This command appends the text of the module file
    # to the beginning of the shell script
    if 'module_filename' in resource_dict:
        module_filename = resource_dict['module_filename']
        
        if module_filename:
            module_file = open(module_filename)
            outfile.write(module_file.read())
    
    run_name = resource_dict['run_name']
    ref_fa = resource_dict['genome_fa']
    fastq_1 = resource_dict['fastq_1']
    fastq_2 = resource_dict['fastq_2']
    
    if args.filter_gff:
        resource_dict['filter_gff'] = args.filter_gff
        
    ref_line = ('name={name}\n'
                'reffa={ref_fa}\n'
                'fastq1={fastq_1}\n'
                'fastq2={fastq_2}\n'
                'conffile={config_file_name}\n').format(
                   ref_fa = ref_fa,
                   fastq_1 = fastq_1,
                   fastq_2 = fastq_2,
                   name = run_name,
                   config_file_name = config_file_name)
                   
    outfile.write(ref_line)
    
    outline = ('echo "Starting run..."\n'
               '\techo "Running -config"\n'
               'python cvish.py -config ${conffile} -fa $reffa -fastq_1 $fastq1 -fastq_2 $fastq2 -run_name $name\n')
    
    outfile.write(outline)
           
    outline = ('\techo "Running -load"\n'
               'python cvish.py -load -run_name $name\n')
    outfile.write(outline)
        
    outline = ('\techo "Running -depth"\n'
               'python cvish.py -depth -run_name $name\n')
    outfile.write(outline)
                
    outline = ('\techo "Running -peaks"\n'
               'python cvish.py -peaks -run_name $name\n')
    outfile.write(outline)
    
    if (resource_dict['filter_gff'] or 
        resource_dict['filter_bed'] or 
        resource_dict['filter_object']):
        outline = ('\techo "Running -filter"\n'
                   'python cvish.py -filter -run_name $name\n')
            
    outline = ('\techo "Running -map"\n'
               'python cvish.py -map -run_name $name\n')
    outfile.write(outline)
        
    outline = ('\techo "Running -map"\n'
               'python cvish.py -map -run_name $name\n')
    outfile.write(outline)
        
    outline = ('\techo "Running -call"\n'
               'python cvish.py -call -run_name $name\n')
    outfile.write(outline)
    
    outfile.close()
        
    if not args.no_run:
        if sbatch_file_name:
            outline = ('Beginning sbatch run: {}').format(outfile_name)
            print (outline)
            bash_command = ('sbatch {output_file}').format(output_file=outfile_name)
            subprocess.run([bash_command],stderr=subprocess.STDOUT,shell=True)
            
        else:
            outline = ('Beginning run: {}').format(outfile_name)
            print (outline)
            bash_command = ('bash {output_file}').format(output_file=outfile_name)
            subprocess.run([bash_command],stderr=subprocess.STDOUT,shell=True)

def pickle_loader(file_name, runmode='undefined'):
    if file_name:
        if runmode == 'undefined':
            print('please define runmode' + file_name)
            1/0
        
        if runmode == 'dict':
            try:
                file = open(file_name,'rb')
                object_file = pickle.load(file)
                file.close()
            except:
                print('error in loading file ' + file_name)
                1/0
        
        if runmode == 'df':
            try:
                print(file_name)
                object_file = pd.read_pickle(file_name)
            except:
                print('error in loading file ' + file_name)
                1/0
                
        return(object_file)

    else:
        return(False)

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
            
def io_overwrite(resource_dict, field, new_val):
    pre_existing_dict = {}
    
    resource_file_name = ('{}/resource.tab').format(final_output_dir)    
    resource_file = open(resource_file_name)
    
    for line in resource_file:
        line = line.strip()
        cat, val = line.split('\t')
        if cat != field:
            pre_existing_dict[cat]=val
        else:
            pre_existing_dict[cat]=new_val
    resource_file.close()

    resource_file_name = ('{}/resource.tab').format(final_output_dir)    
    resource_file = open(resource_file_name, 'w')
    
    for cat, val in pre_existing_dict.items():
        outline = ('{}\t{}\n').format(cat,val)
        resource_file.write(outline)
    resource_file.close()
    
    resource_pickle_name = ('{}/resource.p').format(final_output_dir)         
    with open(resource_pickle_name, 'wb') as file:
        pickle.dump(resource_dict, file)
        
def io_load(runmode='resource'):
    if runmode == 'resource':
        resource_pickle_name = ('{}/resource.p').format(final_output_dir)
        
        resource_dict = pickle_loader(resource_pickle_name, 'dict')
        
        return(resource_dict) 
        
    if runmode == 'fq':
        pickle_in = ("{}/uniuid_to_seq_dict.p").format(pickles_dir)
        
        uniuid_to_seq_dict = pickle_loader(pickle_in, 'dict')
                
        return(uniuid_to_seq_dict) 
    
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
    bash_command = ('samtools view -Sb {output_file_cs}.sam > {output_file_cs}_u.bam').format(output_file_cs=outfile_name)
    subprocess.run([bash_command],stderr=subprocess.STDOUT,shell=True)

    bash_command = ('samtools sort -T tmp_sort -o {output_file_cs}.bam {output_file_cs}_u.bam').format(output_file_cs=outfile_name)
    subprocess.run([bash_command],stderr=subprocess.STDOUT,shell=True)
        
    bash_command = ('samtools index {output_file_cs}.bam').format(output_file_cs=outfile_name)
    subprocess.run([bash_command],stderr=subprocess.STDOUT,shell=True)
        
    bash_command = ('rm {output_file_cs}.sam {output_file_cs}_u.bam').format(output_file_cs=outfile_name)
    subprocess.run([bash_command],stderr=subprocess.STDOUT,shell=True)
        
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
        
def make_fasta(temp_split_seq_file_name, qname_list, uniuid_to_seq_dict):
    temp_split_seq_file = open(temp_split_seq_file_name, 'w')     
    
    # preid_set = set()
    
    for qname in qname_list:
        qname = qname.split('_')[1]
                    
        if qname in uniuid_to_seq_dict:
            
            seq = uniuid_to_seq_dict[qname]
            outline = ('>{}\n{}\n').format(qname,seq)
            temp_split_seq_file.write(outline) 
        else:
            print('missing qname: ', qname)
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
        
def parse_fastq(fastq_full_name, series_number, qname_lookup, uniuid_to_seq_dict):
    
    outline = ('Parsing {}, using {}').format(fastq_full_name, series_number)
    print(outline)
    
    unfiltered_fastq_dict = {}
    
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
    line_ct = 0
    
    read_ct = 0
    hit_ct = 0
    
    for line in fq_file:
        line = line.strip()
        
        if (line[0] == '@') and (line_ct == 0):
            process = False
            read_ct += 1
            
            if (read_ct % 1000000) == 0:
                outline = ('{} reads processed').format(read_ct)
                print(outline)
            
            inline_qname = line.split('@')[1].split(' ')[0]
            
            if inline_qname in qname_lookup:
                unfiltered_fastq_dict[inline_qname] = {} 
                process = True
        
        if process:                                                           
            if line_ct==1:
                unfiltered_fastq_dict[inline_qname]['seq'] = line
                
            if line_ct==3:
                unfiltered_fastq_dict[inline_qname]['phred'] = line

        line_ct += 1
        
        if line_ct >= 4:
            line_ct = 0
            
    fq_file.close()
    
    for inline_qname in qname_lookup:
        if inline_qname in unfiltered_fastq_dict:
        
            if ('seq' not in unfiltered_fastq_dict[inline_qname]):
                print(inline_qname, 'seq failure')
                1/0
                
            if ('phred' not in unfiltered_fastq_dict[inline_qname]):
                print(inline_qname, 'phred failure')
                1/0
            #
            query_qname = inline_qname + series_number
            uniuid_to_seq_dict[query_qname] = unfiltered_fastq_dict[inline_qname]['seq']
                        
            hit_ct+=1
            
            if (hit_ct % 1000) == 0:
                 outline = ('{} discordant/split reads processed').format(hit_ct)
                 print(outline)
            
    return(uniuid_to_seq_dict)

def prep_qname(qname, hypo, outline, complete_qname_dict, qname_lookup_dict, locus_lookup_dict):
     
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
        
    return(complete_qname_dict, qname_lookup_dict, locus_lookup_dict)

def parse_brks(break_file_name, complete_qname_dict, qname_lookup_dict, locus_lookup_dict):
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
                    complete_qname_dict, qname_lookup_dict, locus_lookup_dict = prep_qname(qname, hypo, outline, complete_qname_dict, qname_lookup_dict, locus_lookup_dict)
            else:
                complete_qname_dict, qname_lookup_dict, locus_lookup_dict = prep_qname(qname, hypo, outline, complete_qname_dict, qname_lookup_dict, locus_lookup_dict)    
                
    brks.close()
    
    return(complete_qname_dict, qname_lookup_dict, locus_lookup_dict)
        
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

def make_probable_depth(sub_mean, fg_mean, std, ploidy):    
    int_depth = round((sub_mean*ploidy)/fg_mean, 3)

    z_score = (sub_mean - fg_mean)/std
       
    p_value = stats.norm.pdf(abs(z_score))*2
        
    return(int_depth, p_value)

def gff_to_feature_map(feature_map, name, chromo, start, stop):
    if chromo not in feature_map:
        feature_map[chromo]={}
        
    for nt in range(start, stop +1):
        if nt not in feature_map[chromo]:
            feature_map[chromo][nt] = set()
        
        feature_map[chromo][nt].add(name)
        
    return(feature_map)

def make_gene_by_gene_copy_number(each_sample, ref_gff_filename, fg_median, fg_mean, fg_std, mpileup_df):
    global resource_dict
    
    feature = resource_dict['ref_gff_feature']
    
    gff_feature_name = resource_dict['gff_feature_name']
    ploidy = resource_dict['expected_ploidy']
    
    cn_dict = {}
    feature_map = {}
        
    feature_depth_filename = ('{}/{}_feature_copy_number.tsv').format(final_output_dir, each_sample)
    
    feature_depth_file = open(feature_depth_filename, 'w')
    
    header = ('ID\tchromo\tstart\tstop\tsign\tsub_median\tsub_mean\tsub_std\t'
              'rel_median\trel_mean\tcum_std\tint_depth\tp_value\tdeets\n')
    feature_depth_file.write(header)
            
    ref_gff = open(ref_gff_filename)
    
    for line in ref_gff:
        if line[0]!='#':
            if line.count('\t')>7:
                #I	sgd	CDS	335	649	.	+	0	ID=CDS:YAL069W;Parent=transcript:YAL069W_mRNA;protein_id=YAL069W
                region_type = line.split('\t')[2]
                
                if feature == region_type:
                    line = line.strip()
                    chromo, _s, _r, start, stop, _dot1, sign, _dot2, deets = line.split('\t')
                    
                    start = int(start)
                    stop = int(stop)
                    
                    if gff_feature_name:
                        if gff_feature_name in deets:
                            name = deets.rsplit(gff_feature_name,1)[1].split(';')[0]
                        else:
                            name = deets
                    else:
                        name = deets
                        
                    feature_map = gff_to_feature_map(feature_map, name, chromo, start, stop)
                
                    sub_df = mpileup_df[
                        (mpileup_df["chromo"] == chromo) & (
                                   (mpileup_df["nuc"] >= int(start)) &
                                   (mpileup_df["nuc"] <= int(stop))                                   
                                   )]
                    
                    sub_median = round(sub_df["ct"].median(),3)
                    sub_mean = round(sub_df["ct"].mean(),3)
                    sub_std = round(sub_df["ct"].std(),3) 
                    
                    int_depth, p_value = make_probable_depth(sub_mean, fg_mean, (fg_std+sub_std), ploidy)
                                        
                    outline = ('{name}\t{chromo}\t{start}\t{stop}\t{sign}\t'
                               '{sub_median}\t{sub_mean}\t{sub_std}\t'
                               '{rel_median}\t{rel_mean}\t{cum_std}\t'
                               '{int_depth}\t{p_value}\t'
                               '{deets}\n').format(
                                   name = name,
                                   chromo = chromo, start = start, stop = stop, sign = sign,
                                   sub_median = sub_median, sub_mean = sub_mean, sub_std = sub_std,
                                   rel_median = round(sub_median/fg_median,3),
                                   rel_mean = round(sub_mean/fg_mean,3),
                                   cum_std = round(fg_std + sub_std, 3),
                                   int_depth = int_depth, p_value = p_value,
                                   deets = deets)
                    
                    if p_value <= 0.05:                                                  
                        print(outline)
                        
                    feature_depth_file.write(outline)
                    
                    if name not in cn_dict:
                        cn_dict[name] = {}
                        
                    locus = ('{},{},{}').format(chromo, start, stop)
                    cn_dict[name][locus] = {'cn':int_depth, 'pval':p_value}
                                                     
    feature_depth_file.close()
    
    return(cn_dict, feature_map)
    
def run_bedgraph(bam_file, each_sample):
    print('\tRunning bedtools genomecov on '+ str(bam_file) + '...')

    bashCommand = ('bedtools genomecov -bg -ibam {} > {}{}.bedgraph').format(bam_file, final_output_dir, each_sample)
    print(bashCommand)
    subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)

def run_mpileup(infile_name, each_sample, runmode, filter_by_bed, filter_bed):
    print('\tRunning mpileup genomecov on '+ str(infile_name) + '...')
     
    if runmode == 'RD':
        print('\tRunning samtools to determine sample depth...')
        if filter_by_bed:
            bashCommand = ('samtools mpileup -B -f {} {} -l {} -a -o {}/{}_{}.mpileup').format(fa_file, infile_name, filter_bed, pickles_dir, each_sample, runmode)        
        else:
            bashCommand = ('samtools mpileup -B -f {} {} -a -o {}/{}_{}.mpileup').format(fa_file, infile_name, pickles_dir, each_sample, runmode)                        
        print(bashCommand)
        subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)

    if runmode == 'discordant' or runmode == 'disco':
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
    if gff_file.rsplit('.', 1)[1] == '.gff':
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
    else:
        outline = ('Expecting .gff format for filter regions. Recieved: \n{}').format(gff_file)
        print(outline)
        
def parse_tab(tab_file, anc_QC_list):
    global resource_dict
    split_score = resource_dict['split_score']
    disco_score = resource_dict['disco_score']
    sw = resource_dict['split_weight']
    dw = resource_dict['disco_weight']
    
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

        split_weighted_score = (sw*s_ct)
        disco_weighted_score = (dw*d_ct)        
        
        if (s_ct >= 1) and (
                (split_weighted_score >= split_score) or 
                (disco_weighted_score >= disco_score)
                ):
            anc_QC_list.append(breakpoint_name)
            populate_filter_dict(chromo, start, stop, breakpoint_name)
        
        if (s_ct == 0) and (disco_weighted_score >= disco_score):
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
    
def get_prefilter_sequence(prefilter_object_name):      
    process = False

    temp_filter_object = open(prefilter_object_name)
    
    #print('prefilter_object_name', prefilter_object_name)
    seq_str = ''
    
    for line in temp_filter_object:
        
        if line[0] == '>' and not process:
            process=True
                    
        if line[0] != '>' and process:
            seq = line.strip()
            seq = seq.upper()
            #print('seq, ', seq)
            seq_str+=(seq)
                        
    return(seq_str)   

def check_coverage(source, prefilter_object_name, contig_seq_dict):      
    process = False

    temp_filter_object = open(prefilter_object_name)
    
    #print('prefilter_object_name', prefilter_object_name)
    seq_str = ''
    
    for line in temp_filter_object:
        
        if line[0] == '>' and not process:
            process=True
                    
        if line[0] != '>' and process:
            seq = line.strip()
            seq = seq.upper()
            #print('seq, ', seq)
            seq_str+=(seq)
            
    temp_filter_object.close()

    if len(seq_str)>0:
        contig_seq_dict[source]=seq_str
            
    return(contig_seq_dict)    

def long_walk(nt, long_range, ct_dict):
    miss_ct = 0
    depth = []
    
    lwgap = 10

    for s_step in long_range:
        if s_step in ct_dict:
            depth.append(float(ct_dict[s_step]))
            nt=s_step
            miss_ct = 0
            
        else:
            miss_ct += 1

        if miss_ct >= 3*lwgap:
            return(nt,depth)
            
    return(nt,depth)

def local_walk(short_range, ct_dict):
    global min_length
    
    hit_ct = 0

    for s_step in short_range:
        if s_step in ct_dict:
            hit_ct += 1
            
        if hit_ct >= (purity*min_length):
            return(True)
    else:
        return(False)
              
def build_trace(subz_df, suby_df, each_type, c_median, g_median):
    global min_length
    
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
        print('get_seq(seq_file_name)')
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

def get_otherside_dict(query_deets_dict, split_assignment, anchor):
    process = False
    best_chromo = False
    best_hf = False
    best_ht = False
    max_bit_score = 0
    
    for otherside_type in split_assignment:
        for next_anchor in split_assignment[otherside_type]:
            if next_anchor != anchor:
                chromo = query_deets_dict[next_anchor]['chromo']
                next_bit_score = query_deets_dict[next_anchor]['bit_score']
                next_hf = min(query_deets_dict[next_anchor]['hit_from'], query_deets_dict[next_anchor]['hit_to'])
                next_ht = max(query_deets_dict[next_anchor]['hit_from'], query_deets_dict[next_anchor]['hit_to'])
                    
                if next_bit_score > max_bit_score:
                    process = True
                    best_chromo = chromo
                    best_hf = next_hf
                    best_ht = next_ht
                    max_bit_score = next_bit_score
                
    return(process, best_chromo, best_hf, best_ht, max_bit_score)

def build_otherside_dict(query_deets_dict, split_assignment, anchor, otherside_type, otherside_dict):
    best_chromo = False
    best_hf = False
    best_ht = False
    max_bit_score = 0
    
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
                
            if next_bit_score > max_bit_score:
                best_chromo = chromo
                best_hf = next_hf
                best_ht = next_ht
                max_bit_score = next_bit_score
                
    return(otherside_dict, best_chromo, best_hf, best_ht, max_bit_score)

def parse_json_for_eval(json_file_name):
    '''
    Parameters
    ----------
    json_file_name : str
        json_file_name points to the json file output from the blastn results
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
    global max_eval
        
    try:
        data = json.load(open(json_file_name))
        process_data = True
    except:
        print('error in loading json file ' + json_file_name)
        process_data = False
    
    if process_data:
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
                                                                    evalue_score = float(each_hsps["evalue"])
                                                                
                                                                if evalue_score < max_eval:
                                                                    return(True)
    return(False)

def parse_json(json_file_name, seq, hypothesis_dict, anchor_contig_dict, max_eval):
    global global_hypothesis_dict
    #global global_contig_dict
    '''
    Parameters
    ----------
    json_file_name : str
        json_file_name points to the json file output from the blastn results
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
        process_data = True
    except:
        print('error in loading json file ' + json_file_name)
        process_data = False
    
    if process_data:
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
            bit_score = query_deets_dict[anchor]['bit_score']
            
            process, other_chromo, next_start, next_stop, next_bit_score = get_otherside_dict(query_deets_dict, split_assignment, anchor)
            
            if process:
                uid = len(global_hypothesis_dict)
                score = bit_score + next_bit_score
                if uid in global_hypothesis_dict:
                    print('error if uid in global_hypothesis_dict:')
                    1/0
                    
                global_hypothesis_dict[uid] = {"anchor_chromo":chromo, 
                                                       "anchor_start":hf,
                                                       "anchor_stop":ht,
                                                       "other_chromo":other_chromo,
                                                       "other_start":next_start,
                                                       "other_stop":next_stop,
                                                       "score":score,
                                                       "contig":set([seq]),
                                                       "checked":False,
                                                       "combined_with":False,
                                                       "replaced_by":False
                                                       }
            
            if chromo not in hypothesis_dict:
                hypothesis_dict[chromo] = {}
                
            #build otherside
            for nt in range(hf, ht+1):
                if nt not in hypothesis_dict[chromo]:
                    otherside_dict = {}
                else:
                    otherside_dict = hypothesis_dict[chromo][nt]
                    
                for otherside_type in split_assignment:
                    otherside_dict, other_chromo, next_start, next_stop, next_bit_score = build_otherside_dict(query_deets_dict, split_assignment, anchor, otherside_type, otherside_dict)
                    
                hypothesis_dict[chromo][nt] = otherside_dict
                 
            if chromo not in anchor_contig_dict:
                anchor_contig_dict[chromo] = {}
                
                for nt in range(hf, ht+1):
                    if nt not in anchor_contig_dict[chromo]:
                        anchor_contig_dict[chromo][nt] = set()
                        
                    anchor_contig_dict[chromo][nt].add(seq)
                                
    return(hypothesis_dict, anchor_contig_dict)

def make_anchor_regions(nt_set, gap):
    if len(nt_set) > 1:
        nt_list = list(nt_set)
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
    
    else:
        return(region_dict)

# def make_otherside_regions(chromo, region_dict, chromo_hypothesis_dict, chromo_anchor_contig_dict, gap):
#     '''
#     The key premise here is that there are two types of reads 
#     1. will be a read with unique matches on both sides (2 anchors), this is a trivial case
#     2. a read with only one unique match (1 anchors) one variable.
    
#     In the second case it can be informative to ask two questions:
#         a. what are all the variable loci the anchor maps to?
#         b. what is the relative read depth of those loci?
        
#     When a read is aligned the anchor(s) nucelotides (both read and reference) should be identified as should the read depth weighted variable.
#     The anchor and variables should then be read out. 
#     ''' 
#     global global_hypothesis_dict

#     for region_number in region_dict:
#         start = region_dict[region_number]['start']
#         stop = region_dict[region_number]['stop']
        
#         contig_set = set()
#         for nt in range(start, stop+1):
#             if nt in chromo_anchor_contig_dict:
#                 for contig in chromo_anchor_contig_dict[nt]:
#                     contig_set.add(contig)
        
#         if len(contig_set) < 1:
#             print(contig_set)
#             print(chromo, start, stop)
#             1/0
        
#         for nt in range(start, stop+1):
#             if nt in chromo_hypothesis_dict:                
#                 for otherside_chromo in chromo_hypothesis_dict[nt]:
#                     next_nt_list = []
#                     for next_nt in chromo_hypothesis_dict[nt][otherside_chromo]:
#                         next_nt_list.append(next_nt)
                        
#                     next_region_dict = make_anchor_regions(next_nt_list, gap)
                    
#                     for next_region in next_region_dict:                       
#                         next_bs_list = []
#                         next_start = next_region_dict[next_region]['start']
#                         next_stop = next_region_dict[next_region]['stop']
                        
#                         for next_nt in range(next_start, next_stop+1):
#                             #because of gaps not all nt may be present for scoring, these will be counted as a 0
#                             if next_nt in chromo_hypothesis_dict[nt][otherside_chromo]:
#                                 next_bs_list.append(chromo_hypothesis_dict[nt][otherside_chromo][next_nt])
#                             else:
#                                 next_bs_list.append(0)
                        
#                         score = round(np.median(next_bs_list))
                                                
#                         uid = len(global_hypothesis_dict)
#                         global_hypothesis_dict[uid] = {"anchor_chromo":chromo, 
#                                                                "anchor_start":start,
#                                                                "anchor_stop":stop,
#                                                                "other_chromo":otherside_chromo,
#                                                                "other_start":next_start,
#                                                                "other_stop":next_stop,
#                                                                "score":score,
#                                                                "contig":contig_set,
#                                                                "checked":False,
#                                                                "combined_with":False,
#                                                                "replaced_by":False
#                                                                }
                        
                        
                                
#     return() 

def fourway_test(uid_start, next_start, uid_stop, next_stop, gap):
    #print("uid_start, next_start, uid_stop, next_stop, gap")
    #print(uid_start, next_start, uid_stop, next_stop, gap)
    
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

def collapse_select_set(uid, next_uid):
    global resgap, global_hypothesis_dict
    '''
    This function collapses based on alignment coordinates.
    Sequences of hypotheses are stored as a set

    Parameters
    ----------
    chromoside_dict : dict
        DESCRIPTION.

    Returns
    -------
    dict with coordinates, score, and sequence set.

    '''
    sites_dict = {'anchor':'other',
                  'other':'anchor'}
    
    if (uid != next_uid):
        print('uid', uid)
        print(global_hypothesis_dict[uid])
        print()
        print('next_uid', next_uid)
        print(global_hypothesis_dict[next_uid])
        print()
        
        for site in sites_dict:
            site_uid_chromo = ('{}_chromo').format(site)
            site_uid_o_chromo = ('{}_chromo').format(sites_dict[site])
            site_uid_start = ('{}_start').format(site)
            site_uid_o_start = ('{}_start').format(sites_dict[site])
            site_uid_stop = ('{}_stop').format(site)
            site_uid_o_stop = ('{}_stop').format(sites_dict[site])
            
            for next_site in sites_dict:
                site_next_chromo = ('{}_chromo').format(next_site)
                site_next_o_chromo = ('{}_chromo').format(sites_dict[next_site])
                site_next_start = ('{}_start').format(next_site)
                site_next_o_start = ('{}_start').format(sites_dict[next_site])
                site_next_stop = ('{}_stop').format(next_site)
                site_next_o_stop = ('{}_stop').format(sites_dict[next_site])
                                        
                uid_chromo = global_hypothesis_dict[uid][site_uid_chromo]
                next_chromo = global_hypothesis_dict[next_uid][site_next_chromo]
                uid_o_chromo = global_hypothesis_dict[uid][site_uid_o_chromo]
                next_o_chromo = global_hypothesis_dict[next_uid][site_next_o_chromo]
                
                if (uid_chromo == next_chromo) and (uid_o_chromo == next_o_chromo):
                    print('pass chromo check')
                    
                    uid_start=global_hypothesis_dict[uid][site_uid_start]
                    uid_stop=global_hypothesis_dict[uid][site_uid_stop]
                    next_start=global_hypothesis_dict[next_uid][site_next_start]
                    next_stop=global_hypothesis_dict[next_uid][site_next_stop]
                    
                    print("uid_start, uid_stop", uid_start, uid_stop)
                    print("next_start, next_stop", next_start, next_stop)
                    
                    if fourway_test(uid_start, next_start, uid_stop, next_stop, resgap):
                        print('pass n fourway check')
                        
                        uid_o_start=global_hypothesis_dict[uid][site_uid_o_start]
                        uid_o_stop=global_hypothesis_dict[uid][site_uid_o_stop]
                        next_o_start=global_hypothesis_dict[next_uid][site_next_o_start]
                        next_o_stop=global_hypothesis_dict[next_uid][site_next_o_stop] 
                                                
                        print("uid_o_start, uid_o_stop", uid_o_start, uid_o_stop)
                        print("next_o_start, next_o_stop", next_o_start, next_o_stop)
                                                
                        if fourway_test(uid_o_start, next_o_start, uid_o_stop, next_o_stop, resgap):
                            print('pass o fourway check')
                                
                            return(True)
                        
    return(False)

def nested_collapse():
    #print('in nested_collapse')
    global resgap, global_hypothesis_dict, resource_dict
    
    min_confidence_score = resource_dict["min_confidence_score"]
    with_disco = resource_dict['with_disco']
    
    '''
    This function collapses based on alignment coordinates.
    Sequences of hypotheses are stored as a set

    Parameters
    ----------
    chromoside_dict : dict
        DESCRIPTION.

    Returns
    -------
    dict with coordinates, score, and sequence set.

    '''
    to_add_dict = {}
    sites_dict = {'anchor':'other',
                  'other':'anchor'}
    
    for uid in global_hypothesis_dict:
        
        if (
            (not global_hypothesis_dict[uid]['replaced_by']) or 
            (not global_hypothesis_dict[uid]['combined_with']) or 
            (not global_hypothesis_dict[uid]['checked'])
            ):
            
            if with_disco:
                uid_score = global_hypothesis_dict[uid]['score']
    
                if uid_score < min_confidence_score:
                    global_hypothesis_dict[uid]['replaced_by'] = 'filtered'
                    global_hypothesis_dict[uid]['combined_with'] = 'filtered'
                    global_hypothesis_dict[uid]['checked'] = 'filtered'
    
                    outline = ('{uid} filtered for low score: {uid_score}').format(
                        uid = uid, uid_score = uid_score)
                    print(outline)
                    
                    return(True)
            
            for next_uid in global_hypothesis_dict:
                if (
                    (not global_hypothesis_dict[next_uid]['replaced_by']) or 
                    (not global_hypothesis_dict[next_uid]['combined_with']) or 
                    (not global_hypothesis_dict[next_uid]['checked'])
                    ):
                    
                    if with_disco:
                        uid_score = global_hypothesis_dict[next_uid]['score']
    
                        if uid_score < min_confidence_score:
                            global_hypothesis_dict[next_uid]['replaced_by'] = 'filtered'
                            global_hypothesis_dict[next_uid]['combined_with'] = 'filtered'
                            global_hypothesis_dict[next_uid]['checked'] = 'filtered'
    
                            outline = ('{uid} filtered for low score: {uid_score}').format(
                                uid = next_uid, uid_score = uid_score)
                            print(outline)
                            
                            return(True)
                    
                    # outline = ("Testing {} {}").format(uid, next_uid)
                    # print(outline)
                    
                    if (uid != next_uid):
                        
                        # outline = ("Testing {} {}").format(uid, next_uid)
                        # print(outline)
                        
                        for site in sites_dict:
                            site_uid_chromo = ('{}_chromo').format(site)
                            site_uid_o_chromo = ('{}_chromo').format(sites_dict[site])
                            site_uid_start = ('{}_start').format(site)
                            site_uid_o_start = ('{}_start').format(sites_dict[site])
                            site_uid_stop = ('{}_stop').format(site)
                            site_uid_o_stop = ('{}_stop').format(sites_dict[site])
                            
                            for next_site in sites_dict:
                                site_next_chromo = ('{}_chromo').format(next_site)
                                site_next_o_chromo = ('{}_chromo').format(sites_dict[next_site])
                                site_next_start = ('{}_start').format(next_site)
                                site_next_o_start = ('{}_start').format(sites_dict[next_site])
                                site_next_stop = ('{}_stop').format(next_site)
                                site_next_o_stop = ('{}_stop').format(sites_dict[next_site])
                                                        
                                uid_chromo = global_hypothesis_dict[uid][site_uid_chromo]
                                next_chromo = global_hypothesis_dict[next_uid][site_next_chromo]
                                uid_o_chromo = global_hypothesis_dict[uid][site_uid_o_chromo]
                                next_o_chromo = global_hypothesis_dict[next_uid][site_next_o_chromo]
                                
                                if (uid_chromo == next_chromo) and (uid_o_chromo == next_o_chromo):
                                    uid_start=global_hypothesis_dict[uid][site_uid_start]
                                    uid_stop=global_hypothesis_dict[uid][site_uid_stop]
                                    next_start=global_hypothesis_dict[next_uid][site_next_start]
                                    next_stop=global_hypothesis_dict[next_uid][site_next_stop]
                                    
                                    if fourway_test(uid_start, next_start, uid_stop, next_stop, resgap):
                                        uid_o_start=global_hypothesis_dict[uid][site_uid_o_start]
                                        uid_o_stop=global_hypothesis_dict[uid][site_uid_o_stop]
                                        next_o_start=global_hypothesis_dict[next_uid][site_next_o_start]
                                        next_o_stop=global_hypothesis_dict[next_uid][site_next_o_stop] 
                                        
                                        if fourway_test(uid_o_start, next_o_start, uid_o_stop, next_o_stop, resgap):
                                            new_score = sum([global_hypothesis_dict[uid]['score'], global_hypothesis_dict[next_uid]['score']])
                                                              
                                            #TODO fix contig update to take the combined form
                                            new_contig_set = global_hypothesis_dict[uid]['contig']
                                            
                                            if not isinstance(new_contig_set, set):
                                                new_contig_set = set([new_contig_set])
                                            
                                            for contig in global_hypothesis_dict[next_uid]['contig']:
                                                new_contig_set.add(contig)
                                                print()
                                                                            
                                            to_add_dict = {"anchor_chromo":uid_chromo, 
                                                            "anchor_start":min(int(uid_start), int(next_start)),
                                                            "anchor_stop":max(int(uid_stop), int(next_stop)),
                                                            "other_chromo":uid_o_chromo,
                                                            "other_start":min(int(uid_o_start), int(next_o_start)),
                                                            "other_stop":max(int(uid_o_stop), int(next_o_stop)),
                                                            "score":new_score,
                                                            "contig":new_contig_set,
                                                            "checked":False,
                                                            "combined_with":False,
                                                            "replaced_by":False}
                                            
                                            new_uid = max(list(global_hypothesis_dict.keys()))+1
                                            global_hypothesis_dict[new_uid] = to_add_dict
                                            
                                            global_hypothesis_dict[uid]['combined_with'] = next_uid
                                            global_hypothesis_dict[uid]['replaced_by'] = new_uid
                                            global_hypothesis_dict[uid]['checked'] = True
                                            
                                            global_hypothesis_dict[next_uid]['combined_with'] = new_uid
                                            global_hypothesis_dict[next_uid]['replaced_by'] = new_uid
                                            global_hypothesis_dict[next_uid]['checked'] = True
                                            
                                            outline = ('{uid} combined with {next_uid} and replaced by {new_uid}\n\tcontig set: {new_contig_set}').format(
                                                uid = uid, next_uid = next_uid, new_uid = new_uid, new_contig_set = new_contig_set)
                                            print(outline)
                                            
                                            return(True)
                                        
                        
    return(False)

def last_pass():
    global global_hypothesis_dict
        
    temp_set = set()
    for uid in global_hypothesis_dict:
        if (
            (not global_hypothesis_dict[uid]['replaced_by']) or 
            (not global_hypothesis_dict[uid]['combined_with']) or 
            (not global_hypothesis_dict[uid]['checked'])
            ):
            temp_set.add(uid)
                            
    print('temp_set', len(temp_set))
             
    reduced = nested_collapse()               
    
    if reduced:
        return(True)
    
    #no longer 
    temp_set = set()                
    for uid in global_hypothesis_dict:
        if (
            (not global_hypothesis_dict[uid]['replaced_by']) or 
            (not global_hypothesis_dict[uid]['combined_with']) or 
            (not global_hypothesis_dict[uid]['checked'])
            ):
            print(uid, global_hypothesis_dict[uid])
            temp_set.add(uid)
                
    for uid in temp_set:
        for next_uid in temp_set:
            if uid != next_uid:
                select_set_results = collapse_select_set(uid, next_uid)
                print('uid, next_uid', uid, next_uid)
                print('select_set_results', select_set_results)
                
                if select_set_results:
                    reduced = nested_collapse()               
                    
                    if reduced:
                        return(True)
    
    return(False)

def collapse_hypotheses():
    print ("in collapse_hypotheses")
    global global_hypothesis_dict
    
    # for chromo in otherside_dict:
    #     for uid in otherside_dict[chromo]:
    #         if uid not in global_hypothesis_dict:
    #             global_hypothesis_dict[uid] = otherside_dict[chromo][uid]
    #             global_hypothesis_dict[uid]['checked'] = False
    #         else:
    #             print(uid, chromo, global_hypothesis_dict[uid])
    #             print('present', otherside_dict[chromo][uid])
    #             1/0
                
    #new_uid = original_size+1
    reduced = True
    
    while reduced:
        reduced = nested_collapse()
    
    reduced = True

    while reduced:
        reduced = last_pass()
                        
    return()
                              
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

def get_details_on_each_side(chromo, start, stop, cn_dict, feature_map):
    global resgap
    feature_set = set()
    
    if chromo in feature_map:
        for nt in range(start, stop+1):
            if nt in feature_map[chromo]:
                for name in feature_map[chromo][nt]:
                    feature_set.add(name)
                
    if len(feature_set) > 0:
        details = ('')            
        for feature_name in feature_set:
            if feature_name in cn_dict:
                for locus in cn_dict[feature_name]:
                    locus_chromo, locus_start, locus_stop = locus.split(',')
                    locus_start = int(locus_start)
                    locus_stop = int(locus_stop)
                
                    if chromo == locus_chromo:
                        if fourway_test(start, locus_start, stop, locus_stop, resgap):
                            cn = cn_dict[feature_name][locus]['cn']
                            pval = cn_dict[feature_name][locus]['pval']
                            outline = ('{feature_name};cn={cn};pval={pval}').format(
                                feature_name = feature_name,
                                cn = cn, pval = pval)
                            details+=outline
                        
        return(details)
    else:
        return('no_annotated_feature')
    
def get_details(anchor_chromo, anchor_start, anchor_stop,
                      other_chromo, other_start, other_stop,
                      cn_dict, feature_map):
    
    if not feature_map:
        return('no_ref_gff')
    
    anchor_details = get_details_on_each_side(anchor_chromo, anchor_start, anchor_stop, cn_dict, feature_map)
    other_details = get_details_on_each_side(other_chromo, other_start, other_stop, cn_dict, feature_map)
    
    details = ('anchor_details={};other_details={}').format(anchor_details, other_details)
    
    return(details)

def return_higher_cn(m_chromo_df, nuc, cnv_min_length, c_median):
    if c_median == 0:
        return(1)
    
    m_chromo_left = m_chromo_df.loc[((m_chromo_df["nuc"] >= (nuc - cnv_min_length)) & 
                                           (m_chromo_df["nuc"] <= nuc))]
    
    m_chromo_left_median = abs(np.log2(m_chromo_left["ct"].median()/c_median))
    
    m_chromo_right = m_chromo_df.loc[((m_chromo_df["nuc"] >= nuc) & 
                                           (m_chromo_df["nuc"] <= (nuc + cnv_min_length)))]
                                      
    m_chromo_right_median = abs(np.log2(m_chromo_right["ct"].median()/c_median))
    
    max_cn = max(m_chromo_left_median, m_chromo_right_median)
    
    return(max_cn)

def cn_weight(mpileup_df, chromo, start, stop, score, cnv_min_length):

    m_chromo_df = mpileup_df.loc[mpileup_df["chromo"] == chromo]
    
    c_median = m_chromo_df["ct"].median()
    
    cn_start = return_higher_cn(m_chromo_df, start, cnv_min_length, c_median)
    
    cn_stop = return_higher_cn(m_chromo_df, stop, cnv_min_length, c_median)
    
    max_cn = max(cn_start, cn_stop)

    max_cn = min(max_cn, 1e6)
    
    return(score*max_cn)

def sv_weight(score, anchor_chromo, anchor_start, anchor_stop, other_chromo, other_start, other_stop, chromo_size_dict):
    if anchor_chromo != other_chromo:        
        return(score*1e2)
    
    else:
        chromo_size = int(chromo_size_dict[anchor_chromo])
        left = min(anchor_start, other_start)
        right = max(anchor_stop, other_stop)
        
        weight = 1+(abs(left-right)/chromo_size)
                
        return(score*weight)

def summarize_hypotheses(hypothesis_dict, anchor_contig_dict, gap, resource_dict): 
    '''
    Summarize hypotheses 
    '''
    
    cnv_min_length = resource_dict['cnv_min_length']
    
    pickle_name = resource_dict['mpileup_df_RD']
    mpileup_df = pickle_loader(pickle_name, 'df')   
    
    if 'original_genome_fa' in resource_dict:
        fa_file = resource_dict['original_genome_fa']
    else:
        fa_file = resource_dict['genome_fa']
        
    bam_file = resource_dict['bam_file']
    min_confidence_score = resource_dict['min_confidence_score']
    
    if 'cn_dict' in resource_dict:
        if resource_dict['cn_dict']:
            cn_dict_p = resource_dict['cn_dict']
            cn_dict = pickle_loader(cn_dict_p, 'dict')

        else:
            cn_dict = False
    else:
        cn_dict = False
            
    if 'feature_map' in resource_dict:
        if resource_dict['feature_map']:
            feature_map_p = resource_dict['feature_map']
            feature_map = pickle_loader(feature_map_p, 'dict')
            
            # with open(feature_map_p, 'rb') as fp:
            #     feature_map = pickle.load(fp)
        else:
            feature_map = False
    else:
        feature_map = False
    
    #resource_dict['bam_file']
    
    chromo_size_p = resource_dict['chromo_size']
    chromo_size_dict = pickle_loader(chromo_size_p, 'dict')
    # with open(chromo_size_p, 'rb') as fp:
    #     chromo_size_dict = pickle.load(fp)
        
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
                  '##source=cvish\n'
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
                  '##source=cvish\n'
                  '##reference={fa_file}\n').format(
                      today=today, fa_file=fa_file)
    
    #while the vcf records are still handled by the vcf_set
    vcf_set = set()
    gff_set = set()                  

    print("Starting hypothesis reduction ...")
    good_ct = set()
    bad_ct = set()
    for uid in global_hypothesis_dict:
        if global_hypothesis_dict[uid]["contig"] == set():
            bad_ct.add(uid)
        else:
            good_ct.add(uid)
    print('good_ct', len(good_ct))
    print('bad_ct', len(bad_ct))
        
    collapse_hypotheses()
    
    outfile_name = ('global_hypothesis.csv')
    global_hypothesis_df = pd.DataFrame.from_dict(global_hypothesis_dict, orient='index')

    global_hypothesis_df.to_csv(path_or_buf=outfile_name, na_rep = np.nan)
        
    for uid in global_hypothesis_dict:
        if ((not global_hypothesis_dict[uid]['replaced_by']) or 
            (not global_hypothesis_dict[uid]['combined_with'])):
            
            anchor_chromo = global_hypothesis_dict[uid]['anchor_chromo']
            anchor_start = global_hypothesis_dict[uid]['anchor_start']
            anchor_stop = global_hypothesis_dict[uid]['anchor_stop']
            other_chromo = global_hypothesis_dict[uid]['other_chromo']
            other_start = global_hypothesis_dict[uid]['other_start']
            other_stop = global_hypothesis_dict[uid]['other_stop']
            score = global_hypothesis_dict[uid]['score']
            contig_set = global_hypothesis_dict[uid]['contig']
            
            score = cn_weight(mpileup_df, anchor_chromo, anchor_start, anchor_stop, score, cnv_min_length)
            score = cn_weight(mpileup_df, other_chromo, other_start, other_stop, score, cnv_min_length)
            score = sv_weight(score, anchor_chromo, anchor_start, anchor_stop, other_chromo, other_start, other_stop, chromo_size_dict)
            score = score/1000
                            
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
                
            if score < min_confidence_score:
                outline = ('Filtering hypothesis {uid} for low confidence score.'
                           ' Score:{score}, minimum: {min_confidence_score}').format(
                               uid = uid,
                               score = score,
                               min_confidence_score = min_confidence_score)
                print(outline)
                    
            if score > 0:
                #print('Adding breakpoint ...')

                details = get_details(anchor_chromo, anchor_start, anchor_stop,
                                      other_chromo, other_start, other_stop, 
                                      cn_dict, feature_map)
                
                filter_char = ''
                filter_deet = 'PASS'
                if score < min_confidence_score:
                    filter_char = '#'
                    filter_deet = 'FILTER'
                    
                    outline = ('Filtering hypothesis {uid} for low confidence score.'
                               ' Score:{score}, minimum: {min_confidence_score}').format(
                                   uid = uid,
                                   score = score,
                                   min_confidence_score = min_confidence_score)
                    print(outline)
            
                gff_line = ('{filter_char}{chromo}\tcvish\t{uid}_anchor_split'
                            '\t{start}\t{stop}\t.\t.\t{score}'
                            '\tnode_uid={uid};filter={filter_deet};otherside={otherside}_breeze;{details};contig={contig}\n').format(
                                filter_char = filter_char,
                                chromo = anchor_chromo,
                                uid = uid,
                                start = anchor_start,
                                stop = anchor_stop,
                                score = int(score),
                                filter_deet = filter_deet,
                                otherside = otherside,
                                details = details,
                                contig = contig)
                                
                gff_set.add(gff_line)
                
                vcf_line = ('{filter_char}{chromo}\t{pos}\t{uid}'
                            '\tN\t{alt}\t60\t{filter_deet}\tSVMETHOD=cvish;OTHERSIDE={info}'
                            '\tGT:GQ:DP\t.:.:{dp}\n').format(
                                filter_char = filter_char,
                                chromo = anchor_chromo,
                                pos = anchor_start,
                                uid = uid,
                                info = otherside,
                                alt = contig,
                                filter_deet = filter_deet,
                                dp = int(score))
                                
                vcf_set.add(vcf_line)
                
                anchorside = ('{chromo}:{start}-{stop}').format(
                    chromo = anchor_chromo,
                    start = anchor_start,
                    stop = anchor_stop)
        
                rev_gff_line = ('{filter_char}{filter_char}{chromo}\tcvish\t{uid}_breeze_split'
                            '\t{start}\t{stop}\t.\t.\t{score}'
                            '\tnode_uid={uid};filter={filter_deet};otherside={anchorside};{details};contig={contig}\n').format(
                                filter_char = filter_char,
                                chromo = other_chromo,
                                uid = uid,
                                start = other_start,
                                stop = other_stop,
                                score = int(score),
                                filter_deet = filter_deet,
                                anchorside = anchorside,
                                details = details,
                                contig = contig)
                
                gff_set.add(rev_gff_line)
                
                vcf_line = ('{filter_char}{chromo}\t{pos}\t{uid}'
                            '\tN\t{alt}\t60\t{filter_deet}\tSVMETHOD=cvish;OTHERSIDE={info}'
                            '\tGT:GQ:DP\t.:.:{dp}\n').format(
                                filter_char = filter_char,
                                chromo = other_chromo,
                                pos = other_start,
                                uid = uid,
                                info = anchorside,
                                alt = contig,
                                filter_deet = filter_deet,
                                dp = int(score))
                                
                vcf_set.add(vcf_line)
            
    return(gff_set, gff_header, vcf_set, vcf_header)

def parse_msa_file(filename, hash_to_name, msa_dict):
    '''
    Parse multiple sequence alignment file to dictionaries

    Parameters
    ----------
    filename : string
        multiple sequence alignment file name

    Returns
    -------
    msa_dict : dictionary containing name (key) to sequence

    '''
    infile = open(filename)
    seq = ''
    
    for line in infile:
        if line[0] == '>':            
            hname = (line.split('>')[1].split('_')[0].strip())

            if seq != '':
                msa_dict[hname] = seq
                
            if hname not in hash_to_name:
                hash_to_name[hname] = set()

            seq = ''
        else:
            seq += line.strip()
            
    #picks up any trailers
    msa_dict[hname] = seq
    
    infile.close()
    return(msa_dict)
    
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

def check_for_alignment(seq1, seq2):
    if seq1 == seq2:
        return(True)
    
    json_file_name = ('{temp}/temp_blast_align.json').format(temp=temp_dir)
    
    q_file_name = ('{temp}/temp_query.fa').format(temp=temp_dir)
    
    q_file = open(q_file_name, 'w')
    outline = ('>query_seq\n{}').format(seq1)
    q_file.write(outline)
    q_file.close()
    
    s_file_name = ('{temp}/temp_subject.fa').format(temp=temp_dir)
    
    s_file = open(s_file_name, 'w')
    outline = ('>subject_seq\n{}').format(seq2)
    s_file.write(outline)
    s_file.close()
                                    
    bashCommand = ('blastn -task "blastn-short" -query {q_file_name} -subject {s_file_name} -outfmt 15 -out {json_file_name}').format(
        q_file_name = q_file_name, json_file_name = json_file_name, s_file_name=s_file_name)
    
    subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)
    #print(bashCommand)
    
    results = parse_json_for_eval(json_file_name)

    return(results)

def build_subclusters(msa_dict, combined_set, check_list, prefix_dict):
    global global_contig_dict
    
    cluster = str(len(msa_dict))
        
    for name1 in msa_dict:
        name_to_cluster_set = set()
        
        if name1 not in combined_set:
            for name2 in msa_dict:
                if name2 not in combined_set:
                    check_name = ('{}_{}').format(min([name1, name2]),
                                                  max([name1, name2]))
                    
                    if (name1 != name2) and (check_name not in check_list):
                        check_list.add(check_name)

                        if check_for_alignment(msa_dict[name1], msa_dict[name2]):
                            name_to_cluster_set.add(name2)
                            
                                                            
            if len(name_to_cluster_set) > 1:
                
                combined_set.add(name1)
                
                node_name = ('{hypo}_{cluster}').format(hypo=hypo, cluster=cluster)
                cluster_seq_name = ('{temp}/{node_name}_cluster_seq.fa').format(node_name=node_name, temp=temp_dir)
                
                cluster_seq = open(cluster_seq_name, 'w')
                outline = ('>{}\n{}\n').format(name1, msa_dict[name1])
                print('outline: ', outline)
                cluster_seq.write(outline)

                for each_name2 in name_to_cluster_set:
                    combined_set.add(each_name2)
                    
                    outline = ('>{}\n{}\n').format(name2, msa_dict[name2])
                    #print('outline: ', outline)
                    cluster_seq.write(outline)
                    
                cluster_seq.close()                    
                      
                cluster_align_name = ('{temp}/{node_name}_cluster_aligned.msf').format(node_name=node_name, temp=temp_dir)
                bashCommand = ('mafft --auto --adjustdirection --globalpair --quiet {cluster_seq_name} > {cluster_align_name}').format(
                    cluster_seq_name = cluster_seq_name, 
                    cluster_align_name = cluster_align_name)
                
                subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)
                print(bashCommand)
                               
                bashCommand = ('cons -name {node_name}_{len_hnames} -plurality 1 {cluster_align_name} -outseq {cluster_seq_name}').format(
                    node_name = node_name, 
                    len_hnames = len(name_to_cluster_set), 
                    cluster_align_name = cluster_align_name, 
                    cluster_seq_name = cluster_seq_name)
                
                subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)
                print(bashCommand)
                    
                seq_str = get_prefilter_sequence(cluster_seq_name)
                
                 
                
                if len(seq_str) > 25:
                    prefix = ('{temp}/{node_name}').format(temp=temp_dir, node_name=node_name)
                    prefix_dict[cluster] = prefix
                    
                    return(True, cluster, msa_dict, combined_set, check_list, seq_str, prefix_dict)
    
    return(False, cluster, msa_dict, combined_set, check_list, '', prefix_dict)

def build_clusters(hypo, msa_dict, contig_seq_dict):
    prefix_set = set()
    prefix_dict = {}
    combined_set = set()
    check_list = set()
    
    initial_set = set(list(msa_dict.keys()))
    
    perform_subclustering = True
    while perform_subclustering:
        perform_subclustering, cluster, msa_dict, combined_set, check_list, seq, prefix_dict = build_subclusters(msa_dict, combined_set, check_list, prefix_dict)
        
        if perform_subclustering:
            msa_dict[cluster] = seq
    
    for cluster in msa_dict:
        if (cluster not in combined_set):
            if (cluster not in initial_set):
                node_name = ('{hypo}_{cluster}').format(hypo=hypo, cluster=cluster)
                #cluster_seq_name = ('{temp}/{node_name}_cluster_seq.fa').format(node_name=node_name, temp=temp_dir)
                prefix = prefix_dict[cluster]
                prefix_set.add(prefix)
                
                contig_seq_dict[node_name] = msa_dict[cluster]
            else:
                print('unclustered read ', cluster, msa_dict[cluster])
                    
    return(prefix_set, contig_seq_dict)

def correct_size(chromo, start, stop, chromo_size_dict):
    global resource_dict
    flank_nt = resource_dict['filter_flanking']
    
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
    for is_nt in (range(start, stop+1)):
        region_name = str(chromo + '_' + str(is_nt))

        if region_name in locus_to_hypo_lookup:
            #hypo_id = locus_to_hypo_lookup[region_name]
            return(locus_to_hypo_lookup[region_name])
    
    #In the event the hypo falls within a filtered region find_hypo returns a negative
    return(False)
        
def filter_regions(chromo, start, stop, region_filter_dict, ancestor_filter):
    global resource_dict
    flank_nt = resource_dict['filter_flanking']
    
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
    
def load_peaks(peaks_file_name, region_filter_dict):
    global resource_dict
    flank_nt = resource_dict['filter_flanking']
    ancestor_filter = 0
    
    peaks_file = open(peaks_file_name)
    
    peaks_dict = {}
    for line in peaks_file:
        #'{chromo}\tRead_Depth\tCNV\t{start}\t{stop}\t{rel_c}\t.\t.\tID={chromo}:{start}-{stop};rel_chromosome_RD={rel_c};rel_genome_RD={rel_g};sample={sample}\n'
        chromo = line.split('\t')[0]
        start = int(line.split('\t')[3])
        stop = int(line.split('\t')[4])
        
        process_brks = True
        if region_filter_dict:
            process_brks, ancestor_filter = filter_regions(chromo, start, stop, region_filter_dict, ancestor_filter)

        #section handles flanking and chromosome ends
        if process_brks:
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
    
def map_alignments(uid, uid_align, alignment_map, hypothesis_map, ancestor_filter, region_filter_dict):            
    process_brks = True
    chromo, start, stop, qscore, cigar, instance = uid_align
    uni_uid_instance = ('{}~{}').format(uni_uid, instance)
                                    
    if region_filter_dict:
        process_brks, ancestor_filter = filter_regions(chromo, start, stop, region_filter_dict, ancestor_filter)

    if process_brks:
        if uid not in hypothesis_map:
            hypothesis_map[uid] = [uni_uid_instance]
        else:
            hypothesis_map[uid] += [uni_uid_instance]
        
        if uni_uid_instance not in alignment_map:
            alignment_map[uni_uid_instance] = [chromo, start, stop, qscore, cigar]
        
        #small sanity check here
        else:
            p_chromo = uid_align[0]
            p_start = uid_align[1]
            p_stop = uid_align[2]
            
            if (p_chromo != chromo) or (p_start != start) or (p_stop != stop):
                print('uni_uid_instance collision', uni_uid_instance, uid_align, alignment_map[uni_uid_instance])
    
    return(alignment_map, hypothesis_map, ancestor_filter)

def refine_alignment(uni_uid, alignment, alignment_map, runmode, locus_to_hypo_lookup, refined_map):
       
    if alignment:
        q_uid = ('{}_{}').format(runmode, uni_uid)
        chromo, start, stop, qscore, cigar = alignment_map[alignment]
        # print('refine_alignment')
        # print(alignment)
        # print(alignment_map[alignment])
        breakpoint_hypothesis = find_hypo(locus_to_hypo_lookup, chromo, start, stop)

        if breakpoint_hypothesis:
            if breakpoint_hypothesis not in refined_map:
                refined_map[breakpoint_hypothesis] = set()
            
            refined_map[breakpoint_hypothesis].add(q_uid)
                                               
    return(refined_map)

""" Step Zero """
if args.load_configuration_file:
    config_file_name = args.load_configuration_file
    resource_dict = load_configuration_file(config_file_name)
        
    outline = ('Saving the following configuration parameters to run name {}\n').format(
        resource_dict['run_name'])
    print(outline)
    
# TODO making make_template_file
if args.make_template_file:
    outfile_name = args.make_template_file
    #initialize resources object and file:
    resource_dict = {'run_name':'REQUIRED: USER DEFINED',
                     'genome_fa':'REQUIRED: USER DEFINED',
                     'fastq_1':'REQUIRED: USER DEFINED',
                     'fastq_2':'REQUIRED: USER DEFINED',
                     'verbose':args.verbose,
                     'with_disco':args.with_disco,
                     'depth_region_filter':args.depth_filter_bed,
                     'ref_gff':args.ref_gff_file,
                     'ref_gff_feature':args.ref_gff_feature,
                     'gff_feature_name':args.gff_feature_name,
                     'module_filename':args.module_filename,
                     'max_gap':args.max_gap,
                     'min_purity':args.min_purity,
                     'filter_gff':args.filter_gff,
                     'filter_bed':args.filter_bed,
                     'filter_object':args.filter_object,
                     'resolution_gap':args.resolution_gap,
                     'max_eval':args.max_eval,
                     'mapq_val':args.mapq_val,
                     'filter_flanking':args.filter_flanking,
                     'split_score':args.split_score,
                     'disco_score':args.disco_score,
                     'split_weight':args.split_weight,
                     'disco_weight':args.disco_weight,
                     'sbatch_filename':args.sbatch_filename,
                     'min_confidence_score':args.min_confidence_score,
                     'cnv_min_length':args.cnv_min_length,
                     'filter_chromosome':args.filter_chromosome,
                     'no_post_run_clean_up':args.no_post_run_clean_up,
                     'expected_ploidy':args.expected_ploidy,
                     'high_sensitivity_mode':args.high_sensitivity_mode,
                     'low_sensitivity_mode':args.low_sensitivity_mode,
                     'skip_fasta_index':args.skip_fasta_index
                     }
        
    resource_file_name = io_make('template', outfile_name, resource_dict)
    
    print('Template generated at ', resource_file_name)

if args.configuration:
    """ This command sets the run configuration and experimental parameters.
    These are stored in the resource object. A resource file is created to 
    record run associated information. 
    
    This is located in 
        /idx/resources.tab
        
    """
    _resource_dict = generate_config_file()
    
    #print(resource_dict)
    
""" Step One """
if args.load_sequences:
    """ This command uses bwa mem to align the paired end (required) fastq files to the reference genome
    These are passed to samblaster which seperates out discordant and split reads
    and then samtools for conversion to bam and bam index generation.
    
    """
    #get parameters:
    resource_dict = io_load()
    #print(resource_dict)
    
    fa_file = resource_dict['genome_fa']
    ref_gff_file = resource_dict['ref_gff']
    gff_feature_name = resource_dict['gff_feature_name']
    fastq_1 = resource_dict['fastq_1']
    fastq_2 = resource_dict['fastq_2']
    
    chromo_size_dict = {}
    
    command_file_name=('{}_command_file.sh').format(s_path)
    command_file = open(command_file_name,'w')
    
    if resource_dict['module_filename']:
        module_file = open(resource_dict['module_filename'])
        command_file.write(module_file.read())
        
    if not resource_dict['skip_fasta_index']:
        outline = ('#Original reference fasta: {fa_file}\n').format(fa_file = fa_file)
        temp_fa_file = ('{temp}/temp_reference_fasta.fa').format(temp = temp_dir)
        command_file.write(outline)
        
        outline = ('cp {fa_file} {temp_fa_file}\n').format(fa_file = fa_file,
                                                           temp_fa_file = temp_fa_file)
        command_file.write(outline)
        
        resource_dict['original_genome_fa'] = fa_file
        io_overwrite(resource_dict, 'genome_fa', temp_fa_file)
                
        fa_file = temp_fa_file
        
        outline = ('bwa index {fa_file}\n').format(fa_file = fa_file)
        command_file.write(outline)
            
    ref_fa=('genome_fa={}\n').format(fa_file)
    ref_gff=('\tgenome_gff={}\n').format(ref_gff_file)
    ref_gff_feature_name=('\tgff_feature_name={}\n').format(gff_feature_name)
         
    infastq_1=('fastq_1={}\n').format(fastq_1)
    infastq_2=('fastq_2={}\n').format(fastq_2)
        
    output_filename = ('output_file={bam_dir}/{output_file}\n'
                     'mkdir -p {bam_dir}\n'
                     'echo bam_dir $PWD + {bam_dir}\n').format(
        bam_dir=bam_dir, 
        output_file=output_file)
                         
    command_file.write(ref_fa + ref_gff + ref_gff_feature_name + infastq_1 + infastq_2 + output_filename)
    outline = ('bwa mem -M -t 16 ${genome_fa} ${fastq_1} ${fastq_2} | samblaster -M -e -d ${output_file}_discordant.sam -s ${output_file}_split.sam | samtools sort -@16 -o ${output_file}.bam\n')
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
    
    header_file_name = ('{}/{}.sam').format(bam_dir, output_file)
    read_file = open(header_file_name)
    
    for line in read_file:
        if line.split('\t')[0] == '@SQ':
            chromo_name = line.split('\t')[1].split(':')[1]
            chromo_length = line.split('\t')[2].split(':')[1]
            chromo_size_dict[chromo_name] = chromo_length
            
    pickle_name = ('{}/{}_chromo_size.p').format(pickles_dir, output_file)
    pickle.dump(chromo_size_dict, open(pickle_name, 'wb'))
    
    resource_dict['standard_path'] = s_path
    resource_dict['bam_dir'] = bam_dir
    resource_dict['final_output_dir'] = final_output_dir
    resource_dict['temp_dir'] = temp_dir
    resource_dict['pickles_dir'] = pickles_dir
    
    resource_dict['chromo_size.p'] = pickle_name
     
    resource_dict['chromo_set'] = set(chromo_size_dict.keys())
            
    resource_dict['command_file_name'] = command_file_name
    
    resource_dict['bam_file'] = (bam_dir+output_file+str('.bam'))
    resource_dict['discordant_file'] = (bam_dir+output_file+str('_discordant.bam'))
    resource_dict['split_file'] = (bam_dir+output_file+str('_split.bam'))    
    
    resource_dict['sam_file'] = (bam_dir+output_file+str('.sam'))
    resource_dict['discordant_sam_file'] = (bam_dir+output_file+str('_discordant.sam'))
    resource_dict['split_sam_file'] = (bam_dir+output_file+str('_split.sam'))
    
    resource_dict['chromo_size']=(pickle_name)
    
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
            
    fa_file = resource_dict['genome_fa']
    each_sample = resource_dict['run_name']
    bam_file = resource_dict['bam_file']
    sam_file = resource_dict['sam_file']
    read_type_list = resource_dict['read_types']
    depth_filter_bed = resource_dict['depth_region_filter']
    
    filter_chromosome_set = resource_dict['filter_chromosome']        
    ref_gff = resource_dict['ref_gff']

    #Build boolean check for depth_filter_bed 
    filter_by_bed = convert_bool(depth_filter_bed)
    if not isinstance(filter_by_bed, bool):
        filter_by_bed = False
                                                               
    for each_type in read_type_list:
        if each_type == 'RD':
            infile_name = resource_dict['bam_file']
        if each_type == 'split':
            infile_name = resource_dict['split_file']
        if each_type == 'discordant' or each_type == 'disco':
            infile_name = resource_dict['discordant_file']
        
        outline = ('\nProcessing sample {}...').format(each_sample)
        print(outline)
        
        run_bedgraph(bam_file, each_sample)
        
        print('\tParsing mpileup into chromosomes...')
        run_mpileup(infile_name, each_sample, each_type, filter_by_bed, depth_filter_bed)

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
            
        if filter_chromosome_set:
            print(filter_chromosome_set)
            mpileup_df = mpileup_df[mpileup_df["chromo"].isin(filter_chromosome_set)]
            
        pickle_name = ('{}/mpileup_{}_{}.p').format(pickles_dir, each_sample, each_type)
        mpileup_df.to_pickle(pickle_name)
        
        df_name = ('mpileup_df_{}').format(each_type)
        resource_dict[df_name] = pickle_name

        fg_median = mpileup_df["ct"].median()
        fg_mean = mpileup_df["ct"].mean()
        fg_std = mpileup_df["ct"].std()
                
        zscore_read_type = ('{}_zscore_median').format(each_type)
        resource_dict[zscore_read_type] = fg_median
        zscore_read_type = ('{}_zscore_mean').format(each_type)
        resource_dict[zscore_read_type] = fg_mean
        zscore_read_type = ('{}_zscore_std').format(each_type)
        resource_dict[zscore_read_type] = fg_std
        
        if each_type == 'RD' and ref_gff:
            if (fg_median != 0) and (fg_mean != 0) and (fg_std != 0):
                cn_dict, feature_map = make_gene_by_gene_copy_number(each_sample, ref_gff, fg_median, fg_mean, fg_std, mpileup_df)
                
                pickle_name = ('{}/cn_dict_{}.p').format(pickles_dir, each_sample)
                with open(pickle_name, "wb") as fp:
                    pickle.dump(cn_dict, fp)        
                resource_dict['cn_dict'] = pickle_name
                
                pickle_name = ('{}/feature_map_{}.p').format(pickles_dir, each_sample)
                with open(pickle_name, "wb") as fp:
                    pickle.dump(feature_map, fp)       
                resource_dict['feature_map'] = pickle_name  
                
            else:
                print('Global read depths mean, median, standard deviation include zeros, skipping relative depth analysis.')
        
    io_append(resource_dict)
        
""" Step Three """
if args.peaks:
    """
    Peak detector
    1. load 'mpileup' files into dataframe pickles.    additional feature, filter regions of chromosome with large deletions
        
    """
    print('loading pickles')
    resource_dict = io_load()
    read_type_list = resource_dict['read_types']
    each_sample = resource_dict['run_name']
    chromo_size_p = pickle_loader(resource_dict['chromo_size'], 'dict')
    chromo_set = resource_dict['chromo_set']
    gap = resource_dict['max_gap']
    purity = resource_dict['min_purity']
    min_length = resource_dict['cnv_min_length']
    
    pickles_dir = resource_dict["pickles_dir"]
    print(pickles_dir)
                                
    file_name_lookup = []
    depth_dict = {}
    
    for each_type in read_type_list:
        print('Processing ', each_type)
        
        depth_gff_file_name = ('{}/{}_{}_depth.gff').format(final_output_dir, each_sample, each_type)
        depth_gff_outfile = open(depth_gff_file_name,'w')
        
        chromo_gff_file_name = ('{}/{}_{}_chromosome.gff').format(final_output_dir, each_sample, each_type)
        chromo_gff_outfile = open(chromo_gff_file_name,'w')
            
        df_name = ('mpileup_df_{}').format(each_type)
        pickle_name = resource_dict[df_name]
        mpileup_df = pickle_loader(pickle_name, 'df')
                
        zscore_read_type = ('{}_zscore_mean').format(each_type)
        g_mean = resource_dict[zscore_read_type]
        zscore_read_type = ('{}_zscore_std').format(each_type)
        g_std = resource_dict[zscore_read_type]
                        
        outline = ('Global Median:\t{}\nGlobal Std:\t{}\n').format(g_mean,g_std)
        print(outline)
         
        for every_chromo in chromo_set:
            m_chromo_df = mpileup_df.loc[mpileup_df["chromo"] == every_chromo]
            
            c_mean = m_chromo_df["ct"].mean()
            c_std = m_chromo_df["ct"].std()
            chromo_annu = c_mean+c_std
            
            start, stop = 1, len(m_chromo_df)
            outline = ('{chromo}\t{each_type}_depth\tChromosome\t{start}\t{stop}\t{c_mean}\t.\t.\tID={chromo}:{start}-{stop};Median_RD={c_mean};Std_RD={c_std};rel_genome_RD={rel_g};sample={sample}\n').format(chromo=every_chromo, each_type=each_type, start=start, stop=stop, c_mean=c_mean, c_std=c_std, rel_g=round(c_mean/float(g_mean),2), sample=each_sample)
            chromo_gff_outfile.write(outline)
        
        for every_chromo in chromo_set:
            
            m_chromo_df = mpileup_df.loc[mpileup_df["chromo"] == every_chromo]
    
            if each_type == 'RD':
                fzero = m_chromo_df.replace(0, np.NaN)
                c_median = fzero["ct"].mean()
                c_std = fzero["ct"].std()
                
            else:
                c_median = m_chromo_df["ct"].mean()
                c_std = m_chromo_df["ct"].std()              
            
            #this selects nts with more than some zscore of depth
            if each_type == 'split':
                subz_df = m_chromo_df.loc[(m_chromo_df["ct"] >= (c_median + 3*c_std))]
                suby_df = m_chromo_df.loc[m_chromo_df["ct"] >= (c_median + c_std)]
                
                build_trace(subz_df, suby_df, each_type, c_mean, g_mean)
                
            if each_type == 'disco' or each_type == 'discordant':
                subz_df = m_chromo_df.loc[(m_chromo_df["ct"] >= (c_median + 3*c_std))]
                suby_df = m_chromo_df.loc[m_chromo_df["ct"] >= (c_median + c_std)]
            
                build_trace(subz_df, suby_df, each_type, c_mean, g_mean)
            
            # this detects deletions in the DNA depth, independant of discordance
            if each_type == 'RD':
                subz_df = m_chromo_df.loc[(m_chromo_df["ct"]) == 0]
                suby_df = m_chromo_df.loc[(m_chromo_df["ct"]) == 0]
                
                build_trace(subz_df, suby_df, each_type, c_mean, g_mean)
                                
        depth_gff_outfile.close()
        chromo_gff_outfile.close()
        
        depth_pickle_name = ('{}/peaks_{}_{}.p').format(pickles_dir, each_sample, each_type)         
        with open(depth_pickle_name, 'wb') as file:
            pickle.dump(depth_dict, file)
            
    io_append(resource_dict)
        
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
    
    mapq_val = resource_dict['mapq_val']
    
    filter_object = resource_dict['filter_object']
    
    region_filter_dict = pickle_loader(filter_object, 'dict')
    
    with_disco = resource_dict['with_disco']
    
    if with_disco:
        read_type_list = set(['split', 'disco'])
    else:
        read_type_list = set(['split'])
        
    sw = resource_dict['split_weight']
    dw = resource_dict['disco_weight']
        
    # with open(chromo_size_pickle, 'rb') as fp:
    #     chromo_size_dict = pickle.load(fp)
    chromo_size_dict = pickle_loader(chromo_size_pickle, 'dict')
    
    derived_score = (resource_dict['split_zscore_median'] + resource_dict['split_zscore_std']) * sw
    resource_dict['split_derived_score'] = derived_score
    
    outline = ('Using split z-score: {score}\nSplit median: {median}\t'
               'Split standard deviation: {std}\nSplit weight {weight}').format(
        score = derived_score, 
        median = resource_dict['split_zscore_median'],
        std = resource_dict['split_zscore_std'],
        weight = sw)
                   
    print(outline)
                
    if 'discordant_zscore_median' in resource_dict:
        derived_score = (resource_dict['discordant_zscore_median'] + resource_dict['discordant_zscore_std']) * dw
        resource_dict['discordant_derived_score'] = derived_score
        
        outline = ('Using discordant z-score: {score}\nDisco median: {median}\t'
                   'Disco standard deviation: {std}\nDisco weight {weight}').format(
            score = derived_score, 
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

    alignment_map = {}
    hypothesis_map = {}
    chromo_map = {}
                      
    for runmode in read_type_list:        
        if runmode == 'discordant' or runmode == 'disco':
            print('Parsing discordant reads...')
            read_file_name = discordant_sam_file        
            peaks_file_name = ('{}/{}_discordant_depth.gff').format(final_output_dir, output_file)
        
        if runmode == 'split':
            print('Parsing split reads...')
            read_file_name = split_sam_file
            peaks_file_name = ('{}/{}_split_depth.gff').format(final_output_dir, output_file)
        
        #filter_by_peaks:
        #keys: chromo, values: sets of nucleotides
        peaks_dict = load_peaks(peaks_file_name, region_filter_dict)
            
        ct = 0
        uni_uid_set = set()
        process_list = []
        read_file = open(read_file_name)

        for line in read_file:
            process = False
            ct+=1
            
            if line[0]!='@':                
                uid = line.split('\t')[0]
                if runmode == 'split':
                    uid = uid.rsplit('_',1)[0]
                    
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
                    
                    #check if in a region of interest
                    if chromo in peaks_dict:
                        loci_set = peaks_dict[chromo]
                        
                        for nt in range(start,stop):
                            if nt in loci_set:
                                process = True
                    
                if process:
                    """ here we make the origin of the fastq file explicit, first asking is the reads are paired,
                    then asking if they are from the first ('forward') or second ('reverse') fastq files.
                    """                    
                    if np.unpackbits(np.array([flag], dtype=np.uint8))[-1]==1:
                        if np.unpackbits(np.array([flag], dtype=np.uint8))[-8]==1:
                            fastq_2_list.append(uid)
                            side = '+2'
                        else:
                            fastq_1_list.append(uid)
                            side = '+1'
                    else:
                        side = '+1'
                        
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
        
        if runmode == 'split' or runmode == 'discordant' or runmode == 'disco':     
            #this makes sure I'm adding to both sides
            uni_uid_set = uid_source_dict[runmode]
            
            for uni_uid in uni_uid_set:
                uid = uni_uid.split('+')[0]
                                               
                if (uid+'+1' not in uid_deets_dict):
                    pull_from_rd_1.add(uid)
                    fastq_1_list.append(uid)
                    
                if (uid+'+2' not in uid_deets_dict):
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
                        uni_uid = uid+'+1'
                
                if uid in pull_from_rd_2:
                    if np.unpackbits(np.array([flag], dtype=np.uint8))[-8]==1:
                        process_2 = True
                        uni_uid = uid+'+2'
                        
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
        uniuid_to_seq_dict = parse_fastq(fastq_1_filename, '+1', fastq_1_list, uniuid_to_seq_dict)

    if len(fastq_2_list) > 0:
        fastq_2_file = open(fastq_2_filename)
        uniuid_to_seq_dict = parse_fastq(fastq_2_filename, '+2', fastq_2_list, uniuid_to_seq_dict)
                      
    pickle_out = ("{}/uniuid_to_seq_dict.p").format(pickles_dir)
    pickle.dump(uniuid_to_seq_dict, open(pickle_out, 'wb'))
        
    pickle_out = ("{}/uid_source_dict.p").format(pickles_dir)
    pickle.dump(uid_source_dict, open(pickle_out, 'wb'))
    
    pickle_out = ("{}/uid_deets_dict.p").format(pickles_dir)
    pickle.dump(uid_deets_dict, open(pickle_out, 'wb'))
    
    print('Mapping regions...')
    """ This component converts the two sides of the split reads into a pair of addresses
    (ie: chromo_1, start_1, stop_1 and chromo_2, start_2, stop_2) if these are not filtered 
    by the filter regions command, they are stored in the hypothesis_map dictionary as a pair
    """
    ancestor_filter = 0
    
    for runmode in read_type_list:
        if runmode in uid_source_dict:
            uni_uid_set = uid_source_dict[runmode]
            
            if runmode == 'discordant' or runmode == 'disco':
                print('Running discordant read mode on ' + str(len(uni_uid_set)) + ' reads ...')                  
                
            if runmode == 'split':
                print('Running split read mode on ' + str(len(uni_uid_set)) + ' reads ...')
                
            for uni_uid in uni_uid_set:
                deets = uid_deets_dict[uni_uid]
                    
                if len(deets) > 1:
                    for uid_align in deets:
                        alignment_map, hypothesis_map, ancestor_filter = map_alignments(uni_uid, uid_align, alignment_map, hypothesis_map, ancestor_filter, region_filter_dict)
                else:
                    uid_align = deets[0]
                    alignment_map, hypothesis_map, ancestor_filter = map_alignments(uni_uid, uid_align, alignment_map, hypothesis_map, ancestor_filter, region_filter_dict)
        
    pickle_out = ("{}/hypothesis_map.p").format(pickles_dir)
    pickle.dump(hypothesis_map, open(pickle_out, 'wb'))
    
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
        # if each_chromo not in locus_to_hypo_lookup:
        #     locus_to_hypo_lookup[each_chromo] = {}
             
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
            
    pickle_out = ("{}/locus_to_hypo_lookup.p").format(pickles_dir)
    pickle.dump(locus_to_hypo_lookup, open(pickle_out, 'wb'))
    
    pickle_out = ("{}/hypo_to_locus_lookup.p").format(pickles_dir)
    pickle.dump(hypo_to_locus_lookup, open(pickle_out, 'wb'))

    print('Refining hypotheses ...')                 
    refined_map = {}
    uids_that_were_filtered = set()
    
    for runmode in read_type_list:
        print(runmode)
        processed_alignments = set()
        prev_ct = 0
        uni_ct = 0 
        if runmode in uid_source_dict:
            uni_uid_set = uid_source_dict[runmode]
            
                    
            for uid in uni_uid_set:
                #print(uni_uid)
                #1/0
                #uid = uni_uid.rsplit('.',1)[0]
                
                uni_ct += 1
                
                if uid in hypothesis_map:
                    alignments_set = set(hypothesis_map[uid])
                    
                    if alignments_set != processed_alignments:
                    
                        for alignment in alignments_set:
                            if alignment not in processed_alignments:
                                processed_alignments.add(alignment)
                                
                                #before = len(refined_map)
                                refined_map = refine_alignment(uid, alignment, alignment_map, runmode, locus_to_hypo_lookup, refined_map)
                                #after = len(refined_map)
                                
                                #if before != after:
                                #    print('before, after', before, after)
                    
                else:
                    uids_that_were_filtered.add(uid)
                
                pct_ct = round(100*uni_ct/len(uni_uid_set))
                if pct_ct % 10 == 0:
                    if pct_ct != prev_ct:
                        prev_ct = pct_ct
                        print(pct_ct)
                        print('refined_map size: ', len(refined_map))
                    
    print('uids_that_were_filtered: ', uids_that_were_filtered)
    
    pickle_out = ("{}/refined_map.p").format(pickles_dir)
    pickle.dump(refined_map, open(pickle_out, 'wb'))
                                
    #define outputs
    break_tab = ('{}/break_bt.tab').format(temp_dir)
    break_bed = ('{}/break_bd.bed').format(temp_dir)
    #
    break_tab_file = open(break_tab,'w')
    break_bed_file = open(break_bed,'w')
    
    refined_ct = 0
    resource_dict['break_tab'] = False
    
    for breakpoint_hypothesis, uid_list in refined_map.items():

        #print(breakpoint_hypothesis)
        process_brks_1 = True
        
        hypo_side = int(breakpoint_hypothesis)
                
        chromo_l, start_l, stop_l = hypo_to_locus_lookup[hypo_side]
        
        if filter_object:
            process_brks_1, ancestor_filter = filter_regions(chromo_l, start_l, stop_l, region_filter_dict, ancestor_filter)
        
        process_brks_1 = True
        if process_brks_1:
            
            
            if len(uid_list) > 1:
                #1/0
                split_ct = 0
                disco_ct = 0

                processed_list = []
                uid_str = ''
                
                temp_uid_dict = {}

                for each in uid_list:
                    # if each.count('_') == 1:
                    runmode, uid = each.split('_',1)
                        
                    # else:
                    #     print(each)
                    #     1/0
                    #     runmode, uid = each.split('_', 1)
                        
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
                #print('temp_uid_dict', temp_uid_dict)
                
                if len(temp_uid_dict) >= 2: 
                    s_weight = sw * split_ct
                    d_weight = dw * disco_ct
                    
                    if with_disco:
                        weight = (d_weight + s_weight)
                        score = resource_dict['discordant_derived_score']
                    else:
                        weight = s_weight
                        score = resource_dict['discordant_derived_score']
                                            
                    """generate two file sets
                    one is a text tab file the other is a bed file that can be viewed in IGV
                    """
                    if (split_ct > 1) and (weight > score):
                        refined_ct += 1
                            
                        #else:
                        outline = ('{}_both\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n').format(breakpoint_hypothesis, chromo_l, start_l, stop_l, chromo_l, start_l, stop_l, weight, uid_str)
                        break_tab_file.write(outline)
                        
                        outline = ('{}\t{}\t{}\t{}_both\t{}\t.\n').format(chromo_l, start_l, stop_l, breakpoint_hypothesis, weight)
                        break_bed_file.write(outline)
                        
                        if not resource_dict['break_tab']:
                            resource_dict['break_tab'] = break_tab

    break_tab_file.close()
    break_bed_file.close() 

    outline = ('\t{} significant hypotheses have been resolved to {} read supported breakpoints.').format(len(hypothesis_map), refined_ct)
    print(outline)

    outline = ('\tPotential breakpoints and supporting reads: \n{}').format(break_tab)
    print(outline)        

    outline = ('\tBedfile format of potential breakpoints: \n{}').format(break_bed)
    print(outline)

    io_append(resource_dict)
            
""" Step Five """
if args.make_filter:
    
    resource_dict = io_load()
    run_name = resource_dict['run_name']
    filter_gff = resource_dict['filter_gff']
    filter_bed = resource_dict['filter_bed']
    filter_object = resource_dict['filter_object']
    pickles_dir = resource_dict['pickles_dir']
    
    if filter_gff or filter_bed:        
        if filter_gff:
            outline = ("Generating filter using:\n{}").format(filter_gff)
            print(outline)
            
            gff_file = set([filter_gff])

            for each_gff_file in filter_gff:
                parse_filter_gff(each_gff_file)
                                
        if filter_bed:
            outline = ("Generating filter using:\n{}").format(filter_bed)
            print(outline)
            
            for each in set([filter_bed]):
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
                    
        if not filter_object:
            filter_object = ('{}/{}_filter_object.p').format(pickles_dir, run_name)
            resource_dict['filter_object'] = filter_object
            
        pickle_out = ("{}").format(filter_object)
        pickle.dump(filter_region_dict, open(pickle_out, 'wb'))
        
        io_append(resource_dict)
        
    else:
        outline = ("No filter files defined. Please use -filter_gff or -filter_bed to specify.")
        
""" Step Six """
if args.call_breakpoints:        
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
    resgap = resource_dict['resolution_gap']
    max_eval = resource_dict['max_eval']
    no_post_run_clean_up = resource_dict['no_post_run_clean_up']
    
    uniuid_to_seq_dict = io_load('fq')
    
    # TODO use resource
    break_tab = ('{}/break_bt.tab').format(temp_dir)
    #resource_dict['max_eval']

    """parse_brks:
    This recovers the loci from the hypothesized breakpoints identified in earlier
    steps.
    """
    if 'break_tab' not in resource_dict:
        print("No potential breakpoints found in the -map step.")
        
    if not resource_dict['break_tab']:
        print("No potential breakpoints found in the -map step.")
    
    print('Loading potential breakpoints...')
    # collects the names (complete_qname_dict) of all the reads associated with a breakpoint
    complete_qname_dict, qname_lookup_dict, locus_lookup_dict = parse_brks(break_tab, complete_qname_dict, qname_lookup_dict, locus_lookup_dict)

    # hypothesis_dict stores anchors and othersides by nt
    hypothesis_dict = {}
    # each anchor has a single contig, these contigs are collected by nt
    anchor_contig_dict = {}
    
    #to look up hypo from qname
    qname_to_hypo = {}
    
    cluster = 0
    qname_ct = 0
    
    print('Starting breakpoint testing...')
    #using the names of reads (complete_qname_dict) collects the sequences from (uniuid_to_seq_dict)
    uid_to_seq_lookup = {}
    
    for hypo, qname_list in complete_qname_dict.items():        
        qname_ct +=1
        if qname_ct % 10 == 0:
            outline = ('Breakpoint Hypotheses Tested: {} \t Number remaining {} \t Percent Complete {}\n').format(qname_ct, len(complete_qname_dict)-qname_ct, round(100*qname_ct/float(len(complete_qname_dict)),2))
            print(outline)
            print(hypo, qname_list)
            
        check_qname_set = set()
        
        for qname in qname_list:
            check_qname_set.add(qname.split('+')[0])
            
        if len(check_qname_set) > 1:
            msa_dict = {}
            hash_to_name = {}
            
            #hypo_uniuid_with_unique_seq = set()
            #make sequence file for all reads associated with a hypothesis
            temp_seq_file_name = ('{}/{}_temp_seq.fa').format(temp_dir, hypo)
            
            make_fasta(temp_seq_file_name, qname_list, uniuid_to_seq_dict)
            
            #turn that into a dictionary
            msa_dict = parse_msa_file(temp_seq_file_name, hash_to_name, msa_dict)
            
            print('locus_lookup_dict[hypo]', locus_lookup_dict[hypo])
            print(temp_seq_file_name)
            print(len(msa_dict))
            
            prefix_set, contig_seq_dict = build_clusters(hypo, msa_dict, contig_seq_dict)
                    
            for prefix in prefix_set:
                json_short_file_name = ('{}_short_blastn.json').format(prefix)
                json_long_file_name = ('{}_long_blastn.json').format(prefix)
                cluster_seq_name = ('{}_cluster_seq.fa').format(prefix)
                
                seq = get_seq(cluster_seq_name)
                
                bashCommand = ('blastn -task "blastn-short" -query {} -subject {fa_file} -outfmt 15 -out {}').format(cluster_seq_name, json_short_file_name, fa_file=fa_file)
                subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)
                hypothesis_dict, anchor_contig_dict = parse_json(json_short_file_name, seq, hypothesis_dict, anchor_contig_dict, max_eval)
          
                bashCommand = ('blastn -query {} -subject {fa_file} -outfmt 15 -out {}').format(cluster_seq_name, json_long_file_name, fa_file=fa_file)
                subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)
                hypothesis_dict, anchor_contig_dict = parse_json(json_long_file_name, seq, hypothesis_dict, anchor_contig_dict, max_eval)

    gff_set, gff_header, vcf_set, vcf_header = summarize_hypotheses(hypothesis_dict, anchor_contig_dict, resgap, resource_dict)

    hypothesis_pickle_name = ('{}/hypothesis_dict.p').format(pickles_dir)
    with open(hypothesis_pickle_name, 'wb') as file:
         pickle.dump(hypothesis_dict, file)
         
    resource_dict['hypothesis_dict.p'] = hypothesis_pickle_name
            
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
        
    outline = ('\t{} nodes realigned to {} regions.\n\n\t'
               'For a complete representation refer to:\n\t\t{} or {}').format(
                   len(complete_qname_dict), len(gff_set), gff_file_name, vcf_file_name)
    print(outline)
    
    outline = ('\tFor associated bam files, bedgraphs, and statistics check:'
               '\n\n\t{}').format(final_output_dir)
    print(outline)
    
    bam_dir = resource_dict['bam_dir']
    #temp_dir = resource_dict['temp_dir']
    pickles_dir = resource_dict['pickles_dir']
    
    # vcf_file_name = ('{}/{}_SV_CNV.vcf').format(final_output_dir, output_file)
    # vcf_file = open(vcf_file_name,'w')
    
    bashCommand = ('mv {bam_dir}/*.bam* {final_output_dir}/').format(
        bam_dir = bam_dir, final_output_dir = final_output_dir)
    print(bashCommand)       
    subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)
    
    io_append(resource_dict)
    
    if not no_post_run_clean_up:
        bam_dir = resource_dict['bam_dir']
        temp_dir = resource_dict['temp_dir']
        pickles_dir = resource_dict['pickles_dir']
        
        outline = ('\t***CLEANING INTERMEDIARY FILES***'
                   '\n\n\t{}').format(final_output_dir)
        print(outline)
        
        bashCommand = ('rm -rf {}').format(bam_dir)
        print(bashCommand)       
        subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)

        bashCommand = ('rm -rf {}').format(temp_dir)
        print(bashCommand)       
        subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)
        
        bashCommand = ('rm -rf {}').format(pickles_dir)
        print(bashCommand)       
        subprocess.run([bashCommand],stderr=subprocess.STDOUT,shell=True)        
        
        outline = ('\t***TO PREVENT THIS USE THE -no_clean FLAG ***'
                   '\n\n')
        print(outline)
    
if args.run:
    resource_dict = generate_run_file()
    
