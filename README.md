# CVish
A **C**opy-number **V**ariant and Structural Variant finder for short read sequencing that uses split and discordant reads as well read depth to identify potential breakpoints. 

_Important note for NYU HPC users:
 You can load all necessary modules via:_
 ```
 source demo/module_load.sh
 ```

## Quick Start:
 ### Install:
 ```
 git clone https://github.com/pspealman/CVish.git
 cd CVish
 ```
 ### Run whole analysis with __-run__ command
 * This generates the required split and discordant sam files from fastq files. Should be ran once for the ancestor strain and once for each evolved strain.
 ```
   python cvish.py -run -fa <reference_genome_fasta_file> -fastq_1 <read 1 of 2 PE> -fastq_2 <read 2 of 2 PE> -config <path to config file> -run_name <name of run>
 ```
 * Example run on ancestor strain:
 ```
   python cvish.py -run -fa demo/demo.fna -fastq_1 demo/n01_ancestor.fastq.gz -fastq_2 demo/n02_ancestor.fastq.gz -config demo/demo_config.tsv -run_name demo_anc
   nano results/demo_anc/output/demo_anc_SV_CNV.gff
 ```
  * Example run on evolved strain: 
 ```
  python cvish.py -run -fa demo/demo.fna -fastq_1 demo/n01_evolved.fastq.gz -fastq_2 demo/n02_evolved.fastq.gz -config demo/demo_evo_config.tsv -filter_gff results/demo_anc/output/demo_anc_SV_CNV.gff -run_name demo_evo
  nano results/demo_evo/output/demo_evo_SV_CNV.gff
 ```
 ### Important notes on select configuration file parameters
 Many parameters can and should be modified for optimal performance, these are enumerated in full below in the **Configuration section** but several important components are discussed here:
  #### Analysis parameters
  Several parameters can have a profound effect on the analysis. ```min_confidence_score``` sets the lower limit for potential breakpoints to be reported. Breakpoints scoring beneath these are out
 
  #### _-run_ command optional parameters
  The run command can be used incongunction with modules and SLURM's sbatch job scheduler using the ```module_filename``` and ```sbatch_filename``` parameters. To generate a shell executable, but not run it, use the ```-no_run``` argument from the command line. 
  
  #### Filters for improved performance:
  Using filters can improve the speed and performance of CVish. These filters include:
   * Read depth filters ```depth_filter_bed``` which limit read depth calculations to regions of interest and avoid read depth measures over regions such as the rDNA locus. An example of this is located in: ```filter_files/saccharomyces_cerevisiae_chromosome_NCBI_rDNA_filter.bed``` 
   * Region filters ```filter_bed``` and ```filter_gff``` prevent breakpoint predicitions from occuring within the defined regions. These can be used to filter problematic regions, such as low complexity regions as in this telomere example:
```filter_files/saccharomyces_cerevisiae_chromosome_NCBI_filter_telomeres.bed```
Or regions that are repetitive throughout the genome, such as in this transposon example:
```filter_files/saccharomyces_cerevisiae_chromosome_transposon_GAP1.gff```
   * Ancestor filters ```filter_bed``` and ```filter_gff``` are also used to prevent breakpoint predictions from occurring in regions already identified in the ancestor. An example of this can be seen in the __-run__ section above:
```results/demo_anc/output/demo_anc_SV_CNV.gff```  

## Output Analysis 
After the run has completed several resulting files should be analyzed to help determine likely breakpoints.

All files are located in the ```results/<sample_name>/``` directory.

 ### Bams 
 Discordant and split reads are subset in their own bam files, viewable in the ```/results/<sample_name>/bam/``` directory. Files that may be of interest are:
 ```
 <sample_name>.bam - Complete, aligned reads bam
 <sample_name>.bedGraph - bedgraph of complete, aligned bam
 <sample_name>_discordant.bam - Discordant reads bam
 <sample_name>_split.bam - Split reads bam
 ```

 ### Output files
 Many files intended to help in the analysis of predicted breakpoints are in the ```/results/<sample_name>/output/``` directory. Files that may be of interest are:

 ```
 <sample_name>_SV_CNV.gff - Resolved high confidence breakpoint, with scores and breakpoint sequences, reported as features.
 <sample_name>_SV_CNV.vcf - Resolved high confidence breakpoint, with scores and breakpoint sequences, reported as variants.
 <sample_name>_feature_copy_number.tsv - [optional] - if a reference gff file is provided during the run, the selected features will have their read depth statistics reported here. 
 ```

 #### How to make sense of the results in the _SV_CNV.gff_ file
 The _SV_CNV.gff_ file follows a simple convention for highlighting breakpoints. Each breakpoint has at least two components: an anchor and a breezepoint. The anchor, by definition, contains a uniquely mapping sequence. Conversely, the breezepoint can map to numerous places with no restriction, this allows for potential breezepoints to map to repeat or low-complexity regions, such as gene homologs and transposon elements.

 When possible the assemdbled contig that generated the breakpoint is included.

 #### Example _SV_CNV.gff_

 ```
 [chromosome]    [source]    [unique_id]   [start]   [stop]    [dot] [dot]   [score]   [details]
 I   cvish   7458_anchor_split	39821    39905   .   .   7    node_uid=7458;otherside=NC_001147.6:1091237-1091291_breeze;contig=CACACACACCACACCCACACACCCACACACCACACCCACACACTCTCTCACATCTACCTCTACTCTCGCTGTCAT
 IV   cvish   7458_breeze_split   1091237   1091291   .   .   7    node_uid=7458_anchor;otherside=NC_001144.5:30-111;contig=CACACACACCACACCCACACACCCACACACCACACCCACACACTCTCTCACATCTACCTCTACTCTCGCTGTCAT
 ```
 In this example the proposed breakpoint "7458" spans from chrI:30-111 (anchor) to chrIV:1091237-1091291 (breezepoint). It has a score of 7 and would not be reported using the default _min_score_ of 10.

 ## Note on performance
 CVish performs best with breakpoints spanning unique sequences and at sufficient depth (~30x depth). For lower depths of sequencing or breakpoints that occur in low-complexity or non-unique sequences this performance will suffer. 

# Detailed manual:
 ## Requirements
 CVish intends to be a requirement lite tool and as such is contained entirely in single python script. It requires the following programs to be installed in the environment:

  * python	(tested on version: 3.8.6)
  * bwa		(tested on version: 0.7.17)
  * samtools	(tested on version: 1.14)
  * bedtools	(tested on version: 2.29.2)
  * blast	(tested on version: 2.11.0)
  * samblaster	(tested on version: 0.1.26)
  * mafft	(tested on version: 7.475)
  * emboss	(tested on version: 6.6.0)

## Built in demo and test 
 For demonstration use:
  ```
  python cvish.py -demo
  ```
 To run an install test using defaults, use:
  ```
  python cvish.py -test
  ```
