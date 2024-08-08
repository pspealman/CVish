# CVish
A **C**opy-number **V**ariant and Structural Variant finder for short read sequencing that uses split and discordant reads as well read depth to identify potential breakpoints. 
(_version 2.0_)

_Important note for NYU HPC users:
 You can load all necessary modules via:_
 ```
 source demo/module_load.sh
 ```

## Quick Start:
CVish can be run using a single Illumina paired-end sequencing library (.fastq) along with a reference genome (.fa), it will produce variant calls in both a GFF3 (.gff) and VCFv4.2 (.vcf) format by just using the command line as.
``` 
python cvish.py -run -fa <reference_genome.fa> -fastq_1 <n01_fastq.gz> -fastq_2 <n02.fastq.gz> -run_name <unique_sample_name>
# output is stored in results/<unique_sample_name>/output/<unique_sample_name>_SV_CNV.gff
```

## Example run:
 ### Install and run preparation:
 ```
 git clone https://github.com/pspealman/CVish.git
 cd CVish
 wget https://ftp.ensembl.org/pub/release-112/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz -O S288C_Ensembl.fa.gz
 ```
 ### Run whole analysis with __-run__ command
 This generates the required split and discordant sam files from fastq files. 
 ```
  python cvish.py -run -fa S288C_Ensembl.fa.gz -fastq_1 demo/n01_ancestor.fastq.gz -fastq_2 demo/n02_ancestor.fastq.gz -run_name ancestor_demo
  head results/ancestor_demo/output/ancestor_demo_SV_CNV.gff
 ```
 ### OPTIONAL Use first results as filter for second sample 
 In instances where an ancestor and evolved strain are sequenced the results of one run can be used to filter .
 ```
  python cvish.py -run -fa S288C_Ensembl.fa.gz -fastq_1 demo/n01_ancestor.fastq.gz -fastq_2 demo/n02_ancestor.fastq.gz -run_name ancestor_demo
  python cvish.py -run -fa S288C_Ensembl.fa.gz -fastq_1 demo/n01_evolved.fastq.gz -fastq_2 demo/n02_evolved.fastq.gz -exclude results/ancestor_demo/output/ancestor_demo_SV_CNV.gff -run_name evolved_demo
 ```

## Basic Analysis Tutorial:
 ### Step One: Preparation
 Before the analysis we need to download a reference genome file - this can be done using the wget command
 ```
   wget https://ftp.ensembl.org/pub/release-112/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz -O S288C_Ensembl.fa.gz
   gunzip S288C_Ensembl.fa.gz
   head S288C_Ensembl.fa
 ```
 Note that _Saccharomyces cereivisiae_'s genome in Ensembl uses the 'Roman Numeral' (ie. 'I', 'II', 'XI') chromosome naming convention with the mitochondrial genome named 'Mito'.
 * We can download the matching GFF from the same source as well.
 ```
   wget https://ftp.ensembl.org/pub/release-112/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.112.gff3.gz -O S288C_Ensembl.gff.gz
   gunzip S288C_Ensembl.gff.gz
   head -n 30 S288C_Ensembl.gff
 ```
 * We can verify that the gff is also using the 'Roman Numeral' (ie. 'I', 'II', 'XI') chromosome naming convention.
 * It is also worth noting that 
 * Since we want to exclude regions like the rDNA locus from analysis we want to make sure that we select the appropriate 'excluded regions' file.
 ```
   head filter_files/S288C_Ensembl_exclude_regions.bed
 ```
 * We find the BED file uses the same naming convention, so we can proceed with **Step Two**
   
 ### Step Two: Run Command on Ancestor Strain
 * We're first going to analyze the sequencing data for the Ancestor, using the files 
 ```
   python cvish.py -run -fa S288C_Ensembl.fa.gz -fastq_1 <read 1 of 2 PE> -fastq_2 <read 2 of 2 PE> -run_name <name of run>
 ```
 * Example run on ancestor strain:
 ```
   
   python cvish.py -run -fa demo/S288C_R64_demo.fa.gz -fastq_1 demo/n01_ancestor.fastq.gz -fastq_2 demo/n02_ancestor.fastq.gz -exclude saccharomyces_cerevisiae_chromosome_Ensembl_rDNA_exclude.bed demo/ -run_name demo_anc
   nano results/demo_anc/output/demo_anc_SV_CNV.gff
 ```
  * Example run on evolved strain: 
 ```
  python cvish.py -run -fa demo/demo.fna -fastq_1 demo/n01_evolved.fastq.gz -fastq_2 demo/n02_evolved.fastq.gz -config  -filter_gff results/demo_anc/output/demo_anc_SV_CNV.gff -run_name demo_evo
  nano results/demo_evo/output/demo_evo_SV_CNV.gff
 ```
 ### Important notes on select configuration file parameters
 Many parameters can and should be modified for optimal performance, these will be enumerated in full below in the **Configuration section** but several important components are discussed here:
  #### Make an example configuration file
  You can always make a template configuration file by using the ```make_template_file``` command:
```
  python cvish.py -template example.tsv
```
  #### Ploidy
  CVish runs as haploid by defualt, this can be adjusted by setting the ```expected_ploidy``` value to some non-zero integer.
  
  #### Analysis parameters
  Several parameters can have a profound effect on the analysis. ```min_confidence_score``` sets the lower limit for potential breakpoints to be reported. Breakpoints scoring beneath these are output as **FILTERED**. To limit analysis to specifc chromosomes use the ```filter_chromosome``` parameter followed by a space or comma seperated list of permitted chromosomes. 
  Finally, breakpoint predictions can be generated by both split-reads and discordant reads. Utilizing both is more informative, and increases specificity and accuracy, but can increase run times 10-20% above what is seen if split-reads alone are used. To force only split reads to be used set the ```with_disco``` parameter to false.
  Preloaded alternative defaults are available in the configuration file, ```high_sensitivity_mode``` has lower thresholds for inclusion and analysis and should detect weaker signals, at the risk of increasing FDR. While ```low_sensitivity_mode``` has even higher thresholds than the default and should exclude more, and have lower FDR.  
 
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

## Using reference GFF files
   One optional feature of CVish is that it can calculate expected copy number over any region. 
   Using this option requires supplying a GFF that has the same coordinate system as the reference genome FASTA file (so don't use the Ensembl genome GFF with the NCBI genome FASTA, for example). This can be set in the config file at the ```ref_gff``` line as such: 
```ref_gff /demo/demo.gff```

   If we open this GFF file we'll find the following lines:
```
NC_001140.6	RefSeq	region	1	5505	.	-	.	ID=id3045;Dbxref=SGD:S000028891;Note=TEL08L%3B Telomeric region on the left arm of Chromosome VIII%3B composed of an X element core sequence%2C an X element combinatorial repeat%2C a short Y' element%2C and a short terminal stretch of telomeric repeats;gbkey=telomere
NC_001140.6	RefSeq	gene	445	3311	.	-	.	ID=gene2517;Dbxref=GeneID:856335;Name=YHL050C;gbkey=Gene;locus_tag=YHL050C
NC_001140.6	RefSeq	mRNA	445	3311	.	-	.	ID=rna2505;Parent=gene2517;Dbxref=Genbank:NM_001179130.1;Name=NM_001179130.1;gbkey=mRNA;product=hypothetical protein;transcript_id=NM_001179130.1
NC_001140.6	RefSeq	exon	2671	3311	.	-	.	ID=id3046;Parent=rna2505;Dbxref=Genbank:NM_001179130.1;gbkey=mRNA;product=hypothetical protein;transcript_id=NM_001179130.1
NC_001140.6	RefSeq	exon	445	1897	.	-	.	ID=id3047;Parent=rna2505;Dbxref=Genbank:NM_001179130.1;gbkey=mRNA;product=hypothetical protein;transcript_id=NM_001179130.1
NC_001140.6	RefSeq	CDS	2671	3311	.	-	0	ID=cds2355;Parent=rna2505;Dbxref=SGD:S000001042,Genbank:NP_011813.1;Name=NP_011813.1;Note=hypothetical protein%3B potential Cdc28p substrate;gbkey=CDS;product=hypothetical protein;protein_id=NP_011813.1
NC_001140.6	RefSeq	CDS	445	1897	.	-	1	ID=cds2355;Parent=rna2505;Dbxref=SGD:S000001042,Genbank:NP_011813.1;Name=NP_011813.1;Note=hypothetical protein%3B potential Cdc28p substrate;gbkey=CDS;product=hypothetical protein;protein_id=NP_011813.1
``` 
   Note that the **Feature** column (3rd column), describes the GFF feature and has options like *region, gene, mRNA, exon, CDS* etc. Each record also has an **Attribute** column (9th column) that contains descriptions of the feature. Here we're interested in the "Systematic Name" which is stored in the "locus_tag" attribute of the "gene" feature. To use this we'll need to define these in the config file by editing the ```ref_gff_feature``` and ```gff_feature_name``` parameters:
```
ref_gff_feature	gene
gff_feature_name	locus_tag=
```
   Finally, if you are using a non-haploid organism you will need to specif the expected ploidy number by editing the ```expected_ploidy``` parameter:
```
expected_ploidy 2
```

## Output Analysis 
After the run has completed several resulting files should be analyzed to help determine likely breakpoints.

 ### Output files
 Many files intended to help in the analysis of predicted breakpoints are in the ```/results/<sample_name>/output/``` directory. Files that may be of interest are:

 ```
 <sample_name>_SV_CNV.gff - Resolved high confidence breakpoint, with scores and breakpoint sequences, reported as features.
 <sample_name>_SV_CNV.vcf - Resolved high confidence breakpoint, with scores and breakpoint sequences, reported as variants.
 <sample_name>_feature_copy_number.tsv - [optional] - if a reference gff file is provided during the run, the selected features will have their read depth statistics reported here. 
 <sample_name>.bam - Complete, aligned reads bam
 <sample_name>.bedGraph - Bedgraph of complete, aligned bam
 <sample_name>_discordant.bam - Discordant reads bam
 <sample_name>_split.bam - Split reads bam
 ```

 #### How to make sense of the results in the _SV_CNV.gff_ file
 The _SV_CNV.gff_ file follows a simple convention for highlighting breakpoints. Each breakpoint has at least two components: an anchor and a breezepoint. The anchor, by definition, contains a uniquely mapping sequence. Conversely, the breezepoint can map to numerous places with no restriction, this allows for potential breezepoints to map to repeat or low-complexity regions, such as gene homologs and transposon elements. When possible the assemdbled contig that generated the breakpoint is included.

 #### Example _SV_CNV.gff_
 ```
 [chromosome]    [source]    [unique_id]   [start]   [stop]    [dot] [dot]   [score]   [details]
VI	cvish	16_anchor_split	55630	55685	.	.	25	node_uid=16;filter=PASS;otherside=XI:513836-513945_breeze;anchor_details=name=YFL038C_cn=1.56_pval=0.7328444325612059;other_details=no_annotated_feature;contig=ATGATATTTCAACATCATGTTTTGCTTAGTAGACTCTTGCGGGCGTTCCATCCGTGTGAAATACATCATTTACACCTCGCTCTGGGTCAAGTAATCAAAAAATACCTCGTGTGGCTGCAAGAGATTGATCGGTATGCAACCTCAACAGTGTtgaagctattggtactgtctcttatactcatctgacgctgccg;
XI	cvish	16_breeze_split	513836	513945	.	.	25	node_uid=16_anchor;filter=PASS;otherside=VI:55630-55685;anchor_details=name=YFL038C_cn=1.56_pval=0.7328444325612059;other_details=no_annotated_feature;contig=ATGATATTTCAACATCATGTTTTGCTTAGTAGACTCTTGCGGGCGTTCCATCCGTGTGAAATACATCATTTACACCTCGCTCTGGGTCAAGTAATCAAAAAATACCTCGTGTGGCTGCAAGAGATTGATCGGTATGCAACCTCAACAGTGTtgaagctattggtactgtctcttatactcatctgacgctgccg;
 ```
 In this example the proposed breakpoint "16" spans from chrVI:55630-55685 (anchor) to chrXI:513836-513945 (breezepoint). It has a score of 25. The contig generated by the reads associated with this proposed breakpoint is ```ATGATATTTCAACATCATGTTTTGCTTAGTAGACTCTTGCGGGCGTTCCATCCGTGTGAAATACATCATTTACACCTCGCTCTGGGTCAAGTAATCAAAAAATACCTCGTGTGGCTGCAAGAGATTGATCGGTATGCAACCTCAACAGTGTtgaagctattggtactgtctcttatactcatctgacgctgccg```.

 ## Note on performance
 CVish performs best with breakpoints spanning unique sequences and at sufficient depth (~30x depth). For lower depths of sequencing or breakpoints that occur in low-complexity or non-unique sequences this performance will suffer. 

# Example analysis:
## Data 
The following example uses data from [DGY1657](https://doi.org/10.1128/mra.00729). 

## Run
The [DGY1657_config.tsv](https://github.com/pspealman/CVish/blob/master/demo/demo_config.tsv) is available in the demo folder. This was ran with the following command:
```
python cvish.py -run -config demo/DGY1657_config.tsv -run_name Anc
```

## Interpretation
### Predictions
This resulted in the following _CNV_SV.gff_ file:
```
###fileDate=20_05_2023
##source=cvish
##reference=/scratch/work/cgsb/genomes/Public/Fungi/Saccharomyces_cerevisiae/Ensembl/R64-1-1/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa
VI	cvish	16_anchor_split	55630	55685	.	.	25	node_uid=16;filter=PASS;otherside=XI:513836-513945_breeze;anchor_details=name=YFL038C_cn=1.56_pval=0.7328444325612059;other_details=no_annotated_feature
X	cvish	9_breeze_split	136913	136939	.	.	3	node_uid=9_anchor;filter=PASS;otherside=VI:226057-226209;anchor_details=name=YFR034W-A_cn=1.148_pval=0.7822299621511224name=YFR035C_cn=0.997_pval=0.7978788778605969;other_details=no_annotated_feature;contig=gtGTATGTGTATGTCTGTGTGCAAGTATTTCCTATGCTGCAAGTGCGATTTTCTCGTTTTCTAtTTTTTTTTTTTTTTTGCCTCGCCTAATATGTGGTAGGCGAAAGGCTGACCCGGCCGCTCGCACGGAAATATTTGGCAAATGAGTCTTGAC;
VI	cvish	9_anchor_split	226057	226209	.	.	3	node_uid=9;filter=PASS;otherside=X:136913-136939_breeze;anchor_details=name=YFR034W-A_cn=1.148_pval=0.7822299621511224name=YFR035C_cn=0.997_pval=0.7978788778605969;other_details=no_annotated_feature;contig=gtGTATGTGTATGTCTGTGTGCAAGTATTTCCTATGCTGCAAGTGCGATTTTCTCGTTTTCTAtTTTTTTTTTTTTTTTGCCTCGCCTAATATGTGGTAGGCGAAAGGCTGACCCGGCCGCTCGCACGGAAATATTTGGCAAATGAGTCTTGAC;
XI	cvish	16_breeze_split	513836	513945	.	.	25	node_uid=16_anchor;filter=PASS;otherside=VI:55630-55685;anchor_details=name=YFL038C_cn=1.56_pval=0.7328444325612059;other_details=no_annotated_feature;contig=ATGATATTTCAACATCATGTTTTGCTTAGTAGACTCTTGCGGGCGTTCCATCCGTGTGAAATACATCATTTACACCTCGCTCTGGGTCAAGTAATCAAAAAATACCTCGTGTGGCTGCAAGAGATTGATCGGTATGCAACCTCAACAGTGTtgaagctattggtactgtctcttatactcatctgacgctgccg;
XI	cvish	0_breeze_split	524082	524116	.	.	3	node_uid=0_anchor;filter=PASS;otherside=XI:523865-524009;anchor_details=name=YKR045C_cn=1.063_pval=0.7964352276922085;other_details=no_annotated_feature;contig=GAGACAtacgttgctgtttgagcactatATAAAAGTGAATATACATGGTGCCAAGGGAAAGTCCTAGCATGTCAAACAGTCATCATACATCCCAGGGAAGAAGAAACAAACTTTCAGTGTGGGTTAAAAAAATAATAAACACCACTACCACAACTAACGCATgaCGGTCTCGAGC;
XI	cvish	0_anchor_split	523865	524009	.	.	3	node_uid=0;filter=PASS;otherside=XI:524082-524116_breeze;anchor_details=name=YKR045C_cn=1.063_pval=0.7964352276922085;other_details=no_annotated_feature;contig=GAGACAtacgttgctgtttgagcactatATAAAAGTGAATATACATGGTGCCAAGGGAAAGTCCTAGCATGTCAAACAGTCATCATACATCCCAGGGAAGAAGAAACAAACTTTCAGTGTGGGTTAAAAAAATAATAAACACCACTACCACAACTAACGCATgaCGGTCTCGAGC;
```
Here, we see three hypotheses: "16","9", and "0" with scores of 25, 3, and 3, respectively. **Hypothesis 16** maps from chrVI:55630-55685 on one side to chrXI:513836-513945 on the other - this agrees with what we know about this strain, which is that we used the ACT1 promotor to drive a reporter inserted upstream of the GAP1 locus. The other two, 9 and 0, are low scoring and probably not accurate predictions. To get a better sense we can also load the visualization files located in the output directory.  

![Figure with representative RD, discordant, and split-reads with breakpoint predictions](https://github.com/pspealman/CVish/blob/master/demo/DGY1657_cvish.png)

This figure shows reads distributions over the _GAP1_ neighborhood on chromosome XI in S. cerevisiae. Illumina short reads (light blue), ONT long reads (dark blue) - note the very large variation in read depth seen in the short reads is not observed for the long reads. The variation in read depth in the short reads makes the accurate calling of CNVs by read depth alone difficult. 

The figure also shows the discordant (orange) and split reads (red) extracted from the Illumina short reads. Here see some shallow distribution of discordant reads across the whole region. Conversely, we find a distinct peak of split reads at three sites. From the right to the left, the first and largest peak (14 reads deep) corresponds with **Hypothesis 16**, it also has 5 discordant reads that support this hypothesis as well. And it is a known true positive as well, as it is the known insertion site of the ACT1 promoter containing construct. The next peak has no corresponding Hypothesis but does correspond to the homopolymer region spanning chrXI:518102-518142. Because this region has such low sequence complexity these reads are filtered from our analysis. Finally, there is the third peak, corresponding to both sides of **Hypothesis 0**, this is a low scoring hypothesis and indeed it is only supported by two split reads and no discordant reads.

# Detailed manual:
 ## Requirements
 CVish intends to be a requirement lite tool and as such is contained entirely in single python script. It requires the following programs to be installed in the environment:

  * python	(tested on version: 3.8.6)
  * bwa		(tested on version: 0.7.17) **NB** bwa is required and is not compatible with bwa-mem2 
  * samtools	(tested on version: 1.14)
  * bedtools	(tested on version: 2.29.2)
  * blast	(tested on version: 2.11.0)
  * samblaster	(tested on version: 0.1.26)
  * mafft	(tested on version: 7.475)
  * emboss	(tested on version: 6.6.0)

<!--
## Built in demo and test 
 For demonstration of command syntax use:
  ```
  python cvish.py -demo
  ```
 To run an install test using defaults, use:
  ```
  python cvish.py -test
  ```
  -->
