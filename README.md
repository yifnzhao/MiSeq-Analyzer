# MiSeq-Analyzer

* The goal is to identify mutations that occur in ovarian cancer genes in mice treated with CRISPR and to visualize the data in an informative way. The genes of interest are: ***Trp53, PTEN2, PTEN3, BRCA1.***
* There are 5 mutated mice; from each mice, 5 samples were taken from left and right ovarian tumors and 3 metastasis sites (LOT, ROT, Met#1,2,3). Each sample underwent amplicon sequencing of the 4 genes of interst. There are 100 fastq files in total. 

Note: The reference reads of the 4 genes, primers and sgRNA designed to perform the CRISPR are in the ***MiSeq Primers and sgRNAs (new).docx*** document. (Note the changes in primer sequences for p53) 

### 2019-03-29
- updated MutationFinder (v5) and inputValidator (v2): added sample #26
- please download 03292019_data_summary_update.zip
- PTEN2 #26 : empty input txt file...need to investigate more about this

### 2019-03-25
- (note to self) remember to add sample # 26 (omitted previously)


### 2019-03-11
- uploaded v5_MutationFinder:
A self-contained version of MutationFinder that first generates a ±20 bp window around reference sequence (including PAM). This version is supposed to also have (almost?) correct number for the site of mutation.



Some notes about MiSeq:
(reference : https://web.uri.edu/gsc/illumina-miseq-next-generation-sequencer/ )

- performs clonal amplification, genomic DNA sequencing, and data analysis with base calling, alignment, variant calling, and reporting
- utilises a double-sided, single-lane flow cell and reagent cartridge supplied in kit form
- sequencing is performed by recordinf the synthesis of DNA strands in clusters of sample templates attached to the flow cell.
- each newly attached base liberates a fluorescent dye that is excited by diode lasers and imaged using two digita camera
- up to 96 samples may be sequenced in a single run with DNA libraries prepared with indexed or bar-coded adpatoers

- fastq only generates intermediate analysis files in FASTQ format enabling the use of third-parth tools to analyze sequencing data


### 2019-03-04
### ***lab meeting notes:***

#### next steps:
1. increase query window (20 bp upstream AND downstream of PAM)
1.  clean up and define mutations (pay attention to + and -: upstream of PAM marked as neg) (frameshift vs in-frame mutation) (missense mutation)
1.  fix BRCA1 PAM site misalignment problem
1.  filter mutations that only occurred once, or not? Is this manipulation proper statistically?
1. piechart - meaningful color codes and labels


#### other notes:
+ notion of "reverse" in comp sci vs "reverse complementary" in biology (5'-------NGG 3') 
+ frameshift vs in-frame mutation
+ proper ways to present a mutation around PAM site
+ Katie noticed  some discrepancy in the analysis result of p53 sample # 10 (T sub only detected 5 times in the software she used). We thought my result is still valid after checking with the original unvalidated reads...Dr. Yamanaka suggested it might be due to inclusion criteria...I think the software might be running some algorithm similar to cas-offinder and excluded the reads early on at the validation step. 
+ p53 insufficient reads - due to massive deletion (Katie)


#### related ideas:
+ random selection of 100 reads : is it redundant to generate more than 1 set of 100 reads per sample?
+ right now we have 100 piecharts... We probably will need to look for patterns and extract features in some way
+ could be interesting to employ machine learning approaches / principal component analysis / independent component analysis

#### (done before lab meeting today)
+ fixed some bugs in v4 (some problems with shortening query and reversion check point)
+ p53-8: there is one sequence that cannot be aligned - (search for keyword "Error" in the log)

### 2019-03-03

+ v4: PTEN3, p53 (no reverse operation)
+ v3: PTEN2, BRCA1 (need reversion)

+ not sure about BRCA1 though (it has two PAM sequences and I only used the first one (AGG, first PAM after seq reversion) because it's more convenient

+ no longer using cas-offinder for reads validation (for many reasons)
+ use v2_inputValidator instead (validates reads by matching them with right primer sequence, order shuffled, first 100 reads taken)
+ see the logs for backing-tracking mutation types, number of reads (validated / full 100?), etc
+ also note that in local folder named '<mygenename>/fastq', there are original raw fastq files as well as validated reads (all)


### 2019-02-19
+ processed PTEN2 samples 6-10
+ forgot to include wt in pie chart

### 2019-02-11
#### feedback from Katie:
+ PAM off: some of the instances might just be alignment/sequencing error? For now, just mark and keep the sequence info for further investigation
+ use samples #6-10 (i.e., p53 6-10, BRCA1 6-10, PTEN2 6-10, PTEN3 6-10)
+ for now: use sample from ***PTEN2 6-10, do not merge as they are from different locations, use R1 version (different from R2 in direction of reads) because it's more accurate***

#### a few thoughts:
+ each time we run mutationFinder, we need to:
1. read (to a dataframe) and write to a ./mutations.csv - this file maps mutation name and site to mutated seq, should contain no repeat. If new mutations are identified, write to this file.
1. write to a ./output_parsed_<samplename>.csv - this file contains occurrence counts for each mutation present in this sample (for this particular gene and location)

### 2019-02-10
#### MutationFinder
+ automatic alignment
+ allows manual curation
+ generates a csv file using pandas dataframe


### 2019-02-04
#### talked with Katie and Dr. Yamanaka
+ look for mutation directly upstream of PAM sequence (e.g. first query has 2 bp del)
+ now I am thinking that I can run cas-offinder twice per file (first filter out other lab's sequence, second reduce the window size for comparison) - no it does not work due to the same-length constraint, this program only helps filtering, have to write my own code to identify the mutations precisely
+ categorize mutation types (annotate with mutation as well)
+ how much time does it take to run the whole dataset
+ learned about compound heterozygous, homozygous, heterozygous 
+ reference: figure 1 (Katie's email)

### 2019-01-17
##### cas-offinder
+ This page provides very useful info: https://github.com/snugel/cas-offinder
+ from terminal:
```
Yifans-MBP:MiSeq yifan$ ./cas-offinder
-bash: ./cas-offinder: Permission denied
Yifans-MBP:MiSeq yifan$ chmod u+x ./cas-offinder
Yifans-MBP:MiSeq yifan$ ./cas-offinder
Cas-OFFinder v2.4 (Aug 16 2016)

Copyright (c) 2013 Jeongbin Park and Sangsu Bae
Website: http://github.com/snugel/cas-offinder

Usage: cas-offinder {input_file} {C|G|A}[device_id(s)] {output_file}
(C: using CPUs, G: using GPUs, A: using accelerators)

Example input file:
/var/chromosomes/human_hg19
NNNNNNNNNNNNNNNNNNNNNRG
GGCCGACCTGTCGCTGACGCNNN 5
CGCCAGCGTCAGCGACAGGTNNN 5
ACGGCGCCAGCGTCAGCGACNNN 5
GTCGCTGACGCTGGCGCCGTNNN 5

Available device list:
Type: CPU, ID: 0, <Intel(R) Core(TM) i7-6820HQ CPU @ 2.70GHz> on <Apple>
Type: GPU, ID: 0, <Intel(R) HD Graphics 530> on <Apple>
Type: GPU, ID: 1, <AMD Radeon Pro 455 Compute Engine> on <Apple>
Yifans-MBP:MiSeq yifan$ 
```
#### Next Steps:
+ parse the output txt file - done
+ extract info about mismatched bases - done
+ visualize - csv, can do better if necessary - done
+ streamline for multiple files and for all 4 genes
+ note to self: GPU is more robust in performing this task, use G not C
+ to organize the data files: use regular expression
```
$ mkdir ./R1fastq
$ mv *R1.fastq ./R1fastq/
```
#### Goals:
1. compare the reads - done
1. categorize the type of mutation (insertion, deletion, substitution; also, # of bp changed) for each read - done
1. visualize the result as smth like pie chart - done
1. streamline for the entire dataset - done
