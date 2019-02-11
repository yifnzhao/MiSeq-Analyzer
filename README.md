# MiSeq

* The goal is to identify mutations that occur in ovarian cancer genes in mice treated with CRISPR and to visualize the data in an informative way. The genes of interest are: ***Trp53, PTEN2, PTEN3, BRCA1.***
* There are 5 mutated mice; from each mice, 5 samples were taken from left and right ovarian tumors and 3 metastasis sites (LOT, ROT, Met#1,2,3). Each sample underwent amplicon sequencing of the 4 genes of interst. There are 100 reads in total. 

Note: The reference reads of the 4 genes, primers and sgRNA designed to perform the CRISPR are in the ***MiSeq Primers and sgRNAs.docx*** document. 


## 2019-02-10
#### MutationFinder
+ automatic alignment
+ allows manual curation
+ generates a csv file using pandas dataframe


## 2019-02-04
#### talked with Katie and Dr. Yamanaka
+ look for mutation directly upstream of PAM sequence (e.g. first query has 2 bp del)
+ now I am thinking that I can run cas-offinder twice per file (first filter out other lab's sequence, second reduce the window size for comparison)
+ categorize mutation types (annotate with mutation as well)
+ how much time does it take to run the whole dataset
+ learned about compound heterozygous, homozygous, heterozygous 
+ Read figure 1 (Katie's email)

## 2019-01-17
#### use cas-offinder
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
### my next step
+ parse the output txt file - done
+ extract info about mismatched bases - done
+ visualize - csv, can do better if necessary
+ streamline for multiple files and for all 4 genes

+ note to self: GPU is more robust in performing this task, use G not C
### goal:
1. compare the reads - done
1. categorize the type of mutation (insertion, deletion, substitution; also, # of bp changed) for each read - done
1. visualize the result as smth like pie chart
1. streamline for the entire dataset
