# mice_seq
* This repo contains my codes for mice sequence alignment. 
* The goal is to identify mutations that occur in ovarian cancer genes in mice treated with CRISPR and to visualize the data in an informative way. The genes of interest are: ***Trp53, PTEN2, PTEN3, BRCA1.***
* There are 5 mutated mice; from each mice, 5 samples were taken from left and right ovarian tumors and 3 metastasis sites. Each sample underwent amplicon sequencing of the 4 genes of interst. There are 100 reads in total. 

Note: The reference reads of the 4 genes, primers and sgRNA designed to perform the CRISPR are in the ***MiSeq Primers and sgRNAs.docx*** document. 






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


+ right now I am not sure what each of the lines in input file should correspond.
+ According to the cas-offinder authors, 
" The first line of the input file gives directory path containing chromosomes FASTA files,
The second line indicates the desired pattern including PAM site,
...and following lines are the query sequences and maximum mistmatch numbers, seperated by spaces. (The length of the desired pattern and the query sequences should be the same!)" 

Maybe Katie can provide some insight to this.

After this problem is solved:
+ parse the output file (txt)
+ make comparison (look for lowercase letters, symbols etc coded by the cas-offinder authors
+ visualize
_______________________________

## 2019-01-14
* Use PTEN2 data for now
### to-do:
1. understand what each bio.seq.seq class contains
2. finish compare(myList) function
### goal:
1. compare the reads
1. categorize the type of mutation (insertion, deletion, substitution; also, # of bp changed) for each read
1. visualize the result as smth like pie chart
