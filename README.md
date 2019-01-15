# mice_seq
* This repo includes my codes for mice sequence alignment. 
* The goal is to identify mutations that occur in ovarian cancer genes in mice treated with CRISPR and to visualize the data in an informative way. The genes of interest are: ***Trp53, PTEN2, PTEN3, BRCA1.***
* There are 5 mutated mice; from each mice, 5 samples were taken from left and right ovarian tumors and 3 metastasis sites. Each sample underwent amplicon sequencing of the 4 genes of interst. There are 100 reads in total. 

Note: The reference reads of the 4 genes, primers and sgRNA designed to perform the CRISPR are in the ***MiSeq Primers and sgRNAs.docx*** document. 



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

