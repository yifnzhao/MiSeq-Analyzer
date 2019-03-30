#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 11:06:53 2019

@author: yifan
"""

#inputValidator

''' 
Unfortunately cas-offinder is not very good at validating sequences
so I had to write this module instead :(
This module is to find all reads in fastq file that have primer sequence
and put them in a txt file
'''


import os
from  random import shuffle
import re
import sys 



def parse(fileName):
    """
    this function parses the file and obtain the reads as a dictionary of {id:seq}
    input: fileName (name of the file that you want to parse)
    output: a dictionary whose keys are ids of reads and values are the correpsonding sequences    
    """
    from Bio import SeqIO
    dictIDSeq={}
    for record in SeqIO.parse(fileName, "fastq"):
        dictIDSeq[record.id]=record.seq
    return dictIDSeq

def dict2list(myDict):
    '''
    this function converts a dicitonary of {id:seq} to a list of sequenced reads
    '''
    return [value for key,value in myDict.items()]


def validated_seq_only(list_all_seq, primer):
    list_validated_seq=[]
    for read in list_all_seq:
        str_read = str(read)
        
        if re.search(primer,str_read):
            list_validated_seq.append(str_read)  
    #print(list_validated_seq)        
    return list_validated_seq
            


def random100(list_validated_seq):
   shuffle(list_validated_seq)
   if len(list_validated_seq)>100:
       mylist=list_validated_seq[0:100]
   else:
       mylist=list_validated_seq
   print("\n")
   print("----------------------------------------------------------------------")
   print("The total number of reads in the input txt file below is: " + str(len(mylist)))
   return mylist




def writeInput(filename, primer):
    
#######RANDOM100##################
#    inputfilename=str(filename)+'_input.txt'
#    with open(inputfilename,'w') as f:
#        mydictIDSeq = parse(filename)
#        allSeq = dict2list(mydictIDSeq)
#        all_validated_seq = validated_seq_only(allSeq, primer)
#        my100random_list_seq = random100(all_validated_seq)
#        for item in my100random_list_seq:
#            f.write("%s\n" %item)
    
    
########ALL VALIDATED READS########
#    inputfilename=str(filename)+'_input_all.txt'
#    with open(inputfilename,'w') as f:
#        mydictIDSeq = parse(filename)
#        allSeq = dict2list(mydictIDSeq)
#        all_validated_seq = validated_seq_only(allSeq, primer)
#        
#        for item in all_validated_seq:
#            f.write("%s\n" %item)
#            
#        print("\n")
#        print("----------------------------------------------------------------------")
#        print("The total number of reads in the input txt file below is: " + str(len(all_validated_seq)))    
#    
    
    
##########ALL READS (NO VALIDATION)#######
    inputfilename=str(filename)+'_input_all.txt'
    with open(inputfilename,'w') as f:
        mydictIDSeq = parse(filename)
        allSeq = dict2list(mydictIDSeq)
        
        for item in allSeq:
            f.write("%s\n" %item)
            
        print("\n")
        print("----------------------------------------------------------------------")
        print("The total number of reads in the input txt file below is: " + str(len(allSeq)))    
    
    
    

if __name__=="__main__":
    
    
    
    stdoutOrigin=sys.stdout 
    sys.stdout = open("log20190330_NOT_validated_inputs.txt", "w")

    PTEN2_primer='GTTGCACAGTATCCTTTTGAAGA'
    PTEN3_primer='AAGAAGTCCTTACATGGGTTGG'
    p53_primer='CTGTGCAGTTGTGGGTCAG' #new left primer
    BRCA1_primer='TTGTGAGCGTTTGAATGAGG'
#    
#    #p53
#    os.chdir('/Users/yifan/MiSeq/data/MiSeq20190218/R1/p53/fastq')
#    p53=[
#'MI.M03992_0353.001.FLD0001.p53-1_R1.fastq',
#'MI.M03992_0353.001.FLD0002.p53-13_R1.fastq',
#'MI.M03992_0353.001.FLD0003.p53-25_R1.fastq',
#'MI.M03992_0353.001.FLD0009.p53-2_R1.fastq',
#'MI.M03992_0353.001.FLD0010.p53-14_R1.fastq',
#'MI.M03992_0353.001.FLD0017.p53-3_R1.fastq',
#'MI.M03992_0353.001.FLD0018.p53-15_R1.fastq',
#'MI.M03992_0353.001.FLD0025.p53-4_R1.fastq',
#'MI.M03992_0353.001.FLD0026.p53-16_R1.fastq',
#'MI.M03992_0353.001.FLD0033.p53-5_R1.fastq',
#'MI.M03992_0353.001.FLD0034.p53-17_R1.fastq',
#'MI.M03992_0353.001.FLD0041.p53-6_R1.fastq',
#'MI.M03992_0353.001.FLD0042.p53-18_R1.fastq',
#'MI.M03992_0353.001.FLD0049.p53-7_R1.fastq',
#'MI.M03992_0353.001.FLD0050.p53-19_R1.fastq',
#'MI.M03992_0353.001.FLD0057.p53-8_R1.fastq',
#'MI.M03992_0353.001.FLD0058.p53-20_R1.fastq',
#'MI.M03992_0353.001.FLD0065.p53-9_R1.fastq',
#'MI.M03992_0353.001.FLD0066.p53-21_R1.fastq',
#'MI.M03992_0353.001.FLD0073.p53-10_R1.fastq',
#'MI.M03992_0353.001.FLD0074.p53-22_R1.fastq',
#'MI.M03992_0353.001.FLD0081.p53-11_R1.fastq',
#'MI.M03992_0353.001.FLD0082.p53-23_R1.fastq',
#'MI.M03992_0353.001.FLD0089.p53-12_R1.fastq',
#'MI.M03992_0353.001.FLD0090.p53-24_R1.fastq',
#'MI.M03992_0353.001.FLD0328.p53-26_FACs_R1.fastq'
#
#]    
#    for input in p53:
#        writeInput(input,p53_primer)
#        print('This input file has been created: m' + input+'_input.txt')        
#              
#   #BRCA1            
#    os.chdir('/Users/yifan/MiSeq/data/MiSeq20190218/R1/BRCA1/fastq')
#    BRCA1=['MI.M03992_0353.001.FLD0008.BRCA1-10_R1.fastq',
#'MI.M03992_0353.001.FLD0016.BRCA1-11_R1.fastq',
#'MI.M03992_0353.001.FLD0024.BRCA1-12_R1.fastq',
#'MI.M03992_0353.001.FLD0031.BRCA1-1_R1.fastq',
#'MI.M03992_0353.001.FLD0032.BRCA1-13_R1.fastq',
#'MI.M03992_0353.001.FLD0039.BRCA1-2_R1.fastq',
#'MI.M03992_0353.001.FLD0040.BRCA1-14_R1.fastq',
#'MI.M03992_0353.001.FLD0047.BRCA1-3_R1.fastq',
#'MI.M03992_0353.001.FLD0048.BRCA1-15_R1.fastq',
#'MI.M03992_0353.001.FLD0055.BRCA1-4_R1.fastq',
#'MI.M03992_0353.001.FLD0056.BRCA1-16_R1.fastq',
#'MI.M03992_0353.001.FLD0063.BRCA1-5_R1.fastq',
#'MI.M03992_0353.001.FLD0064.BRCA1-17_R1.fastq',
#'MI.M03992_0353.001.FLD0071.BRCA1-6_R1.fastq',
#'MI.M03992_0353.001.FLD0072.BRCA1-18_R1.fastq',
#'MI.M03992_0353.001.FLD0079.BRCA1-7_R1.fastq',
#'MI.M03992_0353.001.FLD0080.BRCA1-19_R1.fastq',
#'MI.M03992_0353.001.FLD0087.BRCA1-8_R1.fastq',
#'MI.M03992_0353.001.FLD0088.BRCA1-20_R1.fastq',
#'MI.M03992_0353.001.FLD0095.BRCA1-9_R1.fastq',
#'MI.M03992_0353.001.FLD0096.BRCA1-21_R1.fastq',
#'MI.M03992_0353.001.FLD0296.BRCA1-22_R1.fastq',
#'MI.M03992_0353.001.FLD0304.BRCA1-23_R1.fastq',
#'MI.M03992_0353.001.FLD0312.BRCA1-24_R1.fastq',
#'MI.M03992_0353.001.FLD0320.BRCA1-25_R1.fastq',
#'MI.M03992_0353.001.FLD0352.BRCA1-26_FACs_R1.fastq'
#
#]    
#    for input in BRCA1:
#        writeInput(input,BRCA1_primer)
#        print('This input file has been created: m' + input+'_input.txt')    
#        
     
    #PTEN2   
    os.chdir('/Users/yifan/MiSeq/data/MiSeq20190218/R1/PTEN2/fastq')
    PTEN2 = ['MI.M03992_0353.001.FLD0336.PTEN2-26_FACs_R1.fastq']
    
    
#    PTEN2=['MI.M03992_0353.001.FLD0004.PTEN2-12_R1.fastq',
#'MI.M03992_0353.001.FLD0005.PTEN2-24_R1.fastq',
#'MI.M03992_0353.001.FLD0011.PTEN2-1_R1.fastq',
#'MI.M03992_0353.001.FLD0012.PTEN2-13_R1.fastq',
#'MI.M03992_0353.001.FLD0013.PTEN2-25_R1.fastq',
#'MI.M03992_0353.001.FLD0019.PTEN2-2_R1.fastq',
#'MI.M03992_0353.001.FLD0020.PTEN2-14_R1.fastq',
#'MI.M03992_0353.001.FLD0027.PTEN2-3_R1.fastq',
#'MI.M03992_0353.001.FLD0028.PTEN2-15_R1.fastq',
#'MI.M03992_0353.001.FLD0035.PTEN2-4_R1.fastq',
#'MI.M03992_0353.001.FLD0036.PTEN2-16_R1.fastq',
#'MI.M03992_0353.001.FLD0043.PTEN2-5_R1.fastq',
#'MI.M03992_0353.001.FLD0044.PTEN2-17_R1.fastq',
#'MI.M03992_0353.001.FLD0051.PTEN2-6_R1.fastq',
#'MI.M03992_0353.001.FLD0052.PTEN2-18_R1.fastq',
#'MI.M03992_0353.001.FLD0059.PTEN2-7_R1.fastq',
#'MI.M03992_0353.001.FLD0060.PTEN2-19_R1.fastq',
#'MI.M03992_0353.001.FLD0067.PTEN2-8_R1.fastq',
#'MI.M03992_0353.001.FLD0068.PTEN2-20_R1.fastq',
#'MI.M03992_0353.001.FLD0075.PTEN2-9_R1.fastq',
#'MI.M03992_0353.001.FLD0076.PTEN2-21_R1.fastq',
#'MI.M03992_0353.001.FLD0083.PTEN2-10_R1.fastq',
#'MI.M03992_0353.001.FLD0084.PTEN2-22_R1.fastq',
#'MI.M03992_0353.001.FLD0091.PTEN2-11_R1.fastq',
#'MI.M03992_0353.001.FLD0092.PTEN2-23_R1.fastq',
#'MI.M03992_0353.001.FLD0336.PTEN2-26_FACs_R1.fastq'
#
#]    
    for input in PTEN2:
        writeInput(input,PTEN2_primer)
        print('This input file has been created: m' + input+'_input.txt')        
#    #PTEN3
#    os.chdir('/Users/yifan/MiSeq/data/MiSeq20190218/R1/PTEN3/fastq')
#    PTEN3=['MI.M03992_0353.001.FLD0006.PTEN3-11_R1.fastq',
#'MI.M03992_0353.001.FLD0007.PTEN3-23_R1.fastq',
#'MI.M03992_0353.001.FLD0014.PTEN3-12_R1.fastq',
#'MI.M03992_0353.001.FLD0015.PTEN3-24_R1.fastq',
#'MI.M03992_0353.001.FLD0021.PTEN3-1_R1.fastq',
#'MI.M03992_0353.001.FLD0022.PTEN3-13_R1.fastq',
#'MI.M03992_0353.001.FLD0023.PTEN3-25_R1.fastq',
#'MI.M03992_0353.001.FLD0029.PTEN3-2_R1.fastq',
#'MI.M03992_0353.001.FLD0030.PTEN3-14_R1.fastq',
#'MI.M03992_0353.001.FLD0037.PTEN3-3_R1.fastq',
#'MI.M03992_0353.001.FLD0038.PTEN3-15_R1.fastq',
#'MI.M03992_0353.001.FLD0045.PTEN3-4_R1.fastq',
#'MI.M03992_0353.001.FLD0046.PTEN3-16_R1.fastq',
#'MI.M03992_0353.001.FLD0053.PTEN3-5_R1.fastq',
#'MI.M03992_0353.001.FLD0054.PTEN3-17_R1.fastq',
#'MI.M03992_0353.001.FLD0061.PTEN3-6_R1.fastq',
#'MI.M03992_0353.001.FLD0062.PTEN3-18_R1.fastq',
#'MI.M03992_0353.001.FLD0069.PTEN3-7_R1.fastq',
#'MI.M03992_0353.001.FLD0070.PTEN3-19_R1.fastq',
#'MI.M03992_0353.001.FLD0077.PTEN3-8_R1.fastq',
#'MI.M03992_0353.001.FLD0078.PTEN3-20_R1.fastq',
#'MI.M03992_0353.001.FLD0085.PTEN3-9_R1.fastq',
#'MI.M03992_0353.001.FLD0086.PTEN3-21_R1.fastq',
#'MI.M03992_0353.001.FLD0093.PTEN3-10_R1.fastq',
#'MI.M03992_0353.001.FLD0094.PTEN3-22_R1.fastq',
#'MI.M03992_0353.001.FLD0344.PTEN3-26_FACs_R1.fastq'
#]    
#    for input in PTEN3:
#        writeInput(input,PTEN3_primer)
#        print('This input file has been created: m' + input+'_input.txt')
# 
    sys.stdout.close()
    sys.stdout=stdoutOrigin
