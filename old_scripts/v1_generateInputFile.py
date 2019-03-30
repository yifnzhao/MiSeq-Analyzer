#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 10 10:26:53 2019

@author: yifan
"""

# generateInputFile

'''
This module generates an input file for cas-offinder
'''
import os
from  random import shuffle
import pandas as pd

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
    """
    this function converts a dicitonary of {id:seq} to a list of sequenced reads
    """
    return [value for key,value in myDict.items()]


def writeInput(filename):
    inputfilename=str(filename)+'_input.txt'
    with open(inputfilename,'w') as f:
        f.write('/Users/yifan/MiSeq/myreference/\n')
        f.write('NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n')
        mydictIDSeq = parse(filename)
        myListSeq = dict2list(mydictIDSeq)
        for item in myListSeq:
            f.write("%s 100\n" % item)

def writeInput100random(filename, myListSeq):
    inputfilename='./100random/'+str(filename)+'_100random_input.txt'
    with open(inputfilename,'w') as f:
        f.write('/Users/yifan/MiSeq/myreference/\n')
        f.write('NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n')    
        for item in myListSeq:
            seq=item[0]
            f.write("%s\n" % seq)
            
            
            
    
    

if __name__=="__main__":
#    os.chdir('/Users/yifan/MiSeq/data/MiSeq20190218/R1/p53')
#    allinputs=[
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
#'MI.M03992_0353.001.FLD0090.p53-24_R1.fastq'
#
#]
#
#    
#    for input in allinputs:
#        writeInput(input)
#        print('This input file has been created:' + input)
    
    
#
#    directory = '/Users/yifan/MiSeq/data/MiSeq20190218/R1/PTEN3/txt/all/'
#    os.chdir(directory)
#
#    #get 100 reads randomly
#    for filename in os.listdir(directory):
#        if filename.endswith("input.txt"):
#            
#            allseq = pd.read_csv(filename, sep = "\n" ).values.tolist()
#            del allseq[0]
#            shuffle(allseq)
#            myListSeq100=allseq[0:100]
#            print(len(myListSeq100))
#            #print(myListSeq100)
#            writeInput100random(filename, myListSeq100)
#            print('The input file for 100 random sequences has been created: ' + filename)
#                
        
        
        