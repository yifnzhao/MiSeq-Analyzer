#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 10 10:26:53 2019

@author: yifan
"""

# MutationFinder

from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import re
from termcolor import colored, cprint
import pandas


def smaller_query(up5bp, length_ref, querySeq):
    ''' returns a shorter query sequence ending at PAM site and is 7 bp longer than reference sequence'''
    myregex=up5bp+"\w+"
    up2end_index = re.search(myregex, querySeq).span()
    up_index=up2end_index[0]
    end_index= up2end_index[1]
    small=querySeq[up_index:end_index]
    smaller=small[:length_ref+7]
    return smaller
def reverse(seq):
    return seq[::-1]

def get_query_from_PAM(up5bp,reversed_refseq,queryseq):
    length_ref=len(reversed_refseq)
    smallseq=smaller_query(up5bp,length_ref,queryseq)
    reversed_smallseq=reverse(smallseq)
    whereisPAM=reversed_smallseq.find(PAM+"GC")
    if whereisPAM in [0,1,2,3,4,5]:
        seq_from_PAM=reversed_smallseq[whereisPAM:]
    else:
        seq_from_PAM=reversed_smallseq + "  (requires further investigation, PAM sequence is off)"
    return seq_from_PAM

def find_mismatch_position(reversed_refseq,seq_from_PAM):
    if seq_from_PAM.find(reversed_refseq)!=-1:
        return -1 # sequences without mutation will not be printed
    for index in range(len(reversed_refseq)):
            if reversed_refseq[index]== seq_from_PAM[index]:
                continue
            else:
                mismatch_site=index-3
                return mismatch_site 
def already_in_mutations_list(mutations_list, mutation_site,seq_from_PAM):
    exist=False
    if mutations_list!=[]:
        for mutation in mutations_list:
            if mutation[1] == mutation_site:
                    if mutation[2]==seq_from_PAM:
                        mutation[3]+=1
                        print("This mutation has been already identified as \""+ mutation[0]+ "\" at position "+str(mutation[1]))
                        print("-------------------------------")
                        return True
    return exist
                                
                




if __name__=="__main__":
    print("This is the start of comparison")
    #PAM=input("Enter the PAM sequence: ")
    #up5bp=input("Enter 5 bp upstream of sgRNA: ")
    #refseq=input("Enter the reference sequence, including PAM: ")
    PAM="GGG"
    up5bp="TGATT"
    refseq="TGTGCATATTTATTGCATCGGGG"
    f=open("/Users/yifan/MiSeq/test/output.txt","r")
    
    reversed_refseq=reverse(refseq)
    length_ref=len(refseq)
    
    all_lines=f.readlines()
    mutations=[]
    for str_line in all_lines:
        line=str_line.split() 
        queryseq=line[0]
        gene_name=line[1]
        refseqline=line[3]
        seq_from_PAM=get_query_from_PAM(up5bp,reversed_refseq,queryseq) #this shorter sequence is reversed
        alignments = pairwise2.align.globalxx(reversed_refseq,seq_from_PAM)
        print("---01234567890123456789")
        print(format_alignment(*alignments[0]))
        mutation_site=find_mismatch_position(reversed_refseq,seq_from_PAM)
        if mutation_site==-1:        
            print("No mutation found :)  Next query:")
            print("    ")
            continue
        else:
            exist=already_in_mutations_list(mutations, mutation_site,seq_from_PAM)
            if exist==True:
                continue
            else:
                print("The first mismatch occurs at position " + str(mutation_site))
                mutation_name=input("What is the name of this mutation? (e.g. AC in, G sub, etc.) ")
                mutations.append([mutation_name,mutation_site,seq_from_PAM, 0]) # note that seq_from_PAM is a reversed seq
                print(mutations)
                print("We found a new mutation! Next Query:")
        
        
        #generate lists for writing to CSV file
        
        mut_names=[]
        mut_sites=[]
        mut_seq=[]
        occurrence_count=[]
        for mutation in mutations:
            mut_names.append(mutation[0])
            mut_sites.append(mutation[1])
            mut_seq.append(mutation[2])
            occurrence_count.append(mutation[3])
    df = pandas.DataFrame(data={"Mutation Name":mut_names, "Mutation Site": mut_sites, "Reversed Shortened Mutated Seq": mut_seq, "Occurence Count -1": occurrence_count})
    df.to_csv('./outputparsed.csv', sep=',')  
        
        
        
        
        
        
        
        