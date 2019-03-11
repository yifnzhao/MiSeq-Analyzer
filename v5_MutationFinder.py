#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 10 18:19:06 2019

@author: yifan
"""

#v5_MutationFinder
# Compared to previous versions, this version increases the query window
# and includes some minor fixes


from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import re
import pandas
import matplotlib
import matplotlib.pyplot as plt
import sys
import os

    
def window(seq, left, right)
    ''' 
    query window (20 bp upstream of PAM AND downstream of the reference sequence)
    seq: query sequence
    left: the position of leftmost bp of interest
        (either the start of reference seq or the end of PAM
    right: the position of rightmost bp of interest
    '''
    leftmost = left - 20
    rightmost = right + 20
    window = seq[leftmost-1: rightmost-1] #b/c str pos starts at 0
    return window


def flip(seq):
    #use this for PTEN2 and BRCA1 (both query and ref)
    return seq[::-1]

def matchPAM(query,PAM, up5, moreup5):
    whereisPAM = query.find(PAM) # returns the position of the first occurrence of PAM
    five_left=query.find(up5)
    ten_left=query.find(moreup5)
    if whereisPAM != -1:
        myquery= query[whereisPAM:] + " P"
    elif five_left != -1:
        myquery = query[five_left:] + " 5"
    elif ten_left != -1:
        myquery = query[tenleft:] + " 10"
    else: myquery = query + " unable to match with PAM sequence: " + PAM
    return myquery

def compare(ref, query):   
    # to be used in find_mismatch_position function below
    for index in range(len(refseq)):
            if refseq[index]== seq_from_PAM[index]:
                continue
            else:
                mismatch_site=index #starting from PAM
                return mismatch_site 
            
def find_mismatch_position(ref,query,up5,moreup5):
    if query.find(refseq)!=-1:
        return None # no mutation found
    if query[-1]=="5": #starts at up5
        myref = up5+ref
        pos = 5-(compare(myref, query)) +1 #pos indicates downstream of PAM
    elif query[-1]=="0":
        myref=moreup5+up5+ref
        pos = 10-(compare(myref, query))+1
    elif query[-1]=="P":
        myref=ref
        pos = -(compare(myref, query))+1 #neg indicates upstream of PAM
    return pos

def already_in_mutations_list(mutations_list, mutation_site,seq_from_PAM):
    exist=False
    if mutations_list!=[]:
        for mutation in mutations_list:
            if mutation[1] == mutation_site:
                    if mutation[2]==seq_from_PAM:
                        mutation[3]+=1
                        print("In this sample, this mutation has been already identified as \""+ str(mutation[0])+ "\" at position "+str(mutation[1]))
                        print("\n\n")
                        return True
    return exist
                                
def already_in_all_mutations_list(sample_mutations_list, all_mut_list, mutation_site,seq_from_PAM):
   
    exist = False
    for mutation in all_mut_list:
        if mutation[2]==mutation_site:
            if mutation[3]==seq_from_PAM:
                mutation[4]+=1
                sample_mutation=mutation[1:]
                sample_mutation[-2]=0 #set occurrence count back to 0

                print("In sample "+mutation[-1]+" , this mutation has been already identified as \""+ str(mutation[1])+ "\" at position "+str(mutation[2]))
                print("\n\n")
                return sample_mutation
    return exist

def firsttime_all_mutations(mutation_n, my_txt_dir, PAM, up5, moreup5, ref, left, right):
    ''' generates a file with all mutations, more mutations can be added later on upon discovery'''
    print("Start analyzing (first time): "+ my_txt_dir)
    print("We have so far identified "+str(mutation_n)+ " mutations")
    

    

    len_ref=len(refseq)
    
    
    f=open(my_txt_dir,"r") #directory of txt containing my validated reads
    all_lines=f.readlines()
    mutations=[]
    for str_line in all_lines:
        line=str_line.split() 
        queryseq=line[0]
        small_query=window(queryseq,left, right)
        
        #for PTEN2 and BRCA1 only: flip!!!
        if (PAM == "GGG") or (PAM == "GGA") # this holds because the PAM seq for the 4 genes are different
            ref = flip(ref)
            small_query = flip(small_query)

        smaller_query=matchPAM(small_query, PAM, up5, moreup5)
        alignments = pairwise2.align.globalxx(refseq,smaller_query)
        print("12345678901234567890")
        if alignments == []:
            print("Unable to align the query with reference sequence: "+line[0])
            print("Error...")
            continue
        print(format_alignment(*alignments[0]))
        mutation_site=find_mismatch_position(ref,smaller_query, up5, moreup5)
        if mutation_site==None:        
            print("No mutation found :)")
            print("\n\n")
            print("Next query:\n")
            continue
        else:
            exist=already_in_mutations_list(mutations, mutation_site,seq_from_PAM)
            if exist==True:
                continue
            else:
                print("The first mismatch occurs at position " + str(mutation_site))
                mutation_n+=1
                mutation_name = str(mutation_n)                
                mutations.append([mutation_name,mutation_site,seq_from_PAM, 0, my_txt_dir]) # note that seq_from_PAM is a reversed seq
                print(mutations)
                print("We found a new mutation! Next Query:")
        #writing to CSV file
        mutations=[[None, None,None,None,None]]
        df =pandas.DataFrame(mutations, columns=['Mutation Name', 'Mutation Site','Reversed Shortened Mutated Seq', 'Occurrence Count -1','First identified in:'])
        df.to_csv('PTEN3_all_mutations.csv', sep=',')
    return mutation_n        

def mutations_in_a_sample(mutation_n, inputtxt_dir, all_mutations_csv_dir,output_csv_dir, PAM, up5, moreup5, ref, left, right):
    
    print("Start analyzing: "+ inputtxt_dir)
    print("We have so far identified "+str(mutation_n)+ " mutations")
    
    all_mut_df= pandas.read_csv(all_mutations_csv_dir)
    all_mut=all_mut_df.values.tolist()
    

    
    
    len_ref=len(refseq)

    
    f=open(inputtxt_dir,"r") #output file that needs to be parsed (i.e. find the mutations) 
    
    all_lines=f.readlines()
    mutations=[]
    for str_line in all_lines:
        line=str_line.split() 
        queryseq=line[0]
        small_query=window(queryseq, left, right)
        
         
        #for PTEN2 and BRCA1 only: flip!!!
        ref = flip(ref)
        small_query = flip(small_query)
        # comment the two lines above for PTEN3 and p53!!!
        
        
        samller_query=matchPAM(small_query, PAM, up5, moreup5)
        alignments = pairwise2.align.globalxx(refseq,seq_from_PAM)
        
        print("12345678901234567890")
        #print(alignments)
        
        if alignments == []:
            print("Unable to align the query with reference sequence: "+line[0])
            print("Error...")
            continue
        print(format_alignment(*alignments[0]))
        mutation_site=find_mismatch_position(red,samller_query, up5, moreup5)
        if mutation_site==None:        
           print("No mutation found :)")
           print("\n\n")
           print("Next query:\n")
        else:
            exist_in_sample_mutations_list=already_in_mutations_list(mutations, mutation_site,seq_from_PAM)
            
            if exist_in_sample_mutations_list==True:
                continue
            else:
                exist_in_all_mutations_list=already_in_all_mutations_list(mutations, all_mut, mutation_site, seq_from_PAM)
                if  exist_in_all_mutations_list != False:
                    mutations.append(exist_in_all_mutations_list)
                else:
                    print("The first mismatch occurs at position " + str(mutation_site))
                    mutation_n+=1
                    mutation_name = str(mutation_n)
                    mutations.append([mutation_name,mutation_site,seq_from_PAM, 0, inputtxt_dir])
                    # note that seq_from_PAM is a reversed seq
                    print("We found a new mutation!")
                    print("\n\n")
                    print("Next query:\n")
        
    #print(mutations)  

    #writing to CSV file
  
    df =pandas.DataFrame(mutations, columns=['Mutation Name', 'Mutation Site',
                                             'Reversed Shortened Mutated Seq', 
                                             'Occurrence Count -1','First identified in:'
                                             ])
    
    df.to_csv(output_csv_dir, sep=',') 
    
    for mut in all_mut:
        del mut[0]
    for mymutation in mutations:
        all_mut.append(mymutation)

    updated_df_allmut=pandas.DataFrame(all_mut, columns=['Mutation Name', 
                                                         'Mutation Site',
                                                         'Reversed Shortened Mutated Seq', 
                                                         'Occurrence Count -1',
                                                         'First identified in:'])
    updated_df_allmut.to_csv(all_mutations_csv_dir, sep=',')
    return mutation_n


def piechart(csv_dir, piechart_dir):
    ''' this method generates a pie chart for a given csv file containing the following info:    
    mutation name
    mutation site
    reversed shortened mutated seq
    occurrence count-1
    
    reference: https://matplotlib.org/gallery/pie_and_polar_charts/pie_features.html#sphx-glr-gallery-pie-and-polar-charts-pie-features-py
    
    To generate a pie chart, we will need labels and sizes
    '''
    #load csv

    df=pandas.read_csv(csv_dir)
    filelist=df.values.tolist()
    
    labels=[]
    counts=[]
    #total count
    #total = sum(1 for line in open(txt_dir))
    total=100
    wt = total
    for entry in filelist:
        # 0:index, 1: name, 2:site, 3:seq, 4: count-1
        if entry[4]==0:
            continue
        entry[4]+=1 #all counts add 1
        labels.append(str(int(entry[1]))+" at site "+str(int(entry[2])))
        counts.append(int(entry[4])/total)
        wt-= entry[4]
    labels.append("WT")
    counts.append(wt/total)
    
    matplotlib.rcParams['font.size'] = 7.0    
    fig1, ax1 = plt.subplots()
    ax1.pie(counts, explode = None, labels = labels, autopct ='%1.1f%%',
            shadow = False, startangle=90)
    ax1.axis('equal') # Equal aspect ratio ensures that pie is drawn as a circle
    plt.savefig(piechart_dir, dpi=1000)
    plt.show()
    
def BRCA1():
        # BRCA1:
    PAM="GGA" #flipped
    up5="CCTTC" #flipped
    moreup5="GGACC" #flipped
    ref="GCGTCGATCATCCAGAGCGTGGGCTACCGGAACCGTGTCAGAAGG"
    left=83
    right=128
    mutation_n = 0  
    os.chdir('/Users/yifan/MiSeq/data/MiSeq20190218/R1/BRCA1/my_txt')
   
   
    stdoutOrigin=sys.stdout 
    sys.stdout = open("BRCA1_log20190311.txt", "w")
   
    firsttime_all_mutations(mutation_n,'BRCA1-1_R1.fastq_input.txt',PAM, up5, moreup5, ref, left, right)
   
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-1_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-1.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-2_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-2.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-3_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-3.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-4_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-4.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-5_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-5.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-6_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-6.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-7_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-7.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-8_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-8.csv',PAM, up5, moreup5, ref, left, right)        
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-9_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-9.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-10_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-10.csv',PAM, up5, moreup5, ref, left, right)  
    
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-11_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-11.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-12_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-12.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-13_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-13.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-14_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-14.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-15_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-15.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-16_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-16.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-17_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-17.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-18_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-18.csv',PAM, up5, moreup5, ref, left, right)        
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-19_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-19.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-20_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-20.csv',PAM, up5, moreup5, ref, left, right)  
    
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-21_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-21.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-22_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-22.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-23_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-23.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-24_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-24.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-25_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-25.csv',PAM, up5, moreup5, ref, left, right)  
    
    piechart('BRCA1-1.csv','BRCA1-1.png');
    piechart('BRCA1-2.csv','BRCA1-2.png');
    piechart('BRCA1-3.csv','BRCA1-3.png');
    piechart('BRCA1-4.csv','BRCA1-4.png');
    piechart('BRCA1-5.csv','BRCA1-5.png');
    piechart('BRCA1-6.csv','BRCA1-6.png');
    piechart('BRCA1-7.csv','BRCA1-7.png');
    piechart('BRCA1-8.csv','BRCA1-8.png');
    piechart('BRCA1-9.csv','BRCA1-9.png');
    piechart('BRCA1-10.csv','BRCA1-10.png');
    piechart('BRCA1-11.csv','BRCA1-11.png');
    piechart('BRCA1-12.csv','BRCA1-12.png');
    piechart('BRCA1-13.csv','BRCA1-13.png');
    piechart('BRCA1-14.csv','BRCA1-14.png');
    piechart('BRCA1-15.csv','BRCA1-15.png');
    piechart('BRCA1-16.csv','BRCA1-16.png');
    piechart('BRCA1-17.csv','BRCA1-17.png');
    piechart('BRCA1-18.csv','BRCA1-18.png');
    piechart('BRCA1-19.csv','BRCA1-19.png');
    piechart('BRCA1-20.csv','BRCA1-20.png');
    piechart('BRCA1-21.csv','BRCA1-21.png');
    piechart('BRCA1-22.csv','BRCA1-22.png');
    piechart('BRCA1-23.csv','BRCA1-23.png');
    piechart('BRCA1-24.csv','BRCA1-24.png');
    piechart('BRCA1-25.csv','BRCA1-25.png');

    sys.stdout.close()
    sys.stdout=stdoutOrigin
    print("BRCA1 finished")
    
def PTEN3():     
        # PTEN3:
    PAM="CCT"
    up5="CAGAT"
    moreup5="CCACA"
    ref="CCTCAGTTTGTGGTCTGCCAGCT"
    left=189
    right=212
       
    mutation_n = 0  
    os.chdir('/Users/yifan/MiSeq/data/MiSeq20190218/R1/PTEN3/my_txt')
   
   
    stdoutOrigin=sys.stdout 
    sys.stdout = open("PTEN3_log20190311.txt", "w")
   
    firsttime_all_mutations(mutation_n,'PTEN3-1_R1.fastq_input.txt',PAM, up5, moreup5, ref, left, right)
   
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-1_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-1.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-2_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-2.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-3_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-3.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-4_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-4.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-5_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-5.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-6_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-6.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-7_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-7.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-8_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-8.csv',PAM, up5, moreup5, ref, left, right)        
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-9_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-9.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-10_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-10.csv',PAM, up5, moreup5, ref, left, right)  
    
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-11_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-11.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-12_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-12.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-13_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-13.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-14_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-14.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-15_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-15.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-16_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-16.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-17_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-17.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-18_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-18.csv',PAM, up5, moreup5, ref, left, right)        
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-19_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-19.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-20_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-20.csv',PAM, up5, moreup5, ref, left, right)  
    
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-21_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-21.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-22_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-22.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-23_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-23.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-24_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-24.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-25_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-25.csv',PAM, up5, moreup5, ref, left, right)  
    
    piechart('PTEN3-1.csv','PTEN3-1.png');
    piechart('PTEN3-2.csv','PTEN3-2.png');
    piechart('PTEN3-3.csv','PTEN3-3.png');
    piechart('PTEN3-4.csv','PTEN3-4.png');
    piechart('PTEN3-5.csv','PTEN3-5.png');
    piechart('PTEN3-6.csv','PTEN3-6.png');
    piechart('PTEN3-7.csv','PTEN3-7.png');
    piechart('PTEN3-8.csv','PTEN3-8.png');
    piechart('PTEN3-9.csv','PTEN3-9.png');
    piechart('PTEN3-10.csv','PTEN3-10.png');
    piechart('PTEN3-11.csv','PTEN3-11.png');
    piechart('PTEN3-12.csv','PTEN3-12.png');
    piechart('PTEN3-13.csv','PTEN3-13.png');
    piechart('PTEN3-14.csv','PTEN3-14.png');
    piechart('PTEN3-15.csv','PTEN3-15.png');
    piechart('PTEN3-16.csv','PTEN3-16.png');
    piechart('PTEN3-17.csv','PTEN3-17.png');
    piechart('PTEN3-18.csv','PTEN3-18.png');
    piechart('PTEN3-19.csv','PTEN3-19.png');
    piechart('PTEN3-20.csv','PTEN3-20.png');
    piechart('PTEN3-21.csv','PTEN3-21.png');
    piechart('PTEN3-22.csv','PTEN3-22.png');
    piechart('PTEN3-23.csv','PTEN3-23.png');
    piechart('PTEN3-24.csv','PTEN3-24.png');
    piechart('PTEN3-25.csv','PTEN3-25.png');

    sys.stdout.close()
    sys.stdout=stdoutOrigin
    print("PTEN3 finished")
    
def p53():     

         # p53:
    PAM="CCA" #shortened for consistency
    up5="TGATT"
    moreup5="TGTAA"
    ref="CCCCCACCATGAGCGCTGCTCCGATG"
    left=107
    right=130
       
    
    mutation_n = 0  
    os.chdir('/Users/yifan/MiSeq/data/MiSeq20190218/R1/p53/my_txt')
   
   
    stdoutOrigin=sys.stdout 
    sys.stdout = open("p53_log20190311.txt", "w")
    firsttime_all_mutations(mutation_n,'p53-1_R1.fastq_input.txt',PAM, up5, moreup5, ref, left, right)
   
    
    mutation_n=mutations_in_a_sample(mutation_n, 'p53-1_R1.fastq_input.txt','p53_all_mutations.csv','p53-1.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'p53-2_R1.fastq_input.txt','p53_all_mutations.csv','p53-2.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'p53-3_R1.fastq_input.txt','p53_all_mutations.csv','p53-3.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'p53-4_R1.fastq_input.txt','p53_all_mutations.csv','p53-4.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'p53-5_R1.fastq_input.txt','p53_all_mutations.csv','p53-5.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'p53-6_R1.fastq_input.txt','p53_all_mutations.csv','p53-6.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'p53-7_R1.fastq_input.txt','p53_all_mutations.csv','p53-7.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'p53-8_R1.fastq_input.txt','p53_all_mutations.csv','p53-8.csv',PAM, up5, moreup5, ref, left, right)        
    mutation_n=mutations_in_a_sample(mutation_n, 'p53-9_R1.fastq_input.txt','p53_all_mutations.csv','p53-9.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'p53-10_R1.fastq_input.txt','p53_all_mutations.csv','p53-10.csv',PAM, up5, moreup5, ref, left, right)  
    
    mutation_n=mutations_in_a_sample(mutation_n, 'p53-11_R1.fastq_input.txt','p53_all_mutations.csv','p53-11.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'p53-12_R1.fastq_input.txt','p53_all_mutations.csv','p53-12.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'p53-13_R1.fastq_input.txt','p53_all_mutations.csv','p53-13.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'p53-14_R1.fastq_input.txt','p53_all_mutations.csv','p53-14.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'p53-15_R1.fastq_input.txt','p53_all_mutations.csv','p53-15.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'p53-16_R1.fastq_input.txt','p53_all_mutations.csv','p53-16.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'p53-17_R1.fastq_input.txt','p53_all_mutations.csv','p53-17.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'p53-18_R1.fastq_input.txt','p53_all_mutations.csv','p53-18.csv',PAM, up5, moreup5, ref, left, right)        
    mutation_n=mutations_in_a_sample(mutation_n, 'p53-19_R1.fastq_input.txt','p53_all_mutations.csv','p53-19.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'p53-20_R1.fastq_input.txt','p53_all_mutations.csv','p53-20.csv',PAM, up5, moreup5, ref, left, right)  
    
    mutation_n=mutations_in_a_sample(mutation_n, 'p53-21_R1.fastq_input.txt','p53_all_mutations.csv','p53-21.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'p53-22_R1.fastq_input.txt','p53_all_mutations.csv','p53-22.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'p53-23_R1.fastq_input.txt','p53_all_mutations.csv','p53-23.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'p53-24_R1.fastq_input.txt','p53_all_mutations.csv','p53-24.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'p53-25_R1.fastq_input.txt','p53_all_mutations.csv','p53-25.csv',PAM, up5, moreup5, ref, left, right)  

    piechart('p53-1.csv','p53-1.png');
    piechart('p53-2.csv','p53-2.png');
    piechart('p53-3.csv','p53-3.png');
    piechart('p53-4.csv','p53-4.png');
    piechart('p53-5.csv','p53-5.png');
    piechart('p53-6.csv','p53-6.png');
    piechart('p53-7.csv','p53-7.png');
    piechart('p53-8.csv','p53-8.png');
    piechart('p53-9.csv','p53-9.png');
    piechart('p53-10.csv','p53-10.png');
    piechart('p53-11.csv','p53-11.png');
    piechart('p53-12.csv','p53-12.png');
    piechart('p53-13.csv','p53-13.png');
    piechart('p53-14.csv','p53-14.png');
    piechart('p53-15.csv','p53-15.png');
    piechart('p53-16.csv','p53-16.png');
    piechart('p53-17.csv','p53-17.png');
    piechart('p53-18.csv','p53-18.png');
    piechart('p53-19.csv','p53-19.png');
    piechart('p53-20.csv','p53-20.png');
    piechart('p53-21.csv','p53-21.png');
    piechart('p53-22.csv','p53-22.png');
    piechart('p53-23.csv','p53-23.png');
    piechart('p53-24.csv','p53-24.png');
    piechart('p53-25.csv','p53-25.png');

    sys.stdout.close()
    sys.stdout=stdoutOrigin
    print("p53 finished")
    
    
def PTEN2():  
     # PTEN2:
    PAM="GGG"
    up5="TAAAC" #flipped
    moreup5="ATTTT" #flipped
    ref="TGTGCATATTTATTGCATCGGGG"
    left=154
    right=177
    
    mutation_n = 0  
    os.chdir('/Users/yifan/MiSeq/data/MiSeq20190218/R1/PTEN2/my_txt')
    stdoutOrigin=sys.stdout 
    sys.stdout = open("PTEN2_log2019011.txt", "w")
    
    
    firsttime_all_mutations(mutation_n,'MI.M03992_0353.001.FLD0011.PTEN2-1_R1.fastq_input.txt',PAM, up5, moreup5, ref, left, right)
   
    #1-5
    mutation_n=mutations_in_a_sample(mutation_n, 'MI.M03992_0353.001.FLD0011.PTEN2-1_R1.fastq_input.txt','PTEN2_all_mutations.csv','PTEN2-1.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'MI.M03992_0353.001.FLD0019.PTEN2-2_R1.fastq_input.txt','PTEN2_all_mutations.csv','PTEN2-2.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'MI.M03992_0353.001.FLD0027.PTEN2-3_R1.fastq_input.txt','PTEN2_all_mutations.csv','PTEN2-3.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'MI.M03992_0353.001.FLD0035.PTEN2-4_R1.fastq_input.txt','PTEN2_all_mutations.csv','PTEN2-4.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'MI.M03992_0353.001.FLD0043.PTEN2-5_R1.fastq_input.txt','PTEN2_all_mutations.csv','PTEN2-5.csv',PAM, up5, moreup5, ref, left, right)  

    
    #6-10
    mutation_n=mutations_in_a_sample(mutation_n, 'MI.M03992_0353.001.FLD0051.PTEN2-6_R1.fastq_input.txt','PTEN2_all_mutations.csv','PTEN2-6.csv',PAM, up5, moreup5, ref, left, right)
    mutation_n=mutations_in_a_sample(mutation_n, 'MI.M03992_0353.001.FLD0059.PTEN2-7_R1.fastq_input.txt','PTEN2_all_mutations.csv','PTEN2-7.csv',PAM, up5, moreup5, ref, left, right)
    mutation_n=mutations_in_a_sample(mutation_n, 'MI.M03992_0353.001.FLD0067.PTEN2-8_R1.fastq_input.txt','PTEN2_all_mutations.csv','PTEN2-8.csv',PAM, up5, moreup5, ref, left, right)    
    mutation_n=mutations_in_a_sample(mutation_n, 'MI.M03992_0353.001.FLD0075.PTEN2-9_R1.fastq_input.txt','PTEN2_all_mutations.csv','PTEN2-9.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'MI.M03992_0353.001.FLD0083.PTEN2-10_R1.fastq_input.txt','PTEN2_all_mutations.csv','PTEN2-10.csv',PAM, up5, moreup5, ref, left, right)  
    
    #11-15
    mutation_n=mutations_in_a_sample(mutation_n, 'MI.M03992_0353.001.FLD0091.PTEN2-11_R1.fastq_input.txt','PTEN2_all_mutations.csv','PTEN2-11.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'MI.M03992_0353.001.FLD0004.PTEN2-12_R1.fastq_input.txt','PTEN2_all_mutations.csv','PTEN2-12.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'MI.M03992_0353.001.FLD0012.PTEN2-13_R1.fastq_input.txt','PTEN2_all_mutations.csv','PTEN2-13.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'MI.M03992_0353.001.FLD0020.PTEN2-14_R1.fastq_input.txt','PTEN2_all_mutations.csv','PTEN2-14.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'MI.M03992_0353.001.FLD0028.PTEN2-15_R1.fastq_input.txt','PTEN2_all_mutations.csv','PTEN2-15.csv',PAM, up5, moreup5, ref, left, right)  
    
    #16-20
    mutation_n=mutations_in_a_sample(mutation_n, 'MI.M03992_0353.001.FLD0036.PTEN2-16_R1.fastq_input.txt','PTEN2_all_mutations.csv','PTEN2-16.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'MI.M03992_0353.001.FLD0044.PTEN2-17_R1.fastq_input.txt','PTEN2_all_mutations.csv','PTEN2-17.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'MI.M03992_0353.001.FLD0052.PTEN2-18_R1.fastq_input.txt','PTEN2_all_mutations.csv','PTEN2-18.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'MI.M03992_0353.001.FLD0060.PTEN2-19_R1.fastq_input.txt','PTEN2_all_mutations.csv','PTEN2-19.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'MI.M03992_0353.001.FLD0068.PTEN2-20_R1.fastq_input.txt','PTEN2_all_mutations.csv','PTEN2-20.csv',PAM, up5, moreup5, ref, left, right)  
    
    #21-25
    mutation_n=mutations_in_a_sample(mutation_n, 'MI.M03992_0353.001.FLD0076.PTEN2-21_R1.fastq_input.txt','PTEN2_all_mutations.csv','PTEN2-21.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'MI.M03992_0353.001.FLD0084.PTEN2-22_R1.fastq_input.txt','PTEN2_all_mutations.csv','PTEN2-22.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'MI.M03992_0353.001.FLD0092.PTEN2-23_R1.fastq_input.txt','PTEN2_all_mutations.csv','PTEN2-23.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'MI.M03992_0353.001.FLD0005.PTEN2-24_R1.fastq_input.txt','PTEN2_all_mutations.csv','PTEN2-24.csv',PAM, up5, moreup5, ref, left, right)  
    mutation_n=mutations_in_a_sample(mutation_n, 'MI.M03992_0353.001.FLD0013.PTEN2-25_R1.fastq_input.txt','PTEN2_all_mutations.csv','PTEN2-25.csv',PAM, up5, moreup5, ref, left, right)  
    
    
    #mutation_n=mutations_in_a_sample(mutation_n, '','PTEN2_all_mutations.csv','PTEN2-3.csv')  
    
    piechart('PTEN2-1.csv','PTEN2-1.png');
    piechart('PTEN2-2.csv','PTEN2-2.png');
    piechart('PTEN2-3.csv','PTEN2-3.png');
    piechart('PTEN2-4.csv','PTEN2-4.png');
    piechart('PTEN2-5.csv','PTEN2-5.png');
    piechart('PTEN2-6.csv','PTEN2-6.png');
    piechart('PTEN2-7.csv','PTEN2-7.png');
    piechart('PTEN2-8.csv','PTEN2-8.png');
    piechart('PTEN2-9.csv','PTEN2-9.png');
    piechart('PTEN2-10.csv','PTEN2-10.png');
    piechart('PTEN2-11.csv','PTEN2-11.png');
    piechart('PTEN2-12.csv','PTEN2-12.png');
    piechart('PTEN2-13.csv','PTEN2-13.png');
    piechart('PTEN2-14.csv','PTEN2-14.png');
    piechart('PTEN2-15.csv','PTEN2-15.png');
    piechart('PTEN2-16.csv','PTEN2-16.png');
    piechart('PTEN2-17.csv','PTEN2-17.png');
    piechart('PTEN2-18.csv','PTEN2-18.png');
    piechart('PTEN2-19.csv','PTEN2-19.png');
    piechart('PTEN2-20.csv','PTEN2-20.png');
    piechart('PTEN2-21.csv','PTEN2-21.png');
    piechart('PTEN2-22.csv','PTEN2-22.png');
    piechart('PTEN2-23.csv','PTEN2-23.png');
    piechart('PTEN2-24.csv','PTEN2-24.png');
    piechart('PTEN2-25.csv','PTEN2-25.png');

    sys.stdout.close()
    sys.stdout=stdoutOrigin
    
    print("PTEN2 finished")

if name=="__main__":
    BRCA1()
    p53()
    PTEN3()
    PTEN2()
    print("Mutation Finder has exited")
    
 