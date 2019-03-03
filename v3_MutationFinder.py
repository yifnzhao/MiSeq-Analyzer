#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 18:01:50 2019

@author: yifan
"""

#MutationFinder v3


''' 
an update from version3, can be used for more genes, including
PTEN3, BRCA1, p53
and more if desired



right now analyzing BRCA1
'''

from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import re
import pandas
import matplotlib
import matplotlib.pyplot as plt
import sys
import os

def smaller_query(up5bp, length_ref, querySeq,moreup):
    ''' returns a shorter query sequence ending at PAM site and is 7 bp longer than reference sequence'''
    myregex=up5bp+"\w+"
    if re.search(myregex, querySeq) == None:
        myregex=moreup+"\w+"
    if re.search(myregex, querySeq) == None:
        print("unable to identify the target region for:")
        print(querySeq) 
        return querySeq
    up2end_index = re.search(myregex, querySeq).span()    
    up_index=up2end_index[0]
    end_index= up2end_index[1]
    small=querySeq[up_index:end_index]
    smaller=small[:length_ref+7]
    return smaller

def reverse(seq):
    return seq[::-1]

def get_query_from_PAM(up5bp,reversed_refseq,queryseq,PAM,moreup):
    length_ref=len(reversed_refseq)
    
    smallseq=smaller_query(up5bp,length_ref,queryseq,moreup)
    reversed_smallseq=reverse(smallseq)
    whereisPAM=reversed_smallseq.find(PAM+"AG") #change GC for different genes
    if whereisPAM in [0,1,2,3,4,5]:
        seq_from_PAM=reversed_smallseq[whereisPAM:]
    else:
        seq_from_PAM=reversed_smallseq + "  (requires further investigation, PAM sequence is off)"
    return seq_from_PAM

def find_mismatch_position(reversed_refseq,seq_from_PAM):
    if seq_from_PAM.find(reversed_refseq)!=-1:
        return -1 # no mutation found
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
    
def firsttime_all_mutations(mutation_n, my_txt_dir):
    ''' generates a file with all mutations, more mutations can be added later on upon discovery'''
    print("Start analyzing (first time): "+ my_txt_dir)
    print("We have so far identified "+str(mutation_n)+ " mutations")
    
    
#    # PTEN2:
#    PAM="GGG"
#    up5bp="TGATT"
#    moreup="TGTAA"
#    refseq="TGTGCATATTTATTGCATCGGGG"
#    
    
     # BRCA1:
    PAM="AGG"
    up5bp="AGGAG"
    moreup="GAATG"
    refseq="GCGTCGATCATCCAGAGCGTGGGCTACCGGAACCGTGTCAGAAGG"
    

    
    f=open(my_txt_dir,"r") #directory of txt containing my validated reads    
    reversed_refseq=reverse(refseq)
    all_lines=f.readlines()
    mutations=[]
    for str_line in all_lines:
        line=str_line.split() 
        queryseq=line[0]
        seq_from_PAM=get_query_from_PAM(up5bp,reversed_refseq,queryseq,PAM, moreup) #this shorter sequence is reversed
        alignments = pairwise2.align.globalxx(reversed_refseq,seq_from_PAM)
        print("---01234567890123456789")
        print(format_alignment(*alignments[0]))
        mutation_site=find_mismatch_position(reversed_refseq,seq_from_PAM)
        if mutation_site==-1:        
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
        df.to_csv('BRCA1_all_mutations.csv', sep=',')
    return mutation_n        

def mutations_in_a_sample(mutation_n, inputtxt_dir, all_mutations_csv_dir,output_csv_dir):
    print("Start analyzing: "+ inputtxt_dir)
    print("We have so far identified "+str(mutation_n)+ " mutations")
    
    all_mut_df= pandas.read_csv(all_mutations_csv_dir)
    all_mut=all_mut_df.values.tolist()
    
    #mut_name: 1; mut_pos: 2; mut_seq:3; mut_count:4; mut_first_identified_file:5
#    #for PTEN2:
#    PAM="GGG"
#    up5bp="TGATT"
#    moreup="TGTAA"
#    refseq="TGTGCATATTTATTGCATCGGGG"
    
     # BRCA1:
    PAM="AGG"
    up5bp="AGGAG"
    moreup="GAATG"
    refseq="GCGTCGATCATCCAGAGCGTGGGCTACCGGAACCGTGTCAGAAGG"
    
    
    
    
    
    f=open(inputtxt_dir,"r") #output file that needs to be parsed (i.e. find the mutations) 
    reversed_refseq=reverse(refseq)
    all_lines=f.readlines()
    mutations=[]
    for str_line in all_lines:
        line=str_line.split() 
        queryseq=line[0]
        seq_from_PAM=get_query_from_PAM(up5bp,reversed_refseq,queryseq,PAM, moreup) #this shorter sequence is reversed
        alignments = pairwise2.align.globalxx(reversed_refseq,seq_from_PAM)
        print("---01234567890123456789")
        print(format_alignment(*alignments[0]))
        mutation_site=find_mismatch_position(reversed_refseq,seq_from_PAM)
        if mutation_site==-1:        
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
        
if __name__=="__main__":
    mutation_n = 0  
    os.chdir('/Users/yifan/MiSeq/data/MiSeq20190218/R1/BRCA1/my_txt')
   
   
    stdoutOrigin=sys.stdout 
    sys.stdout = open("BRCA1_log20190228.txt", "w")
   
    firsttime_all_mutations(mutation_n,'BRCA1-1_R1.fastq_input.txt')
   
    
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-1_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-1.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-2_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-2.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-3_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-3.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-4_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-4.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-5_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-5.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-6_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-6.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-7_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-7.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-8_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-8.csv')        
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-9_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-9.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-10_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-10.csv')  
    
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-11_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-11.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-12_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-12.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-13_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-13.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-14_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-14.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-15_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-15.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-16_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-16.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-17_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-17.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-18_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-18.csv')        
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-19_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-19.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-20_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-20.csv')  
    
    
    
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-21_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-21.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-22_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-22.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-23_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-23.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-24_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-24.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'BRCA1-25_R1.fastq_input.txt','BRCA1_all_mutations.csv','BRCA1-25.csv')  
    
  
    
    
    
    
  
    
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
    
    
    