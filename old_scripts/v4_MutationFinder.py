#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  3 17:44:21 2019

@author: yifan
"""

#v4_MutationFinder
'''this version does not reverse the reference sequence (for PTEN3 and p53)'''



from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import re
import pandas
import matplotlib
import matplotlib.pyplot as plt
import sys
import os

def query_from_PAM(PAM, length_ref, querySeq):
    ''' returns a shorter query sequence starting at PAM site and is 7 bp longer than reference sequence'''
    
    PAM_pos = querySeq.find(PAM+"CA") # change this for diff genes
#p53
#    if PAM_pos == -1:
#        PAM_pos_13up = querySeq.find("TCGTG")
#        small_seq=querySeq[PAM_pos_13up + 13:]
#    else:
#        small_seq=querySeq[PAM_pos:]
    
 #PTEN3
    if PAM_pos == -1:
        PAM_pos_22up = querySeq.find("ATATGTATTT")
        small_seq=querySeq[PAM_pos_22up + 22:]
    else:
        small_seq=querySeq[PAM_pos:]   
    

    smaller_seq=small_seq[:length_ref+7]
    return smaller_seq




def find_mismatch_position(refseq,seq_from_PAM):
    if seq_from_PAM.find(refseq)!=-1:
        return -1 # no mutation found
    for index in range(len(refseq)):
            if refseq[index]== seq_from_PAM[index]:
              
                continue
            else:
                mismatch_site=index #starting from PAM
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
    
     # PTEN3:
    PAM="CAGATCCT"
    
    refseq="CCTCAGTTTGTGGTCTGCCAGCT"
    
#    
#     # p53:
#    PAM="CCCCCA"
#    refseq="CCCCCACCATGAGCGCTGCTCCGATG"
    len_ref=len(refseq)
    
    
    
    
    
    f=open(my_txt_dir,"r") #directory of txt containing my validated reads    
#    reversed_refseq=reverse(refseq)
    all_lines=f.readlines()
    mutations=[]
    for str_line in all_lines:
        line=str_line.split() 
        queryseq=line[0]
        seq_from_PAM=query_from_PAM(PAM,len_ref, queryseq) #this shorter sequence is  NOT reversed
        alignments = pairwise2.align.globalxx(refseq,seq_from_PAM)
       # print("---01234567890123456789")
        print(format_alignment(*alignments[0]))
        mutation_site=find_mismatch_position(refseq,seq_from_PAM)
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
        df.to_csv('PTEN3_all_mutations.csv', sep=',')
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
    
      # PTEN3:
    PAM="CAGATCCT"
    
    refseq="CCTCAGTTTGTGGTCTGCCAGCT"
#     # p53:
#    PAM="CCCCCA"
#    
#    refseq="CCCCCACCATGAGCGCTGCTCCGATG"
    len_ref=len(refseq)

    
    f=open(inputtxt_dir,"r") #output file that needs to be parsed (i.e. find the mutations) 
#    reversed_refseq=reverse(refseq)
    all_lines=f.readlines()
    mutations=[]
    for str_line in all_lines:
        line=str_line.split() 
        queryseq=line[0]
        seq_from_PAM=query_from_PAM(PAM,len_ref, queryseq) #this shorter sequence is  NOT reversed
        alignments = pairwise2.align.globalxx(refseq,seq_from_PAM)
        #print("---01234567890123456789")
        #print(alignments)
        if alignments == []:
            print("Unable to align the query with reference sequence: "+line[0])
            print("Error...")
            continue
        print(format_alignment(*alignments[0]))
        mutation_site=find_mismatch_position(refseq,seq_from_PAM)
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
    os.chdir('/Users/yifan/MiSeq/data/MiSeq20190218/R1/PTEN3/my_txt')
   
   
    stdoutOrigin=sys.stdout 
    sys.stdout = open("PTEN3_log20190304.txt", "w")
   
    firsttime_all_mutations(mutation_n,'PTEN3-1_R1.fastq_input.txt')
   
    
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-1_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-1.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-2_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-2.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-3_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-3.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-4_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-4.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-5_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-5.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-6_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-6.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-7_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-7.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-8_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-8.csv')        
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-9_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-9.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-10_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-10.csv')  
    
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-11_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-11.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-12_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-12.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-13_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-13.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-14_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-14.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-15_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-15.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-16_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-16.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-17_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-17.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-18_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-18.csv')        
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-19_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-19.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-20_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-20.csv')  
    
    
    
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-21_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-21.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-22_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-22.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-23_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-23.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-24_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-24.csv')  
    mutation_n=mutations_in_a_sample(mutation_n, 'PTEN3-25_R1.fastq_input.txt','PTEN3_all_mutations.csv','PTEN3-25.csv')  
  
    
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


# 
#    firsttime_all_mutations(mutation_n,'p53-1_R1.fastq_input.txt')
#   
#    
#    mutation_n=mutations_in_a_sample(mutation_n, 'p53-1_R1.fastq_input.txt','p53_all_mutations.csv','p53-1.csv')  
#    mutation_n=mutations_in_a_sample(mutation_n, 'p53-2_R1.fastq_input.txt','p53_all_mutations.csv','p53-2.csv')  
#    mutation_n=mutations_in_a_sample(mutation_n, 'p53-3_R1.fastq_input.txt','p53_all_mutations.csv','p53-3.csv')  
#    mutation_n=mutations_in_a_sample(mutation_n, 'p53-4_R1.fastq_input.txt','p53_all_mutations.csv','p53-4.csv')  
#    mutation_n=mutations_in_a_sample(mutation_n, 'p53-5_R1.fastq_input.txt','p53_all_mutations.csv','p53-5.csv')  
#    mutation_n=mutations_in_a_sample(mutation_n, 'p53-6_R1.fastq_input.txt','p53_all_mutations.csv','p53-6.csv')  
#    mutation_n=mutations_in_a_sample(mutation_n, 'p53-7_R1.fastq_input.txt','p53_all_mutations.csv','p53-7.csv')  
#    mutation_n=mutations_in_a_sample(mutation_n, 'p53-8_R1.fastq_input.txt','p53_all_mutations.csv','p53-8.csv')        
#    mutation_n=mutations_in_a_sample(mutation_n, 'p53-9_R1.fastq_input.txt','p53_all_mutations.csv','p53-9.csv')  
#    mutation_n=mutations_in_a_sample(mutation_n, 'p53-10_R1.fastq_input.txt','p53_all_mutations.csv','p53-10.csv')  
#    
#    mutation_n=mutations_in_a_sample(mutation_n, 'p53-11_R1.fastq_input.txt','p53_all_mutations.csv','p53-11.csv')  
#    mutation_n=mutations_in_a_sample(mutation_n, 'p53-12_R1.fastq_input.txt','p53_all_mutations.csv','p53-12.csv')  
#    mutation_n=mutations_in_a_sample(mutation_n, 'p53-13_R1.fastq_input.txt','p53_all_mutations.csv','p53-13.csv')  
#    mutation_n=mutations_in_a_sample(mutation_n, 'p53-14_R1.fastq_input.txt','p53_all_mutations.csv','p53-14.csv')  
#    mutation_n=mutations_in_a_sample(mutation_n, 'p53-15_R1.fastq_input.txt','p53_all_mutations.csv','p53-15.csv')  
#    mutation_n=mutations_in_a_sample(mutation_n, 'p53-16_R1.fastq_input.txt','p53_all_mutations.csv','p53-16.csv')  
#    mutation_n=mutations_in_a_sample(mutation_n, 'p53-17_R1.fastq_input.txt','p53_all_mutations.csv','p53-17.csv')  
#    mutation_n=mutations_in_a_sample(mutation_n, 'p53-18_R1.fastq_input.txt','p53_all_mutations.csv','p53-18.csv')        
#    mutation_n=mutations_in_a_sample(mutation_n, 'p53-19_R1.fastq_input.txt','p53_all_mutations.csv','p53-19.csv')  
#    mutation_n=mutations_in_a_sample(mutation_n, 'p53-20_R1.fastq_input.txt','p53_all_mutations.csv','p53-20.csv')  
#    
#    
#    
#    mutation_n=mutations_in_a_sample(mutation_n, 'p53-21_R1.fastq_input.txt','p53_all_mutations.csv','p53-21.csv')  
#    mutation_n=mutations_in_a_sample(mutation_n, 'p53-22_R1.fastq_input.txt','p53_all_mutations.csv','p53-22.csv')  
#    mutation_n=mutations_in_a_sample(mutation_n, 'p53-23_R1.fastq_input.txt','p53_all_mutations.csv','p53-23.csv')  
#    mutation_n=mutations_in_a_sample(mutation_n, 'p53-24_R1.fastq_input.txt','p53_all_mutations.csv','p53-24.csv')  
#    mutation_n=mutations_in_a_sample(mutation_n, 'p53-25_R1.fastq_input.txt','p53_all_mutations.csv','p53-25.csv')  
#  
#    
#    piechart('p53-1.csv','p53-1.png');
#    piechart('p53-2.csv','p53-2.png');
#    piechart('p53-3.csv','p53-3.png');
#    piechart('p53-4.csv','p53-4.png');
#    piechart('p53-5.csv','p53-5.png');
#    piechart('p53-6.csv','p53-6.png');
#    piechart('p53-7.csv','p53-7.png');
#    piechart('p53-8.csv','p53-8.png');
#    piechart('p53-9.csv','p53-9.png');
#    piechart('p53-10.csv','p53-10.png');
#    piechart('p53-11.csv','p53-11.png');
#    piechart('p53-12.csv','p53-12.png');
#    piechart('p53-13.csv','p53-13.png');
#    piechart('p53-14.csv','p53-14.png');
#    piechart('p53-15.csv','p53-15.png');
#    piechart('p53-16.csv','p53-16.png');
#    piechart('p53-17.csv','p53-17.png');
#    piechart('p53-18.csv','p53-18.png');
#    piechart('p53-19.csv','p53-19.png');
#    piechart('p53-20.csv','p53-20.png');
#    piechart('p53-21.csv','p53-21.png');
#    piechart('p53-22.csv','p53-22.png');
#    piechart('p53-23.csv','p53-23.png');
#    piechart('p53-24.csv','p53-24.png');
#    piechart('p53-25.csv','p53-25.png');
#
#
#    
    
     
    sys.stdout.close()
    sys.stdout=stdoutOrigin
    
    
    