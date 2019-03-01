#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 10 10:26:53 2019

@author: yifan
"""

# MutationFinder

''' This version is used for input from cas-offinder validation algorithm'''

from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import re
import pandas
import matplotlib
import matplotlib.pyplot as plt





def smaller_query(up5bp, length_ref, querySeq,moreup):
    ''' returns a shorter query sequence ending at PAM site and is 7 bp longer than reference sequence'''
    myregex=up5bp+"\w+"
    
    
    if re.search(myregex, querySeq) == None:
        myregex=moreup+"\w+"
    #print(myregex)
    #print(querySeq)    
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
                                
                
def already_in_all_mutations_list(sample_mutations_list, all_mut_list, mutation_site,seq_from_PAM):
    exist=False
    for mutation in all_mut_list:
        if mutation[2]==mutation_site:
            if mutation[3]==seq_from_PAM:
                mutation[4]+=1
                sample_mutation=mutation[1:4]
                sample_mutation.append(0) #occurence count
                sample_mutations_list.append(sample_mutation)                                      
                print("This mutation has been already identified as \""+ mutation[1]+ "\" at position "+str(mutation[2])+"in previous samples")
                print("-------------------------------")
                
                return True
    
    return exist
    
def firsttime_all_mutations(casoutput_dir):
    ''' generates a file with all mutations, more mutations can be added later on upon discovery'''
    print("This is the start of comparison")
    #PAM=input("Enter the PAM sequence: ")
    #up5bp=input("Enter 5 bp upstream of sgRNA: ")
    #refseq=input("Enter the reference sequence, including PAM: ")
    
    #for PTEN2:
    PAM="GGG"
    up5bp="TGATT"
    moreup="TGTAA"
    refseq="TGTGCATATTTATTGCATCGGGG"
    
    
    f=open(casoutput_dir,"r") #output file that needs to be parsed (i.e. find the mutations)
    
    
    
    reversed_refseq=reverse(refseq)
    length_ref=len(refseq)
    
    all_lines=f.readlines()
    mutations=[]
    for str_line in all_lines:
        line=str_line.split() 
        queryseq=line[0]
        gene_name=line[1]
        refseqline=line[3]
        seq_from_PAM=get_query_from_PAM(up5bp,reversed_refseq,queryseq,PAM, moreup) #this shorter sequence is reversed
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
                #mutation_name=input("What is the name of this mutation? (e.g. AC in, G sub, etc.) ")
                #mutation_name = 
                
                mutations.append([mutation_name,mutation_site,seq_from_PAM, 0]) # note that seq_from_PAM is a reversed seq
                print(mutations)
                print("We found a new mutation! Next Query:")
        
        
        #writing to CSV file
        df =pandas.DataFrame(mutations, columns=['Mutation Name', 'Mutation Site','Reversed Shortened Mutated Seq', 'Occurence Count -1'])
        df.to_csv('./PTEN2_all_mutations.csv', sep=',')         
#    mut_names=[]
#    mut_sites=[]
#    mut_seq=[]
#    occurrence_count=[]
#    for mutation in mutations:
#        mut_names.append(mutation[0])
#        mut_sites.append(mutation[1])
#        mut_seq.append(mutation[2])
#        occurrence_count.append(mutation[3])
#    df = pandas.DataFrame(data={"Mutation Name":mut_names, "Mutation Site": mut_sites, "Reversed Shortened Mutated Seq": mut_seq, "Occurence Count -1": occurrence_count})
#    df.to_csv('./PTEN2_all_mutations.csv', sep=',')  
#   

     

def mutations_in_a_sample(casoutput_dir, all_mutations_csv_dir,output_csv_dir):
    all_mut_df= pandas.read_csv(all_mutations_csv_dir)
    all_mut=all_mut_df.values.tolist()
    #mut_name: 1; mut_pos: 2; mut_seq:3; mut_count:4
    
    print("This is the start of comparison")
    #PAM=input("Enter the PAM sequence: ")
    #up5bp=input("Enter 5 bp upstream of sgRNA: ")
    #refseq=input("Enter the reference sequence, including PAM: ")
    
    #for PTEN2:
    PAM="GGG"
    up5bp="TGATT"
    moreup="TGTAA"
    refseq="TGTGCATATTTATTGCATCGGGG"
    
    
    f=open(casoutput_dir,"r") #output file that needs to be parsed (i.e. find the mutations)
    
    
    reversed_refseq=reverse(refseq)
    length_ref=len(refseq)
    
    all_lines=f.readlines()
    mutations=[]
    for str_line in all_lines:
        line=str_line.split() 
        queryseq=line[0]
        gene_name=line[1]
        refseqline=line[3]
        seq_from_PAM=get_query_from_PAM(up5bp,reversed_refseq,queryseq,PAM, moreup) #this shorter sequence is reversed
        alignments = pairwise2.align.globalxx(reversed_refseq,seq_from_PAM)
        print("---01234567890123456789")
        print(format_alignment(*alignments[0]))
        mutation_site=find_mismatch_position(reversed_refseq,seq_from_PAM)
        if mutation_site==-1:        
            print("No mutation found :)  Next query:")
            print("    ")
            continue
        else:
            exist_in_sample_mutations_list=already_in_mutations_list(mutations, mutation_site,seq_from_PAM)
            
            if exist_in_sample_mutations_list==True:
                continue
            else:
                exist_in_all_mutations_list=already_in_all_mutations_list(mutations, all_mut, mutation_site, seq_from_PAM)
                if  exist_in_all_mutations_list==True:
                    continue
            
            
                               
            
                print("The first mismatch occurs at position " + str(mutation_site))
                mutation_name=input("What is the name of this mutation? (e.g. AC in, G sub, etc.) ")
                mutations.append([mutation_name,mutation_site,seq_from_PAM, 0]) # note that seq_from_PAM is a reversed seq
                print(mutations)
                print("We found a new mutation! Next Query:")
        
        
    #print(mutations)          
    #writing to CSV file
    df =pandas.DataFrame(mutations, columns=['Mutation Name', 'Mutation Site','Reversed Shortened Mutated Seq', 'Occurence Count -1'])
    
    df.to_csv(output_csv_dir, sep=',') 
    #print(all_mut)
    #update all_mutations csv
    for mut in all_mut:
        del mut[0]
        
    for mymutation in mutations:
        exist=False
        for mut in all_mut:
            if mymutation[2] == mut[2]:
                exist=True            
        if exist==False:        
            all_mut.append(mymutation)
            print("A new mutation has been appended")
            print(mymutation)

        
    updated_df_allmut=pandas.DataFrame(all_mut, columns=['Mutation Name', 'Mutation Site','Reversed Shortened Mutated Seq', 'Occurence Count(ignore this for now)'])
    updated_df_allmut.to_csv(all_mutations_csv_dir, sep=',')
    print('done')




def piechart(csv_dir, txt_dir, piechart_dir):
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
    total = sum(1 for line in open(txt_dir))
    wt = total
    
    
    for entry in filelist:
        # 0:index, 1: name, 2:site, 3:seq, 4: count-1
        entry[4]+=1 #all counts add 1
        labels.append(str(entry[1])+" at site "+str(entry[2]))
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
    print(total)


#def RT_LT_files():
    #for i in range(1, 26):
        #if i%5 == 1:
            
        
    
    
        
if __name__=="__main__":
    
    
    
    #firsttime_all_mutations("/Users/yifan/MiSeq/test/PTEN2-6_output.txt")
   
    #mutations_in_a_sample("/Users/yifan/MiSeq/test/PTEN2-6_output.txt",'/Users/yifan/MiSeq/test/PTEN2_all_mutations.csv','./PTEN2-6.csv')
    #print("All mutations in PTEN2-6 has been identified") 
    #mutations_in_a_sample("/Users/yifan/MiSeq/test/PTEN2-7_output.txt",'/Users/yifan/MiSeq/test/PTEN2_all_mutations.csv','./PTEN2-7.csv')
    #print("All mutations in PTEN2-7 has been identified") 
    #mutations_in_a_sample("/Users/yifan/MiSeq/test/PTEN2-8_output.txt",'/Users/yifan/MiSeq/test/PTEN2_all_mutations.csv','./PTEN2-8.csv')
    #print("All mutations in PTEN2-8 has been identified") 
    #mutations_in_a_sample("/Users/yifan/MiSeq/test/PTEN2-9_output.txt",'/Users/yifan/MiSeq/test/PTEN2_all_mutations.csv','./PTEN2-9.csv')
    #print("All mutations in PTEN2-9 has been identified") 
    #mutations_in_a_sample("/Users/yifan/MiSeq/test/PTEN2-10_output.txt",'/Users/yifan/MiSeq/test/PTEN2_all_mutations.csv','./PTEN2-10.csv')
    #print("All mutations in PTEN2-10 has been identified")    
        
     #piechart('/Users/yifan/MiSeq/test/PTEN2_6-10/csv_output/PTEN2-6.csv', '/Users/yifan/MiSeq/test/PTEN2_6-10/cas_output/PTEN2-6_output.txt','/Users/yifan/MiSeq/test/PTEN2_6-10/piechart_output/PTEN2-6.png');
     piechart('/Users/yifan/MiSeq/test/PTEN2_6-10/csv_output/PTEN2-7.csv', '/Users/yifan/MiSeq/test/PTEN2_6-10/cas_output/PTEN2-7_output.txt','/Users/yifan/MiSeq/test/PTEN2_6-10/piechart_output/PTEN2-7.png');
     piechart('/Users/yifan/MiSeq/test/PTEN2_6-10/csv_output/PTEN2-8.csv', '/Users/yifan/MiSeq/test/PTEN2_6-10/cas_output/PTEN2-8_output.txt','/Users/yifan/MiSeq/test/PTEN2_6-10/piechart_output/PTEN2-8.png');
     piechart('/Users/yifan/MiSeq/test/PTEN2_6-10/csv_output/PTEN2-9.csv', '/Users/yifan/MiSeq/test/PTEN2_6-10/cas_output/PTEN2-9_output.txt','/Users/yifan/MiSeq/test/PTEN2_6-10/piechart_output/PTEN2-9.png');
     piechart('/Users/yifan/MiSeq/test/PTEN2_6-10/csv_output/PTEN2-10.csv', '/Users/yifan/MiSeq/test/PTEN2_6-10/cas_output/PTEN2-10_output.txt','/Users/yifan/MiSeq/test/PTEN2_6-10/piechart_output/PTEN2-10.png');