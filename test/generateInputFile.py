'''
This module generates an input file for cas-offinder
'''

def parse(fileName):
    """
    this function parses the file and obtain the reads as a dictionary of {id:seq}
    input: fileName (name of the file that you want to parse)
    output: a dictionary whose keys are ids of reads and values are the correpsonding sequences    
    """
    from Bio import SeqIO
    dictIDSeq={}
    for record in SeqIO.parse(fileName, "fastq"):
        print("%s %s" % (record, record.seq))
        #print("%s" % (record.seq))
        dictIDSeq[record.id]=record.seq
    #print (dictIDSeq)
    #for key, value in dictIDSeq.items():
        #print(len(value)) about 250?
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

if __name__=="__main__":
    writeInput("/Users/yifan/MiSeq/test/MI.M03555_0362.001.FLD0083.PTEN2-10_R1.fastq")
    print('The input file has been created.')
