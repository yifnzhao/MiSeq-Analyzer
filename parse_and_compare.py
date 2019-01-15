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
        
def compare(myList):
    """
    this function compares a list of reads to reference reads
    
    output: tbd
    """
    

    
    
#__main__
mydictIDSeq = parse("MI.M03555_0362.001.FLD0004.PTEN2-12_R1.fastq")
myListSeq = dict2list(mydictIDSeq)
