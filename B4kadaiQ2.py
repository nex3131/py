filename = "hs_ref_GRCh38_chr21.fa" ##input("fastafilename:")
gff_file = "ref_GRCh38_scaffolds.gff3" ##input("gfffilename:")
#LOC102724951geneid =  "LOC102724951" 
geneid = input("geneqid:")


#reading fastafilewhih
from Bio import SeqIO
def read_fasta_file(filename):
    sequences = []
    for record in SeqIO.parse(filename, "fasta"):
        sequences.append(record)
    return sequences
sequences = read_fasta_file(filename)


array = []

#reading gfffile
with open(gff_file) as file:
#loop for every line in file
    for line in file:
        if line[0] != "#":
            if geneid in line:
                array.append(line.split('\t'))
file.close()



parent = "gene="+geneid
hairetsu =[]
rhairetsu =[]
j = -1
k = -1

#screening sequences
for seq in sequences:
    for i in range(len(array)):
        if array[i][0] in seq.id:
            if parent in array[i][8]:
                if array[i][2] == "CDS":            
                    a = int(array[i][3])
                    b = int(array[i][4])
                    
                    
                    if array[i][6] == "-":
                        rseq = seq.seq[a-1:b].reverse_complement()
                        #print(rseq)
                        j += 1
                        rhairetsu.append([a,b,rseq])
                        rhairetsu[j].append(array[i][8].split(";"))
                        
                        """
                        c = len(rseq)
                        base = -1
                        while base < c-1:
                            base +=1
                            if rseq[base] == "A":
                                print("T",end = "")
                            if rseq[base] == "T":
                                print("A",end = "")
                            if rseq[base] == "G":
                                print("C",end = "")
                            if rseq[base] == "C":
                                print("G",end = "")
                        """
                    else :
                        #print(seq.seq[a-1:b], end=" ")
                        k += 1
                        hairetsu.append([a,b,seq.seq[a-1:b]])
                        hairetsu[k].append(array[i][8].split(";"))
                        #print(seq)
                        #print(seq.translate())





"""
#output merging sequences
mrna = ""
for seq in sequences:
    for i in range(len(array)):
        if array[i][0] in seq.id:
            if parent in array[i][8]:
                if array[i][2] == "CDS":            
                    a = int(array[i][3])
                    b = int(array[i][4])
                    
                    if array[i][6] == "-":
                        rseq = seq.seq[a-1:b][::-1]
                        
                        print(rseq.translate())
                        
                        c = len(rseq)
                        base = -1
                        while base < c-1:
                            base +=1
                            if rseq[base] == "A":
                                print("T",end = "")
                            if rseq[base] == "T":
                                print("A",end = "")
                            if rseq[base] == "G":
                                print("C",end = "")
                            if rseq[base] == "C":
                                print("G",end = "")
                        
                    else :
                        #print(seq.seq[a-1:b], end=" ")
                        mrna += seq.seq[a-1:b]
                        #print(seq)
                        #print(seq.translate())
                        
"""                
###print (hairetsu)   
sortarray=(sorted(hairetsu, key=lambda x: x[0]))
rsortarray=(sorted(rhairetsu, key=lambda x: x[0], reverse=True))
sortarray=(sorted(sortarray, key=lambda x: x[3][0]))
rsortarray=(sorted(rsortarray, key=lambda x: x[3][0]))


#print(sortarray)
#print (len(sortarray[0][2]))


import itertools

cds = []
for k , g in itertools.groupby(sortarray,lambda x: x[3][0]):
    cds = list(g)
    q = 0
    amino = ""
    #print(k,len(sss))
    while not q == len(cds):
        #print (k,sss[q][2])
        amino += cds[q][2]
        q +=1
    #print(k,amino)
    print(k,amino.translate())
    

import itertools
cds = []
for k , g in itertools.groupby(rsortarray,lambda x: x[3][0]):
    cds = list(g)
    q = 0
    ramino = ""
    #print(k,len(sss))
    while not q == len(cds):
        #print (k,sss[q][2])
        ramino += cds[q][2]
        q +=1
    #print(k,ramino)
    print(k,ramino.translate())