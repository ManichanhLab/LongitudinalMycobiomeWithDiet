# -*- coding: utf-8 -*-
import os
import sys    

def concatenate(genomes):

    fd = open(genomes)
    final = open("gut_fungi_genomes.fasta","a+")
    for genom in fd:
        fd2 = open(genom.strip())
        seq = "" #variable which we'll use to construct the seq which is divided into various lines
        seq += ">%s"%(genom)
        print(seq)
        for line2 in fd2: 
            if (line2[0] == ">"):
                continue
            else: #if we're still looking at the sequence of the seq we concatenate it to the last part of the seq
                seq+=line2.strip()
        fd2.close()        
        print(seq,file=final)
          
    fd.close()
    final.close()


if __name__ == "__main__":
    concatenate("gut_fungi_all.txt")
