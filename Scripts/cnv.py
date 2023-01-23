# -*- coding: utf-8 -*-
import os
import sys
import re
import subprocess
import pandas as pd
import argparse



def obtain_cnv(genome_name):
    gen_file = open(genome_name)
    for line in gen_file:
        #we read each genome name, to access their folder and results hmm file
        CNV = 0
        line = line.rstrip('\n')
        hmm_file = f"/mnt/synology/ALEIX/busco_mapping/hmm/{line}/{line}_filtered_output_final.tbl"
        hmm_results = open(hmm_file)
        contig_dic = {}
        counter = 0
        for align in hmm_results:
            #we read each line in the hmm results file and save the results in a dictionary which has as keys the different contigs inside a genome and all the matches in said contig
            align = align.rstrip('\n')
            alignment = align.split(" ")
            if alignment[0] in contig_dic.keys():
                contig_dic[alignment[0]].append(alignment)

            else:
                contig_dic[alignment[0]] = [alignment]
        for contigs in contig_dic.values():
            #we look at each contig list , with its matches and count them
            for region in contigs:
                counter +=1
            #if the number of matches is even, we have to make sure that we have pairs of LSU-SSU, whereas if the counter is odd we have to determine whether the odd number is because it is an incomplete CNV or it is a false positive
            if (counter % 2) == 0:
                LSU = 0
                SSU = 0
                for regions in contigs:
                    if regions[1] == "LSU_rRNA_eukarya":
                        LSU +=1
                    else:
                        SSU +=1

                    if LSU == SSU:
                        CNV += counter/2
                    else: #if the number of LSU and SSU is not the same, it means that a pair is not correct and we must check which ones are not correct and see whether we can count them or not
                        paired = []
                        for region in contigs:
                            #we will compare LSU with SSU to determine the pairs and see which is the odd match
                            if region[1] == "SSU_rRNA_eukarya":
                                continue
                            for region2 in contigs:
                            #if we are comparing the same match or if they're both LSU or if they're not on the same strand, we go to the next comparison
                                if (region == region2) or (region2[1] == "LSU_rRNA_eukarya") or (region2[4] != region[4]):
                                    continue
                                #we check whether the distance between comparison is as expected
                                if region[4] == "+":
                                    distance = int(region[2]) - int(region2[3])
                                    if distance > 400 and distance < 800:
                                        paired.append(region)
                                        paired.append(region2)
                                else:
                                    distance = int(region2[3]) - int(region[2])
                                    if distance > 400 and distance < 800:
                                        paired.append(region)
                                        paired.append(region2)
                        #the match that does not have a pair is the problematic one, so we must see where it lies inside the genome (beginning or end, to see whether it is incomplete)
                        unmatch = [i for i in contigs if i not in paired]
                        for unmatched in unmatch:
                            #if the unmatched genome is found at the very beginning and therefore the ITS region is cut, we will take it into account
                            if (unmatched[4] == "+" and int(unmatched[2]) < 800 and unmatched[1] == "LSU_rRNA_eukarya") or  (unmatched[4] == "-" and int(unmatched[3]) < 800 and unmatched[1] == "SSU_rRNA_eukarya"):
                                CNV += (counter+2)/2
                        #if it is not in the beginning, we will look if it is at the end
                        else:
                    
                            genome_file = f"/mnt/synology/ALEIX/FUNGAL_GENOMES/NCBI/{line}"
                            fd = open(genome_file)
                            seq = "" #variable which we'll use to construct the seq which is divided into various lines
                            seq_dic = {}
                            for line2 in fd: 
                                if (line2[0] == ">"): #if the line we're looking at is the header of a sequence, we're finished with the last seq, therefore we can calculate the frequency of the residue and if the seq passes the threshold, we'll add it to the counter
                                    if len(seq) > 0:
                                        seq_dic[prev_seq] = seq
                                    prev_seq = line2.strip()
                                    seq="" #we restart the new seq
                                else: #if we're still looking at the sequence of the seq we concatenate it to the last part of the seq
                                    seq+=line2.strip()
                            if len(seq) > 0:
                                seq_dic[prev_seq] = seq    
                            fd.close()
                            for unmatched in unmatch:        
                                p = re.compile(unmatched[0])
                                s = [i for i in seq_dic.keys() if p.search(i) != None][0]
                                if (unmatched[4] == "+" and (int(unmatched[3])+2500) > len(seq_dic[s]) and unmatched[1] == "SSU_rRNA_eukarya") or  (unmatched[4] == "-" and (int(unmatched[2])+2500) > len(seq_dic[s]) and unmatched[1] == "LSU_rRNA_eukarya"):
                                    CNV += (counter+2)/2
                                else:
                                    CNV += (counter-2)/2

            else:
                paired = []
                for region in contigs:
                    #we will compare LSU with SSU to determine the pairs and see which is the odd match
                    if region[1] == "SSU_rRNA_eukarya":
                        continue
                    for region2 in contigs:
                        #if we are comparing the same match or if they're both LSU or if they're not on the same strand, we go to the next comparison
                        if (region == region2) or (region2[1] == "LSU_rRNA_eukarya") or (region2[4] != region[4]):
                            continue
                        #we check whether the distance between comparison is as expected
                        if region[4] == "+":
                            distance = int(region[2]) - int(region2[3])
                            if distance > 400 and distance < 800:
                                paired.append(region)
                                paired.append(region2)
                        else:
                            distance = int(region2[3]) - int(region[2])
                            if distance > 400 and distance < 800:
                                paired.append(region)
                                paired.append(region2)
                #the match that does not have a pair is the problematic one, so we must see where it lies inside the genome (beginning or end, to see whether it is incomplete)
                try:
                    unmatch = [i for i in contigs if i not in paired][0]
                except:
                    continue
                #if the unmatched genome is found at the very beginning and therefore the ITS region is cut, we will take it into account
                if (unmatch[4] == "+" and int(unmatch[2]) < 800 and unmatch[1] == "LSU_rRNA_eukarya") or  (unmatch[4] == "-" and int(unmatch[3]) < 800 and unmatch[1] == "SSU_rRNA_eukarya"):
                    CNV += (counter+1)/2
                #if it is not in the beginning, we will look if it is at the end
                else:
                    
                    genome_file = f"/mnt/synology/ALEIX/FUNGAL_GENOMES/NCBI/{line}"
                    fd = open(genome_file)
                    seq = "" #variable which we'll use to construct the seq which is divided into various lines
                    seq_dic = {}
                    for line2 in fd: 
                        if (line2[0] == ">"): #if the line we're looking at is the header of a sequence, we're finished with the last seq, therefore we can calculate the frequency of the residue and if the seq passes the threshold, we'll add it to the counter
                            if len(seq) > 0:
                                seq_dic[prev_seq] = seq
                            prev_seq = line2.strip()
                            seq="" #we restart the new seq
                        else: #if we're still looking at the sequence of the seq we concatenate it to the last part of the seq
                            seq+=line2.strip()
                    if len(seq) > 0:
                        seq_dic[prev_seq] = seq    
                    fd.close()
                                        
                    p = re.compile(unmatch[0])
                    s = [i for i in seq_dic.keys() if p.search(i) != None][0]
                    if (unmatch[4] == "+" and (int(unmatch[3])+2500) > len(seq_dic[s]) and unmatch[1] == "SSU_rRNA_eukarya") or  (unmatch[4] == "-" and (int(unmatch[2])+2500) > len(seq_dic[s]) and unmatch[1] == "LSU_rRNA_eukarya"):
                        CNV += (counter+1)/2
                    else:
                        CNV += (counter-1)/2
            counter = 0
        print(CNV)
        if int(CNV) == 0:
            continue
        #now we must create a bed file in order to extract the ITS sequence from the whole genome file
        extract_seq = []
        for contigs in contig_dic.values():
            for region in contigs:
            #we will compare LSU with SSU to determine the pairs and see which is the odd match
                if region[1] == "SSU_rRNA_eukarya":
                    continue
                for region2 in contigs:
                #if we are comparing the same match or if they're both LSU or if they're not on the same strand, we go to the next comparison
                    if (region == region2) or (region2[1] == "LSU_rRNA_eukarya") or (region2[4] != region[4]):
                        continue
                    #we check whether the distance between comparison is as expected
                    if region[4] == "+":
                        distance = int(region[2]) - int(region2[3])
                        if distance > 300 and distance < 800:
                            extract_seq = [region[0],region2[3],region[2]]
                            print(extract_seq)
                            break
                    else:
                        distance = int(region2[3]) - int(region[2])
                        if distance > 300 and distance < 800:
                            extract_seq = [region[0],region[2],region2[3]]
                            print(extract_seq)
                            break
        outf = open(f"/mnt/synology/ALEIX/busco_mapping/hmm/{line}/{line}.bed","w")
        print(extract_seq)
        print(*extract_seq,sep="\t",file=outf)
        outf.close()
        seqfile = open(f"/mnt/synology/ALEIX/busco_mapping/hmm/{line}/{line}_itsseq.txt","w+")
        subprocess.call(["bedtools","getfasta","-fi",f"/mnt/synology/ALEIX/FUNGAL_GENOMES/NCBI/{line}","-bed",f"/mnt/synology/ALEIX/busco_mapping/hmm/{line}/{line}.bed"], stdout=seqfile)
        seqfile.close()
        itsfile = open(f"/mnt/synology/ALEIX/busco_mapping/hmm/{line}/{line}_itsseq.txt","r")
        for its in itsfile:
            if its.startswith(">"):
                continue
            else:
                its_seq = its.rstrip()
        #seqfile.close()
        m = re.compile("\.[1|2|3|4]")
        x = m.search(line)
        #try:
        #    line2 = line[0:x.end()]
        #
        #except:
        #    line2 = line
        seqfile.close()
        n = re.compile("/")
        y = n.split(genome_name)
        final = open(f"ITS_CNV_DDBB_{y[-1]}.tsv","a+")
        try:
            final_table = [line,int(CNV),its_seq]
            print(*final_table,sep="\t",file=final)
        except:
            continue


if __name__=="__main__":
    obtain_cnv(sys.argv[1])