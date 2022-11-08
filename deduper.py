#!/usr/bin/env python
from lib2to3.pgen2.token import MINUS
import math
import argparse
import re

#Functions

def get_args():
    parser = argparse.ArgumentParser(description="specifies file inputs, index cutoff")
    parser.add_argument("-u", "--umi", help="file of known UMIs", \
        required=True, type=str)
    parser.add_argument("-f", "--file", help="input: sorted SAM file", \
        required=True, type=str)
    parser.add_argument("-o", "--outfile", help="deduplicated and sorted SAM file", \
        required=True, type=str)
    return parser.parse_args()


#create short variable names from argparse inputs to reference later in script
clinput = get_args()
f = clinput.file
u = clinput.umi
o = clinput.outfile


#initialize empty set to hold list of known UMIs
umi_set = set()
i = 0
#loop through umi file and stores each umi as an object in a list
with open(u, "r") as umis: 
    for line in umis:
        i +=1 
        umi = line.strip("\n")
        umi_set.add(umi) #
    #print(umi_set)


#record-grab gets the current record and needed variables using file handle input
def record_grab(line):
    #takes a line in a sam file, splits it by tab into fields, saves umi, flag, rname, position, cigar ainto variables and returns them as a list. gets umi from splitting the qname by colon
    rec = line.strip('\n') #strip new line character
    #print(line)
    lxl = rec.split("\t") #split line into fields by tab
    qname = lxl[0].split(":")
    umi = str(qname[-1])
    flag = int(lxl[1])
    rname = str(lxl[2])
    pos = int(lxl[3])
    cigar = str(lxl[5])
    #rec_list = umi, flag, rname, pos, cigar
    # print(rec,umi,flag, rname, pos, cigar)
    return umi, flag, rname, pos, cigar


def strander(flag):
#takes an interger variable (flag) and uses bitwise flag code magic to determine if the flag corrsponds to a forward or reverse strand
    if((flag & 16) == 16):
        strand = "reverse"
    else:
        strand = "forward" 
    return strand
        
#print(strander(0)) 
# should be forward
#print(strander(83))
# should be reverse

def cig_parse(pos,cigar,strand):
    #takes left most starting position from sam file and readjusts to the actual starting position using the cigar string and strand
    if strand == "forward":
        S = re.findall(r'^(\d+)S', cigar)    
        if S: S = int(S[0])                 
        if S: pos -= S 
    elif strand == "reverse":
        S = re.findall(r'(\d+)S$', cigar)    # Find all S
        if S: S = int(S[1])                 # Convert it to int and also drop any other S occurances. RIGHT clip.
        if S: pos += S                 # Add S to our start position 
        M = re.findall(r'(\d+)M', cigar)    # Find all M
        if M: M = [int(i) for i in M]       # Convert list.str to list.int
        if M: pos += sum(M)            # Add all Ms to position
        D = re.findall(r'(\d+)D', cigar)    # Find all D
        if D: D = [int(i) for i in D]             
        if D: pos += sum(D)  
        N = re.findall(r'(\d+)N', cigar)    # Find all D
        if N: N = [int(i) for i in N]             
        if N: pos += sum(N)  
    return pos

#test of forward clip-adjust portion (it works)
# print(cig_parse(40,"10S40M10S","forward"))
#print((cig_parse(40,"40M10S","forward")))
# print(cig_parse(40,"40M","forward"))
#test of reverse (it works)
# print(cig_parse(10,"10S10D10M10N20M10S", "reverse"))
# print(cig_parse(10,"10S10M10N20M10S", "reverse"))
# print(cig_parse(10,"10S10D10M20M10S", "reverse"))
# print(cig_parse(10,"10S10D10M10N20M", "reverse"))

unk_umi_set = set()
check_set = set()

with open(o, "w") as fo:
    with open(f, "r") as fh:
        for line in fh:
            rec = line.strip("\n")
            if rec.startswith("@"):
                #print(rec)
                fo.write(line)
            else:
                umi, flag, rname, pos, cigar = record_grab(rec)
                # print(umi, flag, rname, pos, cigar)
                # if umi not in umi_set:
                #     continue
                chrom = rname
                strand = strander(flag)
                pos = (cig_parse(pos, strand, cigar))
                # print(pos)
                if umi in umi_set:
                    check_list = (umi, chrom, strand, pos)
                    # print(check_list)
                    if check_list in check_set:
                        continue
                    else:
                        check_set.add(check_list)
                        fo.write(line)
                else:
                    unk_umi_set.add(umi)
fo.close()
# print(len(unk_umi_set))
# print(len(dup_dict))    
           

            

        


                










    