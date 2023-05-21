#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 12:27:26 2021

@author: mobilab5
"""
import csv
import sys
import os
from os import listdir
from os.path import isfile, join
import re
import subprocess

def cigar_substr_list(cigar):
    MDI_list = ["M", "D", "I"]
    cigar_substr = list()
    pos_list = list()
    
    for MDI in MDI_list:
        for index in re.finditer(MDI,cigar):
            pos_list.append(index.end(0))
    pos_list.sort()
    
    cigar_substr.append(cigar[:pos_list[0]])    # first element    
    for i in range(len(pos_list)-1):            # following elements
        start = pos_list[i]
        end = pos_list[i+1]
        cigar_substr.append(cigar[start:end])
    # print(cigar_substr)    
    return cigar_substr

def SeqBasedCigar(ori_seq, cigar_list, ori_qual, start_pos):
    # print(cigar_list)
    if ("1I" in cigar_list) or ("1D" in cigar_list):
        seq = ori_seq
        qual = ori_qual
        # A. fixed cigar_list with 1D
        count_d = 0
        for i in range(len(cigar_list)):   
            if cigar_list[i] == "1D":    
                count_d += 1
        # print("original_len:",len(seq))
        # print("count_d",count_d)
        
        for i in range(count_d):
            # 1. Separate cigar list to num_list & char_list
            num_list = [int(re.findall(r'\d+', i)[0]) for i in cigar_list]      # Get only num part in cigar_list
            char_list = [i[-1] for i in cigar_list]
            
            # 2. Find left-most 1D index
            j = cigar_list.index("1D")
            
            # 3. Fixed sequence
            alt_pos = sum([num_list[j] for j in range(len(char_list[:j+1])) if char_list[j] == "M" or char_list[j] == "I"]) 
            ref_pos = sum([num_list[j] for j in range(len(char_list[:j+1])) if char_list[j] == "M" or char_list[j] == "D"]) + start_pos -2
            temp_seq = list(seq)
            temp_seq.insert(alt_pos,ref[ref_pos])
            seq =  "".join(temp_seq)
            # print("temp_seq_len",len(seq))
            
            # 4. Fixed qual
            temp_qual = list(qual)
            temp_qual.insert(alt_pos,"G")
            qual =  "".join(temp_qual)
            
            # 5. Fixed cigar
            alt_cigar = str(num_list[j-1] + num_list[j+1] + 1) +"M"
            del cigar_list[j-1:j+2]             
            cigar_list.insert(j-1,alt_cigar)
        # print(cigar_list)     
        
        # B. fixed cigar_list with 1I   
        count_i = 0
        for i in range(len(cigar_list)):   
            if cigar_list[i] == "1I":    
                count_i += 1
        
        # print("count_i",count_i) 
        for i in range(count_i):
            # 1. Separate cigar list to num_list & char_list
            num_list = [int(re.findall(r'\d+', i)[0]) for i in cigar_list]      # Get only num part in cigar_list
            char_list = [i[-1] for i in cigar_list]
            
            # 2. Find left-most 1I index
            j = cigar_list.index("1I")
            
            # 3. Fixed sequence
            alt_pos = sum([num_list[j] for j in range(len(char_list[:j+1])) if char_list[j] == "M" or char_list[j] == "I"])            
            temp_seq = list(seq)
            del temp_seq[alt_pos]
            seq =  "".join(temp_seq)
            # print("temp_seq_len",len(seq))
            
            # 4. Fixed qual
            temp_qual = list(qual)
            del temp_qual[alt_pos]
            qual =  "".join(temp_qual)
            
            # 5. Fixed cigar
            alt_cigar = str(num_list[j-1] + num_list[j+1]) +"M"
            del cigar_list[j-1:j+2]             
            cigar_list.insert(j-1,alt_cigar)
        # print(cigar_list)
        # Final new_cigar
        new_seq = seq
        new_qual = qual
        new_cigar = "".join(cigar_list)
        # print("new_len:",len(new_seq))
        # print(new_cigar)
    
    else:
        new_seq = ori_seq
        new_qual = ori_qual
        new_cigar = "".join(cigar_list)
        
        # print("new_len:",len(new_seq))
        # print(new_cigar)
    return new_seq,new_qual,new_cigar


# RUN: create "$sample_fixed_ngs_error.sam" file
sample = sys.argv[1]
wd_path = sys.argv[2] + "/" + sample + "/2_mapping"
os.chdir(wd_path)
sw_path = sys.argv[4] + "/script" 
sys.path.append(sw_path)
from ExtendedHiFiBR import SeqInfo as si

ref = si.get_refSeq(sys.argv[3])
reader = csv.reader(open(sample+'_FilteredAligned.sam',"r"), dialect="excel-tab")

for line in reader:
    if line[0][0] == "@":
        pass    
    else:          
        cigar_list = cigar_substr_list(line[5])
        char_list = [i[-1] for i in cigar_list]
        if (char_list[0] == "M") and (char_list[-1] == "M"):
            # A. Check if there is only 1D/1I in cigar
            if ("I" in char_list) or ("D" in char_list):
                indel_cigar_i = [i for i,s in enumerate(char_list) if (s == "I") or (s == "D")]
                indel_cigar_L = [cigar_list[i] for i in indel_cigar_i]
                indel_num_L = [i[:-1] for i in indel_cigar_L]
                cond = all([i == "1" for i in indel_num_L])
                #print(indel_cigar_L)
                #print(cond)
                
            if cond:
                start_pos = int(line[3].split(":")[-1])
                new_res = SeqBasedCigar(line[9], cigar_list, line[10], start_pos) # line[9]: seq; line[10]: qual; line[5]: cigar
                new_seq = new_res[0]
                new_qual = new_res[1]
                new_cigar = new_res[2]
                line_list = [line[0], line[1], line[2], line[3], line[4], new_cigar, 
                             line[6], line[7], line[8], new_seq, new_qual, line[11]]
            else:
                line_list = [line[0], line[1], line[2], line[3], line[4], line[5], 
                             line[6], line[7], line[8], line[9], line[10], line[11]]
            
            # B. Write file    
            with open(sample+"_NoNgsError.sam", "a+") as padded_seq:
                padded_seq.write("\t".join(line_list))
                padded_seq.write('\n')
                
        else: 
            with open(sample+"_Undealt_NGSerror.txt","a+") as undealt:
                undealt.write(line[0]+"\n")
    # with open(sample+"_NoNgsError.sam", "a+") as padded_seq:
    #     padded_seq.write('\n')

cmd = "samtools view -H " + sample + '_FilteredAligned.sam'
proc=subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, )
header=proc.communicate()[0].decode("utf-8")

f = open(sample+"_NoNgsError.sam",'r+')
lines = f.readlines()   # read old content
f.seek(0)               # go back to the beginning of the file
f.write(header)         # write new content at the beginning
for line in lines:      # write old content after new
    if not line.isspace():  # remove blank lines
        f.write(line)
f.close()


# RUN2: fixed error line in sam file. 
# NOT pattern: 180M-1D-89M (D/I between 2 M) -> Ex: 528M-37D-1I-2M-4D-437M
# Result: 528M-39M-4D-437M (2 M adjacent)


del_line = 0
while True:
    try:
        cmd = "samtools view -c " + sample + "_NoNgsError.sam"
        proc=subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        error=proc.communicate()[1].decode("utf-8")
    
        del_info = error.split("\n")[1]
        del_num = [int(s) for s in del_info.split() if s.isdigit()][-1]
        #print(del_num)
        cmd = 'sed -i "' + str(del_num) + 'd" ' + sample + "_NoNgsError.sam"
        os.system(cmd)
        del_line += 1
    except IndexError:
        pass
    if "different length" not in error:
        break
    
print("--- REPORT --- \nDeleted lines " + sample, str(del_line))    
