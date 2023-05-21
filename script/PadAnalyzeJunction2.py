#!/usr/bin/env python

import csv
import sys
import os
from os import listdir
from os.path import isfile, join

def cigar_substring_list(cigar_string):
    MDI_list = ["M", "D", "I"]
    cigar_substring_list = []
    cigar_len = len(cigar_string)
    i = 0
    while i < cigar_len:
        sub_cigar = ""
        while cigar_string[i] not in MDI_list:
            sub_cigar = sub_cigar + cigar_string[i]
            i += 1
        sub_cigar = sub_cigar + cigar_string[i]
        cigar_substring_list.append(sub_cigar)
        i += 1
    return cigar_substring_list

# RUN
sample = sys.argv[1]
wd_path = sys.argv[2] + "/" + sample + "/2_mapping"
os.chdir(wd_path) 
sw_path = sys.argv[4] + "/script" 
sys.path.append(sw_path)
from ExtendedHiFiBR import SeqInfo as si

ref = si.get_refSeq(sys.argv[3])     
 
try:
    reader = csv.reader(open(sample+'_AutoAnalyse2.sam',"r"), dialect="excel-tab")
    cmd = f"samtools view -H {sample}_AutoAnalyse2.sam > {sample}_Padded.sam"
    os.system(cmd)
    
except:
    reader = csv.reader(open(sample+'_AutoAnalyse.sam',"r"), dialect="excel-tab")
    cmd = f"samtools view -H {sample}_AutoAnalyse.sam > {sample}_Padded.sam"
    os.system(cmd)
    
padded_seq = open(f"{sample}_Padded.sam", "a+")

for line in reader:
    if line[0][0] == "@":
        pass
    else:
        cigar_list = cigar_substring_list(line[5])
        char_cigar = [i[-1] for i in cigar_list]
        if (char_cigar[0] != "M") or (char_cigar[-1] != "M"):
            continue
        # A. Change seq & qual
        left_index = int(line[3]) - 1    # left-most start position
        
        for j in line:
            if j[:2] == "YS": right_index = int(j.split(":")[-1]) + 1  # right-most position
        quality = line[10]    # qual score
        seq = line[9]         # seq
        
        l_pad = ref[:left_index]
        r_pad = ref[right_index-1:]
        new_seq = l_pad + seq + r_pad
        new_quality = "G" *(len(l_pad)) + quality + "G"*(len(r_pad))
        
        # B. Change cigar        
        
        if len(cigar_list) == 1:
            new_cigar = str(len(l_pad) + int(cigar_list[0][:-1]) + len(r_pad)) + "M"
        else:
            new_M_l = str(len(l_pad) + int(cigar_list[0][:-1])) + "M"
            new_M_r = str(len(r_pad) + int(cigar_list[-1][:-1])) + "M"
            new_cigar = new_M_l + "".join(cigar_list[1:-1]) + new_M_r
        
        # C. Add Tag
        old_cigar = "OC:Z:" + line[5]
        old_seq = "OS:Z:" + line[9]
        
        # D. Change pos -> new_line                
        line_list = [line[0], line[1], line[2], "1", line[4], new_cigar, line[6],
                     line[7], line[8], new_seq, new_quality, line[-1], old_cigar, old_seq]
        line = "\t".join(line_list)
        padded_seq.write(line + "\n")
        
