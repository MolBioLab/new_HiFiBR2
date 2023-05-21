#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 23 23:04:12 2023

@author: mobilab5
"""

import sys
import os

# DEF
def check_match(sub_readSeq, sub_readQual, sub_refSeq):
    mismatch = 0
    for i in range(len(sub_readSeq)):
        # Check if any mismatch (with Q >=30) between left side of read and ref
        # Expect mismatch = 0
        if sub_readSeq[i] != sub_refSeq[i]:
            # Check if mismtach is actually mismatch (Q>30)
            qual_num = ord(sub_readQual[i]) - 33
            
            if qual_num >= 30:
                mismatch += 1
             
    if mismatch == 0:        
        return True
    else:
        return False  
    
    return


# RUN
sample = sys.argv[1]
wd_path = sys.argv[2]
os.chdir(f"{wd_path}/{sample}/2_mapping")
sw_path = sys.argv[4] + "/script" 
sys.path.append(sw_path)

from ExtendedHiFiBR import SeqInfo as si
refSeq = si.get_refSeq(sys.argv[3])

cutPos1 = int(sys.argv[5])
cutPos2 = int(sys.argv[6])

full_info = []
print(sample)
# Filter seq with 2 side totally match
with open(f"{sample}_len.sam", "r") as sam:
    lines = sam.readlines()
    
    for line in lines:
        line = line.strip()
        if line[0] != "@": # Get info -> Filter
            full_info.append(line)
            
        else: # Write sam header
            with open(f"{sample}_AutoAnalyse.sam", "+a") as out:
                out.write(line + "\n")
    
for info in full_info:
    info_L = info.split('\t')
    startPos = int(info_L[3])-1
    
    for j in info_L:
        if j[:2] == "YS": endPos = int(j.split(":")[-1])
        
    cigar = info_L[5]
    readSeq = info_L[9]
    qual = info_L[10]
    
    # Check if left-side, right-side totally match with ref
    lSide_ref = refSeq[startPos : startPos + 11] 
    rSide_ref = refSeq[endPos - 11 : endPos]
    
    l_check = check_match(readSeq[:11], qual[:11], lSide_ref)
    r_check = check_match(readSeq[-11:], qual[-11:], rSide_ref)
    
    # Check if left_side actually left of cutPos1 & right_side actually right of cutPos2
    check1 = startPos + 10 < cutPos1
    check2 = endPos - 10 > cutPos2
    
    # Filter
    if l_check == True and r_check == True and check1 == True and check2 == True:
        with open(f"{sample}_AutoAnalyse.sam", "+a") as out:
            out.write(info + "\n")