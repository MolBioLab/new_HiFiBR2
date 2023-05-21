#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  7 23:39:43 2022

@author: mobilab5
"""
import os, sys
import pysam 
import subprocess

# FUNCTION
def getCigar_info(cigar_T):
    cigar_char = []
    cigar_num = []
    for charID, num in cigar_T:
        cigar_char.append(cigar_info[charID])
        cigar_num.append(num)
    return cigar_char, cigar_num

def sum_SpecificCigar(interestedRange, interestedChar, cigar_char, cigar_num):
    cigar_numSpecific_L = []
    for i in interestedRange: 
        if cigar_char[i] in interestedChar:
            cigar_numSpecific_L.append(cigar_num[i])
    cigar_sumSpecific = sum(cigar_numSpecific_L)
    return cigar_sumSpecific

# RUN
# path = "/mnt/data/TAnh_K17/nhej_long/exp/02_ValidateTool/res/simu03/complex-uncut-1-side"
# sample = "MN00778:127:000H55J52:1:11101:2542:164"
# wd_path = path + "/" + sample + "/2_mapping"
# os.chdir(wd_path)
# sw_path = "/mnt/data/TAnh_K17/nhej_long/lib/new_HiFiBR/script" 
# sys.path.append(sw_path)

sample = sys.argv[1]
wd_path = sys.argv[2] + "/" + sample + "/2_mapping"
os.chdir(wd_path)
sw_path = sys.argv[3] + "/script" 
sys.path.append(sw_path)

cigar_info = {0:"M", 1:"I", 2:"D", 3:"N",
              4:"S", 5:"H", 6:"P", 7:"=",
              8:"X", 9:"B"}
    
for align_i in [1, 2]:
    insam = f"{sample}_secondaryAlignment{align_i}.sam"
    outsam = f"{sample}_fixed_a_align{align_i}.sam"
    
    # A. Grouping uniq & multiple alignment reads
    cmd = f'./step3_SecondaryAlignment.sh {sample} {align_i} {wd_path}'
    process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    
    
    aligned_info = {}
    samfile = pysam.AlignmentFile(insam, "rb")
        # Each read get info: {key:[startPos, endPos, cigar]}
    for read in samfile.fetch():
        info_L = [read.flag,
                  read.reference_name,
                  read.reference_start,
                  read.mapping_quality,
                  read.cigartuples,
                  "*",
                  0,
                  0,
                  read.seq,
                  read.qual,
                  "RG:Z:Unpaired_reads_assembled_against_pCOH_CD4_alt_PCR",
                  read.reference_end]
        
        if read.query_name not in aligned_info: 
            aligned_info[read.query_name] = [info_L]
        else:
            aligned_info[read.query_name].append(info_L)
    
    for query_name, info in aligned_info.items():
        if len(info) == 2:
            startPos1, endPos1, cigar_T1 = [info[0][i] for i in [2, -1, 4]]
            startPos2, endPos2, cigar_T2 = [info[1][i] for i in [2, -1, 4]]
            
            if startPos1 < startPos2:
                start_l, end_l, cigar_Tl = startPos1, endPos1, cigar_T1
                start_r, end_r, cigar_Tr = startPos2, endPos2, cigar_T2
            else:
                start_l, end_l, cigar_Tl = startPos2, endPos2, cigar_T2
                start_r, end_r, cigar_Tr = startPos1, endPos1, cigar_T1
            
            # Only execute alignment pair with pattern
            cigar_char_l, cigar_num_l = getCigar_info(cigar_Tl)
            cigar_char_r, cigar_num_r = getCigar_info(cigar_Tr)
            
            # Check-point 1 
            # Ex:   151M 1D 10M 3S -> Sum all M & I = 151 + 10 = 161 (same as below 161S)
            #       161S 3M
            range_l = range(len(cigar_num_l)-1) # exclude left-side soft-clipped (last position)
            cigar_sumMI_l = sum_SpecificCigar(range_l, ["M","I"], cigar_char_l, cigar_num_l)
            cond1 = cigar_sumMI_l == cigar_num_r[0]
            
            # Check-point 2
            # Ex:   151M 1D 10M 3S
            #       161S 3M -> Sum all M & I = 3 (same as above 3S)
            range_r = range(1,len(cigar_num_r)) # exclude right-side soft-clipped (1st position)
            cigar_sumMI_r = sum_SpecificCigar(range_r, ["M","I"], cigar_char_r, cigar_num_r)
            cond2 = cigar_sumMI_r == cigar_num_l[-1]
            
            if cond1 and cond2:
                # Get all info from left-side, exclude soft-clipped from right-side (1st position)
                new_cigar_l = "".join([str(cigar_num_l[i]) + cigar_char_l[i] for i in range(len(cigar_num_l)-1)])
     
                # Get all info from right-side, exclude soft-clipped from left-side (last position)
                new_cigar_r = "".join([str(cigar_num_r[i]) + cigar_char_r[i] for i in range(1,len(cigar_num_r))]) 
                    
                # Del length
                range_l = range(len(cigar_num_l))
                cigar_sumMD_l = sum_SpecificCigar(range_l, ["M","D"], cigar_char_l, cigar_num_l)
                
                range_r = range(len(cigar_num_r))
                cigar_sumMD_r = sum_SpecificCigar(range_r, ["M","D"], cigar_char_r, cigar_num_r)
                
                deletion = end_r - cigar_sumMD_l - cigar_sumMD_r - start_l 
                new_cigar = new_cigar_l + str(deletion) + "D" + new_cigar_r    

                info[0][2] = start_l + 1
                info[0][4] = new_cigar
                new_line = [query_name] + list(map(str,info[0][:-1]))
                new_line = "\t".join(new_line) # exclude endPos at the end
                
                with open(outsam,"+a") as out:
                    out.write("\n" + new_line)
            
            else:
                with open("log.txt", "+a") as log:
                    log.write(query_name)
                
    cmd = f"sed -i '/^[[:space:]]*$/d' {outsam}"
    os.system(cmd)
    # cmd = f"rm {insam}"
    # os.system(cmd)
