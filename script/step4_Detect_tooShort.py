#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 19:17:16 2023

@author: mobilab5
"""

import os 
import difflib
import sys
import pandas as pd
import numpy as np
import csv
from tqdm import tqdm
# FUNCTION
def get_info(f_name):
    name_L, cigar_L, seq_L, start_pos, end_pos, qual_L = [], [], [], [], [], []
    reader = csv.reader(open(f_name,"r"), dialect="excel-tab") # f_name+'.sam'
    repre_info = []
    # ["auto1_name","auto1_cigar","auto1_seq","auto1_start","auto1_end","auto1_qual"]
    for line in reader:
        if line[0][0] == "@":
            pass    
        else:
            if "N" in line[9]:
                continue
            
            for j in line:
                if j[:2] == "YS": right_index = int(j.split(":")[-1])
            
            left_index = int(line[3])-1
            
            info_L = [line[0],
                      line[5],
                      line[9],
                      left_index,
                      right_index,
                      line[10]]
            repre_info.append(info_L)
    return repre_info

def get_overlap(s1, s2):
    s = difflib.SequenceMatcher(None, s1, s2, autojunk=False)
    pos_a, pos_b, size = s.find_longest_match(0, len(s1), 0, len(s2)) 
    
    s1_start, s1_end = pos_a, pos_a+size
    s2_start, s2_end = pos_b, pos_b+size
    
    if s1_start == 0 and s1_end == 0: # s1 embed in s2
        return True
    elif s1_end == len(s1) and s2_start == 0: # s1-right overlap s2-left
        return True
    elif s1_start == 0 and s2_end == len(s2): # s2-right overlap s1-left
        return True
    else:
        return False

# RUN
# A. INPUT
sample = sys.argv[1]
wd_path = sys.argv[2]
os.chdir(f"{wd_path}/3_HiFiBR") 

sw_path = sys.argv[4] + "/script" 
sys.path.append(sw_path)

from ExtendedHiFiBR import SeqInfo as si
refSeq = si.get_refSeq(sys.argv[3]) 

df = pd.read_excel(f"{sample}_Padded_Final_temp2.xlsx", engine='openpyxl') # f"{sample}_padded_vs_altref_Final2.xlsx"
sub_df = df[["Len of matching left", "Len of matching right", "Represented read"]]
sub_D = sub_df.set_index('Represented read').T.to_dict('list') # {Represented read : [Len of matching left, Len of matching right]}

# B. GET reanalyse sam file
reanalyse_L = df["Represented read"].to_list()
with open("reanalyse.txt","a+") as f:
    f.write("\n".join(reanalyse_L))

cmd = f""" samtools view -H {wd_path}/2_mapping/{sample}_AutoAnalyse.sam > {wd_path}/3_HiFiBR/{sample}_reanalyse_Auto.sam
samtools view {wd_path}/2_mapping/{sample}_AutoAnalyse.sam | rg -f reanalyse.txt >> {wd_path}/3_HiFiBR/{sample}_reanalyse_Auto.sam
"""
os.system(cmd)

# C. ANALYZE maybe same event between junctions
auto1_repre = auto2_repre = get_info(f"{sample}_reanalyse_Auto.sam")
check_all_D = {} # record all overlapped-pair


for in1 in tqdm ( range(len(auto1_repre)) , desc=" outer", position=0):
    auto1_name,auto1_cigar,auto1_seq,auto1_start,auto1_end,auto1_qual = auto1_repre[in1]
    
    for in2 in tqdm( range(in1,len(auto2_repre)) , desc=" inner loop", position=1, leave=False):
        check_L = []
        auto2_name,auto2_cigar,auto2_seq,auto2_start,auto2_end,auto2_qual = auto2_repre[in2]
        
        # 1. Find overlapped-pair and record its overlap
        if auto1_name != auto2_name:            
            check = get_overlap(auto1_seq, auto2_seq)
            
            if check == False:
                continue
            elif check == True:
        
                # 2. Find which junction is "Maybe same event" with others
                auto1_Ml = int(sub_D[auto1_name][0])
                auto2_Ml = int(sub_D[auto2_name][0])
                
                auto1_Mr = int(sub_D[auto1_name][1])
                auto2_Mr = int(sub_D[auto2_name][1])
                
                delta_left = auto1_Ml - auto2_Ml
                delta_right = auto1_Mr - auto2_Mr
                l_gt_r = abs(delta_left) > abs(delta_right)
                r_gt_l = abs(delta_right) > abs(delta_left)
                
                if l_gt_r == True: # check left-side if delta-left greater than delta-right
                    if auto1_Ml > auto2_Ml:
                        check_all_D[auto1_name] = auto2_name 
                    elif auto2_Ml > auto1_Ml:
                        check_all_D[auto2_name] = auto1_name
                elif r_gt_l == True:
                    if auto1_Mr > auto2_Mr:
                        check_all_D[auto1_name] = auto2_name
                    elif auto2_Mr > auto1_Mr:
                        check_all_D[auto2_name] = auto1_name
                
                break
    
print(check_all_D)

# Write to new xlsx
def write_xlsx(row):
    for key,value in check_all_D.items():
        if row["Times of event (file)"] ==1 and row["Represented read"] == key:
            row["Maybe same event"] = value
    return row

df2 = df.apply(write_xlsx,axis=1)
df2 = df2[["Ref seq name", "Cigar", "Len of read", "Recontructed CIGAR",	
              "Len of matching left", "Len of matching right", 
              "Distance: break to left matching", "Distance: break to right matching",
              "Nu del left side", "Nu del right side", "Total del nu", 
              "Start of ins", "End of ins", "Len of ins", "Seq of ins",	
              "Repair event", "Reconstructed seq", "Times of event (popu)", 
              "Microhomo seq", "Len of microhomo", "Match_mismatch MH",
              "Times of event (file)", "% read contain mismatch","Represented read",	
              "No. Time of event (file)", "Maybe same event"]]
df2.to_excel(f"{sample}_Padded_Final3.xlsx",index=False)
