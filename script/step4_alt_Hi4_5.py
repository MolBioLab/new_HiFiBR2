#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  6 23:44:33 2021

@author: mobilab5

UPDATE:
    - Multi-thread with lock
    - Add parameter "line" to function check_non_temp_ins()
"""
# LIBRARY

import sys
import os
from os import listdir, getcwd
from os.path import isfile, join, isdir

import csv
import pandas as pd
import numpy as np
import math
from sklearn.neighbors import KernelDensity
from sklearn.utils.fixes import parse_version
from scipy.signal import argrelextrema
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks

import re
from itertools import repeat
from diff_match_patch import diff_match_patch

import subprocess
import concurrent.futures as ccf
from tqdm import tqdm
from multiprocessing import Manager

from xlsxwriter.workbook import Workbook

# =============================================================================
# FUNCTIONS
# =============================================================================

# SYS FUNCTION----------------------------------------------------------------
def increase_field_size():
    maxInt = sys.maxsize
    decrease = True
    while decrease:
        decrease = False
        try:
            csv.field_size_limit(maxInt)
        except OverflowError:
            maxInt = int(maxInt/10)
            decrease = True
    
def get_files(directory):   # check if executing file is .sam file 
    files = [directory + '/' + f for f in listdir(directory) 
             if (isfile(join(directory, f)) & (f.find(".sam") != -1))]
    return files

# JUNC ANALYSIS FUNCTION------------------------------------------------------- 
def retrieve_ref_seq(ref):
    reader = csv.reader(open(ref, "r"), delimiter = '\t')
    header = next(reader)
    ref_seq = next(reader)[0].strip()
    return ref_seq

def prepare_line(new_line_list):
    new_line_list = [str(s) for s in new_line_list]
    new_line = '\t'.join(new_line_list)
    return new_line

def write_fasta_line(name_list, out_sequence, out_file):
    fasta_name = '>' + '_'.join(name_list)
    out_file.write(fasta_name + '\n')
    out_file.write(out_sequence + '\n')

def make_subcigar_L(cigar):    
    MDI_list = ["M", "D", "I"]
    subcigar_L = []
    cigar_len = len(cigar)
    i = 0
    while i < cigar_len:
        sub_cigar = ""
        while cigar[i] not in MDI_list:
            sub_cigar = sub_cigar + cigar[i]
            i += 1
        sub_cigar = sub_cigar + cigar[i]
        subcigar_L.append(sub_cigar)
        i += 1
    return subcigar_L

def check_by_nu(num1, num2, seq, ins_seq, side):    # diff
    i = 0 if side == "l" else 1                 
    try:
        if side == "l":                       
            while seq[i] == ins_seq[i]:
                i += 1
        else:
            while seq[-i] == ins_seq[-i]:
                i += 1
    except IndexError:
        pass
       
    if side == "l":
        num1 = num1 + i
        return num1, i
    else:
        num2 = num2 + (i-1)
        return num2
    
def report_indel_info(subcigar_L,read_seq,read_len):
    M_l = int(subcigar_L[0][:-1])
    M_r = int(subcigar_L[-1][:-1])
    dist2brk_l = M_l - frag_l
    dist2brk_r = M_r - frag_r
    del_l = dist2brk_l if dist2brk_l < 0 else 0
    del_r = dist2brk_r if dist2brk_r < 0 else 0

    if del_l < 0:
        ins_start = M_l
    else:
        ins_start = frag_l
    if del_r < 0:
        ins_end = M_r
    else:
        ins_end = frag_r
        
    ins_len = read_len - (ins_start + ins_end)
    ins_len = 0 if ins_len < 0 else ins_len
    ins_seq = read_seq[ins_start:-ins_end]
    del_total = del_l + del_r
    uncut_seq_l = ""
    uncut_seq_r = ""
    
    if ins_len > 0: # diff
        if dist2brk_l > 0:
            ins_start = check_by_nu(frag_l, 0, cut_reg_seq, ins_seq, "l")[0]
            uncut_seq_l = ref_seq[frag_l:ins_start]            
            
        if dist2brk_r > 0:
            ins_end = check_by_nu(0, frag_r, cut_reg_seq, ins_seq, "r")
            uncut_seq_r = ref_seq[-ins_end:-frag_r]
            
        ins_seq = read_seq[ins_start:-ins_end] 
        ins_len = len(ins_seq)
    
    return M_l, M_r, dist2brk_l, dist2brk_r, del_l, del_r, del_total, ins_start, ins_end, ins_len, ins_seq, uncut_seq_l, uncut_seq_r

def check_fake_ins(subcigar_L,read_seq,read_len):
    M_l, M_r, dist2brk_l, dist2brk_r, del_l, del_r, del_total, ins_start, ins_end, ins_len, ins_seq, uncut_seq_l, uncut_seq_r = report_indel_info(subcigar_L,read_seq,read_len)
    sub_l = subcigar_L[0]
    sub_r = subcigar_L[-1]
    
    # Only check cigar which has M at both side & has del + ins    
    if sub_l.find("M") != -1 and sub_r.find("M") != -1 and del_total < 0 and ins_len > 0:
        del_start = frag_l + del_l
        del_end = frag_r + del_r
        del_seq = ref_seq[del_start:-del_end]
    
        if 0 > del_l: 
             temp = check_by_nu(int(sub_l[:-1]), 0, del_seq, ins_seq, "l")                                              
             sub_l = str(temp[0]) + "M"  
             ins_seq = ins_seq[temp[1]:] # Useful in case: fake ins nu belong to both r/l side -> Choose 1 side
                      
        if 0 > del_r:              
             sub_r = str(check_by_nu(0, int(sub_r[:-1]), del_seq, ins_seq ,"r")) + "M"                       
        subcigar_L = [sub_l] + subcigar_L[1:-1] + [sub_r]
        return subcigar_L
             
    else:
        return subcigar_L

def match_mis(read,ref):
    mismatch_L = []
    for i in range(len(read)):
        if read[i] == ref[i]:
            mismatch_L.append('0')
        else:
            mismatch_L.append('1')
    return mismatch_L

def window_info(mismatch_pos,direction,M_l_i):
    X_ori = np.array(mismatch_pos)
    X_data = X_ori[:, np.newaxis]
    X_plot = np.linspace(min(X_ori)-10, max(X_ori)+10, 1000)[:, np.newaxis]
    
    # Gaussian KDE Plot        
    kde = gaussian_kde(X_ori)
    f = kde.covariance_factor()
    bw = f * X_ori.std()
    kde = KernelDensity(kernel='gaussian', bandwidth=bw).fit(X_data)
    log_dens = kde.score_samples(X_plot)
    
    x = X_plot[:, 0]
    y = np.exp(log_dens)
    peaks, _ = find_peaks(y, height=0.05)
   
    # Get cut position (window size):     
    if direction == "l":                       
        if len(peaks) == 0 or len(peaks) > 0:
            cut_pos = mismatch_pos[0]
        elif len(peaks) == 1:
            min_local_i = list(argrelextrema(y, np.less)[0])
            if len(min_local_i) == 0:
                cut_pos = mismatch_pos[0]
            else:
                cut_pos = int(np.floor(x[min(min_local_i)]))          
        window = M_l_i - cut_pos
        
    elif direction == "r":
        if len(peaks) == 0 or len(peaks) > 0:
            cut_pos = mismatch_pos[-1]
        elif len(peaks) == 1:
            min_local_i = list(argrelextrema(y, np.less)[0])
            if len(min_local_i) == 0:
                cut_pos = mismatch_pos[-1]
            else:
                cut_pos = int(np.floor(x[max(min_local_i)]))
        window = cut_pos + 1
    
    return window

def get_mismatch(read,ref,qual,direction,M_l_i):
    mismatch = match_mis(read,ref)
    mismatch_pos_temp = [i for i, elem in enumerate(mismatch) if "1" in elem] 
    
    mismatch_qual = [qual[i] for i in mismatch_pos_temp]
    cond_qual = [ord(i)-33 >= 30 for i in mismatch_qual]
    mismatch_pos = [x for x, y in zip(mismatch_pos_temp, cond_qual) if y == True]
    if len(mismatch_pos) == 0:
        window = 0
    elif len(mismatch_pos) == 1:
        if direction == "l":
            window = M_l_i - mismatch_pos[0]
        elif  direction == "r":   
            window = mismatch_pos[0] + 1
    else:
        window = window_info(mismatch_pos,direction,M_l_i)
    return window 

def check_non_temp_ins(line, subcigar_L, read_seq): #fixed
    sub_l = M_l = subcigar_L[0]
    sub_r = M_r = subcigar_L[-1]
    
    qual = line[10]
    M_l_i = int(M_l[:-1])
    M_r_i = int(M_r[:-1])
    read_l = read_seq[:M_l_i]
    read_r = read_seq[-M_r_i:]  
    ref_l = ref_seq[:M_l_i]
    ref_r = ref_seq[-M_r_i:]
    qual_l = qual[:M_l_i]
    qual_r = qual[-M_r_i:]  
    
    # A. Left side window
    thres_l = get_mismatch(read_l,ref_l,qual_l,"l",M_l_i) + 1
    # B. Right side window
    thres_r = get_mismatch(read_r,ref_r,qual_r,"r",0)
    
    # C. 
    exact_ticker = 0
    if len(subcigar_L) == 1:
        if len(ref_seq) > len(read_seq):
            return subcigar_L
        exact_ticker = 1
        sub_l = M_l = str(frag_l) + "M"
        sub_r = M_r = str(frag_r) + "M"
        
    if sub_l.find("M") != -1 and sub_r.find("M") != -1:
        # 1. Left
        match_index_l = int(sub_l[:-1])
        count_l = 1
        mismatch_L_l = []
        while count_l < thres_l:
            if read_seq[match_index_l - count_l] != ref_seq[match_index_l - count_l]:
                mismatch_index = match_index_l - count_l + 1
                mismatch_L_l.append(mismatch_index)
            count_l += 1
        if len(mismatch_L_l) > 0:
            min_mismatch_index = min(mismatch_L_l)
            M_l = str(min_mismatch_index - 1) + "M"
        
        # 2. Right
        match_index_r = int(sub_r[:-1])
        count_r = 0
        mismatch_L_r = []
        while count_r < thres_r:
            if read_seq[-(match_index_r - count_r)] != ref_seq[-(match_index_r - count_r)]:
                mismatch_index = match_index_r - count_r
                mismatch_L_r.append(mismatch_index)
            count_r += 1
        if len(mismatch_L_r) > 0:   
            min_mismatch_index = min(mismatch_L_r)
            M_r = str(min_mismatch_index - 1) + "M"

        if int(M_l[:-1]) + int(M_r[:-1]) == len(ref_seq) and exact_ticker == 1:
            return subcigar_L
        
        # 3. Results
        new_subcigar_L = [M_l] + subcigar_L[1:-1] + [M_r]
        return new_subcigar_L
    else:
        return subcigar_L
            
def report_class(dels, ins_len, read, ref, check):
    # A. For cigar_len = 1 and read_len = ref_len. # diff
    if check == 1:       
        # 1. match_mismatch list
        mismatch_L = []
        cut_reg_read = read[frag_l:-frag_r]
        for i in range(len(cut_reg_read)):
            if cut_reg_read[i] == cut_reg_seq[i]:
                mismatch_L.append('0')
            else:
                mismatch_L.append('1')
        # 2. Check if there's any mismatches in cut_reg (Label as "abnormal")
        num_mis = mismatch_L.count("1")            
        if num_mis > 0:                               
            return "abnormal" 
        else:
            return "exact"
    # B. For other cases
    else: 
        if dels == 0:
            if ins_len == 0:
                return 'exact'
            else:
                return 'insertion'
        else:
            if ins_len == 0:
                return 'deletion'
            else:
                return 'complex'
                
def get_working_sequence(read_len, seq_class, recon_read, read_seq):
    if read_len != ref_len and seq_class != "exact":
        working_seq = recon_read
    if read_len != ref_len and seq_class == "exact":
        working_seq = read_seq
    if read_len == ref_len and seq_class != "exact":
        working_seq = recon_read
    if read_len == ref_len and seq_class == "exact":
        working_seq = ref_seq
    return working_seq

def report_seq_micro(seq, read_index):
    try:
        read_micro = seq[read_index]
        return read_micro
    except IndexError:
        read_micro = ''
        return read_micro
    
def report_MH(item, ref_seq, frag_r, frag_l):
    side_L = ["l", "r"]
    for side in side_L:
        
        # 1. Get 1st pos that may have MH event
        if side == "l":
            read_index = read_index_2 = int(item[4]) - 1   #left item[4] = M_l
            ref_index = ref_index_2 = -(int(item[5]) + 1)  #right item[5] = M_r
            if_statement = ref_index_2 < -(frag_r)
            i = -1  # used as pointer, move 1bp each time to the left (both read & ref)
            
        if side == "r":
            read_index = read_index_2 = -int(item[5]) #right 
            ref_index = ref_index_2 = int(item[4]) #left
            if_statement = ref_index_2 > frag_l -1
            i = 1 # used as pointer, move 1bp each time to the right (both read & ref)

        read_micro = report_seq_micro(item[-1], read_index) # item[-1] = working_seq
        ref_micro = report_seq_micro(ref_seq, ref_index)
        
        count_all = 1   # already counted 1st read = ref MH
        # 2. Allow 20% mismatch per MH len        
        if read_micro == ref_micro:            
            count_mis = 0
            mismatch_L = [] # "M" match, "u" unmatch
            while True:  
                if read_micro != ref_micro:
                    count_mis += 1
                    mismatch_L.append("u")
                else:
                    mismatch_L.append("M") 
                count_all += 1
                read_index_2 += i
                ref_index_2 += i
                read_micro = report_seq_micro(item[-1], read_index_2)
                ref_micro = report_seq_micro(ref_seq, ref_index_2)
                percent = math.ceil((count_mis / count_all)*100)
                if if_statement or percent > 20:
                    break          
            mismatch_S = "".join(mismatch_L)
            
            # Remove unmatch nu in MH at the end of MH string
            for w in mismatch_S[::-1]:  # Execute from the end to the beginning of string
                if w == "u":
                    mismatch_S = mismatch_S[:-1]
                    read_index_2 -= i
                    ref_index_2 -= i
                else:
                    break
        else:
            mismatch_S = ""    
        # 3. Get MH seq & mismatch list    
        if side == "l":
            read_MH_l = item[-1][read_index_2+1 : read_index+1]
            # ref_MH_r = ref_seq[ref_index_2+1 : ref_index+1] # no need
            mismatch_S_l = mismatch_S[::-1]    
        if side == "r":
            read_MH_r = item[-1][read_index : read_index_2] # no need
            # ref_MH_l = ref_seq[ref_index : ref_index_2]
            mismatch_S_r = mismatch_S
    read_MH_total = read_MH_l + read_MH_r
    mismatch_S_total = mismatch_S_l + mismatch_S_r

    return read_MH_total, mismatch_S_total

def execute_each_line(line, recon_L, no_dup_L, read_err_D, read_err_distri_D, lock):
    with lock: # fixed
        cigar = line[5]              
            
        # 1.1.1. If cigar includes "S", "*", "H" tag -> Eliminate read 
        blackL = ["S", "*", "H"] #S appears in CLC, * in Bowtie2, H in BWA
        blackL_ticker = 0
        for item in blackL:
            if cigar.find(item) != -1:
                blackL_ticker = 1
                break
            
        if blackL_ticker != 0:
            with open(f[:-4] + "_Drop.csv", "a+") as drop:
                drop.write("black_list \t" + line[0] + "\n")
            return None
        
        # 1.1.2. Check & Alter cigar string & Fill info (del_total, ins_len, seq_class)
        read_seq = line[9]
        read_len = len(read_seq)
        subcigar_L = make_subcigar_L(cigar) 
        
        char_L = [i[-1] for i in subcigar_L]    # For insertion only, 'cause cant analyse with common pipeline
        if char_L == ["M","I","M"]:
            M_l = int(subcigar_L[0][:-1])
            M_r = int(subcigar_L[-1][:-1])
            dist2brk_l = M_l - frag_l
            dist2brk_r = M_r - frag_r
            
            ins_seq = read_seq[M_l:-M_r]
            ins_len = len(ins_seq)
            seq_class = "insertion"
            del_l, del_r, del_total, ins_start, ins_end, uncut_seq_l, uncut_seq_r, recon_cigar = 0, 0, 0, 0, 0, "", "", ""
            
        elif (len(subcigar_L) == 1) and (int(subcigar_L[0][:-1]) == ref_len):    # diff
            test_L = M_l, M_r, dist2brk_l, dist2brk_r, del_l, del_r, del_total, ins_start, ins_end, ins_len, ins_seq, uncut_seq_l, uncut_seq_r, recon_cigar = read_len, read_len, 0, 0, 0, 0, 0, 0, 0, 0, "", "", "", ""             
            seq_class = report_class(del_total, ins_len, read_seq, ref_seq, 1)
            check_uncut = "Maybe" 
            
        else:               
            subcigar_L = check_non_temp_ins(line, subcigar_L, read_seq)     # diff: non_temp 1st, fake_ins 2nd                       
            info = report_indel_info(subcigar_L, read_seq, read_len)
            # print("info before", info)
            
            subcigar_L = check_fake_ins(subcigar_L,read_seq,read_len)                                    
            info = report_indel_info(subcigar_L, read_seq, read_len)
            # print("info after", info)
            test_L = M_l, M_r, dist2brk_l, dist2brk_r, del_l, del_r, del_total, ins_start, ins_end, ins_len, ins_seq, uncut_seq_l, uncut_seq_r = report_indel_info(subcigar_L,read_seq,read_len)                
            seq_class = report_class(del_total, ins_len, "", "", 0)            
            recon_cigar = "".join(subcigar_L)
            check_uncut = "True" # This mode report True "exact" class (just cut out cut_reg, not already repair junction)
            
        # 1.1.3. Drop abnormal
        if seq_class == "abnormal": # diff
            with open(f[:-4] + "_Drop.csv", "a+") as drop:
                drop.write("abnormal \t" + line[0] + "\n")
            return None
        
            # Note: del_total not be used for further analyse -> Fixed del_total
        if len(uncut_seq_l) + len(uncut_seq_r) == 34:
            del_total = -del_total              
        elif dist2brk_l > 0 or dist2brk_r > 0:  # uncut 1 side/ 2 sides
            del_total = ref_len - M_l - M_r               
        else:
            del_total = -del_total
        
        # 1.1.3. Reconstructed reads: NOT for insertion only
        l_end_index = frag_l + del_l
        r_start_index = frag_r + del_r
    
        half_ref_l = ref_seq[:l_end_index]
        half_ref_r = ref_seq[-r_start_index:]  
        
        recon_read = half_ref_l + uncut_seq_l + ins_seq + uncut_seq_r + half_ref_r # diff
            #  Note: Fix seq_class later, 'cause it is used as "working seq" parameter recently
        working_seq = get_working_sequence(read_len, seq_class, recon_read, read_seq) 
            
        if len(read_seq) != len(working_seq) and (seq_class != "insertion"):    # same as old Hi-FiBR        
            info = [M_l, M_r, dist2brk_l, dist2brk_r, del_l, del_r, del_total, ins_start, ins_end, ins_len, ins_seq, uncut_seq_l, uncut_seq_r]
            info2 = [str(i) for i in info]
            
            with open(f[:-4] + "_Drop.csv", "a+") as drop:
                drop.write("diff len \t" + line[0] +"\t" + str(len(read_seq)) + "\t" + str(len(working_seq)) + "\t" + "\t".join(info2) + "\n")
            return None #This should not occur often; it is an additional precaution for junk reads
            
        recon_L.append(working_seq)
        
        # 1.1.4. Error distribution           
        if working_seq not in read_err_D:
            read_err_D[working_seq] = 0
            read_err_distri_D[working_seq] = {}
            bp = 0
            while bp < len(working_seq):
                read_err_distri_D[working_seq][bp + 1] = 0
                bp += 1
        if read_seq != working_seq:
            read_err_D[working_seq] += 1
            bp_pos = 0
            while bp_pos < len(working_seq):
                read_bp = read_seq[bp_pos]
                recon_bp = working_seq[bp_pos]
                if read_bp != recon_bp:
                    read_err_distri_D[working_seq][bp_pos + 1] += 1
                bp_pos += 1
                
        # 1.1.5. Append all info (in a line)
            # Note: Fix seq_class info 
        if seq_class == "insertion":
            if M_l + M_r == ref_len - len(cut_reg_seq):
                seq_class = "insertion"
            else:
                seq_class = "insertion uncut 1 side"
                 
        else:            
            if (seq_class == "exact") and (check_uncut == "Maybe"):
                seq_class = "uncut 2 sides"
                recon_cigar = str(read_len) + "M"
                
            elif (dist2brk_l > 0) or (dist2brk_r > 0):
                if ins_len >0:
                    seq_class = "complex uncut 1 side"
                else:
                    seq_class = "deletion uncut 1 side"
                    
            elif (dist2brk_l > 0) and (dist2brk_r > 0):
                if ins_len >0:
                    seq_class = "complex uncut 2 sides"
                else:
                    seq_class = "deletion uncut 2 sides" 
            
        additional_line_L = [read_len, recon_cigar, M_l, M_r, dist2brk_l, dist2brk_r, 
                             del_l, del_r, del_total, ins_start, ins_end, ins_len, ins_seq, 
                             seq_class, working_seq]
        new_line_L = line.copy()
        new_line_L.extend(additional_line_L)
        new_line = prepare_line(new_line_L)
        with open(f[:-4] + "_Extended.tsv", "a+") as out:
            out.write(new_line + '\n')
    
        no_dup_line = [line[2], line[5]]        # [ref name, cigar]
        no_dup_line.extend(additional_line_L)   # Add features
        no_dup_line = prepare_line(no_dup_line)
        no_dup_L.append(no_dup_line)            # List of lines (1 line - all info)
        return

def compare2seq(seq1, seq2):
    dmp = diff_match_patch()
    diff_L = dmp.diff_main(seq1, seq2)
    stat_L = [k[0] for k in diff_L] # del = -1, ins = 1, equal = 0 of seq(considering) compare to seq(key)
        
    # 1. seq(considering) totally inside/exact the same seq(key)
    if stat_L == [-1, 0, -1] or stat_L == [-1, 0] or stat_L == [0, -1] or stat_L == [0]:
        return "done1"
    # 2. seq(considering) with longer left-side
    elif stat_L == [1, 0, -1] or stat_L == [1, 0]:
        return "done2"                             
    # 3. seq(considering) with longer right-side
    elif stat_L == [-1, 0, 1] or stat_L == [0, 1]:
        return "done3"
    # 4. seq(considering) with longer both-side
    elif stat_L == [1, 0, 1]:
        return "done4"
    else:
        return "no"
        
def compareToEli(seq, eli_info):
    [eli_group, eli_side, eli_seq] = eli_info
    res = compare2seq(eli_seq, seq) # eli_seq & seq both are old seq (unpadded)
    if eli_side == "left": 
        if res == "done1" or res == "done3": # seq with shorter left-side than eli_seq -> be dropped
            return eli_info
        else:
            return "correct"
    elif eli_side == "right":   
        if res == "done1" or res == "done2": # seq with shorter right-side than eli_seq -> be dropped
            return eli_info
        else:
            return "correct"

def compareToOther(seq, other_info):
    [other_group, other_seq] = other_info # other_group = other_reconSeq. other_seq is old seq (unpadded)
    res = compare2seq(other_seq, seq)
    if res == "done4":
        return other_info   # row is same as other_group
    else:
        return "correct"    # row is uniq

def add_df_value(df, cond, add):
    i = df[df["Reconstructed seq"] == cond].index.values[0]
    v1 = list(df.loc[df.index[i], "Maybe same event"])
    v1.append(add)
    df.loc[df["Reconstructed seq"] == cond, "Maybe same event"] = pd.Series([v1]*df.shape[0])
    return df    

def maybe_same_event(sub_seq, compare_info, func, text, recon, df_info, count_pair = 0):    
    with ccf.ProcessPoolExecutor() as executor:
       results = list(executor.map(func, repeat(sub_seq[0]), compare_info))
    check = [k == "correct" for k in results] # means: smallest seq doesnt belong to any eli groups = True  
    
    if all(check):
        if text == "eli":
            return df_info  
        elif text == "other":
            return df_info, count_pair
        
    else:
        res_short = [k for k in results if k != "correct"][0] # if smallest seq belongs to many compare groups (maybe impossible) -> just take 1 of them
        with ccf.ProcessPoolExecutor() as executor:
            results2 = list(executor.map(func, sub_seq, repeat(res_short)))
        count_drop = len([k for k in results2 if k != "correct"])
        
            # fill info "Maybe same event" to df_info if there's any drop seq in row
        if text == "eli":
            print("drop lines ELI", count_drop)
            df_info = add_df_value(df_info, recon, res_short[0] + " " + str(count_drop))
            return df_info
        
        elif text == "other":
            count_pair += 1
            print("drop lines OTHER", count_drop)
            df_info = add_df_value(df_info, recon, "same" + str(count_pair) + " " + str(count_drop))
            df_info = add_df_value(df_info, res_short[0], "ref" + str(count_pair))
            return df_info, count_pair

# FUNC for extra insertion -----------------------------------------------------------
def check_left(seed_pos, read_seq, read_qual):
    ref_l_side = ref_seq[seed_pos : cutPos1-1]
    read_l_side = read_seq[:len(ref_l_side)]
    read_qual_l_side = read_qual[:len(ref_l_side)]
    
    if len(ref_l_side) != len(read_l_side):
        return "no"
    else:
        mismatch = 0
        # Check if any mismatch (with Q >=30) between left side of read and ref
        # Expect mismatch = 0
        for i in range(len(ref_l_side)):
            if read_l_side[i] != ref_l_side[i]:
                
                # Check if mismtach is actually mismatch (Q>30)
                qual_num = ord(read_qual_l_side[i]) - 33
                if qual_num >= 30:
                    mismatch += 1
                    
        if mismatch == 0:        
            return ref_l_side
        else:
            return "no"

def check_right(seed_pos, read_seq, read_qual):
    ref_r_side = ref_seq[cutPos2 : seed_pos+11]
    read_r_side = read_seq[len(read_seq) - len(ref_r_side):]
    read_qual_r_side = read_qual[len(read_qual) - len(ref_r_side):]
    
    if len(ref_r_side) != len(read_r_side):
        return "no"
    else:
        mismatch = 0
        # Check if any mismatch (with Q >=30) between left side of read and ref
        # Expect mismatch = 0
        for i in range(len(ref_r_side)):
            if read_r_side[i] != ref_r_side[i]:
                
                # Check if mismtach is actually mismatch (Q>30)
                qual_num = ord(read_qual_r_side[i]) - 33
                if qual_num >= 30:
                    mismatch += 1
        
        if mismatch == 0:        
            return ref_r_side
        else:
            return "no"
def extra_insertion(full_SeqInfo, path):
    ins_info = {}   
    for SeqInfo in full_SeqInfo:
        header, seq, qual = SeqInfo
        
        # 1. Find 11 nu for seeding
        seed_l = seq[:11]
        seed_pos_l_L = [m.start() for m in re.finditer(seed_l, ref_seq)]
        
        seed_r = seq[len(seq)-11:]
        seed_pos_r_L = [m.start() for m in re.finditer(seed_r, ref_seq)]
                
        # 2. Execute if both left-side and right-side have seeds
        if len(seed_pos_l_L) > 0 and len(seed_pos_r_L) > 0:
            check_l_L = [ check_left(seed_pos_l, seq, qual) for seed_pos_l in seed_pos_l_L ] # Extending seed left side
            check_l_L2 = []
            for check_l in check_l_L:
                if len(seq) - len(check_l) <= 11: # seq = left-side + insSeq + right-side
                    continue
                if check_l in ["","no"]:
                    continue
                check_l_L2.append(check_l)                
        else:
            continue
    
        # 3. Extending seed right side, if left side of read match ref
        if len(check_l_L2) != 0:
            check_r_L = [ check_right(seed_pos_r, seq, qual) for seed_pos_r in seed_pos_r_L ]
            check_r_L2 = [i for i in check_r_L if i not in ["","no"] ]
        else:
            continue
        
        # 4. Write sam file if both left side & right side of read match ref
        if len(check_r_L2) != 0:
            for check_l in check_l_L2:
                for check_r in check_r_L2:                    
                    Ml = len(check_l)
                    Mr = len(check_r)
                    insSeq = seq[Ml : -Mr]
                    insLen = len(insSeq)
                    
                    if insLen != 0:
                        if insSeq not in ins_info.keys():
                            ins_info[insSeq] = [header]
                        else:
                            ins_info[insSeq].append(header)
                        
                        insSeq_L = list(ins_info.keys())
                        write_lines = [header, check_l, insSeq, check_r, f"{Ml}M{insLen}I{Mr}M"]
                        with open(f"{path}/2_mapping/getInsertion_log.txt", "+a") as f_out:
                            f_out.write("\t".join(write_lines) + "\n")
            
    ins_info2 = {}
    for k, v in ins_info.items():
        times_of_event = len(v)
        represented_read = v[0]
        ins_info2[k] = [times_of_event, represented_read]
    
    print("ins info:", ins_info2) 
    return ins_info2          

# =============================================================================
# RUN 
# =============================================================================

sample = sys.argv[1]
path = f"{sys.argv[2]}/{sample}"
os.chdir(path)

# A. INPUT for analyzing     
ref = sys.argv[3]
ref_seq = retrieve_ref_seq(ref)
ref_len= len(ref_seq)

cutPos1 = int(sys.argv[5])
cutPos2 = int(sys.argv[6])
frag_l = cutPos1 - 1          # length of left fragment
frag_r = ref_len - cutPos2

cut_reg = frag_r - frag_l + 1
cut_reg_seq = ref_seq[frag_l:-frag_r]

directory = "3_HiFiBR"
input_files = get_files(directory)

thread = int(sys.argv[7])

# B. READ each sam file
for f in input_files:
    print(f"alt-HiFiBR AUTO ANALYZE {f}...")
    manager = Manager()
    reader = csv.reader(open(f,"r"), dialect="excel-tab")
    recon_L, no_dup_L, read_err_D, read_err_distri_D = manager.list(), manager.list(), manager.dict(), manager.dict()
    
    
    # 1. READ each line in sam file:
    line_L = []
    for line in reader:
        if line[0][0] == '@':
            continue 
        line_L.append(line)
        
    # Execute each 10 lines (~ 10 cores, could changed by argument)    
    for i in tqdm(range(0,len(line_L), thread)):
        sub_lines = line_L[i:i+thread]
        with ccf.ProcessPoolExecutor() as executor:
            lock = manager.Lock() # fixed
            res = [executor.submit(execute_each_line, line, recon_L, no_dup_L, 
                            read_err_D, read_err_distri_D, lock) for line in sub_lines] # fixed       
        
        
    # 2. EXECUTE ALL READ: Count time of events 
        # 2.1. Remove dup lines (features spilted by "tab")
    no_dup_L = list(set(no_dup_L))
    no_dup_L = list(filter(None, no_dup_L))
    print("no dup ", len(no_dup_L))
    
    split_L = []
    for item in no_dup_L:
        split_L.append(item.split("\t"))
    no_dup_L = split_L     
    
        # 2.2. Fill Time of events & MH features info to EACH LINE
    control_L = []
    sorting_L = []
    for item in no_dup_L:
        w = 1
        line_count = recon_L.count(item[-1])    # item[-1] = working_seq
        control_L.append(item[-1])

        # a. Only del type has MH features
        if "deletion" in item[-2]:              # item[-2] = seq_class
            read_MH_total, mismatch_S_total = report_MH(item, ref_seq, frag_r, frag_l)
            len_micro_total = len(read_MH_total)
            
            end_line_L = [line_count, read_MH_total, len_micro_total, mismatch_S_total]
            full_line_L = item.copy()
            full_line_L.extend(end_line_L)  # Add line_count + 3 more features to item 
                                            # full_line_L: ref_name, cigar, M_l,..., working_seq, time of events,3 MH features 

        # b. Other types has no MH features
        if "deletion" not in item[-2]:
            end_line_L = [line_count, "", "", ""]
            full_line_L = item.copy()
            full_line_L.extend(end_line_L)
        sorting_L.append(full_line_L)
    print("control ", len(control_L))
    print("sort ", len(sorting_L))
        
    # 3. Get all info and write to files
        # 3.1. Writing MH info to "Count_Columns.sam" file
    micro_out = open(f[:-4] + "_Count_Columns.sam", "w")
    sorting_L.sort(key=lambda x:x[-5], reverse=True)    # sorted by working_seq
    
    for sorted_item in sorting_L:
        full_line = prepare_line(sorted_item)
        micro_out.write(full_line + '\n')
    micro_out.close()
    
        # 3.2. Writing .fasta files
    micro_in_reader = csv.reader(open(f[:-4] + "_Count_Columns.sam", "r"), delimiter = '\t')
    complex_out = open(f[:-4] + "_Complex_Seqs.fasta", "w")
    ins_out = open(f[:-4] + "_Insert_Seqs.fasta", "w")
    fasta_out = open(f[:-4] + "_Fasta_Reads.fasta", "w")
    read_err_dist_out = open(f[:-4] + "_read_error_distribution.txt", "w")

    for line in micro_in_reader:
        
        ins_seq, recon_read, read_count, cigar_str, read_class = line[-7], line[-5], str(line[-4]), line[1], line[-6]   # diff
        name_list = [read_count, cigar_str]
        
        if read_class == "complex":
            write_fasta_line(name_list, ins_seq, complex_out)
        if read_class == "insertion":
            write_fasta_line(name_list, ins_seq, ins_out)
        write_fasta_line(name_list, recon_read, fasta_out)
        
        control_count = str(control_L.count(line[-5]))
        read_error = str((read_err_D[recon_read]/int(read_count)) * 100)
        final_line = '\t'.join(line) + '\t' + control_count + '\t' + read_error

        with open(f[:-4] + "_Final.tsv", "a+") as final_out:
            final_out.write(final_line + "\n")

        read_err_D_line = [recon_read, read_count, read_error]
        err_D_key = 1
        for key in read_err_distri_D[recon_read]:
            bp_count = read_err_distri_D[recon_read][err_D_key]
            read_err_D_line.append(bp_count)
            err_D_key += 1
        read_err_out_line = [str(s) for s in read_err_D_line]
        read_err_out_line = '\t'.join(read_err_out_line)
        read_err_dist_out.write(read_err_out_line + '\n')
        
    complex_out.close()
    ins_out.close()
    fasta_out.close()
    read_err_dist_out.close()

    # 4. Change tsv to xlsx
    tsv_reader = csv.reader(open(f"{f[:-4]}_Final.tsv", "rt", encoding= "utf-8"), delimiter='\t')  
    workbook = Workbook(f"{f[:-4]}_Final_temp.xlsx")
    worksheet = workbook.add_worksheet()
    
    label = ["Ref seq name", "Original CIGAR", "Len of read", "Recontructed CIGAR",	
              "Len of matching left", "Len of matching right", 
              "Distance: break to left matching", "Distance: break to right matching",
              "Nu del left side", "Nu del right side", "Total del nu", 
              "Start of ins", "End of ins", "Len of ins", "Seq of ins",	
              "Repair event", "Reconstructed seq", "Times of event (popu)", 
              "Microhomo seq", "Len of microhomo", "Match_mismatch MH",
              "Times of event (file)", "% read contain mismatch"]
    for i in range(len(label)):
        worksheet.write(0,i,label[i])

    for row, data in enumerate(tsv_reader):
        worksheet.write_row(row + 1, 0, data)
    
    # Close the XLSX file.
    workbook.close()
    os.remove(f"{f[:-4]}_Final.tsv")

# D. ADD INFO TO FINAL.XLSX
    print("alt-HiFiBR ADD INFO...")
    # 1. Add "represented read"
    os.chdir(path)
    df_final = pd.read_excel(f"{f[:-4]}_Final_temp.xlsx", engine='openpyxl')
    df_final.rename(columns = {"Original CIGAR":"Cigar"},inplace = True)
    
    col = ["Full Name", "Flag",	"Ref", "Start Pos", "Mapq",	"Cigar", "Pnext", "Rnext",	
           "Tlen", "Seq", "Qual", "End Pos", "Old Cigar", "Old Seq", "Len of read", "Recontructed CIGAR",	
           "Len of matching left", "Len of matching right", "Distance: break to left matching",	
           "Distance: break to right matching",	"Nu del left side", "Nu del right side",	
           "Total del nu", "Start of ins", "End of ins", "Len of ins", "Seq of ins",	
           "Repair event", "Reconstructed seq"]
    
    df_full = pd.read_csv(f"{f[:-4]}_Extended.tsv",sep = "\t",header = None)
    df_full.columns = col
    
    col2 = ["Ref","Cigar"] + col[14:]
    df_get_repre = df_full.drop_duplicates(subset = col2, keep='first')
    
    col3 = col[1:14]
    col3.remove("Cigar")                             
    df_info = df_final.merge(df_get_repre, how="left", on = col2[1:])
    df_info.drop(col3, axis = "columns",inplace=True)
    df_info.rename(columns={"Full Name": "Represented read"},inplace = True)
    
    # 2. Add "Time of events (file) id" & Get "Time of events (file)" sam file
    df_filt = df_info[df_info["Times of event (file)"] >= 2]
    recon_L = list(set(list(df_filt["Reconstructed seq"])))
    df_info["No. Time of event (file)"] = ""

    count = 0
    for recon in recon_L:
        count += 1
        index_L = df_info.index[df_info["Reconstructed seq"]==recon].tolist()
        
        for i in range(len(index_L)):
            df_info["No. Time of event (file)"].iat[index_L[i]] = count     
    
    # 3. Find extra insertion junction from Unmapped
    print("alt-HiFiBR EXTRA INSERTION...")
        # Find insertion
    unaligned2_sam = f"{path}/2_mapping/{sample}_Unaligned2.sam"
    if os.path.exists(unaligned2_sam):
        cmd = f"samtools fastq {unaligned2_sam} > {path}/2_mapping/{sample}_Unaligned2.fq"
        os.system(cmd)
        
        f_input = f"{path}/2_mapping/{sample}_Unaligned2.fq" # "Seq3_insertion_0D-0D-147_droppedByFilter.fq"
        with open(f_input, "r") as f_in:
            lines = f_in.readlines()
            lines = list(map(str.strip, lines))
            
            startLine_L = list(range(0, len(lines), 4))
            full_SeqInfo = [[ lines[startLine], lines[startLine+1], lines[startLine+3] ]
                            for startLine in startLine_L] # [ [Header1, Seq1, Qual1] , [...], [...]] 
        ins_info2 = extra_insertion(full_SeqInfo, path) 
        
            # Add rows to df_info with colName correspond to label variable line 835
        try:
            for k, v in ins_info2.items():
                insSeq = k
                timeofevent, repre = v
                
                # if extra insertion junc already analyzed above
                if insSeq in df_info["Seq of ins"].to_list():
                    df_boolean = ( df_info['Seq of ins'] == insSeq ) & ( df_info['Repair event'] == "insertion" )
                    need_modify_index = df_info[df_boolean].index.tolist()
                    
                    for index_i in need_modify_index:
                        old_timeofevent = int(df_info.loc[df_info.index[index_i], 'Times of event (popu)'])
                        df_info.loc[df_info.index[index_i], 'Times of event (popu)'] = old_timeofevent + timeofevent
                        
                else:
                    Ml = cutPos1 - 1
                    Mr = len(ref_seq) - cutPos2 + 1
                    insLen = len(insSeq)
                    insCigar = f"{Ml}M{insLen}I{Mr}M"
                    readLen = len(ref_seq) - (cutPos2 - cutPos1 + 1) + insLen
                    reconSeq = ref_seq[:cutPos1] + insSeq + ref_seq[cutPos2:]
                    
                    add_ins_info = ["pCOH_CD4_alt_PCR", insCigar, readLen, insCigar,
                                    Ml, Mr, 0, 0, 0, 0, 0,
                                    cutPos1, cutPos1 + insLen, insLen, insSeq,
                                    "insertion", reconSeq, timeofevent,
                                    "", 0, "", 1, "", repre, ""]
                    df_info.loc[len(df_info.index)] = add_ins_info
        except:
            print("no insertion")
    
        # Write fq file of insertion junction (in order to check ins juction by multiple alignment)
    ins_name_L = df_info[df_info["Repair event"] == "insertion"]["Represented read"].to_list()
    if len(ins_name_L) > 0:
        for ins_name in ins_name_L:
            cmd = f"cat {path}/2_mapping/{sample}_Unaligned2.fq | rg -A 3 '{ins_name}' >> {path}/3_HiFiBR/{sample}_insertion.fq"
            os.system(cmd)


    # 4. Filter junctions with "minTimeofEvent" 
    print("alt-HiFiBR FILTER minTimeofEvent...")
    minTimeofEvent = int(sys.argv[8])
    df_info["Times of event (popu)"] = df_info["Times of event (popu)"].astype("int")
    df_info = df_info[df_info["Times of event (popu)"] >= minTimeofEvent]
    df_info.to_excel(f"{f[:-4]}_Final_temp2.xlsx",index=False)

    # 5. Predict tooShort Junc
    print("alt-HiFiBR PREDICT TOOSHORT JUNC...")
    os.chdir(f"{path}/3_HiFiBR")
    cmd = f"python step4_Detect_tooShort.py {sample} {sys.argv[2]} {ref} {sys.argv[4]}"
    os.system(cmd)
    
    os.remove(f"{sample}_Final_temp.xlsx", f"{sample}_Final_temp2.xlsx")