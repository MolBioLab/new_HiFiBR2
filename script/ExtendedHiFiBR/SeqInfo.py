# -*- coding: utf-8 -*-
"""
Created on Sat Nov  6 00:49:12 2021

@author: DELL
"""
import csv
import sys
import warnings

from os import listdir, getcwd
from os.path import isfile, join, isdir

#import tkinter as tk
#from tkinter import *

#from tqdm.tk import trange, tqdm
import time



# SYS FUNCTION ----------------------------------------------------------------
def increase_field_size():
    """ 
    Increase data structure (strings, lists,...) size. 
    
    """
    maxInt = sys.maxsize
    decrease = True
    while decrease:
        decrease = False
        try:
            csv.field_size_limit(maxInt)
        except OverflowError:
            maxInt = int(maxInt/10)
            decrease = True

# GET INPUT FILE --------------------------------------------------------------
def get_refSeq(ref_path):
    """ 
    Get reference sequence from .fa file in working directory. 
    
    Args:
        ref (str): path to reference sequence file name with extension ".fa"
        
    Returns: 
        reference sequence (str)
    """
    reader = csv.reader(open(ref_path, "r"), delimiter = '\t')
    next(reader)
    ref_seq = next(reader)[0].strip()
    return ref_seq