# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 14:12:25 2021

@author: mobilab5
"""

#!/usr/bin/env python3

import pysam
import sys

insam= sys.argv[1]
samfile = pysam.AlignmentFile(insam, "rb")
outfile = pysam.AlignmentFile("-", "wb", template=samfile)

for aln in samfile:
    ys= aln.reference_end
    if not ys:
        ys= -1
    aln.setTag('YS', ys)
    outfile.write(aln)

samfile.close()
outfile.close()
sys.exit()
