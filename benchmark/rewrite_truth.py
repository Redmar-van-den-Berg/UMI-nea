#!/usr/bin/env python3

from sim_to_fastq import fix_length
import sys

fname=sys.argv[1]
length=int(sys.argv[2])

with open(fname) as fin:
    for line in fin:
        umi, id_ = line.strip("\n").split(" ")
        print(fix_length(umi, length), id_)
