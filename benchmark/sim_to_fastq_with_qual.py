#!/usr/bin/env python3

import sys
import random

fname = sys.argv[1]
seed = int(sys.argv[2])

random.seed(seed)

FASTQ_QUAL=r'!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~'

if __name__ == "__main__":
    with open(fname) as fin:
        for i, line in enumerate(fin, 1):
            child, founder, error = line.strip("\n").split("\t")

            # For each base in the child, choose a random quality score (with replacement)
            quality = ''.join(random.choices(FASTQ_QUAL, k=len(child)))

            print(f"@read{i}",sep='')
            print(child)
            print("+")
            print(quality)

