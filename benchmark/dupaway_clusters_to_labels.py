#!/usr/bin/env python3

import sys

input_file=sys.argv[1]
cluster_file=sys.argv[2]

def read_clusters(fname):
    """
    Generator which yields groups of reads which belong together
    """
    with open(fname) as fin:
        read = next(fin).strip("\n")
        reads = [read]
        read = next(fin).strip("\n")
        while True:
            # If we find the next cluster
            if read.startswith("@"):
                yield reads
                reads = [read]
            else:
                reads.append(read[2:])

            try:
                read = next(fin).strip("\n")
            except StopIteration:
                yield reads
                return

clusters = dict()
for i, reads in enumerate(read_clusters(cluster_file),1):
    clusters.update({read: i for read in reads})

# Go through the FastQ file, printing the read and the cluster assignment
with open(input_file) as fin:
    while True:
        try:
            name = next(fin).strip("\n")
        except StopIteration:
            break
        sequence = next(fin).strip("\n")
        # Skip second header and qual
        next(fin);next(fin)
        print(sequence, clusters[name])
