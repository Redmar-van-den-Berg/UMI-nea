#!/usr/bin/env python3

import sys

fname = sys.argv[1]
length = int(sys.argv[2])

def fix_length(umi, lenght):
    """Fix UMI to the specified length

    If the UMI is longer than the specified lenght, cut off the end
    If the UMI is shorter, padd the end with A's
    """
    modified = umi + 'A'*(lenght-len(umi))
    return modified[:lenght]

if __name__ == "__main__":
    with open(fname) as fin:
        for i, line in enumerate(fin, 1):
            child, founder, error = line.strip("\n").split("\t")

            # Fix the length of the UMI
            child = fix_length(child, length)

            print("@read",i, "_", child, sep='')
            print("A")
            print("+")
            print("?")

