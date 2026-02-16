#!/usr/bin/env python

import sys

truth_file = sys.argv[1]
classificiation_file = sys.argv[2]


def read_labels(fname):
    groups = list()
    with open(fname) as fin:
        umis = set()
        umi, prev_label = next(fin).strip("\n").split(" ")
        umis.add(umi)
        for line in fin:
            umi, label =  line.strip("\n").split(" ")

            if label == prev_label:
                umis.add(umi)
            else:
                print(f"Adding label {prev_label}: {umis}")
                groups.append(umis)

                umis=set()
                umis.add(umi)
                prev_label = label
        else:
            print(f"Adding label {prev_label}: {umis}")
            groups.append(umis)
    return groups


def print_discordant(truth, classification):
    """Print the discordant cases between truth and classification"""

    # While there are groups in truth, we find the most similar group in
    # classification, and then remove the group from both
    while truth:
        t = truth[0]

        diff = 9999
        best_fit=set()

        for group in classification:
            d = len(t.symmetric_difference(group))
            if d < diff:
                diff = d
                best_fit=group
        if diff:
            uniq_in_truth = t.difference(best_fit)
            uniq_in_classificiation = best_fit.difference(t)

            # Print if there are difference between the sets
            if t.symmetric_difference(best_fit):
                print(f"Unique in truth: {uniq_in_truth}, Unique in classification: {uniq_in_classificiation}")

            classification.remove(best_fit)
            truth.remove(t)
        else:
            # print(f"No difference, removing {truth}")
            truth.remove(t)
            classification.remove(t)
    if classification:
        print("Clusters not assigned to any true cluster:")
        for c in classification:
            print("\t", c)

print("Reading truth")
truth = read_labels(truth_file)

print("Reading classification")
classification = read_labels(classificiation_file)


print("truth:")
for t in truth:
    print("\t", t)

print("classification")
for c in classification:
    print("\t", c)

print("Discordant classifications")
print_discordant(truth, classification)
