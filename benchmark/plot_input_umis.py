#!/usr/bin/env python3

import argparse
import itertools
import editdistance

def main(input, distance_threshold):
    nodes = dict()
    mutations = list()

    with open(input) as fin:
        for line in fin:
            umi, founder, error = line.strip("\n").split("\t")

            if umi not in nodes:
                nodes[umi] = 1
                if error:
                    mutations.append((founder, error, umi))
            else:
                nodes[umi]+=1

            if founder not in nodes:
                nodes[founder] = 1

    if distance_threshold:
        distances = list()
        for umi1, umi2 in itertools.combinations(nodes.keys(), 2):
            dist = editdistance.eval(umi1, umi2)

            if dist and dist <= distance_threshold:
                distances.append((umi1, dist, umi2))
    else:
        distances = list()

    print_graph(nodes, mutations, distances)

def print_graph(nodes, mutations, distances):
    # Header
    print("digraph UMIs {")
    # Print the nodes
    for i, n in enumerate(nodes):
        count = nodes[n]
        name = f"{n} ({count})"
        print(f'{n} [label="{name}"];')

    # Print the mutations
    for founder, error, umi in mutations:
        # print(f'{founder} -> {umi}')# [label="{error}]";')
        print(f'{founder} -> {umi} [label="{error}" color="red"];')

    # Print the distances
    for founder, distance, umi in distances:
        print(f'{founder} -> {umi} [label="{distance}" color="blue"];')
    # Closing bracket
    print("}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="File containing simulated UMI's")
    parser.add_argument("--distance", help="Show the edit distance between UMI's", type=int, default=0)

    args = parser.parse_args()

    main(args.input, args.distance)
