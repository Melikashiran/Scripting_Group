#! /usr/bin/env python3


"""Combines the count files from the different samples into a single file, 
and performs some basic analysis on the data."""

import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def parse_arguments():
    parser = argparse.ArgumentParser(description="Combine count files and perform analysis.")
    parser.add_argument("-i", "--input_dir", required=True, help="Directory containing count files.")
    parser.add_argument("-o", "--output_file", required=True, help="Output file for combined counts.")
    parser.add_argument("-p", "--plot_file", required=True, help="Output file for the plot.")
    parser.add_argument("-g", "--gene", required=True, default=None, help="Select the gene to analyze.")
    return parser.parse_args()

def combine_counts(input_dir, output_file):
    combined_counts = {}
    
    for filename in os.listdir(input_dir):
        if filename.endswith(".txt"):
            sample_name = filename.split(".")[0]
            with open(os.path.join(input_dir, filename), 'r') as f:
                for line in f:
                    chrom, start, end, count = line.strip().split()
                    key = (chrom, start, end)
                    if key not in combined_counts:
                        combined_counts[key] = {}
                    combined_counts[key][sample_name] = int(count)
    
    with open(output_file, 'w') as f:
        header = "Chrom\tStart\tEnd\t" + "\t".join(sorted(combined_counts[next(iter(combined_counts))].keys())) + "\n"
        f.write(header)
        for key, counts in combined_counts.items():
            line = f"{key[0]}\t{key[1]}\t{key[2]}"
            for sample in sorted(combined_counts[next(iter(combined_counts))].keys()):
                line += f"\t{counts.get(sample, 0)}"
            line += "\n"
            f.write(line)

    
