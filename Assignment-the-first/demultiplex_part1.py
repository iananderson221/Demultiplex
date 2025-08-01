#!/usr/bin/env python

import argparse
import gzip as gz
import matplotlib.pyplot as plt

def get_args():
    
    parser = argparse.ArgumentParser(description="A program to initialize input and output file name")
    parser.add_argument("-l", "--read_length", help="specified read length", required=True, type=int)
    parser.add_argument("-f", "--file_name", help="specified input filename", required=True, type=str)
    return parser.parse_args()

# Parse command line arguments
args = get_args()

def init_list(lst: list, value: float=0.0) -> list:
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with 101 values of 0.0.'''
    for num in range(args.read_length):
        lst.append(value)
    return lst
    
my_list: list = []
my_list = init_list(my_list) 

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    return ord(letter) - 33

def populate_list(file: str) -> tuple[list, int]:
    '''Reads a gzipped FASTQ file and accumulates quality scores for each base position.

    Returns:
        - quality_sums: A list of summed quality scores for each base.
        - line_count: Total number of lines processed.'''
    quality_sums = init_list([]) 
    
    line_count = 0

    with gz.open(args.file_name,"rt") as fh:
        for line in fh:
            if line_count % 4 == 3:
                qual = line.strip()
                for i, score in enumerate(qual):
                    quality_sums[i] += convert_phred(score)
                    
            line_count += 1
    return quality_sums, line_count

my_list, num_lines = populate_list(args.file_name)

for i, quality_sum in enumerate(my_list):
    mean = quality_sum / (num_lines / 4)
    my_list[i] = mean
    print(f"{i}\t{my_list[i]}")

x = list(range(len(my_list)))
y = my_list
plt.figure(figsize=(12, 6))
plt.plot(x, y)
plt.grid()
plt.title("Mean Phred Quality Score per Base Position")
plt.xlabel("Base Position")
plt.ylabel("Mean Phred Quality Score")

plot_name = "R4_mean_quality_scores.png"
plt.savefig(plot_name)
print(plot_name)
plt.show()