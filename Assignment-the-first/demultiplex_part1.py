#!/usr/bin/env python

# Import required libraries
import argparse             # For parsing command line arguments
import gzip as gz           # For reading gzipped FASTQ files
import matplotlib.pyplot as plt  # For plotting

def get_args():
    """
    Parse command line arguments for input FASTQ file and read length.
    """
    parser = argparse.ArgumentParser(description="Calculate and plot per-base mean quality scores from a gzipped FASTQ file.")
    parser.add_argument("-l", "--read_length", help="Specified read length", required=True, type=int)
    parser.add_argument("-f", "--file_name", help="Input gzipped FASTQ filename", required=True, type=str)
    return parser.parse_args()

# Parse arguments from the command line
args = get_args()

def init_list(lst: list, value: float = 0.0) -> list:
    """
    Initialize a list with a specified number of elements (based on read length),
    each initialized to a given value (default 0.0).
    """
    for _ in range(args.read_length):
        lst.append(value)
    return lst

# Initialize list to store summed quality scores at each base position
my_list: list = []
my_list = init_list(my_list)

def convert_phred(letter: str) -> int:
    """
    Convert a single ASCII character to its corresponding Phred score.
    Phred+33 encoding is assumed.
    """
    return ord(letter) - 33

def populate_list(file: str) -> tuple[list, int]:
    """
    Reads a gzipped FASTQ file and accumulates quality scores per base position.
    
    Args:
        file: Path to the gzipped FASTQ file.

    Returns:
        quality_sums: List of summed quality scores for each base position.
        line_count: Total number of lines read from the file (should be a multiple of 4).
    """
    quality_sums = init_list([])  # Create a fresh list for quality score sums
    line_count = 0

    # Open and read the gzipped FASTQ file
    with gz.open(args.file_name, "rt") as fh:
        for line in fh:
            if line_count % 4 == 3:  # Every 4th line is the quality string
                qual = line.strip()
                for i, score in enumerate(qual):
                    quality_sums[i] += convert_phred(score)
            line_count += 1

    return quality_sums, line_count

# Process the file and get per-base quality sums and total line count
my_list, num_lines = populate_list(args.file_name)

# Calculate and store mean quality score at each base position
for i, quality_sum in enumerate(my_list):
    mean = quality_sum / (num_lines / 4)  # Divide by number of reads
    my_list[i] = mean
    print(f"{i}\t{my_list[i]}")  # Print position and mean quality score

# Prepare data for plotting
x = list(range(len(my_list)))  # Base positions
y = my_list                   # Corresponding mean quality scores

# Create the plot
plt.figure(figsize=(12, 6))
plt.plot(x, y)
plt.grid()
plt.title("Mean Phred Quality Score per Base Position")
plt.xlabel("Base Position")
plt.ylabel("Mean Phred Quality Score")

# Save the plot to a file
plot_name = "R4_mean_quality_scores.png"
plt.savefig(plot_name)

# Output the name of the saved plot file
print(plot_name)

# Display the plot on screen
plt.show()
