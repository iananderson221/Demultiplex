# Assignment the First

## Part 1
1. Be sure to upload your Python script. Provide a link to it here:

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | Read 1 | 101 | Phred+33 |
| 1294_S1_L008_R2_001.fastq.gz | Index 1 | 8 | Phred+33 |
| 1294_S1_L008_R3_001.fastq.gz | Index 2 | 8 | Phred+33 |
| 1294_S1_L008_R4_001.fastq.gz | Read 2 | 101 | Phred+33 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
    2. **YOUR ANSWER HERE**
    3. **YOUR ANSWER HERE**
    
## Part 2
1. Define the problem

Create an algorithm to demultiplex paired-end sequencing reads using two index FASTQ files and a list of 24 known index sequences. The algorithm must:

    Sort read pairs into:

        Correctly matched index-pairs (R1 + R2 FASTQ file for each unique index-pair)

        Index-hopped pairs (wrong but known index combinations)

        Reads with unknown indexes (either not in the known list or containing ambiguous bases like N, or below quality threshold)

    Annotate read headers with the index-pair used (e.g. @SEQ_ID AAAAAAAA-TTTTTTTT).

    Report statistics:

        Read-pair count per matched index-pair

        Total count of index-hopped reads

        Total count of unknown index(es)

2. Describe output
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
4. Pseudocode
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
