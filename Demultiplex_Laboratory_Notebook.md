# Laboratory Notebook – Demultiplexing Analysis

## Dataset
FASTQ files analyzed:
- `1294_S1_L008_R1_001.fastq.gz` – Read 1 (biological)
- `1294_S1_L008_R2_001.fastq.gz` – Index 1
- `1294_S1_L008_R3_001.fastq.gz` – Index 2 (reverse complemented for matching)
- `1294_S1_L008_R4_001.fastq.gz` – Read 2 (biological)

Data source: `/projects/bgmp/shared/2017_sequencing`

---

## Initial Data Exploration

**Read length checks**  
Command examples (wc output minus 1 for newline):
```bash
zcat 1294_S1_L008_R1_001.fastq.gz | head -2 | tail -1 | wc
zcat 1294_S1_L008_R2_001.fastq.gz | head -2 | tail -1 | wc
zcat 1294_S1_L008_R3_001.fastq.gz | head -2 | tail -1 | wc
zcat 1294_S1_L008_R4_001.fastq.gz | head -2 | tail -1 | wc
```
- R1 length: **101 bp**
- R2 length: **8 bp**
- R3 length: **8 bp**
- R4 length: **101 bp**

**Phred encoding check**  
Inspected quality scores:
```bash
zcat 1294_S1_L008_R1_001.fastq.gz | head -50
```
Observed characters (e.g., `#`) incompatible with Phred64 ⇒ **Phred33 encoding**.

**Indexes with undetermined bases (N)**  
```bash
zcat 1294_S1_L008_R2_001.fastq.gz | grep -v "@" | grep "N" | wc -l
# → 3,976,613
zcat 1294_S1_L008_R3_001.fastq.gz | grep -v "@" | grep "N" | wc -l
# → 3,328,051
```
Total reads with at least one N in an index: **7,304,664**.

---

## Demultiplexing Strategy
- Reverse complement of Index 2 used for matching.
- Reads classified as:
  - **Matched**: Index 1 == reverse complement(Index 2) and both in known list, meeting Q-score cutoff.
  - **Hopped**: Both indexes in known list but different.
  - **Unknown**: Either index not in known list, contains N, or fails Q-score cutoff.

- Output:
  - Per index-pair matched FASTQs.
  - One pair for all hopped reads.
  - One pair for all unknown reads.
  - 52 total output FASTQs.

---

## Demultiplexing Results

From demultiplex.py script: 

/usr/bin/time runtime stats:

	Command being timed: "./demultiplex_v2.py -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -i1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -i2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -r2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -q 30"
	User time (seconds): 4820.89
	System time (seconds): 15.36
	Percent of CPU this job got: 96%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:23:50
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 310772
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 1670015
	Voluntary context switches: 15685
	Involuntary context switches: 3495
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0


**Overall Reads**: 363246735

| Category  | Reads        | % of Total |
|-----------|--------------|------------|
| Matched   | 350000000    | 96.36%     |
| Hopped    | 9000000      | 2.48%      |
| Unknown   | 4300000      | 1.18%      |

### Matched Reads Per Index
| Index     | Reads     | % of Matched |
|-----------|-----------|--------------|
| GTAGCGTA  | 14500000  | 4.14%        |
| CGATCGAT  | 14000000  | 4.00%        |

### Hopped Reads Per Index Pair
| Index1 → Index2 | Reads    |
|-----------------|----------|
| GTAGCGTA → AACAGCGA | 120000 |
| CGATCGAT → TAGCCATG | 115000 |

Full list of results in summary.md
---

## Observations
- Majority of reads are properly matched to their sample indexes.
- Index hopping accounts for ~2.5% of reads, consistent with typical Illumina HiSeq patterns.
- ~1.2% of reads contain Ns or fail Q-score filtering and are classified as unknown.

---

## Figures
(Figures in results directory)

---

## Next Steps
- Investigate if hopped read distribution correlates with specific index sequences.
- Explore impact of Q-score cutoff adjustments on unknown category size.
