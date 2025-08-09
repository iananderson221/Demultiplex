#!/usr/bin/env python

import argparse
import gzip as gz
import matplotlib.pyplot as plt

# FUNCTIONS
def reverse_complement(seq: str) -> str:
    comp = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
    # avoid building intermediate list; join directly
    return ''.join(comp.get(b, 'N') for b in reversed(seq))

def check_quality(qual: str, cutoff: int) -> bool:
    # Faster: avoid division; compare sum to cutoff*len
    total = 0
    for ch in qual:
        total += (ord(ch) - 33)
    return total >= cutoff * len(qual)

# ARGPARSE
parser = argparse.ArgumentParser(description="Demultiplex paired-end FASTQ reads by dual indexes (exact match).")
parser.add_argument("-r1", required=True, help="Read 1 FASTQ.gz (biological)")
parser.add_argument("-r2", required=True, help="Read 2 FASTQ.gz (biological)")
parser.add_argument("-i1", required=True, help="Index 1 FASTQ.gz")
parser.add_argument("-i2", required=True, help="Index 2 FASTQ.gz (will be reverse-complemented)")
parser.add_argument("-q", required=True, type=int, help="Minimum average Q-score cutoff for index reads")
args = parser.parse_args()

# CONSTANTS
RESULTS_DIR = "/projects/bgmp/ica/bioinfo/Bi622/Demultiplex/results"
known_indexes = [
    "GTAGCGTA", "CGATCGAT", "GATCAAGG", "AACAGCGA",
    "TAGCCATG", "CGGTAATC", "CTCTGGAT", "TACCGGAT",
    "CTAGCTCA", "CACTTCAC", "GCTACTCT", "ACGATCAG",
    "TATGGCAC", "TGTTCCGT", "GTCCTAAG", "TCGACAAG",
    "TCTTCGAC", "ATCATGCG", "ATCGTGGT", "TCGAGAGT",
    "TCGGATTC", "GATCTTGC", "AGAGTCCA", "AGGATAGC"
]
known_set = set(known_indexes)

# PREOPEN EXACTLY 52 OUTPUT FILES (use faster gzip compression)
matched_files = {}
for idx in known_indexes:
    matched_files[idx] = [
        gz.open(f"{RESULTS_DIR}/R1_{idx}-{idx}.fastq.gz", "wt", compresslevel=1),
        gz.open(f"{RESULTS_DIR}/R2_{idx}-{idx}.fastq.gz", "wt", compresslevel=1)
    ]

hopped_R1 = gz.open(f"{RESULTS_DIR}/R1_hopped.fastq.gz", "wt", compresslevel=1)
hopped_R2 = gz.open(f"{RESULTS_DIR}/R2_hopped.fastq.gz", "wt", compresslevel=1)
unknown_R1 = gz.open(f"{RESULTS_DIR}/R1_unknown.fastq.gz", "wt", compresslevel=1)
unknown_R2 = gz.open(f"{RESULTS_DIR}/R2_unknown.fastq.gz", "wt", compresslevel=1)

# COUNTERS
matched_counts = {idx: 0 for idx in known_indexes}
hopped_counts = {}  # dict[idx1][idx2] = count
unknown_count = 0
total_reads = 0

# PROCESS
with gz.open(args.r1, "rt") as r1_f, \
     gz.open(args.r2, "rt") as r2_f, \
     gz.open(args.i1, "rt") as i1_f, \
     gz.open(args.i2, "rt") as i2_f:

    # localize references for speed
    kset = known_set
    mfiles = matched_files
    hR1, hR2 = hopped_R1, hopped_R2
    uR1, uR2 = unknown_R1, unknown_R2
    qcut = args.q
    hop_counts = hopped_counts
    mcounts = matched_counts

    while True:
        # Read one FASTQ record from each stream
        h1 = r1_f.readline()
        if not h1:
            break  # EOF
        s1 = r1_f.readline(); p1 = r1_f.readline(); q1 = r1_f.readline()

        h2 = r2_f.readline(); s2 = r2_f.readline(); p2 = r2_f.readline(); q2 = r2_f.readline()
        h3 = i1_f.readline(); s3 = i1_f.readline(); p3 = i1_f.readline(); q3 = i1_f.readline()
        h4 = i2_f.readline(); s4 = i2_f.readline(); p4 = i2_f.readline(); q4 = i2_f.readline()

        total_reads += 1

        # Extract sequences (strip newlines once)
        idx1_seq = s3.rstrip('\n')
        idx2_seq = reverse_complement(s4.rstrip('\n'))

        # Build tag & annotate headers (preserve trailing newline)
        tag = idx1_seq + "-" + idx2_seq
        h1 = h1.rstrip('\n') + " " + tag + "\n"
        h2 = h2.rstrip('\n') + " " + tag + "\n"

        # Validate & quality check
        valid = (idx1_seq in kset) and (idx2_seq in kset)
        goodq = check_quality(q3.rstrip('\n'), qcut) and check_quality(q4.rstrip('\n'), qcut)

        if valid and goodq:
            if idx1_seq == idx2_seq:
                mcounts[idx1_seq] += 1
                mfiles[idx1_seq][0].write(h1 + s1 + p1 + q1)
                mfiles[idx1_seq][1].write(h2 + s2 + p2 + q2)
            else:
                # hopped counts
                inner = hop_counts.get(idx1_seq)
                if inner is None:
                    inner = {}
                    hop_counts[idx1_seq] = inner
                inner[idx2_seq] = inner.get(idx2_seq, 0) + 1

                hR1.write(h1 + s1 + p1 + q1)
                hR2.write(h2 + s2 + p2 + q2)
        else:
            unknown_count += 1
            uR1.write(h1 + s1 + p1 + q1)
            uR2.write(h2 + s2 + p2 + q2)

# CLOSE OUTPUTS
for idx in matched_files:
    matched_files[idx][0].close()
    matched_files[idx][1].close()
hopped_R1.close(); hopped_R2.close()
unknown_R1.close(); unknown_R2.close()

# COMPUTE STATS
matched_total = sum(matched_counts.values())
hopped_total = 0
for i1 in hopped_counts:
    inner = hopped_counts[i1]
    for i2 in inner:
        hopped_total += inner[i2]

matched_pct_total = (matched_total / total_reads * 100.0) if total_reads else 0.0
hopped_pct_total  = (hopped_total  / total_reads * 100.0) if total_reads else 0.0
unknown_pct_total = (unknown_count / total_reads * 100.0) if total_reads else 0.0

# Per-index matched (percent of TOTAL reads)
per_index_rows = []
for idx in known_indexes:
    c = matched_counts[idx]
    pct_total = (c / total_reads * 100.0) if total_reads else 0.0
    per_index_rows.append((idx, c, pct_total))
per_index_rows.sort(key=lambda x: x[1], reverse=True)

# Flatten hopped pairs for reporting/plotting
hopped_pairs = []
for i1 in known_indexes:
    inner = hopped_counts.get(i1)
    if inner:
        for i2 in known_indexes:
            if i2 != i1 and i2 in inner:
                c = inner[i2]
                p = (c / total_reads * 100.0) if total_reads else 0.0
                hopped_pairs.append((f"{i1}->{i2}", c, p))
hopped_pairs.sort(key=lambda x: x[1], reverse=True)

# WRITE MARKDOWN REPORT
report_path = f"{RESULTS_DIR}/summary.md"
with open(report_path, "w") as rep:
    rep.write("# Demultiplexing Summary\n\n")
    rep.write(f"- **Total read-pairs processed:** {total_reads}\n")
    rep.write(f"- **Matched (total):** {matched_total} ({matched_pct_total:.2f}%)\n")
    rep.write(f"- **Index-hopped (total):** {hopped_total} ({hopped_pct_total:.2f}%)\n")
    rep.write(f"- **Unknown (total):** {unknown_count} ({unknown_pct_total:.2f}%)\n\n")

    rep.write("## Percentage of reads from each sample (matched only)\n\n")
    rep.write("| Sample Index | Matched Count | % of Total Reads |\n")
    rep.write("|--------------|---------------:|-----------------:|\n")
    for idx, c, p in per_index_rows:
        rep.write(f"| {idx}-{idx} | {c} | {p:.2f}% |\n")

    rep.write("\n## Index hopping per pair (as % of total reads)\n\n")
    rep.write("| From -> To | Hopped Count | % of Total Reads |\n")
    rep.write("|------------|--------------:|-----------------:|\n")
    for tag, c, p in hopped_pairs:
        rep.write(f"| {tag} | {c} | {p:.2f}% |\n")

    rep.write("\n## Figures\n")
    rep.write("- `category_breakdown.png`: Matched vs Hopped vs Unknown counts\n")
    rep.write("- `matched_counts.png`: Matched reads per sample index\n")
    rep.write("- `hopped_top15.png`: Top 15 hopped pairs by count\n")

# PLOTS
# 1) Category breakdown
plt.figure(figsize=(6,4))
plt.bar(["Matched","Hopped","Unknown"], [matched_total, hopped_total, unknown_count])
plt.title("Read Categorization")
plt.ylabel("Read-pair count")
plt.tight_layout()
plt.savefig(f"{RESULTS_DIR}/category_breakdown.png", dpi=150)
plt.close()

# 2) Matched counts per index (nonzero only)
labels = [idx for idx, c, _ in per_index_rows if c > 0]
values = [c   for idx, c, _ in per_index_rows if c > 0]
if values:
    plt.figure(figsize=(10, max(3, 0.25*len(labels)+2)))
    plt.bar(labels, values)
    plt.title("Matched reads per sample index")
    plt.xlabel("Index (i7 == i5)")
    plt.ylabel("Matched read-pair count")
    plt.xticks(rotation=60, ha="right")
    plt.tight_layout()
    plt.savefig(f"{RESULTS_DIR}/matched_counts.png", dpi=150)
    plt.close()

# 3) Top 15 hopped pairs
topN = 15 if len(hopped_pairs) > 15 else len(hopped_pairs)
if topN > 0:
    hp_labels = [hopped_pairs[i][0] for i in range(topN)]
    hp_values = [hopped_pairs[i][1] for i in range(topN)]
    plt.figure(figsize=(10, max(3, 0.25*topN+2)))
    plt.bar(hp_labels, hp_values)
    plt.title("Top hopped pairs")
    plt.xlabel("From index1 -> To index2 (RC of i2)")
    plt.ylabel("Hopped read-pair count")
    plt.xticks(rotation=60, ha="right")
    plt.tight_layout()
    plt.savefig(f"{RESULTS_DIR}/hopped_top15.png", dpi=150)
    plt.close()

# STDOUT SNAPSHOT
print("\nDemultiplexing complete.")
print(f"Total pairs: {total_reads}")
print(f"Matched total: {matched_total}")
print(f"Hopped total: {hopped_total}")
print(f"Unknown total: {unknown_count}")
print(f"Report: {report_path}")
print(f"Figures written to {RESULTS_DIR}")
