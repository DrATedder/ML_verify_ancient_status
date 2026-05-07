import os
import random
import gzip
import glob
from Bio import SeqIO

# ==== USER CONFIGURATION ====
num_datasets = 5                   # Number of datasets to create
reads_per_dataset = 10000          # Total reads per dataset (R1+R2)
ancient_ratio = 0.3                # Proportion of ancient reads (0.0–1.0)
input_dir = "/storage02/or-microbio/meta_sim_test_25"        # Directory with FASTQ files
output_dir = "/storage02/or-microbio/meta_sim_test_25/ML_data_sets"  # Output directory for FASTA files

ancient_glob = "aUoBsim*_RNT_R*.fastq*"  # Glob pattern for ancient samples
modern_glob = "UoBsim*_R*.fastq*"        # Glob pattern for modern samples
# =============================

def paired_fastqs(file_list):
    paired = {}
    for f in file_list:
        fname = os.path.basename(f)
        if "_R1" in fname:
            base = fname.split("_R1")[0]
            paired.setdefault(base, {})["1"] = fname
        elif "_R2" in fname:
            base = fname.split("_R2")[0]
            paired.setdefault(base, {})["2"] = fname
    return [v for v in paired.values() if "1" in v and "2" in v]

def open_fastq(filepath):
    return gzip.open(filepath, "rt") if filepath.endswith(".gz") else open(filepath, "r")

def sample_paired_reads(r1_path, r2_path, num_reads):
    with open_fastq(r1_path) as f1, open_fastq(r2_path) as f2:
        r1_records = list(SeqIO.parse(f1, "fastq"))
        r2_records = list(SeqIO.parse(f2, "fastq"))

    if len(r1_records) != len(r2_records):
        raise ValueError(f"Paired files {r1_path} and {r2_path} have unequal reads.")

    indices = random.sample(range(len(r1_records)), min(num_reads, len(r1_records)))
    return [r1_records[i] for i in indices], [r2_records[i] for i in indices]

def write_fasta(records, out_path):
    with open(out_path, "w") as f:
        for r in records:
            f.write(f">{r.id}\n{r.seq}\n")

# Prepare output directory
os.makedirs(output_dir, exist_ok=True)

# List FASTQ file pairs
ancient_files = glob.glob(os.path.join(input_dir, ancient_glob))
modern_files = glob.glob(os.path.join(input_dir, modern_glob))
ancient_pairs = paired_fastqs(ancient_files)
modern_pairs = paired_fastqs(modern_files)

if not ancient_pairs or not modern_pairs:
    raise ValueError(f"Missing input FASTQ pairs.\n"
                     f"Ancient found: {len(ancient_pairs)} pairs\n"
                     f"Modern found: {len(modern_pairs)} pairs\n"
                     f"Check glob patterns and input_dir.")

# Generate datasets
for i in range(1, num_datasets + 1):
    a_reads = int(reads_per_dataset * ancient_ratio / 2)
    m_reads = int(reads_per_dataset * (1 - ancient_ratio) / 2)

    # Sample ancient
    a_pair = random.choice(ancient_pairs)
    a_r1, a_r2 = sample_paired_reads(os.path.join(input_dir, a_pair["1"]),
                                     os.path.join(input_dir, a_pair["2"]),
                                     a_reads)

    # Sample modern
    m_pair = random.choice(modern_pairs)
    m_r1, m_r2 = sample_paired_reads(os.path.join(input_dir, m_pair["1"]),
                                     os.path.join(input_dir, m_pair["2"]),
                                     m_reads)

    # Combine and shuffle
    paired_combined = list(zip(a_r1 + m_r1, a_r2 + m_r2))
    random.shuffle(paired_combined)
    dataset_r1, dataset_r2 = zip(*paired_combined)

    # Write to FASTA
    write_fasta(dataset_r1, os.path.join(output_dir, f"dataset_{i}_R1.fasta"))
    write_fasta(dataset_r2, os.path.join(output_dir, f"dataset_{i}_R2.fasta"))
    print(f"Created dataset {i}: {len(dataset_r1)} paired reads")
