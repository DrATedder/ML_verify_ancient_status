import os
import numpy as np
from Bio import SeqIO

# === USER CONFIGURATION ===
input_dir = "/storage02/or-microbio/meta_sim_test_25/ML_data_sets"
output_dir = "/storage02/or-microbio/meta_sim_test_25/ML_data_sets/encoded_output"
# ===========================

# One-hot encoding map
BASES = "ACGTN"
base_to_onehot = {
    "A": [1, 0, 0, 0, 0],
    "C": [0, 1, 0, 0, 0],
    "G": [0, 0, 1, 0, 0],
    "T": [0, 0, 0, 1, 0],
    "N": [0, 0, 0, 0, 1]
}

def one_hot_encode_sequence(seq):
    return np.array([base_to_onehot.get(base.upper(), [0, 0, 0, 0, 1]) for base in seq])

def load_fasta_sequences(fasta_path):
    return list(SeqIO.parse(fasta_path, "fasta"))

def load_metadata(metadata_path, dataset_name):
    with open(metadata_path, "r") as f:
        header = f.readline().strip().split(",")[1:]  # skip 'dataset' column
        for line in f:
            parts = line.strip().split(",")
            if parts[0] == dataset_name:
                return list(map(int, parts[1:]))
    raise ValueError(f"No metadata found for {dataset_name} in {metadata_path}")

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Process each paired fasta set
def one_hot_encode_sequence(seq):
    return np.array([base_to_onehot.get(base.upper(), [0, 0, 0, 0, 1]) for base in seq])

def load_fasta_sequences(fasta_path):
    return list(SeqIO.parse(fasta_path, "fasta"))

def load_metadata(metadata_path, dataset_name):
    with open(metadata_path, "r") as f:
        header = f.readline().strip().split(",")[1:]  # skip 'dataset' column
        for line in f:
            parts = line.strip().split(",")
            if parts[0] == dataset_name:
                return list(map(int, parts[1:]))
    raise ValueError(f"No metadata found for {dataset_name} in {metadata_path}")

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Process each paired fasta set
for r1_file in sorted(f for f in os.listdir(input_dir) if f.endswith("_R1.fasta")):
    dataset_name = r1_file.replace("_R1.fasta", "")
    r2_file = dataset_name + "_R2.fasta"
    metadata_file = dataset_name + "_metadata.csv"

    r1_path = os.path.join(input_dir, r1_file)
    r2_path = os.path.join(input_dir, r2_file)
    metadata_path = os.path.join(input_dir, metadata_file)

    if not os.path.exists(r2_path):
        print(f"Missing R2 file for {dataset_name}, skipping.")
        continue
    if not os.path.exists(metadata_path):
        print(f"Missing metadata file {metadata_file}, skipping.")
        continue

    print(f"Processing dataset: {dataset_name}")

    r1_reads = load_fasta_sequences(r1_path)
    r2_reads = load_fasta_sequences(r2_path)

    if len(r1_reads) != len(r2_reads):
        print(f"Read count mismatch in {dataset_name}, skipping.")
        continue

    labels = load_metadata(metadata_path, dataset_name)
    if len(labels) != len(r1_reads):
        print(f"Metadata size mismatch in {dataset_name}, skipping.")
        continue

    # Encode each read pair and combine
    encoded = []
    for r1, r2 in zip(r1_reads, r2_reads):
        r1_encoded = one_hot_encode_sequence(str(r1.seq))
        r2_encoded = one_hot_encode_sequence(str(r2.seq))
        combined = np.concatenate((r1_encoded, r2_encoded), axis=0)
        encoded.append(combined)

    # Save as .npy
    X = np.stack(encoded)
    y = np.array(labels)

    np.save(os.path.join(output_dir, f"X_data_{dataset_name}.npy"), X)
    np.save(os.path.join(output_dir, f"y_labels_{dataset_name}.npy"), y)

    print(f"Saved {X.shape[0]} samples for {dataset_name}.")
	
