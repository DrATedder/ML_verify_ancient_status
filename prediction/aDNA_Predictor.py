# ============================================================
# ANCIENT DNA RANDOM FOREST PREDICTOR
# HIGH-PERFORMANCE MULTIPROCESSING VERSION
# ============================================================
#
# FEATURES
# --------
#
# - FASTA / FASTQ support
# - GZIP FASTQ support
# - Parallel encoding
# - Chunked processing
# - Low-memory operation
# - Multi-core CPU utilisation
# - Designed for HPC / SLURM environments
# - Streaming FASTQ reading
#
# ============================================================

import os
import re
import glob
import gzip
import joblib
import multiprocessing as mp

import numpy as np
import pandas as pd

from Bio import SeqIO

# ============================================================
# USER CONFIG
# ============================================================

input_dir = ""

model_path = (
    "~/prediction/"
    "random_forest_combined.joblib"
)

output_dir = os.path.join(
    input_dir,
    "predictions"
)

os.makedirs(output_dir, exist_ok=True)

# ------------------------------------------------------------
# PERFORMANCE SETTINGS
# ------------------------------------------------------------

N_CORES = 32
CHUNK_SIZE = 50000

# ============================================================
# ONE-HOT ENCODING
# ============================================================

base_to_onehot = {

    "A": [1, 0, 0, 0, 0],
    "C": [0, 1, 0, 0, 0],
    "G": [0, 0, 1, 0, 0],
    "T": [0, 0, 0, 1, 0],
    "N": [0, 0, 0, 0, 1]
}

# ============================================================
# FORMAT DETECTION
# ============================================================

def detect_file_format(filepath):

    filename = filepath.lower()

    if filename.endswith((
        ".fastq",
        ".fq",
        ".fastq.gz",
        ".fq.gz"
    )):
        return "fastq"

    elif filename.endswith((
        ".fasta",
        ".fa",
        ".fna"
    )):
        return "fasta"

    else:
        raise ValueError(
            f"Unsupported format:\n{filepath}"
        )

# ============================================================
# FILE OPENING
# ============================================================

def open_sequence_file(filepath):

    if filepath.endswith(".gz"):

        return gzip.open(filepath, "rt")

    else:

        return open(filepath, "r")

# ============================================================
# ENCODING FUNCTION
# ============================================================

def encode_pair(pair_data):

    r1_seq, r2_seq, read_id = pair_data

    r1_encoded = np.array([
        base_to_onehot.get(
            b.upper(),
            [0, 0, 0, 0, 1]
        )
        for b in r1_seq
    ])

    r2_encoded = np.array([
        base_to_onehot.get(
            b.upper(),
            [0, 0, 0, 0, 1]
        )
        for b in r2_seq
    ])

    combined = np.concatenate(
        (r1_encoded, r2_encoded),
        axis=0
    )

    combined = combined.reshape(-1)

    return combined, read_id

# ============================================================
# CHUNK GENERATOR
# ============================================================

def generate_read_chunks(
    r1_path,
    r2_path,
    chunk_size=50000
):

    format1 = detect_file_format(r1_path)
    format2 = detect_file_format(r2_path)

    with open_sequence_file(r1_path) as h1, \
         open_sequence_file(r2_path) as h2:

        r1_iter = SeqIO.parse(h1, format1)
        r2_iter = SeqIO.parse(h2, format2)

        chunk = []

        for r1, r2 in zip(r1_iter, r2_iter):

            chunk.append((
                str(r1.seq),
                str(r2.seq),
                r1.id
            ))

            if len(chunk) >= chunk_size:

                yield chunk
                chunk = []

        if chunk:
            yield chunk

# ============================================================
# PREDICTION FUNCTION
# ============================================================

def process_dataset(
    r1_path,
    r2_path,
    model,
    output_csv
):

    print("\n===================================")
    print(f"PROCESSING:")
    print(os.path.basename(r1_path))
    print("===================================")

    all_results = []

    pool = mp.Pool(N_CORES)

    total_reads = 0

    for i, chunk in enumerate(

        generate_read_chunks(
            r1_path,
            r2_path,
            CHUNK_SIZE
        )

    ):

        print(
            f"\nChunk {i+1}"
            f" | Reads: {len(chunk)}"
        )

        # ----------------------------------------------------
        # PARALLEL ENCODING
        # ----------------------------------------------------

        encoded_results = pool.map(
            encode_pair,
            chunk
        )

        X_chunk = np.array([
            x[0]
            for x in encoded_results
        ])

        read_ids = [
            x[1]
            for x in encoded_results
        ]

        # ----------------------------------------------------
        # PREDICTION
        # ----------------------------------------------------

        predictions = model.predict(X_chunk)

        probabilities = model.predict_proba(
            X_chunk
        )

        ancient_prob = probabilities[:, 1]
        modern_prob = probabilities[:, 0]

        labels = np.where(
            predictions == 1,
            "Ancient",
            "Modern"
        )

        chunk_df = pd.DataFrame({

            "Read_ID":
                read_ids,

            "Ancient_Probability":
                ancient_prob,

            "Modern_Probability":
                modern_prob,

            "Predicted_Class":
                labels
        })

        all_results.append(chunk_df)

        total_reads += len(chunk_df)

        print(
            f"Processed total reads:"
            f" {total_reads:,}"
        )

    pool.close()
    pool.join()

    # --------------------------------------------------------
    # CONCATENATE RESULTS
    # --------------------------------------------------------

    results_df = pd.concat(
        all_results,
        ignore_index=True
    )

    results_df = results_df.sort_values(
        "Ancient_Probability",
        ascending=False
    )

    # --------------------------------------------------------
    # SAVE
    # --------------------------------------------------------

    results_df.to_csv(
        output_csv,
        index=False
    )

    # --------------------------------------------------------
    # SUMMARY
    # --------------------------------------------------------

    n_ancient = np.sum(
        results_df["Predicted_Class"]
        == "Ancient"
    )

    n_modern = np.sum(
        results_df["Predicted_Class"]
        == "Modern"
    )

    print("\n===================================")
    print("SUMMARY")
    print("===================================")

    print(f"Total Reads : {len(results_df):,}")
    print(f"Ancient     : {n_ancient:,}")
    print(f"Modern      : {n_modern:,}")

    print(
        f"\nSaved:\n{output_csv}"
    )

# ============================================================
# FIND INPUT FILES
# ============================================================

def find_r1_files(directory):

    patterns = [

        "*_R1.fastq",
        "*_R1.fastq.gz",

        "*_R1.fq",
        "*_R1.fq.gz",

        "*_R1.fasta",
        "*_R1.fa",
        "*_R1.fna"
    ]

    files = []

    for pattern in patterns:

        files.extend(
            glob.glob(
                os.path.join(
                    directory,
                    pattern
                )
            )
        )

    return sorted(files)

# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    print("\n===================================")
    print("LOADING MODEL")
    print("===================================")

    rf_model = joblib.load(model_path)

    print(f"\nLoaded:\n{model_path}")

    r1_files = find_r1_files(input_dir)

    if len(r1_files) == 0:

        raise ValueError(
            "No R1 files found."
        )

    print("\n===================================")
    print("FOUND DATASETS")
    print("===================================")

    for f in r1_files:
        print(os.path.basename(f))

    # --------------------------------------------------------
    # PROCESS EACH DATASET
    # --------------------------------------------------------

    for r1_path in r1_files:

        r2_path = re.sub(
            r"_R1",
            "_R2",
            r1_path
        )

        if not os.path.exists(r2_path):

            print(
                f"\nMissing R2 for:\n"
                f"{r1_path}"
            )

            continue

        dataset_name = re.sub(
            r"_R1.*",
            "",
            os.path.basename(r1_path)
        )

        output_csv = os.path.join(
            output_dir,
            f"{dataset_name}_predictions.csv"
        )

        process_dataset(
            r1_path,
            r2_path,
            rf_model,
            output_csv
        )

    print("\n===================================")
    print("ALL DATASETS COMPLETE")
    print("===================================")
