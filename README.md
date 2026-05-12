# ML_verify_ancient_status
This project aims to develop a ML model to identify modern contaminants in ancient metagenomic datasets (or to validate the 'ancient' status of each read, if you prefer). The model is currently a pretty basic RandomForest model (I know,  this isn't the best idea, but as a benchmark...), but the early stages scripts collected here will take you from simulating reads (using [InSilicoSeq](https://github.com/HadrienG/InSilicoSeq)), simulating aDNA damage patterns (using [simulate_damage](https://github.com/DrATedder/simulate_damage)), drawing datasets from your simulated reads. The final pre-model step is one-hot encoding your data sets.
There is a run order here, but this is **far** from a pipeline. It is worth noting that the `sbatch` scripts assume you have a `slurm` job submission system set up liket the one i've been using. They also assume you have the correct pre-requisites installed.

## Pre-requisite software requirements (probably not exhaustive...)
**Note**. I've given some token install instructions here. Best always to read the actual install instructions though, and adjust for your own system.
* `InSilicoSeq` - `conda install -c bioconda insilicoseq`
* `simulate_damage` - `pip3 install simulate_damage`
* `python3` - Best to find a suitable installer from [python.org](https://www.python.org/downloads/)

## Python module requirements (again, probably not exhaustive...)
```
numpy
pandas
matplotlib
sklearn.ensemble
sklearn.model_selection
sklearn.metrics
IPython.display
Bio
random
gzip
```
## Script run order

1. Simulate some reads - `simulate_reads_ncbi.sbatch`
2. Apply deamination damage to reads - `deamination_pipeline.sbatch`
3. Adjust read labels for 'ancient' reads - `tweak_read_names.py`
4. Create datasets with user defined proportion of 'ancient' and 'modern' reads - `data_Set_creation.py`
5. Pull metadata file for each pairwise dataset - `create_metadata.py`
6. One-hot encode the datasets - `convert_data_numpy.py`
7. Train the model, and guage precision and accuracy - `RF_ML_model_ancient_status.ipynb`

**Note** A `streamlit` app version of `RF_ML_model_ancient_status.ipynb` is available [here](https://github.com/DrATedder/RF_ancient_validation_streamlit).
