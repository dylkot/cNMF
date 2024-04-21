# Consensus NMF (cNMF)

<img src="https://storage.googleapis.com/sabeti-public/dkotliar/elife-cNMF-fig1.jpg" style="height: 800px;" />

cNMF is a pipeline for inferring gene expression programs from scRNA-Seq

It takes a count matrix (N cells X G genes) as input and produces a (K x G) matrix of gene expression programs (GEPs) and a (N x K) matrix specifying the usage of each program for each cell in the data.

Read more about the method in the [publication](https://elifesciences.org/articles/43803) and check out examples on [simulated data](Tutorials/analyze_simulated_example_data.ipynb) and [PBMCs](Tutorials/analyze_pbmc_example_data.ipynb).

# Installation
cNMF has been tested with Python 3.7 and requires scikit-learn>=1.0 and scanpy>=1.8

You can install with [pip](https://pypi.org/):

```bash
pip install cnmf
```

We also suggest creating a [conda](https://conda.io/miniconda.html) environment to manage the dependencies as follows:

```bash
conda update -yn base conda
conda create -n cnmf_env --yes --channel bioconda --channel conda-forge --channel defaults python=3.7 fastcluster matplotlib numpy palettable pandas scipy 'scikit-learn>=1.0' pyyaml 'scanpy>=1.8' && conda clean --yes --all # Create environment, cnmf_env, containing required packages
conda activate cnmf_env # Activate cnmf_env - necessary before running cnmf
pip install cnmf # install the actual cnmf package
    
## Only needed to load the example notebook in jupyterlab but not needed for non-interactive runs ## 
conda install --yes jupyterlab && conda clean --yes --all
```

Now you can run cNMF as long as the cnmv_enf environment is activated with the `conda activate cnmf_env` command

# Running cNMF

cNMF can be run from the command line without any parallelization using the example commands below:

```bash
cnmf prepare --output-dir ./example_data --name example_cNMF -c ./example_data/counts_prefiltered.txt -k 5 6 7 8 9 10 11 12 13 --n-iter 100 --seed 14

cnmf factorize --output-dir ./example_data --name example_cNMF --worker-index 0 --total-workers 1

cnmf combine --output-dir ./example_data --name example_cNMF

cnmf k_selection_plot --output-dir ./example_data --name example_cNMF

cnmf consensus --output-dir ./example_data --name example_cNMF --components 10 --local-density-threshold 0.01 --show-clustering
```

Or alternatively, the same steps can be run from within a Python environment using the commands below:

```python
from cnmf import cNMF
cnmf_obj = cNMF(output_dir="./example_data", name="example_cNMF")
cnmf_obj.prepare(counts_fn="./example_data/counts_prefiltered.txt", components=np.arange(5,14), n_iter=100, seed=14)
cnmf_obj.factorize(worker_i=0, total_workers=1)
cnmf_obj.combine()
cnmf_obj.k_selection_plot()
cnmf_obj.consensus(k=10, density_threshold=0.01)
usage, spectra_scores, spectra_tpm, top_genes = cnmf_obj.load_results(K=10, density_threshold=0.01)
```

For the Python environment approach, `usage` will contain the usage matrix with each cell normalized to sum to 1. `spectra_scores` contains the gene_spectra_scores output (aka Z-score unit gene expression matrix), `spectra_tpm` contains the GEP spectra in units of TPM and `top_genes` contains an ordered list of the top 100 associated genes for each program.

Output data files will all be available in the ./example_data/example_cNMF directory including:

   - Z-score unit gene expression program matrix - `example_data/example_cNMF/example_cNMF.gene_spectra_score.k_10.dt_0_01.txt`
   - TPM unit gene expression program matrix - `example_data/example_cNMF/example_cNMF.gene_spectra_tpm.k_10.dt_0_01.txt`
   - Usage matrix - `example_data/example_cNMF/example_cNMF.usages.k_10.dt_0_01.consensus.txt`
   - K selection plot - `example_data/example_cNMF/example_cNMF.k_selection.png`
   - Clustergram diagnostic plot - `example_data/example_cNMF/example_cNMF.clustering.k_10.dt_0_01.pdf`

Some usage notes:
 - __Parallelization__: The factorize step can be parallelized with the --total-workers flag and then submitting multiple jobs, one per worker, indexed starting by 0. For example:
  ```
  cnmf factorize --output-dir ./example_data --name example_cNMF --worker-index 0 --total-workers 3 &
  cnmf factorize --output-dir ./example_data --name example_cNMF --worker-index 1 --total-workers 3 &
  cnmf factorize --output-dir ./example_data --name example_cNMF --worker-index 2 --total-workers 3 &
  ```
  would break the factorization jobs up into 3 batches and submit them independently. This can be used with compute clusters to run the factorizations in parallel (see tutorials for example).
 - __Input data__: Input data can be provided in 2 ways:
    - 1. as a raw tab-delimited text file containing row labels with cell IDs (barcodes) and column labels as gene IDs
    - 2. as a scanpy file ending in .h5ad containg counts as the data feature. See the PBMC dataset tutorial for an example of how to generate the Scanpy object from the data provided by 10X. Because Scanpy uses sparse matrices by default, the .h5ad data structure can take up much less memory than the raw counts matrix and can be much faster to load. 
    
See the tutorials or Stepwise_Guide.md for more details

# Change log

### New in version 1.5
- Fixed bug in detecting and printing cells with 0 counts of overdispersed genes
- Added option in load_results() to return normalized or unnormalized usage.
- Added a Preprocess class to batch correct data prior to cNMF. See the added Tutorial analyze_batcheffectcorrect_BaronEtAl.ipynb to illustrate its basic usage.

### New in version 1.4
- Usage is re-fit a final time from gene_spectra_tpm which increases accuracy in simulations
- Use cnmf_obj.load_results(K=_, density_threshold=_) to obtain usage, spectra_scores, spectra_tpm, and top_genes matrices
- cnmf_obj.combine() now has a skip_missing_files=True/False option to skip incomplete factorize iterations
- GEPs are now ordered by maximum total usage
- Now detects and errors when 0 counts of HVGs with interpretable error message

### New in version 1.3
- Installation via pip
- Object oriented interface for Python users and command line script option via `cnmf`  

### New in version 1.2
 - Increased the threshold for ignoring genes with low mean expression for determining high-variance genes from a TPM of 0.01 to 0.5. Some users were identifying uninterpretable programs with very low usage except in a tiny number of cells. We suspect that this was due to including genes as high-variance that are detected in a small number of cells. This change in the default parameter will help offset that problem in most cases.
 - Updated import of NMF for compatibility with scikit-learn versions >22
 - Colorbar for heatmaps included with consensus matrix plot

### New in version 1.1
 - Now operates by default on sparse matrices. Use --densify option in prepare step if data is not sparse
 - Now takes Scanpy AnnData object files (.h5ad) as input
 - Now has option to use KL divergence beta_loss instead of Frobenius. Frobenius is the default because it is much faster.
 - Includes a Docker file for creating a Docker container to run cNMF in parallel with cloud compute
 - Includes a tutorial on a simple PBMC dataset
 - Other minor fixes
