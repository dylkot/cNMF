# cNMF code and example data

> Running cNMF proceeds through a serious of steps starting from a filtered count matrix and ending with a set of estimated gene expression program spectra and their usages for each cell in the dataset. This can be a powerful approach to analyzing single-cell RNA-Seq data where some cells express activity programs in addition to their identity program. For now, if you use this code, please cite [the manuscript](https://www.biorxiv.org/content/early/2018/04/30/310599) on bioarxiv. In addition to this code for running cNMF, all of the analysis in the manuscript available for exploration and re-execution on [Code Ocean](https://codeocean.com/2018/11/20/identifying-gene-expression-programs-of-cell-type-identity-and-cellular-activity-with-single-cell-rna-seq/code). This includes examples of running cNMF on real data and downstream analyses.


# Make sure you have the appropriate packages installed.
> We reccomend [conda](https://conda.io/miniconda.html) as a package management system to install the necessary python packages. After installing and configuring conda, you can create an environment for this workflow as follows
```
conda update -n base conda # Make sure conda is up to date
conda create -n cnmf_env python=3.6
conda activate cnmf_env
conda install --yes --channel bioconda --channel conda-forge --channel defaults cython==0.28.2 fastcluster==1.1.24 jupyter==1.0.0 jupyterlab==0.32.1 matplotlib==2.2.2 numpy==1.12.1 palettable==3.1.1 pandas==0.23.0 scipy==1.1.0 && conda clean --yes --all
pip install git+https://github.com/scikit-learn/scikit-learn.git@3b24a4855db8564ddccb011228c9f4a51a519c9a
pip install --upgrade --no-cache-dir --upgrade-strategy=only-if-needed bhtsne==0.1.9 #(only needed to generate the tsne in the example)
```
   Now your just need to activate cnmf_env each time before running cNMF.

# Step by step guide 
> You can see all possible command line options by running
```
python cnmf.py --help
```

> And also see analyze_example_data.ipynb for a step by step walkthrough with example data. But we also discuss the individual steps below.

### Step 1 prepare the run parameters
> Initializes the run by creating a run parameters file and a normalized input matrix
  Example command:
  ```
   python ./cnmf.py prepare --output-dir ./example_data --name example_cNMF -c ./example_data/counts_prefiltered.txt -k 5 6 7 8 9 10 11 12 13 --n-iter 100 --total-workers 8 --seed 14 --numgenes 2000
  ```
  - --output-dir - the output directory into which all results will be placed
  - --name - a subdirectory output_dir/name will be created and all output files will have name as their prefix
  - -c - path to the counts file that will be used. This should be a raw counts matrix (I.e. not log-transformed or normalized per cell library size)
  - -k - list of K values that will be tested for cNMF
  - --n-iter - specifies the number of NMF iterations to run for each K
  - --total-workers - specifies how many workers (e.g. cores on a machine or nodes on a compute farm) can be used in parallel
  - --seed - the master seed that will be used to generate the individual seed for each NMF replicate
  - --numgenes - the number of higest variance genes that will be used for running the factorization. Removing low variance genes helps amplify the signal and is an important factor in correctly inferring programs in the data. However, don't worry, at the end the spectra is re-fit to include estimates for all genes, even those that weren't included in the high-variance set.
  
**Please note that the counts matrix must not contain any genes or cells with 0 counts. Those must be filtered out prior to running cNMF**

### Step 2 factorize the matrix using the parameters specified previously
> Next you run NMF for all of the replicates specified previously. The basic idea is that the NMF tasks have been divied up amongst the number of workers specified in step 1. You just need to specify the index of a given worker (from 0 to total-workers-1) to start its jobs like in the following:
  ```
   python ./cnmf.py factorize --output-dir ./example_data --name example_cNMF --worker-index 0 
  ```
This would run all the factorization steps for the first worker. See the analyze_example_data.ipynb notebook for examples of how you could submit all of the workers to run in parallel either      using [GNU parralel](https://www.gnu.org/software/parallel/) or an [UGER scheduler](http://www.univa.com/resources/files/univa_user_guide_univa__grid_engine_854.pdf). Tip: The implementation of
NMF in scikit-learn by default will use half of the available cores on a machine. Therefore, if you are using GNU parallel on a large machine, you should use 2 workers to get the best performance.
  
### Step 3 combine the individual spectra results files for each K into a merged file
> Since a separate file has been created for each replicate for each K, we combine the replicates for each K as below:
Example command:
  ```
  python ./cnmf.py combine --output-dir ./example_data --name example_cNMF
  ```
After this, you can optionally delete the individual spectra files like so:
  ```
  rm ./example_data/example_cNMF/cnmf_tmp/example_cNMF.spectra.k_*.iter_*.df.npz
  ```
  
### Step 4 select an optimal K by considering the trade-off between stability and error
> This will iterate through all of the values of K that have been run and will calculate the stability and error.
It then outputs a PDF file plotting this relationship into the output_dir/name directory
Example command:
```
  python ./cnmf.py k_selection_plot --output-dir ./example_data --name example_cNMF
```
This outputs a K selection plot to example_data/example_cNMF/example_cNMF.k_selection.pdf. There is no universally definitive criteria for choosing K but we will typically use the largest value that is reasonably stable and/or a local maximum in stability.



### Step 5 obtain consensus estimates for the programs and their usages at the desired value of K
> This will cluster the spectra after first optionally filtering out ouliers. It then collapses the clusters down to their median.
It finally recomputes 2 estimates of the programs (one in units of TPM, the other in units of TPM Z-scores, thus reflecting whether
having a higher usage of a program would be expected to decrease or increase gene expression), and an estimate of the usage. It also
outputs a diagnostic plot of the clustergram, and the histogram of distances between each spectra and its K nearest neighbors that you
can use for setting a threshold to filter outliers. We recommend inspecting the histogram to determine a threshold to use for filtering
the outliers. By default it determines the number of neighbors to use for this filtering as 30% of the number of iterations done. But this
can be modified from the command line. The --local-density-threshold determines the threshold on average distance to K nearest neighbors to
use. 2.0 or above means that nothing will be filtered out. 
Example command:
```
    python ./cnmf.py consensus --output-dir ./example_data --name example_cNMF --components 10 --local-density-threshold 0.01
```
By the end of this step, you should have the following files in your directory:
   - TPM Z-score unit spectra matrix - example_data/example_cNMF/example_cNMF.gene_spectra_score.k_10.dt_0_01.txt
   - TPM unit spectra matrix - example_data/example_cNMF/example_cNMF.gene_spectra_tpm.k_10.dt_0_01.txt
   - Usage matrix example_data/example_cNMF/example_cNMF.usages.k_10.dt_0_01.consensus.txt
   - Diagnostic plot - example_data/example_cNMF/example_cNMF.clustering.k_10.dt_0_01.pdf
    
    
