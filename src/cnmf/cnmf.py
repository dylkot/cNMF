#!/usr/bin/env python

import numpy as np
import pandas as pd
import os, errno
import datetime
import uuid
import itertools
import yaml
import subprocess
import scipy.sparse as sp

from scipy.spatial.distance import squareform
from sklearn.decomposition import non_negative_factorization
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.utils import sparsefuncs

from fastcluster import linkage
from scipy.cluster.hierarchy import leaves_list

import matplotlib.pyplot as plt

import scanpy as sc

def save_df_to_npz(obj, filename):
    np.savez_compressed(filename, data=obj.values, index=obj.index.values, columns=obj.columns.values)

def save_df_to_text(obj, filename):
    obj.to_csv(filename, sep='\t')

def load_df_from_npz(filename):
    with np.load(filename, allow_pickle=True) as f:
        obj = pd.DataFrame(**f)
    return obj

def check_dir_exists(path):
    """
    Checks if directory already exists or not and creates it if it doesn't
    """
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def worker_filter(iterable, worker_index, total_workers):
    return (p for i,p in enumerate(iterable) if (i-worker_index)%total_workers==0)

def fast_euclidean(mat):
    D = mat.dot(mat.T)
    squared_norms = np.diag(D).copy()
    D *= -2.0
    D += squared_norms.reshape((-1,1))
    D += squared_norms.reshape((1,-1))
    D = np.sqrt(D)
    D[D < 0] = 0
    return squareform(D, checks=False)

def fast_ols_all_cols(X, Y):
    pinv = np.linalg.pinv(X)
    beta = np.dot(pinv, Y)
    return(beta)

def fast_ols_all_cols_df(X,Y):
    beta = fast_ols_all_cols(X, Y)
    beta = pd.DataFrame(beta, index=X.columns, columns=Y.columns)
    return(beta)

def var_sparse_matrix(X):
    mean = np.array(X.mean(axis=0)).reshape(-1)
    Xcopy = X.copy()
    Xcopy.data **= 2
    var = np.array(Xcopy.mean(axis=0)).reshape(-1) - (mean**2)
    return(var)


def get_highvar_genes_sparse(expression, expected_fano_threshold=None,
                       minimal_mean=0.5, numgenes=None):
    # Find high variance genes within those cells
    gene_mean = np.array(expression.mean(axis=0)).astype(float).reshape(-1)
    E2 = expression.copy(); E2.data **= 2; gene2_mean = np.array(E2.mean(axis=0)).reshape(-1)
    gene_var = pd.Series(gene2_mean - (gene_mean**2))
    del(E2)
    gene_mean = pd.Series(gene_mean)
    gene_fano = gene_var / gene_mean

    # Find parameters for expected fano line
    top_genes = gene_mean.sort_values(ascending=False)[:20].index
    A = (np.sqrt(gene_var)/gene_mean)[top_genes].min()
    
    w_mean_low, w_mean_high = gene_mean.quantile([0.10, 0.90])
    w_fano_low, w_fano_high = gene_fano.quantile([0.10, 0.90])
    winsor_box = ((gene_fano > w_fano_low) &
                    (gene_fano < w_fano_high) &
                    (gene_mean > w_mean_low) &
                    (gene_mean < w_mean_high))
    fano_median = gene_fano[winsor_box].median()
    B = np.sqrt(fano_median)

    gene_expected_fano = (A**2)*gene_mean + (B**2)
    fano_ratio = (gene_fano/gene_expected_fano)

    # Identify high var genes
    if numgenes is not None:
        highvargenes = fano_ratio.sort_values(ascending=False).index[:numgenes]
        high_var_genes_ind = fano_ratio.index.isin(highvargenes)
        T=None


    else:
        if not expected_fano_threshold:
            T = (1. + gene_fano[winsor_box].std())
        else:
            T = expected_fano_threshold

        high_var_genes_ind = (fano_ratio > T) & (gene_mean > minimal_mean)

    gene_counts_stats = pd.DataFrame({
        'mean': gene_mean,
        'var': gene_var,
        'fano': gene_fano,
        'expected_fano': gene_expected_fano,
        'high_var': high_var_genes_ind,
        'fano_ratio': fano_ratio
        })
    gene_fano_parameters = {
            'A': A, 'B': B, 'T':T, 'minimal_mean': minimal_mean,
        }
    return(gene_counts_stats, gene_fano_parameters)



def get_highvar_genes(input_counts, expected_fano_threshold=None,
                       minimal_mean=0.5, numgenes=None):
    # Find high variance genes within those cells
    gene_counts_mean = pd.Series(input_counts.mean(axis=0).astype(float))
    gene_counts_var = pd.Series(input_counts.var(ddof=0, axis=0).astype(float))
    gene_counts_fano = pd.Series(gene_counts_var/gene_counts_mean)

    # Find parameters for expected fano line
    top_genes = gene_counts_mean.sort_values(ascending=False)[:20].index
    A = (np.sqrt(gene_counts_var)/gene_counts_mean)[top_genes].min()

    w_mean_low, w_mean_high = gene_counts_mean.quantile([0.10, 0.90])
    w_fano_low, w_fano_high = gene_counts_fano.quantile([0.10, 0.90])
    winsor_box = ((gene_counts_fano > w_fano_low) &
                    (gene_counts_fano < w_fano_high) &
                    (gene_counts_mean > w_mean_low) &
                    (gene_counts_mean < w_mean_high))
    fano_median = gene_counts_fano[winsor_box].median()
    B = np.sqrt(fano_median)

    gene_expected_fano = (A**2)*gene_counts_mean + (B**2)

    fano_ratio = (gene_counts_fano/gene_expected_fano)

    # Identify high var genes
    if numgenes is not None:
        highvargenes = fano_ratio.sort_values(ascending=False).index[:numgenes]
        high_var_genes_ind = fano_ratio.index.isin(highvargenes)
        T=None


    else:
        if not expected_fano_threshold:
            T = (1. + gene_counts_fano[winsor_box].std())
        else:
            T = expected_fano_threshold

        high_var_genes_ind = (fano_ratio > T) & (gene_counts_mean > minimal_mean)

    gene_counts_stats = pd.DataFrame({
        'mean': gene_counts_mean,
        'var': gene_counts_var,
        'fano': gene_counts_fano,
        'expected_fano': gene_expected_fano,
        'high_var': high_var_genes_ind,
        'fano_ratio': fano_ratio
        })
    gene_fano_parameters = {
            'A': A, 'B': B, 'T':T, 'minimal_mean': minimal_mean,
        }
    return(gene_counts_stats, gene_fano_parameters)


def compute_tpm(input_counts):
    """
    Default TPM normalization
    """
    tpm = input_counts.copy()
    sc.pp.normalize_per_cell(tpm, counts_per_cell_after=1e6)
    return(tpm)


class cNMF():


    def __init__(self, output_dir=".", name=None):
        """
        Parameters
        ----------

        output_dir : path, optional (default=".")
            Output directory for analysis files.

        name : string, optional (default=None)
            A name for this analysis. Will be prefixed to all output files.
            If set to None, will be automatically generated from date (and random string).
        """

        self.output_dir = output_dir
        if name is None:
            now = datetime.datetime.now()
            rand_hash =  uuid.uuid4().hex[:6]
            name = '%s_%s' % (now.strftime("%Y_%m_%d"), rand_hash)
        self.name = name
        self.paths = None
        self._initialize_dirs()


    def _initialize_dirs(self):
        if self.paths is None:
            # Check that output directory exists, create it if needed.
            check_dir_exists(self.output_dir)
            check_dir_exists(os.path.join(self.output_dir, self.name))
            check_dir_exists(os.path.join(self.output_dir, self.name, 'cnmf_tmp'))

            self.paths = {
                'normalized_counts' : os.path.join(self.output_dir, self.name, 'cnmf_tmp', self.name+'.norm_counts.h5ad'),
                'nmf_replicate_parameters' :  os.path.join(self.output_dir, self.name, 'cnmf_tmp', self.name+'.nmf_params.df.npz'),
                'nmf_run_parameters' :  os.path.join(self.output_dir, self.name, 'cnmf_tmp', self.name+'.nmf_idvrun_params.yaml'),
                'nmf_genes_list' :  os.path.join(self.output_dir, self.name, self.name+'.overdispersed_genes.txt'),

                'tpm' :  os.path.join(self.output_dir, self.name, 'cnmf_tmp', self.name+'.tpm.h5ad'),
                'tpm_stats' :  os.path.join(self.output_dir, self.name, 'cnmf_tmp', self.name+'.tpm_stats.df.npz'),

                'iter_spectra' :  os.path.join(self.output_dir, self.name, 'cnmf_tmp', self.name+'.spectra.k_%d.iter_%d.df.npz'),
                'iter_usages' :  os.path.join(self.output_dir, self.name, 'cnmf_tmp', self.name+'.usages.k_%d.iter_%d.df.npz'),
                'merged_spectra': os.path.join(self.output_dir, self.name, 'cnmf_tmp', self.name+'.spectra.k_%d.merged.df.npz'),

                'local_density_cache': os.path.join(self.output_dir, self.name, 'cnmf_tmp', self.name+'.local_density_cache.k_%d.merged.df.npz'),
                'consensus_spectra': os.path.join(self.output_dir, self.name, 'cnmf_tmp', self.name+'.spectra.k_%d.dt_%s.consensus.df.npz'),
                'consensus_spectra__txt': os.path.join(self.output_dir, self.name, self.name+'.spectra.k_%d.dt_%s.consensus.txt'),
                'consensus_usages': os.path.join(self.output_dir, self.name, 'cnmf_tmp',self.name+'.usages.k_%d.dt_%s.consensus.df.npz'),
                'consensus_usages__txt': os.path.join(self.output_dir, self.name, self.name+'.usages.k_%d.dt_%s.consensus.txt'),

                'consensus_stats': os.path.join(self.output_dir, self.name, 'cnmf_tmp', self.name+'.stats.k_%d.dt_%s.df.npz'),

                'clustering_plot': os.path.join(self.output_dir, self.name, self.name+'.clustering.k_%d.dt_%s.png'),
                'gene_spectra_score': os.path.join(self.output_dir, self.name, 'cnmf_tmp', self.name+'.gene_spectra_score.k_%d.dt_%s.df.npz'),
                'gene_spectra_score__txt': os.path.join(self.output_dir, self.name, self.name+'.gene_spectra_score.k_%d.dt_%s.txt'),
                'gene_spectra_tpm': os.path.join(self.output_dir, self.name, 'cnmf_tmp', self.name+'.gene_spectra_tpm.k_%d.dt_%s.df.npz'),
                'gene_spectra_tpm__txt': os.path.join(self.output_dir, self.name, self.name+'.gene_spectra_tpm.k_%d.dt_%s.txt'),

                'k_selection_plot' :  os.path.join(self.output_dir, self.name, self.name+'.k_selection.png'),
                'k_selection_stats' :  os.path.join(self.output_dir, self.name, self.name+'.k_selection_stats.df.npz'),
            }


    def prepare(self, counts_fn, components, n_iter = 100, densify=False, tpm_fn=None, seed=None,
                         beta_loss='frobenius',num_highvar_genes=2000, genes_file=None):
        """
        Load input counts, reduce to high-variance genes, and variance normalize genes.
        Subsequently prepare file for distributing jobs over workers.


        Parameters
        ----------
        counts_fn : str
            Path to input counts matrix

        components : list or numpy array
            Values of K to run NMF for
            
        n_iter : integer, optional (defailt=100)
            Number of iterations for factorization. If several ``k`` are specified, this many
            iterations will be run for each value of ``k``.

        densify : boolean, optional (default=False)
            Convert sparse data to dense

        tpm_fn : str or None, optional (default=None)
            If provided, load tpm data from file. Otherwise will compute it from the counts file
            
        seed : int or None, optional (default=None)
            Seed for sklearn random state.
            
        beta_loss : str or None, optional (default='frobenius')

        num_highvar_genes : int or None, optional (default=2000)
            If provided and genes_file is None, will compute this many highvar genes to use for factorization
        
        genes_file : str or None, optional (default=None)
            If provided will load high-variance genes from a list of these genes
        """
        
        
        if counts_fn.endswith('.h5ad'):
            input_counts = sc.read(counts_fn)
        else:
            ## Load txt or compressed dataframe and convert to scanpy object
            if counts_fn.endswith('.npz'):
                input_counts = load_df_from_npz(counts_fn)
            else:
                input_counts = pd.read_csv(counts_fn, sep='\t', index_col=0)
                
            if densify:
                input_counts = sc.AnnData(X=input_counts.values,
                                       obs=pd.DataFrame(index=input_counts.index),
                                       var=pd.DataFrame(index=input_counts.columns))
            else:
                input_counts = sc.AnnData(X=sp.csr_matrix(input_counts.values),
                                       obs=pd.DataFrame(index=input_counts.index),
                                       var=pd.DataFrame(index=input_counts.columns))

                
        if sp.issparse(input_counts.X) & densify:
            input_counts.X = np.array(input_counts.X.todense())
 
        if tpm_fn is None:
            tpm = compute_tpm(input_counts)
            sc.write(self.paths['tpm'], tpm)
        elif tpm_fn.endswith('.h5ad'):
            subprocess.call('cp %s %s' % (tpm_fn, self.paths['tpm']), shell=True)
            tpm = sc.read(self.paths['tpm'])
        else:
            if tpm_fn.endswith('.npz'):
                tpm = load_df_from_npz(tpm_fn)
            else:
                tpm = pd.read_csv(tpm_fn, sep='\t', index_col=0)
            
            if densify:
                tpm = sc.AnnData(X=tpm.values,
                            obs=pd.DataFrame(index=tpm.index),
                            var=pd.DataFrame(index=tpm.columns)) 
            else:
                tpm = sc.AnnData(X=sp.csr_matrix(tpm.values),
                            obs=pd.DataFrame(index=tpm.index),
                            var=pd.DataFrame(index=tpm.columns)) 

            sc.write(self.paths['tpm'], tpm)
        
        if sp.issparse(tpm.X):
            gene_tpm_mean = np.array(tpm.X.mean(axis=0)).reshape(-1)
            gene_tpm_stddev = var_sparse_matrix(tpm.X)**.5
        else:
            gene_tpm_mean = np.array(tpm.X.mean(axis=0)).reshape(-1)
            gene_tpm_stddev = np.array(tpm.X.std(axis=0, ddof=0)).reshape(-1)
            
            
        input_tpm_stats = pd.DataFrame([gene_tpm_mean, gene_tpm_stddev],
             index = ['__mean', '__std']).T
        save_df_to_npz(input_tpm_stats, self.paths['tpm_stats'])
        
        if genes_file is not None:
            highvargenes = open(genes_file).read().rstrip().split('\n')
        else:
            highvargenes = None

        norm_counts = self.get_norm_counts(input_counts, tpm, num_highvar_genes=num_highvar_genes,
                                               high_variance_genes_filter=highvargenes)
        
        
        if norm_counts.X.dtype != np.float64:
            norm_counts.X = norm_counts.X.astype(np.float64)

        self.save_norm_counts(norm_counts)
        (replicate_params, run_params) = self.get_nmf_iter_params(ks=components, n_iter=n_iter, random_state_seed=seed, beta_loss=beta_loss)
        self.save_nmf_iter_params(replicate_params, run_params)
        
    
    def combine(self, components=None):
        run_params = load_df_from_npz(self.paths['nmf_replicate_parameters'])

        if type(components) is int:
            ks = [components]
        elif components is None:
            ks = sorted(set(run_params.n_components))
        else:
            ks = components

        for k in ks:
            self.combine_nmf(k)    
    
    
    
    def get_norm_counts(self, counts, tpm,
                         high_variance_genes_filter = None,
                         num_highvar_genes = None
                         ):
        """
        Parameters
        ----------

        counts : anndata.AnnData
            Scanpy AnnData object (cells x genes) containing raw counts. Filtered such that
            no genes or cells with 0 counts
        
        tpm : anndata.AnnData
            Scanpy AnnData object (cells x genes) containing tpm normalized data matching
            counts

        high_variance_genes_filter : np.array, optional (default=None)
            A pre-specified list of genes considered to be high-variance.
            Only these genes will be used during factorization of the counts matrix.
            Must match the .var index of counts and tpm.
            If set to None, high-variance genes will be automatically computed, using the
            parameters below.

        num_highvar_genes : int, optional (default=None)
            Instead of providing an array of high-variance genes, identify this many most overdispersed genes
            for filtering

        Returns
        -------

        normcounts : anndata.AnnData, shape (cells, num_highvar_genes)
            A counts matrix containing only the high variance genes and with columns (genes)normalized to unit
            variance

        """

        if high_variance_genes_filter is None:
            ## Get list of high-var genes if one wasn't provided
            if sp.issparse(tpm.X):
                (gene_counts_stats, gene_fano_params) = get_highvar_genes_sparse(tpm.X, numgenes=num_highvar_genes)  
            else:
                (gene_counts_stats, gene_fano_params) = get_highvar_genes(np.array(tpm.X), numgenes=num_highvar_genes)
                
            high_variance_genes_filter = list(tpm.var.index[gene_counts_stats.high_var.values])
                
        ## Subset out high-variance genes
        norm_counts = counts[:, high_variance_genes_filter]

        ## Scale genes to unit variance
        if sp.issparse(tpm.X):
            sc.pp.scale(norm_counts, zero_center=False)
            if np.isnan(norm_counts.X.data).sum() > 0:
                print('Warning NaNs in normalized counts matrix')                       
        else:
            norm_counts.X /= norm_counts.X.std(axis=0, ddof=1)
            if np.isnan(norm_counts.X).sum().sum() > 0:
                print('Warning NaNs in normalized counts matrix')                    
        
        ## Save a \n-delimited list of the high-variance genes used for factorization
        open(self.paths['nmf_genes_list'], 'w').write('\n'.join(high_variance_genes_filter))

        ## Check for any cells that have 0 counts of the overdispersed genes
        zerocells = norm_counts.X.sum(axis=1)==0
        if zerocells.sum()>0:
            examples = norm_counts.obs.index[zerocells]
            print('Warning: %d cells have zero counts of overdispersed genes. E.g. %s' % (zerocells.sum(), examples[0]))
            print('Consensus step may not run when this is the case')
        
        return(norm_counts)

    
    def save_norm_counts(self, norm_counts):
        self._initialize_dirs()
        sc.write(self.paths['normalized_counts'], norm_counts)

        
    def get_nmf_iter_params(self, ks, n_iter = 100,
                               random_state_seed = None,
                               beta_loss = 'kullback-leibler'):
        """
        Create a DataFrame with parameters for NMF iterations.


        Parameters
        ----------
        ks : integer, or list-like.
            Number of topics (components) for factorization.
            Several values can be specified at the same time, which will be run independently.

        n_iter : integer, optional (defailt=100)
            Number of iterations for factorization. If several ``k`` are specified, this many
            iterations will be run for each value of ``k``.

        random_state_seed : int or None, optional (default=None)
            Seed for sklearn random state.

        """

        if type(ks) is int:
            ks = [ks]

        # Remove any repeated k values, and order.
        k_list = sorted(set(list(ks)))

        n_runs = len(ks)* n_iter

        np.random.seed(seed=random_state_seed)
        nmf_seeds = np.random.randint(low=1, high=(2**32)-1, size=n_runs)

        replicate_params = []
        for i, (k, r) in enumerate(itertools.product(k_list, range(n_iter))):
            replicate_params.append([k, r, nmf_seeds[i]])
        replicate_params = pd.DataFrame(replicate_params, columns = ['n_components', 'iter', 'nmf_seed'])

        _nmf_kwargs = dict(
                        alpha_W=0.0,
                        alpha_H='same',
                        l1_ratio=0.0,
                        beta_loss=beta_loss,
                        solver='mu',
                        tol=1e-4,
                        max_iter=1000,
                        init='random'
                        )
        
        ## Coordinate descent is faster than multiplicative update but only works for frobenius
        if beta_loss == 'frobenius':
            _nmf_kwargs['solver'] = 'cd'

        return(replicate_params, _nmf_kwargs)


    def save_nmf_iter_params(self, replicate_params, run_params):
        self._initialize_dirs()
        save_df_to_npz(replicate_params, self.paths['nmf_replicate_parameters'])
        with open(self.paths['nmf_run_parameters'], 'w') as F:
            yaml.dump(run_params, F)


    def _nmf(self, X, nmf_kwargs):
        """
        Parameters
        ----------
        X : pandas.DataFrame,
            Normalized counts dataFrame to be factorized.

        nmf_kwargs : dict,
            Arguments to be passed to ``non_negative_factorization``

        """
        (usages, spectra, niter) = non_negative_factorization(X, **nmf_kwargs)

        return(spectra, usages)


    def factorize(self,
                worker_i=0, total_workers=1,
                ):
        """
        Iteratively run NMF with prespecified parameters.

        Use the `worker_i` and `total_workers` parameters for parallelization.

        Generic kwargs for NMF are loaded from self.paths['nmf_run_parameters'], defaults below::

            ``non_negative_factorization`` default arguments:
                alpha=0.0
                l1_ratio=0.0
                beta_loss='kullback-leibler'
                solver='mu'
                tol=1e-4,
                max_iter=200
                regularization=None
                init='random'
                random_state, n_components are both set by the prespecified self.paths['nmf_replicate_parameters'].


        Parameters
        ----------
        norm_counts : pandas.DataFrame,
            Normalized counts dataFrame to be factorized.
            (Output of ``normalize_counts``)

        run_params : pandas.DataFrame,
            Parameters for NMF iterations.
            (Output of ``prepare_nmf_iter_params``)

        """
        run_params = load_df_from_npz(self.paths['nmf_replicate_parameters'])
        norm_counts = sc.read(self.paths['normalized_counts'])
        _nmf_kwargs = yaml.load(open(self.paths['nmf_run_parameters']), Loader=yaml.FullLoader)

        jobs_for_this_worker = worker_filter(range(len(run_params)), worker_i, total_workers)
        for idx in jobs_for_this_worker:

            p = run_params.iloc[idx, :]
            print('[Worker %d]. Starting task %d.' % (worker_i, idx))
            _nmf_kwargs['random_state'] = p['nmf_seed']
            _nmf_kwargs['n_components'] = p['n_components']

            (spectra, usages) = self._nmf(norm_counts.X, _nmf_kwargs)
            spectra = pd.DataFrame(spectra,
                                   index=np.arange(1, _nmf_kwargs['n_components']+1),
                                   columns=norm_counts.var.index)
            save_df_to_npz(spectra, self.paths['iter_spectra'] % (p['n_components'], p['iter']))


    def combine_nmf(self, k, remove_individual_iterations=False):
        run_params = load_df_from_npz(self.paths['nmf_replicate_parameters'])
        print('Combining factorizations for k=%d.'%k)

        combined_spectra = None
        n_iter = sum(run_params.n_components==k)

        run_params_subset = run_params[run_params.n_components==k].sort_values('iter')
        spectra_labels = []

        for i,p in run_params_subset.iterrows():

            spectra = load_df_from_npz(self.paths['iter_spectra'] % (p['n_components'], p['iter']))
            if combined_spectra is None:
                combined_spectra = np.zeros((n_iter, k, spectra.shape[1]))
            combined_spectra[p['iter'], :, :] = spectra.values

            for t in range(k):
                spectra_labels.append('iter%d_topic%d'%(p['iter'], t+1))

        combined_spectra = combined_spectra.reshape(-1, combined_spectra.shape[-1])
        combined_spectra = pd.DataFrame(combined_spectra, columns=spectra.columns, index=spectra_labels)

        save_df_to_npz(combined_spectra, self.paths['merged_spectra']%k)
        return combined_spectra


    def consensus(self, k, density_threshold=0.5, local_neighborhood_size = 0.30,show_clustering = True,
                  skip_density_and_return_after_stats = False, close_clustergram_fig=False):
        merged_spectra = load_df_from_npz(self.paths['merged_spectra']%k)
        norm_counts = sc.read(self.paths['normalized_counts'])

        density_threshold_str = str(density_threshold)
        if skip_density_and_return_after_stats:
            density_threshold_str = '2'
        density_threshold_repl = density_threshold_str.replace('.', '_')
        n_neighbors = int(local_neighborhood_size * merged_spectra.shape[0]/k)

        # Rescale topics such to length of 1.
        l2_spectra = (merged_spectra.T/np.sqrt((merged_spectra**2).sum(axis=1))).T


        if not skip_density_and_return_after_stats:
            # Compute the local density matrix (if not previously cached)
            topics_dist = None
            if os.path.isfile(self.paths['local_density_cache'] % k):
                local_density = load_df_from_npz(self.paths['local_density_cache'] % k)
            else:
                #   first find the full distance matrix
                topics_dist = squareform(fast_euclidean(l2_spectra.values))
                #   partition based on the first n neighbors
                partitioning_order  = np.argpartition(topics_dist, n_neighbors+1)[:, :n_neighbors+1]
                #   find the mean over those n_neighbors (excluding self, which has a distance of 0)
                distance_to_nearest_neighbors = topics_dist[np.arange(topics_dist.shape[0])[:, None], partitioning_order]
                local_density = pd.DataFrame(distance_to_nearest_neighbors.sum(1)/(n_neighbors),
                                             columns=['local_density'],
                                             index=l2_spectra.index)
                save_df_to_npz(local_density, self.paths['local_density_cache'] % k)
                del(partitioning_order)
                del(distance_to_nearest_neighbors)

            density_filter = local_density.iloc[:, 0] < density_threshold
            l2_spectra = l2_spectra.loc[density_filter, :]

        kmeans_model = KMeans(n_clusters=k, n_init=10, random_state=1)
        kmeans_model.fit(l2_spectra)
        kmeans_cluster_labels = pd.Series(kmeans_model.labels_+1, index=l2_spectra.index)

        # Find median usage for each gene across cluster
        median_spectra = l2_spectra.groupby(kmeans_cluster_labels).median()

        # Normalize median spectra to probability distributions.
        median_spectra = (median_spectra.T/median_spectra.sum(1)).T

        # Compute the silhouette score
        stability = silhouette_score(l2_spectra.values, kmeans_cluster_labels, metric='euclidean')

        # Obtain the reconstructed count matrix by re-fitting the usage matrix and computing the dot product: usage.dot(spectra)
        refit_nmf_kwargs = yaml.load(open(self.paths['nmf_run_parameters']), Loader=yaml.FullLoader)
        refit_nmf_kwargs.update(dict(
                                    n_components = k,
                                    H = median_spectra.values,
                                    update_H = False
                                    ))
        
        _, rf_usages = self._nmf(norm_counts.X,
                                          nmf_kwargs=refit_nmf_kwargs)
        rf_usages = pd.DataFrame(rf_usages, index=norm_counts.obs.index, columns=median_spectra.index)
        rf_pred_norm_counts = rf_usages.dot(median_spectra)

        # Compute prediction error as a frobenius norm
        if sp.issparse(norm_counts.X):
            prediction_error = ((norm_counts.X.todense() - rf_pred_norm_counts)**2).sum().sum()
        else:
            prediction_error = ((norm_counts.X - rf_pred_norm_counts)**2).sum().sum()
        
        consensus_stats = pd.DataFrame([k, density_threshold, stability, prediction_error],
                    index = ['k', 'local_density_threshold', 'stability', 'prediction_error'],
                    columns = ['stats'])

        if skip_density_and_return_after_stats:
            return consensus_stats
        
        save_df_to_npz(median_spectra, self.paths['consensus_spectra']%(k, density_threshold_repl))
        save_df_to_npz(rf_usages, self.paths['consensus_usages']%(k, density_threshold_repl))
        save_df_to_npz(consensus_stats, self.paths['consensus_stats']%(k, density_threshold_repl))
        save_df_to_text(median_spectra, self.paths['consensus_spectra__txt']%(k, density_threshold_repl))
        save_df_to_text(rf_usages, self.paths['consensus_usages__txt']%(k, density_threshold_repl))

        # Compute gene-scores for each GEP by regressing usage on Z-scores of TPM
        tpm = sc.read(self.paths['tpm'])
        tpm_stats = load_df_from_npz(self.paths['tpm_stats'])
        
        if sp.issparse(tpm.X):
            norm_tpm = (np.array(tpm.X.todense()) - tpm_stats['__mean'].values) / tpm_stats['__std'].values
        else:
            norm_tpm = (tpm.X - tpm_stats['__mean'].values) / tpm_stats['__std'].values
        
        if norm_tpm.dtype != np.float64:
            norm_tpm = norm_tpm.astype(np.float64)
        
        usage_coef = fast_ols_all_cols(rf_usages.values, norm_tpm)
        usage_coef = pd.DataFrame(usage_coef, index=rf_usages.columns, columns=tpm.var.index)

        save_df_to_npz(usage_coef, self.paths['gene_spectra_score']%(k, density_threshold_repl))
        save_df_to_text(usage_coef, self.paths['gene_spectra_score__txt']%(k, density_threshold_repl))

        # Convert spectra to TPM units, and obtain results for all genes by running last step of NMF
        # with usages fixed and TPM as the input matrix
        norm_usages = rf_usages.div(rf_usages.sum(axis=1), axis=0)
        refit_nmf_kwargs.update(dict(
                                    H = norm_usages.T.values,
                                ))
        
        # Needed otherwise _nmf will crash because with inconsistent dtypes
        if tpm.X.dtype != np.float64:
            tpm.X = tpm.X.astype(np.float64)
        
        _, spectra_tpm = self._nmf(tpm.X.T, nmf_kwargs=refit_nmf_kwargs)
        spectra_tpm = pd.DataFrame(spectra_tpm.T, index=rf_usages.columns, columns=tpm.var.index)
        save_df_to_npz(spectra_tpm, self.paths['gene_spectra_tpm']%(k, density_threshold_repl))
        save_df_to_text(spectra_tpm, self.paths['gene_spectra_tpm__txt']%(k, density_threshold_repl))

        if show_clustering:
            if topics_dist is None:
                topics_dist = squareform(fast_euclidean(l2_spectra.values))
                # (l2_spectra was already filtered using the density filter)
            else:
                # (but the previously computed topics_dist was not!)
                topics_dist = topics_dist[density_filter.values, :][:, density_filter.values]


            spectra_order = []
            for cl in sorted(set(kmeans_cluster_labels)):

                cl_filter = kmeans_cluster_labels==cl

                if cl_filter.sum() > 1:
                    cl_dist = squareform(topics_dist[cl_filter, :][:, cl_filter])
                    cl_dist[cl_dist < 0] = 0 #Rarely get floating point arithmetic issues
                    cl_link = linkage(cl_dist, 'average')
                    cl_leaves_order = leaves_list(cl_link)

                    spectra_order += list(np.where(cl_filter)[0][cl_leaves_order])
                else:
                    ## Corner case where a component only has one element
                    spectra_order += list(np.where(cl_filter)[0])


            from matplotlib import gridspec
            import matplotlib.pyplot as plt

            width_ratios = [0.5, 9, 0.5, 4, 1]
            height_ratios = [0.5, 9]
            fig = plt.figure(figsize=(sum(width_ratios), sum(height_ratios)))
            gs = gridspec.GridSpec(len(height_ratios), len(width_ratios), fig,
                                    0.01, 0.01, 0.98, 0.98,
                                   height_ratios=height_ratios,
                                   width_ratios=width_ratios,
                                   wspace=0, hspace=0)

            dist_ax = fig.add_subplot(gs[1,1], xscale='linear', yscale='linear',
                                      xticks=[], yticks=[],xlabel='', ylabel='',
                                      frameon=True)

            D = topics_dist[spectra_order, :][:, spectra_order]
            dist_im = dist_ax.imshow(D, interpolation='none', cmap='viridis', aspect='auto',
                                rasterized=True)

            left_ax = fig.add_subplot(gs[1,0], xscale='linear', yscale='linear', xticks=[], yticks=[],
                xlabel='', ylabel='', frameon=True)
            left_ax.imshow(kmeans_cluster_labels.values[spectra_order].reshape(-1, 1),
                            interpolation='none', cmap='Spectral', aspect='auto',
                            rasterized=True)


            top_ax = fig.add_subplot(gs[0,1], xscale='linear', yscale='linear', xticks=[], yticks=[],
                xlabel='', ylabel='', frameon=True)
            top_ax.imshow(kmeans_cluster_labels.values[spectra_order].reshape(1, -1),
                              interpolation='none', cmap='Spectral', aspect='auto',
                                rasterized=True)


            hist_gs = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs[1, 3],
                                   wspace=0, hspace=0)

            hist_ax = fig.add_subplot(hist_gs[0,0], xscale='linear', yscale='linear',
                xlabel='', ylabel='', frameon=True, title='Local density histogram')
            hist_ax.hist(local_density.values, bins=np.linspace(0, 1, 50))
            hist_ax.yaxis.tick_right()

            xlim = hist_ax.get_xlim()
            ylim = hist_ax.get_ylim()
            if density_threshold < xlim[1]:
                hist_ax.axvline(density_threshold, linestyle='--', color='k')
                hist_ax.text(density_threshold  + 0.02, ylim[1] * 0.95, 'filtering\nthreshold\n\n', va='top')
            hist_ax.set_xlim(xlim)
            hist_ax.set_xlabel('Mean distance to k nearest neighbors\n\n%d/%d (%.0f%%) spectra above threshold\nwere removed prior to clustering'%(sum(~density_filter), len(density_filter), 100*(~density_filter).mean()))
            
            ## Add colorbar
            cbar_gs = gridspec.GridSpecFromSubplotSpec(8, 1, subplot_spec=hist_gs[1, 0],
                                   wspace=0, hspace=0)
            cbar_ax = fig.add_subplot(cbar_gs[4,0], xscale='linear', yscale='linear',
                xlabel='', ylabel='', frameon=True, title='Euclidean Distance')
            vmin = D.min().min()
            vmax = D.max().max()
            fig.colorbar(dist_im, cax=cbar_ax,
            ticks=np.linspace(vmin, vmax, 3),
            orientation='horizontal')
            
            
            #hist_ax.hist(local_density.values, bins=np.linspace(0, 1, 50))
            #hist_ax.yaxis.tick_right()            

            fig.savefig(self.paths['clustering_plot']%(k, density_threshold_repl), dpi=250)
            if close_clustergram_fig:
                plt.close(fig)


    def k_selection_plot(self, close_fig=False):
        '''
        Borrowed from Alexandrov Et Al. 2013 Deciphering Mutational Signatures
        publication in Cell Reports
        '''
        run_params = load_df_from_npz(self.paths['nmf_replicate_parameters'])
        stats = []
        for k in sorted(set(run_params.n_components)):

            stats.append(self.consensus(k, skip_density_and_return_after_stats=True,
                                        show_clustering=False, close_clustergram_fig=True).stats)

        stats = pd.DataFrame(stats)
        stats.reset_index(drop = True, inplace = True)

        save_df_to_npz(stats, self.paths['k_selection_stats'])

        fig = plt.figure(figsize=(6, 4))
        ax1 = fig.add_subplot(111)
        ax2 = ax1.twinx()


        ax1.plot(stats.k, stats.stability, 'o-', color='b')
        ax1.set_ylabel('Stability', color='b', fontsize=15)
        for tl in ax1.get_yticklabels():
            tl.set_color('b')
        #ax1.set_xlabel('K', fontsize=15)

        ax2.plot(stats.k, stats.prediction_error, 'o-', color='r')
        ax2.set_ylabel('Error', color='r', fontsize=15)
        for tl in ax2.get_yticklabels():
            tl.set_color('r')

        ax1.set_xlabel('Number of Components', fontsize=15)
        ax1.grid('on')
        plt.tight_layout()
        fig.savefig(self.paths['k_selection_plot'], dpi=250)
        if close_fig:
            plt.close(fig)


def main():
    """
    Example commands:

        output_dir="./cnmf_test/"


        python cnmf.py prepare --output-dir $output_dir \
           --name test --counts ./cnmf_test/test_data.df.npz \
           -k 6 7 8 9 --n-iter 5

        python cnmf.py factorize  --name test --output-dir $output_dir

        THis can be parallelized as such:

        python cnmf.py factorize  --name test --output-dir $output_dir --total-workers 2 --worker-index WORKER_INDEX (where worker_index starts with 0)

        python cnmf.py combine  --name test --output-dir $output_dir

        python cnmf.py consensus  --name test --output-dir $output_dir

    """

    import sys, argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('command', type=str, choices=['prepare', 'factorize', 'combine', 'consensus', 'k_selection_plot'])
    parser.add_argument('--name', type=str, help='[all] Name for analysis. All output will be placed in [output-dir]/[name]/...', nargs='?', default='cNMF')
    parser.add_argument('--output-dir', type=str, help='[all] Output directory. All output will be placed in [output-dir]/[name]/...', nargs='?', default='.')
    parser.add_argument('-c', '--counts', type=str, help='[prepare] Input (cell x gene) counts matrix as df.npz or tab delimited text file')
    parser.add_argument('-k', '--components', type=int, help='[prepare] Numper of components (k) for matrix factorization. Several can be specified with "-k 8 9 10"', nargs='+')
    parser.add_argument('-n', '--n-iter', type=int, help='[prepare] Numper of factorization replicates', default=100)
    parser.add_argument('--total-workers', type=int, help='[all] Total number of workers to distribute jobs to', default=1)
    parser.add_argument('--seed', type=int, help='[prepare] Seed for pseudorandom number generation', default=None)
    parser.add_argument('--genes-file', type=str, help='[prepare] File containing a list of genes to include, one gene per line. Must match column labels of counts matrix.', default=None)
    parser.add_argument('--numgenes', type=int, help='[prepare] Number of high variance genes to use for matrix factorization.', default=2000)
    parser.add_argument('--tpm', type=str, help='[prepare] Pre-computed (cell x gene) TPM values as df.npz or tab separated txt file. If not provided TPM will be calculated automatically', default=None)
    parser.add_argument('--beta-loss', type=str, choices=['frobenius', 'kullback-leibler', 'itakura-saito'], help='[prepare] Loss function for NMF.', default='frobenius')
    parser.add_argument('--densify', dest='densify', help='[prepare] Treat the input data as non-sparse', action='store_true', default=False) 
    parser.add_argument('--worker-index', type=int, help='[factorize] Index of current worker (the first worker should have index 0)', default=0)
    parser.add_argument('--local-density-threshold', type=float, help='[consensus] Threshold for the local density filtering. This string must convert to a float >0 and <=2', default=0.5)
    parser.add_argument('--local-neighborhood-size', type=float, help='[consensus] Fraction of the number of replicates to use as nearest neighbors for local density filtering', default=0.30)
    parser.add_argument('--show-clustering', dest='show_clustering', help='[consensus] Produce a clustergram figure summarizing the spectra clustering', action='store_true')

    args = parser.parse_args()

    cnmf_obj = cNMF(output_dir=args.output_dir, name=args.name)
    
    if args.command == 'prepare':
        cnmf_obj.prepare(args.counts, components=args.components, n_iter=args.n_iter, densify=args.densify,
                         tpm_fn=args.tpm, seed=args.seed, beta_loss=args.beta_loss,
                         num_highvar_genes=args.numgenes, genes_file=args.genes_file)

    elif args.command == 'factorize':
        cnmf_obj.factorize(worker_i=args.worker_index, total_workers=args.total_workers)

    elif args.command == 'combine':
        cnmf_obj.combine(components=args.components)

    elif args.command == 'consensus':
        run_params = load_df_from_npz(cnmf_obj.paths['nmf_replicate_parameters'])

        if type(args.components) is int:
            ks = [args.components]
        elif args.components is None:
            ks = sorted(set(run_params.n_components))
        else:
            ks = args.components

        for k in ks:
            merged_spectra = load_df_from_npz(cnmf_obj.paths['merged_spectra']%k)
            cnmf_obj.consensus(k, args.local_density_threshold, args.local_neighborhood_size, args.show_clustering,
                               close_clustergram_fig=True)

    elif args.command == 'k_selection_plot':
        cnmf_obj.k_selection_plot(close_fig=True)


if __name__=="__main__":
    main()