import os
import pandas as pd

class cNMFWrapper():

    """
    Wrapper around the cNMF.py command line interface
    for easy use inside python code

    to run it:
    1) prepare()
    2) factorize()
    3) consensus()
    4) load_results()
    """
    def __init__(self, output_dir:str, run_name:str, cnmf_executable:str, n_workers:int, verbose=False):
        """
        Parameters
        ----------

        output_dir : path
            folder to save results to

        run_name : string
            A name for this analysis. Will be prefixed to all output files.
            Creates a subfolder with this name in output_dir, storing some run-specific files

        run_name : string
            The command line call to cnmf.py, something like `python /some/path/cnmf.py`

        n_workers : int
            Number of processes to run in parallel
        """

        self.output_dir = output_dir
        self.run_name = run_name
        self.executable = cnmf_executable

        # make sure the executable exists
        assert self.executable.startswith('python')
        assert os.path.exists(self.executable.split()[1]), f'executable doesnt exist {self.executable.split()[1]}'
        self.n_workers = n_workers
        self.verbose = verbose

    def _create_call_str(self, mode:str):
        """
        forming the call of the python executable cnmf.py <mode> --output-dir <> --name
        based on the objects output_dir and run_name
        """
        s = f'{self.executable} {mode} --output-dir {self.output_dir} --name {self.run_name}'
        return s

    def prepare(self, h5ad_file:str, k_list:list, num_iter:int, num_genes:int):
        """
        From a raw count matrix (AnnData in .h5ad) to properly normalized data for NMF

        Parameters
        ----------

        h5ad_file : anndata.AnnData
            Path to a Scanpy AnnData object (cells x genes) containing raw counts.
            Filtered such that no genes or cells with 0 counts

        k_list : list<int>
            List of numbers of components to calculate NMF for

        num_genes : int
            Identify this many most overdispersed genes for filtering

        num_iter : int
            Number of times to run NMF for a given parameter setting.
            Since NMF is stochastic, results may vary and considering all repeats
            will make the analysis more robust
        """
        assert isinstance(h5ad_file, str), "should be a filename, not and AnnData object"
        K = ' '.join([str(i) for i in k_list])
        call = self._create_call_str('prepare')
        prepare_cmd = f'{call} -c {h5ad_file}  -k {K} --n-iter {num_iter} --total-workers {self.n_workers} --numgenes {num_genes}'
        if self.verbose:
            print(f'Prepare command ({self.n_workers} Workers):\n{prepare_cmd}')
        os.system(prepare_cmd)

    def factorize(self):
        """
        Actual cNMF and combination of results of different repeats/iterations
        """
        call = self._create_call_str('factorize')

        if self.n_workers == 1:
            factorize_cmd = f'{call}  --worker-index 0'
            if self.verbose:
                print('Factorize command for 1 worker :\n%s' % factorize_cmd)
        else:
            worker_index = ' '.join([str(x) for x in range(self.n_workers)])
            factorize_cmd = f'nohup parallel {call} --worker-index {{}} ::: {worker_index}'
            if self.verbose:
                print(f'Factorize command for {self.n_workers} workers:\n{factorize_cmd}')

        os.system(factorize_cmd)


        combine_cmd = self._create_call_str('combine')
        if self.verbose:
            print('Combine command:\n%s' % combine_cmd)
        os.system(combine_cmd)

    def k_selection_plot(self):
        """
        For each number of topics (k), the error (averaged over repeats) and the
        robustness (across repeats) is plotted
        """
        plot_K_selection_cmd = self._create_call_str('k_selection_plot')
        if self.verbose:
            print('Plot K tradeoff command:\n%s' % plot_K_selection_cmd)
        os.system(plot_K_selection_cmd)

        img_name = f"{self.output_dir}/{self.run_name}/{self.run_name}.k_selection.png"
        from IPython.display import Image
        I = Image(filename=img_name, width=1000, height=1000)
        return I

        # import Image
        # image = Image.open(img_name)
        # image.show()

    def consensus(self, selected_K:int, local_density_threshold=2):
        """
        for a specific k (number of topics), merge the different repeats into
        consensus topics/usages

        Parameters
        ----------

        selected_K : int
            number of topics/components

        local_density_threshold : float
            Cutoff to detect outlier topics/components across repeats.
            See plot produced by this function.
        """
        call = self._create_call_str('consensus')

        local_density_threshold = f'{local_density_threshold:.2f}'
        consensus_cmd = f'{call} --local-density-threshold {local_density_threshold} --components {selected_K} --show-clustering'
        if self.verbose:
            print('Consensus command for K=%d:\n%s' % (selected_K, consensus_cmd))
        os.system(consensus_cmd)

        # load the image created
        dt = f'{local_density_threshold}'.replace('.', '_')
        img_name = f"{self.output_dir}/{self.run_name}/{self.run_name}.clustering.k_{selected_K}.dt_{dt}.png"
        from IPython.display import Image
        I = Image(filename=img_name, width=1000, height=1000)
        return I

    def load_results(self, selected_K:int, local_density_threshold:float):
        """
        Get the results of the cNMF, i.e the consensus components/topics and
        the consensus usages (topic composition of each datapoint)

        Parameters
        ----------

        selected_K : int
            number of topics/components

        local_density_threshold : float
            Cutoff to detect outlier topics/components across repeats.
            See plot produced by this function.
        """
        dt = f'{local_density_threshold:.2f}'.replace('.', '_')
        # f1 = f'{self.output_dir}/{self.run_name}/{self.run_name}.gene_spectra.k_{selected_K}.dt_{dt}.consensus.txt'
        # f2 = f'{self.output_dir}/{self.run_name}/{self.run_name}.gene_spectra_tpm.k_{selected_K}.dt_{dt}.consensus.txt'
        # f3 = f'{self.output_dir}/{self.run_name}/{self.run_name}.usages.k_{selected_K}.dt_{dt}.consensus.txt'


        fn_nmf_components = f'{self.output_dir}/{self.run_name}/{self.run_name}.spectra.k_{selected_K}.dt_{dt}.consensus.txt'
        fn_nmf_usage = f'{self.output_dir}/{self.run_name}/{self.run_name}.usages.k_{selected_K}.dt_{dt}.consensus.txt'

        nmf_components = pd.read_csv(fn_nmf_components, sep='\t', index_col=0)
        nmf_usage = pd.read_csv(fn_nmf_usage, sep='\t', index_col=0)

        return nmf_components, nmf_usage
