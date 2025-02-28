import pytest
import numpy as np
import pandas as pd
import scanpy as sc
import os
import scipy.sparse as sp
from cnmf import cNMF, save_df_to_npz, load_df_from_npz

# Global parameters for data simulation
NUM_CELLS = 100
NUM_GENES = 500
BINOM_N = 100
BINOM_P = 0.01
SEED = 42

@pytest.fixture
def mock_cnmf(tmp_path):
    return cNMF(output_dir=str(tmp_path), name="test")

def generate_counts_file(tmp_path, file_format, dtype=np.int64, zero_count=False):
    """
    Generates a synthetic single-cell RNA-seq counts file in various formats.
    
    Args:
        tmp_path (Path): Temporary path for storing the file.
        file_format (str): One of ['txt', 'npz', 'h5ad'].
        dtype (numpy dtype, optional): The data type to store (int64, float32, etc.).
        zero_count (bool, optional): If True, makes the first cell have zero counts.

    Returns:
        str: Path to the generated counts file.
    """
    np.random.seed(SEED)
    data = np.random.binomial(n=BINOM_N, p=BINOM_P, size=(NUM_CELLS, NUM_GENES)).astype(dtype)
    
    if zero_count:
        data[0, :] = 0  # Introduce zero-count cells

    if file_format == "txt":
        df = pd.DataFrame(data, columns=[f"gene{i}" for i in range(NUM_GENES)],
                          index=[f"cell{i}" for i in range(NUM_CELLS)])
        counts_fn = tmp_path / f"counts_{dtype.__name__}.txt"
        df.to_csv(counts_fn, sep='\t')

    elif file_format == "npz":
        df = pd.DataFrame(data, columns=[f"gene{i}" for i in range(NUM_GENES)],
                          index=[f"cell{i}" for i in range(NUM_CELLS)])
        counts_fn = tmp_path / f"counts_{dtype.__name__}.npz"
        save_df_to_npz(df, counts_fn)

    elif file_format == "h5ad":
        adata = sc.AnnData(X=sp.csr_matrix(data))
        counts_fn = tmp_path / f"counts_{dtype.__name__}.h5ad"
        adata.write_h5ad(counts_fn)

    else:
        raise ValueError("Unsupported file format. Choose from ['txt', 'npz', 'h5ad'].")

    return str(counts_fn)

@pytest.mark.parametrize("file_format", ["txt", "npz", "h5ad"])
@pytest.mark.parametrize("dtype", [np.int64, np.float32, np.float64])
@pytest.mark.parametrize("densify", [True, False])
def test_prepare(mock_cnmf, file_format, dtype, densify, tmp_path):
    counts_fn = generate_counts_file(tmp_path, file_format, dtype)
    
    output_dir = tmp_path / "output"
    os.makedirs(output_dir, exist_ok=True)
    
    mock_cnmf.prepare(counts_fn, components=[5, 10], n_iter=10, densify=densify)
    
    # Check if output files were created
    expected_files = [
        mock_cnmf.paths['normalized_counts'],
        mock_cnmf.paths['nmf_replicate_parameters'],
        mock_cnmf.paths['nmf_run_parameters'],
        mock_cnmf.paths['nmf_genes_list'],
        mock_cnmf.paths['tpm'],
        mock_cnmf.paths['tpm_stats']
    ]
    
    for file in expected_files:
        assert os.path.exists(file), f"Expected output file {file} not found."
    
    # Clean up after test
    for file in expected_files:
        os.remove(file)

@pytest.mark.parametrize("file_format", ["txt", "npz", "h5ad"])
@pytest.mark.parametrize("dtype", [np.int64, np.float32, np.float64])
@pytest.mark.parametrize("densify", [True, False])
def test_prepare_raises_on_zero_count_cells(mock_cnmf, file_format, dtype, densify, tmp_path):
    counts_fn = generate_counts_file(tmp_path, file_format, dtype, zero_count=True)

    with pytest.raises(Exception, match="Error: .* cells have zero counts of overdispersed genes.*"):
        mock_cnmf.prepare(counts_fn, components=[5, 10], n_iter=10, densify=densify)
