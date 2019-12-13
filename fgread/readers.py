import anndata
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
from .dataset import DataSet
from . import BLOGURL
from tqdm.auto import tqdm


def read_loom_to_anndata(dataset: DataSet):
    """Reads a dataset in the loom format into the AnnData format."""

    adata = anndata.read_loom(dataset.file)
    return adata


def read_seurat_to_anndata(dataset: DataSet):
    """Reads a dataset in the Seurat format into the AnnData format (not implemented)."""

    raise NotImplementedError(f"Reading of Seurat files not implemented.\nSee {BLOGURL} for more information.")


def read_anndata_to_anndata(dataset: DataSet):
    """Reads a dataset in the AnnData format into the AnnData format."""

    adata = anndata.read_h5ad(dataset.file)
    return adata


def read_10xhdf5_to_anndata(dataset: DataSet):
    """Reads a dataset in the 10x hdf5 format into the AnnData format."""

    adata = sc.read_10x_h5(dataset.file)
    return adata


def read_10xmtx_to_anndata(dataset: DataSet):
    """Reads a dataset in the 10x mtx format into the AnnData format."""

    adata = sc.read_10x_mtx(dataset.path)
    return adata


def read_densetsv_to_anndata(dataset: DataSet):
    """Reads a dense text file in tsv format into the AnnData format."""

    return read_densemat_to_anndata(dataset, sep="\t")


def read_densecsv_to_anndata(dataset: DataSet):
    """Reads a dense text file in csv format into the AnnData format."""

    return read_densemat_to_anndata(dataset, sep=",")


def read_densemat_to_anndata(dataset: DataSet, sep=None):
    """Helper function to read dense text files in tsv and csv format.
    The separator (tab or comma) is passed by the corresponding function."""

    file = dataset.file

    genes = pd.read_csv(file, skiprows=1, usecols=(0,), header=None, names=["GeneID"], delimiter=sep,
                        squeeze=True)
    with open(file) as f:
        # Get cell names from first row
        cells = f.readline().strip().replace('"', '').split(sep)
        # Read second row and get the real cell count, as the first row can have different formats
        line = f.readline().strip().replace('"', '')
        n_cells = len(line.split(sep)) - 1
        cells = pd.Series(cells[-n_cells:], name="CellID")

        # Initialize matrix
        lil_mat = sp.lil_matrix((len(cells), len(genes)), dtype=np.float64)

        prog = tqdm(total=len(genes) * len(cells), unit="counts", desc="Reading counts", unit_scale=True)
        gene_idx = 0
        while line:
            prog.update(len(cells))
            expr_lst = line.split(sep)[1:]

            for idx in range(len(expr_lst)):
                if expr_lst[idx] == 0:
                    lil_mat[idx, gene_idx] = expr_lst[idx]

            line = f.readline().strip().replace('"', '')
            gene_idx += 1

        prog.close()

    # Convert to anndata
    obs = pd.DataFrame(index=cells)
    var = pd.DataFrame(index=genes)

    adata = anndata.AnnData(X=lil_mat, obs=obs, var=var, dtype=np.float64)

    return adata
