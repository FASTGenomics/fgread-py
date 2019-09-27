import anndata
import re
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
from .dataset import DataSet


def read_loom_to_anndata(dataset: DataSet):
    """Reads a dataset in the loom format into the AnnData format."""

    adata = anndata.read_loom(dataset.file)
    return adata


def read_seurat_to_anndata(dataset: DataSet):
    """Reads a dataset in the Seurat format into the AnnData format (not implemented)."""

    raise NotImplementedError("Reading of Seurat files not implemented.")


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
    return read_densemat_to_anndata(dataset, sep="\t")


def read_densecsv_to_anndata(dataset: DataSet):
    return read_densemat_to_anndata(dataset, sep=",")


def read_densemat_to_anndata(dataset: DataSet, sep=None):
    """Reads a dataset in the DropSeq format into the AnnData format."""

    file = dataset.file

    with open(file) as f:
        cells = f.readline().replace('"', '').split(sep)
        nextline = f.readline().replace('"', '').split(sep)
        n_cells = len(nextline)-1
        cells = cells[-n_cells:]
        samples = [re.search("(.*)_", c).group(1) for c in cells]

    genes = pd.read_csv(
        file, skiprows=1, usecols=(0,), header=None, names=["GeneID"]
    ).set_index("GeneID")
    X = np.loadtxt(
        file,
        delimiter=sep,
        skiprows=1,
        usecols=range(1, len(cells) + 1),
        dtype=np.float32,
    ).T
    X = sp.csr_matrix(X)

    var = genes
    obs = pd.DataFrame(
        samples, columns=["sample"], index=pd.Series(cells, name="CellID")
    )

    adata = anndata.AnnData(X=X, var=var, obs=obs)
    return adata
