import anndata
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
from pathlib import Path
from .dataset import DataSet
from . import DOCSURL


def read_loom_to_anndata(ds_file: Path):
    """Reads a dataset in the loom format into the AnnData format."""

    adata = anndata.read_loom(ds_file)
    return adata


def read_seurat_to_anndata(ds_file: Path):
    """Reads a dataset in the Seurat format into the AnnData format (not implemented)."""

    raise NotImplementedError(
        f"Reading of Seurat files not implemented.\nSee {DOCSURL} for more information."
    )


def read_anndata_to_anndata(ds_file: Path):
    """Reads a dataset in the AnnData format into the AnnData format."""

    adata = anndata.read_h5ad(ds_file)
    return adata


def read_10xhdf5_to_anndata(ds_file: Path):
    """Reads a dataset in the 10x hdf5 format into the AnnData format."""

    adata = sc.read_10x_h5(ds_file)
    return adata


def read_10xmtx_to_anndata(ds_file: Path):
    """Reads a dataset in the 10x mtx format into the AnnData format."""

    adata = sc.read_10x_mtx(ds_file.parent)
    return adata


def read_densetsv_to_anndata(ds_file: Path):
    """Reads a dense text file in tsv format into the AnnData format."""

    return read_densemat_to_anndata(ds_file, sep="\t")


def read_densecsv_to_anndata(ds_file: Path):
    """Reads a dense text file in csv format into the AnnData format."""

    return read_densemat_to_anndata(ds_file, sep=",")


def read_densemat_to_anndata(ds_file: Path, sep=None):
    """Helper function to read dense text files in tsv and csv format.
    The separator (tab or comma) is passed by the corresponding function."""

    file = ds_file

    with open(file) as f:
        cells = f.readline().replace('"', "").split(sep)
        nextline = f.readline().replace('"', "").split(sep)
        n_cells = len(nextline) - 1
        cells = cells[-n_cells:]

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
    obs = pd.DataFrame(cells, columns=["sample"], index=pd.Series(cells, name="CellID"))

    adata = anndata.AnnData(X=X, var=var, obs=obs)
    return adata
