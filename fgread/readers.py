import anndata
import re
import numpy as np
import pandas as pd
import scipy.sparse as sp
from .scanpy_read_10x import read_10x_h5


def read_loom_to_anndata(dataset):
    """Reads a data set in the loom format."""

    adata = anndata.read_loom(dataset.file)
    return adata


def read_seurat_to_anndata(dataset):
    """Reads a data set in the Seurat format (not implemented)."""

    raise NotImplementedError("Reading of Seurat files not implemented.")


def read_anndata_to_anndata(dataset):
    """Reads a data set in the AnnData format."""

    adata = anndata.read_h5ad(dataset.file)
    return adata


def read_10xhdf5_to_anndata(dataset):
    """Reads a data set in the 10x hdf5 format."""

    # todo replace with anndata.read_10x_h5 once read_10x_h5 is moved to anndata (if
    # ever)
    adata = read_10x_h5(dataset.file)
    return adata


def read_dropseq_to_anndata(dataset):
    """Reads a data set in the DropSeq format."""

    file = dataset.file

    with open(file) as f:
        cells = f.readline().replace('"', "").split("\t")
        samples = [re.search("(.*)_", c).group(1) for c in cells]

    genes = pd.read_csv(
        file, sep="\t", skiprows=1, usecols=(0,), header=None, names=["GeneID"]
    ).set_index("GeneID")
    X = np.loadtxt(
        file,
        delimiter="\t",
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
