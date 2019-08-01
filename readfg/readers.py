import scanpy as sc
import json
import re
import numpy as np
import pandas as pd
import scipy.sparse as sp


def read_manifest(dataset_dir):
    with open(dataset_dir / "manifest.json") as f:
        return json.load(f)


def read_loom(dataset_dir):
    adata = sc.read_loom(dataset_dir / "data.loom")
    adata.uns["manifest"] = read_manifest(dataset_dir)
    return adata


def read_seurat(dataset_dir):
    raise NotImplementedError("Reading of Seurat files not implemented.")


def read_anndata(dataset_dir):
    adata = sc.read_h5ad(dataset_dir / "data.h5ad")
    adata.uns["manifest"] = read_manifest(dataset_dir)
    return adata


def read_10x_hdf5(dataset_dir):
    adata = sc.read_10x_h5(dataset_dir / "data.h5")
    adata.uns["manifest"] = read_manifest(dataset_dir)
    return adata


def read_dropseq(dataset_dir):

    file = dataset_dir / "data.tsv"

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

    adata = sc.AnnData(X=X, var=var, obs=obs)
    adata.uns["manifest"] = read_manifest(dataset_dir)
    return adata
