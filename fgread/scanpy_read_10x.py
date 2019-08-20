"""Taken from scanpy/scanpy/readwrite.py as a temporary fix to the lack of read_10x_h5
in the anndata package.

"""

from anndata import AnnData
from typing import Union, Optional
from pathlib import Path
import numpy as np
import tables


def read_10x_h5(
    filename: Union[str, Path], genome: Optional[str] = None, gex_only: bool = True
) -> AnnData:
    """\
    Read 10x-Genomics-formatted hdf5 file.
    Parameters
    ----------
    filename
        Filename.
    genome
        Filter expression to this genes within this genome. For legacy 10x h5
        files, this must be provided if the data contains more than one genome.
    gex_only
        Only keep 'Gene Expression' data and ignore other feature types,
        e.g. 'Antibody Capture', 'CRISPR Guide Capture', or 'Custom'
    Returns
    -------
    Annotated data matrix, where obsevations/cells are named by their
    barcode and variables/genes by gene name. The data matrix is stored in
    `adata.X`, cell names in `adata.obs_names` and gene names in
    `adata.var_names`. The gene IDs are stored in `adata.var['gene_ids']`.
    The feature types are stored in `adata.var['feature_types']`
    """
    with tables.open_file(str(filename), "r") as f:
        v3 = "/matrix" in f
    if v3:
        adata = _read_v3_10x_h5(filename)
        if genome:
            if genome not in adata.var["genome"].values:
                raise ValueError(
                    f"Could not find data corresponding to genome '{genome}' in '{filename}'. "
                    f'Available genomes are: {list(adata.var["genome"].unique())}.'
                )
            adata = adata[:, list(map(lambda x: x == str(genome), adata.var["genome"]))]
        if gex_only:
            adata = adata[
                :,
                list(map(lambda x: x == "Gene Expression", adata.var["feature_types"])),
            ]
        return adata
    else:
        return _read_legacy_10x_h5(filename, genome=genome)


def _read_legacy_10x_h5(filename, *, genome=None, start=None):
    """
    Read hdf5 file from Cell Ranger v2 or earlier versions.
    """
    with tables.open_file(str(filename), "r") as f:
        try:
            children = [x._v_name for x in f.list_nodes(f.root)]
            if not genome:
                if len(children) > 1:
                    raise ValueError(
                        f"'{filename}' contains more than one genome. For legacy 10x h5 "
                        "files you must specify the genome if more than one is present. "
                        f"Available genomes are: {children}"
                    )
                genome = children[0]
            elif genome not in children:
                raise ValueError(
                    f"Could not find genome '{genome}' in '{filename}'. "
                    f"Available genomes are: {children}"
                )
            dsets = {}
            for node in f.walk_nodes("/" + genome, "Array"):
                dsets[node.name] = node.read()
            # AnnData works with csr matrices
            # 10x stores the transposed data, so we do the transposition right away
            from scipy.sparse import csr_matrix

            M, N = dsets["shape"]
            data = dsets["data"]
            if dsets["data"].dtype == np.dtype("int32"):
                data = dsets["data"].view("float32")
                data[:] = dsets["data"]
            matrix = csr_matrix((data, dsets["indices"], dsets["indptr"]), shape=(N, M))
            # the csc matrix is automatically the transposed csr matrix
            # as scanpy expects it, so, no need for a further transpostion
            adata = AnnData(
                matrix,
                dict(obs_names=dsets["barcodes"].astype(str)),
                dict(
                    var_names=dsets["gene_names"].astype(str),
                    gene_ids=dsets["genes"].astype(str),
                ),
            )
            return adata
        except KeyError:
            raise Exception("File is missing one or more required datasets.")


def _read_v3_10x_h5(filename, *, start=None):
    """
    Read hdf5 file from Cell Ranger v3 or later versions.
    """
    with tables.open_file(str(filename), "r") as f:
        try:
            dsets = {}
            for node in f.walk_nodes("/matrix", "Array"):
                dsets[node.name] = node.read()
            from scipy.sparse import csr_matrix

            M, N = dsets["shape"]
            data = dsets["data"]
            if dsets["data"].dtype == np.dtype("int32"):
                data = dsets["data"].view("float32")
                data[:] = dsets["data"]
            matrix = csr_matrix((data, dsets["indices"], dsets["indptr"]), shape=(N, M))
            adata = AnnData(
                matrix,
                dict(obs_names=dsets["barcodes"].astype(str)),
                dict(
                    var_names=dsets["name"].astype(str),
                    gene_ids=dsets["id"].astype(str),
                    feature_types=dsets["feature_type"].astype(str),
                    genome=dsets["genome"].astype(str),
                ),
            )
            return adata
        except KeyError:
            raise Exception("File is missing one or more required datasets.")
