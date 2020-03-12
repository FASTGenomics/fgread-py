import pytest
import fgread


def test_read_anndata(data_dir, list_datasets):
    adata = fgread.load_data("AnnData dataset", data_dir=data_dir)
    list_data = list_datasets[list_datasets["title"] == "AnnData dataset"]
    n_cells, n_genes = adata.X.shape
    print(adata.uns["ds_metadata"])
    assert n_genes == list_data["numberOfGenes"].values
    assert n_cells == list_data["numberOfCells"].values
    assert list(adata.uns["ds_metadata"]) == list_data["id"].values
    assert adata.obs["fg_id"][0] == list_data["id"].values


def test_seurat_raises(data_dir, list_datasets):
    with pytest.raises(NotImplementedError):
        fgread.load_data("Seurat Object dataset", data_dir=data_dir)
