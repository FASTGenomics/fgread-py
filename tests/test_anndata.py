import pytest
import fgread


def test_read_single_dataset(dset, list_datasets):
    assert dset["id"] in list_datasets
    adata = fgread.read_dataset(list_datasets[dset["id"]])
    assert adata.shape == dset["dim"]
    assert adata.uns["metadata"]["title"] == dset["title"]
    assert adata.uns["metadata"]["format"] == dset["format"]
    assert (adata.obs["fg_id"] == dset["id"]).all()
    assert (adata.obs["fg_title"] == dset["title"]).all()


def test_seurat_raises(dset_seurat, list_datasets):
    with pytest.raises(NotImplementedError):
        fgread.read_dataset(list_datasets[dset_seurat["id"]])
