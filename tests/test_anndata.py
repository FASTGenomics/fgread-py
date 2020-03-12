import pytest
import fgread


def test_read_single_dataset(dset, list_datasets, data_dir):
    print(dset['title'])
    print(list_datasets)
    assert dset["title"] in list_datasets['title']
    adata = fgread.load_data("AnnData dataset", data_dir=data_dir)
    assert adata.shape == dset["dim"]
    assert adata.uns["metadata"]["title"] == dset["title"]
    assert adata.uns["metadata"]["format"] == dset["format"]
    assert (adata.obs["fg_id"] == dset["id"]).all()
    assert (adata.obs["fg_title"] == dset["title"]).all()


def test_seurat_raises(dset_seurat, list_datasets):
    with pytest.raises(NotImplementedError):
        fgread.load_data("Seurat Object dataset")
