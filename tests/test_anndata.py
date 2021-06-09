import pytest
import fgread

load_fail = {
    "Seurat Object dataset": {"type": NotImplementedError},
    "No expression": {"type": TypeError},
    "Other dataset": {"type": KeyError},
    "mtx legacy dataset": {"type": KeyError},
    "mtx v3 dataset": {"type": KeyError},
    "Multifile dataset meta": {"type": TypeError},
}


def test_read_anndata(data_dir, dset):
    title = dset["title"]
    id = dset["id"]
    if title in load_fail:
        with pytest.raises(load_fail[title]["type"]):
            fgread.load_data(title, data_dir=data_dir)
    else:
        adata = fgread.load_data(title, data_dir=data_dir)
        n_cells, n_genes = adata.X.shape

        assert n_genes == dset["numberOfGenes"]
        assert n_cells == dset["numberOfCells"]
        assert adata.uns["ds_metadata"] == {id: {"title": title}}
        assert adata.uns["ds_metadata_raw"] == {id: str(dset.to_dict())}
        assert adata.obs["fg_id"][0] == id
