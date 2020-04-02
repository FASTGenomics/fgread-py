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
    if title in load_fail:
        with pytest.raises(load_fail[title]["type"]):
            fgread.load_data(title, data_dir=data_dir)
    else:
        pass
