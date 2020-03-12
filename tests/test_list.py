import fgread


def test_list(data_dir):
    dsets = fgread.ds_info()
    assert dsets.shape()[0] == 11


def test_representation(list_datasets):
    for dset in list_datasets.items():
        rep = dset.__repr__()
        assert "id:" in rep
        assert "title:" in rep
        assert "format:" in rep
        assert "path:" in rep
