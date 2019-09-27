import fgread


def test_list(data_dir):
    dsets_list = fgread.get_datasets(data_dir)
    assert len(dsets_list) == 11
    assert set(dsets_list.keys()) == set(range(1, 12))


def test_representation(list_datasets):
    for dset in list_datasets.items():
        rep = dset.__repr__()
        assert "id:" in rep
        assert "title:" in rep
        assert "format:" in rep
        assert "path:" in rep
