import fgread


def test_list(data_dir):
    dsets_list = fgread.list_datasets(data_dir)
    assert len(dsets_list) == 5
    assert set(dsets_list.keys()) == set([1, 2, 3, 4, 5])
