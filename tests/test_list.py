import fgread


def test_list(data_dir):
    dsets_list = fgread.list_datasets(data_dir)
    assert len(dsets_list) == 7
    assert set(dsets_list.keys()) == set(range(1, 8))
