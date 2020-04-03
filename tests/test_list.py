import os


# test number of datasets
def test_list(data_dir, list_datasets):
    nfolders = len(os.listdir(data_dir))
    assert list_datasets.shape[0] == nfolders


# test equality of metadata in json and list
def test_representation(json_dset, list_datasets):
    json_data = json_dset.drop("schemaVersion", axis=1)
    title = json_data["title"].values[0]
    list_data = list_datasets.loc[list_datasets["title"] == title]
    for col in json_data.columns:
        assert list_data[col].values == json_data[col].values
