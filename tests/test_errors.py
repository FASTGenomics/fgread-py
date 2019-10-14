import fgread
import pytest


def test_other(dset_other, list_datasets):
    with pytest.raises(
        NotImplementedError, match='Datasets with the "Other" format are unsupported'
    ):
        fgread.read_dataset(list_datasets[dset_other["id"]])


def test_notset(dset_notset, list_datasets):
    with pytest.raises(ValueError, match="The format of the dataset .* was not defined. If you can modify the dataset "\
                "please specify its format in its details page, otherwise ask the dataset owner to do that."):
        fgread.read_dataset(list_datasets[dset_notset["id"]])
