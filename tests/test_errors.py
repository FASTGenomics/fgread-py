import fgread
import pytest


def test_other(dset_other, list_datasets):
    with pytest.raises(
        NotImplementedError, match='Data sets with the "Other" format are unsupported'
    ):
        fgread.read_dataset(list_datasets[dset_other["id"]])


def test_other(dset_notset, list_datasets):
    with pytest.raises(KeyError, match="The format of the data set .* was not defined"):
        fgread.read_dataset(list_datasets[dset_notset["id"]])
