import pytest
import pandas as pd
from pathlib import Path
import os

import fgread

HERE = Path(__file__).parent
DATA_DIR = HERE / "data" / "all_datasets"
DSETS = os.listdir(DATA_DIR)

DATASETS = [
    dict(id=1, dim=(298, 16892), title="Loom dataset", format="Loom"),
    dict(id=3, dim=(10, 20), title="AnnData dataset", format="AnnData"),
    dict(id=4, dim=(1222, 33538), title="10x (hdf5) dataset", format="10x (hdf5)"),
    dict(id=5, dim=(20000, 99), title="tab-separated text dataset", format="tab-separated text"),
    dict(id=8, dim=(30, 1000), title="mtx legacy dataset", format="10x (mtx)"),
    dict(id=9, dim=(30, 1000), title="mtx v3 dataset", format="10x (mtx)"),
    dict(id=10, dim=(499, 99), title="comma-separated text dataset", format="comma-separated text"),
    dict(id=11, dim=(499, 99), title="tab-separated text variant dataset", format="tab-separated text"),
]

# special treatment for the unsupported seurat object
DATASET_SEURAT = dict(
    id=2, dim=(1222, 33538), title="Seurat Object dataset", format="Seurat Object"
)

# special treatment for  "other dataset"
DATASET_OTHER = dict(id=6, title="Other dataset", format="Other")

# special treatment for "not set dataset"
DATASET_NOTSET = dict(id=7, title="Not set dataset", format="Not set")


@pytest.fixture
def data_dir():
    return DATA_DIR


# a list of all datasets
@pytest.fixture()
def list_datasets(data_dir):
    return fgread.ds_info(data_dir=data_dir)


# supported datasets
@pytest.fixture(params=DATASETS)
def dset(request):
    return request.param


# The seurat dataset
@pytest.fixture
def dset_seurat():
    return DATASET_SEURAT


# The notset dataset
@pytest.fixture
def dset_notset():
    return DATASET_NOTSET


# read json from all datasets
@pytest.fixture(params=DSETS)
def json_dset(request):
    series = pd.read_json(str(DATA_DIR) + "/" + request.param + "/dataset_info.json",
                          typ="series")
    df = series.to_frame().transpose()
    return df
