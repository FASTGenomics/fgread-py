import pytest
from pathlib import Path

import fgread

HERE = Path(__file__).parent
DATA_DIR = HERE / "data" / "readers"

DATASETS = [
    dict(id=1, dim=(298, 16892), title="Loom dataset", format="Loom"),
    dict(id=3, dim=(10, 20), title="AnnData dataset", format="AnnData"),
    dict(id=4, dim=(1222, 33538), title="10x (hdf5) dataset", format="10x (hdf5)"),
    dict(id=5, dim=(20000, 99), title="Drop-Seq (tsv) dataset", format="Drop-Seq (tsv)"),
    dict(id=8, dim=(2700, 32738), title="3k PBMCs dataset", format="10x (mtx)"),
]

# special treatment for the unsupported seurat object
DATASET_SEURAT = dict(
    id=2, dim=(1222, 33538), title="Seurat Object dataset", format="Seurat Object"
)

# special treatment for the unsupported seurat object
DATASET_OTHER = dict(id=6, title="Other dataset", format="Other")

# special treatment for the unsupported seurat object
DATASET_NOTSET = dict(id=7, title="Not set dataset", format="Not set")


@pytest.fixture
def data_dir():
    return DATA_DIR


# a list of all datasets
@pytest.fixture()
def list_datasets(data_dir):
    return fgread.get_datasets(data_dir)


# supported datasets
@pytest.fixture(params=DATASETS)
def dset(request):
    return request.param


# The seurat dataset
@pytest.fixture
def dset_seurat():
    return DATASET_SEURAT


# The other dataset
@pytest.fixture
def dset_other():
    return DATASET_OTHER


# The notset dataset
@pytest.fixture
def dset_notset():
    return DATASET_NOTSET
