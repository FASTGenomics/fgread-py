import pytest
from pathlib import Path

import fgread

HERE = Path(__file__).parent
DATA_DIR = HERE / "data" / "readers"

DATASETS = [
    dict(id=1, dim=(298, 16892), title="Loom data set", format="Loom"),
    dict(id=3, dim=(10, 20), title="AnnData data set", format="AnnData"),
    dict(id=4, dim=(1222, 33538), title="10x (hdf5) data set", format="10x (hdf5)"),
    dict(
        id=5, dim=(20000, 99), title="Drop-Seq (tsv) data set", format="Drop-Seq (tsv)"
    ),
]

# special treatment for the unsupported seurat object
DATASET_SEURAT = dict(
    id=2, dim=(1222, 33538), title="Seurat Object data set", format="Seurat Object"
)


@pytest.fixture
def data_dir():
    return DATA_DIR


@pytest.fixture()
def list_datasets(data_dir):
    return fgread.list_datasets(data_dir)


# supported data sets
@pytest.fixture(params=DATASETS)
def dset(request):
    return request.param


# The seurat data set
@pytest.fixture
def dset_seurat():
    return DATASET_SEURAT
