import os
from pathlib import Path

import fgread
import pandas as pd
import pytest

HERE = Path(__file__).parent
DATA_DIR = HERE / "data" / "all_datasets"
DSETS = os.listdir(DATA_DIR)

DATASETS = fgread.ds_info(data_dir=DATA_DIR, output=True, pretty=False)
N_DATASETS = DATASETS.shape[0]


@pytest.fixture
def data_dir():
    return DATA_DIR


@pytest.fixture(params=list(range(N_DATASETS)))
def dset(request):
    metadata_keys_delete = [
        "state",
        "expressionDataFileInfos",
        "metaDataFileInfos",
    ]  # These keys do not go into anndata.uns

    dataset = DATASETS.loc[request.param]
    dataset.drop(metadata_keys_delete, inplace=True, errors="ignore")

    return dataset


# read json from all datasets
@pytest.fixture(params=DSETS)
def json_dset(request):
    series = pd.read_json(
        str(DATA_DIR) + "/" + request.param + "/dataset_info.json", typ="series"
    )
    df = series.to_frame().transpose()
    return df


@pytest.fixture
def list_datasets():
    return fgread.ds_info(data_dir=DATA_DIR, output=True, pretty=False)
