import re
from pathlib import Path
import json
from . import readers

DEFAULT_READERS = {
    "Loom": readers.read_loom_to_anndata,
    "Seurat Object": readers.read_seurat_to_anndata,
    "AnnData": readers.read_anndata_to_anndata,
    "10x (hdf5)": readers.read_10xhdf5_to_anndata,
    "Drop-Seq (tsv)": readers.read_dropseq_to_anndata,
}

DATA_DIR = "/fastgenomics/data"
DATASET_INFO_FILE = "dataset_info.json"


class DataSet(object):
    """Represents a data set on FASTGenomics, including the relative location and the
contents of the metadata.json file.

    """

    def __init__(self, path):
        self.path = path

        if not self.path.exists():
            raise FileNotFoundError(filename=self.path)

        self.metadata = self.read_metadata()
        self.format = self.metadata["format"]
        self.title = self.metadata["title"]
        self.file = self.path / self.metadata["file"]
        self.id = int(self.path.name.split("_")[-1])

    def read_metadata(self):
        with open(self.path / DATASET_INFO_FILE) as f:
            return json.load(f)

    def __repr__(self):
        return f"DataSet: {self.title} [{self.format}]"


def read_dataset(dataset: DataSet, additional_readers={}):
    """Reads a single data set.  Dispatches to specific readers based on the contents of the
`dataset.format`.

    """

    format = dataset.format
    title = dataset.title
    path = dataset.path

    readers = {**DEFAULT_READERS, **additional_readers}

    if format in readers:
        print(
            f'Loading data set "{title}" in format "{format}" from directory "{path}".'
        )
        adata = readers[format](dataset)
        adata.uns["metadata"] = dataset.metadata
        adata.var["fg_title"] = dataset.title
        adata.var["fg_id"] = dataset.id
        return adata
    else:
        raise KeyError(f'Unsupported format "{format}", use one of {readers}')


def list_datasets(data_dir=DATA_DIR):
    """Lists available data sets."""

    data_dir = Path(data_dir)
    paths = [
        f
        for f in data_dir.iterdir()
        if f.is_dir() and re.match(r"^dataset_\d{4}$", f.name)
    ]
    return {dataset.id: dataset for dataset in map(DataSet, paths)}


def read_datasets(datasets=None, additional_readers={}, data_dir=DATA_DIR):
    """Reads all data sets."""

    datasets = datasets or list_datasets(data_dir)
    return {
        id: read_dataset(dataset, additional_readers=additional_readers)
        for id, dataset in datasets.items()
    }
