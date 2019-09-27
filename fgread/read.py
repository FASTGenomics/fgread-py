import re
from pathlib import Path
from . import readers
from .dataset import DataSet
from collections import OrderedDict

DEFAULT_READERS = {
    "Loom": readers.read_loom_to_anndata,
    "Seurat Object": readers.read_seurat_to_anndata,
    "AnnData": readers.read_anndata_to_anndata,
    "10x (hdf5)": readers.read_10xhdf5_to_anndata,
    "10x (mtx)": readers.read_10xmtx_to_anndata,
    "tab-separated text": readers.read_densetsv_to_anndata,
    "comma-separated text": readers.read_densecsv_to_anndata,
}

DATA_DIR = "/fastgenomics/data"


def read_dataset(dataset: DataSet, additional_readers={}):
    """Reads a single dataset.  Dispatches to specific readers based on the value of
    the ``dataset.format``.

    :param dataset: Object of class :py:class:`~.dataset.DataSet` to be read.
    :param additional_readers: Used to specify your own readers for the specific data
        set format.  Highly experimental and not tested.

    :returns: AnnData object containing the loaded dataset.
    """

    format = dataset.format
    title = dataset.title
    path = dataset.path

    readers = {**DEFAULT_READERS, **additional_readers}

    if format == "Other":
        raise NotImplementedError(
            f'The format of the dataset "{title}" is "{format}".  Datasets with the "{format}" format are unsupported "\
            "by this module and have to be loaded manually.'
        )
    elif format == "Not set":
        raise KeyError(
            f'The format of the dataset "{title}" was not defined.  If you can modify the dataset please specify its "\
            "format in its Details page, otherwise ask the dataset owner to do that.'
        )
    elif format in readers:
        print(
            f'Loading dataset "{title}" in format "{format}" from directory "{path}"...'
        )
        adata = readers[format](dataset)
        adata.uns["metadata"] = dataset.metadata
        adata.obs["fg_title"] = dataset.title
        adata.obs["fg_id"] = dataset.id
        n_genes = adata.shape[0]
        n_cells = adata.shape[1]
        print(
            f'Loaded dataset "{title}" with {n_genes} genes and {n_cells} cells\n'
        )
        return adata
    else:
        raise KeyError(f'Unsupported format "{format}", use one of {readers}')


def get_datasets(data_dir=DATA_DIR):
    """Lists all available datasets.  This is a convenience function used to gather all
    information specified in the FASTGenomics environment.  The returned value can be
    either used to manually load datasets or passed to the :py:func:`read_dataset` or
    :py:func:`read_datasets` functions.

    :param data_dir: Specify the main data directory.  Useful for testing the module,
        defaults to the FASTGenomics path ``/fastgenomics/data``.

    :returns: A dictionary where keys are dataset ids (the ``xxxx`` part of
              ``/fastgenomics/data/dataset_xxxx``) and values are the corresponding
              :py:class:`~dataset.DataSet` objects.
    """

    data_dir = Path(data_dir)
    paths = [
        subdir
        for subdir in sorted(data_dir.iterdir())
        if subdir.is_dir() and re.match(r"^dataset_\d{4}$", subdir.name)
    ]
    return OrderedDict({dataset.id: dataset for dataset in map(DataSet, paths)})


def print_datasets(data_dir=DATA_DIR):
    """prints the list of available datasets

    :param data_dir: Specify the main data directory.  Useful for testing the module,
        defaults to the FASTGenomics path ``/fastgenomics/data``.
    """

    datasets = get_datasets(data_dir=data_dir)

    for index, ds in datasets.items():
        print(f"Dataset: {index}:", ds)
        print()


def read_datasets(datasets=None, additional_readers={}, data_dir=DATA_DIR):
    """Reads all datasets and returns them as AnnData objects.  Internally uses
    :py:func:`read_dataset` to read the datasets.

    :param datasets: If specified, read the datasets from this dictionary.  Can be
        useful for e.g. filtering some dataset types.
    :param additional_readers: Used to specify your own readers for the specific data
        set format.  Highly experimental and not tested.
    :param data_dir: Specify the main data directory.  Only used when
        ``datasets==None``.  Useful for testing the module, defaults to the FASTGenomics
        path ``/fastgenomics/data``.

    :returns: A dictionary of dataset objects, where the keys are dataset ids and the
              values are the corresponding AnnData objects.
    """

    datasets = datasets or get_datasets(data_dir)
    return {
        dataset.id: read_dataset(dataset, additional_readers=additional_readers)
        for dataset in datasets.values()
    }
