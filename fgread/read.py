import re
from pathlib import Path
from . import readers
from .dataset import DataSet

DEFAULT_READERS = {
    "Loom": readers.read_loom_to_anndata,
    "Seurat Object": readers.read_seurat_to_anndata,
    "AnnData": readers.read_anndata_to_anndata,
    "10x (hdf5)": readers.read_10xhdf5_to_anndata,
    "Drop-Seq (tsv)": readers.read_dropseqtsv_to_anndata,
}

DATA_DIR = "/fastgenomics/data"


def read_dataset(dataset: DataSet, additional_readers={}):
    """Reads a single data set.  Dispatches to specific readers based on the value of
    the ``dataset.format``.

    :param dataset: Object of class :py:class:`~.dataset.DataSet` to be read.
    :param additional_readers: Used to specify your own readers for the specific data
        set format.  Highly experimental and not tested.

    :returns: AnnData object containing the loaded data set.
    """

    format = dataset.format
    title = dataset.title
    path = dataset.path

    readers = {**DEFAULT_READERS, **additional_readers}

    if format == "Other":
        raise NotImplementedError(
            f'The format of the data set "{title}" is "{format}".  Data sets with the "{format}" format are unsupported by this module and have to be loaded manually.'
        )
    elif format == "Not set":
        raise KeyError(
            f'The format of the data set "{title}" was not defined.  If you can modify the data set please specify its format in its Details page, otherwise ask the data set owner to do that.'
        )
    elif format in readers:
        print(
            f'Loading data set "{title}" in format "{format}" from directory "{path}".'
        )
        adata = readers[format](dataset)
        adata.uns["metadata"] = dataset.metadata
        adata.obs["fg_title"] = dataset.title
        adata.obs["fg_id"] = dataset.id
        return adata
    else:
        raise KeyError(f'Unsupported format "{format}", use one of {readers}')


def list_datasets(data_dir=DATA_DIR):
    """Lists all available data sets.  This is a convenience function used to gather all
    information specified in the FASTGenomics environment.  The returned value can be
    either used to manually load data sets or passed to the :py:func:`read_dataset` or
    :py:func:`read_datasets` functions.

    :param data_dir: Specify the main data directory.  Useful for testing the module,
        defaults to the FASTGenomics path ``/fastgenomics/data``.

    :returns: A dictionary where keys are data set ids (the ``xxxx`` part of
              ``/fastgenomics/data/dataset_xxxx``) and values are the corresponding
              :py:class:`~dataset.DataSet` objects.
    """

    data_dir = Path(data_dir)
    paths = [
        f
        for f in data_dir.iterdir()
        if f.is_dir() and re.match(r"^dataset_\d{4}$", f.name)
    ]
    return {dataset.id: dataset for dataset in map(DataSet, paths)}


def read_datasets(datasets=None, additional_readers={}, data_dir=DATA_DIR):
    """Reads all data sets and returns them as AnnData objects.  Internally uses
    :py:func:`read_dataset` to read the datasets.

    :param datasets: If specified, read the datasets from this dictionary.  Can be
        useful for e.g. filtering some data set types.
    :param additional_readers: Used to specify your own readers for the specific data
        set format.  Highly experimental and not tested.
    :param data_dir: Specify the main data directory.  Only used when
        ``datasets==None``.  Useful for testing the module, defaults to the FASTGenomics
        path ``/fastgenomics/data``.

    :returns: A dictionary of data set objects, where the keys are data set ids and the
              values are the corresponding AnnData objects.
    """

    datasets = datasets or list_datasets(data_dir)
    return {
        dataset.id: read_dataset(dataset, additional_readers=additional_readers)
        for dataset in datasets.values()
    }
