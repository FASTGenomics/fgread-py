from . import readers, readers_old, BLOGURL, DS_URL_PREFIX
from .dataset import DataSet, DatasetDict
import re
from pathlib import Path
import pandas as pd
import json
from typing import Optional, Union
import logging


# configure logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
logger.addHandler(ch)


DEFAULT_READERS_OLD = {
    "Loom": readers_old.read_loom_to_anndata,
    "Seurat Object": readers_old.read_seurat_to_anndata,
    "AnnData": readers_old.read_anndata_to_anndata,
    "10x (hdf5)": readers_old.read_10xhdf5_to_anndata,
    "10x (mtx)": readers_old.read_10xmtx_to_anndata,
    "tab-separated text": readers_old.read_densetsv_to_anndata,
    "comma-separated text": readers_old.read_densecsv_to_anndata,
}

DEFAULT_READERS = {
    "Loom": readers.read_loom_to_anndata,
    "Seurat Object": readers.read_seurat_to_anndata,
    "AnnData": readers.read_anndata_to_anndata,
    "10x (hdf5)": readers.read_10xhdf5_to_anndata,
    "10x (mtx)": readers.read_10xmtx_to_anndata,
    "tab-separated text": readers.read_densetsv_to_anndata,
    "comma-separated text": readers.read_densecsv_to_anndata,
}

DATA_DIR = Path("/fastgenomics/data")


def get_ds_paths(data_dir: Union[str, Path] = DATA_DIR) -> list:
    """Gets available datasets for this analysis from path.
    
    Parameters
    ----------
    data_dir : Union[str,Path], optional
        Directory containing the datasets, e.g. "fastgenomics/data", by default DATA_DIR
    
    Returns
    -------
    list
        A list of dataset paths
    """
    data_dir = Path(data_dir)
    paths = [
        subdir
        for subdir in sorted(data_dir.iterdir())
        if subdir.is_dir() and re.match(r"^dataset_\d{4}$", subdir.name)
    ]
    assert paths != [], "There is no data available in this analysis."
    return paths


def ds_info(
    data_dir: Path = DATA_DIR,
    pretty: bool = True,
    id: Optional[str] = None,
    output: bool = True,
) -> pd.DataFrame:
    """[summary]
    
    Parameters
    ----------
    data_dir : Path, optional
        Directory containing the datasets, e.g. "fastgenomics/data", by default DATA_DIR
    pretty : bool, optional
        Whether to display some nice output, by default True
    id : Optional[str], optional
        A single dataset ID. If set, only this dataset will be displayed. Recommended to use with `pretty`, by default None
    output : bool, optional
        Whether to return a DataFrame or not, by default True
    
    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing all, or a single dataset (depends on `id`)
    """

    if not pretty and not output:
        logger.warning(
            '"You have set "pretty" and "output" to false. Hence, this function will do/return nothing.'
        )
    ds_paths = get_ds_paths(data_dir=data_dir)
    ds_df = pd.DataFrame()
    for ds_path in ds_paths:
        with open(ds_path / "dataset_info.json") as f:
            ds_info = json.load(f)
            ds_info["path"] = ds_path
            _ = ds_info.pop("schemaVersion", None)
        ds_df = ds_df.append(ds_info, ignore_index=True)

    # sort colnames
    sort_order = [
        "title",
        "id",
        "format",
        "organism",
        "tissue",
        "numberOfCells",
        "numberOfGenes",
    ]
    col_names = ds_df.columns.values.tolist()
    col_names_sorted = [name for name in sort_order if name in col_names]
    [col_names.remove(name) for name in sort_order]
    col_names_sorted.extend(col_names)
    ds_df = ds_df[col_names_sorted]

    def add_url(title, id):
        return f'<a href="{DS_URL_PREFIX}{id}" target="_blank">{title}</a>'

    ds_df["title"] = ds_df.apply(lambda x: add_url(x.title, x.id), axis=1)

    def disp_pretty_df(df, index=True, header=True):
        try:
            from IPython.display import display, Markdown

            df_html = df.to_html(
                render_links=True, escape=False, header=header, index=index
            )
            display(Markdown(df_html))
        except:
            logger.warning(
                "IPython not available. Pretty printing only works in Jupyter Notebooks."
            )

    if id:
        single_dict = ds_df.loc[ds_df["id"] == id].squeeze().to_dict()
        single_df = pd.DataFrame(columns=["key", "value"])
        for key, value in single_dict.items():
            section = pd.DataFrame([{"key": f"<b>{key}</b>", "value": value}])
            single_df = single_df.append(section, ignore_index=True)

        if pretty:
            disp_pretty_df(single_df, header=False, index=False)

        if output:
            return ds_df.loc[ds_df["id"] == id]

    else:
        if pretty:
            disp_pretty_df(ds_df)

        if output:
            return ds_df


def load_data(
    id: Optional[str] = None, data_dir: Path = DATA_DIR, additional_readers: dict = {}
):
    """Docstring goes here
    """
    ds_df = ds_info(data_dir=data_dir)
    readers = {**DEFAULT_READERS, **additional_readers}

    if not id:
        assert (
            len(ds_df) == 1
        ), "There is more than one dataset available. Please select."

        id = ds_df.iloc[0]["id"]

    format = ds_df.loc[ds_df["id"] == id, ["format"]]
    title = ds_df.loc[ds_df["id"] == id, ["title"]]
    path = ds_df.loc[ds_df["id"] == id, ["path"]]
    file = ds_df.loc[ds_df["id"] == id, ["file"]]

    if format in readers:
        logger.info(
            f'Loading dataset "{title}" in format "{format}" from directory "{path}"...\n'
        )
        adata = readers[format](Path(path) / file)
        adata.uns["metadata"] = {id: ds_df.loc[ds_df["id"] == id].to_dict()}
        adata.obs["fg_id"] = id
        n_genes = adata.shape[1]
        n_cells = adata.shape[0]
        logger.info(
            f'Loaded dataset "{title}" with {n_cells} cells and {n_genes} genes.\n'
            f"==================================================================\n"
        )
        return adata

    elif format == "Other":
        raise NotImplementedError(
            f'The format of the dataset "{title}" is "{format}".  Datasets with the "{format}" format are '
            f"unsupported by this module and have to be loaded manually.\nSee {BLOGURL} for more information."
        )
    elif format == "Not set":
        raise ValueError(
            f'The format of the dataset "{title}" was not defined. If you can modify the dataset please specify '
            f"its format in its details page, otherwise ask the dataset owner to do that.\nSee {BLOGURL} for more information."
        )
    else:
        raise KeyError(
            f'Unsupported format "{format}", use one of {list(readers)} or implement your '
            f"own reading function.\nSee {BLOGURL} for more information."
        )


def get_datasets(data_dir=DATA_DIR):
    """Gets all available datasets.  This is a convenience function used to gather all
    information specified in the FASTGenomics environment. The returned value can be
    either used to manually load datasets or passed to the :py:func:`read_dataset` or
    :py:func:`read_datasets` functions.

    :param data_dir: Specify the main data directory.  Useful for testing the module,
        defaults to the FASTGenomics path ``/fastgenomics/data``.

    :returns: A dictionary where keys are dataset ids (the ``xxxx`` part of
              ``/fastgenomics/data/dataset_xxxx``) and values are the corresponding
              :py:class:`~dataset.DataSet` objects.
    """

    data_dir = Path(data_dir)
    paths = get_ds_paths(data_dir=data_dir)
    datasets = DatasetDict({dataset.id: dataset for dataset in map(DataSet, paths)})

    return datasets


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

    readers = {**DEFAULT_READERS_OLD, **additional_readers}

    if format in readers:
        logger.info(
            f'Loading dataset "{title}" in format "{format}" from directory "{path}"...\n'
        )
        adata = readers[format](dataset)
        adata.uns["metadata"] = dataset.metadata
        adata.obs["fg_title"] = dataset.title
        adata.obs["fg_id"] = dataset.id
        n_genes = adata.shape[1]
        n_cells = adata.shape[0]
        logger.info(
            f'Loaded dataset "{title}" with {n_cells} cells and {n_genes} genes.\n'
            f"==================================================================\n"
        )
        return adata

    elif format == "Other":
        raise NotImplementedError(
            f'The format of the dataset "{title}" is "{format}".  Datasets with the "{format}" format are '
            f"unsupported by this module and have to be loaded manually.\nSee {BLOGURL} for more information."
        )
    elif format == "Not set":
        raise ValueError(
            f'The format of the dataset "{title}" was not defined. If you can modify the dataset please specify '
            f"its format in its details page, otherwise ask the dataset owner to do that.\nSee {BLOGURL} for more information."
        )
    else:
        raise KeyError(
            f'Unsupported format "{format}", use one of {list(readers)} or implement your '
            f"own reading function.\nSee {BLOGURL} for more information."
        )


def read_datasets(datasets=None, additional_readers={}, data_dir=DATA_DIR):
    """Reads all specified datasets and returns them as AnnData objects.  Internally uses
    :py:func:`read_dataset` to read the datasets.

    :param datasets: If specified, read the datasets from this dictionary or list.  Can be
        useful for e.g. filtering some dataset types.
    :param additional_readers: Used to specify your own readers for the specific data
        set format.  Highly experimental and not tested.
    :param data_dir: Specify the main data directory.  Only used when
        ``datasets==None``.  Useful for testing the module, defaults to the FASTGenomics
        path ``/fastgenomics/data``.

    :returns: If multiple datasets are given returns a dictionary of dataset objects,
        where the keys are dataset ids and the values are the corresponding AnnData objects.
        If a single dataset is passed, an AnnData object is returned.
    """

    datasets = datasets or get_datasets(data_dir)

    if isinstance(datasets, DatasetDict):
        return DatasetDict(
            {
                dataset_id: read_dataset(
                    datasets[dataset_id], additional_readers=additional_readers
                )
                for dataset_id in sorted(datasets.keys())
            }
        )
    elif isinstance(datasets, DataSet):
        return read_dataset(datasets, additional_readers=additional_readers)
    else:
        raise TypeError(
            f'The type of "datasets" has to be a DatasetDict or a single DataSet. Use "fgread.get_datasets()" '
            f"to create it.\nSee {BLOGURL} for more information."
        )
