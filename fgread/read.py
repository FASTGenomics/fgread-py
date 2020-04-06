import json
import logging
import re
from pathlib import Path
from typing import Optional, Union

import pandas as pd
from deprecated.sphinx import deprecated

from . import DOCSURL, DS_URL_PREFIX, readers, readers_old
from .dataset import DataSet, DatasetDict

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
    "loom": readers.read_loom_to_anndata,
    "rds": readers.read_seurat_to_anndata,
    "h5ad": readers.read_anndata_to_anndata,
    "hdf5": readers.read_10xhdf5_to_anndata,
    "h5": readers.read_10xhdf5_to_anndata,
    "tsv": readers.read_densetsv_to_anndata,
    "csv": readers.read_densecsv_to_anndata,
}


DATA_DIR = Path("/fastgenomics/data")


def ds_info(
    ds: Optional[str] = None,
    pretty: bool = None,
    output: bool = None,
    data_dir: Path = DATA_DIR,
) -> pd.DataFrame:
    """Get information on all available datasets in this analysis.
    
    Parameters
    ----------
    ds : Optional[str], optional
        A single dataset ID or dataset title. If set, only this dataset will be displayed. Recommended to use with ``pretty``, by default None
    pretty : bool, optional
        Whether to display some nicely formatted output, by default True
    output : bool, optional
        Whether to return a DataFrame or not, by default True
    data_dir : Path, optional
        Directory containing the datasets, e.g. ``fastgenomics/data``, by default DATA_DIR
    
    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing all, or a single dataset (depends on ``ds``)
    """

    if pretty is None:
        pretty = ds is not None
    if output is None:
        output = ds is None

    if not pretty and not output:
        logger.warning(
            'You have set "pretty" and "output" to false. Hence, this function will do/return nothing.'
        )
    ds_paths = get_ds_paths(data_dir=data_dir)
    ds_df = pd.DataFrame()
    for ds_path in ds_paths:
        with open(ds_path / "dataset_info.json") as f:
            info_df = json.load(f)
            info_df["path"] = ds_path
            info_df["numberOfExpressionDataFiles"] = len(
                info_df["expressionDataFileInfos"]
            )
            info_df["numberOfMetaDataFiles"] = len(info_df["metaDataFileInfos"])
            _ = info_df.pop("schemaVersion", None)
        ds_df = ds_df.append(info_df, ignore_index=True)

    # sort colnames
    sort_order = [
        "title",
        "id",
        "organism",
        "tissue",
        "numberOfCells",
        "numberOfGenes",
        "path",
        "numberOfExpressionDataFiles",
        "expressionDataFileNames",
        "numberOfMetaDataFiles",
        "metaDataFileNames",
        "expressionDataFileInfos",
        "metaDataFileInfos",
    ]
    col_names = ds_df.columns.values.tolist()
    col_names_sorted = [name for name in sort_order if name in col_names]
    [col_names.remove(name) for name in sort_order if name in col_names]
    col_names_sorted.extend(col_names)
    ds_df = ds_df[col_names_sorted]

    ds_df = ds_df.astype(
        {
            "numberOfCells": "int32",
            "numberOfGenes": "int32",
            "numberOfExpressionDataFiles": "int32",
            "numberOfMetaDataFiles": "int32",
        }
    )

    def add_url(title, id):
        return f'<a href="{DS_URL_PREFIX}{id}" target="_blank">{title}</a>'

    def disp_pretty_df(df, index=True, header=True):
        try:
            from IPython.display import display, Markdown

            df_html = df.to_html(
                render_links=True,
                escape=False,
                header=header,
                index=index,
                justify="center",
            )
            display(Markdown(df_html))
        except:
            logger.warning(
                "IPython not available. Pretty printing only works in Jupyter Notebooks."
            )

    if ds:
        if ds_df.empty:
            raise ValueError("There are no datasets in your analysis")

        single_ds_df = select_ds_id(ds, df=ds_df)

        single_ds_df["expressionDataFileNames"] = ", ".join(
            [expr["name"] for expr in single_ds_df.loc[0, "expressionDataFileInfos"]]
        )

        single_ds_df["metaDataFileNames"] = ", ".join(
            [expr["name"] for expr in single_ds_df.loc[0, "metaDataFileInfos"]]
        )

        # Sort columns
        single_col_names = single_ds_df.columns.values.tolist()
        single_col_names_sorted = [
            name for name in sort_order if name in single_col_names
        ]
        [
            single_col_names.remove(name)
            for name in sort_order
            if name in single_col_names
        ]
        single_col_names_sorted.extend(single_col_names)
        single_ds_df = single_ds_df[single_col_names_sorted]

        if pretty:
            pretty_df = single_ds_df

            pretty_df["expressionDataFileNames"] = "<br>".join(
                [expr["name"] for expr in pretty_df.loc[0, "expressionDataFileInfos"]]
            )

            pretty_df["metaDataFileNames"] = ", ".join(
                [expr["name"] for expr in pretty_df.loc[0, "metaDataFileInfos"]]
            )

            empty_cols = [
                col for col in pretty_df.columns if pretty_df.loc[0, col] == ""
            ]
            pretty_df = pretty_df.drop(
                labels=["expressionDataFileInfos", "metaDataFileInfos"] + empty_cols,
                axis=1,
                errors="ignore",
            )

            pretty_df.loc[0, "title"] = pretty_df.apply(
                lambda x: add_url(x.title, x.id), axis=1
            ).squeeze()
            disp_pretty_df(pretty_df.T, header=False)

        if output:
            return single_ds_df

    else:
        if pretty and not ds_df.empty:
            pretty_df = ds_df.drop(
                labels=[
                    "description",
                    "license",
                    "preprocessing",
                    "citation",
                    "webLink",
                    "file",
                    "expressionDataFileInfos",
                    "metaDataFileInfos",
                ],
                axis=1,
                errors="ignore",
            )
            pretty_df["title"] = pretty_df.apply(
                lambda x: add_url(x.title, x.id), axis=1
            )
            disp_pretty_df(pretty_df)

        if output:
            return ds_df


def load_data(
    ds: Optional[str] = None, data_dir: Path = DATA_DIR, additional_readers: dict = {}
):
    """This function loads a single dataset into an AnnData object.
    If there are multiple datasets available you need to specify one by setting
    ``ds`` to a dataset `id` or dataset `title`.
    To get an overview of availabe dataset use :py:func:`ds_info`

    Parameters
    ----------
    ds : Optional[str], optional
        A single dataset ID or dataset title to select a dataset to be loaded.
        If only one dataset is available you do not need to set this parameter, by default None
    data_dir : Path, optional
        Directory containing the datasets, e.g. ``fastgenomics/data``, by default DATA_DIR
    additional_readers : dict, optional
        Used to specify your own readers for the specific data set format.
        Dict key needs to be file extension (e.g., h5ad), dict value a function.
        Still experimental, by default {}

    Returns
    -------
    AnnData Object
        A single AnnData object with dataset id in `obs` and all dataset metadata in `uns`
    
    Examples
    --------
    To use a custom reader for files with the extension ".fg", you have to define a function first:

    >>> def my_loader(file):
    ...     anndata = magic_file_loading(file)
    ...     return anndata

    You can then use this reader like this:
    
    >>> fgread.load_data("my_dataset", additional_readers={"fg": my_loader})

    """
    readers = {**DEFAULT_READERS, **additional_readers}

    if ds:
        single_df = select_ds_id(ds, df=ds_info(data_dir=data_dir, pretty=False))
    else:
        single_df = ds_info(data_dir=data_dir, pretty=False)
        assert (
            len(single_df) == 1
        ), f"There is more than one dataset available. Please select one by its ID or title."

    exp_count = single_df.loc[0, "numberOfExpressionDataFiles"]
    meta_count = single_df.loc[0, "numberOfMetaDataFiles"]
    if exp_count == 0:
        raise TypeError(
            f"There is no expression data available in this data set.\n"
            f"Metadata files: {meta_count}."
        )
    elif exp_count >= 2:
        raise TypeError(
            f"There are {exp_count} expression data files and {meta_count} metadata files in this dataset. "
            "Currently we only provide reading functionality for one expression data file and ignore meta data. "
            "Please load the required data yourself."
        )

    title = single_df.loc[0, "title"]
    ds_id = single_df.loc[0, "id"]
    file = single_df.loc[0, "expressionDataFileInfos"][0]["name"]
    path = single_df.loc[0, "path"]

    try:
        _, format = file.rsplit(".", 1)
        logger.info(f'Expression file "{file}" with format "{format}".')
    except ValueError as e:
        raise ValueError(
            f'The expression file "{file}" has no valid file ending.'
        ).with_traceback(e.__traceback__)

    if format in readers:
        if meta_count != 0:
            logger.info(
                f"There are {meta_count} metadata files in this dataset. "
                "This data will not be integrated into the anndata object."
            )
        logger.info(
            f'Loading dataset "{title}" in format "{format}" from directory "{path}"...\n'
        )
        adata = readers[format](Path(path) / file)
        adata.uns["ds_metadata"] = {ds_id: single_df.loc[0].to_dict()}
        adata.obs["fg_id"] = ds_id
        n_genes = adata.shape[1]
        n_cells = adata.shape[0]
        logger.info(
            f'Loaded dataset "{title}" with {n_cells} cells and {n_genes} genes.\n'
            f"==================================================================\n"
        )
        return adata
    else:
        raise KeyError(
            f'Unsupported file format "{format}", use one of {list(readers)} or implement your '
            f"own reading function. See {DOCSURL} for more information."
        )


def select_ds_id(ds: str, df: pd.DataFrame = None) -> pd.DataFrame:
    """Select a single dataset from a pandas DataFrame by its ID or title
    
    Parameters
    ----------
    ds : str
        A single dataset ID or dataset title for selection
    df : pd.DataFrame, optional
        A pandas DataFrame from which a single entry is selected, by default None
    
    Returns
    -------
    pd.DataFrame
        A pandas DataFrame with only the selected dataset.
    """
    single_df = df.loc[(df["id"] == ds) | (df["title"] == ds)].reset_index(drop=True)
    len_df = len(single_df)
    if len_df == 1:
        return single_df.copy()
    else:
        if len_df > 1:
            display(single_df)
        raise KeyError(
            f"Your selection matches {len_df} datasets. Please make sure to select exactly one."
        )


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
    if not data_dir.exists():
        logger.warning("There are no datasets attached to this analysis.")
        return []

    paths = [
        subdir
        for subdir in sorted(data_dir.iterdir())
        if subdir.is_dir() and re.match(r"^dataset_\d{4}$", subdir.name)
    ]
    assert paths != [], "There is no data available in this analysis."
    return paths


@deprecated(
    version="0.4.0",
    reason="Please use the function `ds_info` or `load_data` instead. This function will be removed in the future.",
    category=FutureWarning,
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


@deprecated(
    version="0.4.0",
    reason="Please use the function `load_data` instead. This function will be removed in the future.",
    category=FutureWarning,
)
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
            f'The format of the dataset "{title}" is "{format}". Datasets with the "{format}" format are '
            f"unsupported by this module and have to be loaded manually.\nSee {DOCSURL} for more information."
        )
    elif format == "Not set":
        raise ValueError(
            f'The format of the dataset "{title}" was not defined. If you can modify the dataset please specify '
            f"its format in its details page, otherwise ask the dataset owner to do that.\nSee {DOCSURL} for more information."
        )
    else:
        raise KeyError(
            f'Unsupported format "{format}", use one of {list(readers)} or implement your '
            f"own reading function.\nSee {DOCSURL} for more information."
        )


@deprecated(
    version="0.4.0",
    reason="Please use the function `load_data` instead. This function will be removed in the future.",
    category=FutureWarning,
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
            f"to create it.\nSee {DOCSURL} for more information."
        )
