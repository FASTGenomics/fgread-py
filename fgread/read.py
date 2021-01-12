import json
import logging
import re
from pathlib import Path
from typing import Optional, Union

import pandas as pd

from . import DOCSURL, DS_URL_PREFIX, readers

# configure logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
logger.addHandler(ch)

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
DF_SORT_ORDER = [
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


def get_datasets_df(data_dir: Path = DATA_DIR) -> pd.DataFrame:
    """Constructs a :py:func:`pandas.DataFrame` from all available datasets.

    Parameters
    ----------
    data_dir : Path, optional
        Directory containing the datasets, e.g. ``fastgenomics/data``, by default DATA_DIR

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing all available datasets
    """

    ds_paths = get_ds_paths(data_dir=data_dir)

    ds_df = pd.DataFrame()
    for ds_path in ds_paths:
        with open(ds_path / "dataset_info.json") as f:
            info_df = json.load(f)
            info_df["path"] = str(ds_path)
            info_df["numberOfExpressionDataFiles"] = len(
                info_df["expressionDataFileInfos"]
            )
            info_df["numberOfMetaDataFiles"] = len(info_df["metaDataFileInfos"])
            _ = info_df.pop("schemaVersion", None)
        ds_df = ds_df.append(info_df, ignore_index=True)

    # sort colnames

    col_names = ds_df.columns.values.tolist()
    col_names_sorted = [name for name in DF_SORT_ORDER if name in col_names]
    [col_names.remove(name) for name in DF_SORT_ORDER if name in col_names]
    col_names_sorted.extend(col_names)
    ds_df = ds_df[col_names_sorted]

    # Format types
    ds_df = ds_df.astype(
        {
            "numberOfCells": "int32",
            "numberOfGenes": "int32",
            "numberOfExpressionDataFiles": "int32",
            "numberOfMetaDataFiles": "int32",
        }
    )

    return ds_df


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
        return

    try:
        ds_df = get_datasets_df(data_dir=data_dir)
    except NoDatasetsError as err:
        logger.warning(err)
        return pd.DataFrame()

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
            name for name in DF_SORT_ORDER if name in single_col_names
        ]
        [
            single_col_names.remove(name)
            for name in DF_SORT_ORDER
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
        if pretty:
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
    ds: Optional[str] = None,
    data_dir: Path = DATA_DIR,
    additional_readers: dict = {},
    expression_file: Optional[str] = None,
    as_format: Optional[str] = None,
):
    """This function loads a single dataset into an AnnData object.
    If there are multiple datasets available you need to specify one by setting
    ``ds`` to a dataset `id` or dataset `title`.
    To get an overview of availabe dataset use :py:func:`ds_info`

    Parameters
    ----------
    ds : str, optional
        A single dataset ID or dataset title to select a dataset to be loaded.
        If only one dataset is available you do not need to set this parameter, by default None
    data_dir : Path, optional
        Directory containing the datasets, e.g. ``fastgenomics/data``, by default DATA_DIR
    additional_readers : dict, optional
        Used to specify your own readers for the specific data set format.
        Dict key needs to be file extension (e.g., h5ad), dict value a function.
        Still experimental, by default {}
    expression_file: str, Optional
        The name of the expression file to load.
        Only needed when there are multiple expression files in a dataset.
    as_format: str, optional
        Specifies which reader should be uses for this dataset. Overwrites the auto-detection
        of the format. Possible parameters are the file extensions of our supported data
        formats: ``h5ad``, ``h5``, ``hdf5``, ``loom``, ``rds``, ``csv``, ``tsv``.

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
        single_df = select_ds_id(ds, df=get_datasets_df(data_dir=data_dir))
    else:
        single_df = get_datasets_df(data_dir=data_dir)
        if len(single_df) > 1:
            raise RuntimeError(
                "There is more than one dataset available in this analysis. "
                "Please select one by its ID or title. "
                'You can list available datasets by using "fgread.ds_info()".'
            )

    exp_count = single_df.loc[0, "numberOfExpressionDataFiles"]
    meta_count = single_df.loc[0, "numberOfMetaDataFiles"]

    if exp_count == 0:
        raise TypeError(
            f"There is no expression data available in this data set.\n"
            f"Metadata files: {meta_count}."
        )

    exp_files = [exp["name"] for exp in single_df.loc[0, "expressionDataFileInfos"]]

    if expression_file:
        if expression_file in exp_files:
            file = expression_file
        else:
            raise KeyError(
                f'Expression file "{expression_file}" not found in dataset. '
                f"Available expression files are: {exp_files}."
            )
    else:
        if exp_count == 1:
            file = single_df.loc[0, "expressionDataFileInfos"][0]["name"]
        else:
            raise TypeError(
                f"There are {exp_count} expression data files in this dataset. "
                'Please specify which one you want to load using the parameter "expression_file". '
                f"Available expression files are: {exp_files}."
            )

    title = single_df.loc[0, "title"]
    ds_id = single_df.loc[0, "id"]
    path = single_df.loc[0, "path"]

    metadata_keys_delete = [
        "state",
        "expressionDataFileInfos",
        "metaDataFileInfos",
    ]  # These keys do not go into anndata.uns
    metadata_dict = single_df.loc[0].to_dict()
    for key in metadata_keys_delete:
        try:
            del metadata_dict[key]
        except KeyError:
            logger.warning(
                f"Key {key} can not removed from metadata dict as it is not present."
            )

    if as_format:
        format = as_format.lower()
    else:
        try:
            format = file.rsplit(".", 1)[1].lower()
            logger.info(f'Expression file "{file}" with format "{format}".')
        except ValueError as e:
            raise ValueError(
                f'The expression file "{file}" has no valid file suffix.'
            ).with_traceback(e.__traceback__)

    if format in readers:
        if meta_count != 0:
            logger.info(
                f"There are {meta_count} metadata files in this dataset. "
                "This data will not be integrated into the anndata object."
            )
        logger.info(
            f'Loading file "{file}" from dataset "{title}" in format "{format}" from directory "{path}"...\n'
        )
        adata = readers[format](Path(path) / file)
        adata.uns["ds_metadata"] = {ds_id: metadata_dict}
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
            f'Unsupported file format "{format}", use one of {list(readers)}. '
            f'You can force the usage of a specific reader by setting "as_format" to a supported format. '
            f"In addition, you can also implement your own reading function. See {DOCSURL} for more information."
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
    elif len_df == 0:
        add_err = ""
        if not ds.startswith("dataset-"):
            add_err = " Please note that dataset titles can be changed by the owner. To be safe, you might want to consider dataset IDs instead."
        raise KeyError("Your selection matches no datasets." + add_err)
    else:
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
        raise NoDatasetsError(
            f'There are no datasets attached to this analysis. Path "{data_dir}" does not exist.'
        )

    paths = [
        Path(subdir)
        for subdir in sorted(data_dir.iterdir())
        if subdir.is_dir() and re.match(r"^dataset_\d{4}$", subdir.name)
    ]

    if not paths:
        raise NoDatasetsError(
            f'There are no datasets attached to this analysis. Path "{data_dir}" is empty.'
        )

    return paths


class NoDatasetsError(Exception):
    """Raised when no datasets are attached"""

    pass
