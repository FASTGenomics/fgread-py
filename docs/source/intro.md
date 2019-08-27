# Intro

This package implements convenience functions for loading data sets in the
[FASTGenomics][fg] [analysis][fg_analysis] environment.  The functions from this package
will let you list and load data sets for which the analysis was defined.

[fg]: https://beta.fastgenomics.org/webclient/
[fg_analysis]: https://beta.fastgenomics.org/webclient/searchPage/analyses

## Supported formats

The following formats are supported by this package
- [AnnData](https://github.com/theislab/anndata)
- [CellRanger (hdf5)](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices)
- [Drop-Seq (tsv)](https://github.com/Hoohm/dropSeqPipe/)
- [Loom](http://loompy.org/)

Currently unsupported
- [Seurat Object](https://satijalab.org/seurat/)

# Usage

Start by importing the module with

``` python
import fgread
```

## Listing data sets

To list the data sets simply call the `fgread.list_datasets` function

``` R
dsets_list = fgread.list_datasets()
```

The `dsets_list` would then contain the information about the location, format, title,
etc. about each data set.

```
{1: id: 1
 title: Loom data set
 format: Loom
 path: ../tests/data/readers/dataset_0001,
 2: id: 2
 title: AnnData data set
 format: AnnData
 path: ../tests/data/readers/dataset_0002
}
```

Note, that `fgread.list_datasets()` does not load any of the data sets.  It's purpose
is to get a list of available data sets, from which you can select the ones you would
like to load.

## Loading a single data set

To load a single data set use `fgread.read_dataset`.  The code below loads the first
data set from the list (the "Loom data set") and returns an [AnnData][anndata] object

``` R
adata = fgread.read_dataset(dsets_list[1])
```

To load the second data set simply run

``` R
adata = fgread.read_dataset(dsets_list[2])
```

The `fgread.read_dataset` function resolves the underlying format of the data set
automatically, based on the `format` attributes contained in the `dsets_list[1]`.

[anndata]: https://anndata.readthedocs.io/en/stable/

## Loading multiple data sets

Similarly, one can load multiple data sets with a single command: `fgread.read_datasets`
(note the `s` at the end).  The command loads all available data sets into _separate_
anndata objects and returns a list of these objects (where the indices correspond to the
indices from `fgread.list_datasets`).

``` R
dsets = fgread.read_datasets(dsets_list)
```
Now the `dsets` is a list containing two anndata objects

```
{1: AnnData object with n_obs × n_vars = 298 × 16892
 obs: 'Area', 'Cell_cluster', 'Cell_id'
 var: 'fg_title', 'fg_id'
 uns: 'metadata',
 2: AnnData object with n_obs × n_vars = 10 × 20
 obs: 'Area', 'Cell_cluster', 'Cell_id'
 var: 'fg_title', 'fg_id'
 uns: 'metadata'
}
```

Used without any arguments `fgread.read_datasets()` loads all data sets

``` R
dsets = fgread.read_datasets()
```


```
{1: AnnData object with n_obs × n_vars = 298 × 16892
 obs: 'Area', 'Cell_cluster', 'Cell_id'
 var: 'fg_title', 'fg_id'
 uns: 'metadata',
 2: AnnData object with n_obs × n_vars = 10 × 20
 obs: 'Area', 'Cell_cluster', 'Cell_id'
 var: 'fg_title', 'fg_id'
 uns: 'metadata'
}
```

# Known issues

Please report the issues through [github][issues].

[issues]: https://github.com/FASTGenomics/fgread-py/issues

# Development and testing

Clone the repository along with the test data by running

``` bash
git clone --recurse-submodules git@github.com:FASTGenomics/fgread-py.git
```

Then enter the `fgread-py` directory and install the dependencies with

``` bash
flit install --deps all
```

To test the package use

``` bash
python3 -m pytest
```
