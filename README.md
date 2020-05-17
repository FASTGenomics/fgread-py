[![Documentation Status](https://readthedocs.org/projects/fgread-py/badge/?version=latest)](https://fgread-py.readthedocs.io/en/latest/?badge=latest)
[![Build Status](https://travis-ci.org/FASTGenomics/fgread-py.svg?branch=master)](https://travis-ci.org/FASTGenomics/fgread-py)
[![PyPI version](https://badge.fury.io/py/fgread.svg)](https://badge.fury.io/py/fgread)
[![PyPI download month](https://img.shields.io/pypi/dm/fgread.svg)](https://pypi.python.org/pypi/fgread/)

# FASTGenomics Reader Module for Python

This package implements convenience functions for loading datasets in the
[FASTGenomics][fg] [analysis][fg_analysis] environment. The functions from this package
will let you list and load datasets for which the analysis was defined.

[fg]: https://beta.fastgenomics.org/webclient/
[fg_analysis]: https://beta.fastgenomics.org/webclient/searchPage/analyses

## Documentation

For the general documentation on how to use the reader, please visit our [FASTGenomics Documentation](https://beta.fastgenomics.org/docs/).

For details on the available functions see the [API Documentation](https://fgread-py.readthedocs.io/en/stable/api.html).

## Known issues

Please report the issues through [github][issues].

[issues]: https://github.com/FASTGenomics/fgread-py/issues

## Development and testing

Clone the repository along with the test data by running

```bash
git clone --recurse-submodules git@github.com:FASTGenomics/fgread-py.git
```

Then enter the `fgread-py` directory and install the dependencies with

```bash
flit install --deps all
```

To test the package use

```bash
python3 -m pytest
```
