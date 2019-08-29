# coding: utf-8

"""Module for reading datasets shared on FASTGenomics"""

from .read import print_datasets, get_datasets, read_dataset, read_datasets

from get_version import get_version

__version__ = get_version(__file__)
__author__ = "FASTGenomics"

del get_version
