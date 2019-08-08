# coding: utf-8

"""Module for reading files shared on FASTGenomics"""

from .read import list_data_sets, read_data_set, read_data_sets

from get_version import get_version

__version__ = get_version(__file__)
__author__ = "Paweł Biernat"

del get_version
