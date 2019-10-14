# coding: utf-8

"""Module for reading datasets shared on FASTGenomics"""
# set blog url for readme
BLOGURL = "https://www.fastgenomics.org/blog_posts/readers/"

from .read import get_datasets, read_dataset, read_datasets
from get_version import get_version


__version__ = get_version(__file__)
__author__ = "FASTGenomics"

del get_version
