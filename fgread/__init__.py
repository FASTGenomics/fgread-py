# coding: utf-8

"""Module for reading datasets shared on FASTGenomics"""
# set blog url for readme
import os

try:
    fgurl = os.environ['FG_URL'].rsplit(":", 1)[0]
except:
    fgurl = "https://beta.fastgenomics.org"
DOCSURL = fgurl + "/docs/"
DS_URL_PREFIX = fgurl + "/webclient/ui/#/datasets/detail-"

from .read import ds_info, load_data, get_datasets, read_dataset, read_datasets
from get_version import get_version

__version__ = get_version(__file__)
__author__ = "FASTGenomics"

del get_version
