# coding: utf-8

"""Module for reading datasets shared on FASTGenomics"""
import os
from .helpers import within_flit
from get_version import get_version

__version__ = get_version(__file__)
__author__ = "FASTGenomics"

del get_version

# set blog url for readme
try:
    fgurl = os.environ["FG_URL"].rsplit(":", 1)[0]
except:
    fgurl = "https://beta.fastgenomics.org"
DOCSURL = fgurl + "/docs/"
DS_URL_PREFIX = fgurl + "/datasets/detail-"

if not within_flit():
    from .read import ds_info, load_data
