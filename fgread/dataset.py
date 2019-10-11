import json
from pathlib import Path
import re

DATASET_INFO_FILE = "dataset_info.json"


class DataSet(object):
    """Represents a dataset on FASTGenomics, including the relative location and the
    contents of the ``metadata.json`` file.

    :param path: absolute path to a dataset folder, for example
        ``/fastgenomics/data/dataset_0001``
    """

    def __init__(self, path: str):
        self.path = Path(path)

        if not self.path.exists():
            raise FileNotFoundError(self.path)

        self.metadata = self.read_metadata()
        self.format = self.metadata["format"]
        self.title = self.metadata["title"]
        self.file = self.path / self.metadata["file"]
        self.id = int(self.path.name.split("_")[-1])

    def read_metadata(self):
        with open(self.path / DATASET_INFO_FILE) as f:
            return json.load(f)

    def __repr__(self):
        return (
            f'id:     {self.id}\n'
            f'title:  {self.title}\n'
            f'format: {self.format}\n'
            f'path:   {self.path}\n'
            f'file:   {self.metadata["file"]}'
        )

class DatasetDict(dict):
    '''
    Represents a dictionary for :py:class:`~DataSet` objects.
    '''

    def __getitem__(self, key):
        if isinstance(key, slice):
            stop = key.stop + 1 if key.stop else 0
            newkey = slice(key.start, stop, key.step)
            keys = list(self.keys())
            keys.append(max(self.keys()) + 1)
            idx = list(range(max(keys))[newkey])
            key = list(set(idx) & set(keys))
        if isinstance(key, list):
            return DatasetDict({k: self[k] for k in key})
        return dict.__getitem__(self, key)

    def __repr__(self):
        ds_list = [f"Dataset: {id}\n{indent_multiline(str(ds))}" for id, ds in self.items()]
        return "\n\n".join(ds_list)


def indent_multiline(ml_str, tabs=1):
    '''

    :param ml_str: A multiline string
    :param tabs: the number of tabs (indents) to add to each line
    :return: An indented multiline string
    '''
    return re.sub(r"^", "\t" * tabs, ml_str, flags=re.M)