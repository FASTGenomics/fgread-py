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
            f"id:     {self.id}\n"
            f"title:  {self.title}\n"
            f"format: {self.format}\n"
            f"path:   {self.path}\n"
            f'file:   {self.metadata["file"]}'
        )


class DatasetDict(dict):
    """
    Represents a dictionary for :py:class:`~DataSet` objects. You can select a single dataset by its ID (DatasetDict[ID]),
    or you can pass a list of IDs (DatasetDict[[ID1, ID3, ID4]]), or you can use slices (DatasetDict[1:3]).
    Note that lower and upper bounds are inclusive and you pass dataset IDs not indices (hence starting with 1).
    """

    def __getitem__(self, select):
        if isinstance(select, slice):
            stop = select.stop + 1 if select.stop else None
            newkey = slice(select.start, stop, select.step)
            keys = list(self.keys())
            keys.append(max(self.keys()) + 1)
            idx = list(range(max(keys))[newkey])
            select = list(set(idx) & set(keys))
        if isinstance(select, list):
            return DatasetDict({f"{sel}": self[sel] for sel in select})
        return dict.__getitem__(self, select)

    def __repr__(self):
        ds_list = [
            f"Dataset: {id}\n{indent_multiline(str(ds))}" for id, ds in self.items()
        ]
        return "\n\n".join(ds_list)


def indent_multiline(ml_str, tabs=1):
    """
    Indents a multiline string.
    :param ml_str: A multiline string
    :param tabs: the number of tabs (indents) to add to each line
    :return: An indented multiline string
    """
    return re.sub(r"^", "\t" * tabs, ml_str, flags=re.M)
