import json
from pathlib import Path

DATASET_INFO_FILE = "dataset_info.json"


class DataSet(object):
    """Represents a data set on FASTGenomics, including the relative location and the
    contents of the ``metadata.json`` file.

    :param path: absolute path to a data set folder, for example
        ``/fastgenomics/data/dataset_0001``
    """

    def __init__(self, path: str):
        self.path = Path(path)

        if not self.path.exists():
            raise FileNotFoundError(filename=self.path)

        self.metadata = self.read_metadata()
        self.format = self.metadata["format"]
        self.title = self.metadata["title"]
        self.file = self.path / self.metadata["file"]
        self.id = int(self.path.name.split("_")[-1])

    def read_metadata(self):
        with open(self.path / DATASET_INFO_FILE) as f:
            return json.load(f)

    def __repr__(self):
        return "\n".join(
            [
                f"id: {self.id}",
                f"title: {self.title}",
                f"format: {self.format}",
                f"path: {self.path}",
            ]
        )
