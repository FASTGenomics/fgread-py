import re
from pathlib import Path
from . import readers

DEFAULT_PARSERS = {
    "Loom": readers.read_loom,
    "Seurat Object": readers.read_seurat,
    "AnnData": readers.read_anndata,
    "10x h5": readers.read_10x_hdf5,
    "Drop-Seq": readers.read_dropseq,
}

DATA_DIR = Path("/fastgenomics/data")


def read_data_set(dataset_dir, parsers={}):
    dataset_dir = Path(dataset_dir)
    manifest = readers.read_manifest(dataset_dir)
    format = manifest["format"]
    title = manifest["title"]

    parsers = {**DEFAULT_PARSERS, **parsers}

    if format in parsers:
        print(
            f'Loading data set "{title}" in format "{format} from directory "{dataset_dir}".'
        )
        return parsers[format](dataset_dir)
    else:
        raise KeyError(f'Unsupported format "{format}", use one of {parsers}')


def list_dirs(data_dir=DATA_DIR):
    data_dir = Path(data_dir)
    return {
        int(f.name.split("_")[1]): f
        for f in data_dir.iterdir()
        if f.is_dir() and re.match(r"^dataset_\d{4}$", f.name)
    }


def list_data_sets(data_dir=DATA_DIR):
    data_dir = Path(data_dir)
    dirs = list_dirs(data_dir)
    return {
        id: dict(path=dir, manifest=readers.read_manifest(dir))
        for id, dir in dirs.items()
    }


def read_data_sets(data_dir=DATA_DIR, data_sets=None, parsers={}):
    data_dir = Path(data_dir)
    data_sets = data_sets or list_data_sets(data_dir)
    return {id: read_data_set(dset["path"]) for id, dset in data_sets.items()}
