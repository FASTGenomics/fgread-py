import pytest
import fgread


def test_other(data_dir):
    with pytest.raises(
            NotImplementedError, match='The format of the dataset "Other dataset" is "Other". '
                                       'Datasets with the "Other" format are unsupported by '
                                       'this module and have to be loaded manually.\nSee '
                                       'https://beta.fastgenomics.org/docs/ for more information.'
    ):
        fgread.load_data("Other dataset",
                         data_dir=data_dir)


def test_notset(data_dir):
    with pytest.raises(
            ValueError, match='The format of the dataset "Not set dataset" was not defined. '
                              'If you can modify the dataset please specify its format '
                              'in its details page, otherwise ask the dataset owner to '
                              'do that.\nSee https://beta.fastgenomics.org/docs/ for '
                              'more information.'
    ):
        fgread.load_data("Not set dataset",
                         data_dir=data_dir)
