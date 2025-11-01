import shnitsel as sh
from shnitsel.io.pyrai2md import parse_pyrai2md
from shnitsel.io import read


class TestPyrai2mdFunctionality:
    """Class to test functionality related to the pyrai2md file format"""

    def test_parse_direct(self):
        parse_pyrai2md('tutorials/test_data/pyrai2md/traj')

    def test_parse_wrapper_single(self):
        assert read('tutorials/test_data/pyrai2md/traj',
                    kind='pyrai2md') is not None

    def test_parse_wrapper_single_guess_kind(self):
        assert read('tutorials/test_data/pyrai2md/traj') is not None

    def test_parse_wrapper_trajectory(self):
        assert read('tutorials/test_data/pyrai2md/',
                    kind='pyrai2md') is not None

    def test_parse_wrapper_trajectory_guess_kind(self):
        assert read('tutorials/test_data/pyrai2md/') is not None
