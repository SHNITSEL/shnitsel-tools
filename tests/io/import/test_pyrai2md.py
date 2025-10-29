import shnitsel as sh
from shnitsel.io.pyrai2md import parse_pyrai2md


class TestPyrai2mdFunctionality:
    """Class to test functionality related to the pyrai2md file format"""

    def test_parse(self):
        parse_pyrai2md('tutorials/test_data/pyrai2md')
