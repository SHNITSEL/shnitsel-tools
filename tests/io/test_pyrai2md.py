import shnitsel as sh


class TestPyrai2mdFunctionality:
    """Class to test functionality related to the pyrai2md file format"""
    
    def test_parse(self):
        sh.parse.pyrai2md.read_traj('tutorials/test_data/pyrai2md')