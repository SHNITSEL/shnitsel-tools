import shnitsel as sh


class TestPyrai2mdFunctionality:
    """Class to test functionality related to the pyrai2md file format"""
    
    def test_parse():
        sh.parse.pyrai2md.read_traj('/traj/pyrai2md')