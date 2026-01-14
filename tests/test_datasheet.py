from pytest import fixture

from shnitsel.io import read
from shnitsel.vis.datasheet import Datasheet


@fixture(
    params=[
        'tutorials/tut_data/traj_I02.nc',
        # 'tutorials/test_data/sharc/traj_I01_v3.0_triplets_nacs_socs',
        # 'tutorials/test_data/newtonx/test_I01_v2.6',
    ]
)
def data(request):
    res = read(request.param)
    return res


@fixture
def sheet(data):
    return Datasheet(data)


@fixture
def page(sheet):
    first_page_name = list(sheet.datasheet_pages)[0]
    return sheet.datasheet_pages[first_page_name]


class TestDatasheetFunctionality:
    """Tests for the Datasheet utility class
    """

    def test_is_data_loaded(self, data):
        assert data is not None

    def test_per_state_histograms(self, page):
        page.plot_per_state_histograms()

    def test_nacs_histograms(self, page):
        page.plot_nacs_histograms()

    def test_timeplots(self, page):
        page.plot_timeplots()

    def test_datasheet_full(self, page):
        page.plot()