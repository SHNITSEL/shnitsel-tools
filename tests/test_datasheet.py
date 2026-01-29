from pytest import fixture

from shnitsel.io import read
from shnitsel.vis.datasheet import Datasheet


@fixture(
    params=[
        ('tutorials/test_data/shnitsel/traj_I02.nc', 1),
        # 'tutorials/test_data/sharc/traj_I01_v3.0_triplets_nacs_socs',, 1),
        # 'tutorials/test_data/newtonx/test_I01_v2.6', 1),
    ]
)
def data(request):
    path, charge = request.param
    res = read(path).set_charge(charge)
    return res


@fixture
def datasheet(data):
    return Datasheet(data)


@fixture
def datasheet_with_structure(data):
    from shnitsel.filtering.structure_selection import StructureSelection

    features = StructureSelection.init_from_dataset(
        next(data.collect_data()), ['bonds', 'dihedrals']
    )
    return Datasheet(data, structure_selection=features)


@fixture
def datasheet_with_states(data):
    from shnitsel.filtering.state_selection import StateSelection

    states = StateSelection.init_from_dataset(
        next(data.collect_data()),
    ).singlets_only()
    return Datasheet(data, state_selection=states)


@fixture
def datasheet_page(datasheet):
    return list(datasheet.pages.values())[0]


class TestDatasheetFunctionality:
    """Tests for the Datasheet utility class"""

    def test_is_data_loaded(self, data):
        assert data is not None

    def test_datasheet_from_file(self):
        Datasheet('tutorials/test_data/shnitsel/traj_I02.nc')

    def test_per_state_histograms(self, datasheet_page):
        datasheet_page.plot_per_state_histograms(
            state_selection=datasheet_page.state_selection
        )

    def test_nacs_histograms(self, datasheet_page):
        datasheet_page.plot_nacs_histograms(
            state_selection=datasheet_page.state_selection
        )

    def test_timeplots(self, datasheet_page):
        datasheet_page.plot_timeplots(state_selection=datasheet_page.state_selection)

    def test_datasheet_page_plot(self, datasheet_page):
        datasheet_page.plot()

    def test_per_state_histograms_full(self, datasheet):
        datasheet.plot_per_state_histograms()

    def test_nacs_histograms_full(self, datasheet):
        datasheet.plot_nacs_histograms()

    def test_timeplots_full(self, datasheet):
        datasheet.plot_timeplots()

    def test_datasheet_full(self, datasheet):
        datasheet.plot()

    def test_datasheet_full_with_states(self, datasheet_with_states):
        datasheet_with_states.plot()

    def test_datasheet_full_with_features(self, datasheet_with_structure):
        datasheet_with_structure.plot()
