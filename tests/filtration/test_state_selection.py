from pytest import fixture

from shnitsel.filtering.state_selection import StateSelection
from shnitsel.io import read


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
def state_selection(data):
    return StateSelection.init_from_dataset(data)


class TestStateSelection:
    """Tests for the Datasheet utility class"""

    def test_is_data_loaded(self, data):
        assert data is not None

    def test_init_from_dataset(self, data):
        StateSelection.init_from_dataset(data)

    def test_init_from_descriptor(self, data):
        ssel = StateSelection.init_from_descriptor(['1, 1->2, 3<>1', '2<>1'])
        assert 1 in ssel.states
        assert 2 in ssel.states
        assert 3 in ssel.states
        assert (1, 2) in ssel.state_combinations
        assert (2, 1) not in ssel.state_combinations
        assert (1, 3) in ssel.state_combinations and (3, 1) in ssel.state_combinations
        assert (1, 2) in ssel.state_combinations and (2, 1) in ssel.state_combinations
        assert (3, 2) not in ssel.state_combinations and (
            2,
            3,
        ) not in ssel.state_combinations

    def test_only_singlets(self, state_selection: StateSelection):
        state_selection.singlets_only()

    def test_only_triplets(self, state_selection: StateSelection):
        state_selection.triplets_only()

    def test_state_charges(self, state_selection: StateSelection):
        state_selection.set_state_charges(1)

    def test_state_charges_2(self, state_selection: StateSelection):
        state_selection.set_state_charges({s: 10 for s in state_selection.states})

    def test_same_multiplicity_transition(self, state_selection: StateSelection):
        state_selection.same_multiplicity_transitions()

    def test_different_multiplicity_transitions(self, state_selection: StateSelection):
        state_selection.different_multiplicity_transitions()

    def test_auto_assign_colors(self, state_selection: StateSelection):
        state_selection.auto_assign_colors()

    def test_combination_info_iterator(self, state_selection: StateSelection):
        list(state_selection.combination_info())

    def test_state_info_iterator(self, state_selection: StateSelection):
        list(state_selection.state_info())

    def test_state_comb_name_assignment(self, state_selection: StateSelection):
        state_selection.set_state_combination_names({(1, 2): "party_time"})

    def test_state_name_assignment(self, state_selection: StateSelection):
        state_selection.set_state_names({1: "carrot"})

    def test_state_type_assignment(self, state_selection: StateSelection):
        state_selection.set_state_types({1: 3})

    def test_filter_states(self, state_selection: StateSelection):
        state_selection.filter_states(
            [2, 1, 3],
            exclude_ids=[3],
            charge=0,
            exclude_charge=1,
            multiplicity=1,
            exclude_multiplicity=3,
            min_states_in_selection=2,
        )

    def test_ground_state_transitions(self, state_selection: StateSelection):
        state_selection.ground_state_transitions()
        state_selection.excited_state_transitions()

    def test_get_state_color(self, state_selection: StateSelection):
        state_selection.get_state_color(1)

    def test_get_state_combination_color(self, state_selection: StateSelection):
        state_selection.get_state_combination_color((1, 2))

    def test_get_state_combination_name(self, state_selection: StateSelection):
        state_selection.get_state_combination_name_or_default((1, 2))
        state_selection.get_state_combination_tex_label((1, 2))

    def test_get_state_name(self, state_selection: StateSelection):
        state_selection.get_state_name_or_default(1)
        state_selection.get_state_tex_label(1)

    def test_get_state_degeneracy(self, state_selection: StateSelection):
        state_selection.get_state_degeneracy(1)
        state_selection.get_state_degeneracy(2)
        state_selection.get_state_degeneracy(3)

    def test_filter_state_combinations(self, state_selection: StateSelection):
        state_selection.filter_state_combinations(ids=[(1, 2), (2, 3)])

    def test_non_degenerate(self, state_selection: StateSelection):
        state_selection.non_degenerate()

    def test_has_state(self, state_selection: StateSelection):
        assert state_selection.has_state(1)
        assert not state_selection.has_state(-2)

    def test_has_state_combination(self, state_selection: StateSelection):
        assert state_selection.has_state_combination((1, 2))
