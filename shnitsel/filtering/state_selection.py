from dataclasses import dataclass
from itertools import combinations, permutations
import logging
import re
from typing import Iterable, Self, Sequence, Literal, TypeAlias, overload

import numpy as np
import xarray as xr

from shnitsel.data.dataset_containers import Frames, Trajectory
from shnitsel.data.dataset_containers.shared import ShnitselDataset
from shnitsel.data.state_helpers import state_name_to_tex_label
from ..core.typedefs import (
    StateCombination,
    StateId,
    StateInfo,
    StateCombInfo,
    MultiplicityLabel,
    MultiplicityLabelValues,
)

StateSelectionDescriptor: TypeAlias = (
    Sequence[StateId | StateCombination] | Sequence[str] | str
)

_state_combs_pattern = re.compile(
    r"(?P<state_from>.+)\s*(?P<rel>(->)|(<>))\s*(?P<state_to>.+)"
)
_state_id_pattern = re.compile(r"(?P<state>\d+)")
_state_name_pattern = re.compile(r"(?P<state_name>.+)")
_state_mult_pattern = re.compile(
    r"(?P<state_mult>" + "|".join(MultiplicityLabelValues) + ")"
)


@dataclass
class StateSelection:
    """Class to keep track of a (sub-)selection of states and state transitions for analysis and plotting."""

    is_directed: bool
    states_base: Sequence[StateId]
    states: Sequence[StateId]
    ground_state_id: StateId
    state_types: dict[StateId, int] | None
    state_names: dict[StateId, str] | None
    state_charges: dict[StateId, int] | None

    state_degeneracy_group: dict[StateId, int] | None
    degeneracy_group_states: dict[int, list[StateId]] | None
    # state_magnetic_nums: dict[StateId, float] | None

    state_combinations_base: list[StateCombination]
    # TODO: FIXME: This needs a property wrapper to return both directions if non-directional
    state_combinations: list[StateCombination]
    state_combination_names: dict[StateCombination, str] | None

    state_colors: dict[StateId, str] | None = None
    state_combination_colors: dict[StateCombination, str] | None = None

    def copy_or_update(
        self,
        *,
        is_directed: bool | None = None,
        ground_state_id: StateId | None = None,
        states: Sequence[StateId] | None = None,
        state_types: dict[StateId, int] | None = None,
        state_names: dict[StateId, str] | None = None,
        state_charges: dict[StateId, int] | None = None,
        state_combinations: list[StateCombination] | None = None,
        state_combination_names: dict[StateCombination, str] | None = None,
        state_colors: dict[StateId, str] | None = None,
        state_combination_colors: dict[StateCombination, str] | None = None,
        state_degeneracy_group: dict[StateId, int] | None = None,
        degeneracy_group_states: dict[int, list[StateId]] | None = None,
        inplace: bool = False,
    ) -> Self:
        """Function to create a copy with replaced member values.

        Meant as a helper for the `Frozen` logic of the selection, i.e. method calls return a new instance
        instead of updating the existing instance.

        Parameters
        ----------
        is_directed : bool, optional
            Potentially new flag whether state combinations are considered to be directed (i.e. (1,3) is different from (3,1)) or not. Defaults to None.
        ground_state_id : StateId or None, optional
            Potentially new ground_state id. Defaults to None.
        states : Sequence[StateId] or None, optional
            Potentially new state ids. Defaults to None.
        state_types : dict[StateId, int] or None, optional
            Potentially new state types/multiplicities. Defaults to None.
        state_names : dict[StateId, str] or None, optional
            Potentially new state names. Defaults to None.
        state_degeneracy_group : dict[StateId, int] or None, optional
            Optional degeneracy group indices for states. Defaults to None.
        degeneracy_group_states : dict[int, list[StateId]] or None, optional
            Optional mapping of degeneracy groups to the states they hold. Defaults to None.
        state_charges : dict[StateId, int] or None, optional
            Potentially new state charges. Defaults to None.
        state_combinations : list[StateCombination] or None, optional)
            Potentially new state combinations. Defaults to None.
        state_combination_names : dict[StateCombination, str] or None, optional
            Potentially new names for state combinations. Defaults to None.
        inplace : bool, optional
            A flag whether the existing instance should be updated or a new one should be created. Defaults to False, i.e. a new instance is created.
        state_colors : dict[StateId, str] or None, optional
            An optional colormap for states. Defaults to None.
        state_combination_colors : dict[StateCombination, str] or None, optional
            An optional colormap for state combinations. Defaults to None.

        Returns
        -------
        StateSelection
            The selection update with the new members set. Can either be a copy if `inplace=False` or the old instance with updated members otherwise.
        """
        if inplace:
            # Update and create
            if is_directed is not None:
                self.is_directed = is_directed
            if ground_state_id is not None:
                self.ground_state_id = ground_state_id
            if states is not None:
                self.states = states
            if state_types is not None:
                self.state_types = state_types
            if state_names is not None:
                self.state_names = state_names
            if state_charges is not None:
                self.state_charges = state_charges
            if state_combinations is not None:
                self.state_combinations = state_combinations
            if state_combination_names is not None:
                self.state_combination_names = state_combination_names
            elif state_names and state_combinations:
                state_combination_names = {}
                for comb in state_combinations:
                    first, second = comb

                    if first in state_names and second in state_names:
                        state_combination_names[comb] = (
                            f"{state_names[first]} - {state_names[second]}"
                        )
                    else:
                        logging.warning(
                            f"Could not assign name to state combination {comb} because of missing state names for {first} or {second}."
                        )
                self.state_combination_names = state_combination_names
            if state_colors is not None:
                self.state_colors = state_colors
            if state_combination_colors is not None:
                self.state_combination_colors = state_combination_colors

            if state_degeneracy_group is not None:
                self.state_degeneracy_group = state_degeneracy_group
            if degeneracy_group_states is not None:
                self.degeneracy_group_states = degeneracy_group_states

            return self
        else:
            if is_directed is None:
                is_directed = self.is_directed
            if ground_state_id is None:
                ground_state_id = self.ground_state_id
            if states is None:
                states = self.states
            if state_types is None:
                state_types = self.state_types
            if state_names is None:
                state_names = self.state_names
            if state_charges is None:
                state_charges = self.state_charges
            if state_combinations is None:
                state_combinations = self.state_combinations
            if state_combination_names is None:
                state_combination_names = self.state_combination_names
            elif state_names and state_combinations:
                state_combination_names = {}
                for comb in state_combinations:
                    first, second = comb

                    if first in state_names and second in state_names:
                        state_combination_names[comb] = (
                            f"{state_names[first]} - {state_names[second]}"
                        )
                    else:
                        logging.warning(
                            f"Could not assign name to state combination {comb} because of missing state names for {first} or {second}."
                        )
            if state_colors is None:
                state_colors = self.state_colors
            if state_combination_colors is None:
                state_combination_colors = self.state_combination_colors
            if state_degeneracy_group is None:
                state_degeneracy_group = self.state_degeneracy_group
            if degeneracy_group_states is None:
                degeneracy_group_states = self.degeneracy_group_states

            return type(self)(
                is_directed=is_directed,
                states_base=self.states_base,
                states=states,
                ground_state_id=ground_state_id,
                state_types=state_types,
                state_names=state_names,
                state_charges=state_charges,
                state_combinations_base=self.state_combinations_base,
                state_combinations=state_combinations,
                state_combination_names=state_combination_names,
                state_colors=state_colors,
                state_combination_colors=state_combination_colors,
                state_degeneracy_group=state_degeneracy_group,
                degeneracy_group_states=degeneracy_group_states,
            )

    @classmethod
    def init_from_dataset(
        cls: type[Self],
        dataset: xr.Dataset | ShnitselDataset,
        is_directed: bool = False,
    ) -> Self:
        """Alternative constructor that creates an initial StateSelection object from a dataset using the entire state information in it.

        Parameters
        ----------
        cls : type[StateSelection]
            The type of this StateSelection so that we can create instances of it.
        dataset : xr.Dataset or ShnitselDataset
            The dataset to extract the state information out of. Must have a `state` dimension and preferrably coordinates `state`, `state_names`, `state_types`, `state_charges`, and `statecomb` set.
            If `state` is not set as a coordinate, a potential dimension size of `state` is taken and states are enumerates `1` through `1+dataset.sizes['state']`.
            If `statecomb` is not set as a coordinate, all unordered pairs of states will be used as a default value for `state_combinations`.
        is_directed : bool, default=False
            Flag whether state combinatons should be assumed different, i.e. (1,3) should be considered different from (3,1).


        Returns
        -------
        StateSelection
            A state selection object initially covering all states (and state combinations) present in the dataset.


        Raises
        ------
        ValueError
            If no `state` information could be extracted from the dataset


        """
        assert 'state' in dataset.sizes, (
            "No state information on the provided dataset. Cannot initialize state selection."
        )

        if 'states' in dataset.coords:
            states = list(int(n) for n in dataset.coords['states'].values)
        elif 'state' in dataset.sizes:
            states = list(
                int(n) for n in np.arange(1, 1 + dataset.sizes['state'], dtype=StateId)
            )
        else:
            raise ValueError(
                "No sufficient state information on the provided dataset. Cannot initialize state selection."
            )
        ground_state_id = np.min(states)

        if 'state_types' in dataset.coords:
            state_types = {
                state_id: type_val
                for (state_id, type_val) in zip(
                    states, dataset.coords['state_types'].values
                )
            }
        else:
            logging.warning(
                "No state types available on the dataset. Please assign them yourself."
            )
            state_types = None

        if 'state_names' in dataset.coords:
            state_names = {
                state_id: name_val
                for (state_id, name_val) in zip(
                    states, dataset.coords['state_names'].values
                )
            }
        else:
            logging.warning(
                "No state names available on the dataset. Please assign them yourself."
            )
            state_names = None

        if 'state_charges' in dataset.coords:
            state_charges = {
                state_id: charge_val
                for (state_id, charge_val) in zip(
                    states, dataset.coords['state_charges'].values
                )
            }
        else:
            logging.info(
                "No state charges available on the dataset. Please assign them yourself."
            )
            state_charges = None

        if not is_directed:
            if 'statecomb' in dataset.coords:
                state_combinations = list(dataset.coords['statecomb'].values)
            else:
                state_combinations = list(combinations(states, 2))
        else:
            if 'full_statecomb' in dataset.coords:
                state_combinations = list(dataset.coords['full_statecomb'].values)
            elif 'statecomb' in dataset.coords:
                # Gather both directions of pairs initially.
                state_combinations = list(
                    set(dataset.coords['statecomb'].values).union(
                        (y, x) for x, y in dataset.coords['statecomb'].values
                    )
                )
            else:
                state_combinations = list(permutations(states, 2))

        if state_names is not None:
            state_combination_names = {
                (a, b): f"{state_names[a]} - {state_names[b]}"
                for (a, b) in state_combinations
            }
        else:
            state_combination_names = None

        state_degeneracy_group = {}
        degeneracy_group_states = {}
        if 'state_degeneracy_group' in dataset.coords:
            # print('In ds:', dataset['state_degeneracy_group'])
            # print("State degeneracy data from dataset")
            degeneracy_info: list[tuple[StateId, int]] = list(
                zip(states, dataset.state_degeneracy_group.values)
            )
            for state, deg_group in degeneracy_info:
                state_degeneracy_group[state] = int(deg_group)
                if int(deg_group) not in degeneracy_group_states:
                    degeneracy_group_states[int(deg_group)] = []
                degeneracy_group_states[int(deg_group)].append(state)
        else:
            # TODO: FIXME: Get degeneracy, magn number from the energy, nac, soc, dip_trans, etc.
            state_degeneracy_group = None
            degeneracy_group_states = None

        # print(state_degeneracy_group, degeneracy_group_states)

        # print('Init:', state_degeneracy_group)
        # print('Init:', degeneracy_group_states)
        # Create an initial state selection
        return cls(
            is_directed=is_directed,
            states_base=states,
            states=states,
            ground_state_id=ground_state_id,
            state_types=state_types,
            state_charges=state_charges,
            state_names=state_names,
            state_combinations_base=state_combinations,
            state_combinations=state_combinations,
            state_combination_names=state_combination_names,
            state_degeneracy_group=state_degeneracy_group,
            degeneracy_group_states=degeneracy_group_states,
        )

    @classmethod
    def init_from_descriptor(
        cls: type[Self],
        spec: StateSelectionDescriptor,
        is_directed: bool | None = None,
    ) -> Self:
        """Build a (rather rudimentary) state selection and
        state combination selection from descriptors with no support for determination of
        multiplicity groups and others but to simplify the process of providing a state selection
        to function calls.

        Parameters
        ----------
        cls : type[Self]
            StateSelection class to use for construction.
        spec : StateSelectionDescriptor
            Either a single spec string, a Sequence of state ids or state id pairs
            or a sequence of spec strings.
            A selection of ``[(1, 2), (2, 1), (3, 1)]`` means
            to select only transitions between states 1 and 2 as well as from
            3 to 1 (but not from 1 to 3).
            Alternatively, combinations may be specified as a single string
            in the following style: ``'1<>2, 3->1'`` -- this specification
            selects the same combinations as in the previous example, with ``<>``
            selecting transitions in either direction and ``->`` being
            one-directional.
        is_directed : bool, optional
            Flag whether state combinatons should be assumed different, i.e. (1,3) should be considered different from (3,1).
            If not provided, will be set depending on whether there is a directed transition in the descriptors, i.e. `i -> j`.

        Returns
        -------
        Self
            A StateSelection built from the state specification.
        """
        states_coll: set[StateId] = set()
        state_combs_coll: set[StateCombination] = set()

        inputs_directed: bool = False

        if isinstance(spec, str):
            states, combs, has_directed_comb = StateSelection._standard_state_comb_spec(
                spec
            )
            states_coll.update(states)
            state_combs_coll.update(combs)
            inputs_directed |= has_directed_comb
        else:
            for ispec in spec:
                if isinstance(ispec, str):
                    states, combs, has_directed_comb = (
                        StateSelection._standard_state_comb_spec(ispec)
                    )
                    states_coll.update(states)
                    state_combs_coll.update(combs)
                    inputs_directed |= has_directed_comb
                elif isinstance(ispec, tuple):
                    states_coll.update(ispec)
                    if len(ispec) == 2:
                        state_combs_coll.add(ispec)
                else:
                    states_coll.add(int(ispec))

        ground_state_id = np.min(list(states_coll))

        if is_directed is None:
            is_directed = inputs_directed

        return cls(
            is_directed=is_directed,
            states_base=list(states_coll),
            states=list(states_coll),
            ground_state_id=ground_state_id,
            state_types=None,
            state_charges=None,
            state_names=None,
            state_combinations_base=list(state_combs_coll),
            state_combinations=list(state_combs_coll),
            state_combination_names=None,
            state_degeneracy_group=None,
            degeneracy_group_states=None,
        )

    @staticmethod
    def _standard_state_comb_spec(
        spec: str,
    ) -> tuple[list[StateId], list[StateCombination], bool]:
        """Support extracting states and state combinations from strings.

        Parameters
        ----------
        spec : str
            The spec string with a certain pattern.
            A state selection holding states or state transitions that should be used
            in analysis, e.g.:
            A selection of ``[(1, 2), (2, 1), (3, 1)]`` means
            to select only transitions between states 1 and 2 as well as from
            3 to 1 (but not from 1 to 3).
            Alternatively, combinations may be specified as a single string
            in the following style: ``'1<>2, 3->1'`` -- this specification
            selects the same hops as in the previous example, with ``<>``
            selecting hops in either direction and ``->`` being one-
            directional.

        Returns
        -------
        tuple[list[StateId], list[StateCombination], bool]
            First the list of StateIds mentioned in the selection descriptor, then the list of state combinations listed in
            the state combination mentioned in the state selection descriptor.
            Finally a flag whether there was at least one directed state combination specifier.
        """
        if not isinstance(spec, str):
            return spec

        sub_specs = re.split(r"\s*,\s*", spec)

        res_state: set[StateId] = set()
        res_state_comb: set[StateCombination] = set()
        comb_directed: bool = False
        for spec in sub_specs:
            found = _state_combs_pattern.match(spec)
            if found:
                # TODO: FIXME: Deal with S, T, D, S0, and other state names that could be here
                state_from = int(found.group("state_from"))
                state_to = int(found.group("state_to"))
                rel = found.group("rel")
                res_state.add(state_from)
                res_state.add(state_to)
                if rel == "->":
                    res_state_comb.add((state_from, state_to))
                    comb_directed = True
                else:
                    res_state_comb.update(
                        [(state_from, state_to), (state_to, state_from)]
                    )
            else:
                found = _state_id_pattern.match(spec)
                if found:
                    state_id = found.group("state")
                    # TODO: FIXME: Deal with S, T, D, S0, and other state names that could be here
                    state = int(state_id)
                    res_state.add(state)

        return list(res_state), list(res_state_comb), comb_directed

    @staticmethod
    def _abstract_state_comb_spec(
        spec: str, states_from_sc_statements: bool = False
    ) -> tuple[
        set[StateId],
        set[StateCombination],
        set[MultiplicityLabel | str],
        set[
            tuple[
                MultiplicityLabel | str | StateId,
                MultiplicityLabel | str | StateId,
            ]
        ],
        bool,
    ]:
        """Support for abstract state combination description, allowing for inputs like
        `T->S` to specify transitions from a triplet to a singlet state in addition to the extraction
        of specific state ids in `._standard_state_comb_spec()`

        Parameters
        ----------
        spec : str
            The spec string with a certain pattern.
            A state selection (optionally including patterns) holding states or state transitions that should be used
            in analysis, e.g.:
            A selection of ``[(1, 2), (2, 1), (3, 1)]`` means
            to select only transitions between states 1 and 2 as well as from
            3 to 1 (but not from 1 to 3).
            Alternatively, combinations may be specified as a single string
            in the following style: ``'1<>2, 3->1'`` -- this specification
            selects the same hops as in the previous example, with ``<>``
            selecting hops in either direction and ``->`` being one-
            directional.
            States can also be described using their names (if configured, e.g. `'S0'`)
            or a label describing their multiplicities (if configured, e.g. `'T'` for all triplet states)
        states_from_sc_statements : bool, default=False
            Flag whether states mentioned in state combination statements should be included in the
            resulting sets for explicit states and state patterns.
            E.g. if the selection is `1->2`, normally, the state selection would be empty,
            but with `states_from_sc_statements=True`, the state selection would be `[1,2]`.

        Returns
        -------
        tuple[
            set[StateId],
            set[StateCombination],
            set[MultiplicityLabel | str],
            set[
                tuple[
                    MultiplicityLabel | str | StateId,
                    MultiplicityLabel | str | StateId,
                ]
            ],
            bool]
            First the set of explicit `StateIds` mentioned in the selection descriptor,
            then the set of explicit state combinations listed in
            the state combination mentioned in the state selection descriptor.
            Third the set of patterns for states, i.e. the set of multiplicity labels or state names,
            Fourth the set of patterns for state types, i.e. tuples involving at least one pattern for states.
            Finally a flag whether there was at least one directed state combination specifier.
        """
        if not isinstance(spec, str):
            return spec

        sub_specs = re.split(r"\s*,\s*", spec)

        res_state: set[StateId] = set()
        res_state_comb: set[StateCombination] = set()
        res_state_patterns: set[MultiplicityLabel | str] = set()
        res_state_comb_patterns: set[
            tuple[
                MultiplicityLabel | str | StateId,
                MultiplicityLabel | str | StateId,
            ]
        ] = set()

        comb_directed: bool = False
        for spec in sub_specs:
            spec = spec.strip()
            found = _state_combs_pattern.match(spec)
            if found:
                # TODO: FIXME: Deal with S, T, D, S0, and other state names that could be here
                state_from = found.group("state_from")
                state_to = found.group("state_to")
                rel = found.group("rel")
                comb_is_pattern = False

                from_id = _state_id_pattern.match(state_from)
                if from_id:
                    res_from = int(state_from)
                    if states_from_sc_statements:
                        res_state.add(res_from)
                else:
                    comb_is_pattern = True
                    from_mult = _state_mult_pattern.match(state_from)
                    from_name = _state_name_pattern.match(state_from)
                    if from_mult or from_name:
                        res_from = state_from
                    else:
                        continue

                    if states_from_sc_statements:
                        res_state_patterns.add(res_from)
                to_id = _state_id_pattern.match(state_to)
                if to_id:
                    res_to = int(state_to)
                    if states_from_sc_statements:
                        res_state.add(res_to)
                else:
                    comb_is_pattern = True
                    to_mult = _state_mult_pattern.match(state_to)
                    to_name = _state_name_pattern.match(state_to)
                    if to_mult or to_name:
                        res_to = state_to
                    else:
                        continue

                    if states_from_sc_statements:
                        res_state_patterns.add(res_to)
                if rel == "->":
                    if comb_is_pattern:
                        res_state_comb_patterns.add((res_from, res_to))
                    else:
                        res_state_comb.add((res_from, res_to))  # type: ignore
                    comb_directed = True
                else:
                    # Bidirectional matches
                    if comb_is_pattern:
                        res_state_comb_patterns.add((res_from, res_to))
                        res_state_comb_patterns.add((res_to, res_from))
                    else:
                        res_state_comb.add((res_to, res_from))  # type: ignore
                        res_state_comb.add((res_from, res_to))  # type: ignore
            else:
                found = _state_id_pattern.match(spec)
                if found:
                    state_id = found.group("state")
                    state = int(state_id)
                    res_state.add(state)
                else:
                    found_mult = _state_mult_pattern.match(spec)
                    found_name = _state_name_pattern.match(spec)
                    if found_mult or found_name:
                        res_state_patterns.add(spec)

        return (
            res_state,
            res_state_comb,
            res_state_patterns,
            res_state_comb_patterns,
            comb_directed,
        )

    def as_directed_selection(self) -> Self:
        """Helper method to turn an undirected selection into a directed selection.

        If the selection is already directed, it will be returned unmodified.

        Returns
        -------
        Self
            Either the already directed selection or a copy with all mirrored transitions also initially included.
        """
        if self.is_directed:
            return self
        else:
            new_combs_base = list(self.state_combinations_base) + list(
                (y, x) for (x, y) in self.state_combinations_base
            )
            new_combs = list(self.state_combinations) + list(
                (y, x) for (x, y) in self.state_combinations
            )
            if self.state_combination_names is not None:
                new_comb_names = dict(self.state_combination_names)
                new_comb_names.update(
                    {(y, x): v for (x, y), v in self.state_combination_names.items()}
                )
            else:
                new_comb_names = None

            if self.state_combination_colors is not None:
                new_comb_colors = dict(self.state_combination_colors)
                new_comb_colors.update(
                    {(y, x): v for (x, y), v in self.state_combination_colors.items()}
                )
            else:
                new_comb_colors = None

        res = self.copy_or_update(
            is_directed=True,
            state_combinations=new_combs,
            state_combination_colors=new_comb_colors,
            state_combination_names=new_comb_names,
        )
        res.state_combinations_base = new_combs_base
        return res

    def as_undirected_selection(self) -> Self:
        """Helper method to turn a directed selection into an undirected selection.

        If the selection is already undirected, it will be returned unmodified.

        Returns
        -------
        Self
            Either the already undirected selection or a copy with all transitions reduced to those with canonical order.
        """
        if not self.is_directed:
            return self
        else:
            new_combs_base = sorted(
                list(
                    set(
                        StateSelection._state_comb_canonicalized(sc, is_directed=False)
                        for sc in self.state_combinations_base
                    )
                )
            )
            new_combs = sorted(
                list(
                    set(
                        StateSelection._state_comb_canonicalized(sc, is_directed=False)
                        for sc in self.state_combinations
                    )
                )
            )
            # We need to add swapped names and colors if they have been set.
            if self.state_combination_names is not None:
                new_comb_names = dict(self.state_combination_names)
                for sc in new_comb_names:
                    canonic_sc = StateSelection._state_comb_canonicalized(
                        sc, is_directed=False
                    )
                    if canonic_sc not in new_comb_names:
                        new_comb_names[canonic_sc] = new_comb_names[sc]
            else:
                new_comb_names = None

            if self.state_combination_colors is not None:
                new_comb_colors = dict(self.state_combination_colors)
                new_comb_colors.update(
                    {(y, x): v for (x, y), v in self.state_combination_colors.items()}
                )
                new_comb_colors = dict(self.state_combination_colors)
                for sc in new_comb_colors:
                    canonic_sc = StateSelection._state_comb_canonicalized(
                        sc, is_directed=False
                    )
                    if canonic_sc not in new_comb_colors:
                        new_comb_colors[canonic_sc] = new_comb_colors[sc]
            else:
                new_comb_colors = None

        res = self.copy_or_update(
            is_directed=False,
            state_combinations=new_combs,
            state_combination_colors=new_comb_colors,
            state_combination_names=new_comb_names,
        )
        res.state_combinations_base = new_combs_base
        return res

    @staticmethod
    def _state_comb_canonicalized(
        comb: StateCombination, is_directed: bool
    ) -> StateCombination:
        """Helper to turn transitions into a canonic order if the selection is not directed.

        If the selection is directed, the combination is returned as-is.

        Parameters
        ----------
        comb : StateCombination
            The combination tuple to turn into a canonical order
        is_directed : bool
            Flag whether the order should be canonical in a set with this `is_directed` flag

        Returns
        -------
        StateCombination
            The canonicalized combination tuple.
        """
        if is_directed:
            return comb
        else:
            return (min(comb[0], comb[1]), max(comb[0], comb[1]))

    def _state_id_matches_pattern(
        self, state: StateId, pattern: MultiplicityLabel | str | StateId
    ) -> bool:
        """Helper function to check whether a state Id matches a certain string pattern provided by a user.

        Parameters
        ----------
        state : StateId
            The state id to check for a match
        pattern : MultiplicityLabel | str | StateId
            The pattern to compare the state to.
            Can be a multiplicity label or a state name.
            If the values for multiplicity labels or state names are not set, this may result in an exception
            being raised

        Returns
        -------
        bool
            Whether the state matches the pattern

        Raises
        ------
        RuntimeError
            If matching for multiplicity or name is requested and type information or name data is missin.
        """
        if isinstance(pattern, StateId):
            return state == pattern
        elif pattern in MultiplicityLabelValues:
            if self.state_types is None or state not in self.state_types:
                raise RuntimeError(
                    "State multiplicities are not configured on this state selection. Cannot match for multiplicity labels like {pattern}"
                )

            state_mult = self.state_types[state]
            return state_mult == self._mult_label_transl(pattern)
        else:
            if self.state_names is None or state not in self.state_names:
                logging.warning(
                    "Matching against a state name without state names being set. Using default state names."
                )

            return self.get_state_name_or_default(state) == pattern

    def _state_ids_match_pattern(
        self, base_selection: Iterable[StateId], pattern: MultiplicityLabel | str
    ) -> set[StateId]:
        """Helper function to check which states out of a collection matches a certain string pattern provided by a user.

        Parameters
        ----------
        base_selection : Iterable[StateId]s
            The state id to check for a match
        pattern : MultiplicityLabel | str
            The pattern to compare the state to.
            Can be a multiplicity label or a state name.
            If the values for multiplicity labels or state names are not set, this may result in an error.

        Returns
        -------
        set[StateId]
            The set of state ids from the selection that adhere to the pattern.
        """
        return set(
            x for x in base_selection if self._state_id_matches_pattern(x, pattern)
        )

    def _state_combs_matches_pattern(
        self,
        state_comb: StateCombination,
        pattern: tuple[
            MultiplicityLabel | str | StateId, MultiplicityLabel | str | StateId
        ],
    ) -> bool:
        """Helper function to check whether a specific state combinations matches a certain string pattern provided by a user.

        Parameters
        ----------
        state_comb : StateCombination
            The state combination to check for a match
        pattern : tuple[MultiplicityLabel | str | StateId, MultiplicityLabel | str | StateId]
            The pattern to compare the state to.
            Each entry can be a multiplicity label or a state name and of of both entries can be an explicit state id..
            If the values for multiplicity labels or state names are not set, this may result in an error.

        Returns
        -------
        boole
            Boolean flag whether the state combinations matches the pattern
        """
        return self._state_id_matches_pattern(
            state_comb[0], pattern[0]
        ) and self._state_id_matches_pattern(state_comb[1], pattern[1])

    def _state_combs_match_pattern(
        self,
        base_selection: Iterable[StateCombination],
        pattern: tuple[
            MultiplicityLabel | str | StateId, MultiplicityLabel | str | StateId
        ],
    ) -> set[StateCombination]:
        """Helper function to check which state combinations out of a collection matches a certain string pattern provided by a user.

        Parameters
        ----------
        base_selection : Iterable[StateCombination]
            The state combinations to check for a match
        pattern : tuple[MultiplicityLabel | str, MultiplicityLabel | str]
            The pattern to compare the state to.
            Each entry can be a multiplicity label or a state name.
            If the values for multiplicity labels or state names are not set, this may result in an error.

        Returns
        -------
        set[StateCombination]
            The set of state combination identifiers from the selection that adhere to the pattern.
        """
        return set(
            x for x in base_selection if self._state_combs_matches_pattern(x, pattern)
        )

    def select_states(
        self,
        selectors: Iterable[StateId | StateSelectionDescriptor]
        | StateId
        | StateSelectionDescriptor
        | None = None,
        *,
        exclude_ids: Iterable[StateId] | StateId | None = None,
        charge: Iterable[int] | int | None = None,
        exclude_charge: Iterable[int] | int | None = None,
        multiplicity: Iterable[int | MultiplicityLabel]
        | int
        | MultiplicityLabel
        | None = None,
        exclude_multiplicity: Iterable[int | MultiplicityLabel]
        | int
        | MultiplicityLabel
        | None = None,
        min_states_in_selection: Literal[0, 1, 2] = 0,
        states_from_sc: bool = False,
        inplace: bool = False,
    ) -> Self:
        """Method to get a new state selection only retaining the states satisfying the required inclusion criteria and
        not satisfying the exclusion criteria.

        Will return a new StateSelection object with the resulting set of states.

        Parameters
        ----------
        selectors : Iterable[StateId or StateSelectionDescriptor] or StateId or StateSelectionDescriptor or None, optional
            Explicit ids of states to retain. Either a single id or an iterable collection of state ids can be provided. Defaults to None.
        exclude_ids : Iterable[StateId] or StateId or None, optional)
            Explicit ids of states to exclude. Either a single id or an iterable collection of state ids can be provided. Defaults to None.
        charge : Iterable[int] or int or None, optional
            Charges of states to retain. Defaults to None.
        exclude_charge : Iterable[int] or int or None, optional
            Charges of states to exclude. Defaults to None.
        multiplicity : Iterable[int] or int or None, optional
            Multiplicity of states to retain. Defaults to None.
        exclude_multiplicity : Iterable[int] or int or None, optional
            Multiplicity of states to exclude. Defaults to None.
        min_states_in_selection : {0, 1, 2}, optional
            Optional parameter to determine whether state combinations should be kept if states they include are no longer part of the selection.
            A state combination is retained if at least `min_states_in_selection` of their states are still within the state selection. Defaults to 0, meaning all combinations are kept.
        states_from_sc : bool, default=False
            Flag, whether states should be included in the selection based on states mentioned in state combination
            selectors. By default False.
        inplace : bool, optional
            Flag to update the selection in-place. Defaults to False, meaning a modified copy is returned.

        Returns
        -------
        StateSelection
            The resulting selection after applying all of the requested conditions.
        """
        new_states = set(self.states)

        if selectors:
            selector_state_ids: set[StateId] = set()
            # selector_state_combs: set[StateCombination] = set()

            if isinstance(selectors, StateId):
                selector_state_ids.add(selectors)
            # elif (
            #     isinstance(selectors, tuple)
            #     and len(selectors) == 2
            #     and isinstance(selectors[0], StateId)
            #     and isinstance(selectors[1], StateId)
            # ):
            #     selector_state_combs.add(selectors)
            else:
                if isinstance(selectors, str):
                    selectors = [selectors]

                for sel in selectors:
                    if isinstance(sel, StateId):
                        selector_state_ids.add(sel)
                    elif isinstance(sel, str):
                        (
                            expl_state,
                            expl_comb,
                            pattern_state,
                            pattern_combs,
                            is_directed,
                        ) = self._abstract_state_comb_spec(
                            sel, states_from_sc_statements=states_from_sc
                        )

                        selector_state_ids.update(expl_state)

                        for pattern in pattern_state:
                            selector_state_ids.update(
                                self._state_ids_match_pattern(new_states, pattern)
                            )
                    else:
                        # Skip combination descriptors
                        pass

            new_states = selector_state_ids

        if exclude_ids:
            if isinstance(exclude_ids, StateId):
                exclude_ids = [exclude_ids]

            next_states = []
            for old_state in new_states:
                if old_state not in exclude_ids:
                    next_states.append(old_state)

            new_states = next_states

        if charge:
            if isinstance(charge, int):
                charge = [charge]
            if self.state_charges:
                next_states = []

                for old_state in new_states:
                    if (
                        old_state in self.state_charges
                        and self.state_charges[old_state] in charge
                    ):
                        next_states.append(old_state)

                new_states = next_states
            else:
                raise ValueError(
                    "Requested filtering by charges but state charges are unknown. Please set the charges first."
                )

        if exclude_charge:
            if isinstance(exclude_charge, int):
                exclude_charge = [exclude_charge]
            if self.state_charges:
                next_states = []

                for old_state in new_states:
                    if (
                        old_state not in self.state_charges
                        or self.state_charges[old_state] not in exclude_charge
                    ):
                        next_states.append(old_state)

                new_states = next_states
            else:
                raise ValueError(
                    "Requested filtering by charges but state charges are unknown. Please set the charges first."
                )

        if multiplicity:
            if isinstance(multiplicity, int) or isinstance(multiplicity, str):
                multiplicity = [multiplicity]

            trans_mult = self._mult_label_transl(multiplicity)
            if self.state_types:
                next_states = []

                for old_state in new_states:
                    if (
                        old_state in self.state_types
                        and self.state_types[old_state] in trans_mult
                    ):
                        next_states.append(old_state)

                new_states = next_states
            else:
                raise ValueError(
                    "Requested filtering by multiplicities but state multiplicities are unknown. Please set the multiplicities/types first."
                )

        if exclude_multiplicity:
            if isinstance(exclude_multiplicity, int) or isinstance(
                exclude_multiplicity, str
            ):
                exclude_multiplicity = [exclude_multiplicity]

            trans_mult_excl = self._mult_label_transl(exclude_multiplicity)

            if self.state_types:
                next_states = []

                for old_state in new_states:
                    if (
                        old_state in self.state_types
                        and self.state_types[old_state] not in trans_mult_excl
                    ):
                        next_states.append(old_state)

                new_states = next_states
            else:
                raise ValueError(
                    "Requested filtering by multiplicities but state multiplicities are unknown. Please set the multiplicities/types first."
                )

        updated_selection = self.copy_or_update(
            states=list(new_states), inplace=inplace
        ).select_state_combinations(
            min_states_in_selection=min_states_in_selection, inplace=inplace
        )

        return updated_selection

    @overload
    @staticmethod
    def _mult_label_transl(
        multipl: Iterable[int | MultiplicityLabel],
    ) -> set[int]: ...

    @overload
    @staticmethod
    def _mult_label_transl(
        multipl: int | MultiplicityLabel,
    ) -> int: ...

    @staticmethod
    def _mult_label_transl(
        multipl: int | MultiplicityLabel | Iterable[int | MultiplicityLabel],
    ) -> int | set[int]:
        """Function to translate potential string-based multiplicities to integers

        Parameters
        ----------
        multipl : int or MultiplicityLabel or Iterable[int or MultiplicityLabel]
            List of multiplicities, either ints or string labels

        Returns
        -------
        int or set[int]
            A set representation of the numeric multiplicities or the single translated value
        """
        if isinstance(multipl, int):
            return multipl
        elif isinstance(multipl, str) and multipl in MultiplicityLabelValues:
            assert isinstance(multipl, str)

            lower_label = multipl.lower()
            if lower_label.startswith("s"):
                return 1
            elif lower_label.startswith("d"):
                return 2
            elif lower_label.startswith("t"):
                return 3
            else:
                raise ValueError(
                    f"Label `{multipl}` is not a valid multiplicity label."
                )
        else:
            return set(StateSelection._mult_label_transl(x) for x in multipl)  # type: ignore # If the input is appropriate, this should yield appropriate results.

    def select_state_combinations(
        self,
        selectors: StateSelectionDescriptor | None = None,
        *,
        ids: Iterable[StateCombination] | None = None,
        min_states_in_selection: Literal[0, 1, 2] = 0,
        inplace: bool = True,
    ) -> Self:
        """Method to get a new state selection with a potentially reduced set of state combinations.

        Parameters
        ----------
        selectors : StateSelectionDescriptor or Iterable[StateSelectionDescriptor], optional
            A textual or tuple-based description of the
        ids : Iterable[StateCombination] or None, optional
            Explicit state transitions ids to retain. Defaults to None.
        min_states_in_selection : Literal[0, 1, 2], optional
            Minimum number of states involved in the state combination that still need to be within the state selection to keep this combination. Defaults to 0, meaning no check will be performed.
        inplace : bool, optional
            Flag to update the selection in-place. Defaults to False, meaning a modified copy is returned.


        Returns
        -------
        StateSelection
            A new state selection with potentially fewer state combinations considered.
        """
        new_state_combinations = set(self.state_combinations)

        if selectors:
            # selector_state_ids: set[StateId] = set()
            selector_state_combs: set[StateCombination] = set()

            if isinstance(selectors, StateId):
                pass
                # selector_state_combs.add(selectors)

            elif (
                isinstance(selectors, tuple)
                and len(selectors) == 2
                and isinstance(selectors[0], StateId)
                and isinstance(selectors[1], StateId)
            ):
                selector_state_combs.add(selectors)  # type: ignore
            else:
                if isinstance(selectors, str):
                    selectors = [selectors]

                for sel in selectors:
                    if isinstance(selectors, StateId):
                        continue
                        # selector_state_combs.add(selectors)

                    elif (
                        isinstance(selectors, tuple)
                        and len(selectors) == 2
                        and isinstance(selectors[0], StateId)
                        and isinstance(selectors[1], StateId)
                    ):
                        selector_state_combs.add(selectors)  # type: ignore
                    elif isinstance(sel, str):
                        (
                            expl_state,
                            expl_comb,
                            pattern_state,
                            pattern_combs,
                            is_directed,
                        ) = self._abstract_state_comb_spec(sel)

                        if is_directed and not self.is_directed:
                            # If the selection requires directed states, convert to directed state representations.
                            return (
                                self.as_directed_selection().select_state_combinations(
                                    selectors=selectors,
                                    ids=ids,
                                    min_states_in_selection=min_states_in_selection,
                                    inplace=inplace,
                                )
                            )

                        selector_state_combs.update(expl_comb)

                        for pattern in pattern_combs:
                            selector_state_combs.update(
                                self._state_combs_match_pattern(
                                    new_state_combinations, pattern
                                )
                            )
                    else:
                        # Skip combination descriptors
                        pass

            new_state_combinations = set(
                self._state_comb_canonicalized(x, self.is_directed)
                for x in selector_state_combs
            )

        if ids:
            # Standardize selection tuple order
            ids = list(self._state_comb_canonicalized(x, self.is_directed) for x in ids)
            # # Filter explicit states
            # id_state_combinations = [
            #     comb for comb in new_state_combinations if comb in ids
            # ]
            new_state_combinations.update(ids)

        if min_states_in_selection > 0:
            # Check that there are sufficiently many states of the combination still in teh selection
            retained_combs = []
            for comb in new_state_combinations:
                states = set(comb)
                num_selected_states = len(states.intersection(self.states))

                if num_selected_states >= min_states_in_selection:
                    retained_combs.append(comb)

            new_state_combinations = retained_combs

        if selectors is None and ids is None and min_states_in_selection == 0:
            new_state_combinations = list(self.state_combinations_base)

        return self.copy_or_update(
            state_combinations=list(new_state_combinations), inplace=inplace
        )

    def select(
        self,
        selectors: Iterable[StateId | StateCombination | StateSelectionDescriptor]
        | StateId
        | StateCombination
        | StateSelectionDescriptor
        | None = None,
        *,
        exclude: Iterable[StateId | StateSelectionDescriptor]
        | StateId
        | StateSelectionDescriptor
        | None = None,
        min_states_in_selection: Literal[0, 1, 2] = 0,
        states_from_sc: bool = False,
        # inplace: bool = False,
    ) -> Self:
        """Method to select both states and state combinations in one go.

        Internally calls `.select_states()` and `.select_state_combinations()`.
        Additionally, the `exclude` keyword parameter can be used to remove matched entries from the
        overall selection.
        If no parameters are provided, all states and state combinations will be provided.
        May implicitly convert to a directed selection if descriptors include directed transitions.

        Parameters
        ----------
        selectors : Iterable[StateId  |  StateCombination  |  StateSelectionDescriptor] | StateId | StateCombination | StateSelectionDescriptor | None, optional
            The description of states and state combinations supposed to be included within the resulting selection.
            The description can be a state id (`int`), a transition id, (`tuple[int,int]`) or a `str` with comma-separated statements
            denoting either states or state combinations.
            The following statements are supported in selector strings:
                - A `str` representation of a state id.
                - A multiplicity label e.g. `'S'` or `'t'` (capitalization irrelevant, requires state types to be set)
                - A state name, if state names are configured on this selection, e.g, `'S0'`
                - A state combination representation as a `str`, of one of the following forms:
                    - `'<state_a> -> <state_b>'` (for directed transitions a to b)
                    - `'<state_a> <> <state_b>'` (for undirected/bidirectional transitions between a and b)
                  All `<state>` expressions may take any statement of the state representations above.
            By default, no states will be selected.
        exclude : Iterable[StateId  |  StateSelectionDescriptor] | StateId | StateSelectionDescriptor | None, optional
            Same format at `selectors`. If set, the states and state combinations selected by this argument are excluded from the resulting selection.
            By default None, meaning no states will be removed.
        min_states_in_selection : Literal[0, 1, 2], optional
            Optional parameter to denote, how many states of a state combination must be within the selection for the state combination to be included in a result, by default 0
            E.g. if only states `1` and `3` are selected, and the selection would include `1->2`, this parameter
            needs to be `1` or `0` for `(1,2)` to be included in the result.
        states_from_sc : bool, default=False
            Flag, whether states should be included in the selection based on states mentioned in state combination
            selectors. By default False.

        Returns
        -------
        Self
            The resulting selection with states and state combinations selected.
        """
        tmp_res_exclude = None
        if exclude:
            # We create a copy first and then we can modify in place.
            tmp_res_exclude = self.select_states(
                exclude, inplace=False, states_from_sc=states_from_sc
            ).select_state_combinations(selectors=exclude, inplace=True)

        # We create a copy first and then we can modify in place.
        tmp_res_select = self.select_states(
            selectors, inplace=False, states_from_sc=states_from_sc
        ).select_state_combinations(
            selectors=selectors,
            min_states_in_selection=min_states_in_selection,
            inplace=True,
        )

        if tmp_res_exclude:
            return tmp_res_select - tmp_res_exclude
        else:
            return tmp_res_select

    def set_state_names(
        self, names: Sequence[str] | dict[StateId, str], inplace: bool = True
    ) -> Self:
        """Helper function to assign new state names to the selection.

        Will peform some sanity checks first.

        Parameters
        ----------
        names : Sequence[str] or dict[StateId, str]
            Either a list of state names aligned with `self.states` ids or a dictionary mapping state ids to names.
        inplace : bool, optional
            Flag to determine whether this function should update the existing selection sequence or return a modified copy. Defaults to True, meaning the existing instance is updated.


        Returns
        -------
        Self
            Either the existing selection with updated names or a new instance with modified names.

        Raises
        ------
        ValueError
            If a Sequence is provided that does not have enough values
        ValueError
            If a dict is  provided that does not have mapping for all state ids in `self.states`

        """
        new_state_names = None
        if isinstance(names, dict):
            state_set = set(self.states)
            if state_set.issubset(names.keys()):
                new_state_names = names
            else:
                raise ValueError(
                    f"Provided `names` dict does not have names assigned for all states. It is missing {state_set.difference(names.keys())}."
                )
        else:
            num_names = len(names)
            if num_names >= len(self.states):
                new_state_names = {
                    state_id: state_name
                    for (state_id, state_name) in zip(self.states, names)
                }
            else:
                raise ValueError(
                    f"Provided `names` sequence does not have enough names for the states in this selection. Provided: {num_names}, Required: {len(self.states)}."
                )
        return self.copy_or_update(state_names=new_state_names, inplace=inplace)

    def set_state_types(
        self, types: Sequence[int] | dict[StateId, int], inplace: bool = True
    ) -> Self:
        """Helper function to assign new state types/multiplicites to the selection.

        Will peform some sanity checks first.

        Parameters
        ----------
        types : Sequence[int] or dict[StateId, int]
            Either a list of state types/multiplicities aligned with `self.states` ids or a dictionary mapping state ids to types.
        inplace : bool, optional
            Flag to determine whether this function should update the existing selection sequence or return a modified copy. Defaults to True, meaning the existing instance is updated.


        Returns
        -------
        Self
            Either the existing selection with updated types or a new instance with modified types.

        Raises
        ------
        ValueError
            If a Sequence is provided that does not have enough values
        ValueError
            If a dict is  provided that does not have mapping for all state ids in `self.states`
        """
        new_state_types = None
        if isinstance(types, dict):
            state_set = set(self.states)
            if state_set.issubset(types.keys()):
                new_state_types = types
            else:
                raise ValueError(
                    f"Provided `types` dict does not have names assigned for all states. It is missing {state_set.difference(types.keys())}."
                )
        else:
            num_values = len(types)
            if num_values >= len(self.states):
                new_state_types = {
                    state_id: state_type
                    for (state_id, state_type) in zip(self.states, types)
                }
            else:
                raise ValueError(
                    f"Provided `types` sequence does not have enough types for the states in this selection. Provided: {num_values}, Required: {len(self.states)}."
                )
        return self.copy_or_update(state_types=new_state_types, inplace=inplace)

    def set_state_charges(
        self, charges: int | Sequence[int] | dict[StateId, int], inplace: bool = True
    ) -> Self:
        """Helper function to assign new state charges to the selection.

        Will peform some sanity checks first.

        Parameters
        ----------
        charges : int or Sequence[int] or dict[StateId, int]
            Either a single charge for all states or a list of state charges aligned with `self.states` ids or a dictionary mapping state ids to charges.
        inplace : bool, optional
            Flag to determine whether this function should update the existing selection sequence or return a modified copy. Defaults to True, meaning the existing instance is updated.


        Returns
        -------
        Self
            Either the existing selection with updated charges or a new instance with modified charges.

        Raises
        ------
        ValueError
            If a Sequence is provided that does not have enough charges
        ValueError
            If a dict is provided that does not have mapping for all state ids in `self.states`


        """
        new_state_charges = None
        if isinstance(charges, int):
            new_state_charges = {state_id: charges for state_id in self.states}
        elif isinstance(charges, dict):
            state_set = set(self.states)
            if state_set.issubset(charges.keys()):
                new_state_charges = charges
            else:
                raise ValueError(
                    f"Provided `charges` dict does not have names assigned for all states. It is missing {state_set.difference(charges.keys())}."
                )
        else:
            num_values = len(charges)
            if num_values >= len(self.states):
                new_state_charges = {
                    state_id: state_type
                    for (state_id, state_type) in zip(self.states, charges)
                }
            else:
                raise ValueError(
                    f"Provided `charges` sequence does not have enough charges for the states in this selection. Provided: {num_values}, Required: {len(self.states)}."
                )
        return self.copy_or_update(state_charges=new_state_charges, inplace=inplace)

    def set_state_combinations(
        self, combinations: Sequence[StateCombination], inplace: bool = True
    ) -> Self:
        """Helper function to assign new state combinations to the selection.

        Will peform some sanity checks first.

        Parameters
        ----------
        combinations : Sequence[StateCombination]
            A list of state combination tuples to set to the selection
        inplace : bool, optional
            Flag to determine whether this function should update the existing selection or return a modified copy. Defaults to True, meaning the existing instance is updated.


        Returns
        -------
        Self
            Either the existing selection with updated combinations or a new instance with modified combinations.

        Raises
        ------
        ValueError
            If an entry in the combinations sequence has a non-positive state entry.

        """
        new_state_combinations = None

        for first, second in combinations:
            if first <= 0:
                raise ValueError(f"State {first} from combinations must be positive")
            if second <= 0:
                raise ValueError(f"State {second} from combinations must be positive")

        return self.copy_or_update(
            state_combinations=new_state_combinations, inplace=inplace
        )

    def set_state_combination_names(
        self, names: Sequence[str] | dict[StateCombination, str], inplace: bool = True
    ) -> Self:
        """Helper function to assign new state combination labels to the selection.

        Will peform some sanity checks first.

        Parameters
        ----------
        names : Sequence[str] | dict[StateCombination
            Either a list of state combination names aligned with `self.state_combinations` or a dictionary mapping state combination ids to names.
        inplace : bool, optional
            Flag to determine whether this function should update the existing selection or return a modified copy. Defaults to True, meaning the existing instance is updated.


        Returns
        -------
        Self
            Either the existing selection with updated names or a new instance with modified names.

        Raises
        ------
        ValueError
            If a Sequence is provided that does not have enough values
        ValueError
            If a dict is  provided that does not have mapping for all state combination ids in `self.state_combinations`

        """
        new_state_combination_names = None
        if isinstance(names, dict):
            state_combinations_set = set(self.state_combinations)
            if state_combinations_set.issubset(names.keys()):
                new_state_combination_names = names
            else:
                raise ValueError(
                    f"Provided `names` dict does not have names assigned for all state combinations. It is missing {state_combinations_set.difference(names.keys())}."
                )
        else:
            num_names = len(names)
            if num_names >= len(self.state_combinations):
                new_state_combination_names = {
                    state_comb_id: state_comb_name
                    for (state_comb_id, state_comb_name) in zip(
                        self.state_combinations, names
                    )
                }
            else:
                raise ValueError(
                    f"Provided `names` sequence does not have enough names for the state combinations in this selection. Provided: {num_names}, Required: {len(self.state_combinations)}."
                )
        return self.copy_or_update(
            state_combination_names=new_state_combination_names, inplace=inplace
        )

    def singlets_only(self, inplace: bool = False) -> Self:
        """Helper function to immediately filter only singlet states. Does not affect state combinations.

        Parameters
        ----------
        inplace : bool, optional
            Flag whether the operation should update the selection in-place. Defaults to False.

        Returns
        -------
        StateSelection
            the updated selection only containing singlet states.

        """
        return self.select_states("S", inplace=inplace)

    def triplets_only(self, inplace: bool = False) -> Self:
        """Helper function to immediately filter only triplet states. Does not affect state combinations.

        Parameters
        ----------
        inplace : bool, optional
            Flag whether the operation should update the selection in-place. Defaults to False.

        Returns
        -------
        StateSelection
            the updated selection only containing triplet states.

        """
        return self.select_states("T", inplace=inplace)

    def same_multiplicity_transitions(self, inplace: bool = False) -> Self:
        """Helper function to only retain combinations between states of the same multiplicities (e.g. for NACs)

        Parameters
        ----------
        inplace : bool, optional
            Flag whether the operation should update the selection in-place. Defaults to False.

        Returns
        -------
        StateSelection
            the updated selection only containing transitions between states of same multiplicity (i.e. singlet-singlet, triplet-tiplet).

        """

        if not self.state_types:
            raise ValueError(
                "Cannot filter transitions by state multiplicity without multiplicities/types being set."
            )

        new_state_combs = []
        for comb in self.state_combinations:
            first, second = comb

            if first in self.state_types and second in self.state_types:
                if self.state_types[first] == self.state_types[second]:
                    new_state_combs.append((first, second))

        return self.copy_or_update(state_combinations=new_state_combs, inplace=inplace)

    def different_multiplicity_transitions(self, inplace: bool = False) -> Self:
        """Helper function to only retain combinations between states of the different multiplicities (e.g. for SOCs)

        Parameters
        ----------
        inplace : bool, optional
            Flag whether the operation should update the selection in-place. Defaults to False.

        Returns
        -------
        StateSelection
            the updated selection only containing transitions between states of different multiplicity (i.e. singlet-triplet).

        """

        if not self.state_types:
            raise ValueError(
                "Cannot filter transitions by state multiplicity without multiplicities/types being set."
            )

        new_state_combs = []
        for comb in self.state_combinations:
            first, second = comb

            if first in self.state_types and second in self.state_types:
                if self.state_types[first] != self.state_types[second]:
                    new_state_combs.append((first, second))

        return self.copy_or_update(state_combinations=new_state_combs, inplace=inplace)

    def ground_state_transitions(
        self, ground_state_id: StateId | None = None, inplace: bool = False
    ) -> Self:
        """Helper function to only retain combinations between states containing the lowest-level state id.

        Parameters
        ----------
        ground_state_id : StateId, optional
            Id of the state to be considered the ground state. Defaults to the lowest id of the selected states.
        inplace : bool, optional
            Flag whether the operation should update the selection in-place. Defaults to False.

        Returns
        -------
        StateSelection
            the updated selection only containing transitions between ground state and other states.

        """

        if ground_state_id is None:
            ground_state_id = self.ground_state_id

        new_state_combs = []
        for comb in self.state_combinations:
            first, second = comb

            if first == ground_state_id or second == ground_state_id:
                new_state_combs.append(comb)

        return self.copy_or_update(state_combinations=new_state_combs, inplace=inplace)

    def excited_state_transitions(
        self, ground_state_id: StateId | None = None, inplace: bool = False
    ) -> Self:
        """Helper function to only retain combinations between states not involving the ground state.

        Parameters
        ----------
        ground_state_id : StateId, optional
            Id of the state to be considered the ground state. Defaults to the lowest id of the selected states.
        inplace : bool, optional
            Flag whether the operation should update the selection in-place. Defaults to False.

        Returns
        -------
        StateSelection
            the updated selection only containing transitions between non-ground states.

        """

        if ground_state_id is None:
            ground_state_id = self.ground_state_id

        new_state_combs = []
        for comb in self.state_combinations:
            first, second = comb

            if first != ground_state_id and second != ground_state_id:
                new_state_combs.append(comb)

        return self.copy_or_update(state_combinations=new_state_combs, inplace=inplace)

    def non_degenerate(self, inplace: bool = False) -> Self:
        """Helper function to remove all degenerate states and combinations identical except for degeneracy.

        Parameters
        ----------
        inplace : bool, optional
            Flag whether the operation should update the selection in-place. Defaults to False.

        Returns
        -------
        StateSelection
            the updated selection only containing non-degenerate states and non-degenerate-equivalent combinations.

        """
        if self.state_degeneracy_group is None:
            # If we do not have degeneracy data, return self, no change needed.
            logging.warning("Skipping without degeneracy data")
            return self

        new_states = []
        new_state_combinations = []
        new_state_comb_degeneracy_groups_encountered = set()
        state_group_representative = {}

        for sc in self.state_combinations:
            deg = self.get_state_combination_degeneracy(sc)
            if deg in new_state_comb_degeneracy_groups_encountered or deg[0] == deg[1]:
                # Skip known combinations and internal transitions.
                continue
            new_state_comb_degeneracy_groups_encountered.add(deg)

            s1, s2 = sc
            deg1, deg2 = deg

            if deg1 in state_group_representative:
                s1 = state_group_representative[deg1]
            else:
                state_group_representative[deg1] = s1
                new_states.append(s1)

            if deg2 in state_group_representative:
                s2 = state_group_representative[deg2]
            else:
                state_group_representative[deg2] = s2
                new_states.append(s2)

            # TODO: Check order
            deg_sc = (min(s1, s2), max(s1, s2))
            new_state_combinations.append(deg_sc)

        for state in self.states:
            deg = self.get_state_degeneracy(state)

            if deg in state_group_representative:
                continue
            else:
                state_group_representative[deg] = state
                new_states.append(state)

        return self.copy_or_update(
            states=new_states,
            state_combinations=new_state_combinations,
            inplace=inplace,
        )

    def state_info(self) -> Iterable[StateInfo]:
        """Get an iterator over the states in this selection.


        Yields
        -------
        Iterable[StateInfo]
            An iterator over the available state info

        """

        for id in self.states:
            if self.state_charges:
                charge = self.state_charges[id] if id in self.state_charges else None
            else:
                charge = None

            name = self.get_state_name_or_default(id)

            if self.state_types:
                multiplicity = self.state_types[id] if id in self.state_types else None
            else:
                multiplicity = None

            yield StateInfo(id, name, multiplicity, charge)

    def get_state_name_or_default(self, id: StateId) -> str:
        """Helper method to either get registered state name or a default string to identify the state

        Parameters
        ----------
        id : StateId
            Id of the state to get the name for.


        Returns
        -------
        str
            Label of the state

        """
        if self.state_names:
            if id in self.state_names:
                return self.state_names[id]

        return f"state{id - 1}"

    def get_state_combination_name_or_default(self, comb: StateCombination) -> str:
        """Helper method to either get registered state combination name or a default string to identify the state combination.

        Parameters
        ----------
        comb : StateCombination
            Id of the state combination to get the name for.


        Returns
        -------
        str
            Label of the state combination

        """
        if self.state_combination_names:
            if comb in self.state_combination_names:
                return self.state_combination_names[comb]

        first, second = comb

        s1 = self.get_state_name_or_default(first)
        s2 = self.get_state_name_or_default(second)
        return f"{s1} - {s2}"

    def get_state_tex_label(self, id: StateId) -> str:
        """Function to get a nice tex-printable label with super- and subscripts for the denoted state.

        Parameters
        ----------
        id : StateId
            Id of the state to get the label for


        Returns
        -------
        str
            Tex-label that needs to be enclosed in a math environment to not cause issues.

        """

        statename = self.get_state_name_or_default(id)
        return state_name_to_tex_label(statename)

    def get_state_combination_tex_label(self, comb: StateCombination) -> str:
        """Function to get a nice tex-printable label with super- and subscripts for a state combination in this selection

        Parameters
        ----------
        comb : StateCombination
            Combination identifier to get the label for

        Returns
        -------
        str
            Tex-label that needs to be enclosed in a math environment to not cause issues.

        """
        first, second = comb

        s1 = self.get_state_tex_label(first)
        s2 = self.get_state_tex_label(second)
        return f"{s1} - {s2}"

    def combination_info(
        self, degeneracy_free: bool = False
    ) -> Iterable[StateCombInfo]:
        """Get an iterator over the state combinations in this selection.

        Parameters
        ----------
        degeneracy_free: bool, optional
            If set to true, combinations with already covered degeneracy-groups will be skipped

        Yields
        -------
        Iterable[StateCombInfo]
            An iterator over the available state combination info

        """
        degen_covered = []

        for comb in self.state_combinations:
            if degeneracy_free:
                degen_class = self.get_state_combination_degeneracy(comb)
                if degen_class in degen_covered:
                    continue
                else:
                    degen_covered.append(degen_class)

            name = self.get_state_combination_name_or_default(comb)
            yield StateCombInfo(comb, name)

    def has_state(self, id: StateId) -> bool:
        """Function to check whether a state is in the selection

        Parameters
        ----------
        id : StateId
            The state id to check whether it has been selected


        Returns
        -------
        bool
            True if in the selection, False otherwise.

        """
        return id in self.states

    def has_state_combination(self, comb: StateCombination) -> bool:
        """Function to check whether a state combination is in the selection

        Parameters
        ----------
        comb : StateCombination
            The combination to check whether it has been selected


        Returns
        -------
        bool
            True if in the selection, False otherwise.

        """
        return comb in self.state_combinations

    def auto_assign_colors(self, inplace: bool = True) -> Self:
        """Function to automatically generate colors for states and state combinations

        Parameters
        ----------
        inplace : bool, optional
            Flag whether the operation should update the selection in-place. Defaults to True because setting colors is not a big issue.


        Returns
        -------
        Self
            Returns the updated instance.

        """
        from shnitsel.vis.colormaps import (
            get_default_state_colormap,
            get_default_interstate_colormap_inter_mult,
            get_default_interstate_colormap_same_mult,
        )

        multiplicities = self.state_types

        if multiplicities is None:
            # Consider all states singlets then
            multiplicities = {s: 1 for s in self.states}

        mult_state_collection: dict[int, set[StateId]] = {}

        for state, mult in multiplicities.items():
            if mult not in mult_state_collection:
                mult_state_collection[mult] = set()
            mult_state_collection[mult].add(state)

        full_state_colormap: dict[StateId, str] = {}
        mult_color_maps: dict[int, dict[StateId, str]] = {}
        for mult, states in mult_state_collection.items():
            state_list = list(states)
            state_list.sort()

            if self.state_degeneracy_group is not None:
                state_deg_group = [self.state_degeneracy_group[s] for s in state_list]
            else:
                state_deg_group = None

            mult_color_maps[mult] = {
                state_id: color
                for state_id, color in zip(
                    state_list,
                    get_default_state_colormap(
                        len(state_list),
                        multiplicity=mult,
                        degeneracy_groups=state_deg_group,
                    ),
                )
            }
            full_state_colormap.update(mult_color_maps[mult])

        full_interstate_colormap: dict[StateCombination, str] = {}

        for mult1, state_colors1 in mult_color_maps.items():
            for mult2, state_colors2 in mult_color_maps.items():
                if mult1 == mult2:
                    full_interstate_colormap.update(
                        get_default_interstate_colormap_same_mult(mult1, state_colors1)
                    )
                else:
                    full_interstate_colormap.update(
                        get_default_interstate_colormap_inter_mult(
                            state_colors1, state_colors2
                        )
                    )
        return self.copy_or_update(
            state_colors=full_state_colormap,
            state_combination_colors=full_interstate_colormap,
            inplace=inplace,
        )

    def get_state_color(self, id: StateId) -> str:
        """Function to get a the state color or a default color value

        Parameters
        ----------
        id : StateId
            Id of the state to get the color for


        Returns
        -------
        str
            Hex-str color code

        """
        from shnitsel.vis.colormaps import st_grey

        if self.state_colors is not None and id in self.state_colors:
            return self.state_colors[id]
        else:
            return st_grey

    def get_state_combination_color(self, comb: StateCombination) -> str:
        """Function to get a the state combination color or a default color value

        Parameters
        ----------
        comb : StateCombination
            Id of the state combination to get the color for

        Returns
        -------
        str
            Hex-str color code

        """
        from shnitsel.vis.colormaps import st_grey

        if (
            self.state_combination_colors is not None
            and comb in self.state_combination_colors
        ):
            return self.state_combination_colors[comb]
        else:
            return st_grey

    def get_state_combination_degeneracy(
        self, comb: StateCombination
    ) -> tuple[int, int]:
        """Function to get the combined degeneracy classes of the two states.

        Helpful for not plotting too degenerate entries.

        Parameters
        ----------
        comb : StateCombination
            Id of the state combination to get the color for

        Returns
        -------
        tuple[int, int]
            Degeneracy groups of either state

        """

        return self.get_state_degeneracy(comb[0]), self.get_state_degeneracy(comb[1])

    def get_state_degeneracy(self, state: StateId) -> int:
        """Function to get the combined degeneracy classes of the two states.

        Helpful for not plotting too degenerate entries.

        Parameters
        ----------
        comb : StateCombination
            Id of the state combination to get the color for


        Returns
        -------
        str
            Hex-str color code

        """

        if (
            self.state_degeneracy_group is not None
            and state in self.state_degeneracy_group
        ):
            return self.state_degeneracy_group[state]

        return state

    # TODO: FIXME: Add print output __str__, __html__ and __repr__

    def __add__(self, other: Self | StateSelectionDescriptor) -> Self:
        """Add the states and state combinations of another state selection into the selection
        represented by this selection.

        For consistency reasons, the other StateSelection (which can be provided as a descriptor instead),
        should be built upon the same base state and state combination ground set.

        Parameters
        ----------
        other : Self | StateSelectionDescriptor
            The states and state combinations to add to this selection, either as another state selection or
            as a description of states and state combinations that can be passed to
            `StateSelection.init_from_descriptor()`.

        Returns
        -------
        Self
            A `StateSelection` object representing the union of the selections
        """
        other_selection: StateSelection
        if not isinstance(other, StateSelection):
            other_selection = StateSelection.init_from_descriptor(other)
        else:
            other_selection = other

        # TODO: FIXME: Add some checks that the ground information of the two selections is the same

        new_state_selected = list(set(self.states).union(other_selection.states))
        new_sc_selected = list(
            set(self.state_combinations).union(other_selection.state_combinations)
        )

        return self.copy_or_update(
            states=new_state_selected, state_combinations=new_sc_selected
        )

    # a|b is meant as an alias for a+b in the set-theory sense.
    __or__ = __add__

    def __sub__(self, other: Self | StateSelectionDescriptor) -> Self:
        """Remove the states and state combinations of another state selection from the selection
        represented by this selection.

        For consistency reasons, the other StateSelection (which can be provided as a descriptor instead),
        should be built upon the same base state and state combination ground set.

        Parameters
        ----------
        other : Self | StateSelectionDescriptor
            The states and state combinations to remove to this selection, either as another state selection or
            as a description of states and state combinations that can be passed to
            `StateSelection.init_from_descriptor()`.

        Returns
        -------
        Self
            A `StateSelection` object representing the difference of the selections
        """
        other_selection: StateSelection
        if not isinstance(other, StateSelection):
            other_selection = StateSelection.init_from_descriptor(other)
        else:
            other_selection = other

        # TODO: FIXME: Add some checks that the ground information of the two selections is the same

        new_state_selected = list(set(self.states).difference(other_selection.states))
        new_sc_selected = list(
            set(self.state_combinations).difference(other_selection.state_combinations)
        )

        return self.copy_or_update(
            states=new_state_selected, state_combinations=new_sc_selected
        )

    def __and__(self, other: Self | StateSelectionDescriptor) -> Self:
        """Get a selection of the states and state combinations shared between
        this and another state selection.

        For consistency reasons, the other StateSelection (which can be provided
        as a descriptor instead), should be built upon the same base state and
        state combination ground set.

        Parameters
        ----------
        other : Self | StateSelectionDescriptor
            The states and state combinations to intersect with this
            this selection, either as another state selection or
            as a description of states and state combinations that
            can be passed to `StateSelection.init_from_descriptor()`.

        Returns
        -------
        Self
            A `StateSelection` object representing the intersection of the selections
        """
        other_selection: StateSelection
        if not isinstance(other, StateSelection):
            other_selection = StateSelection.init_from_descriptor(other)
        else:
            other_selection = other

        # TODO: FIXME: Add some checks that the ground information of the two selections is the same

        new_state_selected = list(set(self.states).intersection(other_selection.states))
        new_sc_selected = list(
            set(self.state_combinations).intersection(
                other_selection.state_combinations
            )
        )

        return self.copy_or_update(
            states=new_state_selected, state_combinations=new_sc_selected
        )

    def __invert__(self) -> Self:
        """Get an inverted selection of the states and state combinations in
        this state selection.

        Warning
        -------

        The result of this operation will only be as expected, if the selection was
        built from a full dataset such that `states_base` and `state_combinations_base`
        have been set correctly.
        If the selection was built from a textual or tuple representation, the
        inverted selection will only consider the states listed in the original
        description.

        Returns
        -------
        Self
            A `StateSelection` object representing the inverted selection in this object
        """

        # TODO: FIXME: Add some checks that the ground information of the two selections is the same

        new_state_selected = list(set(self.states_base).difference(self.states))
        new_sc_selected = list(
            set(self.state_combinations_base).intersection(self.state_combinations)
        )

        return self.copy_or_update(
            states=new_state_selected, state_combinations=new_sc_selected
        )

    union = __add__
    intersect = __and__
    difference = __sub__
    invert = __invert__
