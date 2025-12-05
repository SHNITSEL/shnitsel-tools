from dataclasses import dataclass
from itertools import combinations
import logging
from typing import Iterable, Self, Sequence, Literal

import numpy as np
import xarray as xr

from shnitsel.data.state_helpers import state_name_to_tex_label
from ..core.typedefs import (
    StateCombination,
    StateId,
    StateInfo,
    StateCombInfo,
    MultiplicityLabel,
)


@dataclass
class StateSelection:
    """Class to keep track of a (sub-)selection of states and state transitions for analysis and plotting."""

    states: Sequence[StateId]
    state_types: dict[StateId, int] | None
    state_names: dict[StateId, str] | None
    state_charges: dict[StateId, int] | None

    state_combinations: list[StateCombination]
    state_combination_names: dict[StateCombination, str] | None

    state_colors: dict[StateId, str] | None = None
    state_combination_colors: dict[StateCombination, str] | None = None

    def copy_or_update(
        self,
        states: Sequence[StateId] | None = None,
        state_types: dict[StateId, int] | None = None,
        state_names: dict[StateId, str] | None = None,
        state_charges: dict[StateId, int] | None = None,
        state_combinations: list[StateCombination] | None = None,
        state_combination_names: dict[StateCombination, str] | None = None,
        state_colors: dict[StateId, str] | None = None,
        state_combination_colors: dict[StateCombination, str] | None = None,
        inplace: bool = False,
    ) -> Self:
        """Function to create a copy with replaced member values.

        Meant as a helper for the `Frozen` logic of the selection, i.e. method calls return a new instance
        instead of updating the existing instance.

        Args:
            states (Sequence[StateId] | None, optional): Potentially new state ids. Defaults to None.
            state_types (dict[StateId, int] | None, optional): Potentially new state types/multiplicities. Defaults to None.
            state_names (dict[StateId, str] | None, optional): Potentially new state names. Defaults to None.
            state_charges (dict[StateId, int] | None, optional): Potentially new state charges. Defaults to None.
            state_combinations (list[StateCombination] | None, optional): Potentially new state combinations. Defaults to None.
            state_combination_names (dict[StateCombination, str] | None, optional): Potentially new names for state combinations. Defaults to None.
            inplace (bool, optional): A flag whether the existing instance should be updated or a new one should be created. Defaults to False, i.e. a new instance is created.
            state_colors (dict[StateId, str] | None, optional): An optional colormap for states. Defaults to None.
            state_combination_colors (dict[StateCombination, str] | None, optional): An optional colormap for state combinations. Defaults to None.

        Returns:
            StateSelection: The selection update with the new members set. Can either be a copy if `inplace=False` or the old instance with updated members otherwise.
        """
        if inplace:
            # Update and create
            if states:
                self.states = states
            if state_types:
                self.state_types = state_types
            if state_names:
                self.state_names = state_names
            if state_charges:
                self.state_charges = state_charges
            if state_combinations:
                self.state_combinations = state_combinations
            if state_combination_names:
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
            if state_colors:
                self.state_colors = state_colors
            if state_combination_colors:
                self.state_combination_colors = state_combination_colors

            return self
        else:
            if not states:
                states = self.states
            if not state_types:
                state_types = self.state_types
            if not state_names:
                state_names = self.state_names
            if not state_charges:
                state_charges = self.state_charges
            if not state_combinations:
                state_combinations = self.state_combinations
            if not state_combination_names:
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
            if not state_colors:
                state_colors = self.state_colors
            if not state_combination_colors:
                state_combination_colors = self.state_combination_colors

            return type(self)(
                states=states,
                state_types=state_types,
                state_names=state_names,
                state_charges=state_charges,
                state_combinations=state_combinations,
                state_combination_names=state_combination_names,
                state_colors=state_colors,
                state_combination_colors=state_combination_colors,
            )

    @classmethod
    def init_from_dataset(cls: type[Self], dataset: xr.Dataset) -> Self:
        """Alternative constructor that creates an initial StateSelection object from a dataset using the entire state information in it.

        Args:
            cls (type[StateSelection]): The type of this StateSelection so that we can create instances of it.
            dataset (xr.Dataset): The dataset to extract the state information out of. Must have a `state` dimension and preferrably coordinates `state`, `state_names`, `state_types`, `state_charges`, and `statecomb` set.
            If `state` is not set as a coordinate, a potential dimension size of `state` is taken and states are enumerates `1` through `1+dataset.sizes['state']`.
            If `statecomb` is not set as a coordinate, all unordered pairs of states will be used as a default value for `state_combinations`.

        Raises:
            ValueError: If no `state` information could be extracted from the dataset

        Returns:
            StateSelection: A state selection object initially covering all states (and state combinations) present in the dataset.
        """
        assert 'state' in dataset.sizes, (
            "No state information on the provided dataset. Cannot initialize state selection."
        )

        if 'states' in dataset.coords:
            states = list(dataset.coords['states'].values)
        elif 'state' in dataset.sizes:
            states = list(np.arange(1, 1 + dataset.sizes['state'], dtype=StateId))
        else:
            raise ValueError(
                "No sufficient state information on the provided dataset. Cannot initialize state selection."
            )

        if 'state_types' in dataset.coords:
            state_types = {
                state_id: type_val
                for (state_id, type_val) in zip(
                    states, dataset.coords['state_types'].values
                )
            }
        else:
            logging.warning(
                "No state types vailable on the dataset. Please assign them yourself."
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
                "No state names vailable on the dataset. Please assign them yourself."
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
                "No state charges vailable on the dataset. Please assign them yourself."
            )
            state_charges = None

        if 'statecomb' in dataset.coords:
            state_combinations = list(dataset.coords['statecomb'].values)
        else:
            state_combinations = list(combinations(states, 2))

        if state_names is not None:
            state_combination_names = {
                (a, b): f"{state_names[a]} - {state_names[b]}"
                for (a, b) in state_combinations
            }
        else:
            state_combination_names = None

        # Create an initial state selection
        return cls(
            states=states,
            state_types=state_types,
            state_charges=state_charges,
            state_names=state_names,
            state_combinations=state_combinations,
            state_combination_names=state_combination_names,
        )

    def filter_states(
        self,
        ids: Iterable[StateId] | StateId | None = None,
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
        inplace: bool = False,
    ) -> Self:
        """
        Method to get a new state selection only retaining the states satisfying the required inclusion criteria and
        not satisfying the exclusion criteria.

        Will return a new StateSelection object with the resulting set of states.

        Args:
            ids (Iterable[StateId] | StateId | None, optional): Explicit ids of states to retain. Either a single id or an iterable collection of state ids can be provided. Defaults to None.
            exclude_ids (Iterable[StateId] | StateId | None, optional):  Explicit ids of states to exclude. Either a single id or an iterable collection of state ids can be provided. Defaults to None.
            charge (Iterable[int] | int | None, optional): Charges of states to retain. Defaults to None.
            exclude_charge (Iterable[int] | int | None, optional): Charges of states to exclude. Defaults to None.
            multiplicity (Iterable[int] | int | None, optional): Multiplicity of states to retain. Defaults to None.
            exclude_multiplicity (Iterable[int] | int | None, optional): Multiplicity of states to exclude. Defaults to None.
            min_states_in_selection (Literal[0, 1, 2], optional): Optional parameter to determine whether state combinations should be kept if states they include are no longer part of the selection.
                A state combination is retained if at least `min_states_in_selection` of their states are still within the state selection. Defaults to 0, meaning all combinations are kept.
            inplace (bool, optional): Flag to update the selection in-place. Defaults to False, meaning a modified copy is returned.

        Returns:
            StateSelection: The resulting selection after applying all of the requested conditions.
        """
        new_states = list(self.states)
        if ids:
            if isinstance(ids, StateId):
                ids = [ids]
            next_states = []
            for old_state in new_states:
                if old_state in ids:
                    next_states.append(old_state)

            new_states = next_states

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

        def mult_label_transl(multipl: Iterable[int | MultiplicityLabel]) -> set[int]:
            """Function to translate potential string-based multiplicities to integers

            Args:
                multipl (Iterable[int | MultiplicityLabel]): List of multiplicities, either ints or string labels

            Returns:
                set[int]: A set representation of the numeric multiplicities
            """
            mult_translate = []
            for mult in multipl:
                if isinstance(mult, int):
                    mult_translate.append(mult)

                elif isinstance(mult, str):
                    lower_label = mult.lower()
                    if lower_label.startswith("s"):
                        mult_translate.append(1)
                    elif lower_label.startswith("d"):
                        mult_translate.append(2)
                    elif lower_label.startswith("t"):
                        mult_translate.append(3)
                    else:
                        raise ValueError(
                            f"Label `{mult}` is not a valid multiplicity label."
                        )
                else:
                    raise ValueError(
                        f"Invalid state type {mult} of object type {type(mult)}."
                    )

            res = set(mult_translate)
            return res

        if multiplicity:
            if isinstance(multiplicity, int) or isinstance(multiplicity, str):
                multiplicity = [multiplicity]

            trans_mult = mult_label_transl(multiplicity)

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

            trans_mult_excl = mult_label_transl(exclude_multiplicity)

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

        return self.copy_or_update(
            states=new_states, inplace=inplace
        ).filter_state_combinations(
            min_states_in_selection=min_states_in_selection, inplace=inplace
        )

    def filter_state_combinations(
        self,
        *,
        ids: Iterable[StateCombination] | None = None,
        min_states_in_selection: Literal[0, 1, 2] = 0,
        inplace: bool = True,
    ) -> Self:
        """Method to get a new state selection with a potentially reduced set of state combinations.

        Args:
            ids (Iterable[StateCombination] | None, optional): Explicit state transitions ids to retain. Defaults to None.
            min_states_in_selection (Literal[0, 1, 2], optional): Minimum number of states involved in the state combination that still need to be within the state selection to keep this combination. Defaults to 0, meaning no check will be performed.

        Returns:
            StateSelection: A new state selection with potentially fewer state combinations considered.
        """

        new_state_combinations = self.state_combinations

        if ids:
            # Filter explicit states
            next_state_combinations = []
            for old_comb in new_state_combinations:
                if old_comb in ids:
                    next_state_combinations.append(old_comb)

            new_state_combinations = next_state_combinations

        if min_states_in_selection > 0:
            # Check that there are sufficiently many states of the combination still in teh selection
            retained_combs = []
            for comb in new_state_combinations:
                states = set(comb)
                num_selected_states = len(states.intersection(self.states))

                if num_selected_states >= min_states_in_selection:
                    retained_combs.append(comb)

            new_state_combinations = retained_combs

        return self.copy_or_update(
            state_combinations=new_state_combinations, inplace=inplace
        )

    def set_state_names(
        self, names: Sequence[str] | dict[StateId, str], inplace: bool = True
    ) -> Self:
        """Helper function to assign new state names to the selection.

        Will peform some sanity checks first.

        Args:
            names (Sequence[str] | dict[StateId, str]): Either a list of state names aligned with `self.states` ids or a dictionary mapping state ids to names.
            inplace (bool, optional): Flag to determine whether this function should update the existing selection sequence or return a modified copy. Defaults to True, meaning the existing instance is updated.

        Raises:
            ValueError: If a Sequence is provided that does not have enough values
            ValueError: If a dict is  provided that does not have mapping for all state ids in `self.states`

        Returns:
            Self: Either the existing selection with updated names or a new instance with modified names.
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

        Args:
            types (Sequence[int] | dict[StateId, int]): Either a list of state types/multiplicities aligned with `self.states` ids or a dictionary mapping state ids to types.
            inplace (bool, optional): Flag to determine whether this function should update the existing selection sequence or return a modified copy. Defaults to True, meaning the existing instance is updated.

        Raises:
            ValueError: If a Sequence is provided that does not have enough values
            ValueError: If a dict is  provided that does not have mapping for all state ids in `self.states`

        Returns:
            Self: Either the existing selection with updated types or a new instance with modified types.
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

        Args:
            charges (int| Sequence[int] | dict[StateId, int]): Either a single charge for all states or a list of state charges aligned with `self.states` ids or a dictionary mapping state ids to charges.
            inplace (bool, optional): Flag to determine whether this function should update the existing selection sequence or return a modified copy. Defaults to True, meaning the existing instance is updated.

        Raises:
            ValueError: If a Sequence is provided that does not have enough charges
            ValueError: If a dict is provided that does not have mapping for all state ids in `self.states`

        Returns:
            Self: Either the existing selection with updated charges or a new instance with modified charges.
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

        Args:
            combinations (Sequence[StateCombination]): A list of state combination tuples to set to the selection
            inplace (bool, optional): Flag to determine whether this function should update the existing selection or return a modified copy. Defaults to True, meaning the existing instance is updated.

        Raises:
            ValueError: If an entry in the combinations sequence has a non-positive state entry.

        Returns:
            Self: Either the existing selection with updated combinations or a new instance with modified combinations.
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

        Args:
            names (Sequence[str] | dict[StateCombination, str]): Either a list of state combination names aligned with `self.state_combinations` or a dictionary mapping state combination ids to names.
            inplace (bool, optional): Flag to determine whether this function should update the existing selection or return a modified copy. Defaults to True, meaning the existing instance is updated.

        Raises:
            ValueError: If a Sequence is provided that does not have enough values
            ValueError: If a dict is  provided that does not have mapping for all state combination ids in `self.state_combinations`

        Returns:
            Self: Either the existing selection with updated names or a new instance with modified names.
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

        Args:
            inplace (bool, optional): Flag whether the operation should update the selection in-place. Defaults to False.

        Returns:
            StateSelection: the updated selection only containing singlet states.
        """
        return self.filter_states(multiplicity=1, inplace=inplace)

    def triplets_only(self, inplace: bool = False) -> Self:
        """Helper function to immediately filter only triplet states. Does not affect state combinations.

        Args:
            inplace (bool, optional): Flag whether the operation should update the selection in-place. Defaults to False.

        Returns:
            StateSelection: the updated selection only containing triplet states.
        """
        return self.filter_states(multiplicity=1, inplace=inplace)

    def same_multiplicity_transitions(self, inplace: bool = False) -> Self:
        """Helper function to only retain combinations between states of the same multiplicities (e.g. for NACs)

        Args:
            inplace (bool, optional): Flag whether the operation should update the selection in-place. Defaults to False.

        Returns:
            StateSelection: the updated selection only containing transitions between states of same multiplicity (i.e. singlet-singlet, triplet-tiplet).
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

        Args:
            inplace (bool, optional): Flag whether the operation should update the selection in-place. Defaults to False.

        Returns:
            StateSelection: the updated selection only containing transitions between states of different multiplicity (i.e. singlet-triplet).
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

        Args:
            ground_state_id (StateId, optional): Id of the state to be considered the ground state. Defaults to the lowest id of the selected states.
            inplace (bool, optional): Flag whether the operation should update the selection in-place. Defaults to False.

        Returns:
            StateSelection: the updated selection only containing transitions between ground state and other states.
        """

        if ground_state_id is None:
            ground_state_id = np.min(self.states)

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

        Args:
            ground_state_id (StateId, optional): Id of the state to be considered the ground state. Defaults to the lowest id of the selected states.
            inplace (bool, optional): Flag whether the operation should update the selection in-place. Defaults to False.

        Returns:
            StateSelection: the updated selection only containing transitions between non-ground states.
        """

        if ground_state_id is None:
            ground_state_id = np.min(self.states)

        new_state_combs = []
        for comb in self.state_combinations:
            first, second = comb

            if first != ground_state_id and second != ground_state_id:
                new_state_combs.append(comb)

        return self.copy_or_update(state_combinations=new_state_combs, inplace=inplace)

    def state_info(self) -> Iterable[StateInfo]:
        """Get an iterator over the states in this selection.

        Returns:
            Iterable[StateInfo]: An iterator over the available state info
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

        Args:
            id (StateId): Id of the state to get the name for.

        Returns:
            str: Label of the state
        """
        if self.state_names:
            if id in self.state_names:
                return self.state_names[id]

        return f"state{id - 1}"

    def get_state_combination_name_or_default(self, comb: StateCombination) -> str:
        """Helper method to either get registered state combination name or a default string to identify the state combination.

        Args:
            comb (StateCombination): Id of the state combination to get the name for.

        Returns:
            str: Label of the state combination
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

        Args:
            id (StateId): Id of the state to get the label for

        Returns:
            str: Tex-label that needs to be enclosed in a math environment to not cause issues.
        """

        statename = self.get_state_name_or_default(id)
        return state_name_to_tex_label(statename)

    def get_state_combination_tex_label(self, comb: StateCombination) -> str:
        """Function to get a nice tex-printable label with super- and subscripts for a state combination in this selection

        Args:
            comb (StateCombination): Combination identifier to get the label for

        Returns:
            str: Tex-label that needs to be enclosed in a math environment to not cause issues.
        """
        first, second = comb

        s1 = self.get_state_tex_label(first)
        s2 = self.get_state_tex_label(second)
        return f"{s1} - {s2}"

    def combination_info(self) -> Iterable[StateCombInfo]:
        """Get an iterator over the state combinations in this selection.

        Returns:
            Iterable[StateCombInfo]: An iterator over the available state combination info
        """

        for comb in self.state_combinations:
            name = self.get_state_combination_name_or_default(comb)
            yield StateCombInfo(comb, name)

    def has_state(self, id: StateId) -> bool:
        """Function to check whether a state is in the selection

        Args:
            id (StateId): The state id to check whether it has been selected

        Returns:
            bool: True if in the selection, False otherwise.
        """
        return id in self.states

    def has_state_combination(self, comb: StateCombination) -> bool:
        """Function to check whether a state combination is in the selection

        Args:
            comb (StateCombination): The combination to check whether it has been selected

        Returns:
            bool: True if in the selection, False otherwise.
        """
        return comb in self.state_combinations

    def auto_assign_colors(self, inplace: bool = True) -> Self:
        """Function to automatically generate colors for states and state combinations

        Args:
            inplace (bool, optional): Flag whether the operation should update the selection in-place. Defaults to True because setting colors is not a big issue.

        Returns:
            Self: Returns the updated instance.
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

            mult_color_maps[mult] = {
                state_id: color
                for state_id, color in zip(
                    state_list,
                    get_default_state_colormap(len(state_list), multiplicity=mult),
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

        Args:
            id (StateId): Id of the state to get the color for

        Returns:
            str: Hex-str color code
        """
        from shnitsel.vis.colormaps import st_grey

        if self.state_colors is not None and id in self.state_colors:
            return self.state_colors[id]
        else:
            return st_grey

    def get_state_combination_color(self, comb: StateCombination) -> str:
        """Function to get a the state combination color or a default color value

        Args:
            comb (StateCombination): Id of the state combination to get the color for

        Returns:
            str: Hex-str color code
        """
        from shnitsel.vis.colormaps import st_grey

        if (
            self.state_combination_colors is not None
            and comb in self.state_combination_colors
        ):
            return self.state_combination_colors[comb]
        else:
            return st_grey

    # TODO: FIXME: Add print output __str__, __html__ and __repr__
