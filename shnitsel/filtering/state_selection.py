from dataclasses import dataclass
from itertools import combinations
import logging
from typing import Iterable, Self, Sequence, TypeAlias, Literal

import numpy as np
import xarray as xr

StateId: TypeAlias = int
StateCombination: TypeAlias = tuple[StateId, StateId]

MultiplicityLabel: TypeAlias = Literal[
    's', 'S', 'singlet', 'd', 'D', 'doublet', 't', 'T', 'triplet'
]


@dataclass
class StateSelection:
    """Class to keep track of a (sub-)selection of states and state transitions"""

    states: Sequence[StateId]
    state_types: dict[StateId, int] | None
    state_names: dict[StateId, str] | None
    state_charges: dict[StateId, int] | None

    state_combinations: list[StateCombination]
    state_combination_names: dict[StateCombination, str] | None

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
        exclude_multiplicity: Iterable[int] | int | None = None,
        combinations_min_states_in_selection: Literal[0, 1, 2] = 0,
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
            combinations_min_states_in_selection (Literal[0, 1, 2], optional): Optional parameter to determine whether state combinations should be kept if states they include are no longer part of the selection.
                A state combination is retained if at least `combinations_min_states_in_selection` of their states are still within the state selection. Defaults to 0, meaning all combinations are kept.

        Returns:
            StateSelection: The resulting selection after applying all of the requested conditions.
        """
        # TODO: implement filtering
        return self

    def filter_state_combinations(
        self,
        *,
        ids: Iterable[StateCombination] | StateCombination | None = None,
        min_states_in_selection: Literal[0, 1, 2] = 0,
    ) -> Self:
        """Method to get a new state selection with a potentially reduced set of state combinations.

        Args:
            ids (Iterable[StateCombination] | StateCombination | None, optional): Explicit state transitions ids to retain. Defaults to None.
            min_states_in_selection (Literal[0, 1, 2], optional): Minimum number of states involved in the state combination that still need to be within the state selection to keep this combination. Defaults to 0.

        Returns:
            StateSelection: A new state selection with potentially fewer state combinations considered.
        """
        # TODO: implement filtering
        return self
