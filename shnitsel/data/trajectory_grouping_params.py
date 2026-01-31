from dataclasses import dataclass


@dataclass(unsafe_hash=True)
class TrajectoryGroupingMetadata:
    delta_t_in_fs: float | None
    input_format_name: str | None
    input_format_version: str | None
    est_level: str | None
    theory_basis_set: str | None
    charge_in_e: float | None
    num_states: int | None
