import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
import seaborn as sns

def plot_stack(df_all, trajid_list, time='time', 
        figsize=(7, 8), show_hops=True):
    """
    Plot stacked (gap, BLA, HOOP) panels for multiple trajectories
    in a single combined figure.

    Parameters
    ----------
    df_all : pandas.DataFrame
        Must contain columns:
        ['trajid', 'time', 'gap', 'BLA', 'HOOP', 'dihedral', 'astate']
    trajid_list : list of int
        List of trajectory IDs to include.
    figsize : tuple
        Figure size.
    show_hops : bool
        Whether to draw vertical lines at hopping times.
    """

    # Set global font sizes
    plt.rcParams.update({
        "font.size": 14,
        "axes.titlesize": 14,
        "axes.labelsize": 14,
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
        "legend.fontsize": 12,
    })

    # Marker mapping for electronic state
    marker_map = {3: 's', 2: "o", 1: "D"}

    # Filter trajectories
    df_sel = df_all[df_all["trajid"].isin(trajid_list)].copy()
    if df_sel.empty:
        raise ValueError("No data found for the selected trajid_list.")

    # Color mapping for dihedral
    dihedral_norm = mpl.colors.Normalize(vmin=0, vmax=180)
    cmap = mpl.colormaps.get_cmap("vanimo")

    # --- Create stacked plot ---
    fig, axes = plt.subplots(
        3, 1,
        figsize=figsize,
        sharex=True,
        gridspec_kw={"hspace": 0}
    )
    plt.subplots_adjust(right=0.88)

    # ========= 1) Energy gap =========
    sns.scatterplot(
        data=df_sel, x=time, y="gap",
        hue="dihedral", style="astate",
        hue_norm=dihedral_norm,
        palette="vanimo", markers=marker_map,
        ax=axes[0], s=25, edgecolor=None, legend=False
    )
    axes[0].set_ylabel(r"$\Delta E_{01}$ / eV")

    # ========= 2) BLA =========
    sns.scatterplot(
        data=df_sel, x=time, y="BLA",
        hue="dihedral", style="astate",
        hue_norm=dihedral_norm,
        palette="vanimo", markers=marker_map,
        ax=axes[1], s=25, edgecolor=None, legend=False
    )
    axes[1].set_ylabel("BLA / Bohr")

    # ========= 3) HOOP =========
    sns.scatterplot(
        data=df_sel, x=time, y="HOOP",
        hue="dihedral", style="astate",
        hue_norm=dihedral_norm,
        palette="vanimo", markers=marker_map,
        ax=axes[2], s=25, edgecolor=None, legend=False
    )
    axes[2].set_ylabel("HOOP / °")
    axes[2].set_xlabel("time / fs")

    # ========= hopping times =========
    if show_hops:
        for trajid in trajid_list:
            df_i = df_sel[df_sel["trajid"] == trajid]
            astate_changes = np.where(np.diff(df_i["astate"].values) != 0)[0]
            if len(astate_changes):
                hop_time = df_i[time].values[astate_changes[-1]]
                for ax in axes:
                    ax.axvline(hop_time, color="black", linestyle="--", alpha=0.3)

    # ========= COLORBAR =========
    cbar = fig.colorbar(
        mpl.cm.ScalarMappable(norm=dihedral_norm, cmap=cmap),
        ax=axes, location="right", pad=0.02
    )
    cbar.set_label(r"$\varphi_{CC=CC}$ / °")

    fig.suptitle(
        f"Combined trajectories: {', '.join(map(str, trajid_list))}",
        y=0.92
    )

    plt.show()



def align_trajectories_to_last_hop(
    df_all,
    from_state=2,
    to_state=1,
    time_col="time",
    state_col="astate",
    traj_col="trajid"
):
    """
    Align trajectories such that the LAST hop from from_state -> to_state
    occurs at t = 0.

    Returns
    -------
    df_aligned : pandas.DataFrame
        Copy of df_all with an additional column:
        'time_aligned'
    hop_times : dict
        Dictionary mapping trajid -> last hop time
    """

    df_aligned = []
    hop_times = {}

    for trajid, df_i in df_all.groupby(traj_col):
        df_i = df_i.sort_values(time_col).copy()

        astate = df_i[state_col].values
        times = df_i[time_col].values

        # find ALL from_state -> to_state transitions
        hop_idxs = np.where(
            (astate[:-1] == from_state) & (astate[1:] == to_state)
        )[0]

        if len(hop_idxs) == 0:
            # no hop found → skip trajectory
            continue

        # take the LAST hop
        last_hop_idx = hop_idxs[-1]
        hop_time = times[last_hop_idx]

        hop_times[trajid] = hop_time

        # shift time axis
        df_i["time_aligned"] = df_i[time_col] - hop_time
        df_aligned.append(df_i)

    if not df_aligned:
        raise ValueError("No trajectories with the requested hop were found.")

    df_aligned = pd.concat(df_aligned, ignore_index=True)
    return df_aligned, hop_times


def label_trajectory(df, dihedral_col="dihedral", lower_bound=70, upper_bound=110):
    """
    Add an 'isomerization' column based on the final dihedral per trajectory,
    without using pandas apply.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain 'trajid' and the dihedral column.
    dihedral_col : str
        Column name of the dihedral values.
    lower_bound : float
        Dihedral below this value is labeled 'E-to-Z'.
    upper_bound : float
        Dihedral above this value is labeled 'E-to-E'.

    Returns
    -------
    df_labeled : pd.DataFrame
        Original DataFrame with additional 'isomerization' column.
    """

    df_labeled = df.copy()
    df_labeled["isomerization"] = "undetermined"  # default value

    for trajid, df_i in df_labeled.groupby("trajid"):
        final_dihedral = df_i[dihedral_col].iloc[-1]

        if final_dihedral < lower_bound:
            label = "E-to-Z"
        elif final_dihedral > upper_bound:
            label = "E-to-E"
        else:
            label = "undetermined"

        # assign label to all rows of this trajectory
        df_labeled.loc[df_labeled["trajid"] == trajid, "isomerization"] = label

    return df_labeled


import pandas as pd
import numpy as np

def get_last_hops_from_states(df, state_from=2, state_to=1, traj_col="trajid", state_col="astate"):
    """
    For each trajectory, find the last hop from state_from -> state_to,
    and return the row corresponding to the last state_from before the hop.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain trajectory IDs and state column.
    state_from : int
        Initial state of the hop (default 2).
    state_to : int
        Final state of the hop (default 1).
    traj_col : str
        Column name for trajectory ID.
    state_col : str
        Column name for electronic state.

    Returns
    -------
    df_hop : pd.DataFrame
        One row per trajectory: last state_from before last hop.
    """

    records = []

    for trajid, df_i in df.groupby(traj_col):
        df_i = df_i.sort_values("time_aligned").reset_index(drop=True)
        astate = df_i[state_col].values

        # Find all indices where state_from -> state_to occurs
        hop_idxs = np.where((astate[:-1] == state_from) & (astate[1:] == state_to))[0]

        if len(hop_idxs) == 0:
            continue  # skip trajectories without any hop

        # Last hop
        last_hop_idx = hop_idxs[-1]

        # Keep the row corresponding to the last state_from before the hop
        row_before_hop = df_i.iloc[last_hop_idx].copy()
        records.append(row_before_hop)

    if not records:
        raise ValueError("No trajectories with the requested hop found.")

    df_hop = pd.DataFrame(records).reset_index(drop=True)
    return df_hop


