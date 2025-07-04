import os
import shnitsel as sh
from shnitsel.plot import Datasheet


def test_sharc_iconds():
    # parse icond data
    iconds_butene = sh.parse.sharc_icond.dirs_of_iconds(
        path='./test_data/sharc/iconds_butene'
    )
    # save the parsed data to h5netcdf file
    savepath = os.path.join(os.getcwd(), 'test_data', 'sharc', 'iconds_butene.hdf5')
    iconds_butene.reset_index('statecomb').to_netcdf(savepath, engine='h5netcdf')


def test_sharc_trajs():
    # parse trajectory data from SHARC output files
    traj_frames_butene = sh.parse.read_trajs(
        'test_data/sharc/traj_butene', kind='sharc'
    )


def test_nx():
    # parse trajectory data from Newton-X output files
    traj_frames_chd = sh.parse.read_trajs('test_data/newtonx/', kind='nx')


def test_open():
    sh.open_frames('test_data/butene.nc')


def test_biplot():
    # load trajectory data of A01
    a01 = sh.open_frames('A01_ethene_dynamic.nc')
    # create PCA plot over all trajectories with visualization of the
    # four most important PCA-axis on the molecular structure
    # C=C bond color highlighgting via KDE in PCA
    sh.plot.biplot_kde(
        frames=a01, at1=0, at2=1, geo_filter=[[0, 3], [5, 20]], levels=10
    )
    # C-H bond color highlighting via KDE in PCA
    sh.plot.biplot_kde(
        frames=a01, at1=0, at2=2, geo_filter=[[0, 3], [5, 20]], levels=10
    )


def test_per_state_histograms():
    sheet = sh.plot.Datasheet(path='A01_ethene_dynamic.nc')
    sheet.plot_per_state_histograms()


def test_nacs_histograms():
    sheet_a01 = Datasheet(path='A01_ethene_dynamic.nc')
    sheet_i01 = Datasheet(path='I01_ch2nh2_dynamic.nc')
    # inter-state histograms
    sheet_a01.plot_nacs_histograms()
    sheet_i01.plot_nacs_histograms()


def test_ski_plots():
    # load data
    spectra_data = (
        sh.open_frames(path='A01_ethene_dynamic.nc')
        .sh.get_inter_state()
        .sh.spectra_all_times()
    )
    # plot spectra at different simulation times in one plot with a dahsed line that tracks the maximum
    sh.plot.ski_plots(spectra_data)


def test_timeplots():
    # load data
    sheet_a01 = Datasheet(path='A01_ethene_dynamic.nc')
    sheet_i01 = Datasheet(path='I01_ch2nh2_dynamic.nc')
    # time plots
    sheet_a01.plot_timeplots()
    sheet_i01.plot_timeplots()


def test_datasheet_full():
    # load data
    sheet_a01 = Datasheet(path='A01_ethene_dynamic.nc')
    # automatic generation of datasheet
    sheet_a01.charge = 0
    sheet_a01.plot()