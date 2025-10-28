import os
import shnitsel as sh
from shnitsel.plot import Datasheet



class TestDatasheetFunctionality:
    """Tests for the Datasheet utility class
    """

    def test_per_state_histograms():
        sheet = sh.plot.Datasheet(path='A01_ethene_dynamic.nc')
        sheet.plot_per_state_histograms()


    def test_nacs_histograms():
        sheet_a01 = Datasheet(path='A01_ethene_dynamic.nc')
        sheet_i01 = Datasheet(path='I01_ch2nh2_dynamic.nc')
        # inter-state histograms
        sheet_a01.plot_nacs_histograms()
        sheet_i01.plot_nacs_histograms()

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