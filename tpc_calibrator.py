import ROOT
from hipe4ml.tree_handler import TreeHandler
import pandas as pd
import yaml
import argparse
import numpy as np
from itertools import product
import os

import sys

sys.path.append("utils")
import utils as utils

parser = argparse.ArgumentParser(description="Configure the parameters of the script.")
parser.add_argument(
    "--config-file",
    dest="config_file",
    help="path to the YAML file with configuration.",
    default="",
)
args = parser.parse_args()
if args.config_file == "":
    print("** No config file provided. Exiting. **")
    exit()

config_file = open(args.config_file, "r")
config = yaml.full_load(config_file)

common_input_path = config["common_input_path"]
input_dir_names = config["input_dir_names"]
input_file_name = config["input_file_name"]

output_dir_name = config["output_dir_name"]
output_file_name = config["output_file_name"]

nuclei_tree_name = config["nuclei_tree_name"]
ep_tree_name = config["ep_tree_name"]

mandatory_selections = config["mandatory_selections"]
selection_dict = config["selection_dict"]
selection_list = selection_dict.values()
selections = " and ".join(selection_list)

pdep_selection_dict = config["pdep_selection_dict"]
pdep_selections = pdep_selection_dict["fAvgItsClusSizeCosLambda"]

p_bins = config["p_bins"]
p_bins_array = np.array(config["p_bins"], dtype=np.float64)
n_p_bins = len(p_bins) - 1

default_parameters = config["default_parameters"]

# create output file
if not os.path.exists(output_dir_name):
    os.makedirs(output_dir_name)
output_file = ROOT.TFile(f"{output_dir_name}/{output_file_name}", "recreate")

for dir_name in input_dir_names:

    # getting input tree
    path = common_input_path + "/" + dir_name + "/" + input_file_name
    nuclei_hdl = TreeHandler(path, f"{nuclei_tree_name};", folder_name="DF*")
    nuclei_df = nuclei_hdl._full_data_frame

    # creating output directory
    dataset_dir = output_file.mkdir(dir_name)

    # define new columns
    utils.redefineColumnsLight(nuclei_df)

    # apply mandatory selections
    nuclei_df.query(mandatory_selections + " and " + selections, inplace=True)

    hTPCdEdXvsP = ROOT.TH2D(
        f"hTPCdEdXvsP_{dir_name}",
        r";#it{p}/z (GeV/#it{c}); d#it{E}/d#it{X} (a. u.)",
        25,
        1.0,
        6,
        175,
        0,
        1400,
    )

    hTPCdEdXvsP_toFit = ROOT.TH1D(
        f"hTPCdEdXvsP_toFit_{dir_name}",
        r";#it{p}/z (GeV/#it{c}); d#it{E}/d#it{X} (a. u.)",
        n_p_bins,
        p_bins_array,
    )
    utils.setHistStyle(hTPCdEdXvsP_toFit, ROOT.kBlack)

    for i_p in range(0, n_p_bins):

        p_bin = [p_bins[i_p], p_bins[i_p + 1]]
        p_sel = f"abs(fTPCInnerParam) > {p_bin[0]} and abs(fTPCInnerParam) < {p_bin[1]}"
        p_sel = p_sel + " and " + pdep_selections[i_p]
        print(f"psel: {p_sel}")
        bin_df = nuclei_df.query(p_sel, inplace=False)

        p_label = (
            f"{p_bin[0]:.2f} "
            + r"#leq #it{p}/z < "
            + f"{p_bin[1]:.2f}"
            + r" GeV/#it{c}"
        )

        hTPCdEdX_pbin = ROOT.TH1D(
            f"hTPCdEdX_p{i_p}",
            p_label + r";d#it{E}/d#it{X} (a. u.); counts",
            175,
            0,
            1400,
        )
        for dEdX, p in zip(bin_df["fTPCsignal"], bin_df[f"fTPCInnerParam"]):
            hTPCdEdXvsP.Fill(p, dEdX)
            hTPCdEdX_pbin.Fill(dEdX)

        mean = hTPCdEdX_pbin.GetMean()
        rms = hTPCdEdX_pbin.GetRMS()
        hTPCdEdX_pbin.Fit("gaus", "MQRL+", "", mean - 3 * rms, mean + 3 * rms)
        fitFunc = hTPCdEdX_pbin.GetFunction("gaus")
        hTPCdEdXvsP_toFit.SetBinContent(i_p + 1, fitFunc.GetParameter(1))
        hTPCdEdXvsP_toFit.SetBinError(i_p + 1, fitFunc.GetParameter(2))

        dataset_dir.cd()
        hTPCdEdX_pbin.Write()

    # Drawing default BB function
    func_BB_default = ROOT.TF1("func_BB_default", utils.func_string, 0.5, 6, 5)
    func_BB_default.SetParameters(
        default_parameters[0],
        default_parameters[1],
        default_parameters[2],
        default_parameters[3],
        default_parameters[4],
    )
    func_BB_default.SetLineColor(ROOT.kRed)

    # Defining BB function for fit
    func_BB_fit = ROOT.TF1("func_BB_fit", utils.func_string, 0.5, 6, 5)
    func_BB_fit.SetLineColor(ROOT.kBlue)
    func_BB_fit.SetParameters(
        default_parameters[0],
        default_parameters[1],
        default_parameters[2],
        default_parameters[3],
        default_parameters[4],
    )
    hTPCdEdXvsP_toFit.Fit(func_BB_fit)

    cTPCdEdXvsP = ROOT.TCanvas(
        f"cTPCdEdXvsP_{dir_name}", f"cTPCdEdXvsP_{dir_name}", 800, 600
    )
    cTPCdEdXvsP.DrawFrame(
        0.8, 0.0, 4.0, 1400.0, r";#it{p}/z (GeV/#it{c}); d#it{E}/d#it{X} (a. u.)"
    )
    hTPCdEdXvsP.Draw("colz same")
    func_BB_default.Draw("L same")
    func_BB_fit.Draw("L same")

    dataset_dir.cd()
    hTPCdEdXvsP.Write()
    cTPCdEdXvsP.Write()
    hTPCdEdXvsP_toFit.Write()
