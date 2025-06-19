import ROOT
import pandas as pd
import yaml
import argparse
import os
import re
from itertools import combinations
import numpy as np

import sys

sys.path.append("utils")
import utils as utils

parser = argparse.ArgumentParser(description="Configure the parameters of the script.")
parser.add_argument(
    "--config",
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

input_file_name = config["input_file_name"]
resolution_file_name = config["resolution_file_name"]
output_dir_name = config["output_dir_name"]

nuclei_tree_name = config["nuclei_tree_name"]
ep_tree_name = config["ep_tree_name"]

useSP = True  # config["useSP"]

mandatory_selections = config["mandatory_selections"]
selection_dict = config["selection_dict"]
selection_list = selection_dict.values()
selections = " and ".join(selection_list)

pattern = r"\s*and\s+(abs\(fNsigmaTPC3He\)\s*<\s*(\d+\.?\d*))"

# Search and extract the match without "and"
match = re.search(pattern, selections)
nSigmaTPC_selection_string = match.group(1).strip() if match else ""
n_sigma_selection_from_string = float(match.group(2)) if match else None

selections = re.sub(pattern, "", selections).strip()

# ptdep_selection_dict = config["ptdep_selection_dict"]["fAvgItsClusSizeCosLambda"]

cent_detector_label = config["cent_detector_label"]
reference_flow_detector = config["reference_flow_detector"]
resolution_flow_detectors = config["resolution_flow_detectors"]

centrality_classes = config["centrality_classes"]
pt_bins = config["pt_bins"]
cent_colours = config["cent_colours"]


# create output file
if not os.path.exists(output_dir_name):
    os.makedirs(output_dir_name)
output_file = ROOT.TFile(f"{output_dir_name}/phi_studies.root", "recreate")

# get a unique df from nuclei and ep trees
nuclei_df = utils.get_df_from_tree(input_file_name, nuclei_tree_name)

nucleiflow_df = utils.get_df_from_tree(input_file_name, ep_tree_name)

complete_df = pd.concat([nuclei_df, nucleiflow_df], axis=1, join="inner")

# define new columns
utils.redefineColumns(complete_df, useSP=useSP)

# apply common selections
complete_df.query(f"{mandatory_selections} and {selections}", inplace=True)

for i_cent, cent in enumerate(centrality_classes):
    cent_selection = f"fCentFT0C > {cent[0]} and fCentFT0C < {cent[1]}"
    cent_dir = output_file.mkdir(f"{cent[0]}_{cent[1]}")
    cent_df = complete_df.query(cent_selection, inplace=False)
    vPhi_cent = []
    phi_dir = cent_dir.mkdir("phi")

    pt_bins_arr = np.array(pt_bins[i_cent], dtype=np.float64)
    n_pt_bins = len(pt_bins[i_cent]) - 1

    for i_pt in range(0, n_pt_bins):
        # select the correct pt bin
        pt_sel = f"abs(fPt) > {pt_bins[i_cent][i_pt]} and abs(fPt) < {pt_bins[i_cent][i_pt+1]}"  # and {ptdep_selection_list[i_pt]}"
        bin_df = cent_df.query(pt_sel, inplace=False)

        hPhi = ROOT.TH1F(
            f"hPhi_cent_{cent[0]}_{cent[1]}_pt_{i_pt}", r";#phi; counts", 140, -7.0, 7.0
        )
        utils.setHistStyle(hPhi, ROOT.kRed + 2)

        for phi in bin_df["fPhi"]:
            hPhi.Fill(phi)

        vPhi_cent.append(hPhi)
        phi_dir.cd()
        hPhi.Write()

        hPhiMinusPsi_formatted_title = f";phi - #Psi^{{FT0C}}_{{2}}; counts"
        hPhiMinusPsi = ROOT.TH1F(
            f"hPhiMinusPsiFT0C_cent_{cent[0]}_{cent[1]}_pt_{i_pt}",
            hPhiMinusPsi_formatted_title,
            40,
            0.0,
            ROOT.TMath.Pi(),
        )
        utils.setHistStyle(hPhiMinusPsi, ROOT.kRed + 2)

        for phi, q, psi, v2 in zip(bin_df["fPhi"], bin_df[f"fPsiFT0C"]):
            delta_phi = utils.getPhiInRange(phi - psi)
            sinPhiMinusPsi = ROOT.TMath.Sin(delta_phi)
            hPhiMinusPsi.Fill(delta_phi)

        cent_dir.cd()
        hPhiMinusPsi.Write()
