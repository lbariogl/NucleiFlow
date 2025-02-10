import ROOT
import pandas as pd
import numpy as np

import sys

sys.path.append("utils")
import utils as utils

# 0 - 10% centrality

input_file_name_0_10 = "../theoretical_models/v2_pt_Helium3_5360GeV_PbPb_010.txt"

df_0_10 = pd.read_csv(input_file_name_0_10, index_col=False, sep=" ")

pt_0_10 = df_0_10["pt"].tolist()
v2_0_10 = df_0_10["v2"].tolist()
v2_err_0_10 = df_0_10["v2_err"].tolist()

n_points_0_10 = len(v2_err_0_10)

pt_0_10_array = np.array(pt_0_10, dtype=np.float64)
v2_0_10_array = np.array(v2_0_10, dtype=np.float64)
v2_err_0_10_array = np.array(v2_err_0_10, dtype=np.float64)
zeros_0_10_array = np.zeros(n_points_0_10, dtype=np.float64)

gPredWenbin010 = ROOT.TGraphErrors(
    n_points_0_10, pt_0_10_array, v2_0_10_array, zeros_0_10_array, v2_err_0_10_array
)
gPredWenbin010.SetName("gPredWenbin_0_10")
utils.setHistStyle(gPredWenbin010, ROOT.kCyan + 2)

# 10 - 20% centrality

input_file_name_10_20 = "../theoretical_models/v2_pt_Helium3_5360GeV_PbPb_1020.txt"

df_10_20 = pd.read_csv(input_file_name_10_20, index_col=False, sep=" ")

pt_10_20 = df_10_20["pt"].tolist()
v2_10_20 = df_10_20["v2"].tolist()
v2_err_10_20 = df_10_20["v2_err"].tolist()

n_points_10_20 = len(v2_err_10_20)

pt_10_20_array = np.array(pt_10_20, dtype=np.float64)
v2_10_20_array = np.array(v2_10_20, dtype=np.float64)
v2_err_10_20_array = np.array(v2_err_10_20, dtype=np.float64)
zeros_10_20_array = np.zeros(n_points_10_20, dtype=np.float64)

gPredWenbin1020 = ROOT.TGraphErrors(
    n_points_10_20,
    pt_10_20_array,
    v2_10_20_array,
    zeros_10_20_array,
    v2_err_10_20_array,
)
gPredWenbin1020.SetName("gPredWenbin_10_20")
utils.setHistStyle(gPredWenbin1020, ROOT.kGreen)

# 20 - 30% centrality

input_file_name_20_30 = "../theoretical_models/v2_pt_Helium3_5360GeV_PbPb_2030.txt"

df_20_30 = pd.read_csv(input_file_name_20_30, index_col=False, sep=" ")

pt_20_30 = df_20_30["pt"].tolist()
v2_20_30 = df_20_30["v2"].tolist()
v2_err_20_30 = df_20_30["v2_err"].tolist()

n_points_20_30 = len(v2_err_20_30)

pt_20_30_array = np.array(pt_20_30, dtype=np.float64)
v2_20_30_array = np.array(v2_20_30, dtype=np.float64)
v2_err_20_30_array = np.array(v2_err_20_30, dtype=np.float64)
zeros_20_30_array = np.zeros(n_points_20_30, dtype=np.float64)

gPredWenbin2030 = ROOT.TGraphErrors(
    n_points_20_30,
    pt_20_30_array,
    v2_20_30_array,
    zeros_20_30_array,
    v2_err_20_30_array,
)
gPredWenbin2030.SetName("gPredWenbin_20_30")
utils.setHistStyle(gPredWenbin2030, ROOT.kOrange)

# 30 - 40% centrality

input_file_name_30_40 = "../theoretical_models/v2_pt_Helium3_5360GeV_PbPb_3040.txt"

df_30_40 = pd.read_csv(input_file_name_30_40, index_col=False, sep=" ")

pt_30_40 = df_30_40["pt"].tolist()
v2_30_40 = df_30_40["v2"].tolist()
v2_err_30_40 = df_30_40["v2_err"].tolist()

n_points_30_40 = len(v2_err_30_40)

pt_30_40_array = np.array(pt_30_40, dtype=np.float64)
v2_30_40_array = np.array(v2_30_40, dtype=np.float64)
v2_err_30_40_array = np.array(v2_err_30_40, dtype=np.float64)
zeros_30_40_array = np.zeros(n_points_30_40, dtype=np.float64)

gPredWenbin3040 = ROOT.TGraphErrors(
    n_points_30_40,
    pt_30_40_array,
    v2_30_40_array,
    zeros_30_40_array,
    v2_err_30_40_array,
)
gPredWenbin3040.SetName("gPredWenbin_30_40")
utils.setHistStyle(gPredWenbin3040, ROOT.kRed + 1)

# 40 - 60% centrality

input_file_name_40_60 = "../theoretical_models/v2_pt_Helium3_5360GeV_PbPb_4060.txt"

df_40_60 = pd.read_csv(input_file_name_40_60, index_col=False, sep=" ")

pt_40_60 = df_40_60["pt"].tolist()
v2_40_60 = df_40_60["v2"].tolist()
v2_err_40_60 = df_40_60["v2_err"].tolist()

n_points_40_60 = len(v2_err_40_60)

pt_40_60_array = np.array(pt_40_60, dtype=np.float64)
v2_40_60_array = np.array(v2_40_60, dtype=np.float64)
v2_err_40_60_array = np.array(v2_err_40_60, dtype=np.float64)
zeros_40_60_array = np.zeros(n_points_40_60, dtype=np.float64)

gPredWenbin4060 = ROOT.TGraphErrors(
    n_points_40_60,
    pt_40_60_array,
    v2_40_60_array,
    zeros_40_60_array,
    v2_err_40_60_array,
)
gPredWenbin4060.SetName("gPredWenbin_40_60")
utils.setHistStyle(gPredWenbin4060, ROOT.kViolet + 1)


# saving to file
output_file = ROOT.TFile("../theoretical_models/Predictions.root", "recreate")
gPredWenbin010.Write()
gPredWenbin1020.Write()
gPredWenbin2030.Write()
gPredWenbin3040.Write()
gPredWenbin4060.Write()
