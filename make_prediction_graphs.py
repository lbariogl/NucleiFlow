import ROOT
import pandas as pd
import numpy as np

import sys
sys.path.append('utils')
import utils as utils

# 0 - 20% centrality

input_file_name_0_20 = '../theoretical_models/v2_pt_Helium3_5360GeV_PbPb_020.txt'

df_0_20 = pd.read_csv(input_file_name_0_20, index_col=False, sep= ' ')

print(df_0_20)

pt_0_20 = df_0_20['pt'].tolist()
v2_0_20 = df_0_20['v2'].tolist()
v2_err_0_20 = df_0_20['v2_err'].tolist()

n_points_0_20 = len(v2_err_0_20)

pt_0_20_array = np.array(pt_0_20, dtype=np.float64)
v2_0_20_array = np.array(v2_0_20, dtype=np.float64)
v2_err_0_20_array = np.array(v2_err_0_20, dtype=np.float64)
zeros_0_20_array = np.zeros(n_points_0_20, dtype=np.float64)

gPredWenbin020 = ROOT.TGraphErrors(n_points_0_20, pt_0_20_array, v2_0_20_array, zeros_0_20_array, v2_err_0_20_array)
gPredWenbin020.SetName('gPredWenbin_0_20')
utils.setHistStyle(gPredWenbin020, ROOT.kCyan+2)

# 20 - 40% centrality

input_file_name_20_40 = '../theoretical_models/v2_pt_Helium3_5360GeV_PbPb_2040.txt'

df_20_40 = pd.read_csv(input_file_name_20_40, index_col=False, sep= ' ')
print(df_20_40)

pt_20_40 = df_20_40['pt'].tolist()
v2_20_40 = df_20_40['v2'].tolist()
v2_err_20_40 = df_20_40['v2_err'].tolist()

n_points_20_40 = len(v2_err_20_40)

pt_20_40_array = np.array(pt_20_40, dtype=np.float64)
v2_20_40_array = np.array(v2_20_40, dtype=np.float64)
v2_err_20_40_array = np.array(v2_err_20_40, dtype=np.float64)
zeros_20_40_array = np.zeros(n_points_20_40, dtype=np.float64)

gPredWenbin2040 = ROOT.TGraphErrors(n_points_20_40, pt_20_40_array, v2_20_40_array, zeros_20_40_array, v2_err_20_40_array)
gPredWenbin2040.SetName('gPredWenbin_20_40')
utils.setHistStyle(gPredWenbin2040, ROOT.kViolet+1)

# saving to file
output_file = ROOT.TFile('../theoretical_models/WenbinPredictions.root', 'recreate')
gPredWenbin020.Write()
gPredWenbin2040.Write()



