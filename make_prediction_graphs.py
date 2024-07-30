import ROOT
import pandas as pd
import numpy as np

import sys
sys.path.append('utils')
import utils as utils

# 0 - 20% centrality

input_file_name_0_20 = '../theoretical_models/v2_pt_Helium3_5360GeV_PbPb_020.txt'

df_0_20 = pd.read_csv(input_file_name_0_20, index_col=False, sep= ' ')

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
output_file = ROOT.TFile('../theoretical_models/Predictions.root', 'recreate')
gPredWenbin020.Write()
gPredWenbin2040.Write()

# Sun

sun_df_0_20 = pd.read_excel('../theoretical_models/sun_v2_3He_0_20.xlsx')
sun_df_20_40 = pd.read_excel('../theoretical_models/sun_v2_3He_20_40.xlsx')
sun_df_40_60 = pd.read_excel('../theoretical_models/sun_v2_3He_40_60.xlsx')

pt_array_sun = sun_df_0_20['pt'].to_numpy(dtype=np.float64)
n_points_sun = len(pt_array_sun)
zeros_array_sun = np.zeros(n_points_sun, dtype=np.float64)

v2_0_20_array_sun = sun_df_0_20['v2'].to_numpy(dtype=np.float64)
v2_err_0_20_array_sun = sun_df_0_20['err'].to_numpy(dtype=np.float64)

gPredSun020 = ROOT.TGraphErrors(n_points_sun, pt_array_sun, v2_0_20_array_sun, zeros_array_sun, v2_err_0_20_array_sun)
gPredSun020.SetName('gPredSun_0_20')
utils.setHistStyle(gPredSun020, ROOT.kSpring+5)

v2_20_40_array_sun = sun_df_20_40['v2'].to_numpy(dtype=np.float64)
v2_err_20_40_array_sun = sun_df_20_40['err'].to_numpy(dtype=np.float64)

gPredSun2040 = ROOT.TGraphErrors(n_points_sun, pt_array_sun, v2_20_40_array_sun, zeros_array_sun, v2_err_20_40_array_sun)
gPredSun2040.SetName('gPredSun_20_40')
utils.setHistStyle(gPredSun2040, ROOT.kOrange-3)

v2_40_60_array_sun = sun_df_40_60['v2'].to_numpy(dtype=np.float64)
v2_err_40_60_array_sun = sun_df_40_60['err'].to_numpy(dtype=np.float64)

gPredSun4060 = ROOT.TGraphErrors(n_points_sun, pt_array_sun, v2_40_60_array_sun, zeros_array_sun, v2_err_40_60_array_sun)
gPredSun4060.SetName('gPredSun_40_60')
utils.setHistStyle(gPredSun4060, ROOT.kRed+1)

# saving to file
gPredSun020.Write()
gPredSun2040.Write()
gPredSun4060.Write()
