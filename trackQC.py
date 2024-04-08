import ROOT
from hipe4ml.tree_handler import TreeHandler
import numpy as np
import pandas as pd

import os
import sys
sys.path.append('utils')
import utils as utils


input_file_name = '/data/lbariogl/flow/LHC23_PbPb_pass2_new/AO2D_merged.root'
output_file_name = 'track_qc.root'

output_file = ROOT.TFile(output_file_name, 'recreate')

if not os.path.exists('qc_plots'):
  os.makedirs('qc_plots')

plots_dir = 'qc_plots'

cent_detector_label = 'FT0C'
cent_limits = [30, 50]

nuclei_hdl = TreeHandler(input_file_name, 'O2nucleitable;', folder_name='DF*')
nuclei_df = nuclei_hdl._full_data_frame

nucleiflow_hdl = TreeHandler(input_file_name, 'O2nucleitableflow;', folder_name='DF*')
nucleiflow_df = nucleiflow_hdl._full_data_frame

complete_df = pd.concat([nuclei_df, nucleiflow_df], axis=1, join='inner')

utils.redifineColumns(complete_df)

#BB parameters
p_train = [-321.34, 0.6539, 1.591, 0.8225, 2.363]
resolution_train = 0.09
n_sigma = 3

# complete_df.query(f'fCentFT0C > {cent_limits[0]} and fCentFT0C < {cent_limits[1]}', inplace=True)
selections = 'abs(fEta) < 0.8 and abs(fDCAxy) < 0.1 and abs(fRapidity) < 0.5 and abs(fNsigmaTPC3He) < 3'
complete_df.query(selections, inplace=True)

hEta = ROOT.TH1F('hEta', ';#eta;', 200, -1., 1.)
utils.setHistStyle(hEta, ROOT.kRed+2)
hAvgItsClusSize = ROOT.TH1F('hAvgItsClusSize', r';#LT ITS cluster size #GT;', 20, 0, 20)
utils.setHistStyle(hAvgItsClusSize, ROOT.kRed+2)
hTPCsignalVsPoverZ = ROOT.TH2F('hTPCsignalVsPoverZ', r';#it{p}/z (GeV/#it{c}); d#it{E} / d#it{x} (a.u.)', 600, -6., 6., 1400, 0., 1400.)
hPhi = ROOT.TH1F('hPhi', r';#phi;', 140, -7., 7.)
utils.setHistStyle(hPhi, ROOT.kRed+2)
hPsiFT0C = ROOT.TH1F('hPsiFT0C', r';#Psi_{FTOC};', 140, -7., 7.)
utils.setHistStyle(hPsiFT0C, ROOT.kRed+2)
hPhiMinusPsiFT0C = ROOT.TH1F('hPhiMinusPsiFT0C', r';phi - #Psi_{FTOC};', 140, -7., 7.)
utils.setHistStyle(hPhiMinusPsiFT0C, ROOT.kRed+2)
hV2 = ROOT.TH1F('hV2', r';cos(2(#phi - #Psi_{FTOC}))', 100, -1, 1)
utils.setHistStyle(hV2, ROOT.kRed+2)

print('Filling eta')
for eta in complete_df['fEta']:
  hEta.Fill(eta)

print('Filling cluster size')
for avgClus in complete_df['fAvgItsClusSize']:
  hAvgItsClusSize.Fill(avgClus)

print('Filling specific energy loss')
for rig, sign, signal in zip(complete_df['fTPCInnerParam'], complete_df['fSign'], complete_df['fTPCsignal']):
  hTPCsignalVsPoverZ.Fill(sign*rig, signal)

print('Filling phi')
for phi, psi_ft0c in zip(complete_df['fPhi'], complete_df['fPsiFT0C']):
  hPhi.Fill(phi)
  hPsiFT0C.Fill(psi_ft0c)
  delta_phi = phi - psi_ft0c
  hPhiMinusPsiFT0C.Fill(delta_phi)
  hV2.Fill(ROOT.TMath.Cos(2*delta_phi))


functions = utils.getBBAfunctions(parameters=p_train, resolution=resolution_train)

cTPC = ROOT.TCanvas('cTPC', 'cTPC', 800, 600)
hTPCsignalVsPoverZ.Draw('colz')

for f in functions:
  f.Draw('L SAME')

output_file.cd()
hEta.Write()
hAvgItsClusSize.Write()
hTPCsignalVsPoverZ.Write()
hPhi.Write()
hPsiFT0C.Write()
hPhiMinusPsiFT0C.Write()
cTPC.Write()
hV2.Write()
for f in functions:
  f.Write()

# saving histograms as PDF
utils.saveCanvasAsPDF(hEta, plots_dir)
utils.saveCanvasAsPDF(hAvgItsClusSize, plots_dir)
utils.saveCanvasAsPDF(hAvgItsClusSize, plots_dir, is2D=True)
utils.saveCanvasAsPDF(hPhi, plots_dir)
utils.saveCanvasAsPDF(hPsiFT0C, plots_dir)
utils.saveCanvasAsPDF(hPhiMinusPsiFT0C, plots_dir)
utils.saveCanvasAsPDF(hV2, plots_dir)
cTPC.SaveAs(f'{plots_dir}/cTPC.pdf')

