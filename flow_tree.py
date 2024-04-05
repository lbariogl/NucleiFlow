import ROOT
from hipe4ml.tree_handler import TreeHandler
import pandas as pd

import sys
sys.path.append('utils')
import utils as utils
from flow import FlowMaker


input_file_name = '/data/lbariogl/flow/LHC23_PbPb_pass2_new/AO2D_merged.root'
output_file_name = 'output.root'
output_file = ROOT.TFile(output_file_name, 'recreate')

cent_detector_label = 'FT0C'

nuclei_hdl = TreeHandler(input_file_name, 'O2nucleitable;', folder_name='DF*')
nuclei_df = nuclei_hdl._full_data_frame

nucleiflow_hdl = TreeHandler(input_file_name, 'O2nucleitableflow;', folder_name='DF*')
nucleiflow_df = nucleiflow_hdl._full_data_frame

complete_df = pd.concat([nuclei_df, nucleiflow_df], axis=1, join='inner')

utils.redifineColumns(complete_df)

# other selections
selections = 'fSign < 0 and abs(fEta) < 0.8 and abs(fDCAxy) < 0.1 and fAvgItsClusSize > 4.5 and fTrackedAsHe == True and abs(fRapidity) < 0.5'
complete_df.query(selections, inplace=True)

# QC plots
hEta = ROOT.TH1F('hEta', ';#eta;', 200, -1., 1.)
hAvgItsClusSize = ROOT.TH1F('hAvgItsClusSize', ';<ITS cluster size>;', 20, 0, 20)
hTPCsignalVsPoverZ = ROOT.TH2F('hTPCsignalVsPoverZ', ';#it{p}/z (GeV/#it{c}); d#it{E} / d#it{x} (a.u.)', 600, -6., 6., 1400, 0., 1400.)
hPhi = ROOT.TH1F('hPhi', ';#phi;', 140, -7., 7.)
hPsiFT0C = ROOT.TH1F('hPsiFT0C', r';#Psi_{FTOC};', 140, -7., 7.)
hPhiMinusPsiFT0C = ROOT.TH1F('hPhiMinusPsiFT0C', r';phi - #Psi_{FTOC};', 140, -7., 7.)
hV2 = ROOT.TH1F('hV2', r';cos(2(#phi - #Psi_{FTOC}))', 100, -1, 1)

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

print('Filling v2')
for phi, psi_ft0c, v2 in zip(complete_df['fPhi'], complete_df['fPsiFT0C'], complete_df['fV2FT0C']):
  hPhi.Fill(phi)
  hPsiFT0C.Fill(psi_ft0c)
  delta_phi = phi - psi_ft0c
  hPhiMinusPsiFT0C.Fill(delta_phi)
  hV2.Fill(v2)


qc_dir = output_file.mkdir('QC')
qc_dir.cd()
hEta.Write()
hAvgItsClusSize.Write()
hTPCsignalVsPoverZ.Write()
hPhi.Write()
hV2.Write()

hPsiFT0C.Write()
hPhiMinusPsiFT0C.Write()

# Flow measurement
output_file.cd()

# 0-10%
flow_maker_0_10 = FlowMaker()
flow_maker_0_10.data_df = complete_df.query(f'fCentFT0C > {0} and fCentFT0C < {10}')
flow_maker_0_10.pt_bins = [2, 2.4, 2.8, 3.2, 4., 4.8, 5.6, 6.4, 7.2]
flow_maker_0_10.cent_limits = [0, 10]
flow_maker_0_10.output_file = output_file
flow_maker_0_10.plot_dir = 'plots'

flow_maker_0_10.make_flow()
flow_maker_0_10.dump_to_output_file()

# 10-20%
flow_maker_10_20 = FlowMaker()
flow_maker_10_20.data_df = complete_df.query(f'fCentFT0C > {10} and fCentFT0C < {20}')
flow_maker_10_20.pt_bins = [2, 2.4, 2.8, 3.2, 4., 4.8, 5.6, 6.4, 7.2]
flow_maker_10_20.cent_limits = [10, 20]
flow_maker_10_20.output_file = output_file
flow_maker_10_20.plot_dir = 'plots'

flow_maker_10_20.make_flow()
flow_maker_10_20.dump_to_output_file()

# 20-40%
flow_maker_20_40 = FlowMaker()
flow_maker_20_40.data_df = complete_df.query(f'fCentFT0C > {20} and fCentFT0C < {40}')
flow_maker_20_40.pt_bins = [2, 2.4, 2.8, 3.2, 4., 4.8, 5.6, 6.4, 7.2]
flow_maker_20_40.cent_limits = [20, 40]
flow_maker_20_40.output_file = output_file
flow_maker_20_40.plot_dir = 'plots'

flow_maker_20_40.make_flow()
flow_maker_20_40.dump_to_output_file()

# 40-60%
flow_maker_40_60 = FlowMaker()
flow_maker_40_60.data_df = complete_df.query(f'fCentFT0C > {40} and fCentFT0C < {60}')
flow_maker_40_60.pt_bins = [2., 2.4, 2.8, 3.2, 4., 4.8, 6.]
flow_maker_40_60.cent_limits = [40, 60]
flow_maker_40_60.output_file = output_file
flow_maker_40_60.plot_dir = 'plots'

flow_maker_40_60.make_flow()
flow_maker_40_60.dump_to_output_file()

resolution_file = ROOT.TFile('Resolution_FT0C.root')
hResolution = resolution_file.Get('Resolutuion')
hResolution.SetDirectory(0)

res_0_10 = hResolution.GetBinContent(1)
res_10_20 = hResolution.GetBinContent(2)
res_20_40 = (hResolution.GetBinContent(3) + hResolution.GetBinContent(4)) / 2
res_40_60 = (hResolution.GetBinContent(5) + hResolution.GetBinContent(6)) / 2

resolution = [res_0_10, res_10_20, res_20_40, res_40_60]
uncorr_v2 = [flow_maker_0_10.hV2vsPt, flow_maker_10_20.hV2vsPt, flow_maker_20_40.hV2vsPt, flow_maker_40_60.hV2vsPt]
colors = [ROOT.kAzure+2, ROOT.kGreen+2, ROOT.kOrange-3, ROOT.kRed+1]
labels = [r'0 - 10%', r'10 - 20%', r'20 - 40%', r'40 - 60%']

cV2 = ROOT.TCanvas('cV2', 'cV2', 800, 600)
frame = cV2.DrawFrame(1.7, -0.2, 9, 1., r';#it{p}_{T} (GeV/#it{c}); v_{2}')
cV2.SetBottomMargin(0.13)
cV2.cd()
legend = ROOT.TLegend(0.53, 0.67, 0.87, 0.87, 'FT0C centrality', 'brNDC')
legend.SetBorderSize(0)
legend.SetNColumns(2)

for i, res in enumerate(resolution):
  print(f'{labels[i]} -> resolution: {res}')
  uncorr_v2[i].Scale(1/res)
  utils.setHistStyle(uncorr_v2[i], colors[i])
  legend.AddEntry(uncorr_v2[i], labels[i], 'PF')
  uncorr_v2[i].Draw('PE SAME')

legend.Draw()

output_file.cd()
cV2.Write()
