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
cent_limits = [30, 50]

nuclei_hdl = TreeHandler(input_file_name, 'O2nucleitable', folder_name='DF')
nuclei_df = nuclei_hdl._full_data_frame

nucleiflow_hdl = TreeHandler(input_file_name, 'O2nucleitableflow', folder_name='DF')
nucleiflow_df = nucleiflow_hdl._full_data_frame

complete_df = pd.concat([nuclei_df, nucleiflow_df], axis=1, join='inner')
utils.redifineColumns(complete_df)

# select centrality
complete_df.query(f'fCentFT0C > {cent_limits[0]} and fCentFT0C < {cent_limits[1]}', inplace=True)

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
flow_maker = FlowMaker()
flow_maker.data_df = complete_df
flow_maker.pt_bins = [1, 1.2, 1.4, 1.6, 2., 2.4, 2.8, 3.2, 3.6, 4]
flow_maker.output_dir = output_file

flow_maker.make_flow()
flow_maker.dump_to_output_dir()

