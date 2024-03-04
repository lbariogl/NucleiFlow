import ROOT
from hipe4ml.tree_handler import TreeHandler
import numpy as np
import pandas as pd
import struct

import sys
sys.path.append('utils')
import utils as utils


input_file_name = '/data/lbariogl/HL_data/unmerged/AO2D_out.root'
output_file_name = 'track_qc.root'

output_file = ROOT.TFile(output_file_name, 'recreate')

cent_detector_label = 'FT0C'
cent_limits = [30, 50]

nuclei_hdl = TreeHandler(input_file_name, 'O2nucleitable', folder_name='DF')
nuclei_df = nuclei_hdl._full_data_frame

nucleiflow_hdl = TreeHandler(input_file_name, 'O2nucleitableflow', folder_name='DF')
nucleiflow_df = nucleiflow_hdl._full_data_frame

complete_df = pd.concat([nuclei_df, nucleiflow_df], axis=1, join='inner')

utils.redifineColumns(complete_df)

#BB parameters
p_train = [-321.34, 0.6539, 1.591, 0.8225, 2.363]
resolution_train = 0.09
n_sigma = 5

selections = 'abs(fEta) < 0.8 and abs(fDCAxy) < 0.1 and fAvgItsClusSize > 4.5 and fTrackedAsHe == True'

complete_df.query(selections, inplace=True)

hEta = ROOT.TH1F('hEta', ';#eta;', 200, -1., 1.)
hAvgItsClusSize = ROOT.TH1F('hAvgItsClusSize', ';<ITS cluster size>;', 20, 0, 20)
hTPCsignalVsPoverZ = ROOT.TH2F('hTPCsignalVsPoverZ', ';#it{p}/z (GeV/#it{c}); d#it{E} / d#it{x} (a.u.)', 600, -6., 6., 1400, 0., 1400.)

print('Filling eta')
for eta in complete_df['fEta']:
  hEta.Fill(eta)


print('Filling cluster size')
for avgClus in complete_df['fAvgItsClusSize']:
  hAvgItsClusSize.Fill(avgClus)

print('Filling specific energy loss')
for rig, sig in zip(complete_df['fTPCInnerParam'], complete_df['fTPCsignal']):
  hTPCsignalVsPoverZ.Fill(rig, sig)

print('Done')

functions = utils.getBBAfunctions(parameters=p_train, resolution=resolution_train)

cTPC = ROOT.TCanvas('cTPC', 'cTPC', 800, 600)
hTPCsignalVsPoverZ.Draw('colz')

for f in functions:
  f.Draw('L SAME')

output_file.cd()
hEta.Write()
hAvgItsClusSize.Write()
hTPCsignalVsPoverZ.Write()
cTPC.Write()
for f in functions:
  f.Write()

