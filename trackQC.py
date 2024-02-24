import ROOT
from hipe4ml.tree_handler import TreeHandler
import numpy as np
import struct

import sys
sys.path.append('utils')
import utils as utils


input_file_name = '/data/lbariogl/flow/LHC23zzh_pass2/AO2D.root'
input_dir_name = 'nucleiFlow'
output_file_name = 'track_qc.root'

output_file = ROOT.TFile(output_file_name, 'recreate')

cent_detector_label = 'FT0C'
cent_limits = [30, 50]

nuclei_hdl = TreeHandler(input_file_name, 'O2nucleitable', folder_name='DF')
nuclei_df = nuclei_hdl._full_data_frame

def getITSClSize(itsSizeBitmap):
  sum = 0
  nClus = 0
  for i_layer in range(0,7):
    size = (itsSizeBitmap >> i_layer*4) & 15
    if size > 0:
      nClus = nClus+1
      sum += size
  return sum / nClus

getITSClSize_vectorised = np.vectorize(getITSClSize)

nuclei_df.eval('fAvgItsClusSize = @getITSClSize_vectorised(fITSclusterSizes)', inplace=True)
nuclei_df.drop(columns=['fITSclusterSizes'])

#BB parameters
p0_train =- 321.34
p1_train = 0.6539
p2_train = 1.591
p3_train = 0.8225
p4_train = 2.363
resolution_train = 0.09


def BBA(mom, p0=p0_train, p1=p1_train, p2=p2_train, p3=p3_train, p4=p4_train):

  beta = mom/ np.sqrt(1 + mom * mom)

  aa = np.power(beta, p3)
  bb = np.power(1 / mom, p4)
  bb = np.log(p2 + bb)

  return (p1 - aa - bb) * p0 / aa

def getNsigmaTPC(mom, tpc_signal, resolution_perc=resolution_train):
  exp_signal = BBA(mom)
  resolution = exp_signal * resolution_perc
  return (tpc_signal - exp_signal) / resolution

getNsigmaTPC_vectorised = np.vectorize(getNsigmaTPC)
nuclei_df.eval('fNsigmaTPC3He = @getNsigmaTPC_vectorised(fTPCInnerParam, fTPCsignal)', inplace=True)


selections = 'abs(fEta) < 0.8 and abs(fDCAxy) < 0.1'

nuclei_df.query(selections, inplace=True)

hEta = ROOT.TH1F('hEta', ';#eta;', 200, -1., 1.)
hAvgItsClusSize = ROOT.TH1F('hAvgItsClusSize', ';<ITS cluster size>;', 20, 0, 20)
hTPCsignalVsP = ROOT.TH2F('hTPCsignalVsP', ';#it{p} (GeV/#it{c}); d#{E} / d#it{x} (a.u.)', 600, -6., 6., 0, 1400, 0, 1400)

for eta in nuclei_df['fEta']:
  hEta.Fill(eta)

for avgClus in nuclei_df['fAvgItsClusSize']:
  print('avgClus: ' ,
        avgClus)
  hAvgItsClusSize.Fill(avgClus)

for mom, sig in zip(nuclei_df['fTPCInnerParam'], nuclei_df['fTPCsignal']):
  hTPCsignalVsP.Fill(mom, sig)

output_file.cd()
hEta.Write()
hAvgItsClusSize.Write()
hTPCsignalVsP.Write()




