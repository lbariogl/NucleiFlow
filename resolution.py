import ROOT

import sys
sys.path.append('utils')
import utils as utils

from itertools import combinations

input_file = ROOT.TFile('/data/lbariogl/flow/LHC23zzh_pass2/AnalysisResults.root')

ref_names = ['FT0C', 'FT0A', 'FV0A', 'Bpos', 'Bneg']

hSP_dict = {}
hProfile_dict = {}
cSpProfile_dict = {}
hResolution_dict = {}

def getResolution(det_name1, det_name2, det_name3):

  print(f'det1: {det_name1} det2: {det_name2} det3: {det_name3}')

  # first combination
  det1_det2_name = f'{det_name1}_{det_name2}'
  print(f'det1_det2_name: {det1_det2_name}')
  if det1_det2_name not in hSP_dict.keys():
    hSP_dict[det1_det2_name] = input_file.Get(f'flow-qc/flow/hScalarProduct_{det1_det2_name}')

    hProfile_dict[det1_det2_name] = hSP_dict[det1_det2_name].ProfileX(f'hProfile_{det1_det2_name}')
    utils.setHistStyle(hProfile_dict[det1_det2_name], ROOT.kRed)

    cSpProfile_dict[det1_det2_name] = ROOT.TCanvas(f'cSpProfile_{det1_det2_name}', f'cSpProfile_{det1_det2_name}', 800, 600)
    hSP_dict[det1_det2_name].Draw('colz')
    hProfile_dict[det1_det2_name].Draw('pe same')
  y_title1 = r'#LT ' + hSP_dict[det1_det2_name].GetYaxis().GetTitle() + r' #GT'
  print(f'y_title1: {y_title1}')

  # second combination
  det2_det3_name = f'{det_name2}_{det_name3}'
  print(f'det2_det3_name: {det2_det3_name}')
  if det2_det3_name not in hSP_dict.keys():
    hSP_dict[det2_det3_name] = input_file.Get(f'flow-qc/flow/hScalarProduct_{det2_det3_name}')

    hProfile_dict[det2_det3_name] = hSP_dict[det2_det3_name].ProfileX(f'hProfile_{det2_det3_name}')
    utils.setHistStyle(hProfile_dict[det2_det3_name], ROOT.kRed)

    cSpProfile_dict[det2_det3_name] = ROOT.TCanvas(f'cSpProfile_{det2_det3_name}', f'cSpProfile_{det2_det3_name}', 800, 600)
    hSP_dict[det2_det3_name].Draw('colz')
    hProfile_dict[det2_det3_name].Draw('pe same')
  y_title2 = r'#LT ' + hSP_dict[det2_det3_name].GetYaxis().GetTitle() + r' #GT'
  print(f'y_title2: {y_title2}')

  # third combination
  det1_det3_name = f'{det_name1}_{det_name3}'
  print(f'det1_det3_name: {det1_det3_name}')
  if det1_det3_name not in hSP_dict.keys():
    hSP_dict[det1_det3_name] = input_file.Get(f'flow-qc/flow/hScalarProduct_{det1_det3_name}')

    hProfile_dict[det1_det3_name] = hSP_dict[det1_det3_name].ProfileX(f'hProfile_{det1_det3_name}')
    utils.setHistStyle(hProfile_dict[det1_det3_name], ROOT.kRed)

    cSpProfile_dict[det1_det3_name] = ROOT.TCanvas(f'cSpProfile_{det1_det3_name}', f'cSpProfile_{det1_det3_name}', 800, 600)
    hSP_dict[det1_det3_name].Draw('colz')
    hProfile_dict[det1_det3_name].Draw('pe same')
  y_title3 = r'#LT ' + hSP_dict[det1_det3_name].GetYaxis().GetTitle() + r' #GT'
  print(f'y_title3: {y_title3}')

  # evaluate resolutions

  n_bins = hSP_dict[det1_det2_name].GetNbinsX()
  cent_axis_limits = [hSP_dict[det1_det2_name].GetXaxis().GetBinLowEdge(1), hSP_dict[det1_det2_name].GetXaxis().GetBinUpEdge(n_bins)]
  cent_axis_title = hSP_dict[det1_det2_name].GetXaxis().GetTitle()

  resolution_name1 = f'{det_name1}_{det_name2}_{det_name3}'
  resolution_title1 = f'{y_title1}' + r' #upoint' + f'{y_title2}' + r' / ' + f'{y_title3}'
  hResolution_dict[resolution_name1] = ROOT.TH1F(f'hResolution_{resolution_name1}', f';{cent_axis_title};{resolution_title1}', n_bins, cent_axis_limits[0], cent_axis_limits[1])

  resolution_name2 = f'{det_name2}_{det_name3}_{det_name1}'
  resolution_title2 = f'{y_title2}' + r' #upoint' + f'{y_title3}' + r' / ' + f'{y_title1}'
  hResolution_dict[resolution_name2] = ROOT.TH1F(f'hResolution_{resolution_name2}', f';{cent_axis_title};{resolution_title2}', n_bins, cent_axis_limits[0], cent_axis_limits[1])

  resolution_name3 = f'{det_name3}_{det_name1}_{det_name2}'
  resolution_title3 = f'{y_title3}' + r' #upoint' + f'{y_title1}' + r' / ' + f'{y_title2}'
  hResolution_dict[resolution_name3] = ROOT.TH1F(f'hResolution_{resolution_name3}', f';{cent_axis_title};{resolution_title3}', n_bins, cent_axis_limits[0], cent_axis_limits[1])

det_combos = list(combinations(ref_names, 3))

for det_combo in det_combos:
  getResolution(det_combo[0], det_combo[1], det_combo[2])

output_file = ROOT.TFile('resolution.root', 'recreate')
SP_dir = output_file.mkdir('SP')
Resolution_dir = output_file.mkdir('Resolution')

SP_dir.cd()
for key in hSP_dict.keys():
  hSP_dict[key].Write()
  hProfile_dict[key].Write()
  cSpProfile_dict[key].Write()

Resolution_dir.cd()
for res in hResolution_dict.values():
  res.Write()


