from itertools import combinations
import yaml
import argparse
import numpy as np
import ROOT

import sys
sys.path.append('utils')
import utils as utils

parser = argparse.ArgumentParser(
    description='Configure the parameters of the script.')
parser.add_argument('--config-file', dest='config_file',
                    help="path to the YAML file with configuration.", default='')
args = parser.parse_args()
if args.config_file == "":
    print('** No config file provided. Exiting. **')
    exit()

config_file = open(args.config_file, 'r')
config = yaml.full_load(config_file)

input_file_AR_name = config['input_file_AR_name']
input_file = ROOT.TFile(input_file_AR_name)

ref_names = ['FT0C', 'FT0A', 'TPCl', 'TPCr']

hSP_dict = {}
hProfile_dict = {}
cSpProfile_dict = {}
hResolution_dict = {}

cent_bins = np.array([0., 10., 20., 30., 40., 50., 60.,
                     70., 80., 90., 100.], dtype=np.float64)
n_cent_bins = len(cent_bins) - 1


def getResolution(histo_resolution, det1_det2_name, det2_det3_name, det1_det3_name):
    histo_resolution.GetYaxis().SetRangeUser(0., 1.)
    n_bins = histo_resolution.GetNbinsX()
    for ibin in range(1, n_bins+1):
        val1 = hProfile_dict[det1_det2_name].GetBinContent(ibin)
        val2 = hProfile_dict[det2_det3_name].GetBinContent(ibin)
        val3 = hProfile_dict[det1_det3_name].GetBinContent(ibin)

        err1 = hProfile_dict[det1_det2_name].GetBinError(ibin)
        err2 = hProfile_dict[det2_det3_name].GetBinError(ibin)
        err3 = hProfile_dict[det1_det3_name].GetBinError(ibin)

        if val2 > 0:
            val = val1 * val3 / val2
        else:
            val = -999

        if val2 > 0 and val1 > 0 and val2 > 0:
            err = np.sqrt(err1*err1/val1/val1 + err2*err2 /
                          val2/val2 + err3*err3/val3/val3) * val
        else:
            err = 0
        if val > 0:
            histo_resolution.SetBinContent(ibin, np.sqrt(val))
            histo_resolution.SetBinError(ibin, np.sqrt(err))


def doAllPlots(det_name1, det_name2, det_name3):

    # first combination
    det1_det2_name = f'{det_name1}_{det_name2}'
    if det1_det2_name not in hSP_dict.keys():
        hSP_dict[det1_det2_name] = input_file.Get(
            f'flow-qc/flow_ep/hScalarProduct_{det1_det2_name}_EP')
        hSP_dict[det1_det2_name].RebinX(10)

        hProfile_dict[det1_det2_name] = hSP_dict[det1_det2_name].ProfileX(
            f'hProfile_{det1_det2_name}')
        utils.setHistStyle(hProfile_dict[det1_det2_name], ROOT.kRed)

        cSpProfile_dict[det1_det2_name] = ROOT.TCanvas(
            f'cSpProfile_{det1_det2_name}', f'cSpProfile_{det1_det2_name}', 800, 600)
        hSP_dict[det1_det2_name].Draw('colz')
        hProfile_dict[det1_det2_name].Draw('pe same')
    y_title1 = r'#LT ' + \
        hSP_dict[det1_det2_name].GetYaxis().GetTitle() + r' #GT'

    # second combination
    det2_det3_name = f'{det_name2}_{det_name3}'
    if det2_det3_name not in hSP_dict.keys():
        hSP_dict[det2_det3_name] = input_file.Get(
            f'flow-qc/flow_ep/hScalarProduct_{det2_det3_name}_EP')
        hSP_dict[det2_det3_name].RebinX(10)

        hProfile_dict[det2_det3_name] = hSP_dict[det2_det3_name].ProfileX(
            f'hProfile_{det2_det3_name}')
        utils.setHistStyle(hProfile_dict[det2_det3_name], ROOT.kRed)

        cSpProfile_dict[det2_det3_name] = ROOT.TCanvas(
            f'cSpProfile_{det2_det3_name}', f'cSpProfile_{det2_det3_name}', 800, 600)
        hSP_dict[det2_det3_name].Draw('colz')
        hProfile_dict[det2_det3_name].Draw('pe same')
    y_title2 = r'#LT ' + \
        hSP_dict[det2_det3_name].GetYaxis().GetTitle() + r' #GT'

    # third combination
    det1_det3_name = f'{det_name1}_{det_name3}'
    if det1_det3_name not in hSP_dict.keys():
        hSP_dict[det1_det3_name] = input_file.Get(
            f'flow-qc/flow_ep/hScalarProduct_{det1_det3_name}_EP')
        hSP_dict[det1_det3_name].RebinX(10)

        hProfile_dict[det1_det3_name] = hSP_dict[det1_det3_name].ProfileX(
            f'hProfile_{det1_det3_name}')
        utils.setHistStyle(hProfile_dict[det1_det3_name], ROOT.kRed)

        cSpProfile_dict[det1_det3_name] = ROOT.TCanvas(
            f'cSpProfile_{det1_det3_name}', f'cSpProfile_{det1_det3_name}', 800, 600)
        hSP_dict[det1_det3_name].Draw('colz')
        hProfile_dict[det1_det3_name].Draw('pe same')
    y_title3 = r'#LT ' + \
        hSP_dict[det1_det3_name].GetYaxis().GetTitle() + r' #GT'

    # evaluate resolutions
    n_bins = hSP_dict[det1_det2_name].GetNbinsX()
    cent_axis_limits = [hSP_dict[det1_det2_name].GetXaxis().GetBinLowEdge(
        1), hSP_dict[det1_det2_name].GetXaxis().GetBinUpEdge(n_bins)]
    cent_axis_title = hSP_dict[det1_det2_name].GetXaxis().GetTitle()

    resolution_name1 = f'{det_name1}_{det_name2}_{det_name3}'
    resolution_title1 = r'[' + f'{y_title1}' + r' #upoint' + \
        f'{y_title2}' + r' / ' + f'{y_title3}' + r']^{1/2}'
    hResolution_dict[resolution_name1] = ROOT.TH1F(
        f'hResolution_{resolution_name1}', f';{cent_axis_title};{resolution_title1}', n_bins, cent_axis_limits[0], cent_axis_limits[1])
    getResolution(hResolution_dict[resolution_name1],
                  det1_det2_name, det2_det3_name, det1_det3_name)
    utils.setHistStyle(hResolution_dict[resolution_name1], ROOT.kRed)

    resolution_name2 = f'{det_name2}_{det_name3}_{det_name1}'
    resolution_title2 = r'[' + f'{y_title2}' + r' #upoint' + \
        f'{y_title3}' + r' / ' + f'{y_title1}' + r']^{1/2}'
    hResolution_dict[resolution_name2] = ROOT.TH1F(
        f'hResolution_{resolution_name2}', f';{cent_axis_title};{resolution_title2}', n_bins, cent_axis_limits[0], cent_axis_limits[1])
    getResolution(hResolution_dict[resolution_name2],
                  det2_det3_name, det1_det3_name, det1_det2_name)
    utils.setHistStyle(hResolution_dict[resolution_name2], ROOT.kRed)

    resolution_name3 = f'{det_name3}_{det_name1}_{det_name2}'
    resolution_title3 = r'[' + f'{y_title3}' + r' #upoint' + \
        f'{y_title1}' + r' / ' + f'{y_title2}' + r']^{1/2}'
    hResolution_dict[resolution_name3] = ROOT.TH1F(
        f'hResolution_{resolution_name3}', f';{cent_axis_title};{resolution_title3}', n_bins, cent_axis_limits[0], cent_axis_limits[1])
    getResolution(hResolution_dict[resolution_name3],
                  det1_det2_name, det1_det3_name, det2_det3_name)
    utils.setHistStyle(hResolution_dict[resolution_name3], ROOT.kRed)

# rebin histograms


det_combos = list(combinations(ref_names, 3))

for det_combo in det_combos:
    doAllPlots(det_combo[0], det_combo[1], det_combo[2])

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
