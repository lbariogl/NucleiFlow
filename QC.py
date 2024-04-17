import ROOT
from hipe4ml.tree_handler import TreeHandler
import pandas as pd
import yaml
import argparse

import sys
sys.path.append('utils')
import utils as utils

parser = argparse.ArgumentParser(
    description='Configure the parameters of the script.')
parser.add_argument('--config-file', dest='config_file',
                    help="path to the YAML file with configuration.", default='config/config_analysis.yaml')
args = parser.parse_args()
if args.config_file == "":
    print('** No config file provided. Exiting. **')
    exit()

config_file = open(args.config_file, 'r')
config = yaml.full_load(config_file)

input_file_name = config['input_file_name']
input_file_AR_name = config['input_file_AR_name']
output_dir_name = config['output_dir_name']
output_file_name = config['output_file_name_qc']

nuclei_tree_name = config['nuclei_tree_name']
ep_tree_name = config['ep_tree_name']

mandatory_selections = config['mandatory_selections']
selection_dict = config['selection_dict']
selection_list = selection_dict.values()
selections = " and ".join(selection_list)

cent_detector_label = config['cent_detector_label']

centrality_classes = config['centrality_classes']
pt_bins = config['pt_bins']
cent_colours = config['cent_colours']

do_syst = config['do_syst']

#BB parameters
p_train = config['p_train']
resolution_train = config['resolution_train']
n_sigma_plot = config['n_sigma_plot']

# create output file
output_file = ROOT.TFile(f'{output_dir_name}/{output_file_name}', 'recreate')


# get a unique df from nuclei and ep trees
nuclei_hdl = TreeHandler(input_file_name, f'{nuclei_tree_name};', folder_name='DF*')
nuclei_df = nuclei_hdl._full_data_frame

nucleiflow_hdl = TreeHandler(input_file_name, f'{ep_tree_name};', folder_name='DF*')
nucleiflow_df = nucleiflow_hdl._full_data_frame

complete_df = pd.concat([nuclei_df, nucleiflow_df], axis=1, join='inner')

# set input files for EP qc
input_file_AR = ROOT.TFile(input_file_AR_name)
input_dir_AR = input_file_AR.Get('flow-qc/flow')

# define new columns
utils.redifineColumns(complete_df)

# apply common selections
complete_df.query(f'{mandatory_selections} and {selections}', inplace=True)

# Create QC histograms
hEta = ROOT.TH1F('hEta', ';#eta;', 200, -1., 1.)
utils.setHistStyle(hEta, ROOT.kRed+2)
hAvgItsClusSize = ROOT.TH1F('hAvgItsClusSize', r';#LT ITS cluster size #GT;', 20, 0, 20)
utils.setHistStyle(hAvgItsClusSize, ROOT.kRed+2)

hAvgItsClusSizeCosLambda = ROOT.TH1F('hAvgItsClusSizeCosLambda', r';#LT ITS cluster size #GT #times Cos(#lambda);', 20, 0, 20)
utils.setHistStyle(hAvgItsClusSizeCosLambda, ROOT.kRed+2)
hTPCsignalVsPoverZ = ROOT.TH2F('hTPCsignalVsPoverZ', r';#it{p}/z (GeV/#it{c}); d#it{E} / d#it{x} (a.u.)', 600, -6., 6., 1400, 0., 1400.)
hPhi = ROOT.TH1F('hPhi', r';#phi;', 140, -7., 7.)
utils.setHistStyle(hPhi, ROOT.kRed+2)
hPsiFT0C = ROOT.TH1F('hPsiFT0C', r';#Psi_{FTOC};', 140, -7., 7.)
utils.setHistStyle(hPsiFT0C, ROOT.kRed+2)
hPhiMinusPsiFT0C = ROOT.TH1F('hPhiMinusPsiFT0C', r';phi - #Psi_{FTOC};', 140, -7., 7.)
utils.setHistStyle(hPhiMinusPsiFT0C, ROOT.kRed+2)
hV2 = ROOT.TH1F('hV2', r';cos(2(#phi - #Psi_{FTOC}))', 100, -1, 1)
utils.setHistStyle(hV2, ROOT.kRed+2)

# get histograms for EP qc
hZvtx = input_dir_AR.Get('hRecVtxZData')
utils.setHistStyle(hZvtx, ROOT.kRed+2)
hCentFT0C = input_dir_AR.Get('hCentFT0C')
utils.setHistStyle(hCentFT0C, ROOT.kRed+2)
hDeltaPsi_FT0A_TPCl = input_dir_AR.Get('hDeltaPsi_FT0A_TPCl')
hDeltaPsi_FT0A_TPCr = input_dir_AR.Get('hDeltaPsi_FT0A_TPCr')
hDeltaPsi_FT0C_TPCl = input_dir_AR.Get('hDeltaPsi_FT0C_TPCl')
hDeltaPsi_FT0C_TPCr = input_dir_AR.Get('hDeltaPsi_FT0C_TPCr')
hDeltaPsi_FT0C_FT0A = input_dir_AR.Get('hDeltaPsi_FT0C_FT0A')
hDeltaPsi_TPCl_TPCr = input_dir_AR.Get('hDeltaPsi_TPCl_TPCr')

# Fill QC histograms
print('Filling eta')
for eta in complete_df['fEta']:
  hEta.Fill(eta)

print('Filling cluster size')
for avgClus in complete_df['fAvgItsClusSize']:
  hAvgItsClusSize.Fill(avgClus)

print('Filling cluster size * cos(lambda)')
for avgClus in complete_df['fAvgItsClusSizeCosLambda']:
  hAvgItsClusSizeCosLambda.Fill(avgClus)

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

# Saving histograms to file
qc_dir = output_file.mkdir('QC')
qc_dir.cd()
hEta.Write()
hAvgItsClusSize.Write()
hAvgItsClusSizeCosLambda.Write()
hTPCsignalVsPoverZ.Write()
hPhi.Write()
hPsiFT0C.Write()
hPhiMinusPsiFT0C.Write()
cTPC.Write()
hV2.Write()
for f in functions:
  f.Write()

# Save QC histogram as PDF
utils.saveCanvasAsPDF(hEta, f'{output_dir_name}/qc_plots')
utils.saveCanvasAsPDF(hAvgItsClusSize, f'{output_dir_name}/qc_plots')
utils.saveCanvasAsPDF(hAvgItsClusSizeCosLambda, f'{output_dir_name}/qc_plots')
utils.saveCanvasAsPDF(hTPCsignalVsPoverZ, f'{output_dir_name}/qc_plots', is2D=True)
utils.saveCanvasAsPDF(hPhi, f'{output_dir_name}/qc_plots')
utils.saveCanvasAsPDF(hPsiFT0C, f'{output_dir_name}/qc_plots')
utils.saveCanvasAsPDF(hPhiMinusPsiFT0C, f'{output_dir_name}/qc_plots')
utils.saveCanvasAsPDF(hV2, f'{output_dir_name}/qc_plots')
cTPC.SaveAs(f'{output_dir_name}/qc_plots/cTPC.pdf')
utils.saveCanvasAsPDF(hZvtx, f'{output_dir_name}/qc_plots')
utils.saveCanvasAsPDF(hCentFT0C, f'{output_dir_name}/qc_plots')
utils.saveCanvasAsPDF(hDeltaPsi_FT0A_TPCl, f'{output_dir_name}/qc_plots', is2D=True)
utils.saveCanvasAsPDF(hDeltaPsi_FT0A_TPCr, f'{output_dir_name}/qc_plots', is2D=True)
utils.saveCanvasAsPDF(hDeltaPsi_FT0C_TPCl, f'{output_dir_name}/qc_plots', is2D=True)
utils.saveCanvasAsPDF(hDeltaPsi_FT0C_TPCr, f'{output_dir_name}/qc_plots', is2D=True)
utils.saveCanvasAsPDF(hDeltaPsi_FT0C_FT0A, f'{output_dir_name}/qc_plots', is2D=True)
utils.saveCanvasAsPDF(hDeltaPsi_TPCl_TPCr, f'{output_dir_name}/qc_plots', is2D=True)

# resolution plot
resolution_file = ROOT.TFile('Resolution_FT0C.root')
hResolutionFT0C = resolution_file.Get('Resolutuion')
hResolutionFT0C.SetDirectory(0)
hResolutionFT0C.SetName('hResolutionFT0C')
hResolutionFT0C.SetTitle(r';FT0C percentile (%); R_{2}')
utils.setHistStyle(hResolutionFT0C, ROOT.kRed+2)
cResolutionFT0C = ROOT.TCanvas('cResolutionFT0C', 'cResolutionFT0C', 800, 600)
cResolutionFT0C.SetBottomMargin(0.15)
hResolutionFT0C.Draw('PE')
cResolutionFT0C.SaveAs(f'{output_dir_name}/qc_plots/cResolutionFT0C.pdf')
