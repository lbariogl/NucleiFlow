import ROOT
from hipe4ml.tree_handler import TreeHandler
import pandas as pd
import yaml
import argparse

import sys
sys.path.append('utils')
import utils as utils
from flow import FlowMaker

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

input_file_name = config['input_file_name']
output_dir_name = config['output_dir_name']
output_file_name = config['output_file_name']

nuclei_tree_name = config['nuclei_tree_name']
ep_tree_name = config['ep_tree_name']

selections = config['selections']

cent_detector_label = config['cent_detector_label']

#BB parameters
p_train = config['p_train']
resolution_train = config['resolution_train']
n_sigma_plot = config['n_sigma_plot']

#create output file
output_file = ROOT.TFile(f'{output_dir_name}/{output_file_name}', 'recreate')


# get a unique df from nuclei and ep trees
nuclei_hdl = TreeHandler(input_file_name, f'{nuclei_tree_name};', folder_name='DF*')
nuclei_df = nuclei_hdl._full_data_frame

nucleiflow_hdl = TreeHandler(input_file_name, f'{ep_tree_name};', folder_name='DF*')
nucleiflow_df = nucleiflow_hdl._full_data_frame

complete_df = pd.concat([nuclei_df, nucleiflow_df], axis=1, join='inner')

# define new columns
utils.redifineColumns(complete_df)

# apply common selections
complete_df.query(selections, inplace=True)

selections = 'fSign < 0 and abs(fEta) < 0.8 and abs(fDCAxy) < 0.1 and fAvgItsClusSize > 4.5 and abs(fRapidity) < 0.5'
complete_df.query(selections, inplace=True)

# Create QC histograms
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

# Fill QC histograms
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

# Saving histograms to file
qc_dir = output_file.mkdir('QC')
qc_dir.cd()
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

# Save histogram as PDF
utils.saveCanvasAsPDF(hEta, f'{output_dir_name}/qc_plots')
utils.saveCanvasAsPDF(hAvgItsClusSize, f'{output_dir_name}/qc_plots')
utils.saveCanvasAsPDF(hAvgItsClusSize, f'{output_dir_name}/qc_plots', is2D=True)
utils.saveCanvasAsPDF(hPhi, f'{output_dir_name}/qc_plots')
utils.saveCanvasAsPDF(hPsiFT0C, f'{output_dir_name}/qc_plots')
utils.saveCanvasAsPDF(hPhiMinusPsiFT0C, f'{output_dir_name}/qc_plots')
utils.saveCanvasAsPDF(hV2, f'{output_dir_name}/qc_plots')
cTPC.SaveAs(f'{output_dir_name}/qc_plots/cTPC.pdf')

# Flow measurement
output_file.cd()

# apply common selections
complete_df.query(selections, inplace=True)

# 0-10%
flow_maker_0_10 = FlowMaker()
flow_maker_0_10.data_df = complete_df.query(f'fCentFT0C > {0} and fCentFT0C < {10}')
flow_maker_0_10.pt_bins = [2, 2.4, 2.8, 3.2, 4., 4.8, 5.6, 6.4, 7.2]
flow_maker_0_10.cent_limits = [0, 10]
flow_maker_0_10.output_file = output_file
flow_maker_0_10.plot_dir = '../results/plots/cent_0_10'
flow_maker_0_10.color = ROOT.kRed+1

flow_maker_0_10.make_flow()
flow_maker_0_10.dump_to_output_file()
flow_maker_0_10.dump_to_pdf()

# 10-20%
flow_maker_10_20 = FlowMaker()
flow_maker_10_20.data_df = complete_df.query(f'fCentFT0C > {10} and fCentFT0C < {20}')
flow_maker_10_20.pt_bins = [2, 2.4, 2.8, 3.2, 4., 4.8, 5.6, 6.4, 7.2]
flow_maker_10_20.cent_limits = [10, 20]
flow_maker_10_20.output_file = output_file
flow_maker_10_20.plot_dir = '../results/plots/cent_10_20'
flow_maker_10_20.color = ROOT.kOrange-3
flow_maker_10_20.make_flow()
flow_maker_10_20.dump_to_output_file()
flow_maker_10_20.dump_to_pdf()

# 20-40%
flow_maker_20_40 = FlowMaker()
flow_maker_20_40.data_df = complete_df.query(f'fCentFT0C > {20} and fCentFT0C < {40}')
flow_maker_20_40.pt_bins = [2, 2.4, 2.8, 3.2, 4., 4.8, 5.6, 6.4, 7.2]
flow_maker_20_40.cent_limits = [20, 40]
flow_maker_20_40.output_file = output_file
flow_maker_20_40.plot_dir = '../results/plots/cent_20_40'
flow_maker_20_40.color = ROOT.kGreen+2

flow_maker_20_40.make_flow()
flow_maker_20_40.dump_to_output_file()
flow_maker_20_40.dump_to_pdf()

# 40-60%
flow_maker_40_60 = FlowMaker()
flow_maker_40_60.data_df = complete_df.query(f'fCentFT0C > {40} and fCentFT0C < {60}')
flow_maker_40_60.pt_bins = [2., 2.4, 2.8, 3.2, 4., 4.8, 6.]
flow_maker_40_60.cent_limits = [40, 60]
flow_maker_40_60.output_file = output_file
flow_maker_40_60.plot_dir = '../results/plots/cent_40_60'
flow_maker_40_60.color = ROOT.kAzure+2

flow_maker_40_60.make_flow()
flow_maker_40_60.dump_to_output_file()
flow_maker_40_60.dump_to_pdf()

# 60-80%
flow_maker_60_80 = FlowMaker()
flow_maker_60_80.data_df = complete_df.query(f'fCentFT0C > {60} and fCentFT0C < {80}')
flow_maker_60_80.pt_bins = [2., 2.4, 2.8, 3.2, 4., 4.8, 6.]
flow_maker_60_80.cent_limits = [60, 80]
flow_maker_60_80.output_file = output_file
flow_maker_60_80.plot_dir = '../results/plots/cent_60_80'
flow_maker_60_80.color = ROOT.kViolet+5

flow_maker_60_80.make_flow()
flow_maker_60_80.dump_to_output_file()
flow_maker_60_80.dump_to_pdf()

resolution_file = ROOT.TFile('Resolution_FT0C.root')
hResolution = resolution_file.Get('Resolutuion')
hResolution.SetDirectory(0)

res_0_10 = hResolution.GetBinContent(1)
res_10_20 = hResolution.GetBinContent(2)
res_20_40 = (hResolution.GetBinContent(3) + hResolution.GetBinContent(4)) / 2
res_40_60 = (hResolution.GetBinContent(5) + hResolution.GetBinContent(6)) / 2
res_60_80 = (hResolution.GetBinContent(7) + hResolution.GetBinContent(8)) / 2

resolution = [res_0_10, res_10_20, res_20_40, res_40_60, res_60_80]
uncorr_v2 = [flow_maker_0_10.hV2vsPt, flow_maker_10_20.hV2vsPt, flow_maker_20_40.hV2vsPt, flow_maker_40_60.hV2vsPt, flow_maker_60_80.hV2vsPt]
labels = [r'0 - 10%', r'10 - 20%', r'20 - 40%', r'40 - 60%', r'60 - 80%']

cV2 = ROOT.TCanvas('cV2', 'cV2', 800, 600)
frame = cV2.DrawFrame(1.7, -0.2, 9, 1., r';#it{p}_{T} (GeV/#it{c}); v_{2}')
cV2.SetBottomMargin(0.13)
cV2.cd()
legend = ROOT.TLegend(0.61, 0.58, 0.87,0.81, 'FT0C centrality', 'brNDC')
legend.SetBorderSize(0)
legend.SetNColumns(2)

for i, res in enumerate(resolution):
  print(f'{labels[i]} -> resolution: {res}')
  uncorr_v2[i].Scale(1/res)
  legend.AddEntry(uncorr_v2[i], labels[i], 'PF')
  uncorr_v2[i].Draw('PE SAME')

legend.Draw()

output_file.cd()
cV2.Write()
