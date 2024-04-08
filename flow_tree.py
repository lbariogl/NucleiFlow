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
legend = ROOT.TLegend(0.53, 0.62, 0.87, 0.87, 'FT0C centrality', 'brNDC')
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
