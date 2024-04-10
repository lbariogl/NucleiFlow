import ROOT
from hipe4ml.tree_handler import TreeHandler
import pandas as pd
import yaml
import argparse
import numpy as np
from itertools import product

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

centrality_classes = config['centrality_classes']
pt_bins = config['pt_bins']
cent_colours = config['cent_colours']

do_syst = config['do_syst']

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

# Get resolution from file
resolution_file = ROOT.TFile('Resolution_FT0C.root')
hResolution = resolution_file.Get('Resolutuion')
hResolution.SetDirectory(0)

res_0_10 = hResolution.GetBinContent(1)
res_10_20 = hResolution.GetBinContent(2)
res_20_40 = (hResolution.GetBinContent(3) + hResolution.GetBinContent(4)) / 2
res_40_60 = (hResolution.GetBinContent(5) + hResolution.GetBinContent(6)) / 2
res_60_80 = (hResolution.GetBinContent(7) + hResolution.GetBinContent(8)) / 2

resolutions = [res_0_10, res_10_20, res_20_40, res_40_60, res_60_80]

# Flow measurement
output_file.cd()

n_cent_classes = len(centrality_classes)

flow_makers = []

for i_cent in range(n_cent_classes):
  flow_maker = FlowMaker()
  flow_maker.data_df = complete_df
  flow_maker.selection_string = selections
  flow_maker.pt_bins = pt_bins[i_cent]
  flow_maker.cent_limits = centrality_classes[i_cent]
  flow_maker.resolution = resolutions[i_cent]

  flow_maker.output_file = output_file
  flow_maker.plot_dir = f'../results/plots/cent_{flow_maker.cent_limits[0]}_{flow_maker.cent_limits[1]}'
  flow_maker.color = cent_colours[i_cent]

  flow_maker.make_flow()
  flow_maker.dump_to_output_file()
  flow_maker.dump_to_pdf()

  flow_makers.append(flow_maker)


cV2 = ROOT.TCanvas('cV2', 'cV2', 800, 600)
frame = cV2.DrawFrame(1.7, -0.2, 9, 1., r';#it{p}_{T} (GeV/#it{c}); v_{2}')
cV2.SetBottomMargin(0.13)
cV2.cd()
legend = ROOT.TLegend(0.61, 0.58, 0.87,0.81, 'FT0C centrality', 'brNDC')
legend.SetBorderSize(0)
legend.SetNColumns(2)

for i_cent in range(n_cent_classes):
  cent_label = f'{flow_makers[i_cent].cent_limits[0]} - {flow_makers[i_cent].cent_limits[1]}' + r'%'
  legend.AddEntry(flow_makers[i_cent].hV2vsPt, cent_label,'PF')
  flow_makers[i_cent].hV2vsPt.Draw('PE SAME')

legend.Draw()

output_file.cd()
cV2.Write()
