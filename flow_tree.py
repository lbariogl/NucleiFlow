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

mandatory_selections = config['mandatory_selections']
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

# apply mandatory selections
complete_df.query(mandatory_selections, inplace=True)

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

# Standard selections
flow_makers = []
default_values = []
cent_dirs = []

for i_cent in range(n_cent_classes):
  flow_maker = FlowMaker()
  flow_maker.data_df = complete_df
  flow_maker.selection_string = selections
  flow_maker.pt_bins = pt_bins[i_cent]
  flow_maker.cent_limits = centrality_classes[i_cent]
  flow_maker.resolution = resolutions[i_cent]

  # create output_dir
  cent_dir = output_file.mkdir(f'cent_{flow_maker.cent_limits[0]}_{flow_maker.cent_limits[1]}')
  cent_dirs.append(cent_dir)
  default_dir = cent_dir.mkdir('default')

  flow_maker.output_dir = default_dir
  flow_maker.plot_dir = f'../results/plots/cent_{flow_maker.cent_limits[0]}_{flow_maker.cent_limits[1]}'
  flow_maker.color = cent_colours[i_cent]

  flow_maker.make_flow()
  flow_maker.dump_to_output_file()
  flow_maker.dump_to_pdf()

  flow_makers.append(flow_maker)
  default_values.append(flow_maker.getFlowValues())

# Systematic uncertainties
if do_syst:

  print("** Starting systematic variations **")
  n_trials = config['n_trials']
  # output_dir_syst = output_file.mkdir('trials')

  # List of trial strings to be printed to a text file
  trial_strings = []
  print("----------------------------------")
  print("** Starting systematics analysis **")
  print(f'** {n_trials} trials will be tested **')
  print("----------------------------------")

  cut_dict_syst = config['cut_dict_syst']

  # Create all possible variations (values) fot all the variables of interest (key)
  cut_string_dict = {}
  for var in cut_dict_syst:
    var_dict = cut_dict_syst[var]
    cut_greater = var_dict['cut_greater']
    cut_greater_string = " > " if cut_greater else " < "
    cut_list = var_dict['cut_list']
    cut_arr = np.linspace(cut_list[0], cut_list[1], cut_list[2])
    cut_string_dict[var] = []
    for cut in cut_arr:
        cut_string_dict[var].append(
            var + cut_greater_string + str(cut))

  combos = list(product(*list(cut_string_dict.values())))

  if n_trials < len(combos):
    combo_random_indices = np.random.randint(len(combos), size=n_trials)
  else:
    print(
        f"** Warning: n_trials > n_combinations ({n_trials}, {len(combos)}), taking all the possible combinations **")
    combo_random_indices = np.arange(len(combos))
    np.random.shuffle(combo_random_indices)

  # vary selections for each analysis_object, by centrality class
  for i_cent in range(n_cent_classes):

    histo_v2_syst = []
    n_pt_bins = len(pt_bins[i_cent]) - 1
    for i_pt in range(0, n_pt_bins):
      histo_v2_syst.append(ROOT.TH1F(f'hV2syst_cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}_pt{i_pt}',
                                     ';v_{2}', 100, default_values[i_cent][i_pt][0] - 3*default_values[i_cent][i_pt][1], default_values[i_cent][i_pt][0] + 3*default_values[i_cent][i_pt][1]))

    for i_combo, combo_index in enumerate(combo_random_indices):

      flow_maker_syst = FlowMaker()
      flow_maker_syst.data_df = complete_df
      flow_maker_syst.selection_string = selections
      flow_maker_syst.pt_bins = pt_bins[i_cent]
      flow_maker_syst.cent_limits = centrality_classes[i_cent]
      flow_maker_syst.resolution = resolutions[i_cent]

      combo_suffix = f'_combo{i_combo}'
      trial_strings.append("----------------------------------")
      trial_num_string = f'Trial: {i_combo} / {len(combo_random_indices)}'
      trial_strings.append(trial_num_string)
      print(trial_num_string)
      combo = combos[combo_index]
      sel_string = " & ".join(combo)
      trial_strings.append(str(sel_string))

      flow_maker_syst.selection_string = sel_string
      print(f'selections: {flow_maker_syst.selection_string}')
      flow_maker_syst.suffix = combo_suffix

      flow_maker_syst.make_flow()
      flow_values = flow_maker_syst.getFlowValues()

      # write histograms to file
      trial_dir = cent_dirs[i_cent].mkdir(f'trial_{i_combo}')
      flow_maker_syst.output_dir = trial_dir
      flow_maker_syst.dump_to_output_file()

      print(f'cent: {flow_maker_syst.cent_limits[0]} - {flow_maker_syst.cent_limits[1]}')
      for i_pt in range(0, flow_maker_syst.nPtBins()):
        histo_v2_syst[i_pt].Fill(flow_values[i_pt][0])
        print(f'pt_bin: {i_pt} -> v2: {flow_values[i_pt][0]} +- {flow_values[i_pt][1]}')
      print("----------------------------------")

      del flow_maker_syst
      del flow_values

    output_file.cd(f'cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}')
    for i_pt in range(0, n_pt_bins):
      histo_v2_syst[i_pt].Write()

# Final plots
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
