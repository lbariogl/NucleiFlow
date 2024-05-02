import ROOT
from hipe4ml.tree_handler import TreeHandler
import pandas as pd
import yaml
import argparse
import numpy as np
from itertools import product
import os

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
selection_dict = config['selection_dict']
selection_list = selection_dict.values()
selections = " and ".join(selection_list)

ptdep_selection_dict = config['ptdep_selection_dict']
ptdep_selections = ptdep_selection_dict['fAvgItsClusSizeCosLambda']

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
if not os.path.exists(output_dir_name):
  os.makedirs(output_dir_name)
output_file = ROOT.TFile(f'{output_dir_name}/{output_file_name}', 'recreate')
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

resolutions = [res_0_10, res_10_20, res_20_40, res_40_60]

# Flow measurement
output_file.cd()

n_cent_classes = len(centrality_classes)

# Standard selections
flow_makers = []
default_values = []
cent_dirs = []
cent_plots_dir_names = []

for i_cent in range(n_cent_classes):
  flow_maker = FlowMaker()
  flow_maker.data_df = complete_df
  flow_maker.selection_string = selections
  flow_maker.ptdep_selection_string = ptdep_selections
  flow_maker.pt_bins = pt_bins[i_cent]
  flow_maker.cent_limits = centrality_classes[i_cent]
  flow_maker.resolution = resolutions[i_cent]
  flow_maker.print_frame = True

  # create output_dir
  cent_dir = output_file.mkdir(f'cent_{flow_maker.cent_limits[0]}_{flow_maker.cent_limits[1]}')
  cent_dirs.append(cent_dir)
  default_dir = cent_dir.mkdir('default')

  flow_maker.output_dir = default_dir

  plot_dir_name = f'{output_dir_name}/plots/cent_{flow_maker.cent_limits[0]}_{flow_maker.cent_limits[1]}'
  if not os.path.exists(plot_dir_name):
        os.makedirs(plot_dir_name)


  cent_plots_dir_names.append(plot_dir_name)

  flow_maker.plot_dir = plot_dir_name
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
  ptdep_cut_dict_syst = config['ptdep_cut_dict_syst']

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

  ptdep_cut_string_list = []
  for i_pt, pt_el in enumerate(ptdep_cut_dict_syst["fAvgItsClusSizeCosLambda"]):
    ptdep_cut_greater = pt_el['cut_greater']
    ptdep_cut_greater_string = " > " if ptdep_cut_greater else " < "
    ptdep_cut_list = pt_el['cut_list']
    ptdep_cut_arr = np.linspace(ptdep_cut_list[0], ptdep_cut_list[1], ptdep_cut_list[2])
    ptdep_cut_string_list.append([])
    for ptdep_cut in ptdep_cut_arr:
      ptdep_cut_string_list[i_pt].append(
        "fAvgItsClusSizeCosLambda" + ptdep_cut_greater_string + str(ptdep_cut))

  n_ptdep_variations = len(ptdep_cut_string_list[0])
  n_max_pt_bins = len(pt_bins[0]) - 1

  sorted_ptdep_cut_string_list = []

  for i_var in range(0, n_ptdep_variations):
    pt_var_list = []
    for i_pt in range(0, n_max_pt_bins):
      pt_var_list.append(ptdep_cut_string_list[i_pt][i_var])
    sorted_ptdep_cut_string_list.append(pt_var_list)

  mega_combos = list(product(combos, sorted_ptdep_cut_string_list))

  if n_trials < len(mega_combos):
    mega_combo_random_indices = np.random.randint(len(mega_combos), size=n_trials)
  else:
    print(
        f"** Warning: n_trials > n_combinations ({n_trials}, {len(mega_combos)}), taking all the possible combinations **")
    mega_combo_random_indices = np.arange(len(combos))
    np.random.shuffle(mega_combo_random_indices)

  # vary selections for each analysis_object, by centrality class
  for i_cent in range(n_cent_classes):

    histo_v2_syst = []
    canvas_v2_syst = []
    n_pt_bins = len(pt_bins[i_cent]) - 1
    for i_pt in range(0, n_pt_bins):
      histo_v2_syst.append(ROOT.TH1F(f'hV2syst_cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}_pt{i_pt}',
                                     ';v_{2}', 20, default_values[i_cent][i_pt][0] - 3*default_values[i_cent][i_pt][1], default_values[i_cent][i_pt][0] + 3*default_values[i_cent][i_pt][1]))

    for i_mega_combo, mega_combo_index in enumerate(mega_combo_random_indices):

      flow_maker_syst = FlowMaker()
      flow_maker_syst.data_df = complete_df
      flow_maker_syst.pt_bins = pt_bins[i_cent]
      flow_maker_syst.cent_limits = centrality_classes[i_cent]
      flow_maker_syst.resolution = resolutions[i_cent]

      mega_combo_suffix = f'_combo{i_mega_combo}'
      trial_strings.append("----------------------------------")
      trial_num_string = f'Trial: {i_mega_combo} / {len(mega_combo_random_indices)}'
      trial_strings.append(trial_num_string)
      print(trial_num_string)
      combo = mega_combos[mega_combo_index][0]
      sel_string = " & ".join(combo)
      trial_strings.append(str(sel_string))

      flow_maker_syst.selection_string = sel_string
      flow_maker_syst.ptdep_selection_string = mega_combos[mega_combo_index][1]
      print(f'selections: {flow_maker_syst.selection_string}')
      print(f'pt-dependent selections: {flow_maker_syst.ptdep_selection_string}')
      flow_maker_syst.suffix = mega_combo_suffix

      flow_maker_syst.make_flow()
      flow_values = flow_maker_syst.getFlowValues()

      # write histograms to file
      trial_dir = cent_dirs[i_cent].mkdir(f'trial_{i_mega_combo}')
      flow_maker_syst.output_dir = trial_dir
      flow_maker_syst.dump_to_output_file()

      print(f'cent: {flow_maker_syst.cent_limits[0]} - {flow_maker_syst.cent_limits[1]}')
      for i_pt in range(0, flow_maker_syst.nPtBins()):
        histo_v2_syst[i_pt].Fill(flow_values[i_pt][0])
        # print(f'pt_bin: {i_pt} -> v2: {flow_values[i_pt][0]} +- {flow_values[i_pt][1]}')
      print("----------------------------------")

      del flow_maker_syst
      del flow_values

    output_file.cd(f'cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}')
    for i_pt in range(0, n_pt_bins):

      # create canvas with the central value
      histo_name = histo_v2_syst[i_pt].GetName()
      canvas_syst_name = histo_name.replace('h', 'c', 1)
      canvas_syst = ROOT.TCanvas(canvas_syst_name, canvas_syst_name, 800, 600)
      canvas_syst.SetBottomMargin(0.13)
      canvas_syst.SetLeftMargin(0.13)

      std_line = ROOT.TLine(default_values[i_cent][i_pt][0], 0, default_values[i_cent][i_pt][0], 1.05 * histo_v2_syst[i_pt].GetMaximum())
      std_line.SetLineColor(ROOT.kAzure+2)
      std_line.SetLineWidth(2)
      # create box for statistical uncertainty
      std_errorbox = ROOT.TBox(default_values[i_cent][i_pt][0] - default_values[i_cent][i_pt][1], 0,
                               default_values[i_cent][i_pt][0] + default_values[i_cent][i_pt][1], 1.05 * histo_v2_syst[i_pt].GetMaximum())
      std_errorbox.SetFillColorAlpha(ROOT.kAzure+1, 0.5)
      std_errorbox.SetLineWidth(0)

      info_panel = ROOT.TPaveText(0.6, 0.6, 0.8, 0.82, 'NDC')
      info_panel.SetBorderSize(0)
      info_panel.SetFillStyle(0)
      info_panel.SetTextAlign(12)
      info_panel.SetTextFont(42)
      info_panel.AddText(r'PbPb, #sqrt{#it{s}_{nn}} = 5.36 TeV')
      info_panel.AddText(f'{centrality_classes[i_cent][0]} - {centrality_classes[i_cent][1]} % {cent_detector_label}')
      pt_label = f'{pt_bins[i_cent][i_pt]:.1f}' + r' #leq #it{p}_{T} < ' + \
          f'{pt_bins[i_cent][i_pt+1]:.1f}' + r' GeV/#it{c}'
      info_panel.AddText(pt_label)

      canvas_syst.cd()
      histo_v2_syst[i_pt].Draw()
      std_errorbox.Draw()
      std_line.Draw()
      info_panel.Draw()

      histo_v2_syst[i_pt].Write()
      canvas_syst.Write()
      canvas_syst.SaveAs(f'{cent_plots_dir_names[i_cent]}/{canvas_syst.GetName()}.pdf')


# Final plots

print("Making final plots")
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
cV2.SaveAs(f'{output_dir_name}/plots/{cV2.GetName()}.pdf')

cPurity = ROOT.TCanvas('cPurity', 'cPurity', 800, 600)
frame = cPurity .DrawFrame(1.7, 0.5, 9, 1.1, r';#it{p}_{T} (GeV/#it{c}); purity')
cPurity.SetBottomMargin(0.13)
cPurity.SetLeftMargin(0.13)
cPurity.cd()
legend_purity = ROOT.TLegend(0.23, 0.22, 0.49, 0.57, 'FT0C centrality', 'brNDC')
legend_purity.SetBorderSize(0)
legend_purity.SetNColumns(2)

for i_cent in range(n_cent_classes):
  cent_label = f'{flow_makers[i_cent].cent_limits[0]} - {flow_makers[i_cent].cent_limits[1]}' + r'%'
  legend_purity.AddEntry(flow_makers[i_cent].hPurityVsPt, cent_label,'PF')
  flow_makers[i_cent].hPurityVsPt.Draw('SAME')

legend_purity.Draw()

output_file.cd()
cPurity.Write()
cPurity.SaveAs(f'{output_dir_name}/plots/{cPurity.GetName()}.pdf')
