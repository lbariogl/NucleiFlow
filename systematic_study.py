import ROOT
from hipe4ml.tree_handler import TreeHandler
import pandas as pd
import yaml
import argparse
import numpy as np
import copy

import sys
sys.path.append('utils')
from flow import FlowMaker
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

input_file_name = config['input_file_name']
output_dir_name = config['output_dir_name']
output_file_name = config['output_file_name']
output_file_name = output_file_name[:-5] + '_separated.root'


nuclei_tree_name = config['nuclei_tree_name']
ep_tree_name = config['ep_tree_name']

mandatory_selections = config['mandatory_selections']
selection_dict = config['selection_dict']

cent_detector_label = config['cent_detector_label']

centrality_classes = config['centrality_classes']
pt_bins = config['pt_bins']
cent_colours = config['cent_colours']

# BB parameters
p_train = config['p_train']
resolution_train = config['resolution_train']
n_sigma_plot = config['n_sigma_plot']

# create output file
output_file = ROOT.TFile(f'{output_dir_name}/{output_file_name}', 'recreate')

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

# get Standard Spectrum
standard_file_name = f'{output_dir_name}/' + config['output_file_name']
standard_file = ROOT.TFile(standard_file_name)

n_cent_classes = len(centrality_classes)

default_v2_histos = []
default_v2_values = []
cent_dir_names = []
cent_dirs = []

for i_cent in range(n_cent_classes):
    cent_dir_name = f'cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}'
    cent_dir_names.append(cent_dir_name)
    cent_dir = output_file.mkdir(cent_dir_name)
    cent_dirs.append(cent_dir)
    histo = standard_file.Get(
        f'{cent_dir_name}/default/hV2vsPt_{cent_dir_name}')
    histo.SetDirectory(0)
    default_v2_histos.append(histo)
    default_v2_values.append(utils.getValuesFromHisto(histo))


# get a unique df from nuclei and ep trees
nuclei_hdl = TreeHandler(
    input_file_name, f'{nuclei_tree_name};', folder_name='DF*')
nuclei_df = nuclei_hdl._full_data_frame

nucleiflow_hdl = TreeHandler(
    input_file_name, f'{ep_tree_name};', folder_name='DF*')
nucleiflow_df = nucleiflow_hdl._full_data_frame

complete_df = pd.concat([nuclei_df, nucleiflow_df], axis=1, join='inner')

# define new columns
utils.redifineColumns(complete_df)

# apply mandatory selections
complete_df.query(mandatory_selections, inplace=True)


#########################
#     varied cuts
#########################

print("** Starting systematic variations **")

# create a dictionary with all the possible selections for a specific variable
cut_dict_syst = config['cut_dict_syst']
cut_string_dict = {}
cut_label_dict = {}
for var in cut_dict_syst:
    var_dict = cut_dict_syst[var]
    cut_greater = var_dict['cut_greater']
    cut_greater_string = " > " if cut_greater else " < "

    cut_list = var_dict['cut_list']
    cut_arr = np.linspace(cut_list[0], cut_list[1], cut_list[2])
    cut_string_dict[var] = []
    cut_label_dict[var] = []
    for cut in cut_arr:
        sel_string = var + cut_greater_string + str(cut)
        for var2 in selection_dict:
            if var2 != var:
                sel_string = sel_string + f' and {selection_dict[var2]}'
        cut_string_dict[var].append(sel_string)
        cut_label_dict[var].append(cut_greater_string + f'{cut:.3f}')

print("  ** separated cuts **")

for i_cent in range(n_cent_classes):

    print(
        f'Centrality: {centrality_classes[i_cent][0]} - {centrality_classes[i_cent][1]}')

    v2_dict = {}
    canvas_dict = {}
    legend_dict = {}

    for var, cuts in cut_string_dict.items():
        print(f'var: {var}')
        var_dir = cent_dirs[i_cent].mkdir(var)

        v2_dict[var] = []
        canvas_dict[var] = ROOT.TCanvas(f'c{var}_cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}',
                                        f'c{var}_cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}', 800, 600)
        legend_dict[var] = ROOT.TLegend(0.16, 0.61, 0.86, 0.88, '', 'brNDC')
        legend_dict[var].SetBorderSize(0)

        for i_cut, cut in enumerate(cuts):

            print(f'{var}: {i_cut} / {len(cuts)} ==> {cut}')

            output_dir_varied = var_dir.mkdir(f'{i_cut}')

            flow_maker_syst = FlowMaker()
            flow_maker_syst.data_df = complete_df
            flow_maker_syst.pt_bins = pt_bins[i_cent]
            flow_maker_syst.cent_limits = centrality_classes[i_cent]
            flow_maker_syst.resolution = resolutions[i_cent]

            var_suffix = f'_{var}_{i_cut}'
            flow_maker_syst.selection_string = cut
            flow_maker_syst.suffix = var_suffix

            flow_maker_syst.make_flow()
            flow_values = flow_maker_syst.getFlowValues()
            # write histograms to file
            flow_maker_syst.output_dir = output_dir_varied
            flow_maker_syst.dump_to_output_file()

            histo = copy.deepcopy(flow_maker_syst.hV2vsPt)

            v2_dict[var].append(histo)

            del flow_maker_syst

    # get color paletter
    cols = ROOT.TColor.GetPalette()

    cent_dirs[i_cent].cd()
    cent_dirs[i_cent].mkdir('std')
    default_v2_histos[i_cent].Write()

    for var, histos in v2_dict.items():
        cent_dirs[i_cent].cd(f'{var}')
        canvas_dict[var].cd()
        canvas_title = f'{var}' + r';#it{p}_{T} (GeV/#it{c}); v_{2}'
        canvas_dict[var].DrawFrame(
            1.7, -0.2, 9., 1., canvas_title)
        period = int(cols.GetSize() / len(histos))
        for i_histo, histo in enumerate(histos):
            utils.setHistStyle(histo, cols.At(i_histo*period))
            legend_dict[var].AddEntry(
                histo, f'{cut_label_dict[var][i_histo]}', 'PE')
            histo.Draw('PE SAME')
        legend_dict[var].AddEntry(
            default_v2_histos[i_cent], 'std', 'PE')
        default_v2_histos[i_cent].Draw('PE SAME')
        legend_dict[var].Draw()
        legend_dict[var].SetNColumns(5)
        canvas_dict[var].Write()
