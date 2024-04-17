import ROOT
import yaml
import argparse
import numpy as np

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

output_dir_name = config['output_dir_name']

cent_detector_label = config['cent_detector_label']
centrality_classes = config['centrality_classes']
pt_bins = config['pt_bins']
cent_colours = config['cent_colours']

# retreiving cuts to be tested
cut_dict_syst = config['cut_dict_syst']
# cluster size
its_clus_size_cut_dict = cut_dict_syst['fAvgItsClusSize']
its_clus_size_cut_list = its_clus_size_cut_dict['cut_list']
its_clus_size_cut_arr = np.linspace(
    its_clus_size_cut_list[0], its_clus_size_cut_list[1], its_clus_size_cut_list[2])
n_its_clus_size_cuts = len(its_clus_size_cut_arr)
its_clus_size_cut_greater = its_clus_size_cut_dict['cut_greater']
its_clus_size_cut_greater_string = " > " if its_clus_size_cut_greater else " < "
its_clus_size_cut_label_dict = []
for cut in its_clus_size_cut_arr:
    its_clus_size_cut_label_dict.append(its_clus_size_cut_greater_string + f'{cut:.3f}')

# cluster size * cos(lambda)
its_clus_size_coslambda_cut_dict = cut_dict_syst['fAvgItsClusSizeCosLambda']
its_clus_size_coslambda_cut_list = its_clus_size_cut_dict['cut_list']
its_clus_size_coslambda_cut_arr = np.linspace(
    its_clus_size_coslambda_cut_list[0], its_clus_size_coslambda_cut_list[1], its_clus_size_coslambda_cut_list[2])
n_its_clus_size_coslambda_cuts = len(its_clus_size_coslambda_cut_arr)
its_clus_size_coslambda_cut_greater = its_clus_size_coslambda_cut_dict['cut_greater']
its_clus_size_coslambda_cut_greater_string = " > " if its_clus_size_cut_greater else " < "
its_clus_size_coslambda_cut_label_dict = []
for cut in its_clus_size_coslambda_cut_arr:
    its_clus_size_coslambda_cut_label_dict.append(its_clus_size_coslambda_cut_greater_string + f'{cut:.3f}')


# configurations
input_file = ROOT.TFile('../results_default/flow_separated.root')
input_file_coslambda = ROOT.TFile('../results_coslambda/flow_separated.root')

output_file = ROOT.TFile('../results/its_cluster_size_study.root', 'RECREATE')
n_cent_classes = len(centrality_classes)

# get color paletter
cols = ROOT.TColor.GetPalette()

for i_cent in range(n_cent_classes):

    # create output_dir
    cent_dir = output_file.mkdir(
        f'cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}')

    its_clus_size_dir = cent_dir.mkdir('fAvgItsClusSize')
    its_clus_size_coslambda_dir = cent_dir.mkdir('fAvgItsClusSizeCosLambda')

    for i_pt in range(0, len(pt_bins[i_cent])-1):

        pt_limits = [pt_bins[i_cent][i_pt], pt_bins[i_cent][i_pt + 1]]

        # its cluster size
        canvas_its_clus_size = ROOT.TCanvas(f'cNsigmaItsClustSizeCheck_cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}_pt{i_pt}',
                                            f'cNsigmaItsClustSizeCheck_cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}_pt{i_pt}', 800, 600)
        canvas_its_clus_size.SetBottomMargin(0.13)
        canvas_its_clus_size.SetLeftMargin(0.13)

        legend_its_clus_size = ROOT.TLegend(0.18, 0.66, 0.57, 0.83, 'fAvgItsClusSize', 'brNDC')
        legend_its_clus_size.SetBorderSize(0)
        legend_its_clus_size.SetNColumns(5)

        # set colour span
        period = int(cols.GetSize() / n_its_clus_size_cuts)

        for i_cut in range(0, n_its_clus_size_cuts):

            histo_name = f'cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}/fAvgItsClusSize/{i_cut}/pt_{pt_limits[0]}_{pt_limits[1]}/hNsigma3He_cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}_pt{i_pt}_fAvgItsClusSize_{i_cut}'
            histo = input_file.Get(histo_name)
            histo.SetDirectory(0)
            utils.setHistStyle(histo, cols.At(i_cut*period), linewidth=2)

            canvas_its_clus_size.cd()
            if i_cut == 0 :
                histo.Draw()
            else:
                histo.Draw('SAME')
            legend_its_clus_size.AddEntry(histo, its_clus_size_cut_label_dict[i_cut], 'PE')

        legend_its_clus_size.Draw()
        its_clus_size_dir.cd()
        canvas_its_clus_size.Write()

        # its cluster size * cos(lambda)
        canvas_its_clus_size_coslambda = ROOT.TCanvas(f'cNsigmaItsClustCosLambdaSizeCheck_cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}_pt{i_pt}',
                                            f'cNsigmaItsClustCosLambdaSizeCheck_cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}_pt{i_pt}', 800, 600)
        canvas_its_clus_size_coslambda.SetBottomMargin(0.13)
        canvas_its_clus_size_coslambda.SetLeftMargin(0.13)

        legend_its_clus_size_coslambda = ROOT.TLegend(0.18, 0.66, 0.57, 0.83, 'fAvgItsClusSizeCosLambda', 'brNDC')
        legend_its_clus_size_coslambda.SetBorderSize(0)
        legend_its_clus_size_coslambda.SetNColumns(5)

        # set colour span
        period = int(cols.GetSize() / n_its_clus_size_coslambda_cuts)

        for i_cut in range(0, n_its_clus_size_coslambda_cuts):

            histo_name = f'cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}/fAvgItsClusSizeCosLambda/{i_cut}/pt_{pt_limits[0]}_{pt_limits[1]}/hNsigma3He_cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}_pt{i_pt}_fAvgItsClusSizeCosLambda_{i_cut}'
            histo = input_file_coslambda.Get(histo_name)
            histo.SetDirectory(0)
            utils.setHistStyle(histo, cols.At(i_cut*period),linewidth=2)

            canvas_its_clus_size.cd()
            if i_cut == 0 :
                histo.Draw()
            else:
                histo.Draw('SAME')
            legend_its_clus_size_coslambda.AddEntry(histo, its_clus_size_cut_label_dict[i_cut], 'PE')

        legend_its_clus_size_coslambda.Draw()
        its_clus_size_coslambda_dir.cd()
        canvas_its_clus_size_coslambda.Write()
