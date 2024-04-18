import ROOT
import yaml
import argparse

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

input_file_name = config['output_dir_name'] + config['output_file_name']

cent_detector_label = config['cent_detector_label']

centrality_classes = config['centrality_classes']
n_cent = len(centrality_classes)
pt_bins = config['pt_bins']
cent_colours = config['cent_colours']

input_file = ROOT.TFile(input_file_name)
output_file_name = config['output_dir_name'] + 'final.root'
output_file = ROOT.TFile(output_file_name, 'recreate')

stat_list = []
syst_list = []

for i_cent in range(0, n_cent):
  cent_name = f'cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}'
  h_stat = input_file.Get(f'{cent_name}/default/hV2vsPt_{cent_name}')
  h_stat.SetDirectory(0)
  stat_list.append(h_stat)

  h_syst = h_stat.Clone(f'hV2vsPt_{cent_name}_syst')

  n_pt_bins = len(pt_bins[i_cent]) - 1
  for i_pt in range(0, n_pt_bins):
    hSystDist = input_file.Get(f'{cent_name}/hV2syst_{cent_name}_pt{i_pt}')
    syst_err = hSystDist.GetRMS()
    h_syst.SetBinError(i_pt+1, syst_err)

  syst_list.append(h_syst)

  cV2_cent = ROOT.TCanvas(f'cV2_{cent_name}', 'cV2_{cent_name}', 800, 600)
  frame_cent = cV2_cent.DrawFrame(1.7, -0.2, 9, 1., r';#it{p}_{T} (GeV/#it{c}); v_{2}')
  cV2_cent.SetBottomMargin(0.13)

  info_panel = ROOT.TPaveText(0.6, 0.6, 0.8, 0.82, 'NDC')
  info_panel.SetBorderSize(0)
  info_panel.SetFillStyle(0)
  info_panel.SetTextAlign(12)
  info_panel.SetTextFont(42)
  info_panel.AddText(r'PbPb, #sqrt{#it{s}_{NN}} = 5.36 TeV')
  info_panel.AddText(f'{centrality_classes[i_cent][0]} - {centrality_classes[i_cent][1]} % {cent_detector_label}')

  cV2_cent.cd()
  info_panel.Draw()
  h_stat.Draw('PE SAME')
  h_syst.Draw('PE2 SAME')

  cent_dir = output_file.mkdir(cent_name)
  cent_dir.cd()
  h_stat.Write()
  h_syst.Write()
  cV2_cent.Write()

cV2 = ROOT.TCanvas('cV2', 'cV2', 800, 600)
frame = cV2.DrawFrame(1.7, -0.2, 9, 1., r';#it{p}_{T} (GeV/#it{c}); v_{2}')
cV2.SetBottomMargin(0.13)
cV2.cd()
legend = ROOT.TLegend(0.61, 0.58, 0.87,0.81, 'FT0C centrality', 'brNDC')
legend.SetBorderSize(0)
legend.SetNColumns(2)

for i_cent in range(n_cent):
  cent_label = f'{centrality_classes[i_cent][0]} - {centrality_classes[i_cent][1]}' + r'%'
  legend.AddEntry(stat_list[i_cent], cent_label,'PF')
  stat_list[i_cent].Draw('PE SAME')
  syst_list[i_cent].Draw('PE2 SAME')

legend.Draw()

output_file.cd()
cV2.Write()
