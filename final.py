import ROOT
import yaml
import argparse
import os
import numpy as np

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
input_file_separated_syst = ROOT.TFile(
    input_file_name[:-5] + '_separated.root')
output_file_name = config['output_dir_name'] + 'final.root'
output_file = ROOT.TFile(output_file_name, 'recreate')
output_dir_name = config['output_dir_name']

output_dir_plots_name = output_dir_name + 'final_plots/'

if not os.path.exists(output_dir_plots_name):
    os.makedirs(output_dir_plots_name)


stat_list = []
syst_list = []
separated_abs_syst_dict = {}
syst_colours = {}

# get color paletter
ROOT.gStyle.SetPalette(ROOT.kRainbow)
cols = ROOT.TColor.GetPalette()
n_cols = cols.GetSize()

# get names of cut variables
selection_dict = config['selection_dict']
for var in selection_dict:
    separated_abs_syst_dict[var] = []

ptdep_selection_dict = config['ptdep_selection_dict']
for var in ptdep_selection_dict:
    separated_abs_syst_dict[var] = []

period = int(n_cols / len(separated_abs_syst_dict))

separated_abs_syst_dict['table'] = []
separated_abs_syst_dict['total'] = []

for i_cent in range(0, n_cent):
    cent_name = f'cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}'
    h_stat = input_file.Get(f'{cent_name}/default/hV2vsPt_{cent_name}')
    h_stat.SetDirectory(0)
    stat_list.append(h_stat)

    for i_var, var in enumerate(separated_abs_syst_dict):
        if var == 'total':
            continue
        elif var == 'table':
            continue
        h_abs_syst_separated = input_file_separated_syst.Get(
            f'{cent_name}/{var}/hAbsSystVsPt_{var}')
        h_abs_syst_separated.SetDirectory(0)
        utils.setHistStyle(h_abs_syst_separated,
                           cols.At(n_cols - i_var*period - 1))
        separated_abs_syst_dict[var].append(h_abs_syst_separated)

    # systematic from different tables
    h_table_syst = input_file_separated_syst.Get(
        f'{cent_name}/table_comp/hV2tableDiff_{cent_name}')
    h_table_syst.SetDirectory(0)
    h_table_syst.Scale(0.5)
    separated_abs_syst_dict['table'].append(h_table_syst)

    # summing systematic uncertainties
    h_syst = h_stat.Clone(f'hV2vsPt_{cent_name}_syst')
    h_abs_syst = h_stat.Clone(f'hAbsSystVsPt_{cent_name}')
    h_abs_syst.Reset()
    h_abs_syst.GetYaxis().SetTitle('systematic error')
    h_abs_syst.GetYaxis().SetRangeUser(0., 0.5)
    utils.setHistStyle(h_abs_syst, ROOT.kBlack)

    n_pt_bins = len(pt_bins[i_cent]) - 1

    # check wheter summing systs from single source gives larger uncertainties
    for i_pt in range(0, n_pt_bins):
        total_sum = 0.
        for i_var, var in enumerate(separated_abs_syst_dict):
            if var == 'total':
                continue
            elif var == 'table':
                continue
            var_err = separated_abs_syst_dict[var][i_cent].GetBinContent(i_pt+1)
            total_sum = total_sum + var_err * var_err
        total_sum = np.sqrt(total_sum)

        hSystDist = input_file.Get(f'{cent_name}/hV2syst_{cent_name}_pt{i_pt}')
        syst_err = hSystDist.GetRMS()
        if total_sum > syst_err:
            syst_err = total_sum
        table_err = h_table_syst.GetBinContent(i_pt+1)
        tot_err = np.hypot(syst_err, table_err)
        h_syst.SetBinError(i_pt+1, tot_err)
        h_abs_syst.SetBinContent(i_pt+1, tot_err)
        h_abs_syst.SetBinError(i_pt+1, 0)

    syst_list.append(h_syst)
    separated_abs_syst_dict['total'].append(h_abs_syst)

    cV2_cent = ROOT.TCanvas(f'cV2_{cent_name}', f'cV2_{cent_name}', 800, 600)
    frame_cent = cV2_cent.DrawFrame(
        1.7, -0.2, 9, 1., r';#it{p}_{T} (GeV/#it{c}); v_{2}')
    cV2_cent.SetBottomMargin(0.13)
    cV2_cent.SetBorderSize(0)

    info_panel = ROOT.TPaveText(0.6, 0.6, 0.8, 0.82, 'NDC')
    info_panel.SetBorderSize(0)
    info_panel.SetFillStyle(0)
    info_panel.SetTextAlign(12)
    info_panel.SetTextFont(42)
    info_panel.AddText(r'PbPb, #sqrt{#it{s}_{NN}} = 5.36 TeV')
    info_panel.AddText(
        f'{centrality_classes[i_cent][0]} - {centrality_classes[i_cent][1]} % {cent_detector_label}')

    cV2_cent.cd()
    info_panel.Draw()
    h_stat.Draw('PE SAME')
    h_syst.Draw('PE2 SAME')

    cent_dir = output_file.mkdir(cent_name)
    cent_dir.cd()
    h_stat.Write()
    h_syst.Write()
    h_abs_syst.Write()
    cV2_cent.Write()
    cV2_cent.SaveAs(f'{output_dir_plots_name}{cV2_cent.GetName()}.pdf')

cV2 = ROOT.TCanvas('cV2_tot', 'cV2_tot', 800, 600)
info_panel_total = ROOT.TPaveText(0.22, 0.67, 0.42, 0.83, 'NDC')
info_panel_total.SetBorderSize(0)
info_panel_total.SetFillStyle(0)
info_panel_total.SetTextAlign(12)
info_panel_total.SetTextFont(42)
info_panel_total.SetTextSize(0.04)
info_panel_total.AddText(r'ALICE Preliminary')
info_panel_total.AddText(r'Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.36 TeV')
info_panel_total.AddText(r'{}^{3}#bar{He}, |#eta| < 0.8')
frame = cV2.DrawFrame(1.7, -0.2, 8., 1., r';#it{p}_{T} (GeV/#it{c}); #it{v}_{2}{EP}')
cV2.SetBottomMargin(0.13)
cV2.SetLeftMargin(0.13)
cV2.SetBorderSize(0)
cV2.cd()
legend = ROOT.TLegend(0.56, 0.67, 0.82, 0.83, 'FT0C centrality', 'brNDC')
legend.SetBorderSize(0)
legend.SetNColumns(2)

for i_cent in range(n_cent):
    cent_label = f'{centrality_classes[i_cent][0]}-{centrality_classes[i_cent][1]}' + r'%'
    legend.AddEntry(stat_list[i_cent], cent_label, 'PF')
    stat_list[i_cent].Draw('PEX0 SAME')
    syst_list[i_cent].Draw('PE2 SAME')

legend.Draw()
info_panel_total.Draw()

output_file.cd()
cV2.Write()
cV2.SaveAs(f'{output_dir_plots_name}{cV2.GetName()}.pdf')

# systematic uncertainties

for i_cent in range(n_cent):
    cent_name = f'cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}'
    cent_label = f'{centrality_classes[i_cent][0]} - {centrality_classes[i_cent][1]}' + r'%'
    cAbsSyst = ROOT.TCanvas(
        f'cAbsSystVsPt_{cent_name}', f'cAbsSystVsPt_{cent_name}', 800, 600)
    cAbsSyst.SetBorderSize(0)
    cAbsSyst.SetBottomMargin(0.13)
    cAbsSyst.SetLeftMargin(0.13)
    legend_syst = ROOT.TLegend(
        0.14, 0.64, 0.40, 0.87, 'Systematic uncertainties', 'brNDC')
    legend_syst.SetBorderSize(0)
    cAbsSyst.DrawFrame(
        1.7, 0., 9, 0.025, f'{cent_label}' + r';#it{p}_{T} (GeV/#it{c}); systematic uncertainty')
    for var in separated_abs_syst_dict:
        separated_abs_syst_dict[var][i_cent].Draw('histo same')
        legend_syst.AddEntry(separated_abs_syst_dict[var][i_cent], var, 'PL')
    legend_syst.Draw()
    output_file.cd(cent_name)
    cAbsSyst.Write()
    cAbsSyst.SaveAs(f'{output_dir_plots_name}/{cAbsSyst.GetName()}.pdf')

# comparison with models

cV2comp_0_20 = ROOT.TCanvas('cV2comp_0_20', 'cV2comp_0_20', 800, 600)
framecomp_0_20 = cV2comp_0_20.DrawFrame(
    1.7, -0.07, 8., 0.35, r';#it{p}_{T} (GeV/#it{c}); #it{v}_{2}{EP}')
cV2comp_0_20.SetBottomMargin(0.13)
cV2comp_0_20.SetLeftMargin(0.13)
cV2comp_0_20.SetBorderSize(0)

info_panel_comp_0_20 = ROOT.TPaveText(0.17, 0.67, 0.37, 0.83, 'NDC')
info_panel_comp_0_20.SetBorderSize(0)
info_panel_comp_0_20.SetFillStyle(0)
info_panel_comp_0_20.SetTextAlign(12)
info_panel_comp_0_20.SetTextFont(42)
info_panel_comp_0_20.SetTextSize(0.04)
info_panel_comp_0_20.AddText(r'ALICE Preliminary')
info_panel_comp_0_20.AddText(r'Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.36 TeV')
info_panel_comp_0_20.AddText(r'0-20% FT0C centrality')

legend_comp_0_20 = ROOT.TLegend(0.54, 0.16, 0.9, 0.38, '', 'brNDC')
legend_comp_0_20.SetBorderSize(0)

merged_file = ROOT.TFile('v2_0_20_merged.root')
hV2vsPt_0020 = merged_file.Get('hV2vsPt_0020')
hV2vsPt_0020_syst = merged_file.Get('hV2vsPt_0020_syst')
utils.setHistStyle(hV2vsPt_0020, ROOT.kRed+1)
utils.setHistStyle(hV2vsPt_0020_syst, ROOT.kRed+1)

wenbin_theory_file = ROOT.TFile('../theoretical_models/WenbinPredictions.root')
gPredWenbin020 = wenbin_theory_file.Get('gPredWenbin_0_20')
gPredWenbin020.SetFillStyle(1001)
gPredWenbin020.SetFillColorAlpha(ROOT.kCyan+2, 0.6)
gPredWenbin2040 = wenbin_theory_file.Get('gPredWenbin_20_40')
gPredWenbin2040.SetFillStyle(1001)
gPredWenbin2040.SetFillColorAlpha(ROOT.kViolet+1, 0.6)

BW_0_20_file = ROOT.TFile('../theoretical_models/BWflow_0_20.root')
v2_BW_0_20 = BW_0_20_file.Get('fv23He')
v2_BW_0_20.SetLineColor(ROOT.kCyan+2)
BW_20_40_file = ROOT.TFile('../theoretical_models/BWflow_20_40.root')
v2_BW_20_40 = BW_20_40_file.Get('fv23He')
v2_BW_20_40.SetLineColor(ROOT.kViolet+1)

cV2comp_0_20.cd()
gPredWenbin020.Draw('E3same][')
v2_BW_0_20.Draw('L SAME')
hV2vsPt_0020.Draw('PEX0 SAME')
hV2vsPt_0020_syst.Draw('PE2 SAME')

legend_comp_0_20.AddEntry(hV2vsPt_0020, r'{}^{3}#bar{He}, |#eta| < 0.8', 'PF')
legend_comp_0_20.AddEntry(
    gPredWenbin020, r'#splitline{IP Glasma + MUSIC +}{+ UrQMD + Coalescence}', 'F')
legend_comp_0_20.AddEntry(v2_BW_0_20, r'Blast-wave', 'L')

info_panel_comp_0_20.Draw()
legend_comp_0_20.Draw()

cV2comp_0_20.SaveAs(f'{output_dir_plots_name}{cV2comp_0_20.GetName()}.pdf')

cV2comp_20_40 = ROOT.TCanvas('cV2comp_20_40', 'cV2comp_20_40', 800, 600)
framecomp_20_40 = cV2comp_20_40.DrawFrame(
    1.7, -0.05, 8., 0.6, r';#it{p}_{T} (GeV/#it{c}); #it{v}_{2}{EP}')
cV2comp_20_40.SetBottomMargin(0.13)
cV2comp_20_40.SetLeftMargin(0.13)
cV2comp_20_40.SetBorderSize(0)

cV2comp_20_40.cd()

info_panel_comp_20_40 = ROOT.TPaveText(0.17, 0.67, 0.37, 0.83, 'NDC')
info_panel_comp_20_40.SetBorderSize(0)
info_panel_comp_20_40.SetFillStyle(0)
info_panel_comp_20_40.SetTextAlign(12)
info_panel_comp_20_40.SetTextFont(42)
info_panel_comp_20_40.SetTextSize(0.04)
info_panel_comp_20_40.AddText(r'ALICE Preliminary')
info_panel_comp_20_40.AddText(r'Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.36 TeV')
info_panel_comp_20_40.AddText(r'20-40% FT0C centrality')

legend_comp_20_40 = ROOT.TLegend(0.54, 0.18, 0.90, 0.41, '', 'brNDC')
legend_comp_20_40.SetBorderSize(0)

gPredWenbin2040.Draw('E3same][')
v2_BW_20_40.Draw('L SAME')
stat_list[2].Draw('PEX0 SAME')
syst_list[2].Draw('PE2 SAME')

legend_comp_20_40.AddEntry(stat_list[2], r'{}^{3}#bar{He}, |#eta| < 0.8', 'PF')
legend_comp_20_40.AddEntry(
    gPredWenbin2040, r'#splitline{IP Glasma + MUSIC +}{+ UrQMD + Coalescence}', 'F')
legend_comp_20_40.AddEntry(v2_BW_20_40, r'Blast-wave', 'L')

info_panel_comp_20_40.Draw()
legend_comp_20_40.Draw()

cV2comp_20_40.SaveAs(f'{output_dir_plots_name}{cV2comp_20_40.GetName()}.pdf')

output_file.cd()
cV2comp_0_20.Write()
cV2comp_20_40.Write()
