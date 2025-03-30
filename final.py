import ROOT
import yaml
import argparse
import os
import numpy as np

import sys

sys.path.append("utils")
import utils as utils


parser = argparse.ArgumentParser(description="Configure the parameters of the script.")
parser.add_argument(
    "--config-file",
    dest="config_file",
    help="path to the YAML file with configuration.",
    default="",
)
args = parser.parse_args()
if args.config_file == "":
    print("** No config file provided. Exiting. **")
    exit()

config_file = open(args.config_file, "r")
config = yaml.full_load(config_file)

input_file_name = config["output_dir_name"] + config["output_file_name"]

cent_detector_label = config["cent_detector_label"]

centrality_classes = config["centrality_classes"]
n_cent = len(centrality_classes)
pt_bins = config["pt_bins"]
cent_colours = config["cent_colours"]

input_file = ROOT.TFile(input_file_name)
input_file_separated_syst = ROOT.TFile(input_file_name[:-5] + "_separated.root")
output_file_name = config["output_dir_name"] + "final.root"
output_file = ROOT.TFile(output_file_name, "recreate")
output_dir_name = config["output_dir_name"]

output_dir_plots_name = output_dir_name + "final_plots/"

if not os.path.exists(output_dir_plots_name):
    os.makedirs(output_dir_plots_name)

do_alternative_syst = True


stat_list = []
syst_list = []
separated_abs_syst_dict = {}
syst_colours = {}

# get color paletter
ROOT.gStyle.SetPalette(ROOT.kRainbow)
cols = ROOT.TColor.GetPalette()
n_cols = cols.GetSize()

# get names of cut variables
selection_dict = config["selection_dict"]
for var in selection_dict:
    separated_abs_syst_dict[var] = []

# ptdep_selection_dict = config["ptdep_selection_dict"]
# for var in ptdep_selection_dict:
#     separated_abs_syst_dict[var] = []

period = int(n_cols / len(separated_abs_syst_dict))

if do_alternative_syst:
    separated_abs_syst_dict["table"] = []
separated_abs_syst_dict["total"] = []

for i_cent in range(0, n_cent):
    cent_name = f"cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}"
    h_stat = input_file.Get(f"{cent_name}/default/hV2vsPt_{cent_name}")
    h_stat.SetDirectory(0)
    h_stat.GetYaxis().SetTitle(r"#it{v}_{2}{EP, |#Delta#eta| > 1.3}")
    utils.setHistStyle(h_stat, cent_colours[i_cent])
    stat_list.append(h_stat)

    for i_var, var in enumerate(separated_abs_syst_dict):
        if var == "total":
            continue
        elif var == "table":
            continue
        h_abs_syst_separated = input_file_separated_syst.Get(
            f"{cent_name}/{var}/hAbsSystVsPt_{var}"
        )
        h_abs_syst_separated.SetDirectory(0)
        utils.setHistStyle(h_abs_syst_separated, cols.At(n_cols - i_var * period - 1))
        separated_abs_syst_dict[var].append(h_abs_syst_separated)

    # systematic from different tables
    if do_alternative_syst:
        h_table_syst = input_file_separated_syst.Get(
            f"{cent_name}/table_comp/hV2tableDiff_{cent_name}"
        )
        h_table_syst.SetDirectory(0)
        h_table_syst.Scale(0.5)
        separated_abs_syst_dict["table"].append(h_table_syst)

    # summing systematic uncertainties
    h_syst = h_stat.Clone(f"hV2vsPt_{cent_name}_syst")
    utils.setHistStyle(h_syst, cent_colours[i_cent])
    h_abs_syst = h_stat.Clone(f"hAbsSystVsPt_{cent_name}")
    h_abs_syst.Reset()
    h_abs_syst.GetYaxis().SetTitle("systematic error")
    h_abs_syst.GetYaxis().SetRangeUser(0.0, 0.5)
    utils.setHistStyle(h_abs_syst, ROOT.kBlack)

    n_pt_bins = len(pt_bins[i_cent]) - 1

    # check wheter summing systs from single source gives larger uncertainties
    for i_pt in range(0, n_pt_bins):
        total_sum = 0.0
        for i_var, var in enumerate(separated_abs_syst_dict):
            if var == "total":
                continue
            elif var == "table":
                continue
            var_err = separated_abs_syst_dict[var][i_cent].GetBinContent(i_pt + 1)
            total_sum = total_sum + var_err * var_err
        total_sum = np.sqrt(total_sum)

        hSystDist = input_file.Get(f"{cent_name}/hV2syst_{cent_name}_pt{i_pt}")
        syst_err = hSystDist.GetRMS()
        if total_sum > syst_err:
            syst_err = total_sum
        if do_alternative_syst:
            table_err = h_table_syst.GetBinContent(i_pt + 1)
            tot_err = np.hypot(syst_err, table_err)
        else:
            tot_err = syst_err
        h_syst.SetBinError(i_pt + 1, tot_err)
        h_abs_syst.SetBinContent(i_pt + 1, tot_err)
        h_abs_syst.SetBinError(i_pt + 1, 0)

    syst_list.append(h_syst)
    separated_abs_syst_dict["total"].append(h_abs_syst)

    cV2_cent = ROOT.TCanvas(f"cV2_{cent_name}", f"cV2_{cent_name}", 800, 600)
    frame_cent = cV2_cent.DrawFrame(
        1.7,
        -0.1,
        12.0,
        1.1,
        r";#it{p}_{T} (GeV/#it{c}); #it{v}_{2}{EP, |#Delta#eta| > 1.3}",
    )
    cV2_cent.SetBottomMargin(0.13)
    cV2_cent.SetBorderSize(0)

    info_panel = ROOT.TPaveText(0.6, 0.6, 0.8, 0.82, "NDC")
    info_panel.SetBorderSize(0)
    info_panel.SetFillStyle(0)
    info_panel.SetTextAlign(12)
    info_panel.SetTextFont(42)
    info_panel.AddText(r"PbPb, #sqrt{#it{s}_{NN}} = 5.36 TeV")
    info_panel.AddText(
        f"{centrality_classes[i_cent][0]} - {centrality_classes[i_cent][1]} % {cent_detector_label}"
    )

    cV2_cent.cd()
    info_panel.Draw()
    h_stat.Draw("PE SAME")
    h_syst.Draw("PE2 SAME")

    cent_dir = output_file.mkdir(cent_name)
    cent_dir.cd()
    h_stat.Write()
    h_syst.Write()
    h_abs_syst.Write()
    cV2_cent.Write()
    cV2_cent.SaveAs(f"{output_dir_plots_name}{cV2_cent.GetName()}.pdf")

cV2 = ROOT.TCanvas("cV2_tot", "cV2_tot", 800, 600)
info_panel_total = ROOT.TPaveText(0.16, 0.74, 0.36, 0.86, "NDC")
info_panel_total.SetBorderSize(0)
info_panel_total.SetFillStyle(0)
info_panel_total.SetTextAlign(12)
info_panel_total.SetTextFont(42)
info_panel_total.SetTextSize(0.04)
info_panel_total.AddText(r"ALICE Preliminary")
info_panel_total.AddText(r"Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.36 TeV")
info_panel_total.AddText(r"{}^{3}#bar{He}, |#eta| < 0.8")
frame = cV2.DrawFrame(
    1.7,
    -0.1,
    12.0,
    1.1,
    r";#it{p}_{T} (GeV/#it{c}); #it{v}_{2}{EP, |#Delta#eta| > 1.3}",
)
cV2.SetBottomMargin(0.13)
cV2.SetLeftMargin(0.13)
cV2.SetBorderSize(0)
cV2.cd()
legend = ROOT.TLegend(0.66, 0.63, 0.90, 0.87, "FT0C centrality", "brNDC")
legend.SetBorderSize(0)
legend.SetNColumns(2)

for i_cent in range(n_cent):
    cent_label = (
        f"{centrality_classes[i_cent][0]}-{centrality_classes[i_cent][1]}" + r"%"
    )
    legend.AddEntry(stat_list[i_cent], cent_label, "PF")
    stat_list[i_cent].Draw("PEX0 SAME")
    syst_list[i_cent].Draw("PE2 SAME")

legend.Draw()
info_panel_total.Draw()

output_file.cd()
cV2.Write()
cV2.SaveAs(f"{output_dir_plots_name}{cV2.GetName()}.pdf")

# systematic uncertainties

syst_axis_limits = [0.025, 0.025, 0.025, 0.025, 0.06, 0.07, 0.07]

for i_cent in range(n_cent):
    cent_name = f"cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}"
    cent_label = (
        f"{centrality_classes[i_cent][0]} - {centrality_classes[i_cent][1]}" + r"%"
    )
    cAbsSyst = ROOT.TCanvas(
        f"cAbsSystVsPt_{cent_name}", f"cAbsSystVsPt_{cent_name}", 800, 600
    )
    cAbsSyst.SetBorderSize(0)
    cAbsSyst.SetBottomMargin(0.13)
    cAbsSyst.SetLeftMargin(0.13)
    legend_syst = ROOT.TLegend(
        0.14, 0.64, 0.40, 0.87, "Systematic uncertainties", "brNDC"
    )
    legend_syst.SetBorderSize(0)
    cAbsSyst.DrawFrame(
        1.7,
        0.0,
        12.0,
        syst_axis_limits[i_cent],
        f"{cent_label}" + r";#it{p}_{T} (GeV/#it{c}); systematic uncertainty",
    )
    for var in separated_abs_syst_dict:
        separated_abs_syst_dict[var][i_cent].Draw("histo same")
        legend_syst.AddEntry(separated_abs_syst_dict[var][i_cent], var, "PL")
    legend_syst.Draw()
    output_file.cd(cent_name)
    cAbsSyst.Write()
    cAbsSyst.SaveAs(f"{output_dir_plots_name}/{cAbsSyst.GetName()}.pdf")

# comparison with models

for i_cent in range(n_cent):
    cent_name = f"cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}"

gPredWenbin = []
coalescence_theory_file = ROOT.TFile("../theoretical_models/Predictions.root")
gPredWenbin010 = coalescence_theory_file.Get("gPredWenbin_0_10")
gPredWenbin010.SetFillStyle(1001)
gPredWenbin010.SetFillColorAlpha(ROOT.kCyan + 2, 0.6)
gPredWenbin.append(gPredWenbin010)
gPredWenbin1020 = coalescence_theory_file.Get("gPredWenbin_10_20")
gPredWenbin1020.SetFillStyle(1001)
gPredWenbin1020.SetFillColorAlpha(ROOT.kSpring - 9, 0.6)
gPredWenbin.append(gPredWenbin1020)
gPredWenbin2030 = coalescence_theory_file.Get("gPredWenbin_20_30")
gPredWenbin2030.SetFillStyle(1001)
gPredWenbin2030.SetFillColorAlpha(ROOT.kOrange + 1, 0.6)
gPredWenbin.append(gPredWenbin2030)
gPredWenbin3040 = coalescence_theory_file.Get("gPredWenbin_30_40")
gPredWenbin3040.SetFillStyle(1001)
gPredWenbin3040.SetFillColorAlpha(ROOT.kOrange + 1, 0.6)
gPredWenbin.append(gPredWenbin3040)

gPredWenbin4060 = coalescence_theory_file.Get("gPredWenbin_40_60")
gPredWenbin4060.SetFillStyle(1001)
gPredWenbin4060.SetFillColorAlpha(ROOT.kOrange + 1, 0.6)
gPredWenbin.append(gPredWenbin4060)

cV2comp = []
x_limits = [11.0, 11.0, 11.0, 9.0, 9.0]
y_limits = [0.30, 0.50, 0.60, 0.80, 1.00]

for i_cent in range(5):
    cent_name = f"cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}"
    cV2comp.append(
        ROOT.TCanvas(f"cV2comp_{cent_name}", f"cV2comp_{cent_name}", 800, 600)
    )
    framecomp = cV2comp[i_cent].DrawFrame(
        1.7,
        -0.07,
        x_limits[i_cent],
        y_limits[i_cent],
        r";#it{p}_{T} (GeV/#it{c}); #it{v}_{2}{EP, |#Delta#eta| > 1.3}",
    )
    cV2comp[i_cent].SetBottomMargin(0.13)
    cV2comp[i_cent].SetLeftMargin(0.13)
    cV2comp[i_cent].SetBorderSize(0)

    info_panel_comp = ROOT.TPaveText(0.17, 0.67, 0.37, 0.83, "NDC")
    info_panel_comp.SetBorderSize(0)
    info_panel_comp.SetFillStyle(0)
    info_panel_comp.SetTextAlign(12)
    info_panel_comp.SetTextFont(42)
    info_panel_comp.SetTextSize(0.04)
    info_panel_comp.AddText(r"ALICE Preliminary")
    info_panel_comp.AddText(r"Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.36 TeV")
    info_panel_comp.AddText(
        f"{centrality_classes[i_cent][0]}-{centrality_classes[i_cent][1]}% FT0C centrality"
    )

    legend_comp = ROOT.TLegend(0.55, 0.23, 0.87, 0.4, "", "brNDC")
    legend_comp.SetBorderSize(0)

    gPredWenbin[i_cent].Draw("E3same][")

    stat_list[i_cent].Draw("PEX0 SAME")
    syst_list[i_cent].Draw("PE2 SAME")

    legend_comp.AddEntry(stat_list[i_cent], r"{}^{3}#bar{He}, |#eta| < 0.8", "PF")
    legend_comp.AddEntry(
        gPredWenbin[i_cent],
        r"#splitline{IP Glasma + MUSIC +}{+ UrQMD + Coalescence}",
        "F",
    )

    info_panel_comp.Draw()
    legend_comp.Draw()

    cV2comp[i_cent].SaveAs(
        f"{output_dir_plots_name}{cV2comp[i_cent].GetName()}_comp.pdf"
    )

    output_file.cd()
    cV2comp[i_cent].Write()
