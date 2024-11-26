import ROOT
from hipe4ml.tree_handler import TreeHandler
import pandas as pd
import yaml
import argparse
import numpy as np
import copy
import os

import sys

sys.path.append("utils")
import utils as utils
from flow import FlowMaker

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

input_file_name = config["input_file_name"]
resolution_file_name = config["resolution_file_name"]
output_dir_name = config["output_dir_name"]
output_file_name = config["output_file_name"]
output_file_name = output_file_name[:-5] + "_separated.root"


nuclei_tree_name = config["nuclei_tree_name"]
ep_tree_name = config["ep_tree_name"]

# Tof analysis (True for 4He)
tof_analysis = config["tof_analysis"]

mandatory_selections = config["mandatory_selections"]
selection_dict = config["selection_dict"]
standard_selection_list = selection_dict.values()
standard_selections = " and ".join(standard_selection_list)

standard_ptdep_selection_dict = config["ptdep_selection_dict"][
    "fAvgItsClusSizeCosLambda"
]

cent_detector_label = config["cent_detector_label"]

centrality_classes = config["centrality_classes"]
pt_bins = config["pt_bins"]
cent_colours = config["cent_colours"]

use_Barlow = config["use_Barlow"]

# BB parameters
p_train = config["p_train"]
resolution_train = config["resolution_train"]
n_sigma_plot = config["n_sigma_plot"]

# create output file
if not os.path.exists(output_dir_name):
    os.makedirs(output_dir_name)
output_file = ROOT.TFile(f"{output_dir_name}/{output_file_name}", "recreate")

# Get resolution from file
resolution_file = ROOT.TFile(resolution_file_name)
hResolution = resolution_file.Get("Resolution/hResolution_FT0C_TPCl_TPCr")
hResolution.SetDirectory(0)

res_0_10 = hResolution.GetBinContent(1)
res_10_20 = hResolution.GetBinContent(2)
res_20_40 = (hResolution.GetBinContent(3) + hResolution.GetBinContent(4)) / 2
res_40_60 = (hResolution.GetBinContent(5) + hResolution.GetBinContent(6)) / 2

res_20_30 = hResolution.GetBinContent(3)
res_30_40 = hResolution.GetBinContent(4)
res_40_50 = hResolution.GetBinContent(5)
res_50_60 = hResolution.GetBinContent(6)
res_60_80 = (hResolution.GetBinContent(7) + hResolution.GetBinContent(8)) / 2

# resolutions = [res_0_10, res_10_20, res_20_40, res_40_60]
resolutions = [
    res_0_10,
    res_10_20,
    res_20_30,
    res_30_40,
    res_40_50,
    res_50_60,
    res_60_80,
]

# get Standard Spectrum
standard_file_name = f"{output_dir_name}/" + config["output_file_name"]
standard_file = ROOT.TFile(standard_file_name)

do_alternative_systematic = False

# get alternative flow table
input_file_alternative = ROOT.TFile("../results_pass4_alternative/flow.root")

n_cent_classes = len(centrality_classes)

default_v2_histos = []
alternative_v2_histos = []
default_v2_values = []
cent_dir_names = []
cent_dirs = []

for i_cent in range(n_cent_classes):
    cent_dir_name = (
        f"cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}"
    )
    cent_dir_names.append(cent_dir_name)
    cent_dir = output_file.mkdir(cent_dir_name)
    cent_dirs.append(cent_dir)
    histo = standard_file.Get(f"{cent_dir_name}/default/hV2vsPt_{cent_dir_name}")
    histo.SetDirectory(0)
    default_v2_histos.append(histo)
    default_v2_values.append(utils.getValuesFromHisto(histo))
    histo_varied = input_file_alternative.Get(
        f"{cent_dir_name}/default/hV2vsPt_{cent_dir_name}"
    )
    alternative_v2_histos.append(histo_varied)


# get a unique df from nuclei and ep trees
nuclei_hdl = TreeHandler(input_file_name, f"{nuclei_tree_name};", folder_name="DF*")
nuclei_df = nuclei_hdl._full_data_frame

nucleiflow_hdl = TreeHandler(input_file_name, f"{ep_tree_name};", folder_name="DF*")
nucleiflow_df = nucleiflow_hdl._full_data_frame

complete_df = pd.concat([nuclei_df, nucleiflow_df], axis=1, join="inner")

# define new columns
utils.redefineColumns(complete_df)

# apply mandatory selections
complete_df.query(mandatory_selections, inplace=True)


#########################
#     varied cuts
#########################

print("** Starting systematic variations **")

# create a dictionary with all the possible selections for a specific variable

# pt-independent selections
cut_dict_syst = config["cut_dict_syst"]
cut_string_dict = {}
cut_label_dict = {}
for var in cut_dict_syst:
    var_dict = cut_dict_syst[var]
    cut_greater = var_dict["cut_greater"]
    cut_greater_string = " > " if cut_greater else " < "

    cut_list = var_dict["cut_list"]
    cut_arr = np.linspace(cut_list[0], cut_list[1], cut_list[2])
    cut_string_dict[var] = []
    cut_label_dict[var] = []
    for cut in cut_arr:
        sel_string = var + cut_greater_string + str(cut)
        for var2 in selection_dict:
            if var2 != var:
                sel_string = sel_string + f" and {selection_dict[var2]}"
        cut_string_dict[var].append(sel_string)
        cut_label_dict[var].append(cut_greater_string + f"{cut:.3f}")

# pt-dependent selections
ptdep_cut_dict_syst = config["ptdep_cut_dict_syst"]
n_ptdep_variations_syst = config["n_ptdep_variations_syst"]["fAvgItsClusSizeCosLambda"]

ptdep_cut_string_dict = {}

for pt_key, pt_el in ptdep_cut_dict_syst["fAvgItsClusSizeCosLambda"].items():
    ptdep_cut_greater = pt_el["cut_greater"]
    ptdep_cut_greater_string = " > " if ptdep_cut_greater else " < "
    ptdep_cut_list = pt_el["cut_list"]
    ptdep_cut_arr = np.linspace(
        ptdep_cut_list[0], ptdep_cut_list[1], n_ptdep_variations_syst
    )
    ptdep_cut_string_dict[pt_key] = []
    for ptdep_cut in ptdep_cut_arr:
        ptdep_cut_string_dict[pt_key].append(
            "fAvgItsClusSizeCosLambda" + ptdep_cut_greater_string + str(ptdep_cut)
        )

ptdep_cut_string_dict_list = []
for i_variation in range(n_ptdep_variations_syst):
    variation_dict = {}
    for pt_key in ptdep_cut_dict_syst["fAvgItsClusSizeCosLambda"]:
        variation_dict[pt_key] = ptdep_cut_string_dict[pt_key][i_variation]
    ptdep_cut_string_dict_list.append(variation_dict)


n_ptdep_variations = len(ptdep_cut_string_dict_list)

cut_label_dict["fAvgItsClusSizeCosLambda"] = []
for i_var in range(n_ptdep_variations):

    index = int(i_var - (n_ptdep_variations - 1) / 2)

    ptdep_suffix_string = f"fAvgItsClusSizeCosLambda > def. + {index}*step"
    cut_label_dict["fAvgItsClusSizeCosLambda"].append(ptdep_suffix_string)


print("  ** separated cuts **")

for i_cent in range(n_cent_classes):

    print(
        f"Centrality: {centrality_classes[i_cent][0]} - {centrality_classes[i_cent][1]}"
    )

    v2_dict = {}
    canvas_dict = {}
    legend_dict = {}
    dict_hist_abs_syst = {}

    # pt-independent cuts
    for var, cuts in cut_string_dict.items():
        print(f"var: {var}")
        var_dir = cent_dirs[i_cent].mkdir(var)

        v2_dict[var] = []
        canvas_dict[var] = ROOT.TCanvas(
            f"c{var}_cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}",
            f"c{var}_cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}",
            800,
            600,
        )
        legend_dict[var] = ROOT.TLegend(0.16, 0.61, 0.86, 0.88, "", "brNDC")
        legend_dict[var].SetBorderSize(0)

        histo_v2_syst = []
        n_pt_bins = len(pt_bins[i_cent]) - 1
        for i_pt in range(0, n_pt_bins):
            histo_v2_syst.append(
                ROOT.TH1F(
                    f"hV2syst_{var}_cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}_pt{i_pt}",
                    ";v_{2}",
                    20,
                    default_v2_values[i_cent][i_pt][0]
                    - 3 * default_v2_values[i_cent][i_pt][1],
                    default_v2_values[i_cent][i_pt][0]
                    + 3 * default_v2_values[i_cent][i_pt][1],
                )
            )

        for i_cut, cut in enumerate(cuts):

            print(f"{var}: {i_cut} / {len(cuts)} ==> {cut}")

            output_dir_varied = var_dir.mkdir(f"{i_cut}")

            flow_maker_syst = FlowMaker()
            flow_maker_syst.data_df = complete_df
            flow_maker_syst.tof_analysis = tof_analysis
            flow_maker_syst.pt_bins = pt_bins[i_cent]
            flow_maker_syst.cent_limits = centrality_classes[i_cent]
            flow_maker_syst.resolution = resolutions[i_cent]

            var_suffix = f"_{var}_{i_cut}"
            flow_maker_syst.selection_string = cut
            flow_maker_syst.ptdep_selection_dict = standard_ptdep_selection_dict
            flow_maker_syst.suffix = var_suffix

            flow_maker_syst.create_histograms()
            flow_maker_syst.make_flow()
            flow_values = flow_maker_syst.getFlowValues()
            # write histograms to file
            flow_maker_syst.output_dir = output_dir_varied
            flow_maker_syst.dump_to_output_file()

            for i_pt in range(0, flow_maker_syst.nPtBins()):
                print(f"pt: {i_pt} / {flow_maker_syst.nPtBins()}")
                if use_Barlow:
                    if utils.passBarlow(
                        default_v2_values[i_cent][i_pt][0],
                        flow_values[i_pt][0],
                        default_v2_values[i_cent][i_pt][1],
                        flow_values[i_pt][1],
                    ):
                        histo_v2_syst[i_pt].Fill(flow_values[i_pt][0])
                    else:
                        print("    Rejected for Barlow")

            histo = copy.deepcopy(flow_maker_syst.hV2vsPt)

            v2_dict[var].append(histo)

            del flow_maker_syst

        dict_hist_abs_syst[var] = histo_v2_syst

    # pt-dependent cuts
    print(f"var: fAvgItsClusSizeCosLambda")
    var = "fAvgItsClusSizeCosLambda"
    var_dir = cent_dirs[i_cent].mkdir(var)

    v2_dict[var] = []
    canvas_dict[var] = ROOT.TCanvas(
        f"c{var}_cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}",
        f"c{var}_cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}",
        800,
        600,
    )
    legend_dict[var] = ROOT.TLegend(0.16, 0.61, 0.86, 0.88, "", "brNDC")
    legend_dict[var].SetBorderSize(0)

    histo_v2_syst = []
    n_pt_bins = len(pt_bins[i_cent]) - 1
    for i_pt in range(0, n_pt_bins):
        histo_v2_syst.append(
            ROOT.TH1F(
                f"hV2syst_{var}_cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}_pt{i_pt}",
                ";v_{2}",
                20,
                default_v2_values[i_cent][i_pt][0]
                - 3 * default_v2_values[i_cent][i_pt][1],
                default_v2_values[i_cent][i_pt][0]
                + 3 * default_v2_values[i_cent][i_pt][1],
            )
        )

    for i_cut in range(n_ptdep_variations):

        print(f"{var}: {i_cut} / {n_ptdep_variations}")

        output_dir_varied = var_dir.mkdir(f"{i_cut}")

        flow_maker_syst = FlowMaker()
        flow_maker_syst.data_df = complete_df
        flow_maker_syst.tof_analysis = tof_analysis
        flow_maker_syst.pt_bins = pt_bins[i_cent]
        flow_maker_syst.cent_limits = centrality_classes[i_cent]
        flow_maker_syst.resolution = resolutions[i_cent]

        var_suffix = f"_{var}_{i_cut}"
        flow_maker_syst.selection_string = standard_selections
        flow_maker_syst.ptdep_selection_dict = ptdep_cut_string_dict_list[i_cut]
        flow_maker_syst.suffix = var_suffix

        flow_maker_syst.create_histograms()
        flow_maker_syst.make_flow()
        flow_values = flow_maker_syst.getFlowValues()
        # write histograms to file
        flow_maker_syst.output_dir = output_dir_varied
        flow_maker_syst.dump_to_output_file()

        for i_pt in range(0, flow_maker_syst.nPtBins()):
            histo_v2_syst[i_pt].Fill(flow_values[i_pt][0])

        histo = copy.deepcopy(flow_maker_syst.hV2vsPt)

        v2_dict[var].append(histo)

        del flow_maker_syst

    dict_hist_abs_syst[var] = histo_v2_syst

    # get color palette
    cols = ROOT.TColor.GetPalette()

    cent_dirs[i_cent].cd()
    cent_dirs[i_cent].mkdir("std")
    default_v2_histos[i_cent].Write()

    cent_syst_dir_name = (
        output_dir_name
        + f"/syst_plots/cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}/"
    )
    if not os.path.exists(cent_syst_dir_name):
        os.makedirs(cent_syst_dir_name)

    for var, histos in v2_dict.items():
        cent_dirs[i_cent].cd(f"{var}")
        canvas_dict[var].cd()
        canvas_title = f"{var}" + r";#it{p}_{T} (GeV/#it{c}); v_{2}"
        canvas_dict[var].DrawFrame(1.7, -0.2, 9.0, 1.0, canvas_title)
        period = int(cols.GetSize() / len(histos))
        for i_histo, histo in enumerate(histos):
            utils.setHistStyle(histo, cols.At(i_histo * period))
            legend_dict[var].AddEntry(histo, f"{cut_label_dict[var][i_histo]}", "PE")
            histo.Draw("PE SAME")
        legend_dict[var].AddEntry(default_v2_histos[i_cent], "std", "PE")
        default_v2_histos[i_cent].Draw("PE SAME")
        legend_dict[var].Draw()
        legend_dict[var].SetNColumns(5)
        canvas_dict[var].Write()
        canvas_dict[var].SaveAs(
            f"{cent_syst_dir_name}/{canvas_dict[var].GetName()}.pdf"
        )

        histo_abs_syst = default_v2_histos[i_cent].Clone(f"hAbsSystVsPt_{var}")
        histo_abs_syst.Reset()
        histo_abs_syst.GetYaxis().SetTitle("systematic error")
        histo_abs_syst.GetYaxis().SetRangeUser(0.0, 0.05)

        # getting systematic distributions
        for i_histo, histo in enumerate(dict_hist_abs_syst[var]):

            utils.setHistStyle(histo, ROOT.kRed + 1)
            # create canvas with the central value
            histo_name = histo.GetName()
            canvas_syst_name = histo_name.replace("h", "c", 1)
            canvas_syst = ROOT.TCanvas(canvas_syst_name, canvas_syst_name, 800, 600)
            canvas_syst.SetBottomMargin(0.13)
            canvas_syst.SetLeftMargin(0.13)

            std_line = ROOT.TLine(
                default_v2_values[i_cent][i_histo][0],
                0,
                default_v2_values[i_cent][i_histo][0],
                1.05 * histo.GetMaximum(),
            )
            std_line.SetLineColor(ROOT.kAzure + 2)
            std_line.SetLineWidth(2)
            # create box for statistical uncertainty
            std_errorbox = ROOT.TBox(
                default_v2_values[i_cent][i_histo][0]
                - default_v2_values[i_cent][i_histo][1],
                0,
                default_v2_values[i_cent][i_histo][0]
                + default_v2_values[i_cent][i_histo][1],
                1.05 * histo.GetMaximum(),
            )
            std_errorbox.SetFillColorAlpha(ROOT.kAzure + 1, 0.5)
            std_errorbox.SetLineWidth(0)

            info_panel = ROOT.TPaveText(0.6, 0.6, 0.8, 0.82, "NDC")
            info_panel.SetBorderSize(0)
            info_panel.SetFillStyle(0)
            info_panel.SetTextAlign(12)
            info_panel.SetTextFont(42)
            info_panel.AddText(r"PbPb, #sqrt{#it{s}_{nn}} = 5.36 TeV")
            info_panel.AddText(
                f"{centrality_classes[i_cent][0]} - {centrality_classes[i_cent][1]} % {cent_detector_label}"
            )
            pt_label = (
                f"{pt_bins[i_cent][i_pt]:.1f}"
                + r" #leq #it{p}_{T} < "
                + f"{pt_bins[i_cent][i_pt+1]:.1f}"
                + r" GeV/#it{c}"
            )
            info_panel.AddText(pt_label)

            canvas_syst.cd()
            histo.Draw("histo")
            std_errorbox.Draw()
            std_line.Draw()
            info_panel.Draw()

            cent_dirs[i_cent].cd(f"{var}")
            histo.Write()
            canvas_syst.Write()
            canvas_syst.SaveAs(f"{cent_syst_dir_name}/{canvas_syst.GetName()}.pdf")

            # absolute syst vs pt
            if histo.Integral() > 2:
                histo_abs_syst.SetBinContent(i_histo + 1, histo.GetRMS())
                histo_abs_syst.SetBinError(i_histo + 1, 0)

        cent_dirs[i_cent].cd(f"{var}")
        histo_abs_syst.Smooth(5)
        histo_abs_syst.Write()

if not do_alternative_systematic:
    exit(1)

# evaluate systematic from different table

for i_cent in range(n_cent_classes):

    cent_dir_name = (
        f"cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}"
    )
    table_comp_dir = output_file.mkdir(f"{cent_dir_name}/table_comp")
    cent_syst_dir_name = (
        output_dir_name
        + f"/syst_plots/cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}/"
    )

    histo_ep = default_v2_histos[i_cent]
    histo_qvec = alternative_v2_histos[i_cent]

    comp_table_canvas = ROOT.TCanvas(
        f"cCompTableCanvas_{cent_dir_name}",
        f"cCompTableCanvas_{cent_dir_name}",
        800,
        600,
    )

    utils.setHistStyle(histo_ep, ROOT.kRed)
    utils.setHistStyle(histo_qvec, ROOT.kBlack, marker=21)

    histo_ep.Draw("PE")
    histo_qvec.Draw("PE SAME")

    legend_comp_table = ROOT.TLegend(
        0.15,
        0.62,
        0.41,
        0.85,
        f"FT0C {centrality_classes[i_cent][0]} - {centrality_classes[i_cent][1]}%",
        "brNDC",
    )
    legend_comp_table.SetBorderSize(0)
    legend_comp_table.AddEntry(histo_ep, "EP", "PE")
    legend_comp_table.AddEntry(histo_qvec, "Qvector", "PE")
    legend_comp_table.Draw()

    cent_dirs[i_cent].cd("table_comp")
    comp_table_canvas.Write()

    comp_table_canvas.SaveAs(f"{cent_syst_dir_name}{comp_table_canvas.GetName()}.pdf")

    histo_table_diff = histo_ep.Clone(f"hV2tableDiff_{cent_dir_name}")
    histo_table_diff.Reset()
    histo_table_diff.GetYaxis().SetTitle(r"v_{2}^{EP} - v_{2}^{Qvec}")

    n_pt_bins = len(pt_bins[i_cent]) - 1

    for i_pt in range(1, n_pt_bins + 1):
        diff = histo_ep.GetBinContent(i_pt) - histo_qvec.GetBinContent(i_pt)
        diff_err = np.hypot(histo_ep.GetBinError(i_pt), histo_qvec.GetBinError(i_pt))
        histo_table_diff.SetBinContent(i_pt, abs(diff))
        histo_table_diff.SetBinError(i_pt, diff_err)

    histo_table_diff.Smooth(3)
    cent_dirs[i_cent].cd("table_comp")
    histo_table_diff.Write()

    utils.saveCanvasAsPDF(histo_table_diff, cent_syst_dir_name)
