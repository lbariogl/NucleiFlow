import ROOT
import pandas as pd
import yaml
import argparse
import numpy as np
from itertools import product
import os
import copy
import random

import sys

sys.path.append("utils")
from flow import FlowMaker
import utils as utils

parser = argparse.ArgumentParser(description="Configure the parameters of the script.")
parser.add_argument(
    "--config",
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

nuclei_tree_name = config["nuclei_tree_name"]
ep_tree_name = config["ep_tree_name"]

useSP = config["useSP"]

# n-sigma TPC bins
n_nsigmaTPC_bins = config["n_nsigmaTPC_bins"]
nsigmaTPC_bin_limits = config["nsigmaTPC_bin_limits"]

# Tof analysis (True for 4He)
tof_analysis = config["tof_analysis"]

mandatory_selections = config["mandatory_selections"]
selection_dict = config["selection_dict"]
selection_list = selection_dict.values()
selections = " and ".join(selection_list)

ptdep_selection_dict = config["ptdep_selection_dict"]["fAvgItsClusSizeCosLambda"]

cent_detector_label = config["cent_detector_label"]
reference_flow_detector = config["reference_flow_detector"]

centrality_classes = config["centrality_classes"]
pt_bins = config["pt_bins"]
cent_colours = config["cent_colours"]

do_syst = config["do_syst"]
use_Barlow = config["use_Barlow"]

# BB parameters
p_train = config["p_train"]
resolution_train = config["resolution_train"]
n_sigma_plot = config["n_sigma_plot"]

# create output file
if not os.path.exists(output_dir_name):
    os.makedirs(output_dir_name)
output_file = ROOT.TFile(f"{output_dir_name}/{output_file_name}", "recreate")


# get a unique df from nuclei and ep trees
nuclei_df = utils.get_df_from_tree(input_file_name, nuclei_tree_name)

nucleiflow_df = utils.get_df_from_tree(input_file_name, ep_tree_name)

complete_df = pd.concat([nuclei_df, nucleiflow_df], axis=1, join="inner")

# define new columns
utils.redefineColumns(
    complete_df, mass=utils.mass_helion, parameters=p_train, useSP=useSP
)

# apply mandatory selections
complete_df.query(mandatory_selections, inplace=True)

# Get resolution from file
resolution_file = ROOT.TFile(resolution_file_name)
if useSP:
    res_histo_name = "Resolution_SP/hResolution_FT0C_TPCl_TPCr_SP"
else:
    res_histo_name = "Resolution_EP/hResolution_FT0C_TPCl_TPCr_EP"
hResolution = resolution_file.Get(res_histo_name)
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

res_0_20 = (hResolution.GetBinContent(1) + hResolution.GetBinContent(2)) / 2
res_20_50 = (hResolution.GetBinContent(3) + hResolution.GetBinContent(5)) / 3
res_50_80 = (hResolution.GetBinContent(6) + hResolution.GetBinContent(8)) / 3

# resolutions = [res_0_10, res_10_20, res_20_40, res_40_60]
resolutions_3He = [
    res_0_10,
    res_10_20,
    res_20_30,
    res_30_40,
    res_40_50,
    res_50_60,
    res_60_80,
]

resolutions_4He = [res_0_20, res_20_50, res_50_80]

# Flow measurement
output_file.cd()

n_cent_classes = len(centrality_classes)

# Standard selections
flow_makers = []
default_values = []
cent_dirs = []
cent_plots_dir_names = []

for i_cent in range(n_cent_classes):
    print("*****************************")
    print(
        f" cent class: {centrality_classes[i_cent][0]} - {centrality_classes[i_cent][1]}"
    )
    print("*****************************")
    flow_maker = FlowMaker()
    flow_maker.data_df = complete_df
    flow_maker.selection_string = selections
    flow_maker.tof_analysis = tof_analysis
    flow_maker.n_nsigmaTPC_bins = n_nsigmaTPC_bins
    flow_maker.nsigmaTPC_bin_limits = nsigmaTPC_bin_limits
    flow_maker.pt_bins = pt_bins[i_cent]

    flow_maker.ptdep_selection_dict = ptdep_selection_dict

    flow_maker.cent_limits = centrality_classes[i_cent]
    if tof_analysis:
        flow_maker.resolution = resolutions_4He[i_cent]
    else:
        flow_maker.resolution = resolutions_3He[i_cent]
    flow_maker.print_frame = True
    flow_maker.ref_detector = reference_flow_detector

    # create output_dir
    cent_dir = output_file.mkdir(
        f"cent_{flow_maker.cent_limits[0]}_{flow_maker.cent_limits[1]}"
    )
    cent_dirs.append(cent_dir)
    default_dir = cent_dir.mkdir("default")

    flow_maker.output_dir = default_dir

    plot_dir_name = f"{output_dir_name}/plots/cent_{flow_maker.cent_limits[0]}_{flow_maker.cent_limits[1]}"
    if not os.path.exists(plot_dir_name):
        os.makedirs(plot_dir_name)

    cent_plots_dir_names.append(plot_dir_name)

    flow_maker.plot_dir = plot_dir_name
    flow_maker.color = cent_colours[i_cent]

    flow_maker.create_histograms()
    flow_maker.make_flow()
    flow_maker.dump_to_output_file()
    flow_maker.dump_to_pdf()
    flow_maker.dump_summary_to_pdf()

    flow_makers.append(flow_maker)
    default_values.append(flow_maker.getFlowValues())

# Systematic uncertainties
if do_syst:

    print("** Starting systematic variations **")
    n_trials = config["n_trials"]
    # output_dir_syst = output_file.mkdir('trials')

    # List of trial strings to be printed to a text file
    trial_strings = []
    print("----------------------------------")
    print("** Starting systematics analysis **")
    print(f"** {n_trials} trials will be tested **")
    print("----------------------------------")

    cut_dict_syst = config["cut_dict_syst"]
    ptdep_cut_dict_syst = config["ptdep_cut_dict_syst"]
    n_ptdep_variations_syst = config["n_ptdep_variations_syst"][
        "fAvgItsClusSizeCosLambda"
    ]

    # Create all possible variations (values) fot all the variables of interest (key)
    cut_string_dict = {}
    for var in cut_dict_syst:
        var_dict = cut_dict_syst[var]
        cut_greater = var_dict["cut_greater"]
        cut_greater_string = " > " if cut_greater else " < "
        cut_list = var_dict["cut_list"]
        cut_arr = np.linspace(cut_list[0], cut_list[1], cut_list[2])
        cut_string_dict[var] = []
        for cut in cut_arr:
            cut_string_dict[var].append(var + cut_greater_string + str(cut))

    pt_independent_selections = list(product(*list(cut_string_dict.values())))

    ptdep_cut_string_dict = {}
    n_max_combinations = 0
    n_pt_independent_selections = len(pt_independent_selections)
    n_pt_dependent_selections = n_ptdep_variations_syst

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
    for i_variation in range(n_pt_dependent_selections):
        variation_dict = {}
        for pt_key in ptdep_cut_dict_syst["fAvgItsClusSizeCosLambda"]:
            variation_dict[pt_key] = ptdep_cut_string_dict[pt_key][i_variation]
        ptdep_cut_string_dict_list.append(variation_dict)

    if n_pt_dependent_selections != 0:
        n_max_combinations = n_pt_independent_selections * n_pt_dependent_selections
        all_complete_selection_indices = list(
            product(
                range(n_pt_independent_selections), range(n_pt_dependent_selections)
            )
        )
    else:
        n_max_combinations = n_pt_independent_selections
        all_complete_selection_indices = list(
            product(range(n_pt_independent_selections), [-1])
        )

    complete_selection_random_indices = copy.deepcopy(all_complete_selection_indices)

    if n_trials < n_max_combinations:
        random.shuffle(complete_selection_random_indices)
        complete_selection_random_indices = complete_selection_random_indices[:n_trials]
    else:
        print(
            f"** Warning: n_trials > n_combinations ({n_trials}, {n_max_combinations}), taking all the possible combinations **"
        )

    # vary selections for each analysis_object, by centrality class
    for i_cent in range(n_cent_classes):

        histo_v2_syst = []
        n_pt_bins = len(pt_bins[i_cent]) - 1
        for i_pt in range(0, n_pt_bins):
            histo_v2_syst.append(
                ROOT.TH1F(
                    f"hV2syst_cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}_pt{i_pt}",
                    ";v_{2}",
                    20,
                    default_values[i_cent][i_pt][0]
                    - 3 * default_values[i_cent][i_pt][1],
                    default_values[i_cent][i_pt][0]
                    + 3 * default_values[i_cent][i_pt][1],
                )
            )

        for i_complete_selection, complete_selection_indices in enumerate(
            complete_selection_random_indices
        ):

            flow_maker_syst = FlowMaker()
            flow_maker_syst.silent_mode = True
            flow_maker_syst.data_df = complete_df
            flow_maker_syst.tof_analysis = tof_analysis
            flow_maker_syst.pt_bins = pt_bins[i_cent]
            flow_maker_syst.cent_limits = centrality_classes[i_cent]
            if tof_analysis:
                flow_maker_syst.resolution = resolutions_4He[i_cent]
            else:
                flow_maker_syst.resolution = resolutions_3He[i_cent]
            flow_maker_syst.ref_detector = reference_flow_detector

            complete_selection_suffix = f"_sel{i_complete_selection}"
            trial_strings.append("----------------------------------")
            trial_num_string = f"Trial: {i_complete_selection} / {len(complete_selection_random_indices)}"
            trial_strings.append(trial_num_string)
            print(trial_num_string)
            pt_independent_sel = pt_independent_selections[
                complete_selection_indices[0]
            ]
            sel_string = " and ".join(pt_independent_sel)
            trial_strings.append(str(sel_string))

            flow_maker_syst.selection_string = sel_string
            flow_maker_syst.ptdep_selection_dict = {}
            # ptdep_cut_string_dict_list[
            #     complete_selection_indices[1]
            # ]

            flow_maker_syst.suffix = complete_selection_suffix

            flow_maker_syst.create_histograms()
            flow_maker_syst.make_flow()
            flow_values = flow_maker_syst.getFlowValues()

            # write histograms to file
            trial_dir = cent_dirs[i_cent].mkdir(f"trial_{i_complete_selection}")
            flow_maker_syst.output_dir = trial_dir
            flow_maker_syst.dump_to_output_file()

            print(
                f"cent: {flow_maker_syst.cent_limits[0]} - {flow_maker_syst.cent_limits[1]}"
            )
            for i_pt in range(0, flow_maker_syst.nPtBins()):
                print(f"pt: {i_pt} / {flow_maker_syst.nPtBins()}")
                if use_Barlow:
                    if utils.passBarlow(
                        default_values[i_cent][i_pt][0],
                        flow_values[i_pt][0],
                        default_values[i_cent][i_pt][1],
                        flow_values[i_pt][1],
                    ):
                        histo_v2_syst[i_pt].Fill(flow_values[i_pt][0])
                    else:
                        print("    Rejected for Barlow")
                # print(f'pt_bin: {i_pt} -> v2: {flow_values[i_pt][0]} +- {flow_values[i_pt][1]}')
            print("----------------------------------")

            del flow_maker_syst
            del flow_values

        output_file.cd(
            f"cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}"
        )
        for i_pt in range(0, n_pt_bins):

            utils.setHistStyle(histo_v2_syst[i_pt], ROOT.kRed + 1)
            # create canvas with the central value
            histo_name = histo_v2_syst[i_pt].GetName()
            canvas_syst_name = histo_name.replace("h", "c", 1)
            canvas_syst = ROOT.TCanvas(canvas_syst_name, canvas_syst_name, 800, 600)
            canvas_syst.SetBottomMargin(0.13)
            canvas_syst.SetLeftMargin(0.13)

            std_line = ROOT.TLine(
                default_values[i_cent][i_pt][0],
                0,
                default_values[i_cent][i_pt][0],
                1.05 * histo_v2_syst[i_pt].GetMaximum(),
            )
            std_line.SetLineColor(ROOT.kAzure + 2)
            std_line.SetLineWidth(2)
            # create box for statistical uncertainty
            std_errorbox = ROOT.TBox(
                default_values[i_cent][i_pt][0] - default_values[i_cent][i_pt][1],
                0,
                default_values[i_cent][i_pt][0] + default_values[i_cent][i_pt][1],
                1.05 * histo_v2_syst[i_pt].GetMaximum(),
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
            histo_v2_syst[i_pt].Draw("histo")
            std_errorbox.Draw()
            std_line.Draw()
            info_panel.Draw()

            histo_v2_syst[i_pt].Write()
            canvas_syst.Write()
            canvas_syst.SaveAs(
                f"{cent_plots_dir_names[i_cent]}/{canvas_syst.GetName()}.pdf"
            )


# Final plots

print("Making final plots")
cV2 = ROOT.TCanvas("cV2", "cV2", 800, 600)
frame = cV2.DrawFrame(1.7, -0.1, 12.0, 1.1, r";#it{p}_{T} (GeV/#it{c}); v_{2}")
cV2.cd()
legend = ROOT.TLegend(0.66, 0.62, 0.92, 0.85, "FT0C centrality", "brNDC")
legend.SetBorderSize(0)
legend.SetNColumns(2)

for i_cent in range(n_cent_classes):
    cent_label = (
        f"{flow_makers[i_cent].cent_limits[0]} - {flow_makers[i_cent].cent_limits[1]}"
        + r"%"
    )
    legend.AddEntry(flow_makers[i_cent].hV2vsPt, cent_label, "PF")
    flow_makers[i_cent].hV2vsPt.Draw("PE SAME")

legend.Draw()

output_file.cd()
cV2.Write()
cV2.SaveAs(f"{output_dir_name}/plots/{cV2.GetName()}.pdf")

cPurity = ROOT.TCanvas("cPurity", "cPurity", 800, 600)
frame = cPurity.DrawFrame(1.7, 0.5, 12.0, 1.1, r";#it{p}_{T} (GeV/#it{c}); purity")
cPurity.cd()
legend_purity = ROOT.TLegend(0.23, 0.22, 0.49, 0.57, "FT0C centrality", "brNDC")
legend_purity.SetBorderSize(0)
legend_purity.SetNColumns(2)

for i_cent in range(n_cent_classes):
    cent_label = (
        f"{flow_makers[i_cent].cent_limits[0]} - {flow_makers[i_cent].cent_limits[1]}"
        + r"%"
    )
    legend_purity.AddEntry(flow_makers[i_cent].hPurityVsPt, cent_label, "PF")
    flow_makers[i_cent].hPurityVsPt.Draw("SAME")

legend_purity.Draw()

output_file.cd()
cPurity.Write()
cPurity.SaveAs(f"{output_dir_name}/plots/{cPurity.GetName()}.pdf")

cRawCounts = ROOT.TCanvas("cRawCounts", "cRawCounts", 800, 600)
frame = cRawCounts.DrawFrame(1.7, 0.5, 12.0, 25000, r";#it{p}_{T} (GeV/#it{c}); counts")
cRawCounts.cd()
legend_raw = ROOT.TLegend(0.61, 0.49, 0.87, 0.84, "FT0C centrality", "brNDC")
legend_raw.SetBorderSize(0)
legend_raw.SetNColumns(2)

cent_rebins = []
for cent in centrality_classes:
    cent_rebins.append(cent[0])
cent_rebins.append(centrality_classes[-1][1])
cent_rebins_arr = np.array(cent_rebins, dtype=np.float64)

hRawCounts = ROOT.TH1F(
    "hRawCounts", ";FTOC centrality (%); counts;", n_cent_classes, cent_rebins_arr
)
utils.setHistStyle(hRawCounts, ROOT.kRed)

for i_cent in range(n_cent_classes):
    cent_label = (
        f"{flow_makers[i_cent].cent_limits[0]} - {flow_makers[i_cent].cent_limits[1]}"
        + r"%"
    )
    legend_raw.AddEntry(flow_makers[i_cent].hRawCountsVsPt, cent_label, "PF")
    flow_makers[i_cent].hRawCountsVsPt.Draw("SAME")
    hRawCounts.SetBinContent(i_cent + 1, flow_makers[i_cent].hRawCountsVsPt.Integral())

legend_raw.Draw()

output_file.cd()
cRawCounts.Write()
cRawCounts.SaveAs(f"{output_dir_name}/plots/{cRawCounts.GetName()}.pdf")
hRawCounts.Write()
