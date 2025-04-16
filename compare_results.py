import ROOT
import yaml
import argparse
import numpy as np

import sys

sys.path.append("utils")
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

resolution_file_name = config["resolution_file_name"]

output_dir_name = config["output_dir_name"]

resolition_file = ROOT.TFile(resolution_file_name)

check_igor = False

# different resolution detectors

res1 = resolition_file.Get("Resolution_EP/hResolution_FT0C_TPCl_TPCr_EP")
utils.setHistStyle(res1, ROOT.kBlack)
res1.GetYaxis().SetTitle(r"R_{2} (FT0C)")

res2 = resolition_file.Get("Resolution_EP/hResolution_FT0C_FT0A_TPC_EP")
utils.setHistStyle(res2, ROOT.kRed)
res2.GetYaxis().SetTitle(r"R_{2} (FT0C)")


hRatioRes = res1.Clone("hRatioRes")
hRatioRes.GetYaxis().SetTitle("ratio")
hRatioRes.GetYaxis().SetRangeUser(0.5, 1.5)
hRatioRes.Divide(res2)

line_cent = ROOT.TLine(0.0, 1.0, 100.0, 1.0)
line_cent.SetLineStyle(ROOT.kDashed)

cCompDetRes = ROOT.TCanvas("cCompDetRes", "cCompDetRes", 800, 600)

pad_top_res = ROOT.TPad("pad_top_res", "pad_top_res", 0.0, 0.5, 1.0, 1.0, 0)
pad_top_res.SetLeftMargin(0.15)
pad_top_res.SetBottomMargin(0.0)
pad_top_res.Draw()
pad_top_res.cd()

res1.Draw("PE")
res2.Draw("PE SAME")

legendRes = ROOT.TLegend(0.72, 0.49, 0.89, 0.67, "", "brNDC")
legendRes.AddEntry(res1, "TPCl, TPCr")
legendRes.AddEntry(res2, "FT0A, TPC")
legendRes.Draw()

cCompDetRes.cd()
pad_bottom_res = ROOT.TPad("pad_bottom_res", "pad_bottom_res", 0.0, 0.0, 1.0, 0.5, 0)
pad_bottom_res.SetLeftMargin(0.15)
pad_bottom_res.SetTopMargin(0.0)
pad_bottom_res.SetBottomMargin(0.3)
pad_bottom_res.Draw()
pad_bottom_res.cd()
hRatioRes.Draw("PE")
line_cent.Draw()


output_file = ROOT.TFile(f"{output_dir_name}/v2_comp.root", "RECREATE")
res_dir = output_file.mkdir("Resolution_Detectors")
res_dir.cd()
res1.Write()
res2.Write()
cCompDetRes.Write()
hRatioRes.Write()

cCompDetRes.SaveAs(f"{output_dir_name}/syst_plots/comparisons/v2_comp_resolutions.pdf")

centrality_classes = config["centrality_classes"]
n_cent_classes = len(centrality_classes)
pt_bins = config["pt_bins"]

input_file_name_default = config["output_dir_name"] + config["output_file_name"]
input_file_default = ROOT.TFile(input_file_name_default)

print(f"Input file: {input_file_name_default}")

# Igor's check

if check_igor:
    input_file_name_igor = (
        config["output_dir_name"][:-1] + "_igor/" + config["output_file_name"]
    )
    input_file_igor = ROOT.TFile(input_file_name_igor)

    for i_cent in range(n_cent_classes):
        print(f"Processing centrality class {i_cent}.")
        histo_default = input_file_default.Get(
            f"cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}/default/hV2vsPt_cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}"
        )
        histo_default.SetDirectory(0)

        histo_igor = input_file_igor.Get(
            f"cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}/default/hV2vsPt_cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}"
        )
        histo_igor.SetDirectory(0)
        histo_igor.SetName(
            f"hV2vsPt_cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}_igor"
        )
        utils.setHistStyle(histo_igor, ROOT.kGray + 1, 21)

        hRatioIgor = histo_default.Clone(
            f"hRatioIgor_cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}"
        )
        hRatioIgor.Reset()
        hRatioIgor.GetYaxis().SetTitle("ratio")
        hRatioIgor.GetYaxis().SetRangeUser(0.5, 1.5)
        for i_pt in range(1, len(pt_bins[i_cent])):
            val1 = histo_default.GetBinContent(i_pt)
            err1 = histo_default.GetBinError(i_pt)
            val2 = histo_igor.GetBinContent(i_pt)
            err2 = histo_igor.GetBinError(i_pt)
            ratio_val = val1 / val2
            ratio_err = ratio_val * np.sqrt((err1 / val1) ** 2 + (err2 / val2) ** 2)
            hRatioIgor.SetBinError(i_pt, ratio_err)
            hRatioIgor.SetBinContent(i_pt, ratio_val)

        pt_limits = [pt_bins[i_cent][0], pt_bins[i_cent][-1]]

        line = ROOT.TLine(pt_limits[0], 1.0, pt_limits[1], 1.0)
        line.SetLineStyle(ROOT.kDashed)

        cCompIgor = ROOT.TCanvas(
            f"cCompDetIgor_cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}",
            f"cCompDetIgor_cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}",
            800,
            600,
        )

        pad_top_igor = ROOT.TPad(
            f"pad_top_igor_cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}",
            f"pad_top_igor_cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}",
            0.0,
            0.5,
            1.0,
            1.0,
            0,
        )
        pad_top_igor.SetLeftMargin(0.15)
        pad_top_igor.SetBottomMargin(0.0)
        pad_top_igor.Draw()
        pad_top_igor.cd()

        histo_default.Draw("PE")
        histo_igor.Draw("PE SAME")

        legendIgor = ROOT.TLegend(
            0.20,
            0.58,
            0.37,
            0.76,
            f"{centrality_classes[i_cent][0]} - {centrality_classes[i_cent][1]} #%",
            "brNDC",
        )
        legendIgor.AddEntry(histo_default, "Standard selections")
        legendIgor.AddEntry(histo_igor, "isGoodITSLayersAll")
        legendIgor.Draw()

        cCompIgor.cd()
        pad_bottom_igor = ROOT.TPad(
            "pad_bottom_igor", "pad_bottom_igor", 0.0, 0.0, 1.0, 0.5, 0
        )
        pad_bottom_igor.SetLeftMargin(0.15)
        pad_bottom_igor.SetTopMargin(0.0)
        pad_bottom_igor.SetBottomMargin(0.3)
        pad_bottom_igor.Draw()
        pad_bottom_igor.cd()
        hRatioIgor.Draw("PE")
        line.Draw()

        output_file.cd()
        histo_default.Write()
        histo_default.Write()
        cCompIgor.Write()
        hRatioIgor.Write()

        cCompIgor.SaveAs(
            f"{output_dir_name}/syst_plots/comparisons/v2_comp_igor_cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}.pdf"
        )

# SP check

input_file_name_sp = (
    config["output_dir_name"][:-1] + "_SP/" + config["output_file_name"]
)
input_file_sp = ROOT.TFile(input_file_name_sp)
print(f"Input file: {input_file_name_sp}")
print(f"Output file: {output_dir_name}/v2_comp.root")


for i_cent in range(n_cent_classes):
    print(f"Processing centrality class {i_cent}.")
    histo_default = input_file_default.Get(
        f"cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}/default/hV2vsPt_cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}"
    )
    histo_default.SetDirectory(0)

    histo_sp = input_file_sp.Get(
        f"cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}/default/hV2vsPt_cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}"
    )
    histo_sp.SetDirectory(0)
    histo_sp.SetName(
        f"hV2vsPt_cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}_sp"
    )
    utils.setHistStyle(histo_sp, ROOT.kGray + 1, 21)

    hRatioSP = histo_default.Clone(
        f"hRatioSP_cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}"
    )
    hRatioSP.Reset()
    hRatioSP.GetYaxis().SetTitle("EP / SP")
    hRatioSP.GetYaxis().SetRangeUser(0.12, 1.88)
    for i_pt in range(1, len(pt_bins[i_cent])):
        val1 = histo_default.GetBinContent(i_pt)
        err1 = histo_default.GetBinError(i_pt)
        val2 = histo_sp.GetBinContent(i_pt)
        err2 = histo_sp.GetBinError(i_pt)
        ratio_val = val1 / val2
        ratio_err = ratio_val * np.sqrt((err1 / val1) ** 2 + (err2 / val2) ** 2)
        hRatioSP.SetBinError(i_pt, ratio_err)
        hRatioSP.SetBinContent(i_pt, ratio_val)

    pt_limits = [pt_bins[i_cent][0], pt_bins[i_cent][-1]]

    line_SP = ROOT.TLine(pt_limits[0], 1.0, pt_limits[1], 1.0)
    line_SP.SetLineStyle(ROOT.kDashed)

    hRatioSP.Fit("pol0", "Q")

    cCompSP = ROOT.TCanvas(
        f"cCompSP_cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}",
        f"cCompSP_cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}",
        800,
        600,
    )

    pad_top_sp = ROOT.TPad(
        f"pad_top_sp_cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}",
        f"pad_top_sp_cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}",
        0.0,
        0.5,
        1.0,
        1.0,
        0,
    )
    pad_top_sp.SetLeftMargin(0.15)
    pad_top_sp.SetBottomMargin(0.0)
    pad_top_sp.Draw()
    pad_top_sp.cd()

    histo_sp.Draw("PE")
    histo_default.Draw("PE SAME")

    legend_sp = ROOT.TLegend(
        0.20,
        0.58,
        0.37,
        0.76,
        f"{centrality_classes[i_cent][0]} - {centrality_classes[i_cent][1]} #%",
        "brNDC",
    )
    legend_sp.AddEntry(histo_default, "Event-Plane")
    legend_sp.AddEntry(histo_sp, "Scalar Product")
    legend_sp.Draw()

    cCompSP.cd()
    pad_bottom_sp = ROOT.TPad("pad_bottom_sp", "pad_bottom_sp", 0.0, 0.0, 1.0, 0.5, 0)
    pad_bottom_sp.SetLeftMargin(0.15)
    pad_bottom_sp.SetTopMargin(0.0)
    pad_bottom_sp.SetBottomMargin(0.3)
    pad_bottom_sp.Draw()
    pad_bottom_sp.cd()
    hRatioSP.Draw("PE")
    line_SP.Draw()

    output_file.cd()
    histo_default.Write()
    histo_default.Write()
    cCompSP.Write()
    hRatioSP.Write()

    cCompSP.SaveAs(
        f"{output_dir_name}/syst_plots/comparisons/v2_comp_sp_cent_{centrality_classes[i_cent][0]}_{centrality_classes[i_cent][1]}.pdf"
    )
    if i_cent == 0:
        cCompSP.SaveAs(f"{output_dir_name}/syst_plots/comparisons/v2_comp_sp_all.pdf[")
    cCompSP.SaveAs(f"{output_dir_name}/syst_plots/comparisons/v2_comp_sp_all.pdf")
cCompSP.SaveAs(f"{output_dir_name}/syst_plots/comparisons/v2_comp_sp_all.pdf]")
