import ROOT
import yaml
import os
import numpy as np
import pandas as pd

import sys

sys.path.append("utils")
import utils as utils

# setting configurable for eptables
config_file_eptables = open("config/config_analysis_pass4_new.yaml", "r")

config_eptables = yaml.full_load(config_file_eptables)

# setting configurable for qvectortables
config_file_qvectortables = open(
    "config/config_analysis_alternative_pass4_new.yaml", "r"
)

config_qvectortables = yaml.full_load(config_file_qvectortables)

# get QC file
qc_file_eptables = ROOT.TFile(
    config_eptables["output_dir_name"] + config_eptables["output_file_name_qc"]
)

qc_file_qvectortables = ROOT.TFile(
    config_qvectortables["output_dir_name"]
    + config_qvectortables["output_file_name_qc"]
)

cent_detector_label = config_eptables["cent_detector_label"]
reference_flow_detector = config_eptables["reference_flow_detector"]
resolution_flow_detectors = config_eptables["resolution_flow_detectors"]

centrality_classes = config_eptables["centrality_classes"]
pt_bins = config_eptables["qc_pt_bins"]
n_pt_bins = len(pt_bins) - 1
cent_colours = config_eptables["cent_colours"]

output_file = ROOT.TFile(
    config_eptables["output_dir_name"] + "table_compare_qc.root", "RECREATE"
)

detectors = ["FT0A", "FT0C", "TPCl", "TPCr", "TPC"]

for i_cent, cent in enumerate(centrality_classes):
    cent_dir = output_file.mkdir(f"{cent[0]}_{cent[1]}")
    for det in detectors:
        cent_dir.mkdir(det)

    for i_pt in range(0, n_pt_bins):
        # select the correct pt bin
        pt_label = (
            f"{pt_bins[i_pt]}"
            + r" #leq #it{p}_{T} < "
            + f"{pt_bins[i_pt + 1]}"
            + r" (GeV/#it{c})"
        )

        for det in detectors:

            hPhiMinusPsi_name = f"hPhiMinusPsi_{det}_cent_{cent[0]}_{cent[1]}_pt_{i_pt}"
            hPhiMinusPsi_eptables = qc_file_eptables.Get(
                f"{cent[0]}_{cent[1]}/{det}/{hPhiMinusPsi_name}"
            )
            hPhiMinusPsi_eptables.SetDirectory(0)
            hPhiMinusPsi_eptables.SetName(
                f"hPhiMinusPsi_eptables_{det}_cent_{cent[0]}_{cent[1]}_pt_{i_pt}"
            )
            utils.setHistStyle(hPhiMinusPsi_eptables, ROOT.kBlue + 2)
            hPhiMinusPsi_eptables.Rebin(5)

            hPhiMinusPsi_qvectortables = qc_file_qvectortables.Get(
                f"{cent[0]}_{cent[1]}/{det}/{hPhiMinusPsi_name}"
            )
            hPhiMinusPsi_qvectortables.SetDirectory(0)
            hPhiMinusPsi_qvectortables.SetName(
                f"hPhiMinusPsi_qvectortables_{det}_cent_{cent[0]}_{cent[1]}_pt_{i_pt}"
            )
            utils.setHistStyle(hPhiMinusPsi_qvectortables, ROOT.kRed + 2)
            hPhiMinusPsi_qvectortables.Rebin(5)

            hPhiMinusPsi_eptables_integral = hPhiMinusPsi_eptables.Integral()
            hPhiMinusPsi_qvectortables_integral = hPhiMinusPsi_qvectortables.Integral()
            if (
                hPhiMinusPsi_eptables_integral > 0
                and hPhiMinusPsi_qvectortables_integral > 0
            ):
                hPhiMinusPsi_eptables.Scale(1 / hPhiMinusPsi_eptables_integral)
                hPhiMinusPsi_qvectortables.Scale(
                    1 / hPhiMinusPsi_qvectortables_integral
                )

            yMaxPhiMinusPsi = hPhiMinusPsi_eptables.GetMaximum()
            if hPhiMinusPsi_qvectortables.GetMaximum() > yMaxPhiMinusPsi:
                yMaxPhiMinusPsi = hPhiMinusPsi_qvectortables.GetMaximum()
            hPhiMinusPsi_eptables.GetYaxis().SetRangeUser(0, yMaxPhiMinusPsi * 1.2)
            hPhiMinusPsi_qvectortables.GetYaxis().SetRangeUser(0, yMaxPhiMinusPsi * 1.2)

            cPhiMinusPsi = ROOT.TCanvas(
                f"2PhiMinusPsi_{det}_cent_{cent[0]}_{cent[1]}_pt_{i_pt}",
                f"cPhiMinusPsi_{det}_cent_{cent[0]}_{cent[1]}_pt_{i_pt}",
                800,
                600,
            )
            cPhiMinusPsi.cd()
            hPhiMinusPsi_eptables.Draw("hist")
            hPhiMinusPsi_qvectortables.Draw("hist same")
            lPhiMinusPsi = ROOT.TLegend(
                0.6, 0.7, 0.9, 0.9, f"{det}, {pt_label}", "brNDC"
            )
            lPhiMinusPsi.AddEntry(
                hPhiMinusPsi_eptables,
                f"eptables",
                "l",
            )
            lPhiMinusPsi.AddEntry(
                hPhiMinusPsi_qvectortables,
                f"qvectortables",
                "l",
            )
            lPhiMinusPsi.Draw()
            cent_dir.cd(f"{det}")
            cPhiMinusPsi.Write()
            hPhiMinusPsi_eptables.Write()
            hPhiMinusPsi_qvectortables.Write()

            hQx_name = f"hQx_{det}_cent_{cent[0]}_{cent[1]}_pt_{i_pt}"
            hQx_eptables = qc_file_eptables.Get(f"{cent[0]}_{cent[1]}/{det}/{hQx_name}")
            hQx_eptables.SetDirectory(0)
            hQx_eptables.SetName(
                f"hQx_eptables_{det}_cent_{cent[0]}_{cent[1]}_pt_{i_pt}"
            )
            utils.setHistStyle(hQx_eptables, ROOT.kBlue + 2)
            hQx_eptables.Rebin(5)

            hQx_qvectortables = qc_file_qvectortables.Get(
                f"{cent[0]}_{cent[1]}/{det}/{hQx_name}"
            )
            hQx_qvectortables.SetDirectory(0)
            hQx_qvectortables.SetName(
                f"hQx_qvectortables_{det}_cent_{cent[0]}_{cent[1]}_pt_{i_pt}"
            )
            utils.setHistStyle(hQx_qvectortables, ROOT.kRed + 2)
            hQx_qvectortables.Rebin(5)

            hQx_eptables_integral = hQx_eptables.Integral()
            hQx_qvectortables_integral = hQx_qvectortables.Integral()
            if hQx_eptables_integral > 0 and hQx_qvectortables_integral > 0:
                hQx_eptables.Scale(1 / hQx_eptables_integral)
                hQx_qvectortables.Scale(1 / hQx_qvectortables_integral)

            yMaxQx = hQx_eptables.GetMaximum()
            if hQx_qvectortables.GetMaximum() > yMaxQx:
                yMaxQx = hQx_qvectortables.GetMaximum()
            hQx_eptables.GetYaxis().SetRangeUser(0, yMaxQx * 1.2)
            hQx_qvectortables.GetYaxis().SetRangeUser(0, yMaxQx * 1.2)

            cQx = ROOT.TCanvas(
                f"cQx_{det}_cent_{cent[0]}_{cent[1]}_pt_{i_pt}",
                f"cQx_{det}_cent_{cent[0]}_{cent[1]}_pt_{i_pt}",
                800,
                600,
            )
            cQx.cd()
            hQx_eptables.Draw("hist")
            hQx_qvectortables.Draw("hist same")
            lQx = ROOT.TLegend(0.6, 0.7, 0.9, 0.9, f"{det}, {pt_label}", "brNDC")
            lQx.AddEntry(
                hQx_eptables,
                f"eptables",
                "l",
            )
            lQx.AddEntry(
                hQx_qvectortables,
                f"qvectortables",
                "l",
            )
            lQx.Draw()
            cent_dir.cd(f"{det}")
            cQx.Write()
            hQx_eptables.Write()
            hQx_qvectortables.Write()

            hQy_name = f"hQy_{det}_cent_{cent[0]}_{cent[1]}_pt_{i_pt}"
            hQy_eptables = qc_file_eptables.Get(f"{cent[0]}_{cent[1]}/{det}/{hQy_name}")
            hQy_eptables.SetDirectory(0)
            hQy_eptables.SetName(
                f"hQy_eptables_{det}_cent_{cent[0]}_{cent[1]}_pt_{i_pt}"
            )
            utils.setHistStyle(hQy_eptables, ROOT.kBlue + 2)
            hQy_eptables.Rebin(5)

            hQy_qvectortables = qc_file_qvectortables.Get(
                f"{cent[0]}_{cent[1]}/{det}/{hQy_name}"
            )
            hQy_qvectortables.SetDirectory(0)
            hQy_qvectortables.SetName(
                f"hQy_qvectortables_{det}_cent_{cent[0]}_{cent[1]}_pt_{i_pt}"
            )
            utils.setHistStyle(hQy_qvectortables, ROOT.kRed + 2)
            hQy_qvectortables.Rebin(5)

            hQy_eptables_integral = hQy_eptables.Integral()
            hQy_qvectortables_integral = hQy_qvectortables.Integral()
            if hQy_eptables_integral > 0 and hQy_qvectortables_integral > 0:
                hQy_eptables.Scale(1 / hQy_eptables_integral)
                hQy_qvectortables.Scale(1 / hQy_qvectortables_integral)

            yMaxQy = hQy_eptables.GetMaximum()
            if hQy_qvectortables.GetMaximum() > yMaxQy:
                yMaxQy = hQy_qvectortables.GetMaximum()
            hQy_eptables.GetYaxis().SetRangeUser(0, yMaxQy * 1.2)
            hQy_qvectortables.GetYaxis().SetRangeUser(0, yMaxQy * 1.2)

            cQy = ROOT.TCanvas(
                f"cQy_{det}_cent_{cent[0]}_{cent[1]}_pt_{i_pt}",
                f"cQy_{det}_cent_{cent[0]}_{cent[1]}_pt_{i_pt}",
                800,
                600,
            )
            cQy.cd()
            hQy_eptables.Draw("hist")
            hQy_qvectortables.Draw("hist same")
            lQy = ROOT.TLegend(0.6, 0.7, 0.9, 0.9, f"{det}, {pt_label}", "brNDC")
            lQy.AddEntry(
                hQy_eptables,
                f"eptables",
                "l",
            )
            lQy.AddEntry(
                hQy_qvectortables,
                f"qvectortables",
                "l",
            )
            lQy.Draw()
            cent_dir.cd(f"{det}")
            cQy.Write()
            hQy_eptables.Write()
            hQy_qvectortables.Write()

            hQ_name = f"hQ_{det}_cent_{cent[0]}_{cent[1]}_pt_{i_pt}"
            hQ_eptables = qc_file_eptables.Get(f"{cent[0]}_{cent[1]}/{det}/{hQ_name}")
            hQ_eptables.SetDirectory(0)
            hQ_eptables.SetName(f"hQ_eptables_{det}_cent_{cent[0]}_{cent[1]}_pt_{i_pt}")
            utils.setHistStyle(hQ_eptables, ROOT.kBlue + 2)
            hQ_eptables.Rebin(5)

            hQ_qvectortables = qc_file_qvectortables.Get(
                f"{cent[0]}_{cent[1]}/{det}/{hQ_name}"
            )
            hQ_qvectortables.SetDirectory(0)
            hQ_qvectortables.SetName(
                f"hQ_qvectortables_{det}_cent_{cent[0]}_{cent[1]}_pt_{i_pt}"
            )
            utils.setHistStyle(hQ_qvectortables, ROOT.kRed + 2)
            hQ_eptables.Rebin(5)

            hQ_eptables_integral = hQ_eptables.Integral()
            hQ_qvectortables_integral = hQ_qvectortables.Integral()
            if hQ_eptables_integral > 0 and hQ_qvectortables_integral > 0:
                hQ_eptables.Scale(1 / hQ_eptables_integral)
                hQ_qvectortables.Scale(1 / hQ_qvectortables_integral)

            yMaxQ = hQ_eptables.GetMaximum()
            if hQ_qvectortables.GetMaximum() > yMaxQ:
                yMaxQ = hQ_qvectortables.GetMaximum()
            hQ_eptables.GetYaxis().SetRangeUser(0, yMaxQ * 1.2)
            hQ_qvectortables.GetYaxis().SetRangeUser(0, yMaxQ * 1.2)

            cQ = ROOT.TCanvas(
                f"cQ_{det}_cent_{cent[0]}_{cent[1]}_pt_{i_pt}",
                f"cQ_{det}_cent_{cent[0]}_{cent[1]}_pt_{i_pt}",
                800,
                600,
            )
            cQ.cd()
            hQ_eptables.Draw("hist")
            hQ_qvectortables.Draw("hist same")
            lQ = ROOT.TLegend(0.6, 0.7, 0.9, 0.9, f"{det}, {pt_label}", "brNDC")
            lQ.AddEntry(
                hQ_eptables,
                f"eptables",
                "l",
            )
            lQ.AddEntry(
                hQ_qvectortables,
                f"qvectortables",
                "l",
            )
            lQ.Draw()
            cent_dir.cd(f"{det}")
            cQ.Write()
            hQ_eptables.Write()
            hQ_qvectortables.Write()

output_file.cd()

resolution_file_eptables = ROOT.TFile("../results_pass4_new/resolution_EP.root")

hResolutionEP_eptables = resolution_file_eptables.Get(
    "Resolution_EP/hResolution_FT0C_FT0A_TPC_EP"
)
hResolutionEP_eptables.SetDirectory(0)
hResolutionEP_eptables.SetName("hResolution_FT0C_FT0A_TPC_EP_eptables")

hResolutionSP_eptables = resolution_file_eptables.Get(
    "Resolution_SP/hResolution_FT0C_FT0A_TPC_SP"
)
hResolutionSP_eptables.SetDirectory(0)
hResolutionSP_eptables.SetName("hResolution_FT0C_FT0A_TPC_SP_eptables")

resolution_file_qvectortables = ROOT.TFile("../results_pass4_new/resolution_Qvec.root")
hResolutionEP_qvectortables = resolution_file_qvectortables.Get(
    "Resolution_EP/hResolution_FT0C_FT0A_TPC_EP"
)
hResolutionEP_qvectortables.SetDirectory(0)
hResolutionEP_qvectortables.SetName("hResolution_FT0C_FT0A_TPC_EP_qvectortables")
utils.setHistStyle(hResolutionEP_qvectortables, ROOT.kBlack)

hResolutionSP_qvectortables = resolution_file_qvectortables.Get(
    "Resolution_SP/hResolution_FT0C_FT0A_TPC_SP"
)
hResolutionSP_qvectortables.SetDirectory(0)
hResolutionSP_qvectortables.SetName("hResolution_FT0C_FT0A_TPC_SP_qvectortables")
utils.setHistStyle(hResolutionSP_qvectortables, ROOT.kBlack)

output_file.cd()
hResolutionEP_eptables.Write()
hResolutionEP_qvectortables.Write()
hResolutionSP_eptables.Write()
hResolutionSP_qvectortables.Write()

output_file.cd()
eptables_final_comp_dir = output_file.mkdir("eptables_final_comp")
qvectortables_final_comp_dir = output_file.mkdir("qvectortables_final_comp")

input_eptables_compare_file = ROOT.TFile("../results_pass4_new/v2_comp.root")
input_qvectortables_compare_file = ROOT.TFile("../results_pass4_new/v2_comp.root")

for i_cent, cent in enumerate(centrality_classes):
    eptables_canvas = input_eptables_compare_file.Get(
        f"cCompSP_cent_{cent[0]}_{cent[1]}"
    )
    eptables_final_comp_dir.cd()
    eptables_canvas.Write()
    qvectortables_canvas = input_qvectortables_compare_file.Get(
        f"cCompSP_cent_{cent[0]}_{cent[1]}"
    )
    qvectortables_final_comp_dir.cd()
    qvectortables_canvas.Write()
