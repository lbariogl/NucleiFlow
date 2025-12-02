import ROOT
import pandas as pd
import yaml
import argparse
import os
import re
from itertools import combinations
import numpy as np

import sys

sys.path.append("utils")
import utils as utils

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)

ROOT.ROOT.EnableImplicitMT()
ROOT.RooMsgService.instance().setSilentMode(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

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

centrality_classes = config["centrality_classes"]
pt_bins = config["pt_bins"]
cent_colours = config["cent_colours"]

# Crea la cartella per i PDF se non esiste
output_dir = "../results_phi_studies"
os.makedirs(output_dir, exist_ok=True)

# create output file
output_file = ROOT.TFile(f"{output_dir}/phi_studies/phi_studies.root", "recreate")
print(f"Output file: {output_file.GetName()}")

# get resolutions
res_histo_name = f"Resolution_EP/hResolution_FT0C_FT0A_TPC_EP"

resolution_file = ROOT.TFile(resolution_file_name, "read")
hResolution = resolution_file.Get(res_histo_name)
hResolution.SetDirectory(0)

resolution_v4fromv2_file_name = (
    "../results_pass4_new_cent_SP/resolution_EP_v2_for_v4.root"
)
resolution_v4fromv2_file = ROOT.TFile(resolution_v4fromv2_file_name, "read")
hResolution_v4fromv2 = resolution_v4fromv2_file.Get(res_histo_name)
hResolution_v4fromv2.SetDirectory(0)
utils.setHistStyle(hResolution_v4fromv2, ROOT.kBlue + 2)

res_dir = output_file.mkdir("resolutions")
res_dir.cd()

cent_limits_tens = [0, 10, 20, 30, 40, 50, 60, 70, 80]

# evaluating true resolutions
resolutions = []
hResolution.Write("hResolution_FT0C_FT0A_TPC_EP_standard")
for i_cent, cent in enumerate(centrality_classes):
    cent_left_index = cent_limits_tens.index(cent[0])
    cent_right_index = cent_limits_tens.index(cent[1])
    print(f"cent_limits_tens[{cent_left_index}]: {cent_limits_tens[cent_left_index]}")
    print(f"cent_limits_tens[{cent_right_index}]: {cent_limits_tens[cent_right_index]}")
    if cent_right_index - cent_left_index == 1:
        resol = hResolution.GetBinContent(cent_left_index + 1)
        print(f"cent: {cent[0]} - {cent[1]}, resolution: {resol}")
    else:
        resol = np.mean(
            [
                hResolution.GetBinContent(i)
                for i in range(cent_left_index + 1, cent_right_index + 1)
            ]
        )
        print(f"cent: {cent[0]} - {cent[1]}, resolution: {resol}")
    resolutions.append(resol)

# evaluating v4 from v2 resolutions
resolutions_v4fromv2 = []
hResolution_v4fromv2.Write("hResolution_FT0C_FT0A_TPC_EP_v4fromv2")
for i_cent, cent in enumerate(centrality_classes):
    cent_left_index = cent_limits_tens.index(cent[0])
    cent_right_index = cent_limits_tens.index(cent[1])
    print(f"cent_limits_tens[{cent_left_index}]: {cent_limits_tens[cent_left_index]}")
    print(f"cent_limits_tens[{cent_right_index}]: {cent_limits_tens[cent_right_index]}")
    if cent_right_index - cent_left_index == 1:
        resol = hResolution_v4fromv2.GetBinContent(cent_left_index + 1)
        print(f"cent: {cent[0]} - {cent[1]}, resolution: {resol}")
    else:
        resol = np.mean(
            [
                hResolution_v4fromv2.GetBinContent(i)
                for i in range(cent_left_index + 1, cent_right_index + 1)
            ]
        )
        print(f"cent: {cent[0]} - {cent[1]}, resolution: {resol}")
    resolutions_v4fromv2.append(resol)


# get histograms from QC file
qc_file = ROOT.TFile("../results_phi_studies/qc.root", "read")

# ge v2 histograms from final file
final_file_name = "../results_phi_studies/final.root"
final_file = ROOT.TFile(final_file_name, "read")

vV2stat = []
vV2syst = []
vV4 = []

for i_cent, cent in enumerate(centrality_classes):
    hV2stat = final_file.Get(
        f"cent_{cent[0]}_{cent[1]}/hV2vsPt_cent_{cent[0]}_{cent[1]}"
    )
    hV2stat.SetDirectory(0)
    vV2stat.append(hV2stat)
    hV2syst = final_file.Get(
        f"cent_{cent[0]}_{cent[1]}/hV2vsPt_cent_{cent[0]}_{cent[1]}_syst"
    )
    hV2syst.SetDirectory(0)
    vV2syst.append(hV2syst)

vV2fromFit = []
vV4fromFit = []
vSigmaFromFit = []

# Before the loop (so it persists)
dummy_hist = ROOT.TH1F("dummy_hist", "", 1, 0, 1)
dummy_hist.SetLineColor(ROOT.kBlue)
dummy_hist.SetLineWidth(2)

dummy_hist_v2only = ROOT.TH1F("dummy_hist_v2only", "", 1, 0, 1)
dummy_hist_v2only.SetLineColor(ROOT.kRed)
dummy_hist_v2only.SetLineWidth(2)

cPhiWithPred_list = []

for i_cent, cent in enumerate(centrality_classes):

    cPhiWithPred_list.append([])

    cent_dir = output_file.mkdir(f"cent_{cent[0]}_{cent[1]}")

    pt_bins_arr = np.array(pt_bins[i_cent], dtype=np.float64)
    n_pt_bins = len(pt_bins[i_cent]) - 1

    hV2fromFit = ROOT.TH1F(
        f"hV2FromFit_cent_{cent[0]}_{cent[1]}",
        r"; #it{p}_{T} (GeV/#it{c}); v_{2}",
        n_pt_bins,
        pt_bins_arr,
    )
    utils.setHistStyle(hV2fromFit, cent_colours[i_cent], marker=25)
    vV2fromFit.append(hV2fromFit)

    hV4 = ROOT.TH1F(
        f"hV4_cent_{cent[0]}_{cent[1]}",
        r"; #it{p}_{T} (GeV/#it{c}); v_{4}",
        n_pt_bins,
        pt_bins_arr,
    )
    utils.setHistStyle(hV4, cent_colours[i_cent])
    vV4.append(hV4)

    hV4fromFit = ROOT.TH1F(
        f"hV4FromFit_cent_{cent[0]}_{cent[1]}",
        r"; #it{p}_{T} (GeV/#it{c}); v_{4}",
        n_pt_bins,
        pt_bins_arr,
    )
    utils.setHistStyle(hV4fromFit, cent_colours[i_cent], marker=25)
    vV4fromFit.append(hV4fromFit)

    hSigmaFromFit = ROOT.TH1F(
        f"hSigmaFromFit_cent_{cent[0]}_{cent[1]}",
        r"; #it{p}_{T} (GeV/#it{c}); #sigma_{fit}",
        n_pt_bins,
        pt_bins_arr,
    )
    utils.setHistStyle(hSigmaFromFit, cent_colours[i_cent], marker=25)
    vSigmaFromFit.append(hSigmaFromFit)

    cent_plot_dir = f"../results_phi_studies/phi_studies/cent_{cent[0]}_{cent[1]}"
    os.makedirs(cent_plot_dir, exist_ok=True)

    for i_pt in range(0, n_pt_bins):

        bin_centre = (pt_bins[i_cent][i_pt] + pt_bins[i_cent][i_pt + 1]) / 2
        bin_left = pt_bins[i_cent][i_pt]
        bin_right = pt_bins[i_cent][i_pt + 1]

        hPhiMinusPsi_name = f"{cent[0]}_{cent[1]}/FT0C/hPhiMinusPsi_FT0C_cent_{cent[0]}_{cent[1]}_pt_{i_pt}"
        print(f"Getting histogram: {hPhiMinusPsi_name}")
        hPhiMinusPsi = qc_file.Get(hPhiMinusPsi_name)
        hPhiMinusPsi.SetDirectory(0)
        hPhiMinusPsi.Rebin(4)
        utils.setHistStyle(hPhiMinusPsi, ROOT.kBlack)
        cent_dir.cd()
        hPhiMinusPsi.Write()

        # RooFit convolution fit
        x = ROOT.RooRealVar("x", "x", 0, ROOT.TMath.Pi())
        # Definisci la normalizzazione come RooRealVar (yield)
        norm = ROOT.RooRealVar("norm", "norm", hPhiMinusPsi.Integral(), 0, 1e6)
        v2raw = ROOT.RooRealVar("v2raw", "v2raw", 0.1, -1, 1)
        v4raw = ROOT.RooRealVar("v4raw", "v4raw", 0.01, -1, 1)
        sigma = ROOT.RooRealVar(
            "sigma", "sigma", 0.01, 0.00001, 0.1
        )  # width of Gaussian

        base_formula = "1+2*v2raw*exp(-2*sigma*sigma)*cos(2*x)+2*v4raw*exp(-8*sigma*sigma)*cos(4*x)"
        base_func = ROOT.RooGenericPdf(
            "base_func", base_formula, ROOT.RooArgList(x, v2raw, v4raw, sigma)
        )

        # phi resolution vs pT
        pol2_func = ROOT.TF1("pol2_func", "pol2", 0, 10)
        pol2_func.SetParameter(0, 0.002)
        pol2_func.SetParameter(1, 0)
        pol2_func.SetParameter(2, 0.00047)

        phi_resolution = pol2_func.Eval(bin_centre / 2)
        sigma.setVal(phi_resolution)
        sigma.setConstant(True)

        # Crea il modello esteso: RooAddPdf con la funzione base come unico componente
        model = ROOT.RooAddPdf(
            "model", "model", ROOT.RooArgList(base_func), ROOT.RooArgList(norm)
        )

        dh = ROOT.RooDataHist(
            f"dhPhiMinusPsi_FT0C_cent_{cent[0]}_{cent[1]}_pt_{i_pt}",
            f"dhPhiMinusPsi_FT0C_cent_{cent[0]}_{cent[1]}_pt_{i_pt}",
            ROOT.RooArgList(x),
            hPhiMinusPsi,
        )
        model.fitTo(dh, ROOT.RooFit.Extended(True), ROOT.RooFit.PrintLevel(-1))

        # Extract fit results
        norm_val = norm.getVal()
        norm_err = norm.getError()
        v2raw_val = v2raw.getVal()
        v2raw_err = v2raw.getError()
        v4raw_val = v4raw.getVal()
        v4raw_err = v4raw.getError()
        chi2 = base_func.createChi2(dh).getVal()
        ndf = hPhiMinusPsi.GetNbinsX() - 3  # 3 parameters

        # Create a new norm for the reduced fit
        norm_reduced = ROOT.RooRealVar(
            "norm_reduced", "norm_reduced", hPhiMinusPsi.Integral(), 0, 1e6
        )
        v2raw_reduced = ROOT.RooRealVar("v2raw_reduced", "v2raw_reduced", 0.1, -1, 1)
        sigma_reduced = ROOT.RooRealVar(
            "sigma_reduced", "sigma_reduced", 0.01, 0.00001, 0.1
        )  # width of Gaussian
        sigma_reduced.setVal(phi_resolution)
        sigma_reduced.setConstant(True)
        base_formula_reduced = (
            "1+2*v2raw_reduced*exp(-2*sigma_reduced*sigma_reduced)*cos(2*x)"
        )
        base_func_reduced = ROOT.RooGenericPdf(
            "base_func_reduced",
            base_formula_reduced,
            ROOT.RooArgList(x, v2raw_reduced, sigma_reduced),
        )
        model_reduced = ROOT.RooAddPdf(
            "model_reduced",
            "model_reduced",
            ROOT.RooArgList(base_func_reduced),
            ROOT.RooArgList(norm_reduced),
        )

        model_reduced.fitTo(dh, ROOT.RooFit.Extended(True), ROOT.RooFit.PrintLevel(-1))
        chi2_reduced = base_func_reduced.createChi2(dh).getVal()
        ndf_reduced = hPhiMinusPsi.GetNbinsX() - 2  # 2 parameters

        # For plotting: RooFit
        frame = x.frame()
        dh.plotOn(frame)
        base_func.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kBlue))
        base_func_reduced.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kRed))

        frame.GetXaxis().SetTitle("#phi - #Psi_{2}")
        nbins = hPhiMinusPsi.GetNbinsX()
        frame.GetYaxis().SetTitle(f"Entries / (#pi/{nbins})")

        # Trova min e max dei punti (con errore) nel frame
        # Trova min e max dei punti (con errore) da hPhiMinusPsi
        ymin = float("inf")
        ymax = float("-inf")
        for i in range(1, hPhiMinusPsi.GetNbinsX() + 1):
            bin_content = hPhiMinusPsi.GetBinContent(i)
            bin_error = hPhiMinusPsi.GetBinError(i)
            if bin_content - bin_error < ymin:
                ymin = bin_content - bin_error
            if bin_content + bin_error > ymax:
                ymax = bin_content + bin_error

        # Imposta il range dell'asse y del frame
        frame.SetMinimum(0.90 * ymin)
        frame.SetMaximum(1.05 * ymax)

        # Titolo con centralità e intervallo di pt
        pt_low = pt_bins[i_cent][i_pt]
        pt_high = pt_bins[i_cent][i_pt + 1]
        canvas_title = f"Centrality: {cent[0]}-{cent[1]}%, {pt_low:.2f} #leq #it{{p}}_{{T}} < {pt_high:.2f} GeV/c"

        c = ROOT.TCanvas(
            f"cPhiMinusPsiFT0C_{cent[0]}_{cent[1]}_pt{i_pt}",
            canvas_title,
            800,
            600,
        )
        frame.SetTitle(canvas_title)
        frame.Draw()

        legFit = ROOT.TLegend(0.35, 0.51, 0.65, 0.85)
        legFit.SetFillColor(0)
        legFit.SetFillStyle(0)
        legFit.SetTextSize(0.04)
        legFit.SetTextFont(42)
        legFit.SetTextColor(1)
        legFit.SetBorderSize(0)
        # Add entries to the legend
        legFit.AddEntry(hPhiMinusPsi, r"{}^{3}#bar{He}", "PE")
        legFit.AddEntry(
            dummy_hist_v2only,
            "#it{A} [1 + 2#it{B} e^{-2#sigma^{2}}cos(2#it{#phi})]",
            "L",
        )
        legFit.AddEntry(0, f"#chi^{{2}}/ndf = {chi2_reduced:.1f}/{ndf_reduced}", "")
        legFit.AddEntry(
            dummy_hist,
            "#it{A} [1 + 2#it{B} e^{-2#sigma^{2}}cos(2#it{#phi}) + ",
            "L",
        )
        legFit.AddEntry(
            0,
            "+ 4#it{C} e^{-8#sigma^{2}}cos(4#it{#phi})]",
            "",
        )
        legFit.AddEntry(0, f"#chi^{{2}}/ndf = {chi2:.1f}/{ndf}", "")
        legFit.Draw()

        # Salva il canvas nel file ROOT
        cent_dir.cd()
        c.Write()

        # Salva il canvas come PDF nella cartella
        pdf_path = os.path.join(
            cent_plot_dir, f"cPhiMinusPsi_FT0C_cent_{cent[0]}_{cent[1]}_pt_{i_pt}.pdf"
        )
        c.SaveAs(pdf_path)

        vV2fromFit[i_cent].SetBinContent(i_pt + 1, v2raw_val)
        vV2fromFit[i_cent].SetBinError(i_pt + 1, v2raw_err)
        vV4fromFit[i_cent].SetBinContent(i_pt + 1, v4raw_val)
        vV4fromFit[i_cent].SetBinError(i_pt + 1, v4raw_err)
        vSigmaFromFit[i_cent].SetBinContent(i_pt + 1, sigma.getVal())
        vSigmaFromFit[i_cent].SetBinError(i_pt + 1, sigma.getError())

        # Draw theoretical model prediction as an orange band if available
        pred_file = ROOT.TFile(
            "../theoretical_models/Predictions_august2025.root", "READ"
        )
        pred_graph_name = (
            f"hPredPhi_pt_{bin_left:.1f}_{bin_right:.1f}_{cent[0]}{cent[1]}"
        )
        hPred = pred_file.Get(pred_graph_name)
        hPred.SetDirectory(0)
        hPredClone = hPred.Clone(
            f"hPredClone_{bin_left}_{bin_right}_{cent[0]}{cent[1]}"
        )

        # Fold hPredClone from [0,2pi] to [0,pi]
        n_bins = hPredClone.GetNbinsX()
        x_min = hPredClone.GetXaxis().GetXmin()
        x_max = hPredClone.GetXaxis().GetXmax()
        pi = np.pi
        # Only fold if the range is [0,2pi] or similar
        if x_max > pi + 0.1:
            # Create a new histogram with half the bins, range [0,pi]
            n_bins_folded = int(n_bins // 2)
            hPredFolded = ROOT.TH1D(
                f"hPredFolded_{bin_left}_{bin_right}_{cent[0]}{cent[1]}",
                hPredClone.GetTitle(),
                n_bins_folded,
                0,
                pi,
            )
            orange = ROOT.TColor.GetColor("#ff9900")
            hPredFolded.SetLineColor(orange)
            hPredFolded.SetLineWidth(2)
            hPredFolded.SetFillColor(orange)
            hPredFolded.SetFillStyle(3001)
            for i in range(1, n_bins + 1):
                x = hPredClone.GetBinCenter(i)
                y = hPredClone.GetBinContent(i)
                err = hPredClone.GetBinError(i)
                if x <= pi:
                    bin_fold = hPredFolded.FindBin(x)
                    hPredFolded.SetBinContent(
                        bin_fold, hPredFolded.GetBinContent(bin_fold) + y
                    )
                    hPredFolded.SetBinError(
                        bin_fold,
                        np.sqrt(err**2),
                    )
                else:
                    x_fold = 2 * pi - x
                    bin_fold = hPredFolded.FindBin(x_fold)
                    hPredFolded.SetBinContent(
                        bin_fold, hPredFolded.GetBinContent(bin_fold) + y
                    )
                    hPredFolded.SetBinError(
                        bin_fold,
                        np.sqrt(hPredFolded.GetBinError(bin_fold) ** 2 + err**2),
                    )

        hPredFolded.Scale(1.0 / hPredFolded.Integral("width"))
        hPhiMinusPsiClone = hPhiMinusPsi.Clone(
            f"hPhiMinusPsiClone_{bin_left}_{bin_right}_{cent[0]}{cent[1]}"
        )
        hPhiMinusPsiClone.Scale(1.0 / hPhiMinusPsiClone.Integral("width"))

        # --- New canvas: hPhiMinusPsi (PE) and hPredClone (E3) ---

        cPred = ROOT.TCanvas(
            f"cPhiWithPred_pt_{bin_left}_{bin_right}_{cent[0]}{cent[1]}",
            canvas_title,
            800,
            600,
        )

        # Draw prediction band first for proper layering
        hPredFolded.Draw("E3")
        hPhiMinusPsiClone.Draw("PE SAME")

        # Add TPaveText with ALICE, centrality, and pt info
        pave_pred = ROOT.TPaveText(0.39, 0.73, 0.68, 0.87, "NDC")
        pave_pred.SetFillColor(0)
        pave_pred.SetBorderSize(0)
        pave_pred.SetTextFont(42)
        pave_pred.SetTextSize(0.035)
        pave_pred.AddText("ALICE")
        pave_pred.AddText(f"Pb-Pb {cent[0]}-{cent[1]}%")
        pave_pred.AddText(f"{pt_low:.2f} < #it{{p}}_{{T}} < {pt_high:.2f} GeV/#it{{c}}")
        pave_pred.Draw()

        # Legend for this canvas
        legPred = ROOT.TLegend(0.40, 0.59, 0.63, 0.73)
        legPred.SetFillColor(0)
        legPred.SetFillStyle(0)
        legPred.SetTextSize(0.04)
        legPred.SetTextFont(42)
        legPred.SetTextColor(1)
        legPred.SetBorderSize(0)
        legPred.AddEntry(hPhiMinusPsi, r"{}^{3}#bar{He}", "PE")
        legPred.AddEntry(hPredFolded, "Theory prediction", "F")
        legPred.Draw()

        # Save the canvas in the output file and as PDF
        cent_dir.cd()
        hPredFolded.Write()

        cPred.Write()
        pdf_path_pred = os.path.join(
            cent_plot_dir,
            f"cPhiWithPred_pt_{bin_left}_{bin_right}_{cent[0]}{cent[1]}.pdf",
        )
        cPred.SaveAs(pdf_path_pred)

        if i_cent == 0 and i_pt == 0:
            cPred.SaveAs(f"{output_dir}/cPhiWithPred_all.pdf[")
        cPred.SaveAs(f"{output_dir}/cPhiWithPred_all.pdf")
        if i_cent == len(centrality_classes) - 1 and i_pt == n_pt_bins - 1:
            cPred.SaveAs(f"{output_dir}/cPhiWithPred_all.pdf]")

        pred_file.Close()

        # compute v4 from v2
        hCos4PhiMinusPsi_name = f"{cent[0]}_{cent[1]}/FT0C/hCos4PhiMinusPsi_FT0C_cent_{cent[0]}_{cent[1]}_pt_{i_pt}"
        print(f"Getting histogram: {hCos4PhiMinusPsi_name}")
        hCos4PhiMinusPsi = qc_file.Get(hCos4PhiMinusPsi_name)
        hCos4PhiMinusPsi.SetDirectory(0)
        hCos4PhiMinusPsi.Rebin(4)
        utils.setHistStyle(hCos4PhiMinusPsi, ROOT.kBlack)
        cent_dir.cd()
        hCos4PhiMinusPsi.Write()

        # Crea e salva il canvas per la distribuzione di cos(4(phi-psi))
        cCos4 = ROOT.TCanvas(
            f"cCos4PhiMinusPsiFT0C_{cent[0]}_{cent[1]}_pt{i_pt}",
            f"cCos4PhiMinusPsiFT0C_{cent[0]}_{cent[1]}_pt{i_pt}",
            800,
            600,
        )
        hCos4PhiMinusPsi.SetMaximum(1.2 * hCos4PhiMinusPsi.GetMaximum())
        hCos4PhiMinusPsi.Draw("PE")

        # TPaveText con centralità e range di pt
        pave_cos4 = ROOT.TPaveText(0.55, 0.75, 0.89, 0.89, "NDC")
        pave_cos4.SetFillColor(0)
        pave_cos4.SetBorderSize(0)
        pave_cos4.SetTextFont(42)
        pave_cos4.SetTextSize(0.035)
        pave_cos4.AddText(f"Centrality: {cent[0]}-{cent[1]}%")
        pt_low = pt_bins[i_cent][i_pt]
        pt_high = pt_bins[i_cent][i_pt + 1]
        pave_cos4.AddText(f"{pt_low:.1f} < #it{{p}}_{{T}} < {pt_high:.1f} GeV/c")
        pave_cos4.Draw()

        # Salva il canvas nel file ROOT
        cCos4.Write()

        # Salva il canvas come PDF nella cartella
        pdf_path_cos4 = os.path.join(
            cent_plot_dir,
            f"cCos4PhiMinusPsi_FT0C_cent_{cent[0]}_{cent[1]}_pt_{i_pt}.pdf",
        )
        cCos4.SaveAs(pdf_path_cos4)

        v4_val = hCos4PhiMinusPsi.GetMean()
        v4_err = hCos4PhiMinusPsi.GetMeanError()

        vV4[i_cent].SetBinContent(i_pt + 1, v4_val)
        vV4[i_cent].SetBinError(i_pt + 1, v4_err)

    vV2fromFit[i_cent].Scale(1.0 / resolutions[i_cent])
    vV4[i_cent].Scale(1 / resolutions_v4fromv2[i_cent])
    vV4fromFit[i_cent].Scale(1.0 / resolutions_v4fromv2[i_cent])

    cent_dir.cd()
    vV2fromFit[i_cent].Write()
    vV2stat[i_cent].Write()
    vV2syst[i_cent].Write()
    vV4[i_cent].Write()
    vV4fromFit[i_cent].Write()

    cV2fromFit = ROOT.TCanvas(
        f"cV2FromFit_cent_{cent[0]}_{cent[1]}",
        f"cV2FromFit_cent_{cent[0]}_{cent[1]}",
        800,
        600,
    )

    # Esempio per vV2stat, vV2syst, vV2fromFit
    hists = [vV2stat[i_cent], vV2syst[i_cent], vV2fromFit[i_cent]]
    ymin, ymax = utils.get_min_max_with_errors(hists)
    # Margine extra
    yrange = ymax - ymin
    ymin -= 0.1 * yrange
    ymax += 0.1 * yrange

    frame = cV2fromFit.DrawFrame(
        hists[0].GetXaxis().GetXmin(),
        ymin,
        hists[0].GetXaxis().GetXmax(),
        ymax,
        f";{hists[0].GetXaxis().GetTitle()};{hists[0].GetYaxis().GetTitle()}",
    )
    for hist in hists:
        hist.Draw("PEsame")

    legend = ROOT.TLegend(0.20, 0.70, 0.39, 0.78, "", "brNDC")
    legend.SetBorderSize(0)
    legend.AddEntry(vV2stat[i_cent], "v_{2} standard", "PF")
    legend.AddEntry(vV2fromFit[i_cent], "v_{2} from fit", "PE")
    legend.Draw()

    # Add centrality class panel to cV2fromFit
    cent_pave = ROOT.TPaveText(0.18, 0.80, 0.38, 0.88, "NDC")
    cent_pave.SetFillColor(0)
    cent_pave.SetBorderSize(0)
    cent_pave.SetTextFont(42)
    cent_pave.SetTextSize(0.04)
    cent_pave.AddText(f"Centrality: {cent[0]}-{cent[1]}%")
    cent_pave.Draw()

    cV2fromFit.Write()
    if i_cent == 0:
        cV2fromFit.SaveAs(f"{output_dir}/cV2compVsPt.pdf[")
    cV2fromFit.SaveAs(f"{output_dir}/cV2compVsPt.pdf")

    cV4fromFit = ROOT.TCanvas(
        f"cV4fromFit_cent_{cent[0]}_{cent[1]}",
        f"cV4fromFit_cent_{cent[0]}_{cent[1]}",
        800,
        600,
    )
    hists = [vV4[i_cent], vV4fromFit[i_cent]]
    ymin, ymax = utils.get_min_max_with_errors(hists)
    # Margine extra
    yrange = ymax - ymin
    ymin -= 0.1 * yrange
    ymax += 0.1 * yrange
    frame = cV4fromFit.DrawFrame(
        hists[0].GetXaxis().GetXmin(),
        ymin,
        hists[0].GetXaxis().GetXmax(),
        ymax,
        f";{hists[0].GetXaxis().GetTitle()};{hists[0].GetYaxis().GetTitle()}",
    )
    vV4[i_cent].Draw("PEsame")
    vV4fromFit[i_cent].Draw("PEsame")

    legend = ROOT.TLegend(0.20, 0.70, 0.39, 0.78, "", "brNDC")
    legend.SetBorderSize(0)
    legend.AddEntry(vV4[i_cent], "v_{4} standard wrt #Psi_{2}", "PF")
    legend.AddEntry(vV4fromFit[i_cent], "v_{4} from fit", "PE")
    legend.Draw()

    # Add centrality class panel to cV4fromFit
    cent_pave = ROOT.TPaveText(0.18, 0.80, 0.38, 0.88, "NDC")
    cent_pave.SetFillColor(0)
    cent_pave.SetBorderSize(0)
    cent_pave.SetTextFont(42)
    cent_pave.SetTextSize(0.04)
    cent_pave.AddText(f"Centrality: {cent[0]}-{cent[1]}%")
    cent_pave.Draw()

    cV4fromFit.Write()
    if i_cent == 0:
        cV4fromFit.SaveAs(f"{output_dir}/cV4compVsPt.pdf[")
    cV4fromFit.SaveAs(f"{output_dir}/cV4compVsPt.pdf")

    cSigmaFromFit = ROOT.TCanvas(
        f"cSigmaFromFit_cent_{cent[0]}_{cent[1]}",
        f"cSigmaFromFit_cent_{cent[0]}_{cent[1]}",
        800,
        600,
    )
    hists = [vSigmaFromFit[i_cent]]
    ymin, ymax = utils.get_min_max_with_errors(hists)
    # Margine extra
    yrange = ymax - ymin
    ymin -= 0.1 * yrange
    ymax += 0.1 * yrange
    frame = cV2fromFit.DrawFrame(
        hists[0].GetXaxis().GetXmin(),
        ymin,
        hists[0].GetXaxis().GetXmax(),
        ymax,
        f";{hists[0].GetXaxis().GetTitle()};{hists[0].GetYaxis().GetTitle()}",
    )
    vSigmaFromFit[i_cent].Draw("PE SAME")

    legend = ROOT.TLegend(0.20, 0.70, 0.39, 0.78, "", "brNDC")
    legend.SetBorderSize(0)
    legend.AddEntry(vSigmaFromFit[i_cent], "#sigma_{fit} from convolution", "PE")
    legend.Draw()

    cent_pave = ROOT.TPaveText(0.18, 0.80, 0.38, 0.88, "NDC")
    cent_pave.SetFillColor(0)
    cent_pave.SetBorderSize(0)
    cent_pave.SetTextFont(42)
    cent_pave.SetTextSize(0.04)
    cent_pave.AddText(f"Centrality: {cent[0]}-{cent[1]}%")
    cent_pave.Draw()

    cSigmaFromFit.Write()
    pdf_path_sigma = os.path.join(output_dir, f"cSigmaFromFit.pdf")
    if i_cent == 0:
        cSigmaFromFit.SaveAs(f"{output_dir}/cSigmaFromFit.pdf[")
    cSigmaFromFit.SaveAs(f"{output_dir}/cSigmaFromFit.pdf")


cV2fromFit.SaveAs(f"{output_dir}/cV2compVsPt.pdf]")

cV2allCent = ROOT.TCanvas("cV2allCent", "v2 vs pT for all centrality classes", 800, 600)
cV2allCent.cd()


# Draw all centrality classes on the same canvas
drawn = False
legend_centrality = ROOT.TLegend(
    0.15, 0.62, 0.45, 0.84, "^{3}He, FT0C centrality", "brNDC"
)
legend_centrality.SetBorderSize(0)
legend_centrality.SetTextFont(42)
legend_centrality.SetTextSize(0.035)
legend_centrality.SetFillStyle(0)

for i_cent, cent in enumerate(centrality_classes):
    # Draw hV2stat
    opt = "PE" if not drawn else "PE SAME"
    vV2stat[i_cent].Draw(opt)
    legend_centrality.AddEntry(vV2stat[i_cent], f"{cent[0]}-{cent[1]}%", "FP")
    # Draw hV2syst
    vV2syst[i_cent].Draw("PE2x0 SAME")
    # Draw vV2fromFit
    vV2fromFit[i_cent].SetMarkerStyle(25)
    vV2fromFit[i_cent].Draw("PE SAME")
    drawn = True

# Set axis titles and ranges
first_hist = vV2stat[0]
first_hist.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")
first_hist.GetYaxis().SetTitle("v_{2}")
first_hist.GetYaxis().SetRangeUser(0, 0.8)

# Add a second legend for marker meaning
legend_markers = ROOT.TLegend(0.17, 0.51, 0.37, 0.62, "", "brNDC")
legend_markers.SetBorderSize(0)
legend_markers.SetTextFont(42)
legend_markers.SetTextSize(0.035)
legend_markers.SetFillStyle(0)

# Dummy for v2, SP (same marker as hV2stat)
dummy_v2sp = ROOT.TH1F("dummy_v2sp", "", 1, 0, 1)
dummy_v2sp.SetMarkerStyle(vV2stat[0].GetMarkerStyle())
dummy_v2sp.SetMarkerColor(ROOT.kGray + 2)
dummy_v2sp.SetLineColor(ROOT.kGray + 2)
dummy_v2sp.SetMarkerSize(vV2stat[0].GetMarkerSize())
legend_markers.AddEntry(dummy_v2sp, "v_{2}, Scalar Product", "P")

# Dummy for parameter B from the fit (same marker as vV2fromFit)
dummy_v2fit = ROOT.TH1F("dummy_v2fit", "", 1, 0, 1)
dummy_v2fit.SetMarkerStyle(25)
dummy_v2fit.SetMarkerColor(ROOT.kGray + 2)
dummy_v2fit.SetLineColor(ROOT.kGray + 2)
dummy_v2fit.SetMarkerSize(vV2fromFit[0].GetMarkerSize())
legend_markers.AddEntry(dummy_v2fit, "parameter B from the fit", "P")

# Info panel (top-right corner)
info_panel_comp = ROOT.TPaveText(0.65, 0.75, 0.85, 0.85, "NDC")
info_panel_comp.SetFillColor(0)
info_panel_comp.SetBorderSize(0)
info_panel_comp.SetTextFont(42)
info_panel_comp.SetTextSize(0.04)
info_panel_comp.AddText(r"ALICE")
info_panel_comp.AddText(r"Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.36 TeV")
info_panel_comp.Draw()

legend_centrality.Draw()
legend_markers.Draw()
output_file.cd()
cV2allCent.Write()
cV2allCent.SaveAs(f"{output_dir}/cV2allCent.pdf")

# New canvas for all vV4fromFit
cV4allCent = ROOT.TCanvas(
    "cV4allCent", "B parameter from fit for all centrality classes", 800, 600
)
cV4allCent.cd()

drawn = False
legend_centrality_v4 = ROOT.TLegend(
    0.15, 0.62, 0.45, 0.84, "^{3}He, FT0C centrality", "brNDC"
)
legend_centrality_v4.SetBorderSize(0)
legend_centrality_v4.SetTextFont(42)
legend_centrality_v4.SetTextSize(0.035)
legend_centrality_v4.SetFillStyle(0)

for i_cent, cent in enumerate(centrality_classes):
    opt = "PE" if not drawn else "PE SAME"
    vV4fromFit[i_cent].Draw(opt)
    legend_centrality_v4.AddEntry(vV4fromFit[i_cent], f"{cent[0]}-{cent[1]}%", "FP")
    drawn = True

# Set axis titles and ranges
first_hist_v4 = vV4fromFit[0]
first_hist_v4.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")
first_hist_v4.GetYaxis().SetTitle("C parameter from fit")
first_hist_v4.GetYaxis().SetRangeUser(-0.03, 0.40)

# Add info panel (reuse same as before)
info_panel_v4 = ROOT.TPaveText(0.65, 0.75, 0.85, 0.85, "NDC")
info_panel_v4.SetFillColor(0)
info_panel_v4.SetBorderSize(0)
info_panel_v4.SetTextFont(42)
info_panel_v4.SetTextSize(0.04)
info_panel_v4.AddText(r"ALICE")
info_panel_v4.AddText(r"Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.36 TeV")
info_panel_v4.Draw()

# Add dashed black horizontal line at y=0
line0 = ROOT.TLine(
    first_hist_v4.GetXaxis().GetXmin(), 0, first_hist_v4.GetXaxis().GetXmax(), 0
)
line0.SetLineColor(ROOT.kBlack)
line0.SetLineStyle(2)
line0.SetLineWidth(2)
line0.Draw()

legend_centrality_v4.Draw()
output_file.cd()
cV4allCent.Write()
cV4allCent.SaveAs(f"{output_dir}/cV4allCent.pdf")
cV4fromFit.SaveAs(f"{output_dir}/cV4compVsPt.pdf]")
cSigmaFromFit.SaveAs(f"{output_dir}/cSigmaFromFit.pdf]")
