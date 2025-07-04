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
output_file = ROOT.TFile(f"{output_dir}/phi_studies.root", "recreate")
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

# Crea un file ROOT per i canvas
output_root = ROOT.TFile("phi_studies.root", "RECREATE")

vV2fromFit = []
vV4fromFit = []

for i_cent, cent in enumerate(centrality_classes):

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

    cent_plot_dir = f"../results_phi_studies/phi_studies/cent_{cent[0]}_{cent[1]}"
    os.makedirs(cent_plot_dir, exist_ok=True)

    for i_pt in range(0, n_pt_bins):

        hPhiMinusPsi_name = f"{cent[0]}_{cent[1]}/FT0C/hPhiMinusPsi_FT0C_cent_{cent[0]}_{cent[1]}_pt_{i_pt}"
        print(f"Getting histogram: {hPhiMinusPsi_name}")
        hPhiMinusPsi = qc_file.Get(hPhiMinusPsi_name)
        hPhiMinusPsi.SetDirectory(0)
        hPhiMinusPsi.Rebin(4)
        utils.setHistStyle(hPhiMinusPsi, ROOT.kBlack)
        cent_dir.cd()
        hPhiMinusPsi.Write()

        # Definisci la funzione di fit con nomi dei parametri
        fit_func = ROOT.TF1(
            "fit_func",
            "[0]*(1+2.0*[1]*TMath::Cos(2.0*x)+2.0*[2]*TMath::Cos(4.0*x))",
            0,
            ROOT.TMath.Pi(),
        )
        fit_func.SetParName(0, "Norm")
        fit_func.SetParName(1, "v2raw")
        fit_func.SetParName(2, "v4raw")

        # Imposta valori iniziali ragionevoli
        fit_func.SetParameters(hPhiMinusPsi.GetMaximum(), 0.1, 0.01)

        # Esegui il fit
        hPhiMinusPsi.Fit(fit_func, "R")

        # Parametri e errori
        norm = fit_func.GetParameter(0)
        norm_err = fit_func.GetParError(0)
        v2raw = fit_func.GetParameter(1)
        v2raw_err = fit_func.GetParError(1)
        v4raw = fit_func.GetParameter(2)
        v4raw_err = fit_func.GetParError(2)
        chi2 = fit_func.GetChisquare()
        ndf = fit_func.GetNDF()

        c = ROOT.TCanvas(
            f"cPhiMinusPsiFT0C_{cent[0]}_{cent[1]}_pt{i_pt}",
            f"cPhiMinusPsiFT0C_{cent[0]}_{cent[1]}_pt{i_pt}",
            800,
            600,
        )
        hPhiMinusPsi.SetMaximum(1.2 * hPhiMinusPsi.GetMaximum())
        hPhiMinusPsi.Draw("PE")
        fit_func.Draw("same")

        # TPaveText per i risultati del fit
        pave = ROOT.TPaveText(0.37, 0.63, 0.70, 0.87, "NDC")
        pave.SetFillColor(0)
        pave.SetBorderSize(0)
        pave.AddText(r"Fit: Norm * (1 + 2 v_{2}^{raw}cos(2x) + 2 v_{4}^{raw} cos(4x))")
        pave.AddText(f"Norm = {norm:.3g} #pm {norm_err:.2g}")
        pave.AddText(f"v2raw = {v2raw:.3g} #pm {v2raw_err:.2g}")
        pave.AddText(f"v4raw = {v4raw:.3g} #pm {v4raw_err:.2g}")
        pave.AddText(f"#chi^2/ndf = {chi2:.1f}/{ndf}")
        pave.AddText(f"Centrality: {cent[0]}-{cent[1]}%")  # <--- aggiungi questa riga
        pt_low = pt_bins[i_cent][i_pt]
        pt_high = pt_bins[i_cent][i_pt + 1]
        pave.AddText(f"{pt_low:.1f} < #it{{p}}_{{T}} < {pt_high:.1f} GeV/c")
        pave.Draw()

        # Salva il canvas nel file ROOT
        c.Write()

        # Salva il canvas come PDF nella cartella
        pdf_path = os.path.join(
            cent_plot_dir, f"cPhiMinusPsi_FT0C_cent_{cent[0]}_{cent[1]}_pt_{i_pt}.pdf"
        )
        c.SaveAs(pdf_path)

        vV2fromFit[i_cent].SetBinContent(i_pt + 1, v2raw)
        vV2fromFit[i_cent].SetBinError(i_pt + 1, v2raw_err)

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

        vV4fromFit[i_cent].SetBinContent(i_pt + 1, v4raw)
        vV4fromFit[i_cent].SetBinError(i_pt + 1, v4raw_err)

    vV2fromFit[i_cent].Scale(1.0 / resolutions[i_cent])
    vV4[i_cent].Scale(1.0 / resolutions_v4fromv2[i_cent])
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
    vV2stat[i_cent].Draw("PE")
    vV2syst[i_cent].Draw("PEsame")
    vV2fromFit[i_cent].Draw("PEsame")

    legend = ROOT.TLegend(0.69, 0.22, 0.88, 0.30, "", "brNDC")
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
    vV4[i_cent].Draw("PE")
    vV4fromFit[i_cent].Draw("PEsame")

    legend = ROOT.TLegend(0.69, 0.22, 0.88, 0.30, "", "brNDC")
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


cV2fromFit.SaveAs(f"{output_dir}/cV2compVsPt.pdf]")
cV4fromFit.SaveAs(f"{output_dir}/cV4compVsPt.pdf]")
