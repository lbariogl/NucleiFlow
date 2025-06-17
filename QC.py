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
input_file_AR_name = config["input_file_AR_name"]
resolution_file_name = config["resolution_file_name"]
output_dir_name = config["output_dir_name"]
output_file_name = config["output_file_name_qc"]

nuclei_tree_name = config["nuclei_tree_name"]
ep_tree_name = config["ep_tree_name"]

useSP = True  # config["useSP"]

mandatory_selections = config["mandatory_selections"]
selection_dict = config["selection_dict"]
selection_list = selection_dict.values()
selections = " and ".join(selection_list)

pattern = r"\s*and\s+(abs\(fNsigmaTPC3He\)\s*<\s*(\d+\.?\d*))"

# Search and extract the match without "and"
match = re.search(pattern, selections)
nSigmaTPC_selection_string = match.group(1).strip() if match else ""
n_sigma_selection_from_string = float(match.group(2)) if match else None

selections = re.sub(pattern, "", selections).strip()

# ptdep_selection_dict = config["ptdep_selection_dict"]["fAvgItsClusSizeCosLambda"]

cent_detector_label = config["cent_detector_label"]
reference_flow_detector = config["reference_flow_detector"]
resolution_flow_detectors = config["resolution_flow_detectors"]

centrality_classes = config["centrality_classes"]
pt_bins = config["pt_bins"]
pt_bins_no_cent = config["qc_pt_bins"]
n_pt_bins_no_cent = len(pt_bins_no_cent) - 1
cent_colours = config["cent_colours"]

do_syst = config["do_syst"]

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

# set input files for EP qc
input_file_AR = ROOT.TFile(input_file_AR_name)
input_dir_AR_general = input_file_AR.Get("flow-qc/general")
input_dir_AR_flow_ep = input_file_AR.Get("flow-qc/flow_ep")

# define new columns
utils.redefineColumns(complete_df, useSP=useSP)

# apply common selections
complete_df.query(f"{mandatory_selections} and {selections}", inplace=True)

# fill hNsigmaITSvsP before ITS-cluster size selection
hNsigmaITSvsP = ROOT.TH2F(
    "hNsigmaITSvsP",
    r"; #it{p} (GeV/#it{c}); n#sigma_{ITS}",
    120,
    0.0,
    12.0,
    500,
    -5.0,
    5.0,
)

for p, n_sigma_its in zip(complete_df["fP"], complete_df["fNsigmaITS3He"]):
    hNsigmaITSvsP.Fill(p, n_sigma_its)

its_nsigma_parameters = utils.its_nsigma_parameters
its_nsigma_func_string = utils.its_nsigma_func_string
its_nsigma_func = ROOT.TF1("its_nsigma_func", its_nsigma_func_string, 1.0, 15.0, 3)
its_nsigma_func.SetParameters(*its_nsigma_parameters)
its_nsigma_func.SetLineColor(ROOT.kRed)
its_nsigma_func.SetLineWidth(2)

cNsigmaITSvsP = ROOT.TCanvas("cNsigmaITSvsP", "cNsigmaITSvsP", 800, 600)
hNsigmaITSvsP.Draw("colz")
its_nsigma_func.Draw("SAME")
cNsigmaITSvsP.SaveAs(f"{output_dir_name}/qc_plots/cNsigmaITSvsP.pdf")

output_file.cd()
hNsigmaITSvsP.Write()
its_nsigma_func.Write()
cNsigmaITSvsP.Write()

complete_df.query("fNsigmaITS3HeMinusOffset > 0", inplace=True)


# fill hNsigmaITSvsP before ITS-cluster size selection after selection
hNsigmaITSvsPafterSel = ROOT.TH2F(
    "hNsigmaITSvsPafterSel",
    r"; #it{p} (GeV/#it{c}); n#sigma_{ITS}",
    120,
    0.0,
    12.0,
    500,
    -5.0,
    5.0,
)

for p, n_sigma_its in zip(complete_df["fP"], complete_df["fNsigmaITS3He"]):
    hNsigmaITSvsPafterSel.Fill(p, n_sigma_its)

cNsigmaITSvsPafterSel = ROOT.TCanvas(
    "cNsigmaITSvsPafterSel", "cNsigmaITSvsPafterSel", 800, 600
)
hNsigmaITSvsPafterSel.Draw("colz")
cNsigmaITSvsPafterSel.SaveAs(f"{output_dir_name}/qc_plots/cNsigmaITSvsPafterSel.pdf")

output_file.cd()
hNsigmaITSvsPafterSel.Write()
cNsigmaITSvsPafterSel.Write()

hEta = ROOT.TH1F("hEta", ";#eta;", 200, -1.0, 1.0)
utils.setHistStyle(hEta, ROOT.kRed + 2)

hTPCsignalVsPoverZ = ROOT.TH2F(
    "hTPCsignalVsPoverZ",
    r";#it{p}/z (GeV/#it{c}); d#it{E} / d#it{x} (a.u.)",
    300,
    0.0,
    6.0,
    1400,
    0.0,
    1400.0,
)

mass_bin_shift = int(utils.mass_alpha * utils.mass_alpha * 10) / 10.0

hTOFmassSquaredVsPt = ROOT.TH2F(
    "hTOFmassSquaredVsPt",
    r";#it{p}_{T} (GeV/#it{c}); m^{2}_{TOF} - m^{2}_{{}^{4}He}(Gev^{2}/c^{4})",
    160,
    0.0,
    8.0,
    210,
    4.0 - mass_bin_shift,
    25.0 - mass_bin_shift,
)

hNsigmaVsPt = ROOT.TH2F(
    "hNsigmaVsPt",
    r";#it{p}_{T} (GeV/#it{c}); n_{\sigma_{TPC}}",
    160,
    0.0,
    8.0,
    100,
    -5.0,
    5.0,
)

hNsigmaVsTOFmassSquared = ROOT.TH2F(
    "hNsigmaVsTOFmassSquared",
    r";n_{\sigma_{TPC}}; m^{2}_{TOF} - m^{2}_{{}^{4}He}(Gev^{2}/c^{4})",
    100,
    -5.0,
    5.0,
    210,
    4.0 - mass_bin_shift,
    25.0 - mass_bin_shift,
)

# get histograms for EP qc
hZvtx = input_dir_AR_general.Get("hRecVtxZData")
utils.setHistStyle(hZvtx, ROOT.kRed + 2)
hCentFT0C = input_dir_AR_general.Get("hCentFT0C")
utils.setHistStyle(hCentFT0C, ROOT.kRed + 2)

input_file_AR_small = ROOT.TFile(
    "/data/lbariogl/flow/LHC23zzh_pass4_small_new/AnalysisResults_general.root"
)

input_file_AR_small_igor = ROOT.TFile(
    "/data/lbariogl/flow/LHC23zzh_pass4_small_new/AnalysisResults_igor.root"
)
hEventSelections_normal = input_file_AR_small.Get(
    "nuclei-spectra/spectra/hEventSelections"
)
utils.setHistStyle(hEventSelections_normal, ROOT.kRed + 2)

hEventSelections_igor = input_file_AR_small_igor.Get(
    "nuclei-spectra/spectra/hEventSelections"
)
utils.setHistStyle(hEventSelections_igor, ROOT.kAzure + 1)

detectors = ["FT0A", "FT0C", "TPCl", "TPCr", "TPC"]
combinations_detectors = list(combinations(detectors, 2))
vDeltaPsi = {}

for det1, det2 in combinations_detectors:
    # check if the histograms are in the file
    if not input_dir_AR_flow_ep.Get(f"hDeltaPsi_{det1}_{det2}_EP"):
        print(f"Missing histogram hDeltaPsi_{det1}_{det2}_EP")
        continue
    vDeltaPsi[det1 + " - " + det2] = input_dir_AR_flow_ep.Get(
        f"hDeltaPsi_{det1}_{det2}_EP"
    )

# get histrogram for event-plane angle distributions
vPsi2D = {}
vPsi = {}
for det in detectors:
    hPsi2D = input_dir_AR_flow_ep.Get(f"hPsi_{det}_EP")
    vPsi2D[det] = hPsi2D
    vPsi[det] = []
    yMaxPsi = 0
    legend_psi = ROOT.TLegend(0.40, 0.22, 0.70, 0.45, det, "brNDC")
    legend_psi.SetNColumns(2)
    det_dir = output_file.mkdir(det)
    cPsi = ROOT.TCanvas(f"cPsi{det}", f"cPsi{det}", 800, 600)
    for i_cent, cent in enumerate(centrality_classes):
        left_bin = hPsi2D.GetXaxis().FindBin(cent[0])
        right_bin = hPsi2D.GetXaxis().FindBin(cent[1]) - 1
        hPsi = hPsi2D.ProjectionY(f"h{det}_{cent[0]}_{cent[1]}", left_bin, right_bin)
        utils.setHistStyle(hPsi, cent_colours[i_cent])
        hPsi.Scale(1.0 / hPsi.Integral())
        vPsi[det].append(hPsi)
        loc_max = hPsi.GetMaximum()
        if loc_max > yMaxPsi:
            yMaxPsi = loc_max
        legend_psi.AddEntry(hPsi, f"{cent[0]}-{cent[1]}%", "l")
        det_dir.cd()
        hPsi.Write()
    formatted_title = f";#Psi_{{2}}^{det}; norm. counts"
    cPsi.DrawFrame(-3.5, 0, 3.5, yMaxPsi * 1.1, formatted_title)
    for i_cent, cent in enumerate(centrality_classes):
        vPsi[det][i_cent].Draw("HISTO SAME")

    legend_psi.Draw("SAME")
    cPsi.SaveAs(f"{output_dir_name}/qc_plots/cPsi{det}.pdf")
    det_dir.cd()
    cPsi.Write()

vPhi = []
vPsiTree = {}
vPhiMinusPsi = {}
vSinPhiMinusPsi = {}
vCos2PhiMinusPsi = {}
vSin2PhiMinusPsi = {}
vSinPhiMinusPsiVsPt = {}
vCos2PhiMinusPsiVsPt = {}
vSin2PhiMinusPsiVsPt = {}
cSinPhiMinusPsiVsPt = {}
cCos2PhiMinusPsiVsPt = {}
cSin2PhiMinusPsiVsPt = {}
vQ = {}
vQx = {}
vQy = {}
vV2ScalarProduct = {}
vV2ScalarProduct_check = {}
vV2Tree = {}
for det in detectors:
    vPsiTree[det] = []
    vPhiMinusPsi[det] = []
    vSinPhiMinusPsi[det] = []
    vCos2PhiMinusPsi[det] = []
    vSin2PhiMinusPsi[det] = []
    vQ[det] = []
    vQx[det] = []
    vQy[det] = []
    vV2ScalarProduct[det] = []
    vV2ScalarProduct_check[det] = []
    vV2Tree[det] = []

    vSinPhiMinusPsiVsPt[det] = []
    vCos2PhiMinusPsiVsPt[det] = []
    vSin2PhiMinusPsiVsPt[det] = []

for i_cent, cent in enumerate(centrality_classes):
    cent_seleration = f"fCentFT0C > {cent[0]} and fCentFT0C < {cent[1]}"
    cent_dir = output_file.mkdir(f"{cent[0]}_{cent[1]}")
    cent_df = complete_df.query(cent_seleration, inplace=False)
    vPhi_cent = []
    phi_dir = cent_dir.mkdir("phi")

    pt_bins_arr = np.array(pt_bins[i_cent], dtype=np.float64)
    n_pt_bins = len(pt_bins[i_cent]) - 1

    for det in detectors:
        vPsiTree[det].append([])
        vPhiMinusPsi[det].append([])
        vSinPhiMinusPsi[det].append([])
        vCos2PhiMinusPsi[det].append([])
        vSin2PhiMinusPsi[det].append([])
        vQ[det].append([])
        vQx[det].append([])
        vQy[det].append([])
        vV2ScalarProduct[det].append([])
        vV2ScalarProduct_check[det].append([])
        vV2Tree[det].append([])
        cent_dir.mkdir(det)

        hSinPhiMinusPsiVsPt = ROOT.TH1F(
            f"hSinPhiMinusPsiVsPt_{det}_cent_{cent[0]}_{cent[1]}",
            rf"; #it{{p}}_{{T}} (GeV/#it{{c}}); sin(#phi - #Psi_{{{det}}})",
            n_pt_bins,
            pt_bins_arr,
        )
        utils.setHistStyle(hSinPhiMinusPsiVsPt, cent_colours[i_cent])
        vSinPhiMinusPsiVsPt[det].append(hSinPhiMinusPsiVsPt)

        hCos2PhiMinusPsiVsPt = ROOT.TH1F(
            f"hCos2PhiMinusPsiVsPt_{det}_cent_{cent[0]}_{cent[1]}",
            rf"; #it{{p}}_{{T}} (GeV/#it{{c}}); cos(2(#phi - #Psi_{{{det}}}))",
            n_pt_bins,
            pt_bins_arr,
        )
        utils.setHistStyle(hCos2PhiMinusPsiVsPt, cent_colours[i_cent])
        vCos2PhiMinusPsiVsPt[det].append(hCos2PhiMinusPsiVsPt)

        hSin2PhiMinusPsiVsPt = ROOT.TH1F(
            f"hSin2PhiMinusPsiVsPt_{det}_cent_{cent[0]}_{cent[1]}",
            rf"; #it{{p}}_{{T}} (GeV/#it{{c}}); sin(2(#phi - #Psi_{{{det}}}))",
            n_pt_bins,
            pt_bins_arr,
        )
        utils.setHistStyle(hSin2PhiMinusPsiVsPt, cent_colours[i_cent])
        vSin2PhiMinusPsiVsPt[det].append(hSin2PhiMinusPsiVsPt)

    for i_pt in range(0, n_pt_bins):
        # select the correct pt bin
        pt_sel = f"abs(fPt) > {pt_bins[i_cent][i_pt]} and abs(fPt) < {pt_bins[i_cent][i_pt+1]}"  # and {ptdep_selection_list[i_pt]}"
        bin_df = cent_df.query(pt_sel, inplace=False)

        hPhi = ROOT.TH1F(
            f"hPhi_cent_{cent[0]}_{cent[1]}_pt_{i_pt}", r";#phi; counts", 140, -7.0, 7.0
        )
        utils.setHistStyle(hPhi, ROOT.kRed + 2)

        for phi in bin_df["fPhi"]:
            hPhi.Fill(phi)

        vPhi_cent.append(hPhi)
        phi_dir.cd()
        hPhi.Write()

        for det in detectors:
            hPhiMinusPsi_formatted_title = f";phi - #Psi_{{{det}}}; counts"
            hPhiMinusPsi = ROOT.TH1F(
                f"hPhiMinusPsi_{det}_cent_{cent[0]}_{cent[1]}_pt_{i_pt}",
                hPhiMinusPsi_formatted_title,
                40,
                0.0,
                ROOT.TMath.Pi(),
            )
            utils.setHistStyle(hPhiMinusPsi, ROOT.kRed + 2)

            hSinPhiMinusPsi_formatted_title = rf";sin(#phi - #Psi_{{{det}}}); counts"
            hSinPhiMinusPsi = ROOT.TH1F(
                f"hSinPhiMinusPsi_{det}_cent_{cent[0]}_{cent[1]}_pt_{i_pt}",
                hSinPhiMinusPsi_formatted_title,
                100,
                -1,
                1,
            )
            utils.setHistStyle(hSinPhiMinusPsi, ROOT.kRed + 2)

            hCos2PhiMinusPsi_formatted_title = f";cos(2(#phi - #Psi_{{{det}}}); counts"
            hCos2PhiMinusPsi = ROOT.TH1F(
                f"hCos2PhiMinusPsi_{det}_cent_{cent[0]}_{cent[1]}_pt_{i_pt}",
                hCos2PhiMinusPsi_formatted_title,
                100,
                -1,
                1,
            )
            utils.setHistStyle(hCos2PhiMinusPsi, ROOT.kRed + 2)

            hSin2PhiMinusPsi_formatted_title = f";sin(2(#phi - #Psi_{{{det}}}); counts"
            hSin2PhiMinusPsi = ROOT.TH1F(
                f"hSin2PhiMinusPsi_{det}_cent_{cent[0]}_{cent[1]}_pt_{i_pt}",
                hSin2PhiMinusPsi_formatted_title,
                100,
                -1,
                1,
            )
            utils.setHistStyle(hSin2PhiMinusPsi, ROOT.kRed + 2)

            hQx_formatted_title = f";Q_{{x, {det}}}; counts"
            hQx = ROOT.TH1F(
                f"hQx_{det}_cent_{cent[0]}_{cent[1]}_pt_{i_pt}",
                hQx_formatted_title,
                500,
                -5,
                5,
            )
            utils.setHistStyle(hQx, ROOT.kRed + 2)

            hQy_formatted_title = f";Q_{{y, {det}}}; counts"
            hQy = ROOT.TH1F(
                f"hQy_{det}_cent_{cent[0]}_{cent[1]}_pt_{i_pt}",
                f";Qy_{det}; counts",
                500,
                -5,
                5,
            )
            utils.setHistStyle(hQy, ROOT.kRed + 2)

            hQ_formatted_title = f";Q_{{{det}}}; counts"
            hQ = ROOT.TH1F(
                f"hQ_{det}_cent_{cent[0]}_{cent[1]}_pt_{i_pt}",
                hQ_formatted_title,
                500,
                0,
                5,
            )
            utils.setHistStyle(hQ, ROOT.kRed + 2)

            hV2ScalarProduct_formatted_title = (
                f";Q_{{{det}}} cos(2(#phi - #Psi_{{{det}}}); counts"
            )
            hV2ScalarProduct = ROOT.TH1F(
                f"hV2ScalarProduct_{det}_cent_{cent[0]}_{cent[1]}_pt_{i_pt}",
                hV2ScalarProduct_formatted_title,
                200,
                -4,
                4,
            )
            utils.setHistStyle(hV2ScalarProduct, ROOT.kRed + 2)

            hV2ScalarProduct_check_formatted_title = (
                f";Q_{{x, {det}}} u_{{x}} + Q_{{y, {det}}} u_{{y}}); counts"
            )
            hV2ScalarProduct_check = ROOT.TH1F(
                f"hV2ScalarProduct_check_{det}_cent_{cent[0]}_{cent[1]}_pt_{i_pt}",
                hV2ScalarProduct_check_formatted_title,
                200,
                -4,
                4,
            )
            utils.setHistStyle(hV2ScalarProduct_check, ROOT.kAzure + 1)

            hV2Tree = ROOT.TH1F(
                f"hV2Tree_{det}_cent_{cent[0]}_{cent[1]}_pt_{i_pt}",
                hV2ScalarProduct_formatted_title,
                200,
                -4,
                4,
            )
            utils.setHistStyle(hV2Tree, ROOT.kGreen + 3)

            for phi, q, psi, v2 in zip(
                bin_df["fPhi"],
                bin_df[f"fQ{det}"],
                bin_df[f"fPsi{det}"],
                bin_df[f"fV2{det}"],
            ):
                delta_phi = utils.getPhiInRange(phi - psi)
                sinPhiMinusPsi = ROOT.TMath.Sin(delta_phi)
                cos2PhiMinusPsi = ROOT.TMath.Cos(2 * delta_phi)
                sin2PhiMinusPsi = ROOT.TMath.Sin(2 * delta_phi)
                hPhiMinusPsi.Fill(delta_phi)
                hSinPhiMinusPsi.Fill(sinPhiMinusPsi)
                hCos2PhiMinusPsi.Fill(cos2PhiMinusPsi)
                hSin2PhiMinusPsi.Fill(sin2PhiMinusPsi)
                hQ.Fill(q)
                qx = q * ROOT.TMath.Cos(2 * psi)
                hQx.Fill(qx)
                qy = q * ROOT.TMath.Sin(2 * psi)
                hQy.Fill(qy)
                hV2ScalarProduct.Fill(q * cos2PhiMinusPsi)
                ux = ROOT.TMath.Cos(2 * phi)
                uy = ROOT.TMath.Sin(2 * phi)
                v2_check = qx * ux + qy * uy
                hV2ScalarProduct_check.Fill(v2_check)
                hV2Tree.Fill(v2)

            vPhiMinusPsi[det][i_cent].append(hPhiMinusPsi)
            vSinPhiMinusPsi[det][i_cent].append(hSinPhiMinusPsi)
            vCos2PhiMinusPsi[det][i_cent].append(hCos2PhiMinusPsi)
            vSin2PhiMinusPsi[det][i_cent].append(hSin2PhiMinusPsi)
            vQ[det][i_cent].append(hQ)
            vQx[det][i_cent].append(hQx)
            vQy[det][i_cent].append(hQy)
            vV2ScalarProduct[det][i_cent].append(hV2ScalarProduct)
            vV2ScalarProduct_check[det][i_cent].append(hV2ScalarProduct_check)

            vSinPhiMinusPsiVsPt[det][i_cent].SetBinContent(
                i_pt + 1, hSinPhiMinusPsi.GetMean()
            )
            vSinPhiMinusPsiVsPt[det][i_cent].SetBinError(
                i_pt + 1, hSinPhiMinusPsi.GetMeanError()
            )
            vCos2PhiMinusPsiVsPt[det][i_cent].SetBinContent(
                i_pt + 1, hCos2PhiMinusPsi.GetMean()
            )
            vCos2PhiMinusPsiVsPt[det][i_cent].SetBinError(
                i_pt + 1, hCos2PhiMinusPsi.GetMeanError()
            )
            vSin2PhiMinusPsiVsPt[det][i_cent].SetBinContent(
                i_pt + 1, hSin2PhiMinusPsi.GetMean()
            )
            vSin2PhiMinusPsiVsPt[det][i_cent].SetBinError(
                i_pt + 1, hSin2PhiMinusPsi.GetMeanError()
            )

            cent_dir.cd(det)
            hPhiMinusPsi.Write()
            hSinPhiMinusPsi.Write()
            hCos2PhiMinusPsi.Write()
            hSin2PhiMinusPsi.Write()
            hQ.Write()
            hQx.Write()
            hQy.Write()
            hV2ScalarProduct.Write()
            hV2ScalarProduct_check.Write()
            hV2Tree.Write()

    for det in detectors:
        cent_dir.cd(det)
        vSinPhiMinusPsiVsPt[det][i_cent].Write()
        vCos2PhiMinusPsiVsPt[det][i_cent].Write()
        vSin2PhiMinusPsiVsPt[det][i_cent].Write()

    vPhi.append(vPhi_cent)

for det in detectors:
    cSinPhiMinusPsiVsPt[det] = ROOT.TCanvas(
        f"cCos2PhiMinusPsiVsPt_{det}", f"cCos2PhiMinusPsiVsPt_{det}", 800, 600
    )
    cSinLegend = ROOT.TLegend(0.40, 0.22, 0.70, 0.45, det, "brNDC")
    cSinPhiMinusPsiVsPt[det].DrawFrame(
        0.0,
        -1.2,
        11.0,
        1.2,
        f";#it{{p}}_{{T}} (GeV/#it{{c}}); sin(#phi - #Psi_{{{det}}})",
    )
    for i_cent, cent in enumerate(centrality_classes):
        vSinPhiMinusPsiVsPt[det][i_cent].Draw("PE SAME")
        cSinLegend.AddEntry(
            vSinPhiMinusPsiVsPt[det][i_cent], f"{cent[0]}-{cent[1]}%", "PL"
        )
    cSinLegend.Draw("SAME")
    cSinPhiMinusPsiVsPt[det].SaveAs(
        f"{output_dir_name}/qc_plots/cSinPhiMinusPsiVsPt_{det}.pdf"
    )
    output_file.cd(f"{det}")
    cSinPhiMinusPsiVsPt[det].Write()

    cCos2PhiMinusPsiVsPt[det] = ROOT.TCanvas(
        f"cCos2PhiMinusPsiVsPt_{det}", f"cCos2PhiMinusPsiVsPt_{det}", 800, 600
    )
    cCos2Legend = ROOT.TLegend(0.40, 0.22, 0.70, 0.45, det, "brNDC")
    cCos2PhiMinusPsiVsPt[det].DrawFrame(
        0.0,
        -1.2,
        11.0,
        1.2,
        f";#it{{p}}_{{T}} (GeV/#it{{c}}); cos(2(#phi - #Psi_{{{det}}}))",
    )
    for i_cent, cent in enumerate(centrality_classes):
        vCos2PhiMinusPsiVsPt[det][i_cent].Draw("PE SAME")
        cCos2Legend.AddEntry(
            vCos2PhiMinusPsiVsPt[det][i_cent], f"{cent[0]}-{cent[1]}%", "PL"
        )
    cCos2Legend.Draw("SAME")
    cCos2PhiMinusPsiVsPt[det].SaveAs(
        f"{output_dir_name}/qc_plots/cCos2PhiMinusPsiVsPt_{det}.pdf"
    )
    output_file.cd(f"{det}")
    cCos2PhiMinusPsiVsPt[det].Write()

    cSin2PhiMinusPsiVsPt[det] = ROOT.TCanvas(
        f"cSin2PhiMinusPsiVsPt_{det}", f"cSin2PhiMinusPsiVsPt_{det}", 800, 600
    )
    cSin2Legend = ROOT.TLegend(0.40, 0.22, 0.70, 0.45, det, "brNDC")
    cSin2PhiMinusPsiVsPt[det].DrawFrame(
        0.0,
        -1.2,
        11.0,
        1.2,
        f";#it{{p}}_{{T}} (GeV/#it{{c}}); sin(2(#phi - #Psi_{{{det}}}))",
    )
    for i_cent, cent in enumerate(centrality_classes):
        vSin2PhiMinusPsiVsPt[det][i_cent].Draw("PE SAME")
        cSin2Legend.AddEntry(
            vSin2PhiMinusPsiVsPt[det][i_cent], f"{cent[0]}-{cent[1]}%", "PL"
        )
    cSin2Legend.Draw("SAME")
    cSin2PhiMinusPsiVsPt[det].SaveAs(
        f"{output_dir_name}/qc_plots/cSin2PhiMinusPsiVsPt_{det}.pdf"
    )
    output_file.cd(f"{det}")
    cSin2PhiMinusPsiVsPt[det].Write()

# filling QC-plots after pt-dependent selections
for i_pt in range(0, n_pt_bins_no_cent):
    # select the correct pt bin
    pt_sel = f"abs(fPt) > {pt_bins_no_cent[i_pt]} and abs(fPt) < {pt_bins_no_cent[i_pt+1]}"  # and {ptdep_selection_list[i_pt]}"
    bin_df = complete_df.query(pt_sel, inplace=False)

    # Fill QC histograms
    print("Filling eta")
    for eta in bin_df["fEta"]:
        hEta.Fill(eta)

    print("Filling specific energy loss")
    for rig, signal in zip(bin_df["fTPCInnerParam"], bin_df["fTPCsignal"]):
        hTPCsignalVsPoverZ.Fill(rig, signal)

    print("Filling TOF squared mass and TPC nsigma")
    for pt, tof_mass_squared, n_sigma_tpc in zip(
        bin_df["fPt"], bin_df["fTOFmassSquared"], bin_df["fNsigmaTPC3He"]
    ):
        hTOFmassSquaredVsPt.Fill(
            abs(pt), tof_mass_squared - utils.mass_alpha * utils.mass_alpha
        )
        hNsigmaVsPt.Fill(abs(pt), n_sigma_tpc)
        hNsigmaVsTOFmassSquared.Fill(
            n_sigma_tpc, tof_mass_squared - utils.mass_alpha * utils.mass_alpha
        )

functions = utils.getBBAfunctions(
    parameters=p_train, resolution=resolution_train, n_sigma=2
)

functions_alpha = utils.getBBAfunctions(
    parameters=p_train,
    resolution=resolution_train,
    n_sigma=3,
    color=ROOT.kGreen,
    mass=utils.mass_alpha,
)

cTPC = ROOT.TCanvas("cTPC", "cTPC", 800, 600)
hTPCsignalVsPoverZ.Draw("colz")

for f in functions[3:]:
    f.Draw("L SAME")

# Saving histograms to file
qc_dir = output_file.mkdir("QC")
qc_dir.cd()
hEta.Write()
hNsigmaITSvsP.Write()
its_nsigma_func.Write()
cNsigmaITSvsP.Write()

hTPCsignalVsPoverZ.Write()
hTOFmassSquaredVsPt.Write()
hNsigmaVsPt.Write()
hNsigmaVsTOFmassSquared.Write()
for f in functions:
    f.Write()
hZvtx.Write()
hCentFT0C.Write()
hEventSelections_normal.Write()
hEventSelections_igor.Write()

# Save QC histogram as PDF
plot_dir_name = f"{output_dir_name}/qc_plots"
if not os.path.exists(plot_dir_name):
    os.makedirs(plot_dir_name)

utils.saveCanvasAsPDF(hEta, plot_dir_name)
utils.saveCanvasAsPDF(hTPCsignalVsPoverZ, plot_dir_name, is2D=True)
utils.saveCanvasAsPDF(hTOFmassSquaredVsPt, plot_dir_name, is2D=True)
utils.saveCanvasAsPDF(hNsigmaVsPt, plot_dir_name, is2D=True)
utils.saveCanvasAsPDF(hNsigmaVsTOFmassSquared, plot_dir_name, is2D=True)
cTPC.SaveAs(f"{plot_dir_name}/cTPC.pdf")
# cAvgItsClusSizeCosLambdaPt.SaveAs(f"{plot_dir_name}/cAvgItsClusSizeCosLambdaPt.pdf")
utils.saveCanvasAsPDF(hZvtx, plot_dir_name)
utils.saveCanvasAsPDF(hCentFT0C, plot_dir_name)
utils.saveCanvasAsPDF(hEventSelections_normal, plot_dir_name)
utils.saveCanvasAsPDF(hEventSelections_igor, plot_dir_name)
for histo in vDeltaPsi.values():
    histo.Write()
    utils.saveCanvasAsPDF(histo, plot_dir_name, is2D=True)
# resolution plot
resolution_file = ROOT.TFile(resolution_file_name)
if useSP:
    res_histo_name = f"Resolution_SP/hResolution_{reference_flow_detector}_{resolution_flow_detectors[0]}_{resolution_flow_detectors[1]}_SP"
else:
    res_histo_name = f"Resolution_EP/hResolution_{reference_flow_detector}_{resolution_flow_detectors[0]}_{resolution_flow_detectors[1]}_EP"
hResolution = resolution_file.Get(res_histo_name)
hResolution.SetDirectory(0)
hResolution.SetTitle(r";FT0C percentile (%); R_{2}")
utils.setHistStyle(hResolution, ROOT.kRed + 2)
cResolution = ROOT.TCanvas(
    f"cResolution_{reference_flow_detector}_{resolution_flow_detectors[0]}_{resolution_flow_detectors[1]}",
    f"cResolution_{reference_flow_detector}_{resolution_flow_detectors[0]}_{resolution_flow_detectors[1]}",
    800,
    600,
)
cResolution.SetBottomMargin(0.15)
hResolution.Draw("PE")
cResolution.SaveAs(
    f"{output_dir_name}/qc_plots/cResolution_{reference_flow_detector}_{resolution_flow_detectors[0]}_{resolution_flow_detectors[1]}.pdf"
)
