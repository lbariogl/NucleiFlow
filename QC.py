import ROOT
import pandas as pd
import yaml
import argparse
import os
import re

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

useSP = config["useSP"]

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
pt_bins = config["qc_pt_bins"]
n_pt_bins = len(pt_bins) - 1
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
utils.redefineColumns(complete_df)

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

complete_df.query("fNsigmaITS3He > - 1.5", inplace=True)

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
hPhi = ROOT.TH1F("hPhi", r";#phi; counts", 140, -7.0, 7.0)
utils.setHistStyle(hPhi, ROOT.kRed + 2)
hPsiFT0C = ROOT.TH1F("hPsiFT0C", r";#Psi_{FTOC}; counts", 140, -7.0, 7.0)
utils.setHistStyle(hPsiFT0C, ROOT.kRed + 2)
hPhiMinusPsiFT0C = ROOT.TH1F(
    "hPhiMinusPsiFT0C", r";phi - #Psi_{FTOC}; counts", 140, -7.0, 7.0
)
utils.setHistStyle(hPhiMinusPsiFT0C, ROOT.kRed + 2)
hV2 = ROOT.TH1F("hV2", r";cos(2(#phi - #Psi_{FTOC}))", 100, -1, 1)
utils.setHistStyle(hV2, ROOT.kRed + 2)

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
hDeltaPsi_FT0A_TPCl = input_dir_AR_flow_ep.Get("hDeltaPsi_FT0A_TPCl_EP")
hDeltaPsi_FT0A_TPCr = input_dir_AR_flow_ep.Get("hDeltaPsi_FT0A_TPCr_EP")
hDeltaPsi_FT0C_TPCl = input_dir_AR_flow_ep.Get("hDeltaPsi_FT0C_TPCl_EP")
hDeltaPsi_FT0C_TPCr = input_dir_AR_flow_ep.Get("hDeltaPsi_FT0C_TPCr_EP")
hDeltaPsi_FT0C_FT0A = input_dir_AR_flow_ep.Get("hDeltaPsi_FT0C_FT0A_EP")
hDeltaPsi_TPCl_TPCr = input_dir_AR_flow_ep.Get("hDeltaPsi_TPCl_TPCr_EP")
hDeltaPsi_FT0A_TPC = input_dir_AR_flow_ep.Get("hDeltaPsi_FT0A_TPC_EP")
hDeltaPsi_FT0C_TPC = input_dir_AR_flow_ep.Get("hDeltaPsi_FT0C_TPC_EP")

input_file_AR_small = ROOT.TFile(
    "/data/lbariogl/flow/LHC23zzh_pass4_small_new/AnalysisResults_general.root"
)

input_file_AR_small_igor = ROOT.TFile(
    "/data/lbariogl/flow/LHC23zzh_pass4_small_new/AnalysisResults_igor.root"
)
hEventSelections_normal = input_file_AR_small.Get(
    "nuclei-spectra/spectra/hEventSelections"
)
utils.setHistStyle(hZvtx, ROOT.kRed + 2)
hEventSelections_igor = input_file_AR_small_igor.Get(
    "nuclei-spectra/spectra/hEventSelections"
)
utils.setHistStyle(hZvtx, ROOT.kAzure + 1)


# get histrogram for FT0C distrubution
hPsiFT0C_2D = input_dir_AR_flow_ep.Get("hPsi_FT0C_EP")
vPsiFT0C = []
yMaxPsi = 0
legend_psi_FT0C = ROOT.TLegend(
    0.40, 0.22, 0.70, 0.45, "FTOC centrality classes", "brNDC"
)
legend_psi_FT0C.SetNColumns(2)
FT0C_dir = output_file.mkdir("FT0C")

for i_cent, cent in enumerate(centrality_classes):
    left_bin = hPsiFT0C_2D.GetXaxis().FindBin(cent[0])
    right_bin = hPsiFT0C_2D.GetXaxis().FindBin(cent[1]) - 1
    hPsiFT0C = hPsiFT0C_2D.ProjectionY(
        f"hFT0C_{cent[0]}_{cent[1]}", left_bin, right_bin
    )
    utils.setHistStyle(hPsiFT0C, cent_colours[i_cent])
    hPsiFT0C.Scale(1.0 / hPsiFT0C.Integral())
    vPsiFT0C.append(hPsiFT0C)
    loc_max = hPsiFT0C.GetMaximum()
    if loc_max > yMaxPsi:
        yMaxPsi = loc_max
    legend_psi_FT0C.AddEntry(hPsiFT0C, f"{cent[0]}-{cent[1]}%", "l")

cPsiFT0C = ROOT.TCanvas("cPsiFT0C", "cPsiFT0C", 800, 600)
cPsiFT0C.DrawFrame(-3.5, 0, 3.5, yMaxPsi * 1.1, r";#Psi_{2}^{FT0C}; norm. counts")
for i_cent, hPsiFT0C in enumerate(vPsiFT0C):
    hPsiFT0C.Draw("L SAME")
    FT0C_dir.cd()
    hPsiFT0C.Write()
legend_psi_FT0C.Draw("SAME")
cPsiFT0C.SaveAs(f"{output_dir_name}/qc_plots/cPsiFT0C.pdf")
FT0C_dir.cd()
cPsiFT0C.Write()

# get histrogram for TPCl distrubution
hPsiTPCl_2D = input_dir_AR_flow_ep.Get("hPsi_TPCl_EP")
vPsiTPCl = []
yMaxPsi = 0
legend_psi_TPCl = ROOT.TLegend(
    0.40, 0.22, 0.70, 0.45, "TPCl centrality classes", "brNDC"
)
legend_psi_TPCl.SetNColumns(2)
TPCl_dir = output_file.mkdir("TPCl")

for i_cent, cent in enumerate(centrality_classes):
    left_bin = hPsiTPCl_2D.GetXaxis().FindBin(cent[0])
    right_bin = hPsiTPCl_2D.GetXaxis().FindBin(cent[1]) - 1
    hPsiTPCl = hPsiTPCl_2D.ProjectionY(
        f"hTPCl_{cent[0]}_{cent[1]}", left_bin, right_bin
    )
    utils.setHistStyle(hPsiTPCl, cent_colours[i_cent])
    hPsiTPCl.Scale(1.0 / hPsiTPCl.Integral())
    vPsiTPCl.append(hPsiTPCl)
    loc_max = hPsiTPCl.GetMaximum()
    if loc_max > yMaxPsi:
        yMaxPsi = loc_max
    legend_psi_TPCl.AddEntry(hPsiTPCl, f"{cent[0]}-{cent[1]}%", "l")

cPsiTPCl = ROOT.TCanvas("cPsiTPCl", "cPsiTPCl", 800, 600)
cPsiTPCl.DrawFrame(-3.5, 0, 3.5, yMaxPsi * 1.1, r";#Psi_{2}^{TPCl}; norm. counts")
for i_cent, hPsiTPCl in enumerate(vPsiTPCl):
    hPsiTPCl.Draw("L SAME")
    TPCl_dir.cd()
    hPsiTPCl.Write()
legend_psi_TPCl.Draw("SAME")
cPsiTPCl.SaveAs(f"{output_dir_name}/qc_plots/cPsiTPCl.pdf")
TPCl_dir.cd()
cPsiTPCl.Write()

# get histrogram for TPCr distrubution
hPsiTPCr_2D = input_dir_AR_flow_ep.Get("hPsi_TPCr_EP")
vPsiTPCr = []
yMaxPsi = 0
legend_psi_TPCr = ROOT.TLegend(
    0.40, 0.22, 0.70, 0.45, "TPCr centrality classes", "brNDC"
)
legend_psi_TPCr.SetNColumns(2)
TPCr_dir = output_file.mkdir("TPCr")

for i_cent, cent in enumerate(centrality_classes):
    left_bin = hPsiTPCr_2D.GetXaxis().FindBin(cent[0])
    right_bin = hPsiTPCr_2D.GetXaxis().FindBin(cent[1]) - 1
    hPsiTPCr = hPsiTPCr_2D.ProjectionY(
        f"hTPCr_{cent[0]}_{cent[1]}", left_bin, right_bin
    )
    utils.setHistStyle(hPsiTPCr, cent_colours[i_cent])
    hPsiTPCr.Scale(1.0 / hPsiTPCr.Integral())
    vPsiTPCr.append(hPsiTPCr)
    loc_max = hPsiTPCr.GetMaximum()
    if loc_max > yMaxPsi:
        yMaxPsi = loc_max
    legend_psi_TPCr.AddEntry(hPsiTPCr, f"{cent[0]}-{cent[1]}%", "l")

cPsiTPCr = ROOT.TCanvas("cPsiTPCr", "cPsiTPCr", 800, 600)
cPsiTPCr.DrawFrame(-3.5, 0, 3.5, yMaxPsi * 1.1, r";#Psi_{2}^{TPCr}; norm. counts")
for i_cent, hPsiTPCr in enumerate(vPsiTPCr):
    hPsiTPCr.Draw("L SAME")
    TPCr_dir.cd()
    hPsiTPCr.Write()
legend_psi_TPCr.Draw("SAME")
cPsiTPCr.SaveAs(f"{output_dir_name}/qc_plots/cPsiTPCr.pdf")
TPCr_dir.cd()
cPsiTPCr.Write()

# get histrogram for TPC distrubution
hPsiTPC_2D = input_dir_AR_flow_ep.Get("hPsi_TPC_EP")
vPsiTPC = []
yMaxPsi = 0
legend_psi_TPC = ROOT.TLegend(0.40, 0.22, 0.70, 0.45, "TPC centrality classes", "brNDC")
legend_psi_TPC.SetNColumns(2)
TPC_dir = output_file.mkdir("TPC")

for i_cent, cent in enumerate(centrality_classes):
    left_bin = hPsiTPC_2D.GetXaxis().FindBin(cent[0])
    right_bin = hPsiTPC_2D.GetXaxis().FindBin(cent[1]) - 1
    hPsiTPC = hPsiTPC_2D.ProjectionY(f"hTPC_{cent[0]}_{cent[1]}", left_bin, right_bin)
    utils.setHistStyle(hPsiTPC, cent_colours[i_cent])
    hPsiTPC.Scale(1.0 / hPsiTPC.Integral())
    vPsiTPC.append(hPsiTPC)
    loc_max = hPsiTPC.GetMaximum()
    if loc_max > yMaxPsi:
        yMaxPsi = loc_max
    legend_psi_TPC.AddEntry(hPsiTPC, f"{cent[0]}-{cent[1]}%", "l")

cPsiTPC = ROOT.TCanvas("cPsiTPCr", "cPsiTPCr", 800, 600)
cPsiTPC.DrawFrame(-3.5, 0, 3.5, yMaxPsi * 1.1, r";#Psi_{2}^{TPCr}; norm. counts")
for i_cent, hPsiTPCr in enumerate(vPsiTPCr):
    hPsiTPC.Draw("L SAME")
    TPC_dir.cd()
    hPsiTPC.Write()
legend_psi_TPCr.Draw("SAME")
cPsiTPC.SaveAs(f"{output_dir_name}/qc_plots/cPsiTPC.pdf")
TPC_dir.cd()
cPsiTPC.Write()


# filling QC-plots after pt-dependent selections
for i_pt in range(0, n_pt_bins):
    # select the correct pt bin
    pt_sel = f"abs(fPt) > {pt_bins[i_pt]} and abs(fPt) < {pt_bins[i_pt+1]}"  # and {ptdep_selection_list[i_pt]}"
    bin_df = complete_df.query(pt_sel, inplace=False)

    # Fill QC histograms
    print("Filling eta")
    for eta in bin_df["fEta"]:
        hEta.Fill(eta)

    # print("Filling cluster size")
    # for avgClus in bin_df["fAvgItsClusSize"]:
    #     hAvgItsClusSize.Fill(avgClus)

    # print("Filling cluster size * cos(lambda)")
    # for avgClus in bin_df["fAvgItsClusSizeCosLambda"]:
    #     hAvgItsClusSizeCosLambdaIntegrated.Fill(avgClus)

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

    print("Filling phi")
    for phi, psi_ft0c in zip(bin_df["fPhi"], bin_df["fPsiFT0C"]):
        hPhi.Fill(phi)
        hPsiFT0C.Fill(psi_ft0c)
        delta_phi = phi - psi_ft0c
        hPhiMinusPsiFT0C.Fill(delta_phi)
        hV2.Fill(ROOT.TMath.Cos(2 * delta_phi))

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

# for f in functions_alpha:
#     f.Draw("L SAME")


# Saving histograms to file
qc_dir = output_file.mkdir("QC")
qc_dir.cd()
hEta.Write()
hNsigmaITSvsP.Write()
# hAvgItsClusSize.Write()
# hAvgItsClusSizeCosLambdaIntegrated.Write()
# cAvgItsClusSizeCosLambdaPt.Write()
# for i_pt in range(0, n_pt_bins):
#     hAvgItsClusSizeCosLambda[i_pt].Write()
hTPCsignalVsPoverZ.Write()
hTOFmassSquaredVsPt.Write()
hNsigmaVsPt.Write()
hNsigmaVsTOFmassSquared.Write()
hPhi.Write()
hPsiFT0C.Write()
hPhiMinusPsiFT0C.Write()
cTPC.Write()
hV2.Write()
for f in functions:
    f.Write()
hZvtx.Write()
hCentFT0C.Write()
hEventSelections_normal.Write()
hEventSelections_igor.Write()
hDeltaPsi_FT0A_TPCl.Write()
hDeltaPsi_FT0A_TPCr.Write()
hDeltaPsi_FT0C_TPCl.Write()
hDeltaPsi_FT0C_TPCr.Write()
hDeltaPsi_FT0C_FT0A.Write()
hDeltaPsi_TPCl_TPCr.Write()
hDeltaPsi_FT0A_TPC.Write()
hDeltaPsi_FT0C_TPC.Write()

# Save QC histogram as PDF
plot_dir_name = f"{output_dir_name}/qc_plots"
if not os.path.exists(plot_dir_name):
    os.makedirs(plot_dir_name)

utils.saveCanvasAsPDF(hEta, plot_dir_name)
# utils.saveCanvasAsPDF(hAvgItsClusSize, plot_dir_name)
# utils.saveCanvasAsPDF(hAvgItsClusSizeCosLambdaIntegrated, plot_dir_name)
# for i_pt in range(0, n_pt_bins):
#     utils.saveCanvasAsPDF(hAvgItsClusSizeCosLambda[i_pt], plot_dir_name)
utils.saveCanvasAsPDF(hTPCsignalVsPoverZ, plot_dir_name, is2D=True)
utils.saveCanvasAsPDF(hTOFmassSquaredVsPt, plot_dir_name, is2D=True)
utils.saveCanvasAsPDF(hNsigmaVsPt, plot_dir_name, is2D=True)
utils.saveCanvasAsPDF(hNsigmaITSvsP, plot_dir_name, is2D=True)
utils.saveCanvasAsPDF(hNsigmaVsTOFmassSquared, plot_dir_name, is2D=True)
utils.saveCanvasAsPDF(hPhi, plot_dir_name)
utils.saveCanvasAsPDF(hPsiFT0C, plot_dir_name)
utils.saveCanvasAsPDF(hPhiMinusPsiFT0C, plot_dir_name)
utils.saveCanvasAsPDF(hV2, plot_dir_name)
cTPC.SaveAs(f"{plot_dir_name}/cTPC.pdf")
# cAvgItsClusSizeCosLambdaPt.SaveAs(f"{plot_dir_name}/cAvgItsClusSizeCosLambdaPt.pdf")
utils.saveCanvasAsPDF(hZvtx, plot_dir_name)
utils.saveCanvasAsPDF(hCentFT0C, plot_dir_name)
utils.saveCanvasAsPDF(hEventSelections_normal, plot_dir_name)
utils.saveCanvasAsPDF(hEventSelections_igor, plot_dir_name)
utils.saveCanvasAsPDF(hDeltaPsi_FT0A_TPCl, plot_dir_name, is2D=True)
utils.saveCanvasAsPDF(hDeltaPsi_FT0A_TPCr, plot_dir_name, is2D=True)
utils.saveCanvasAsPDF(hDeltaPsi_FT0C_TPCl, plot_dir_name, is2D=True)
utils.saveCanvasAsPDF(hDeltaPsi_FT0C_TPCr, plot_dir_name, is2D=True)
utils.saveCanvasAsPDF(hDeltaPsi_FT0C_FT0A, plot_dir_name, is2D=True)
utils.saveCanvasAsPDF(hDeltaPsi_TPCl_TPCr, plot_dir_name, is2D=True)
utils.saveCanvasAsPDF(hDeltaPsi_FT0A_TPC, plot_dir_name, is2D=True)
utils.saveCanvasAsPDF(hDeltaPsi_FT0C_TPC, plot_dir_name, is2D=True)

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
