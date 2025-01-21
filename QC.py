import ROOT
from hipe4ml.tree_handler import TreeHandler
import pandas as pd
import yaml
import argparse
import os

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

input_file_name = config["input_file_name"]
input_file_AR_name = config["input_file_AR_name"]
resolution_file_name = config["resolution_file_name"]
output_dir_name = config["output_dir_name"]
output_file_name = config["output_file_name_qc"]

nuclei_tree_name = config["nuclei_tree_name"]
ep_tree_name = config["ep_tree_name"]

mandatory_selections = config["mandatory_selections"]
selection_dict = config["selection_dict"]
selection_list = selection_dict.values()
selections = " and ".join(selection_list)
ptdep_selection_dict = config["ptdep_selection_dict"]["fAvgItsClusSizeCosLambda"]

cent_detector_label = config["cent_detector_label"]

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
nuclei_hdl = TreeHandler(input_file_name, f"{nuclei_tree_name};", folder_name="DF*")
nuclei_df = nuclei_hdl._full_data_frame

nucleiflow_hdl = TreeHandler(input_file_name, f"{ep_tree_name};", folder_name="DF*")
nucleiflow_df = nucleiflow_hdl._full_data_frame

complete_df = pd.concat([nuclei_df, nucleiflow_df], axis=1, join="inner")

# set input files for EP qc
input_file_AR = ROOT.TFile(input_file_AR_name)
input_dir_AR_general = input_file_AR.Get("flow-qc/general")
input_dir_AR_flow_ep = input_file_AR.Get("flow-qc/flow_ep")

# define new columns
utils.redefineColumns(complete_df)

# apply common selections
complete_df.query(f"{mandatory_selections} and {selections}", inplace=True)

# Create QC histograms
hEta = ROOT.TH1F("hEta", ";#eta;", 200, -1.0, 1.0)
utils.setHistStyle(hEta, ROOT.kRed + 2)
hAvgItsClusSize = ROOT.TH1F(
    "hAvgItsClusSize", r";#LT ITS cluster size #GT; counts", 100, 0, 20
)
utils.setHistStyle(hAvgItsClusSize, ROOT.kRed + 2)

hAvgItsClusSizeCosLambdaIntegrated = ROOT.TH1F(
    "hAvgItsClusSizeCosLambdaIntegrated",
    r";#LT ITS cluster size #GT #times Cos(#lambda); counts",
    100,
    0,
    20,
)
utils.setHistStyle(hAvgItsClusSizeCosLambdaIntegrated, ROOT.kRed + 2)
hAvgItsClusSizeCosLambda = []

cols = ROOT.TColor.GetPalette()
period = int(cols.GetSize() / n_pt_bins)
for i_pt in range(0, n_pt_bins):
    pt_label = (
        f"{pt_bins[i_pt]:.1f}"
        + r" #leq #it{p}_{T} < "
        + f"{pt_bins[i_pt+1]:.1f}"
        + r" GeV/#it{c}"
    )
    hAvgItsClusSizeCosLambda_tmp = ROOT.TH1F(
        f"hAvgItsClusSizeCosLambda_pt{i_pt}",
        pt_label + r";#LT ITS cluster size #GT #times Cos(#lambda); counts",
        100,
        0,
        20,
    )
    utils.setHistStyle(
        hAvgItsClusSizeCosLambda_tmp, cols.At(i_pt * period), linewidth=2
    )
    hAvgItsClusSizeCosLambda.append(hAvgItsClusSizeCosLambda_tmp)

hNsigmaITSvsP = ROOT.TH2F(
    "hNsigmaITSvsP",
    r";n#sigma_{ITS}; #it{p} (GeV/#it{c})",
    120,
    0.0,
    12.0,
    500,
    -5.0,
    5.0,
)

for p, n_sigma_its in zip(complete_df["fP"], complete_df["fNsigmaITS3He"]):
    hNsigmaITSvsP.Fill(p, n_sigma_its)

hTPCsignalVsPoverZ = ROOT.TH2F(
    "hTPCsignalVsPoverZ",
    r";#it{p}/z (GeV/#it{c}); d#it{E} / d#it{x} (a.u.)",
    600,
    -6.0,
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
# add line corresponding to expected 4He squared mass


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

cAvgItsClusSizeCosLambdaPt = ROOT.TCanvas(
    "cAvgItsClusSizeCosLambdaPt", "cAvgItsClusSizeCosLambdaPt", 800, 600
)
cAvgItsClusSizeCosLambdaPt.DrawFrame(
    0, 0, 20, 4000, r";#LT ITS cluster size #GT #times Cos(#lambda); counts"
)
legend = ROOT.TLegend(0.6, 0.61, 0.86, 0.88, "", "brNDC")
legend.SetBorderSize(0)

for i_pt in range(0, n_pt_bins):
    # select the correct pt bin
    pt_sel = f"abs(fPt) > {pt_bins[i_pt]} and abs(fPt) < {pt_bins[i_pt+1]}"
    bin_df = complete_df.query(pt_sel, inplace=False)
    for avgClus in bin_df["fAvgItsClusSizeCosLambda"]:
        hAvgItsClusSizeCosLambda[i_pt].Fill(avgClus)
    cAvgItsClusSizeCosLambdaPt.cd()
    hAvgItsClusSizeCosLambda[i_pt].Draw("same")
    legend.AddEntry(
        hAvgItsClusSizeCosLambda[i_pt], hAvgItsClusSizeCosLambda[i_pt].GetTitle(), "PE"
    )
legend.Draw()

# define list of pt-dependent selections
ptdep_selection_list = []
for i in range(0, n_pt_bins):
    bin_centre = (pt_bins[i + 1] + pt_bins[i]) / 2
    condition = utils.get_condition(bin_centre, ptdep_selection_dict)
    ptdep_selection_list.append(condition)

# filling QC-plots after pt-dependent selections
for i_pt in range(0, n_pt_bins):
    # select the correct pt bin
    pt_sel = f"abs(fPt) > {pt_bins[i_pt]} and abs(fPt) < {pt_bins[i_pt+1]} and {ptdep_selection_list[i_pt]}"
    bin_df = complete_df.query(pt_sel, inplace=False)

    # Fill QC histograms
    print("Filling eta")
    for eta in bin_df["fEta"]:
        hEta.Fill(eta)

    print("Filling cluster size")
    for avgClus in bin_df["fAvgItsClusSize"]:
        hAvgItsClusSize.Fill(avgClus)

    print("Filling cluster size * cos(lambda)")
    for avgClus in bin_df["fAvgItsClusSizeCosLambda"]:
        hAvgItsClusSizeCosLambdaIntegrated.Fill(avgClus)

    print("Filling specific energy loss")
    for rig, sign, signal in zip(
        bin_df["fTPCInnerParam"], bin_df["fSign"], bin_df["fTPCsignal"]
    ):
        hTPCsignalVsPoverZ.Fill(sign * rig, signal)

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
    parameters=p_train, resolution=resolution_train, n_sigma=1
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

for f in functions:
    f.Draw("L SAME")

for f in functions_alpha:
    f.Draw("L SAME")


# Saving histograms to file
qc_dir = output_file.mkdir("QC")
qc_dir.cd()
hEta.Write()
hNsigmaITSvsP.Write()
hAvgItsClusSize.Write()
hAvgItsClusSizeCosLambdaIntegrated.Write()
cAvgItsClusSizeCosLambdaPt.Write()
for i_pt in range(0, n_pt_bins):
    hAvgItsClusSizeCosLambda[i_pt].Write()
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
hDeltaPsi_FT0A_TPCl.Write()
hDeltaPsi_FT0A_TPCr.Write()
hDeltaPsi_FT0C_TPCl.Write()
hDeltaPsi_FT0C_TPCr.Write()
hDeltaPsi_FT0C_FT0A.Write()
hDeltaPsi_TPCl_TPCr.Write()

# Save QC histogram as PDF
plot_dir_name = f"{output_dir_name}/qc_plots"
if not os.path.exists(plot_dir_name):
    os.makedirs(plot_dir_name)

utils.saveCanvasAsPDF(hEta, plot_dir_name)
utils.saveCanvasAsPDF(hAvgItsClusSize, plot_dir_name)
utils.saveCanvasAsPDF(hAvgItsClusSizeCosLambdaIntegrated, plot_dir_name)
for i_pt in range(0, n_pt_bins):
    utils.saveCanvasAsPDF(hAvgItsClusSizeCosLambda[i_pt], plot_dir_name)
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
cAvgItsClusSizeCosLambdaPt.SaveAs(f"{plot_dir_name}/cAvgItsClusSizeCosLambdaPt.pdf")
utils.saveCanvasAsPDF(hZvtx, plot_dir_name)
utils.saveCanvasAsPDF(hCentFT0C, plot_dir_name)
utils.saveCanvasAsPDF(hDeltaPsi_FT0A_TPCl, plot_dir_name, is2D=True, logScale=True)
utils.saveCanvasAsPDF(hDeltaPsi_FT0A_TPCr, plot_dir_name, is2D=True, logScale=True)
utils.saveCanvasAsPDF(hDeltaPsi_FT0C_TPCl, plot_dir_name, is2D=True, logScale=True)
utils.saveCanvasAsPDF(hDeltaPsi_FT0C_TPCr, plot_dir_name, is2D=True, logScale=True)
utils.saveCanvasAsPDF(hDeltaPsi_FT0C_FT0A, plot_dir_name, is2D=True, logScale=True)
utils.saveCanvasAsPDF(hDeltaPsi_TPCl_TPCr, plot_dir_name, is2D=True, logScale=True)

# resolution plot
resolution_file = ROOT.TFile(resolution_file_name)
hResolution = resolution_file.Get("Resolution/hResolution_FT0C_TPCl_TPCr")
hResolution.SetDirectory(0)
hResolution.SetTitle(r";FT0C percentile (%); R_{2}")
utils.setHistStyle(hResolution, ROOT.kRed + 2)
cResolutionFT0C = ROOT.TCanvas("cResolutionFT0C", "cResolutionFT0C", 800, 600)
cResolutionFT0C.SetBottomMargin(0.15)
hResolution.Draw("PE")
cResolutionFT0C.SaveAs(f"{output_dir_name}/qc_plots/cResolutionFT0C.pdf")
