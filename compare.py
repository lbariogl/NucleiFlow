import ROOT
import sys
from array import array

sys.path.append("utils")
import utils as utils

old_centrality_classes = [[0, 10], [10, 20], [20, 40], [40, 60]]

old_v2_stat = []
old_v2_syst = []

old_file = ROOT.TFile("../results_pass3_new/final.root")

for cent_limits in old_centrality_classes:
    cent_dir = old_file.Get(f"cent_{cent_limits[0]}_{cent_limits[1]}")
    old_v2_stat.append(cent_dir.Get(f"hV2vsPt_cent_{cent_limits[0]}_{cent_limits[1]}"))
    old_v2_syst.append(
        cent_dir.Get(f"hV2vsPt_cent_{cent_limits[0]}_{cent_limits[1]}_syst")
    )

new_centrality_classes = [
    [0, 10],
    [10, 20],
    [20, 30],
    [30, 40],
    [40, 50],
    [50, 60],
    [60, 80],
]

new_cent_colours = [633, 797, 418, 402, 433, 862, 874]

new_v2_stat = []
new_file = ROOT.TFile("../results_pass4_new/flow.root")

for cent_limits in new_centrality_classes:
    cent_dir = new_file.Get(f"cent_{cent_limits[0]}_{cent_limits[1]}")
    hV2_tmp = cent_dir.Get(f"default/hV2vsPt_cent_{cent_limits[0]}_{cent_limits[1]}")
    hV2_tmp.SetMarkerStyle(25)
    new_v2_stat.append(hV2_tmp)

output_file = ROOT.TFile("../results_pass4_new/comp_pass3_pass4.root", "RECREATE")

cV2 = ROOT.TCanvas("cV2_comp", "cV2_comp", 800, 600)
info_panel_total = ROOT.TPaveText(0.16, 0.71, 0.36, 0.87, "NDC")
info_panel_total.SetBorderSize(0)
info_panel_total.SetFillStyle(0)
info_panel_total.SetTextAlign(12)
info_panel_total.SetTextFont(42)
info_panel_total.SetTextSize(0.04)
info_panel_total.AddText(r"ALICE")
info_panel_total.AddText(r"Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.36 TeV")
info_panel_total.AddText(r"{}^{3}#bar{He}, |#eta| < 0.8")
frame = cV2.DrawFrame(1.7, -0.1, 8.5, 1.1, r";#it{p}_{T} (GeV/#it{c}); #it{v}_{2}{EP}")
cV2.SetBottomMargin(0.13)
cV2.SetLeftMargin(0.13)
cV2.SetBorderSize(0)
cV2.cd()

new_legend = ROOT.TLegend(0.17, 0.60, 0.46, 0.69, "pass4", "brNDC")
new_legend.SetBorderSize(0)
new_legend.SetNColumns(3)

new_legend.AddEntry(new_v2_stat[0], r"0-10 %", "PE")
new_v2_stat[0].Draw("PE SAME")
new_legend.AddEntry(new_v2_stat[4], r"40-50 %", "PE")
new_v2_stat[4].Draw("PE SAME")
new_legend.AddEntry(new_v2_stat[5], r"50-60 %", "PE")
new_v2_stat[5].Draw("PE SAME")

# for i_cent, cent_limits in enumerate(new_centrality_classes):
#     cent_label = f"{cent_limits[0]}-{cent_limits[1]}" + r"%"
#     new_legend.AddEntry(new_v2_stat[i_cent], cent_label, "P")
#     new_v2_stat[i_cent].Draw("PE SAME")

old_legend = ROOT.TLegend(0.17, 0.51, 0.35, 0.60, "pass3", "brNDC")
old_legend.SetBorderSize(0)
old_legend.SetNColumns(2)

# for i_cent, cent_limits in enumerate(old_centrality_classes):
#     cent_label = f"{cent_limits[0]}-{cent_limits[1]}" + r"%"
#     old_legend.AddEntry(old_v2_stat[i_cent], cent_label, "PF")
#     old_v2_stat[i_cent].Draw("PEX0 SAME")
#     old_v2_syst[i_cent].Draw("PE2 SAME")

old_legend.AddEntry(old_v2_stat[0], r"0-10 %", "PF")
utils.setHistStyle(old_v2_stat[0], 418)
utils.setHistStyle(old_v2_syst[0], 418)
old_v2_stat[0].Draw("PEX0 SAME")
old_v2_syst[0].Draw("PE2 SAME")

old_legend.AddEntry(old_v2_stat[3], r"40-60 %", "PF")
utils.setHistStyle(old_v2_stat[3], 797)
utils.setHistStyle(old_v2_syst[3], 797)
old_v2_stat[3].Draw("PEX0 SAME")
old_v2_syst[3].Draw("PE2 SAME")

old_legend.Draw()
new_legend.Draw()
info_panel_total.Draw()

output_file.cd()
cV2.Write()
cV2.SaveAs(f"../results_pass4_new/plots/{cV2.GetName()}.pdf")


def getHistos(file_name, histo_name, color, marker):
    # Open the ROOT file
    input_file = ROOT.TFile.Open(file_name)

    if input_file.IsZombie():
        raise ValueError(f"Failed to open file {file_name}")

    # Automatically find the TDirectory
    key_list = input_file.GetListOfKeys()
    t_directory = None
    for key in key_list:
        obj = key.ReadObj()
        if obj.InheritsFrom("TDirectory"):
            t_directory = obj
            break

    if not t_directory:
        raise ValueError("No TDirectory found in the ROOT file.")

    # Debug: print contents of the TDirectory
    # print("Contents of the TDirectory:")
    # t_directory.ls()

    # Dynamically get the exact names of the histograms
    def find_key(name):
        for key in t_directory.GetListOfKeys():
            if name in key.GetName():
                # print(f"Found histogram key: {key.GetName()}")
                return key.GetName()
        raise ValueError(f"Histogram '{name}' not found in the TDirectory.")

    # Access histograms within the TDirectory
    try:
        hCentralValues = t_directory.Get(find_key("Hist1D_y1"))
        hStatErrors = t_directory.Get(find_key("Hist1D_y1_e1"))
        hSystErrors = t_directory.Get(find_key("Hist1D_y1_e2"))
    except Exception as e:
        raise ValueError(f"Error accessing histograms: {e}")

    # Debug: Check if histograms are loaded correctly
    if not hCentralValues or not hStatErrors or not hSystErrors:
        raise ValueError(
            "One or more required histograms are missing in the TDirectory."
        )

    # Clone and set up the histograms
    n_bins = hCentralValues.GetNbinsX()
    hStat = hCentralValues.Clone(f"{histo_name}_stat")
    hStat.SetDirectory(0)
    utils.setHistStyle(hStat, color, marker)
    hStat.SetTitle(r";#it{p}_{T} (GeV/#it{c}); #it{v}_{2}")
    hSyst = hCentralValues.Clone(f"{histo_name}_syst")
    hSyst.SetTitle(r";#it{p}_{T} (GeV/#it{c}); #it{v}_{2}")
    hSyst.SetDirectory(0)
    utils.setHistStyle(hSyst, color, marker)

    # Loop through bins and set errors
    for i_bin in range(1, n_bins + 1):
        stat_err = hStatErrors.GetBinContent(i_bin)
        hStat.SetBinError(i_bin, stat_err)
        syst_err = hSystErrors.GetBinContent(i_bin)
        hSyst.SetBinError(i_bin, syst_err)

    # Unset directory association
    hCentralValues.SetDirectory(0)
    hStatErrors.SetDirectory(0)
    hSystErrors.SetDirectory(0)

    return hStat, hSyst


proton_1020_file_name = "../results_pass4/v2_proton/proton_1020.root"
hProton_1020_stat, hProton_1020_syst = getHistos(
    proton_1020_file_name, "hProton_1020", color=797, marker=25
)
proton_2030_file_name = "../results_pass4/v2_proton/proton_2030.root"
hProton_2030_stat, hProton_2030_syst = getHistos(
    proton_2030_file_name, "hProton_2030", color=418, marker=25
)
proton_3040_file_name = "../results_pass4/v2_proton/proton_3040.root"
hProton_3040_stat, hProton_3040_syst = getHistos(
    proton_3040_file_name, "hProton_3040", color=402, marker=21
)
proton_4050_file_name = "../results_pass4/v2_proton/proton_4050.root"
hProton_4050_stat, hProton_4050_syst = getHistos(
    proton_4050_file_name, "hProton_4050", color=433, marker=25
)
proton_5060_file_name = "../results_pass4/v2_proton/proton_5060.root"
hProton_5060_stat, hProton_5060_syst = getHistos(
    proton_5060_file_name, "hProton_5060", color=862, marker=25
)

proton_dir = output_file.mkdir("proton_comp")
proton_dir.cd()

hProton_1020_stat.Write()
hProton_1020_syst.Write()
hProton_2030_stat.Write()
hProton_2030_syst.Write()
hProton_3040_stat.Write()
hProton_3040_syst.Write()
hProton_4050_stat.Write()
hProton_4050_syst.Write()
hProton_5060_stat.Write()
hProton_5060_syst.Write()


def rescale_bins(original_hist, scale_factor=3, hist_name_suffix="_scaled"):
    """
    Create a new histogram with bin limits divided by a given scale factor.

    Args:
        original_hist (ROOT.TH1): The original histogram.
        scale_factor (float): The factor by which to divide the bin edges.
        hist_name_suffix (str): Suffix for the name of the new histogram.

    Returns:
        ROOT.TH1: A new histogram with scaled bin edges.
    """
    # Get the original histogram's binning information
    n_bins = original_hist.GetNbinsX()
    x_min = original_hist.GetXaxis().GetXmin()
    x_max = original_hist.GetXaxis().GetXmax()

    # Create new bin edges
    new_bin_edges = [
        x_min / scale_factor + i * ((x_max - x_min) / (scale_factor * n_bins))
        for i in range(n_bins + 1)
    ]

    # Create a new histogram with adjusted bin edges
    hist_name = original_hist.GetName() + hist_name_suffix
    new_hist = ROOT.TH1F(
        hist_name, original_hist.GetTitle(), n_bins, array("d", new_bin_edges)
    )
    new_hist.SetDirectory(0)

    # Fill the new histogram by mapping content from the original
    for i_bin in range(1, n_bins + 1):
        original_content = original_hist.GetBinContent(i_bin)
        original_error = original_hist.GetBinError(i_bin)
        new_hist.SetBinContent(i_bin, original_content)
        new_hist.SetBinError(i_bin, original_error)

    new_hist.Scale(1.0 / 3.0)

    return new_hist


v2_stat_scaled = []
scaled_dir = output_file.mkdir("scaled")

for i_cent, cent_limits in enumerate(new_centrality_classes):
    new_hist = rescale_bins(new_v2_stat[i_cent])
    utils.setHistStyle(new_hist, new_cent_colours[i_cent])
    v2_stat_scaled.append(new_hist)
    scaled_dir.cd()
    v2_stat_scaled[i_cent].Write()


cV2_proton = ROOT.TCanvas("cV2_proton", "cV2_proton", 800, 600)
info_panel_proton = ROOT.TPaveText(0.16, 0.71, 0.36, 0.87, "NDC")
info_panel_proton.SetBorderSize(0)
info_panel_proton.SetFillStyle(0)
info_panel_proton.SetTextAlign(12)
info_panel_proton.SetTextFont(42)
info_panel_proton.SetTextSize(0.04)
info_panel_proton.AddText(r"ALICE")
info_panel_proton.AddText(r"Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.36 TeV")
frame_proton = cV2_proton.DrawFrame(
    0.0, -0.042, 5.0, 0.39, r";#it{p}_{T}/A (GeV/#it{c}); #it{v}_{2} / A"
)
cV2_proton.SetBottomMargin(0.13)
cV2_proton.SetLeftMargin(0.13)
cV2_proton.SetBorderSize(0)
cV2_proton.cd()

legend_proton = ROOT.TLegend(
    0.59, 0.29, 0.93, 0.34, "p, JHEP 09 (2018) 006, 2018", "brNDC"
)
legend_proton.SetBorderSize(0)
legend_proton.SetNColumns(2)

legend_proton.AddEntry(hProton_1020_stat, r"10-20 %", "PE")
hProton_1020_stat.Draw("PE SAME")
hProton_1020_syst.Draw("PE2 SAME")

# legend_proton.AddEntry(hProton_2030_stat, r"20-30 %", "PE")
# hProton_2030_stat.Draw("PE SAME")
# hProton_2030_syst.Draw("PE2 SAME")

# legend_proton.AddEntry(hProton_3040_stat, r"30-40 %", "PE")
# hProton_3040_stat.Draw("PE SAME")
# hProton_3040_syst.Draw("PE2 SAME")

legend_proton.AddEntry(hProton_4050_stat, r"40-50 %", "PE")
hProton_4050_stat.Draw("PE SAME")
hProton_4050_syst.Draw("PE2 SAME")


legend_helium_scaled = ROOT.TLegend(0.60, 0.20, 0.93, 0.27, r"{}^{3}He", "brNDC")
legend_helium_scaled.SetBorderSize(0)
legend_helium_scaled.SetNColumns(2)

legend_helium_scaled.AddEntry(v2_stat_scaled[1], r"10-20 %", "PF")
v2_stat_scaled[1].Draw("PE SAME")

# legend_helium_scaled.AddEntry(v2_stat_scaled[2], r"20-30 %", "PF")
# v2_stat_scaled[2].Draw("PE SAME")

# legend_helium_scaled.AddEntry(v2_stat_scaled[3], r"30-40 %", "PF")
# v2_stat_scaled[3].Draw("PE SAME")

legend_helium_scaled.AddEntry(v2_stat_scaled[4], r"40-50 %", "PF")
v2_stat_scaled[4].Draw("PE SAME")

legend_proton.Draw()
legend_helium_scaled.Draw()
info_panel_proton.Draw()

output_file.cd()
cV2_proton.Write()
cV2_proton.SaveAs(f"../results_pass4_new/plots/{cV2_proton.GetName()}.pdf")
