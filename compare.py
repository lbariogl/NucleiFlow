import ROOT

old_centrality_classes = [[0, 10], [10, 20], [20, 40], [40, 60]]

old_v2_stat = []
old_v2_syst = []

old_file = ROOT.TFile("../results_pass3_new/final.root")
old_file.ls()

for cent_limits in old_centrality_classes:
    print(f"cent_limits: [{cent_limits[0]}, {cent_limits[1]}]")
    cent_dir = old_file.Get(f"cent_{cent_limits[0]}_{cent_limits[1]}")
    old_v2_stat.append(cent_dir.Get(f"hV2vsPt_cent_{cent_limits[0]}_{cent_limits[1]}"))
    old_v2_syst.append(
        cent_dir.Get(f"hV2vsPt_cent_{cent_limits[0]}_{cent_limits[1]}_syst")
    )

print(f"{old_v2_stat=}")
print(f"{old_v2_syst=}")

new_centrality_classes = [
    [0, 10],
    [10, 20],
    [20, 30],
    [30, 40],
    [40, 50],
    [50, 60],
    [60, 80],
]

new_v2_stat = []
new_file = ROOT.TFile("../results_pass4/flow.root")

for cent_limits in new_centrality_classes:
    print(f"{cent_limits[0]} - {cent_limits[1]}")
    cent_dir = new_file.Get(f"cent_{cent_limits[0]}_{cent_limits[1]}")
    hV2_tmp = cent_dir.Get(f"default/hV2vsPt_cent_{cent_limits[0]}_{cent_limits[1]}")
    hV2_tmp.SetMarkerStyle(25)
    new_v2_stat.append(hV2_tmp)

output_file = ROOT.TFile("../results_pass4/comp_pass3_pass4.root", "RECREATE")

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

new_legend = ROOT.TLegend(0.16, 0.56, 0.32, 0.70, "pass4", "brNDC")
new_legend.SetBorderSize(0)
new_legend.SetNColumns(2)

for i_cent, cent_limits in enumerate(new_centrality_classes):
    cent_label = f"{cent_limits[0]}-{cent_limits[1]}" + r"%"
    new_legend.AddEntry(new_v2_stat[i_cent], cent_label, "P")
    new_v2_stat[i_cent].Draw("PE SAME")

old_legend = ROOT.TLegend(0.16, 0.44, 0.32, 0.54, "pass3", "brNDC")
old_legend.SetBorderSize(0)
old_legend.SetNColumns(2)

for i_cent, cent_limits in enumerate(old_centrality_classes):
    cent_label = f"{cent_limits[0]}-{cent_limits[1]}" + r"%"
    old_legend.AddEntry(old_v2_stat[i_cent], cent_label, "PF")
    old_v2_stat[i_cent].Draw("PEX0 SAME")
    old_v2_syst[i_cent].Draw("PE2 SAME")

old_legend.Draw()
new_legend.Draw()
info_panel_total.Draw()

output_file.cd()
cV2.Write()
cV2.SaveAs(f"../results_pass4/plots/{cV2.GetName()}.pdf")
