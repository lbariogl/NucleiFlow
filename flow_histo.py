import ROOT

import sys
sys.path.append('utils')
import utils as utils


input_file_name = './run_task/AnalysisResults.root'
input_dir_name = 'nucleiFlow'
output_file_name = 'output.root'

cent_detector_label = 'FT0C'
cent_limits = [30, 50]

input_file = ROOT.TFile(input_file_name)
output_file = ROOT.TFile(output_file_name, 'RECREATE')

hCentFV0A = input_file.Get(f'{input_dir_name}/hCentFV0A')
hCentFV0A.SetDirectory(0)
hCentFT0A = input_file.Get(f'{input_dir_name}/hCentFT0A')
hCentFT0A.SetDirectory(0)
hCentFT0C = input_file.Get(f'{input_dir_name}/hCentFT0C')
hCentFT0C.SetDirectory(0)
hCentFT0M = input_file.Get(f'{input_dir_name}/hCentFT0M')
hCentFT0M.SetDirectory(0)

hSpFT0AvsNsigmaHe3VsPtvsCent = input_file.Get(
    f'{input_dir_name}/hSpFT0AvsNsigmaHe3VsPtvsCent')
hSpFT0CvsNsigmaHe3VsPtvsCent = input_file.Get(
    f'{input_dir_name}/hSpFT0CvsNsigmaHe3VsPtvsCent')
hSpFV0AvsNsigmaHe3VsPtvsCent = input_file.Get(
    f'{input_dir_name}/hSpFV0AvsNsigmaHe3VsPtvsCent')

pt_axis = hSpFT0AvsNsigmaHe3VsPtvsCent.GetAxis(2)
n_pt_bins = pt_axis.GetNbins()
pt_min = pt_axis.GetBinLowEdge(1)
pt_max = pt_axis.GetBinUpEdge(n_pt_bins)

hSPvsPT_FT0C = ROOT.TH1F('hSPvsPT_FT0C', r';#it{p}_{T} (GeV/#it{c}) ;\langle #hat{u}_{2} #upoint #vec{Q}_{2}^{FT0C} \rangle', n_pt_bins, pt_min, pt_max)
utils.setHistStyle(hSPvsPT_FT0C, ROOT.kAzure+2)
hSPvsPT_FT0A = ROOT.TH1F('hSPvsPT_FT0A', r';#it{p}_{T} (GeV/#it{c}) ;\langle #hat{u}_{2} #upoint #vec{Q}_{2}^{FT0A} \rangle', n_pt_bins, pt_min, pt_max)
utils.setHistStyle(hSPvsPT_FT0A, ROOT.kRed+1)
hSPvsPT_FV0A = ROOT.TH1F('hSPvsPT_FV0A', r';#it{p}_{T} (GeV/#it{c}) ;\langle #hat{u}_{2} #upoint #vec{Q}_{2}^{FV0A} \rangle', n_pt_bins, pt_min, pt_max)
utils.setHistStyle(hSPvsPT_FV0A, ROOT.kGreen+3)


for i_pt in range(1, n_pt_bins + 1):

    pt_dir = output_file.mkdir(f'pt_{i_pt}')
    pt_dir.cd()

    cent_low = cent_limits[0]
    cent_up = cent_limits[1]

    pt_low = pt_axis.GetBinLowEdge(i_pt)
    pt_up = pt_axis.GetBinUpEdge(i_pt)

    utils.getCompleteCanvas(hSpFT0CvsNsigmaHe3VsPtvsCent, cent_low, cent_up, i_pt, i_pt,
                            pt_dir, hSPvsPT_FT0C, qvec_detector_label='FT0C', out_pt_bin=i_pt)

    utils.getCompleteCanvas(hSpFT0AvsNsigmaHe3VsPtvsCent, cent_low, cent_up, i_pt, i_pt,
                            pt_dir, hSPvsPT_FT0A, qvec_detector_label='FT0A', out_pt_bin=i_pt)

    utils.getCompleteCanvas(hSpFV0AvsNsigmaHe3VsPtvsCent, cent_low, cent_up, i_pt, i_pt,
                            pt_dir, hSPvsPT_FV0A, qvec_detector_label='FV0A', out_pt_bin=i_pt)

output_file.cd()
hSPvsPT_FT0C.Write()
hSPvsPT_FT0A.Write()
hSPvsPT_FV0A.Write()

cSPvsPT= ROOT.TCanvas('cSPvsPT', 'cSPvsPT', 800, 600)
cSPvsPT.SetLeftMargin(0.13)
cSPvsPT.SetBottomMargin(0.13)
frame = cSPvsPT.DrawFrame(0., -2., 6., 2., r';#it{p}_{T} (GeV/#it{c}) ;\langle \vec{u}_{2} \cdot \vec{Q} \rangle')
hSPvsPT_FT0C.Draw('PE SAME')
hSPvsPT_FT0A.Draw('PE SAME')
hSPvsPT_FV0A.Draw('PE SAME')
frame.GetYaxis().SetTitleOffset(1.1)
frame.GetXaxis().SetTitleOffset(1.1)
legend = ROOT.TLegend(0.69, 0.7, 0.88, 0.86, '', 'BRNDC')
legend.AddEntry(hSPvsPT_FT0C, 'FT0C', 'PE')
legend.AddEntry(hSPvsPT_FT0A, 'FT0A', 'PE')
legend.AddEntry(hSPvsPT_FV0A, 'FV0A', 'PE')
legend.Draw()
cSPvsPT.Write()
