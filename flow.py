import ROOT

import sys
sys.path.append('utils')
import utils as utils


input_file_name = './run_task/AnalysisResults.root'
input_dir_name = 'nucleiFlow'
output_file_name = 'output.root'

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

for i_pt in range(1, n_pt_bins + 1):
    pt_dir = output_file.mkdir(f'pt_{i_pt}')
    pt_dir.cd()

    hSpVsNsigma, hNsigma = utils.getHistos1D(
        hSpFT0CvsNsigmaHe3VsPtvsCent, i_pt, i_pt, 30, 50, f'hNsigma_{i_pt}', f'hSpVsNsigma_{i_pt}')
    utils.setHistStyle(hNsigma, ROOT.kRed+1, linewidth=2)
    utils.setHistStyle(hSpVsNsigma, ROOT.kAzure+1, linewidth=2)

    info_panel = ROOT.TPaveText(0.7, 0.6, 0.9, 0.82, 'NDC')
    info_panel.SetBorderSize(0)
    info_panel.SetFillStyle(0)
    info_panel.SetTextAlign(12)
    info_panel.SetTextFont(42)
    info_panel.AddText('ALICE')
    info_panel.AddText(r'PbPb, #sqrt{#it{s}_{nn}} = 5.36 TeV')
    info_panel.AddText(r'30 - 50 % FTOC')
    pt_label = f'{pt_axis.GetBinLowEdge(i_pt)}' + r' #leq #it{p}_{T} < ' + \
        f'{pt_axis.GetBinUpEdge(i_pt)}' + r' GeV/#it{c}'
    info_panel.AddText(pt_label)

    canvas = utils.geCanvasWithTwoPanels(
        f'cSpVsNsigma_{i_pt}', hSpVsNsigma, hNsigma, info_panel)

    hSpVsNsigma.Write()
    hNsigma.Write()
    canvas.Write()
