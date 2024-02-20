import ROOT

import sys
sys.path.append('utils')
import utils as utils


input_file_name = './run_task/AnalysisResults_small.root'
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

v2_ft0c_val = []
v2_ft0c_err = []
v2_ft0a_val = []
v2_ft0a_err = []
v2_fv0a_val = []
v2_fv0a_err = []

for i_pt in range(1, n_pt_bins + 1):

    pt_dir = output_file.mkdir(f'pt_{i_pt}')
    pt_dir.cd()

    # system infopanel
    info_panel_ft0a = ROOT.TPaveText(0.6, 0.6, 0.8, 0.82, 'NDC')
    info_panel_ft0a.SetBorderSize(0)
    info_panel_ft0a.SetFillStyle(0)
    info_panel_ft0a.SetTextAlign(12)
    info_panel_ft0a.SetTextFont(42)
    info_panel_ft0a.AddText('ALICE')
    info_panel_ft0a.AddText(r'PbPb, #sqrt{#it{s}_{nn}} = 5.36 TeV')
    info_panel_ft0a.AddText(f'{cent_limits[0]} - {cent_limits[1]} % {cent_detector_label}')
    pt_label = f'{pt_axis.GetBinLowEdge(i_pt):.1f}' + r' #leq #it{p}_{T} < ' + \
        f'{pt_axis.GetBinUpEdge(i_pt):.1
           f}' + r' GeV/#it{c}'
    info_panel_ft0a.AddText(pt_label)

    # FT0C
    hSpVsNsigma_FT0C, hNsigma_FT0C = utils.getHistos1D(
        hSpFT0CvsNsigmaHe3VsPtvsCent, i_pt, i_pt, 30, 50, f'hNsigma_{i_pt}_FT0C', f'hSpVsNsigma_{i_pt}_FT0C', n_rebin=4)
    utils.setHistStyle(hNsigma_FT0C, ROOT.kRed+1, linewidth=2)
    utils.setHistStyle(hSpVsNsigma_FT0C, ROOT.kAzure+1, linewidth=2)
    hSpVsNsigma_FT0C.GetYaxis().SetRangeUser(-2., 2.)

    # fit with a pol0
    fit_FT0C = ROOT.TF1('fit_FT0C', 'pol0', -1, 1)
    hSpVsNsigma_FT0C.Fit(fit_FT0C, 'R')

    val_ft0c = fit_FT0C.GetParameter(0)
    v2_ft0c_val.append(val_ft0c)
    err_ft0c = fit_FT0C.GetParError(0)
    v2_ft0c_err.append(err_ft0c)

    # fit panel ftoc
    fit_panel_ft0a = ROOT.TPaveText(0.6, 0.6, 0.8, 0.82, 'NDC')
    fit_panel_ft0a.SetBorderSize(0)
    fit_panel_ft0a.SetFillStyle(0)
    fit_panel_ft0a.SetTextAlign(12)
    fit_panel_ft0a.SetTextFont(42)
    fit_panel_ft0a.AddText(f'p_{{0}} = {val_ft0c:.3f} #pm {err_ft0c:.3f}')

    canvas_FT0C = utils.geCanvasWithTwoPanels(
        f'cSpVsNsigma_{i_pt}_FT0C', hSpVsNsigma_FT0C, hNsigma_FT0C, top_panel=fit_panel_ft0a, bottom_panel=info_panel_ft0a)

    hSpVsNsigma_FT0C.Write()
    hNsigma_FT0C.Write()
    canvas_FT0C.Write()

    # FT0A
    hSpVsNsigma_FT0A, hNsigma_FT0A = utils.getHistos1D(
        hSpFT0AvsNsigmaHe3VsPtvsCent, i_pt, i_pt, 30, 50, f'hNsigma_{i_pt}_FT0A', f'hSpVsNsigma_{i_pt}_FT0A', n_rebin=4)
    utils.setHistStyle(hNsigma_FT0A, ROOT.kRed+1, linewidth=2)
    utils.setHistStyle(hSpVsNsigma_FT0A, ROOT.kAzure+1, linewidth=2)
    hSpVsNsigma_FT0A.GetYaxis().SetRangeUser(-2., 2.)

    info_panel_ft0a = info_panel_ft0a.Clone('panel_ft0a')

    # fit with a pol0
    fit_FT0A = ROOT.TF1('fit_FT0A', 'pol0', -1, 1)
    hSpVsNsigma_FT0A.Fit(fit_FT0A, 'R')

    val_ft0a = fit_FT0A.GetParameter(0)
    v2_ft0a_val.append(val_ft0a)
    err_ft0a = fit_FT0A.GetParError(0)
    v2_ft0a_err.append(err_ft0a)

    # fit panel ftoa
    fit_panel_ft0a = ROOT.TPaveText(0.6, 0.6, 0.8, 0.82, 'NDC')
    fit_panel_ft0a.SetBorderSize(0)
    fit_panel_ft0a.SetFillStyle(0)
    fit_panel_ft0a.SetTextAlign(12)
    fit_panel_ft0a.SetTextFont(42)
    fit_panel_ft0a.AddText(f'p_{{0}} = {val_ft0a:.3f} #pm {err_ft0a:.3f}')

    canvas_FT0A = utils.geCanvasWithTwoPanels(
        f'cSpVsNsigma_{i_pt}_FT0A', hSpVsNsigma_FT0A, hNsigma_FT0A, top_panel=fit_panel_ft0a, bottom_panel=info_panel_ft0a)

    hSpVsNsigma_FT0A.Write()
    hNsigma_FT0A.Write()
    canvas_FT0A.Write()

    # FV0A
    hSpVsNsigma_FV0A, hNsigma_FV0A = utils.getHistos1D(
        hSpFV0AvsNsigmaHe3VsPtvsCent, i_pt, i_pt, 30, 50, f'hNsigma_{i_pt}_FV0A', f'hSpVsNsigma_{i_pt}_FV0A', n_rebin=4)
    utils.setHistStyle(hNsigma_FV0A, ROOT.kRed+1, linewidth=2)
    utils.setHistStyle(hSpVsNsigma_FV0A, ROOT.kAzure+1, linewidth=2)
    hSpVsNsigma_FV0A.GetYaxis().SetRangeUser(-2., 2.)

    info_panel_fv0a = info_panel_ft0a.Clone('panel_fv0a')

     # fit with a pol0
    fit_FV0A = ROOT.TF1('fit_FV0A', 'pol0', -1, 1)
    hSpVsNsigma_FV0A.Fit(fit_FV0A, 'R')

    val_fv0a = fit_FV0A.GetParameter(0)
    v2_fv0a_val.append(val_fv0a)
    err_fv0a = fit_FV0A.GetParError(0)
    v2_fv0a_err.append(err_fv0a)

    # fit panel fv0a
    fit_panel_fv0a = ROOT.TPaveText(0.6, 0.6, 0.8, 0.82, 'NDC')
    fit_panel_fv0a.SetBorderSize(0)
    fit_panel_fv0a.SetFillStyle(0)
    fit_panel_fv0a.SetTextAlign(12)
    fit_panel_fv0a.SetTextFont(42)
    fit_panel_fv0a.AddText(f'p_{{0}} = {val_fv0a:.3f} #pm {err_fv0a:.3f}')

    canvas_FV0A = utils.geCanvasWithTwoPanels(
        f'cSpVsNsigma_{i_pt}_FV0A', hSpVsNsigma_FV0A, hNsigma_FV0A, top_panel=fit_panel_fv0a, bottom_panel=info_panel_fv0a)

    hSpVsNsigma_FV0A.Write()
    hNsigma_FV0A.Write()
    canvas_FV0A.Write()

