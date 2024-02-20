import ROOT

ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
ROOT.gStyle.SetOptStat(0)


def setHistStyle(hist, color, marker=20, fillstyle=0, linewidth=1):
    hist.SetMarkerColor(color)
    hist.SetLineColor(color)
    hist.SetMarkerStyle(marker)
    hist.SetFillStyle(fillstyle)
    hist.SetLineWidth(linewidth)


def getAverage2D(histo2D, histo_name="histo"):
    x_axis = histo2D.GetXaxis()
    nX_bins = x_axis.GetNbins()
    histo2D_y_title = histo2D.GetYaxis().GetTitle()
    new_title = r'\langle ' + histo2D_y_title + r' \rangle'
    average_histo = histo2D.ProjectionX(histo_name)
    average_histo.Reset()
    average_histo.GetYaxis().SetTitle(new_title)
    for i_bin in range(1, nX_bins + 1):
        histo_tmp = histo2D.ProjectionY('histo_tmp', i_bin, i_bin)
        average_histo.SetBinContent(i_bin, histo_tmp.GetMean())
        average_histo.SetBinContent(i_bin, histo_tmp.GetRMS())
        del histo_tmp
    return average_histo


def getHistos1D(thn_sparse, pt_bin_low, pt_bin_up, cent_low, cent_up, histo1_name='histo1', histo2_name='histo2', n_rebin=1):
    pt_axis = thn_sparse.GetAxis(2)
    pt_axis.SetRange(pt_bin_low, pt_bin_up)
    cent_axis = thn_sparse.GetAxis(3)
    cent_axis.SetRange(cent_axis.FindBin(cent_low),
                       cent_axis.FindBin(cent_up)-1)
    hNsigma = thn_sparse.Projection(1)
    hNsigma.SetName(histo1_name)
    hNsigma.SetTitle('')
    hSPvsNsigma_tmp = thn_sparse.Projection(0, 1)
    if n_rebin > 1:
        hNsigma.RebinX(n_rebin)
        hSPvsNsigma_tmp.RebinX(n_rebin)
    hSPvsNsigma = getAverage2D(hSPvsNsigma_tmp, histo2_name)
    hSPvsNsigma.SetTitle('')
    return hSPvsNsigma, hNsigma


def geCanvasWithTwoPanels(canvas_name, histo_1, histo_2, top_panel=None, bottom_panel=None):
    canvas = ROOT.TCanvas(canvas_name, canvas_name, 800, 600)
    pad_top = ROOT.TPad("pad_top", "pad_top", 0.0, 0.5, 1.0, 1.0, 0)
    pad_top.SetLeftMargin(0.15)
    pad_top.SetBottomMargin(0.)
    pad_top.Draw()
    pad_top.cd()
    histo_1.GetYaxis().SetTitleSize(0.07)
    histo_1.GetYaxis().SetLabelSize(0.05)
    histo_1.Draw('PE')
    if top_panel:
        histo_1.GetListOfFunctions().Add(top_panel)
    canvas.cd()
    pad_bottom = ROOT.TPad("pad_bottom", "pad_bottom", 0.0, 0.0, 1.0, 0.5, 0)
    pad_bottom.SetLeftMargin(0.15)
    pad_bottom.SetTopMargin(0.)
    pad_bottom.SetBottomMargin(0.3)
    pad_bottom.Draw()
    pad_bottom.cd()
    histo_2.GetYaxis().SetTitleSize(0.07)
    histo_2.GetYaxis().SetLabelSize(0.05)
    histo_2.GetXaxis().SetTitleSize(0.07)
    histo_2.GetXaxis().SetLabelSize(0.05)
    histo_2.Draw('PE')
    if bottom_panel:
        histo_2.GetListOfFunctions().Add(bottom_panel)
    return canvas


def getCompleteCanvas(hSpvsNsigmaVsPtvsCent, cent_low, cent_up, pt_bin_low, pt_bin_up, output_dir, arr_fit_val, arr_fit_err, qvec_detector_label='FT0C', cent_detector_label='FTOC', suffix = 0):

    output_dir.cd()

    # get axis info
    pt_axis = hSpvsNsigmaVsPtvsCent.GetAxis(2)
    pt_low = pt_axis.GetBinLowEdge(pt_bin_low)
    pt_up = pt_axis.GetBinUpEdge(pt_bin_up)

    cent_axis = hSpvsNsigmaVsPtvsCent.GetAxis(3)
    cent_bin_low = cent_axis.FindBin(cent_low)
    cent_bin_up = cent_axis.FindBin(cent_up) - 1

    # system infopanel
    info_panel = ROOT.TPaveText(0.6, 0.6, 0.8, 0.82, 'NDC')
    info_panel.SetBorderSize(0)
    info_panel.SetFillStyle(0)
    info_panel.SetTextAlign(12)
    info_panel.SetTextFont(42)
    info_panel.AddText('ALICE')
    info_panel.AddText(r'PbPb, #sqrt{#it{s}_{nn}} = 5.36 TeV')
    info_panel.AddText(f'{cent_low} - {cent_up} % {cent_detector_label}')
    pt_label = f'{pt_low:.1f}' + r' #leq #it{p}_{T} < ' + \
        f'{pt_up:.1f}' + r' GeV/#it{c}'
    info_panel.AddText(pt_label)

    # create 1D histograms
    hSpVsNsigma, hNsigma = getHistos1D(
        hSpvsNsigmaVsPtvsCent, pt_bin_low, pt_bin_up, 30, 50, f'hNsigma_{suffix}_{qvec_detector_label}', f'hSpVsNsigma_{suffix}_{qvec_detector_label}', n_rebin=4)
    setHistStyle(hNsigma, ROOT.kRed+1, linewidth=2)
    setHistStyle(hSpVsNsigma, ROOT.kAzure+1, linewidth=2)
    hSpVsNsigma.GetYaxis().SetRangeUser(-2., 2.)

    # fit with a pol0
    fit = ROOT.TF1('fit', 'pol0', -1, 1)
    hSpVsNsigma.Fit(fit, 'R')

    val = fit.GetParameter(0)
    arr_fit_val.append(val)
    err = fit.GetParError(0)
    arr_fit_err.append(err)

    fit_panel = ROOT.TPaveText(0.6, 0.6, 0.8, 0.82, 'NDC')
    fit_panel.SetBorderSize(0)
    fit_panel.SetFillStyle(0)
    fit_panel.SetTextAlign(12)
    fit_panel.SetTextFont(42)
    fit_panel.AddText(f'p_{{0}} = {val:.3f} #pm {err:.3f}')

    canvas = geCanvasWithTwoPanels(
        f'cSpVsNsigma_{suffix}_{qvec_detector_label}', hSpVsNsigma, hNsigma, top_panel=fit_panel, bottom_panel=info_panel)

    hSpVsNsigma.Write()
    hNsigma.Write()
    canvas.Write()
