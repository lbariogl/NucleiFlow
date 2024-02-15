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
        hNsigma.Rebin(n_rebin)
        hSPvsNsigma_tmp.Rebin(n_rebin)
    hSPvsNsigma = getAverage2D(hSPvsNsigma_tmp, histo2_name)
    hSPvsNsigma.SetTitle('')
    return hSPvsNsigma, hNsigma


def geCanvasWithTwoPanels(canvas_name, histo_1, histo_2, info_panel=None):
    canvas = ROOT.TCanvas(canvas_name, canvas_name, 800, 600)
    pad_top = ROOT.TPad("pad_top", "pad_top", 0.0, 0.3, 1.0, 1.0, 0)
    pad_top.SetLeftMargin(0.15)
    pad_top.SetBottomMargin(0.)
    pad_top.Draw()
    pad_top.cd()
    histo_1.GetYaxis().SetTitleSize(0.07)
    histo_1.Draw('PE')
    canvas.cd()
    pad_bottom = ROOT.TPad("pad_bottom", "pad_bottom", 0.0, 0.0, 1.0, 0.3, 0)
    pad_bottom.SetLeftMargin(0.15)
    pad_bottom.SetTopMargin(0.)
    pad_bottom.SetBottomMargin(0.3)
    pad_bottom.Draw()
    pad_bottom.cd()
    histo_2.GetYaxis().SetTitleSize(0.07)
    histo_2.GetYaxis().SetLabelSize(0.05)
    histo_2.Draw('PE')
    if info_panel:
        histo_2.GetListOfFunctions().Add(info_panel)
    return canvas
