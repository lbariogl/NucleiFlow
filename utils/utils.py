import ROOT
from hipe4ml.tree_handler import TreeHandler
import numpy as np

ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
ROOT.gStyle.SetOptStat(0)


def setHistStyle(hist, color, marker=20, fillstyle=0, linewidth=1):
    hist.SetMarkerColor(color)
    hist.SetLineColor(color)
    hist.SetMarkerStyle(marker)
    hist.SetFillStyle(fillstyle)
    hist.SetLineWidth(linewidth)

# get average ITS cluster size from its size bitmap
def getITSClSize(itsSizeBitmap):
    sum = 0
    nClus = 0
    for i_layer in range(0, 7):
        size = (itsSizeBitmap >> i_layer*4) & 15
        if size > 0:
            nClus = nClus+1
            sum += size
    return sum / nClus

# vectorised version of getITSClSize
getITSClSize_vectorised = np.vectorize(getITSClSize)

# get the sign of the track from the flag
def getSign(flags):
    if (flags & 256):
        return 1
    else:
        return -1

# vectorised version of getSign
getSign_vectorised = np.vectorize(getSign)

# get
def trackedAsHe(flags):
    pid_flag = (flags >> 12) & 15
    if (pid_flag == 7 or pid_flag == 8):
        return True
    else:
        return False

# vectorised version of getSign
trackedAsHe_vectorised = np.vectorize(trackedAsHe)

# Bethe-Bloch-Aleph function as defined in O2
def BBA(x, p):

    bg = x/2.80839160743
    # (TMath::Abs(x) / TMath::Sqrt(1 + x*x))
    beta = abs(bg) / np.sqrt(1 + bg * bg)

    # TMath::Power((TMath::Abs(x) / TMath::Sqrt(1 + x*x)),[3])
    aa = np.power(beta, p[3])
    bb = np.power(1 / abs(bg), p[4])  # TMath::Power(1/TMath::Abs(x),[4])
    # TMath::Log([2] + TMath::Power(1/TMath::Abs(x),[4]))
    bb = np.log(p[2] + bb)

    # '([1] - TMath::Power((TMath::Abs(x) / TMath::Sqrt(1 + x*x)),[3]) - TMath::Log([2] + TMath::Power(1/TMath::Abs(x),[4]))) * [0] / TMath::Power((TMath::Abs(x) / TMath::Sqrt(1 + x*x)),[3])'
    return (p[1] - aa - bb) * p[0] / aa

# Literal form of BBA, to be used inside TF1
func_string = '([1] - TMath::Power((TMath::Abs(2*x/2.80839160743) / TMath::Sqrt(1 + 2*2*x*x/2.80839160743/2.80839160743)),[3]) - TMath::Log([2] + TMath::Power(1/TMath::Abs(2*x/2.80839160743),[4]))) * [0] / TMath::Power((TMath::Abs(2*x/2.80839160743) / TMath::Sqrt(1 + 2*2*x*x/2.80839160743/2.80839160743)),[3])'

default_bb_parameters = p_train = [-321.34, 0.6539, 1.591, 0.8225, 2.363]

# N-sigma TPC at specific rigidity
def getNsigmaTPC(x, tpc_signal, parameters=default_bb_parameters, resolution_perc=0.09):
    exp_signal = BBA(x, parameters)
    resolution = exp_signal * resolution_perc
    return (tpc_signal - exp_signal) / resolution

# vectorised version of N-sigma TPC at specific rigidity
getNsigmaTPC_vectorised = np.vectorize(getNsigmaTPC)

# redifine columns in the complete data-frame
def redifineColumns(complete_df):
    complete_df.eval(
        'fAvgItsClusSize = @getITSClSize_vectorised(fITSclusterSizes)', inplace=True)
    complete_df.drop(columns=['fITSclusterSizes'])
    complete_df.eval('fSign = @getSign_vectorised(fFlags)', inplace=True)
    complete_df.eval(
        'fTrackedAsHe = @trackedAsHe_vectorised(fFlags)', inplace=True)
    complete_df.loc[complete_df['fTrackedAsHe'] == True, 'fTPCInnerParam'] = complete_df['fTPCInnerParam']/2
    complete_df.eval(
        'fNsigmaTPC3He = @getNsigmaTPC_vectorised(2*fTPCInnerParam, fTPCsignal)', inplace=True)
    # ScalarProducts
    complete_df.eval('fSpFT0C = fPt*cos(fPhi) * fXQvecFT0C + fPt*sin(fPhi) * fYQvecFT0C', inplace=True)
    complete_df.eval('fSpFT0A = fPt*cos(fPhi) * fXQvecFT0A + fPt*sin(fPhi) * fYQvecFT0A', inplace=True)
    complete_df.eval('fSpFV0A = fPt*cos(fPhi) * fXQvecFV0A + fPt*sin(fPhi) * fYQvecFV0A', inplace=True)
    complete_df.eval('fSpTPCpos = fPt*cos(fPhi) * fXQvecTPCpos + fPt*sin(fPhi) * fYQvecTPCpos', inplace=True)
    complete_df.eval('fSpTPCneg = fPt*cos(fPhi) * fXQvecTPCneg + fPt*sin(fPhi) * fYQvecTPCneg', inplace=True)

def getBBAfunctions(parameters, resolution, n_sigma=5):
    upper_scale = 1 + resolution * n_sigma
    lower_scale = 1 - resolution * n_sigma
    func_string_up = f'{upper_scale} * ' + func_string
    func_string_down = f'{lower_scale} * ' + func_string

    func_BB_left = ROOT.TF1('func_BB_left', func_string, -6, -0.5, 5)
    func_BB_left.SetParameters(parameters[0],parameters[1],parameters[2], parameters[3], parameters[4])
    func_BB_left.SetLineColor(ROOT.kRed)

    func_BB_left_up = ROOT.TF1('func_BB_left_up', func_string_up, -6, -0.5, 5)
    func_BB_left_up.SetParameters(parameters[0],parameters[1],parameters[2], parameters[3], parameters[4])
    func_BB_left_up.SetLineColor(ROOT.kRed)
    func_BB_left_up.SetLineStyle(ROOT.kDashed)

    func_BB_left_down = ROOT.TF1('func_BB_left_down', func_string_down, -6, -0.5, 5)
    func_BB_left_down.SetParameters(parameters[0],parameters[1],parameters[2], parameters[3], parameters[4])
    func_BB_left_down.SetLineColor(ROOT.kRed)
    func_BB_left_down.SetLineStyle(ROOT.kDashed)

    func_BB_right = ROOT.TF1('func_BB_right', func_string, 0.5, 6., 5)
    func_BB_right.SetParameters(parameters[0],parameters[1],parameters[2], parameters[3], parameters[4])
    func_BB_right.SetLineColor(ROOT.kRed)

    func_BB_right_up = ROOT.TF1('func_BB_right_up', func_string_up, 0.5, 6., 5)
    func_BB_right_up.SetParameters(parameters[0],parameters[1],parameters[2], parameters[3], parameters[4])
    func_BB_right_up.SetLineColor(ROOT.kRed)
    func_BB_right_up.SetLineStyle(ROOT.kDashed)

    func_BB_right_down = ROOT.TF1('func_BB_right_down', func_string_down, 0.5, 6., 5)
    func_BB_right_down.SetParameters(parameters[0],parameters[1],parameters[2], parameters[3], parameters[4])
    func_BB_right_down.SetLineColor(ROOT.kRed)
    func_BB_right_down.SetLineStyle(ROOT.kDashed)

    functions = [func_BB_left, func_BB_left_up, func_BB_left_down, func_BB_right, func_BB_right_up, func_BB_right_down]

    return functions


def getAverage2D(histo2D, histo_name="histo"):
    x_axis = histo2D.GetXaxis()
    nX_bins = x_axis.GetNbins()
    histo2D_y_title = histo2D.GetYaxis().GetTitle()
    new_title = r'#LT ' + histo2D_y_title + r' #GT'
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
    hSPvsNsigma_tmp.SetName(f'{histo2_name}_notAveraged')
    if n_rebin > 1:
        hNsigma.RebinX(n_rebin)
        hSPvsNsigma_tmp.RebinX(n_rebin)
    hSPvsNsigma = getAverage2D(hSPvsNsigma_tmp, histo2_name)
    hSPvsNsigma.SetTitle('')
    return hSPvsNsigma, hNsigma, hSPvsNsigma_tmp


def getCanvasWithTwoPanels(canvas_name, histo_1, histo_2, top_panel=None, bottom_panel=None):
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


def getCompleteCanvas(hSpvsNsigmaVsPtvsCent, cent_low, cent_up, pt_bin_low, pt_bin_up, output_dir, hSPvsPt, qvec_detector_label='FT0C', cent_detector_label='FTOC', out_pt_bin=0):

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
    hSpVsNsigma, hNsigma, hSpVsNsigma_notAveraged = getHistos1D(
        hSpvsNsigmaVsPtvsCent, pt_bin_low, pt_bin_up, cent_bin_low, cent_bin_up, f'hNsigma_{out_pt_bin}_{qvec_detector_label}', f'hSpVsNsigma_{out_pt_bin}_{qvec_detector_label}', n_rebin=4)
    setHistStyle(hNsigma, ROOT.kRed+1, linewidth=2)
    setHistStyle(hSpVsNsigma, ROOT.kAzure+1, linewidth=2)
    hSpVsNsigma.GetYaxis().SetRangeUser(-2., 2.)

    # # fit with a pol0
    # fit = ROOT.TF1('fit', 'pol0', -1, 1)
    # hSpVsNsigma.Fit(fit, 'R')

    # val = fit.GetParameter(0)
    # err = fit.GetParError(0)

    # hSPvsPt.SetBinContent(out_pt_bin, val)
    # hSPvsPt.SetBinError(out_pt_bin, err)

    # fit_panel = ROOT.TPaveText(0.6, 0.6, 0.8, 0.82, 'NDC')
    # fit_panel.SetBorderSize(0)
    # fit_panel.SetFillStyle(0)
    # fit_panel.SetTextAlign(12)
    # fit_panel.SetTextFont(42)
    # fit_panel.AddText(f'p_{{0}} = {val:.3f} #pm {err:.3f}')

    canvas = getCanvasWithTwoPanels(
        f'cSpVsNsigma_{out_pt_bin}_{qvec_detector_label}', hSpVsNsigma, hNsigma, bottom_panel=info_panel)

    hSpVsNsigma.Write()
    hNsigma.Write()
    hSpVsNsigma_notAveraged.Write()
    canvas.Write()
