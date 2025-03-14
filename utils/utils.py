import ROOT
import uproot
import pandas as pd
import numpy as np
from scipy import special

ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
ROOT.gStyle.SetOptStat(0)

ROOT.gStyle.SetFrameLineColor(ROOT.gStyle.GetCanvasColor())

mass_alpha = 3.7273794066
mass_helion = 2.80839160743


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
    max_cluster_size = 0
    for i_layer in range(0, 7):
        size = (itsSizeBitmap >> i_layer * 4) & 15
        max_cluster_size = np.maximum(max_cluster_size, size)
        if size > 0:
            nClus = nClus + 1
            sum += size
    return (sum - max_cluster_size) / (nClus - 1)


# vectorised version of getITSClSize
getITSClSize_vectorised = np.vectorize(getITSClSize)


# get the sign of the track from the flag
def getSign(flags):
    if flags & 256:
        return 1
    else:
        return -1


# vectorised version of getSign
getSign_vectorised = np.vectorize(getSign)


# get pid-in-tracking flag
def trackedAsHe(flags):
    pid_flag = (flags >> 12) & 15
    if pid_flag == 7 or pid_flag == 8:
        return True
    else:
        return False


# vectorised version of getSign
trackedAsHe_vectorised = np.vectorize(trackedAsHe)


# Bethe-Bloch-Aleph function as defined in O2
def BBA(x, p):

    bg = x
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
func_string_to_format = (
    r"([1] - TMath::Power((TMath::Abs({charge}*x/{mass}) / "
    r"TMath::Sqrt(1 + {charge}*{charge}*x*x/{mass}/{mass})),[3]) - "
    r"TMath::Log([2] + TMath::Power(1/TMath::Abs({charge}*x/{mass}),[4]))) * [0] / "
    r"TMath::Power((TMath::Abs({charge}*x/{mass}) / TMath::Sqrt(1 + "
    r"{charge}*{charge}*x*x/{mass}/{mass})),[3])"
)

default_bb_parameters = p_train = [
    -449.90,
    1.1210,
    1.9451,
    0.7336,
    2.0896,
]  # [-321.34, 0.6539, 1.591, 0.8225, 2.363]


# N-sigma TPC at specific rigidity
def getNsigmaTPC(
    x,
    tpc_signal,
    parameters=default_bb_parameters,
    resolution_perc=0.09,
    mass=mass_helion,
):
    exp_signal = BBA(x / mass, parameters)
    resolution = exp_signal * resolution_perc
    return (tpc_signal - exp_signal) / resolution


# vectorised version of N-sigma TPC at specific rigidity
getNsigmaTPC_vectorised = np.vectorize(getNsigmaTPC)

# its n-sigma

its_signal_parameters = [2.3512, 1.8035, 5.1436]
its_resolution_parameters = [0.0874, -1.8280, 0.5064]


def expSignal(p, mass=mass_helion):
    inverseMass = 1.0 / mass
    bg = abs(p) * inverseMass
    return (
        its_signal_parameters[0] / (np.power(bg, its_signal_parameters[1]))
        + its_signal_parameters[2]
    )


def expResolution(p, mass=mass_helion):
    inverseMass = 1.0 / mass
    bg = abs(p) * inverseMass
    relRes = its_resolution_parameters[0] * special.erf(
        (bg - its_resolution_parameters[1]) / its_resolution_parameters[2]
    )
    return relRes


def getNsigmaITS(avgItsClusterSizeCosLambda, p):
    exp = expSignal(p)
    resolution = expResolution(p) * exp
    return (avgItsClusterSizeCosLambda - exp) / resolution


# vectorised version of N-sigma TPC at specific rigidity
getNsigmaITS_vectorised = np.vectorize(getNsigmaITS)

# pt-dependent its n-sigma selection
its_nsigma_parameters = [4.6, 0.5, -4.5]
its_nsigma_func_string = r"([0] / TMath::Power(x, [1])) + [2]"


def nSigmaITSoffset(p):
    return its_nsigma_parameters[2] + (
        its_nsigma_parameters[0] / pow(abs(p), its_nsigma_parameters[1])
    )


def getNsigmaITSminusOffset(nsigma, p):
    offset = nSigmaITSoffset(p)
    diff = nsigma - offset
    return diff


# vectorised version of N-sigma TPC at specific rigidity
getNsigmaITSminusOffset_vectorised = np.vectorize(getNsigmaITSminusOffset)


# get rapidity
def getRapidity(pt, eta, phi, mass=mass_helion):
    vector = ROOT.TLorentzVector()
    vector.SetPtEtaPhiM(pt, eta, phi, mass)
    return vector.Rapidity()


# vectorised version of getRapidity
getRapidity_vectorised = np.vectorize(getRapidity)


# return phi in [-pi, pi]
def getCorrectPhi(phi):
    new_phi = phi
    if phi > np.pi:
        new_phi = new_phi - (2 * np.pi)
    return new_phi


getCorrectPhi_vectorised = np.vectorize(getCorrectPhi)


def get_df_from_tree(input_file_name, tree_name):

    # create empty data-frame
    df = pd.DataFrame()

    # list of folders in the file
    file_folders = uproot.open(input_file_name).keys()

    # filter folders: only folders without tree_name, not containing "/", and ending with ";1"
    file_folders = [folder for folder in file_folders if "/" not in folder]

    # first we sort to have as first one the last cycle
    file_folders.sort(reverse=True)

    # check if there are multiple cycles of the same tree, keep only last one
    file_folders_to_remove = []
    for ifolder, folder in enumerate(file_folders[1:]):
        obj_nocycle = folder.split(";")[0]
        if obj_nocycle in file_folders[ifolder]:
            file_folders_to_remove.append(folder)
    for folder_to_remove in file_folders_to_remove:
        file_folders.remove(folder_to_remove)

    # loop over the folders
    for folder in file_folders:
        df = pd.concat(
            [
                df,
                uproot.open(f"{input_file_name}:{folder}/{tree_name}").arrays(
                    library="pd"
                ),
            ],
            ignore_index=True,
            copy=False,
        )
    return df


# redefine columns in the complete data-frame
def redefineColumns(
    complete_df,
    mass=mass_helion,
    charge=2,
    parameters=default_bb_parameters,
    useSP=False,
):
    print("Redefining columns")
    print("fPt")
    complete_df["fPt"] = charge * complete_df["fPt"]
    print("fP")
    complete_df.eval("fP = fPt * cosh(fEta)", inplace=True)
    print("fPhi")
    complete_df["fPhi"] = getCorrectPhi_vectorised(complete_df["fPhi"])
    print("fCosLambda")
    complete_df.eval("fCosLambda = 1 / cosh(fEta)", inplace=True)
    print("fAvgItsClusSize")
    complete_df.eval(
        "fAvgItsClusSize = @getITSClSize_vectorised(fITSclusterSizes)", inplace=True
    )
    print("fAvgItsClusSize * fCosLambda")
    complete_df.eval(
        "fAvgItsClusSizeCosLambda = fAvgItsClusSize * fCosLambda", inplace=True
    )
    complete_df.drop(columns=["fITSclusterSizes"])
    print("fSign")
    complete_df.eval("fSign = @getSign_vectorised(fFlags)", inplace=True)
    print("fNsigmaTPC3He")
    complete_df["fNsigmaTPC3He"] = complete_df.apply(
        lambda row: getNsigmaTPC(
            charge * row["fTPCInnerParam"],
            row["fTPCsignal"],
            mass=mass,
            parameters=parameters,
        ),
        axis=1,
    )
    print("fNsigmaITS3He")
    complete_df["fNsigmaITS3He"] = complete_df.apply(
        lambda row: getNsigmaITS(row["fAvgItsClusSizeCosLambda"], row["fP"]),
        axis=1,
    )
    print("fNsigmaITS3HeMinusOffset")
    complete_df["fNsigmaITS3HeMinusOffset"] = complete_df.apply(
        lambda row: getNsigmaITSminusOffset(row["fNsigmaITS3He"], row["fP"]),
        axis=1,
    )
    # print('fRapidity')
    # complete_df.eval(
    #     'fRapidity = @getRapidity_vectorised(fPt, fEta, fPhi)', inplace=True)
    print("fTOFmassSquared")
    complete_df.eval(
        f"fTOFmassSquared = {charge} * {charge} * fTPCInnerParam * fTPCInnerParam * (1/(fBeta * fBeta) - 1)",
        inplace=True,
    )
    # v2 with event-plane method
    if useSP:
        print("Using Scalar Product method")
        print("fV2FT0C")
        complete_df.eval("fV2FT0C = fQFT0C * cos(2 * (fPhi-fPsiFT0C))", inplace=True)
        print("fV2FT0A")
        complete_df.eval("fV2FT0A = fQFT0A * cos(2 * (fPhi-fPsiFT0A))", inplace=True)
        print("fV2TPCl")
        complete_df.eval("fV2TPCl = fQTPCl * cos(2 * (fPhi-fPsiTPCl))", inplace=True)
        print("fV2TPCr")
        complete_df.eval("fV2TPCr = fQTPCr * cos(2 * (fPhi-fPsiTPCr))", inplace=True)
        print("fV2TPC")
        complete_df.eval("fV2TPC = fQTPC * cos(2 * (fPhi-fPsiTPC))", inplace=True)
    else:
        print("Using Event Plane method")
        print("fV2FT0C")
        complete_df.eval("fV2FT0C = cos(2 * (fPhi-fPsiFT0C))", inplace=True)
        print("fV2FT0A")
        complete_df.eval("fV2FT0A = cos(2 * (fPhi-fPsiFT0A))", inplace=True)
        print("fV2TPCl")
        complete_df.eval("fV2TPCl = cos(2 * (fPhi-fPsiTPCl))", inplace=True)
        print("fV2TPCr")
        complete_df.eval("fV2TPCr = cos(2 * (fPhi-fPsiTPCr))", inplace=True)
        print("fV2TPC")
        complete_df.eval("fV2TPC = cos(2 * (fPhi-fPsiTPC))", inplace=True)


def redefineColumnsLight(complete_df, charge=2):
    print("Redefining columns")
    print("fPt")
    complete_df["fPt"] = charge * complete_df["fPt"]
    print("fP")
    complete_df.eval("fP = fPt * sinh(fEta)", inplace=True)
    print("fPhi")
    complete_df["fPhi"] = getCorrectPhi_vectorised(complete_df["fPhi"])
    print("fCosLambda")
    complete_df.eval("fCosLambda = 1 / cosh(fEta)", inplace=True)
    print("fAvgItsClusSize")
    complete_df.eval(
        "fAvgItsClusSize = @getITSClSize_vectorised(fITSclusterSizes)", inplace=True
    )
    print("fAvgItsClusSize * fCosLambda")
    complete_df.eval(
        "fAvgItsClusSizeCosLambda = fAvgItsClusSize * fCosLambda", inplace=True
    )
    print("fNsigmaITS3He")
    complete_df["fNsigmaITS3He"] = complete_df.apply(
        lambda row: getNsigmaITS(row["fAvgItsClusSizeCosLambda"], row["fP"]),
        axis=1,
    )
    complete_df.drop(columns=["fITSclusterSizes"])
    complete_df.drop(columns=["fAvgItsClusSize"])
    complete_df.drop(columns=["fAvgItsClusSizeCosLambda"])
    print("fSign")
    complete_df.eval("fSign = @getSign_vectorised(fFlags)", inplace=True)


def getBBAfunctions(
    parameters,
    resolution,
    n_sigma=5,
    mass=mass_helion,
    charge=2,
    color=ROOT.kRed,
    line_width=2,
):
    upper_scale = 1 + resolution * n_sigma
    lower_scale = 1 - resolution * n_sigma

    # Formattazione stringa con i valori corretti
    func_string = func_string_to_format.format(charge=charge, mass=mass)

    func_string_up = f"{upper_scale} * " + func_string
    func_string_down = f"{lower_scale} * " + func_string

    # Creazione delle funzioni
    func_BB_left = ROOT.TF1("func_BB_left", func_string, -6, -0.5, 5)
    func_BB_left.SetParameters(*parameters)
    func_BB_left.SetLineColor(color)
    func_BB_left.SetLineWidth(line_width)

    func_BB_left_up = ROOT.TF1("func_BB_left_up", func_string_up, -6, -0.5, 5)
    func_BB_left_up.SetParameters(*parameters)
    func_BB_left_up.SetLineColor(color)
    func_BB_left_up.SetLineStyle(ROOT.kDashed)
    func_BB_left_up.SetLineWidth(line_width)

    func_BB_left_down = ROOT.TF1("func_BB_left_down", func_string_down, -6, -0.5, 5)
    func_BB_left_down.SetParameters(*parameters)
    func_BB_left_down.SetLineColor(color)
    func_BB_left_down.SetLineStyle(ROOT.kDashed)
    func_BB_left_down.SetLineWidth(line_width)

    func_BB_right = ROOT.TF1("func_BB_right", func_string, 0.5, 6.0, 5)
    func_BB_right.SetParameters(*parameters)
    func_BB_right.SetLineColor(color)
    func_BB_right.SetLineWidth(line_width)

    func_BB_right_up = ROOT.TF1("func_BB_right_up", func_string_up, 0.5, 6.0, 5)
    func_BB_right_up.SetParameters(*parameters)
    func_BB_right_up.SetLineColor(color)
    func_BB_right_up.SetLineStyle(ROOT.kDashed)
    func_BB_right_up.SetLineWidth(line_width)

    func_BB_right_down = ROOT.TF1("func_BB_right_down", func_string_down, 0.5, 6.0, 5)
    func_BB_right_down.SetParameters(*parameters)
    func_BB_right_down.SetLineColor(color)
    func_BB_right_down.SetLineStyle(ROOT.kDashed)
    func_BB_right_down.SetLineWidth(line_width)

    functions = [
        func_BB_left,
        func_BB_left_up,
        func_BB_left_down,
        func_BB_right,
        func_BB_right_up,
        func_BB_right_down,
    ]

    return functions


def getAverage2D(histo2D, histo_name="histo"):
    x_axis = histo2D.GetXaxis()
    nX_bins = x_axis.GetNbins()
    histo2D_y_title = histo2D.GetYaxis().GetTitle()
    new_title = r"#LT " + histo2D_y_title + r" #GT"
    average_histo = histo2D.ProjectionX(histo_name)
    average_histo.Reset()
    average_histo.GetYaxis().SetTitle(new_title)
    for i_bin in range(1, nX_bins + 1):
        histo_tmp = histo2D.ProjectionY("histo_tmp", i_bin, i_bin)
        average_histo.SetBinContent(i_bin, histo_tmp.GetMean())
        average_histo.SetBinError(i_bin, histo_tmp.GetStdDev())
        del histo_tmp
    return average_histo


def getHistos1D(
    thn_sparse,
    pt_bin_low,
    pt_bin_up,
    cent_low,
    cent_up,
    histo1_name="histo1",
    histo2_name="histo2",
    n_rebin=1,
):
    pt_axis = thn_sparse.GetAxis(2)
    pt_axis.SetRange(pt_bin_low, pt_bin_up)
    cent_axis = thn_sparse.GetAxis(3)
    cent_axis.SetRange(cent_axis.FindBin(cent_low), cent_axis.FindBin(cent_up) - 1)
    hNsigma = thn_sparse.Projection(1)
    hNsigma.SetName(histo1_name)
    hNsigma.SetTitle("")
    hSPvsNsigma_tmp = thn_sparse.Projection(0, 1)
    hSPvsNsigma_tmp.SetName(f"{histo2_name}_notAveraged")
    if n_rebin > 1:
        hNsigma.RebinX(n_rebin)
        hSPvsNsigma_tmp.RebinX(n_rebin)
    hSPvsNsigma = getAverage2D(hSPvsNsigma_tmp, histo2_name)
    hSPvsNsigma.SetTitle("")
    return hSPvsNsigma, hNsigma, hSPvsNsigma_tmp


def getCanvasWithTwoPanels(
    canvas_name, histo_1, histo_2, top_panel=None, bottom_panel=None, line=None
):
    canvas = ROOT.TCanvas(canvas_name, canvas_name, 800, 600)
    pad_top = ROOT.TPad("pad_top", "pad_top", 0.0, 0.5, 1.0, 1.0, 0)
    pad_top.SetLeftMargin(0.15)
    pad_top.SetBottomMargin(0.0)
    pad_top.Draw()
    pad_top.cd()
    histo_1.Draw("PE")
    if top_panel:
        histo_1.GetListOfFunctions().Add(top_panel)
    canvas.cd()
    pad_bottom = ROOT.TPad("pad_bottom", "pad_bottom", 0.0, 0.0, 1.0, 0.5, 0)
    pad_bottom.SetLeftMargin(0.15)
    pad_bottom.SetTopMargin(0.0)
    pad_bottom.SetBottomMargin(0.3)
    pad_bottom.Draw()
    pad_bottom.cd()
    histo_2.Draw("PE")
    if bottom_panel:
        histo_2.GetListOfFunctions().Add(bottom_panel)
    if line:
        line.Draw()
    return canvas


def getCompleteCanvas(
    hSpvsNsigmaVsPtvsCent,
    cent_low,
    cent_up,
    pt_bin_low,
    pt_bin_up,
    output_dir,
    hSPvsPt,
    qvec_detector_label="FT0C",
    cent_detector_label="FTOC",
    out_pt_bin=0,
):

    output_dir.cd()

    # get axis info
    pt_axis = hSpvsNsigmaVsPtvsCent.GetAxis(2)
    pt_low = pt_axis.GetBinLowEdge(pt_bin_low)
    pt_up = pt_axis.GetBinUpEdge(pt_bin_up)

    cent_axis = hSpvsNsigmaVsPtvsCent.GetAxis(3)
    cent_bin_low = cent_axis.FindBin(cent_low)
    cent_bin_up = cent_axis.FindBin(cent_up) - 1

    # system infopanel
    info_panel = ROOT.TPaveText(0.6, 0.6, 0.8, 0.82, "NDC")
    info_panel.SetBorderSize(0)
    info_panel.SetFillStyle(0)
    info_panel.SetTextAlign(12)
    info_panel.SetTextFont(42)
    info_panel.AddText("ALICE")
    info_panel.AddText(r"PbPb, #sqrt{#it{s}_{nn}} = 5.36 TeV")
    info_panel.AddText(f"{cent_low} - {cent_up} % {cent_detector_label}")
    pt_label = (
        f"{pt_low:.1f}" + r" #leq #it{p}_{T} < " + f"{pt_up:.1f}" + r" GeV/#it{c}"
    )
    info_panel.AddText(pt_label)

    # create 1D histograms
    hSpVsNsigma, hNsigma, hSpVsNsigma_notAveraged = getHistos1D(
        hSpvsNsigmaVsPtvsCent,
        pt_bin_low,
        pt_bin_up,
        cent_bin_low,
        cent_bin_up,
        f"hNsigma_{out_pt_bin}_{qvec_detector_label}",
        f"hSpVsNsigma_{out_pt_bin}_{qvec_detector_label}",
        n_rebin=4,
    )
    setHistStyle(hNsigma, ROOT.kRed + 1, linewidth=2)
    setHistStyle(hSpVsNsigma, ROOT.kAzure + 1, linewidth=2)
    hSpVsNsigma.GetYaxis().SetRangeUser(-2.0, 2.0)

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
        f"cSpVsNsigma_{out_pt_bin}_{qvec_detector_label}",
        hSpVsNsigma,
        hNsigma,
        bottom_panel=info_panel,
    )

    hSpVsNsigma.Write()
    hNsigma.Write()
    hSpVsNsigma_notAveraged.Write()
    canvas.Write()


def saveCanvasAsPDF(histo, plots_dir, is2D=False, logScale=False):
    histo_name = histo.GetName()
    canvas_name = histo_name.replace("h", "c", 1)
    canvas = ROOT.TCanvas(canvas_name, canvas_name, 800, 600)
    canvas.SetBottomMargin(0.13)
    canvas.SetLeftMargin(0.13)
    if not is2D:
        histo.Draw("histo")
    else:
        histo.Draw("colz")
    if logScale:
        if not is2D:
            canvas.SetLogy()
        else:
            canvas.SetLogz()
    canvas.SaveAs(f"{plots_dir}/{canvas_name}.pdf")


def getValuesFromHisto(histo):
    n_bins = histo.GetXaxis().GetNbins()
    histo_content = []
    for i_bin in range(1, n_bins + 1):
        histo_content.append([histo.GetBinContent(i_bin), histo.GetBinError(i_bin)])
    return histo_content


def passBarlow(def_val, varied_val, def_err, varied_err, n_sigma=1):
    numerator = abs(def_val - varied_val)
    denominator = abs(def_err - varied_err)
    if denominator < 1.0e-5:
        return True
    return (numerator / denominator) < n_sigma


def get_condition(val, conditions):
    """
    @brief returns the first condition for which val < conditions.key

    @param val value to be tested
    @param conditions dictionary with form {right_edge: string}, with keys in ascending order

    """
    # Sort the dictionary keys
    sorted_keys = sorted(conditions.keys())

    # Loop through the keys and find the last key that satisfies val < key
    for key in sorted_keys:
        if val < key:
            return conditions[key]

    # If no condition is found, return None or a default message
    return None
