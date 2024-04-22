import os
import ROOT
import numpy as np

import sys
sys.path.append('utils')
import utils as utils

ROOT.ROOT.EnableImplicitMT()
ROOT.RooMsgService.instance().setSilentMode(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

class FlowMaker:

    def __init__(self):

        # data frame
        self.data_df = None
        self.selection_string = ''
        self.n_sigma_selection = 2

        # variable related members
        self.cent_limits = [30, 50]
        self.pt_bins = []
        self.cent_detector = 'FT0C'
        self.ref_detector = 'FT0C'

        # histograms
        self.n_nsigmaTPC_bins = 50
        self.nsigmaTPC_bin_limits = [-5., 5.]

        self.n_tofMassSquared_bins = 50
        self.tofMassSquared_bin_limits = [-5., 5.]

        self.n_V2_bins = 20
        self.n_V2_bin_limits = [-1., 1.]

        self.hNigma3He = []
        self.frNigma3He = []
        self.hTOFmassSquared = []
        self.hV2vsNsigma3He2D = []
        self.hV2vsNsigma3He = []
        self.cV2vsNsigma3He = []
        self.hV2 = []
        self.hV2vsPt = None
        self.hRawCountsVsPt = None
        self.color = ROOT.kBlack
        self.suffix = ''

        # general
        self.output_dir = None
        self.plot_dir = None

        # resolution
        self.resolution = 1.

    def _check_members(self):

        if self.data_df is None:
            raise ValueError(f'data frame not correctly set')

        if not self.cent_limits:
            raise ValueError(f'centrality limits not set')

        if not self.pt_bins:
            raise ValueError(f'pt bins not set')

    def make_flow(self):

        self._check_members()

        n_pt_bins = len(self.pt_bins) - 1
        pt_bins_arr = np.array(self.pt_bins, dtype=np.float64)
        self.hRawCountsVsPt = ROOT.TH1F(
            f'hRawCountsVsPt_cent_{self.cent_limits[0]}_{self.cent_limits[1]}{self.suffix}', r'; #it{p}_{T} (GeV/#it{c}); counts', n_pt_bins, pt_bins_arr)
        utils.setHistStyle(self.hRawCountsVsPt, self.color)

        for i_pt in range(0, n_pt_bins):
            pt_bin = [self.pt_bins[i_pt], self.pt_bins[i_pt + 1]]
            pt_sel = f'abs(fPt) > {pt_bin[0]} and abs(fPt) < {pt_bin[1]}'
            cent_sel = f'fCentFT0C > {self.cent_limits[0]} and fCentFT0C < {self.cent_limits[1]}'
            if self.selection_string:
                bin_sel = f'{cent_sel} and {pt_sel} and {self.selection_string}'
            else:
                bin_sel = f'{cent_sel} and {pt_sel}'

            pt_label = f'{pt_bin[0]:.2f} ' + r'#leq #it{p}_{T} < ' + \
                f'{pt_bin[1]:.2f}' + r' GeV/#it{c}'

            # select the correct pt bin
            bin_df = self.data_df.query(bin_sel, inplace=False)

            self.hRawCountsVsPt.SetBinContent(i_pt+1, len(bin_df))
            self.hRawCountsVsPt.SetBinError(i_pt+1, np.sqrt(len(bin_df)))

            # crate and fill histograms
            hNsigma3He_tmp = ROOT.TH1F(f'hNsigma3He_cent_{self.cent_limits[0]}_{self.cent_limits[1]}_pt{i_pt}{self.suffix}', pt_label + r';n#sigma^{TPC} (a.u.);',
                                       self.n_nsigmaTPC_bins, self.nsigmaTPC_bin_limits[0], self.nsigmaTPC_bin_limits[1])
            # hNsigma3He_tmp.SetTitle(pt_label)
            utils.setHistStyle(hNsigma3He_tmp, ROOT.kRed+1, linewidth=2)
            hTOFmassSquared_tmp = ROOT.TH1F(f'hTOFmassSquared_cent_{self.cent_limits[0]}_{self.cent_limits[1]}_pt{i_pt}{self.suffix}', pt_label + r';m_{TOF}^{2} - m_{0}^{2} (GeV/#it{c}^{2})^{2};',
                                            self.n_tofMassSquared_bins, self.tofMassSquared_bin_limits[0], self.tofMassSquared_bin_limits[1])
            utils.setHistStyle(hTOFmassSquared_tmp, ROOT.kRed+1, linewidth=2)
            histo_title = r';n#sigma^{TPC} (a.u.); cos(2(#phi - #Psi))'
            hV2vsNsigma3He2D_tmp = ROOT.TH2F(f'hV2{self.ref_detector}vsNsigma3He2D_cent_{self.cent_limits[0]}_{self.cent_limits[1]}_pt{i_pt}{self.suffix}', histo_title,
                                             self.n_nsigmaTPC_bins, self.nsigmaTPC_bin_limits[0], self.nsigmaTPC_bin_limits[1], self.n_V2_bins, self.n_V2_bin_limits[0], self.n_V2_bin_limits[1])

            for nSigmaTPC, v2, m2 in zip(bin_df['fNsigmaTPC3He'], bin_df[f'fV2{self.ref_detector}'], bin_df['fTOFmassSquared']):
                hNsigma3He_tmp.Fill(nSigmaTPC)
                hV2vsNsigma3He2D_tmp.Fill(nSigmaTPC, v2)
                hTOFmassSquared_tmp.Fill(m2 - 2.80839160743 * 2.80839160743)

            hV2vsNsigma3He_tmp = utils.getAverage2D(
                hV2vsNsigma3He2D_tmp, f'hV2{self.ref_detector}vsNsigma3He_cent_{self.cent_limits[0]}_{self.cent_limits[1]}_pt{i_pt}{self.suffix}')
            utils.setHistStyle(hV2vsNsigma3He_tmp, ROOT.kAzure+1, linewidth=2)

            self.hNigma3He.append(hNsigma3He_tmp)
            self.hTOFmassSquared.append(hTOFmassSquared_tmp)
            self.hV2vsNsigma3He2D.append(hV2vsNsigma3He2D_tmp)
            self.hV2vsNsigma3He.append(hV2vsNsigma3He_tmp)

            # system infopanel
            info_panel = ROOT.TPaveText(0.6, 0.6, 0.8, 0.82, 'NDC')
            info_panel.SetBorderSize(0)
            info_panel.SetFillStyle(0)
            info_panel.SetTextAlign(12)
            info_panel.SetTextFont(42)
            info_panel.AddText('ALICE')
            info_panel.AddText(r'PbPb, #sqrt{#it{s}_{nn}} = 5.36 TeV')
            info_panel.AddText(
                f'{self.cent_limits[0]} - {self.cent_limits[1]} % {self.cent_detector}')
            info_panel.AddText(pt_label)

            canvas = utils.getCanvasWithTwoPanels(
                f'cV2{self.ref_detector}vsNsigma_cent_{self.cent_limits[0]}_{self.cent_limits[1]}_pt{i_pt}{self.suffix}', hV2vsNsigma3He_tmp, hNsigma3He_tmp, bottom_panel=info_panel)

            self.cV2vsNsigma3He.append(canvas)

            hV2_tmp = ROOT.TH1F(f'hV2{self.ref_detector}_cent_{self.cent_limits[0]}_{self.cent_limits[1]}_pt{i_pt}{self.suffix}', f'{pt_label}' + r';cos(2(#phi - #Psi))',
                                self.n_V2_bins, self.n_V2_bin_limits[0], self.n_V2_bin_limits[1])
            df_bin_3He = bin_df.query(
                f'abs(fNsigmaTPC3He) < {self.n_sigma_selection}')
            for v2 in df_bin_3He[f'fV2{self.ref_detector}']:
                hV2_tmp.Fill(v2)

            # system infopanel
            info_panel_bis = ROOT.TPaveText(0.4, 0.6, 0.6, 0.82, 'NDC')
            info_panel_bis.SetBorderSize(0)
            info_panel_bis.SetFillStyle(0)
            info_panel_bis.SetTextAlign(12)
            info_panel_bis.SetTextFont(42)
            info_panel_bis.AddText('ALICE')
            info_panel_bis.AddText(r'PbPb, #sqrt{#it{s}_{nn}} = 5.36 TeV')
            info_panel_bis.AddText(
                f'{self.cent_limits[0]} - {self.cent_limits[1]} % {self.cent_detector}')
            info_panel_bis.AddText(pt_label)

            hV2_tmp.GetListOfFunctions().Add(info_panel_bis)
            utils.setHistStyle(hV2_tmp, ROOT.kAzure+1, linewidth=2)

            self.hV2.append(hV2_tmp)

            # fit to the n-sigma TPC distribution
            # Define TPC variable
            ns = ROOT.RooRealVar(r'n#sigma_{p}', r'n#sigma_{TPC}({}^{3}He})',
                                 self.nsigmaTPC_bin_limits[0], self.nsigmaTPC_bin_limits[1], 'a. u.')

            # 3He signal
            nSig = ROOT.RooRealVar(
                r'N_{sig}', 'counts TPC signal', 1e+2, 0, 1e+4)
            muSig = ROOT.RooRealVar(r'#mu_{sig}', 'mu TPC signal', 0, -1, 1)
            sigmaSig = ROOT.RooRealVar(
                r'#sigma_{sig}', 'sigma TPC signal', 0.1, 0.1, 1.5)

            sigModel = ROOT.RooGaussian(
                'signal', 'signal function tpc', ns, muSig, sigmaSig)

           # bkg
            nBkg = ROOT.RooRealVar(
                r'N_{bkg}', 'counts TOF bkg0', 1e+2, 0, 1e+4)
            tauBkg = ROOT.RooRealVar(
                r'#tau_{bkg}', 'tau TOF bkg0', -0.3, -10, -1.e-5)

            bkgModel = ROOT.RooExponential('bkg', 'bkg function', ns, tauBkg)

            totModel = ROOT.RooAddPdf('model', 'tot function', ROOT.RooArgList(
                sigModel, bkgModel), ROOT.RooArgList(nSig, nBkg))

            frame_title = f'{self.cent_limits[0]} - {self.cent_limits[1]} %, {pt_label}'
            histo_name = hNsigma3He_tmp.GetName()
            frame_name = histo_name.replace('h', 'fr', 1)
            frame, _, _ = utils.plotData(ns, hNsigma3He_tmp, totModel, frame_title, self.nsigmaTPC_bin_limits[0], self.nsigmaTPC_bin_limits[1], frame_name)

            self.frNigma3He.append(frame)

        # v2 vs Pt
        self.hV2vsPt = ROOT.TH1F(
            f'hV2vsPt_cent_{self.cent_limits[0]}_{self.cent_limits[1]}{self.suffix}', r'; #it{p}_{T} (GeV/#it{c}); v_{2}', n_pt_bins, pt_bins_arr)
        utils.setHistStyle(self.hV2vsPt, self.color)
        for i_pt in range(0, len(self.pt_bins) - 1):
            self.hV2vsPt.SetBinContent(i_pt+1, self.hV2[i_pt].GetMean())
            self.hV2vsPt.SetBinError(i_pt+1, self.hV2[i_pt].GetMeanError())

        self.hV2vsPt.Scale(1/self.resolution)

    def getFlowValues(self):
        return utils.getValuesFromHisto(self.hV2vsPt)

    def nPtBins(self):
        return len(self.pt_bins) - 1

    def dump_to_output_file(self):
        self.output_dir.cd()
        for i_pt in range(0, len(self.pt_bins) - 1):
            bin = [self.pt_bins[i_pt], self.pt_bins[i_pt + 1]]
            self.output_dir.mkdir(f'pt_{bin[0]}_{bin[1]}')
            self.output_dir.cd(f'pt_{bin[0]}_{bin[1]}')
            self.hNigma3He[i_pt].Write()
            self.frNigma3He[i_pt].Write()
            self.hTOFmassSquared[i_pt].Write()
            self.hV2vsNsigma3He2D[i_pt].Write()
            self.hV2vsNsigma3He[i_pt].Write()
            self.cV2vsNsigma3He[i_pt].Write()
            self.hV2[i_pt].Write()
        self.output_dir.cd()
        self.hV2vsPt.Write()
        self.hRawCountsVsPt.Write()

    def dump_to_pdf(self):
        if not self.plot_dir:
            raise ValueError('plot_dir not set.')
        if not os.path.exists(self.plot_dir):
            os.makedirs(self.plot_dir)
        for i_pt in range(0, len(self.pt_bins) - 1):
            utils.saveCanvasAsPDF(self.hNigma3He[i_pt], self.plot_dir)
            utils.saveCanvasAsPDF(self.hTOFmassSquared[i_pt], self.plot_dir)
            utils.saveCanvasAsPDF(self.hV2[i_pt], self.plot_dir)
            utils.saveFrameAsPDF(self.frNigma3He[i_pt], self.plot_dir)
        utils.saveCanvasAsPDF(self.hV2vsPt, self.plot_dir)
        utils.saveCanvasAsPDF(self.hRawCountsVsPt, self.plot_dir)

    def del_dyn_members(self):
        self.hNigma3He = []
        self.frNigma3He = []
        self.hTOFmassSquared = []
        self.hV2vsNsigma3He2D = []
        self.hV2vsNsigma3He = []
        self.cV2vsNsigma3He = []
        self.hV2 = []
        self.hV2vsPt = None
        self.hRawCountsVsPt = None
