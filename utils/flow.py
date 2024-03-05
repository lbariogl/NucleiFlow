import ROOT
import numpy as np

import sys
sys.path.append('utils')
import utils as utils


class FlowMaker:

    def __init__(self):

        # data frame
        self.data_df = None

        # variable related members
        self.cent_limits = [30, 50]
        self.pt_bins = []
        self.cent_detector = 'FT0C'
        self.ref_detector = 'FT0C'

        # histograms
        self.n_nsigmaTPC_bins = 100
        self.nsigmaTPC_bin_limits = [-5., 5.]

        self.n_SP_bins = 200
        self.n_SP_bin_limits = [-2., 2.]

        self.hNigma3He = []
        self.hSPvsNsigma3He2D = []
        self.hSPvsNsigma3He = []
        self.cSPvsNsigma3He = []

        # general
        self.output_dir = None

    def _check_members(self):

        if self.data_df is None:
            raise ValueError(f'data frame not correctly set')

        if not self.cent_limits:
            raise ValueError(f'centrality limits not set')

        if not self.pt_bins:
            raise ValueError(f'pt bins not set')

    def make_flow(self):

        self._check_members()

        for ibin in range(0, len(self.pt_bins) - 1):
            pt_bin = [self.pt_bins[ibin], self.pt_bins[ibin + 1]]
            bin_sel = f'abs(fPt) > {pt_bin[0]} and abs(fPt) < {pt_bin[1]}'
            pt_label = f'{pt_bin[0]:.2f} ' + r'#leq #it{p}_{T} < ' + \
                f'{pt_bin[1]:.2f}' + r' GeV/#it{c}'

            # select the correct pt bin
            bin_df = self.data_df.query(bin_sel)

            # crate and fill histograms
            hNsigma3He_tmp = ROOT.TH1F(f'hNsigma3He_{ibin}', r';n#sigma^{TPC} (a.u.);',
                                       self.n_nsigmaTPC_bins, self.nsigmaTPC_bin_limits[0], self.nsigmaTPC_bin_limits[1])
            # hNsigma3He_tmp.SetTitle(pt_label)
            utils.setHistStyle(hNsigma3He_tmp, ROOT.kRed+1, linewidth=2)
            histo_title = r';n#sigma^{TPC} (a.u.); #hat{u}_{2} #upoint #vec{Q}_{2}'
            hSPvsNsigma3He2D_tmp = ROOT.TH2F(f'hSp{self.ref_detector}vsNsigma3He2D_{ibin}', histo_title,
                                             self.n_nsigmaTPC_bins, self.nsigmaTPC_bin_limits[0], self.nsigmaTPC_bin_limits[1], self.n_nsigmaTPC_bins, self.nsigmaTPC_bin_limits[0], self.nsigmaTPC_bin_limits[1])

            for nSigmaTPC, sp in zip(bin_df['fNsigmaTPC3He'], bin_df[f'fSp{self.ref_detector}']):
                hNsigma3He_tmp.Fill(nSigmaTPC)
                hSPvsNsigma3He2D_tmp.Fill(nSigmaTPC, sp)

            hSPvsNsigma3He_tmp = utils.getAverage2D(
                hSPvsNsigma3He2D_tmp, f'hSp{self.ref_detector}vsNsigma3He_{ibin}')
            utils.setHistStyle(hSPvsNsigma3He_tmp, ROOT.kAzure+1, linewidth=2)

            self.hNigma3He.append(hNsigma3He_tmp)
            self.hSPvsNsigma3He2D.append(hSPvsNsigma3He2D_tmp)
            self.hSPvsNsigma3He.append(hSPvsNsigma3He_tmp)

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
                f'cSp{self.ref_detector}vsNsigma_{ibin}', hSPvsNsigma3He_tmp, hNsigma3He_tmp, bottom_panel=info_panel)

            self.cSPvsNsigma3He.append(canvas)

    def dump_to_output_dir(self):
        self.output_dir.cd()
        for ibin in range(0, len(self.pt_bins) - 1):
            bin = [self.pt_bins[ibin], self.pt_bins[ibin + 1]]
            self.output_dir.mkdir(f'pt_{bin[0]}_{bin[1]}')
            self.output_dir.cd(f'pt_{bin[0]}_{bin[1]}')
            self.hNigma3He[ibin].Write()
            self.hSPvsNsigma3He2D[ibin].Write()
            self.hSPvsNsigma3He[ibin].Write()
            self.cSPvsNsigma3He[ibin].Write()
