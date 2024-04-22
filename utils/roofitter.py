import os
import ROOT
import numpy as np

import sys
sys.path.append('utils')
import utils as utils

ROOT.ROOT.EnableImplicitMT()
ROOT.RooMsgService.instance().setSilentMode(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)


class RooFitter:
    def __init__(self):
        self.variable = None
        self.variable_range = []

        self.hist = None

        # signal
        self.n_signal = None
        self.n_signal_limits = [1e+2, 0, 1e+4]
        self.mu_signal = None
        self.mu_signal_limits = [0, -1.5, 1.5]
        self.sigma_signal = None
        self.sigma_signal_limits = [0.1, 0.01, 2]

        self.model_signal = None

        # background
        self.n_background = None
        self.n_background_limits = [1e+2, 0, 1e+4]
        self.tau_bakground = None
        self.tau_bakground_limits = [-0.3, -10, -1.e-5]

        self.model_background = None

        # total
        self.total_model = None

        # integral range
        self.integral_range = []

        self.frame = None
        self.canvas = None

        self.pt_label = ''
        self.cent_label = ''

    def _check_members(self):
        if self.hist is None:
            raise ValueError(f'histogram not set')
        if self.variable is None:
            raise ValueError(f'variable not correctly set')
        if not self.variable_range:
            raise ValueError(f'variable range is empty')

        if self.n_signal is None:
            raise ValueError(f'n_signal not correctly set')
        if self.mu_signal is None:
            raise ValueError(f'mu_signal not correctly set')
        if self.sigma_signal is None:
            raise ValueError(f'sigma_signal not correctly set')
        if self.model_signal is None:
            raise ValueError(f'model_signal not correctly set')

        if self.n_background is None:
            raise ValueError(f'n_background not correctly set')
        if self.tau_bakground is None:
            raise ValueError(f'tau_bakground_background not correctly set')
        if self.model_bakground is None:
            raise ValueError(f'model_bakground_background not correctly set')

        if self.total_model is None:
            raise ValueError(f'total_model not correctly set')

        if not self.integral_range:
            raise ValueError(f'integral_range is empty')

    def initialise(self):

        self.variable = ROOT.RooRealVar(r'n#sigma_{TPC}({}^{3}He)', r'n#sigma_{TPC}({}^{3}He)',
                                        self.variable_range[0], self.variable_range[1], 'a. u.')

        # signal
        self.n_signal = ROOT.RooRealVar(
            r'N_{sig}', r'N_{sig}', self.n_signal_limits[0], self.n_signal_limits[1], self.n_signal_limits[2])
        self.mu_signal = ROOT.RooRealVar(
            r'#mu_{sig}', r'#mu_{sig}', self.mu_signal_limits[0], self.mu_signal_limits[1], self.mu_signal_limits[2])
        self.sigma_signal = ROOT.RooRealVar(
            r'#sigma_{sig}',  r'#sigma_{sig}', self.sigma_signal_limits[0], self.sigma_signal_limits[1], self.sigma_signal_limits[2])

        self.model_signal = ROOT.RooGaussian(
            'signal', 'signal', self.variable, self.mu_signal, self.sigma_signal)

        # background
        self.n_background = ROOT.RooRealVar(
            r'N_{bkg}', r'N_{bkg}', self.n_background_limits[0], self.n_background_limits[1], self.n_background_limits[2])
        self.tau_bakground = ROOT.RooRealVar(
            r'#tau_{bkg}', r'#tau_{bkg}', self.tau_bakground_limits[0], self.tau_bakground_limits[1], self.tau_bakground_limits[2])

        self.model_background = ROOT.RooExponential(
            'background', 'background', self.variable, self.tau_bakground)

        # total model
        self.total_model = ROOT.RooAddPdf('model', 'total model', ROOT.RooArgList(
            self.model_signal, self.model_background), ROOT.RooArgList(self.n_signal, self.n_background))

    def fit(self):
        data_hist = ROOT.RooDataHist(f'{self.hist.GetName()}_roofit', 'n-sigma distribution',
                                 ROOT.RooArgList(self.variable), ROOT.RooFit.Import(self.hist))
        fit_results = self.model.fitTo(data_hist, ROOT.RooFit.Extended(), ROOT.RooFit.Verbose(
            False), ROOT.RooFit.PrintEvalErrors(-1), ROOT.RooFit.PrintLevel(-1), ROOT.RooFit.Range(self.variable_range[0], self.variable_range[1]), ROOT.RooFit.Save())
        self.frame = self.variable.frame(self.variable_range[0], self.variable_range[1])
        frame_title = f'{self.cent_label}, {self.pt_label}'
        self.frame.SetTitle(frame_title)
        histo_name = self.hist.GetName()
        frame_name = histo_name.replace('h', 'fr', 1)
        self.frame.SetName(frame_name)
        data_hist.plotOn(self.frame, ROOT.RooFit.Name(
            'data'), ROOT.RooFit.DrawOption('pz'))

        # background line
        self.model.plotOn(self.frame, ROOT.RooFit.Components('bkg'), ROOT.RooFit.LineStyle(
            ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed+2), ROOT.RooFit.Name('bkg'))

        # Total line
        self.model.plotOn(self.frame, ROOT.RooFit.DrawOption(
            ''), ROOT.RooFit.LineColor(ROOT.kAzure+2), ROOT.RooFit.Name('model'))


# # determine purity
#     if bkg_model:
#         # signal histogram
#         tot_integral = hist.Integral(hist.FindBin(
#             integral_range[0]), hist.FindBin(integral_range[1]))
#         tot_integral_err = np.sqrt(tot_integral)

#         # build a copy of bkg function with covariance matrix from fit
#         bkg_par_rooarglist = ROOT.RooArgList(
#             bkg_model_par_list[0], bkg_model_par_list[1])
#         n_bkg = bkg_model_par_list[1]
#         bkg0_cov_matrix = fit_results.reducedCovarianceMatrix(
#             bkg_par_rooarglist).GetMatrixArray()

#         bkg_model_copy = ROOT.RooAddPdf('bkg_model', 'bkg function', ROOT.RooArgList(
#             bkg_model), ROOT.RooArgList(n_bkg))
#         func_bkg = bkg_model_copy.asTF(
#             ROOT.RooArgList(variable), bkg_par_rooarglist)

#         # evaluate background integral
#         bkg_integral_full = func_bkg.Integral(
#             integral_range[0], integral_range[1])
#         bkg_integral = n_bkg.getVal() * func_bkg.Integral(integral_range[0], integral_range[1], 0) / bkg_integral_full
#         bkg_integral_err = n_bkg.getVal() * func_bkg.IntegralError(-3, 3, ROOT.nullptr, bkg0_cov_matrix) / bkg_integral_full
#         bkg_integral_err = np.sqrt(bkg_integral_err)

#         # system infopanel
#         info_panel = ROOT.TPaveText(0.4, 0.6, 0.6, 0.82, 'NDC')
#         info_panel.SetBorderSize(0)
#         info_panel.SetFillStyle(0)
#         info_panel.SetTextAlign(12)
#         info_panel.SetTextFont(42)
#         info_panel.AddText('ALICE')
#         info_panel.AddText(r'PbPb, #sqrt{#it{s}_{nn}} = 5.36 TeV')
#         info_panel.AddText(
#             f'{cent_limits[0]} - {cent_limits[1]} % {cent_detector}')
#         info_panel.AddText(pt_label)
#         info_panel.AddText(f'S+B ({integral_range[0]}' + r'#sigma, ' + f'{integral_range[0]}) = {tot_integral} ' + r'#pm ' + {tot_integral_err})
#         info_panel.AddText(f'B ({integral_range[0]}' + r'#sigma, ' + f'{integral_range[0]}) = {bkg_integral} ' + r'#pm ' + {bkg_integral})
