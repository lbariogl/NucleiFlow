import utils as utils
import os
import ROOT
import numpy as np

import sys
sys.path.append('utils')

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
        self.tau_background = None
        self.tau_background_limits = [-0.3, -10, -1.e-5]

        self.model_background = None

        # total
        self.total_model = None


        self.fit_results = None
        self.tot_integral = 0.
        self.tot_integral_error = 0.
        self.bkg_integral = 0.
        self.bkg_integral_error = 0.
        self.purity = 1.
        self.purity_error = 0.

        # integral range
        self.integral_range = []

        self.frame = None
        self.info_panel = None
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
        if self.tau_background is None:
            raise ValueError(f'tau_background_background not correctly set')
        if self.model_background is None:
            raise ValueError(f'model_background_background not correctly set')

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
        self.tau_background = ROOT.RooRealVar(
            r'#tau_{bkg}', r'#tau_{bkg}', self.tau_background_limits[0], self.tau_background_limits[1], self.tau_background_limits[2])

        self.model_background = ROOT.RooExponential(
            'background', 'background', self.variable, self.tau_background)

        # total model
        self.total_model = ROOT.RooAddPdf('model', 'total model', ROOT.RooArgList(
            self.model_signal, self.model_background), ROOT.RooArgList(self.n_signal, self.n_background))

    def fit(self):
        self._check_members()

        data_hist = ROOT.RooDataHist(f'{self.hist.GetName()}_roofit', 'n-sigma distribution',
                                     ROOT.RooArgList(self.variable), ROOT.RooFit.Import(self.hist))
        self.fit_results = self.total_model.fitTo(data_hist, ROOT.RooFit.Extended(), ROOT.RooFit.Verbose(
            False), ROOT.RooFit.PrintEvalErrors(-1), ROOT.RooFit.PrintLevel(-1), ROOT.RooFit.Range(self.variable_range[0], self.variable_range[1]), ROOT.RooFit.Save())
        self.frame = self.variable.frame(
            self.variable_range[0], self.variable_range[1])
        frame_title = f'{self.cent_label}, {self.pt_label}'
        self.frame.SetTitle(frame_title)
        histo_name = self.hist.GetName()
        frame_name = histo_name.replace('h', 'fr', 1)
        self.frame.SetName(frame_name)
        data_hist.plotOn(self.frame, ROOT.RooFit.Name(
            'data'), ROOT.RooFit.DrawOption('pz'))

        # background line
        self.total_model.plotOn(self.frame, ROOT.RooFit.Components('background'), ROOT.RooFit.LineStyle(
            ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed+2), ROOT.RooFit.Name('background'))

        # Total line
        self.total_model.plotOn(self.frame, ROOT.RooFit.DrawOption(
            ''), ROOT.RooFit.LineColor(ROOT.kAzure+2), ROOT.RooFit.Name('model'))

        # purity determination

        # total integral in the signal region
        self.tot_integral = self.hist.Integral(self.hist.FindBin(self.integral_range[0]), self.hist.FindBin(self.integral_range[1]))
        self.tot_integral_error = np.sqrt(self.tot_integral)

        # background integral in the signal region
        bkg_par_list = ROOT.RooArgList(self.tau_background, self.n_background)
        bkg_model = ROOT.RooAddPdf('bkg_model', 'bkg function TOF', ROOT.RooArgList(self.model_background), ROOT.RooArgList(self.n_background))
        bkg_model.fixCoefNormalization(ROOT.RooArgSet(self.variable))
        func_bkg = bkg_model.asTF(ROOT.RooArgList(self.variable), bkg_par_list)
        bkg_integral_full = func_bkg.Integral(self.variable_range[0], self.variable_range[1])
        self.bkg_integral = self.n_background.getVal() * func_bkg.Integral(self.integral_range[0], self.integral_range[1], 0) / bkg_integral_full
        self.bkg_integral_error = np.sqrt(self.bkg_integral)

        # purity determination
        self.purity = 1 - (self.bkg_integral / self.tot_integral)
        self.purity_error = np.hypot(self.bkg_integral_error/self.bkg_integral, self.tot_integral_error/self.tot_integral)

        # create canvas
        canvas_name = frame_name.replace('fr', 'c', 1) + '_Fit'
        self.canvas = ROOT.TCanvas(canvas_name, canvas_name, 800, 600)
        self.canvas.cd()

        self.info_panel = ROOT.TPaveText(0.7, 0.6, 0.9, 0.6, 'brNDC')
        self.info_panel.SetBorderSize(0)
        self.info_panel.SetFillStyle(0)
        self.info_panel.SetTextAlign(12)
        self.info_panel.SetTextFont(42)
        self.info_panel.AddText(r'#mu = ' + f'{self.mu_signal.getVal():.3f} #pm {self.mu_signal.getError():.3f}')
        self.info_panel.AddText(r'#sigma = ' + f'{self.sigma_signal.getVal():.3f} #pm {self.sigma_signal.getError():.3f}')
        self.info_panel.AddText(r'#tau = ' + f'{self.tau_background.getVal():.3f} #pm {self.tau_background.getError():.3f}')
        self.info_panel.AddText(r'N_{tot}' + f'([{self.integral_range[0]}, {self.integral_range[1]}]) =  {self.tot_integral:.1f} #pm {self.tot_integral_error:.1f}')
        self.info_panel.AddText(r'N_{bkg}' + f'([{self.integral_range[0]}, {self.integral_range[1]}]) = {self.bkg_integral:.1f} #pm {self.bkg_integral_error:.1f}')
        self.info_panel.AddText(r'purity' + f'([{self.integral_range[0]}, {self.integral_range[1]}]) = {self.purity:.2f}')

        self.frame.Draw()
        self.canvas.cd()
        self.canvas.Update()
        self.info_panel.Draw()
        self.canvas.Update()
        self.canvas.SetLogy()

    def saveFrameAsPDF(self, plots_dir):
        self.canvas.SaveAs(f'{plots_dir}/{self.canvas.GetName()}.pdf')
