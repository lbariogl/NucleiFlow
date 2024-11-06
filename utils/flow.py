import roofitter
import utils as utils
import os
import ROOT
import numpy as np
import copy

import sys

sys.path.append("utils")

ROOT.ROOT.EnableImplicitMT()
ROOT.RooMsgService.instance().setSilentMode(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)


class FlowMaker:

    def __init__(self):

        # data frame
        self.data_df = None
        self.selection_string = ""
        self.ptdep_selection_dict = {}
        self.n_sigma_selection = 2
        self.tof_analysis = False

        # defined from self.ptdep_selection_dict
        self.ptdep_selection_string = []

        # variable related members
        self.cent_limits = [30, 50]
        self.pt_bins = []
        self.n_pt_bins = 1

        self.cent_detector = "FT0C"
        self.ref_detector = "FT0C"

        # histograms
        self.n_nsigmaTPC_bins = 50
        self.nsigmaTPC_bin_limits = [-5.0, 5.0]

        self.n_tofMassSquared_bins = 200
        self.tofMassSquared_bin_limits = [-5.0, 5.0]

        self.n_V2_bins = 20
        self.n_V2_bin_limits = [-1.0, 1.0]

        self.hNsigma3He = []
        self.cNsigma3HeFit = []
        self.hTOFmassSquared = []
        self.hV2vsNsigma3He2D = []
        self.hV2vsNsigma3He = []
        self.cV2vsNsigma3He = []
        self.hV2vsTOFmassSquared2D = []
        self.hV2vsTOFmassSquared = []
        self.cV2vsTOFmassSquared = []
        self.hV2 = []
        self.hV2vsPt = None
        self.hRawCountsVsPt = None
        self.color = ROOT.kBlack
        self.suffix = ""

        # general
        self.output_dir = None
        self.plot_dir = None

        # resolution
        self.resolution = 1.0

        # save frame as PDF
        self.print_frame = False

        # purity
        self.purity = []
        self.hPurityVsPt = None

    def _check_members(self):

        if self.data_df is None:
            raise ValueError(f"data frame not correctly set")

        if not self.cent_limits:
            raise ValueError(f"centrality limits not set")

        if not self.pt_bins:
            raise ValueError(f"pt bins not set")

    def create_histograms(self):
        self._check_members()

        self.n_pt_bins = len(self.pt_bins) - 1

        if self.ptdep_selection_dict:
            for i in range(0, self.n_pt_bins):
                bin_centre = (self.pt_bins[i + 1] + self.pt_bins[i]) / 2
                condition = utils.get_condition(bin_centre, self.ptdep_selection_dict)
                self.ptdep_selection_string.append(condition)

        pt_bins_arr = np.array(self.pt_bins, dtype=np.float64)

        # plots vs Pt
        self.hRawCountsVsPt = ROOT.TH1F(
            f"hRawCountsVsPt_cent_{self.cent_limits[0]}_{self.cent_limits[1]}{self.suffix}",
            r"; #it{p}_{T} (GeV/#it{c}); counts",
            self.n_pt_bins,
            pt_bins_arr,
        )
        utils.setHistStyle(self.hRawCountsVsPt, self.color)

        self.hV2vsPt = ROOT.TH1F(
            f"hV2vsPt_cent_{self.cent_limits[0]}_{self.cent_limits[1]}{self.suffix}",
            r"; #it{p}_{T} (GeV/#it{c}); v_{2}",
            self.n_pt_bins,
            pt_bins_arr,
        )
        utils.setHistStyle(self.hV2vsPt, self.color)

        self.hPurityVsPt = ROOT.TH1F(
            f"hPurityVsPt_cent_{self.cent_limits[0]}_{self.cent_limits[1]}{self.suffix}",
            r"; #it{p}_{T} (GeV/#it{c}); purity",
            self.n_pt_bins,
            pt_bins_arr,
        )
        self.hPurityVsPt.GetYaxis().SetRangeUser(0.0, 1.1)
        utils.setHistStyle(self.hPurityVsPt, self.color)

        for i_pt in range(0, self.n_pt_bins):

            pt_bin = [self.pt_bins[i_pt], self.pt_bins[i_pt + 1]]
            pt_label = (
                f"{pt_bin[0]:.2f} "
                + r"#leq #it{p}_{T} < "
                + f"{pt_bin[1]:.2f}"
                + r" GeV/#it{c}"
            )

            # V2
            hV2_tmp = ROOT.TH1F(
                f"hV2{self.ref_detector}_cent_{self.cent_limits[0]}_{self.cent_limits[1]}_pt{i_pt}{self.suffix}",
                f"{pt_label}" + r";cos(2(#phi - #Psi)); counts",
                self.n_V2_bins,
                self.n_V2_bin_limits[0],
                self.n_V2_bin_limits[1],
            )

            # system infopanel
            info_panel_bis = ROOT.TPaveText(0.4, 0.6, 0.6, 0.82, "NDC")
            info_panel_bis.SetBorderSize(0)
            info_panel_bis.SetFillStyle(0)
            info_panel_bis.SetTextAlign(12)
            info_panel_bis.SetTextFont(42)
            info_panel_bis.AddText("ALICE")
            info_panel_bis.AddText(r"PbPb, #sqrt{#it{s}_{NN}} = 5.36 TeV")
            info_panel_bis.AddText(
                f"{self.cent_limits[0]} - {self.cent_limits[1]} % {self.cent_detector}"
            )
            info_panel_bis.AddText(pt_label)

            hV2_tmp.GetListOfFunctions().Add(info_panel_bis)
            utils.setHistStyle(hV2_tmp, ROOT.kAzure + 1, linewidth=2)
            self.hV2.append(hV2_tmp)

            # TPC analysis
            hNsigma3He_tmp = ROOT.TH1F(
                f"hNsigma3He_cent_{self.cent_limits[0]}_{self.cent_limits[1]}_pt{i_pt}{self.suffix}",
                pt_label + r";n#sigma^{TPC} (a.u.); counts",
                self.n_nsigmaTPC_bins,
                self.nsigmaTPC_bin_limits[0],
                self.nsigmaTPC_bin_limits[1],
            )
            utils.setHistStyle(hNsigma3He_tmp, ROOT.kRed + 1, linewidth=2)

            self.hNsigma3He.append(hNsigma3He_tmp)

            histo_title_nsigma_tpc = r";n#sigma^{TPC} (a.u.); cos(2(#phi - #Psi))"
            hV2vsNsigma3He2D_tmp = ROOT.TH2F(
                f"hV2{self.ref_detector}vsNsigma3He2D_cent_{self.cent_limits[0]}_{self.cent_limits[1]}_pt{i_pt}{self.suffix}",
                histo_title_nsigma_tpc,
                self.n_nsigmaTPC_bins,
                self.nsigmaTPC_bin_limits[0],
                self.nsigmaTPC_bin_limits[1],
                self.n_V2_bins,
                self.n_V2_bin_limits[0],
                self.n_V2_bin_limits[1],
            )

            hV2vsNsigma3He_tmp = utils.getAverage2D(
                hV2vsNsigma3He2D_tmp,
                f"hV2{self.ref_detector}vsNsigma3He_cent_{self.cent_limits[0]}_{self.cent_limits[1]}_pt{i_pt}{self.suffix}",
            )
            utils.setHistStyle(hV2vsNsigma3He_tmp, ROOT.kAzure + 1, linewidth=2)

            self.hV2vsNsigma3He2D.append(hV2vsNsigma3He2D_tmp)
            self.hV2vsNsigma3He.append(hV2vsNsigma3He_tmp)

            # TOF analysis
            if self.tof_analysis:
                hTOFmassSquared_tmp = ROOT.TH1F(
                    f"hTOFmassSquared_cent_{self.cent_limits[0]}_{self.cent_limits[1]}_pt{i_pt}{self.suffix}",
                    pt_label
                    + r";m_{TOF}^{2} - m_{{}^{4}He}^{2} (GeV/#it{c}^{2})^{2}; counts",
                    self.n_tofMassSquared_bins,
                    self.tofMassSquared_bin_limits[0],
                    self.tofMassSquared_bin_limits[1],
                )
                utils.setHistStyle(hTOFmassSquared_tmp, ROOT.kRed + 1, linewidth=2)
                self.hTOFmassSquared.append(hTOFmassSquared_tmp)

                histo_title_nsigma_tof = r";m_{TOF}^{2} - m_{{}^{4}He}^{2} (GeV/#it{c}^{2})^{2}; cos(2(#phi - #Psi))"
                hV2vsTOFmassSquared2D_tmp = ROOT.TH2F(
                    f"hV2{self.ref_detector}vsTOFmassSquared2D_cent_{self.cent_limits[0]}_{self.cent_limits[1]}_pt{i_pt}{self.suffix}",
                    histo_title_nsigma_tof,
                    self.n_tofMassSquared_bins,
                    self.tofMassSquared_bin_limits[0],
                    self.tofMassSquared_bin_limits[1],
                    self.n_V2_bins,
                    self.n_V2_bin_limits[0],
                    self.n_V2_bin_limits[1],
                )

                hV2vsTOFmassSquared_tmp = utils.getAverage2D(
                    hV2vsTOFmassSquared_tmp,
                    f"hV2{self.ref_detector}vsTOFmassSquared_cent_{self.cent_limits[0]}_{self.cent_limits[1]}_pt{i_pt}{self.suffix}",
                )
                utils.setHistStyle(
                    hV2vsTOFmassSquared_tmp, ROOT.kAzure + 1, linewidth=2
                )

                self.hV2vsTOFmassSquared2D.append(hV2vsTOFmassSquared2D_tmp)
                self.hV2vsTOFmassSquared.append(hV2vsTOFmassSquared_tmp)

    def make_flow(self):

        self._check_members()

        if self.ptdep_selection_dict:
            for i in range(0, self.n_pt_bins):
                bin_centre = (self.pt_bins[i + 1] + self.pt_bins[i]) / 2
                condition = utils.get_condition(bin_centre, self.ptdep_selection_dict)
                self.ptdep_selection_string.append(condition)

        for i_pt in range(0, self.n_pt_bins):
            pt_bin = [self.pt_bins[i_pt], self.pt_bins[i_pt + 1]]
            pt_centre = (self.pt_bins[i_pt] + self.pt_bins[i_pt + 1]) / 2
            pt_sel = f"abs(fPt) > {pt_bin[0]} and abs(fPt) < {pt_bin[1]}"
            if self.ptdep_selection_string:
                pt_sel = pt_sel + " and " + self.ptdep_selection_string[i_pt]
            cent_sel = f"fCentFT0C > {self.cent_limits[0]} and fCentFT0C < {self.cent_limits[1]}"
            if self.selection_string:
                bin_sel = f"{cent_sel} and {pt_sel} and {self.selection_string}"
            else:
                bin_sel = f"{cent_sel} and {pt_sel}"

            pt_label = (
                f"{pt_bin[0]:.2f} "
                + r"#leq #it{p}_{T} < "
                + f"{pt_bin[1]:.2f}"
                + r" GeV/#it{c}"
            )

            # select the correct pt bin
            bin_df = self.data_df.query(bin_sel, inplace=False)

            for nSigmaTPC, v2, m2 in zip(
                bin_df["fNsigmaTPC3He"],
                bin_df[f"fV2{self.ref_detector}"],
                bin_df["fTOFmassSquared"],
            ):
                self.hNsigma3He[i_pt].Fill(nSigmaTPC)
                self.hV2vsNsigma3He2D[i_pt].Fill(nSigmaTPC, v2)
                if self.tof_analysis:
                    self.hTOFmassSquared[i_pt].Fill(
                        m2 - utils.mass_alpha * utils.mass_alpha
                    )

            self._make_purity(
                self.hNsigma3He[i_pt],
                pt_label=pt_label,
                pt_centre=pt_centre,
                pt_thr=3.0,
            )

            # system infopanel
            info_panel = ROOT.TPaveText(0.6, 0.6, 0.8, 0.82, "NDC")
            info_panel.SetBorderSize(0)
            info_panel.SetFillStyle(0)
            info_panel.SetTextAlign(12)
            info_panel.SetTextFont(42)
            info_panel.AddText("ALICE")
            info_panel.AddText(r"PbPb, #sqrt{#it{s}_{NN}} = 5.36 TeV")
            info_panel.AddText(
                f"{self.cent_limits[0]} - {self.cent_limits[1]} % {self.cent_detector}"
            )
            info_panel.AddText(pt_label)

            canvas = utils.getCanvasWithTwoPanels(
                f"cV2{self.ref_detector}vsNsigma_cent_{self.cent_limits[0]}_{self.cent_limits[1]}_pt{i_pt}{self.suffix}",
                self.hV2vsNsigma3He[i_pt],
                self.hNsigma3He[i_pt],
                bottom_panel=info_panel,
            )

            self.cV2vsNsigma3He.append(canvas)

            df_bin_3He = bin_df.query(f"abs(fNsigmaTPC3He) < {self.n_sigma_selection}")
            if not self.tof_analysis:
                self.hRawCountsVsPt.SetBinContent(i_pt + 1, len(df_bin_3He))
                self.hRawCountsVsPt.SetBinError(i_pt + 1, np.sqrt(len(df_bin_3He)))
            for v2 in df_bin_3He[f"fV2{self.ref_detector}"]:
                self.hV2[i_pt].Fill(v2)

        # v2 vs Pt
        for i_pt in range(0, len(self.pt_bins) - 1):
            self.hV2vsPt.SetBinContent(i_pt + 1, self.hV2[i_pt].GetMean())
            self.hV2vsPt.SetBinError(i_pt + 1, self.hV2[i_pt].GetMeanError())

        self.hV2vsPt.Scale(1 / self.resolution)

        # purity vs pt
        for i_pt in range(0, len(self.pt_bins) - 1):
            self.hPurityVsPt.SetBinContent(i_pt + 1, self.purity[i_pt])

    def _make_purity(self, hist, pt_label, pt_centre, pt_thr):
        # fit n-sigma distribution
        sigma_rootfitter = roofitter.RooFitter()
        sigma_rootfitter.variable_range = self.nsigmaTPC_bin_limits
        sigma_rootfitter.hist = hist
        sigma_rootfitter.integral_range = [
            -1 * self.n_sigma_selection,
            self.n_sigma_selection,
        ]
        sigma_rootfitter.pt_label = pt_label
        sigma_rootfitter.cent_label = (
            f"{self.cent_limits[0]} - {self.cent_limits[1]} % {self.cent_detector}"
        )
        if pt_centre > pt_thr:
            sigma_rootfitter.no_bkg = True
        sigma_rootfitter.initialise()
        sigma_rootfitter.fit()
        if self.print_frame:
            sigma_rootfitter.saveFrameAsPDF(self.plot_dir)

        self.cNsigma3HeFit.append(sigma_rootfitter.canvas)

        self.purity.append(sigma_rootfitter.purity)

    def getFlowValues(self):
        return utils.getValuesFromHisto(self.hV2vsPt)

    def nPtBins(self):
        return self.n_pt_bins

    def dump_to_output_file(self):
        self.output_dir.cd()
        for i_pt in range(0, len(self.pt_bins) - 1):
            bin = [self.pt_bins[i_pt], self.pt_bins[i_pt + 1]]
            self.output_dir.mkdir(f"pt_{bin[0]}_{bin[1]}")
            self.output_dir.cd(f"pt_{bin[0]}_{bin[1]}")
            self.hNsigma3He[i_pt].Write()
            self.cNsigma3HeFit[i_pt].Write()
            if self.tof_analysis:
                self.hTOFmassSquared[i_pt].Write()
            self.hV2vsNsigma3He2D[i_pt].Write()
            self.hV2vsNsigma3He[i_pt].Write()
            self.cV2vsNsigma3He[i_pt].Write()
            self.hV2[i_pt].Write()
        self.output_dir.cd()
        self.hV2vsPt.Write()
        self.hPurityVsPt.Write()
        self.hRawCountsVsPt.Write()

    def dump_to_pdf(self):
        if not self.plot_dir:
            raise ValueError("plot_dir not set.")
        if not os.path.exists(self.plot_dir):
            os.makedirs(self.plot_dir)
        for i_pt in range(0, len(self.pt_bins) - 1):
            utils.saveCanvasAsPDF(self.hNsigma3He[i_pt], self.plot_dir)
            utils.saveCanvasAsPDF(self.hV2[i_pt], self.plot_dir)
        utils.saveCanvasAsPDF(self.hV2vsPt, self.plot_dir)
        utils.saveCanvasAsPDF(self.hPurityVsPt, self.plot_dir)
        utils.saveCanvasAsPDF(self.hRawCountsVsPt, self.plot_dir)

    def dump_summary_to_pdf(self, rows_per_page=3):
        if not self.plot_dir:
            raise ValueError("plot_dir not set.")
        if not os.path.exists(self.plot_dir):
            os.makedirs(self.plot_dir)

        # Define the name of the final PDF
        multipage_pdf = f"{self.plot_dir}/v2_and_nSigma_cent_{self.cent_limits[0]}_{self.cent_limits[1]}.pdf"

        # Create a multipage PDF object, first print command with '(' to initiate PDF
        canvas = ROOT.TCanvas("canvas", "", 595, 842)  # A4 size in portrait
        canvas.Print(f"{multipage_pdf}[", "pdf")

        total_pt_bins = len(self.pt_bins) - 1

        # Loop through pt bins in groups of rows_per_page
        for page_start in range(0, total_pt_bins, rows_per_page):
            # Create a new canvas for each page with A4 size
            canvas.Divide(
                2, rows_per_page
            )  # 2 columns (hNsigma3He and hV2), 3 rows per page

            # Add histograms to the canvas
            for i in range(rows_per_page):
                i_pt = page_start + i
                if i_pt >= total_pt_bins:
                    break  # Exit the loop if we've exceeded the available bins

                # Draw hNsigma3He in the (2*i+1) section of the canvas
                canvas.cd(2 * i + 1)
                self.hNsigma3He[i_pt].SetTitleOffset(1.5)  # Add margin to the histogram
                self.hNsigma3He[i_pt].GetXaxis().SetTitleOffset(
                    1.5
                )  # Add margin to x-axis
                self.hNsigma3He[i_pt].GetYaxis().SetTitleOffset(
                    1.5
                )  # Add margin to y-axis
                self.hNsigma3He[i_pt].Draw()
                # Set margins for the current panel
                ROOT.gPad.SetBottomMargin(
                    0.15
                )  # Set bottom margin for the current panel
                ROOT.gPad.SetLeftMargin(0.15)  # Set left margin for the current panel

                # Draw hV2 in the (2*i+2) section of the canvas
                canvas.cd(2 * i + 2)
                self.hV2[i_pt].SetTitleOffset(1.5)  # Add margin to the histogram
                self.hV2[i_pt].GetXaxis().SetTitleOffset(1.5)  # Add margin to x-axis
                self.hV2[i_pt].GetYaxis().SetTitleOffset(1.5)  # Add margin to y-axis
                self.hV2[i_pt].Draw()
                # Set margins for the current panel
                ROOT.gPad.SetBottomMargin(
                    0.15
                )  # Set bottom margin for the current panel
                ROOT.gPad.SetLeftMargin(0.15)  # Set left margin for the current panel

                # Add the info panel to hV2
                pt_bin = [self.pt_bins[i_pt], self.pt_bins[i_pt + 1]]
                pt_label = (
                    f"{pt_bin[0]:.2f} "
                    + r"#leq #it{p}_{T} < "
                    + f"{pt_bin[1]:.2f}"
                    + r" GeV/#it{c}"
                )
                info_panel_bis = ROOT.TPaveText(0.4, 0.6, 0.6, 0.82, "NDC")
                info_panel_bis.SetBorderSize(0)
                info_panel_bis.SetFillStyle(0)
                info_panel_bis.SetTextAlign(12)
                info_panel_bis.SetTextFont(42)
                info_panel_bis.AddText("ALICE")
                info_panel_bis.AddText(r"PbPb, #sqrt{#it{s}_{NN}} = 5.36 TeV")
                info_panel_bis.AddText(
                    f"{self.cent_limits[0]} - {self.cent_limits[1]} % {self.cent_detector}"
                )
                info_panel_bis.AddText(pt_label)

                self.hV2[i_pt].GetListOfFunctions().Add(info_panel_bis)

            # Add the page to the multipage PDF
            canvas.Print(f"{multipage_pdf}", "pdf")
            canvas.Clear()

        # Add remaining histograms (hV2vsPt, hPurityVsPt, hRawCountsVsPt) in one page
        canvas.Divide(1, 3)  # 1 column, 3 rows

        # Draw hV2vsPt
        canvas.cd(1)
        ROOT.gPad.SetBottomMargin(0.15)  # Set bottom margin for the current panel
        ROOT.gPad.SetLeftMargin(0.15)  # Set left margin for the current panel
        self.hV2vsPt.Draw("PE")
        self.hV2vsPt.GetXaxis().SetTitleOffset(1.5)  # Add margin to x-axis
        self.hV2vsPt.GetYaxis().SetTitleOffset(1.5)  # Add margin to y-axis
        info_panel_v2 = ROOT.TPaveText(0.2, 0.6, 0.4, 0.82, "NDC")
        info_panel_v2.SetBorderSize(0)
        info_panel_v2.SetFillStyle(0)
        info_panel_v2.SetTextAlign(12)
        info_panel_v2.SetTextFont(42)
        info_panel_v2.AddText("ALICE")
        info_panel_v2.AddText(r"PbPb, #sqrt{#it{s}_{NN}} = 5.36 TeV")
        info_panel_v2.AddText(
            f"{self.cent_limits[0]} - {self.cent_limits[1]} % {self.cent_detector}"
        )
        info_panel_v2.AddText("hV2vsPt")
        info_panel_v2.Draw()

        # Draw hPurityVsPt
        canvas.cd(2)
        ROOT.gPad.SetBottomMargin(0.15)  # Set bottom margin for the current panel
        ROOT.gPad.SetLeftMargin(0.15)  # Set left margin for the current panel
        self.hPurityVsPt.Draw("histo")
        self.hPurityVsPt.GetXaxis().SetTitleOffset(1.5)  # Add margin to x-axis
        self.hPurityVsPt.GetYaxis().SetTitleOffset(1.5)  # Add margin to y-axis
        info_panel_purity = ROOT.TPaveText(0.6, 0.5, 0.8, 0.72, "NDC")
        info_panel_purity.SetBorderSize(0)
        info_panel_purity.SetFillStyle(0)
        info_panel_purity.SetTextAlign(12)
        info_panel_purity.SetTextFont(42)
        info_panel_purity.AddText("ALICE")
        info_panel_purity.AddText(r"PbPb, #sqrt{#it{s}_{NN}} = 5.36 TeV")
        info_panel_purity.AddText(
            f"{self.cent_limits[0]} - {self.cent_limits[1]} % {self.cent_detector}"
        )
        info_panel_purity.AddText("hPurityVsPt")
        info_panel_purity.Draw()

        # Draw hRawCountsVsPt
        canvas.cd(3)
        ROOT.gPad.SetBottomMargin(0.15)  # Set bottom margin for the current panel
        ROOT.gPad.SetLeftMargin(0.15)  # Set left margin for the current panel
        self.hRawCountsVsPt.Draw("histo")
        self.hRawCountsVsPt.GetXaxis().SetTitleOffset(1.5)  # Add margin to x-axis
        self.hRawCountsVsPt.GetYaxis().SetTitleOffset(1.5)  # Add margin to y-axis
        info_panel_raw = ROOT.TPaveText(0.6, 0.6, 0.8, 0.82, "NDC")
        info_panel_raw.SetBorderSize(0)
        info_panel_raw.SetFillStyle(0)
        info_panel_raw.SetTextAlign(12)
        info_panel_raw.SetTextFont(42)
        info_panel_raw.AddText("ALICE")
        info_panel_raw.AddText(r"PbPb, #sqrt{#it{s}_{NN}} = 5.36 TeV")
        info_panel_raw.AddText(
            f"{self.cent_limits[0]} - {self.cent_limits[1]} % {self.cent_detector}"
        )
        info_panel_raw.AddText("hRawCountsVsPt")
        info_panel_raw.Draw()

        # Add the page to the PDF
        canvas.Print(f"{multipage_pdf}", "pdf")

        # Close the multipage PDF
        canvas.Print(f"{multipage_pdf}]", "pdf")

    def del_dyn_members(self):
        self.hNsigma3He = []
        self.cNsigma3HeFit = []
        self.hTOFmassSquared = []
        self.hV2vsNsigma3He2D = []
        self.hV2vsNsigma3He = []
        self.cV2vsNsigma3He = []
        self.hV2vsTOFmassSquared2D = []
        self.hV2vsTOFmassSquared = []
        self.cV2vsTOFmassSquared = []
        self.hV2 = []
        self.hV2vsPt = None
        self.hPurityVsP = None
        self.hRawCountsVsPt = None
        self.purity = []
