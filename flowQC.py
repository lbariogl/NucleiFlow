import ROOT
from hipe4ml.tree_handler import TreeHandler
import numpy as np
import pandas as pd
import struct

import sys
sys.path.append('utils')
import utils as utils
import copy
import os


input_file_name = '/data/lbariogl/flow/LHC23_PbPb_pass2_new/AnalysisResults.root'
input_file = ROOT.TFile(input_file_name)

input_dir = input_file.Get('flow-qc/flow')

plots_dir = '../results/ep_plots'

if not os.path.exists(plots_dir):
  os.makedirs(plots_dir)

hZvtx = input_dir.Get('hRecVtxZData')
utils.setHistStyle(hZvtx, ROOT.kRed+2)
utils.saveCanvasAsPDF(hZvtx, plots_dir)

hCentFT0C = input_dir.Get('hCentFT0C')
utils.setHistStyle(hCentFT0C, ROOT.kRed+2)
utils.saveCanvasAsPDF(hCentFT0C, plots_dir)

hDeltaPsi_FT0A_TPCl = input_dir.Get('hDeltaPsi_FT0A_TPCl')
utils.saveCanvasAsPDF(hDeltaPsi_FT0A_TPCl, plots_dir, is2D=True)

hDeltaPsi_FT0A_TPCr = input_dir.Get('hDeltaPsi_FT0A_TPCr')
utils.saveCanvasAsPDF(hDeltaPsi_FT0A_TPCr, plots_dir, is2D=True)

hDeltaPsi_FT0C_TPCl = input_dir.Get('hDeltaPsi_FT0C_TPCl')
utils.saveCanvasAsPDF(hDeltaPsi_FT0C_TPCl, plots_dir, is2D=True)

hDeltaPsi_FT0C_TPCr = input_dir.Get('hDeltaPsi_FT0C_TPCr')
utils.saveCanvasAsPDF(hDeltaPsi_FT0C_TPCr, plots_dir, is2D=True)

hDeltaPsi_FT0C_FT0A = input_dir.Get('hDeltaPsi_FT0C_FT0A')
utils.saveCanvasAsPDF(hDeltaPsi_FT0C_FT0A, plots_dir, is2D=True)

hDeltaPsi_TPCl_TPCr = input_dir.Get('hDeltaPsi_TPCl_TPCr')
utils.saveCanvasAsPDF(hDeltaPsi_TPCl_TPCr, plots_dir, is2D=True)

resolution_file = ROOT.TFile('Resolution_FT0C.root')
hResolutionFT0C = resolution_file.Get('Resolutuion')
hResolutionFT0C.SetDirectory(0)
hResolutionFT0C.SetName('hResolutionFT0C')
hResolutionFT0C.SetTitle(r';FT0C percentile (%); R_{2}')
utils.setHistStyle(hResolutionFT0C, ROOT.kRed+2)
cResolutionFT0C = ROOT.TCanvas('cResolutionFT0C', 'cResolutionFT0C', 800, 600)
cResolutionFT0C.SetBottomMargin(0.15)
hResolutionFT0C.Draw('PE')
cResolutionFT0C.SaveAs(f'{plots_dir}/cResolutionFT0C.pdf')

