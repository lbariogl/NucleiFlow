import ROOT

import sys
sys.path.append('utils')
import utils as utils

hCentFT0C = []

input_files = []
input_files.append(ROOT.TFile('/data/lbariogl/flow/LHC23_PbPb_pass3/544124/AnalysisResults.root'))
input_files.append(ROOT.TFile('/data/lbariogl/flow/LHC23_PbPb_pass3/544123/AnalysisResults.root'))
input_files.append(ROOT.TFile('/data/lbariogl/flow/LHC23_PbPb_pass3/544475/AnalysisResults.root'))
input_files.append(ROOT.TFile('/data/lbariogl/flow/LHC23_PbPb_pass3/544032/AnalysisResults.root'))
input_files.append(ROOT.TFile('/data/lbariogl/flow/LHC23_PbPb_pass3/544477/AnalysisResults.root'))
input_files.append(ROOT.TFile('/data/lbariogl/flow/LHC23_PbPb_pass3/544474/AnalysisResults.root'))
input_files.append(ROOT.TFile('/data/lbariogl/flow/LHC23_PbPb_pass3/544491/AnalysisResults.root'))
input_files.append(ROOT.TFile('/data/lbariogl/flow/LHC23_PbPb_pass3/544492/AnalysisResults.root'))
input_files.append(ROOT.TFile('/data/lbariogl/flow/LHC23_PbPb_pass3/544091/AnalysisResults.root'))
input_files.append(ROOT.TFile('/data/lbariogl/flow/LHC23_PbPb_pass3/544184/AnalysisResults.root'))
input_files.append(ROOT.TFile('/data/lbariogl/flow/LHC23_PbPb_pass3/544390/AnalysisResults.root'))
input_files.append(ROOT.TFile('/data/lbariogl/flow/LHC23_PbPb_pass3/544451/AnalysisResults.root'))
input_files.append(ROOT.TFile('/data/lbariogl/flow/LHC23_PbPb_pass3/544122/AnalysisResults.root'))
input_files.append(ROOT.TFile('/data/lbariogl/flow/LHC23_PbPb_pass3/544116/AnalysisResults.root'))
input_files.append(ROOT.TFile('/data/lbariogl/flow/LHC23_PbPb_pass3/544095/AnalysisResults.root'))
input_files.append(ROOT.TFile('/data/lbariogl/flow/LHC23_PbPb_pass3/544098/AnalysisResults.root'))
input_files.append(ROOT.TFile('/data/lbariogl/flow/LHC23_PbPb_pass3/544490/AnalysisResults.root'))
input_files.append(ROOT.TFile('/data/lbariogl/flow/LHC23_PbPb_pass3/544013/AnalysisResults.root'))
input_files.append(ROOT.TFile('/data/lbariogl/flow/LHC23_PbPb_pass3/544391/AnalysisResults.root'))
input_files.append(ROOT.TFile('/data/lbariogl/flow/LHC23_PbPb_pass3/544389/AnalysisResults.root'))
input_files.append(ROOT.TFile('/data/lbariogl/flow/LHC23_PbPb_pass3/544392/AnalysisResults.root'))
input_files.append(ROOT.TFile('/data/lbariogl/flow/LHC23_PbPb_pass3/544121/AnalysisResults.root'))
input_files.append(ROOT.TFile('/data/lbariogl/flow/LHC23_PbPb_pass3/544510/AnalysisResults.root'))
input_files.append(ROOT.TFile('/data/lbariogl/flow/LHC23_PbPb_pass3/544476/AnalysisResults.root'))
input_files.append(ROOT.TFile('/data/lbariogl/flow/LHC23_PbPb_pass3/544508/AnalysisResults.root'))
input_files.append(ROOT.TFile('/data/lbariogl/flow/LHC23_PbPb_pass3/544028/AnalysisResults.root'))

run_number = [544124, 544123, 544475, 544032, 544477, 544474, 544491, 544492, 544091, 544184, 544390, 544451, 544122, 544116, 544095, 544098, 544490, 544013, 544391, 544389, 544392, 544121, 544510, 544476, 544508, 544028]

for run, input_file in zip(run_number, input_files):
  histo_tmp = input_file.Get('flow-qc/flow/hCentFT0C')
  histo_tmp.SetDirectory(0)
  histo_tmp.SetName(f'{histo_tmp.GetName()}_{run}')
  hCentFT0C.append(histo_tmp)

n_histos = len(hCentFT0C)
cols = ROOT.TColor.GetPalette()
period = int(cols.GetSize() / n_histos)

output_file = ROOT.TFile('../results_pass3/cent_check.root', 'recreate')
canvas = ROOT.TCanvas('cCentCheck', 'cCentCheck', 800, 600)
canvas.DrawFrame(0, 0, 100, 3.5e+6, r';FT0C percentile (%)')
legend = ROOT.TLegend(0.13, 0.75, 0.85, 0.87, 'Run numbers', 'brNDC')
legend.SetBorderSize(0)
legend.SetNColumns(9)

for i, histo in enumerate(hCentFT0C):
  utils.setHistStyle(histo, cols.At(i*period), linewidth=2)
  histo.Draw('same')
  legend.AddEntry(histo, f'{run_number[i]}','PE')

legend.Draw()
canvas.Write()
canvas.SaveAs('../results_pass3/qc_plots/cent_check.pdf')
