import ROOT

input_file = ROOT.TFile("/data/lbariogl/flow/MC/AnalysisResults_Ivan.root")

hMMomResHe3_3D = input_file.Get("nuclei-spectra/spectra/hMMomResHe3")
hAMomResHe3_3D = input_file.Get("nuclei-spectra/spectra/hAMomResHe3")

hResolutionMatter = hMMomResHe3_3D.Project3D("zy")
hResolutionMatter.SetName("hResolutionMatter")

hResolutionAntimatter = hAMomResHe3_3D.Project3D("zy")
hResolutionAntimatter.SetName("hResolutionAntimatter")

# Save the projected histogram into a new ROOT file
output_file = ROOT.TFile("../results_pass4/MC_check.root", "RECREATE")
hResolutionMatter.Write()
hResolutionAntimatter.Write()
