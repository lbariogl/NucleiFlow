import ROOT
import sys

sys.path.append("utils")
import utils as utils

input_file_name = "../phi_resolution/file_resolution.root"
input_file = ROOT.TFile(input_file_name, "READ")

hPhiResolution3D = input_file.Get("efficiency-q-a/phiRes")
if not hPhiResolution3D:  # Check if the histogram was loaded successfully
    raise RuntimeError(
        f"Histogram 'efficiency-q-a/phiRes' not found in {input_file_name}"
    )

hPhiResolution3D.SetDirectory(0)
input_file.Close()  # Close the input file

# Project hPhiResolution3D onto YZ for bin 5 of the X axis
bin_x = 5  # ROOT bins are 1-based

# Get axis info
y_axis = hPhiResolution3D.GetYaxis()
z_axis = hPhiResolution3D.GetZaxis()
ny = y_axis.GetNbins()
nz = z_axis.GetNbins()

# Create a new TH2D for the projection
hPhiResolution = ROOT.TH2D(
    f"hPhiResolution_bin{bin_x}",
    ";%s;%s;" % (y_axis.GetTitle(), z_axis.GetTitle()),
    ny,
    y_axis.GetXmin(),
    y_axis.GetXmax(),
    nz,
    z_axis.GetXmin(),
    z_axis.GetXmax(),
)

# Fill the new TH2D with selection on axis ranges
for iy in range(1, ny + 1):
    y_center = y_axis.GetBinCenter(iy)
    if -0.5 < y_center < 0.5:
        continue  # Skip bins between -0.5 and 0.5 on the x axis (y_axis here)
    for iz in range(1, nz + 1):
        z_center = z_axis.GetBinCenter(iz)
        if not (-ROOT.TMath.Pi() / 2 <= z_center <= ROOT.TMath.Pi() / 2):
            continue  # Skip bins outside [-pi/2, pi/2] on the z axis (z_axis here)
        content = hPhiResolution3D.GetBinContent(bin_x, iy, iz)
        error = hPhiResolution3D.GetBinError(bin_x, iy, iz)
        hPhiResolution.SetBinContent(iy, iz, content)
        hPhiResolution.SetBinContent(iy, iz, error)


hPhiResolution.SetName("hPhiResolution")

hPhiResolution.SetTitle(r";#it{p}_{T}; #phi^{rec} - #phi^{MC} (rad);")

output_file_name = "../results_phi_studies/phi_resolution.root"
output_file = ROOT.TFile(output_file_name, "RECREATE")
hPhiResolution.Write()

# Create a TProfile along the X axis (mean and RMS of Y for each X bin)
profile = hPhiResolution.ProfileX("profile_phi_resolution")
profile.SetLineColor(ROOT.kRed)
profile.SetLineWidth(2)

# Create a TH1D to store the RMS (width) vs pT
n_pt_bins = profile.GetNbinsX()
pt_axis = profile.GetXaxis()
hRMSVsPt = ROOT.TH1D(
    "hRMSVsPt",
    ";#it{p}_{T}; RMS",
    n_pt_bins,
    pt_axis.GetXmin(),
    pt_axis.GetXmax(),
)

for ix in range(1, n_pt_bins + 1):
    rms = profile.GetBinError(ix)  # This is the RMS for the bin
    n = profile.GetBinEntries(ix)
    hRMSVsPt.SetBinContent(ix, rms)
    if n > 1:
        rms_error = rms / (2 * n) ** 0.5  # Error on RMS: RMS / sqrt(2N)
    else:
        rms_error = 0
    hRMSVsPt.SetBinError(ix, rms_error)

# Write the profile and RMS histogram to file
profile.Write()
hRMSVsPt.Write()

# Draw both the 2D histogram and the profile on a canvas
canvas = ROOT.TCanvas("canvas", "Phi Resolution 2D and Profile", 900, 700)
hPhiResolution.Draw("COLZ")
profile.Draw("SAME")

# Save the canvas to file
canvas.SaveAs("../results_phi_studies/phi_resolution_with_profile.pdf")
canvas.Write()

# Fit the RMS vs pT with a second order polynomial (pol2), using only bins with pt > 0
fit_func = ROOT.TF1("fit_pol2", "pol2", 0.5, pt_axis.GetXmax())

# Find first bin with pt > 0
first_bin = 1
while first_bin <= n_pt_bins and pt_axis.GetBinCenter(first_bin) <= 0.5:
    first_bin += 1

# Only fit if there are bins with pt > 0
if first_bin <= n_pt_bins:
    hRMSVsPt.Fit(fit_func, "R", "", pt_axis.GetBinCenter(first_bin), pt_axis.GetXmax())

# Draw the fit on a new canvas
canvas_fit = ROOT.TCanvas("canvas_fit", "RMS vs pT with Fit", 800, 600)
hRMSVsPt.Draw("PE")
fit_func.SetLineColor(ROOT.kRed)
fit_func.Draw("L SAME")
canvas_fit.SaveAs("../results_phi_studies/phi_resolution_rms_fit.pdf")
canvas_fit.Write()
