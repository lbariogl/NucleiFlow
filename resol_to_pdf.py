import ROOT
import sys

def print_histograms_to_pdf(root_filename, pdf_filename):
    # Open the ROOT file
    root_file = ROOT.TFile.Open(root_filename)
    if not root_file or root_file.IsZombie():
        print(f"Error: Cannot open {root_filename}")
        return

    # Get list of keys in the file
    res_EP = root_file.Get("Resolution_EP")
    res_SP = root_file.Get("Resolution_SP")
    # Prepare canvas
    canvas = ROOT.TCanvas("canvas", "Canvas", 800, 600)
    canvas.Print(f"{pdf_filename}[")

    # Loop over all objects in the file
    res_EP.cd()
    for obj in res_EP.GetListOfKeys():
        # Check if the object is a histogram
        obj = res_EP.Get(obj.GetName())
        if obj.InheritsFrom("TH1"):
            info_panel = ROOT.TPaveText(0.6, 0.85, 0.8, 0.65, "NDC")
            info_panel.SetBorderSize(0)
            info_panel.SetFillStyle(0)
            info_panel.SetTextAlign(12)
            info_panel.SetTextFont(42)
            info_panel.AddText("EP Method")

            obj.GetListOfFunctions().Add(info_panel)
            obj.Draw()
            canvas.Print(pdf_filename)
    res_SP.cd()
    for obj in res_SP.GetListOfKeys():
        # Check if the object is a histogram
        obj = res_SP.Get(obj.GetName())
        if obj and obj.InheritsFrom("TH1"):
            info_panel = ROOT.TPaveText(0.6, 0.85, 0.8, 0.65, "NDC")
            info_panel.SetBorderSize(0)
            info_panel.SetFillStyle(0)
            info_panel.SetTextAlign(12)
            info_panel.SetTextFont(42)
            info_panel.AddText("SP Method")

            obj.GetListOfFunctions().Add(info_panel)
            obj.Draw()
            canvas.Print(pdf_filename)

    canvas.Print(f"{pdf_filename}]")
    root_file.Close()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py input.root output.pdf")
        sys.exit(1)
    print_histograms_to_pdf(sys.argv[1], sys.argv[2])