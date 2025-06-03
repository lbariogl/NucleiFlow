import ROOT
import sys

import utils

def combine_root_files(file1_path, file2_path, output_path, harmonic, do_print_pt=True):
    file1 = ROOT.TFile.Open(file1_path, "READ")
    file2 = ROOT.TFile.Open(file2_path, "READ")
    
    output_file = ROOT.TFile(output_path, "RECREATE")
    pdf_output = output_path.replace(".root", ".pdf")
    
    pdf_canvas = ROOT.TCanvas("pdf_canvas", "pdf_canvas", 800, 800)
    pdf_canvas.Print(f"{pdf_output}[")

    centrality = ['0-20%', '20-40%', '40-60%']
    i = 0
    
    for cent_key in file1.GetListOfKeys():
        if cent_key.GetClassName() == "TDirectoryFile":
            cent_name = cent_key.GetName()
            
            cent_dir = output_file.mkdir(cent_name)
            cent_dir.cd()
            
            cent1 = file1.Get(cent_name)
            cent2 = file2.Get(cent_name)
            
            for default_key in cent1.GetListOfKeys():
                if default_key.GetClassName() == "TDirectoryFile":
                    default_name = default_key.GetName()
                    
                    default_dir = cent_dir.mkdir(default_name)
                    default_dir.cd()
                    
                    default1 = cent1.Get(default_name)
                    default2 = cent2.Get(default_name)
                    
                    for pt_key in default1.GetListOfKeys():
                        if pt_key.GetClassName() == "TDirectoryFile":
                            if do_print_pt: continue
                            pt_name = pt_key.GetName()
                            
                            pt_dir = default_dir.mkdir(pt_name)
                            pt_dir.cd()
                            
                            pt1 = default1.Get(pt_name)
                            pt2 = default2.Get(pt_name)
                            
                            for obj_key in pt1.GetListOfKeys():
                                obj_name = obj_key.GetName()
                                print(obj_name)
                                obj_class = obj_key.GetClassName()
                                print(obj_class)
                                
                                obj1 = pt1.Get(obj_name)
                                obj2 = pt2.Get(obj_name)
                                
                                if obj_class.startswith("TH1"):  
                                    obj1.Scale(1. / obj1.Integral())
                                    obj2.Scale(1. / obj2.Integral())
                                    
                                    max_val = max(obj1.GetMaximum(), obj2.GetMaximum()) * 1.1
                                    
                                    # utils.getCanvasWithTwoPanels(f"c_{obj_name}", obj1, obj2)

                                    canvas = ROOT.TCanvas(f"c_{obj_name}", obj_name, 800, 600)
                                    pad1 = ROOT.TPad("pad1", "pad1", 0.0, 0.5, 1.0, 1.0, 0)
                                    pad2 = ROOT.TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.5, 0)
                                    pad1.SetBottomMargin(0.0)
                                    pad1.SetLeftMargin(0.15)
                                    pad2.SetTopMargin(0.0)
                                    pad2.SetBottomMargin(0.3)
                                    pad2.SetLeftMargin(0.15)
                                    pad1.Draw()
                                    pad2.Draw()
                                    
                                    pad1.cd()
                                    frame = pad1.DrawFrame(obj1.GetXaxis().GetXmin(), 0, obj1.GetXaxis().GetXmax(), max_val)
                                    obj1.SetLineColor(ROOT.kRed)
                                    obj1.SetMarkerColor(ROOT.kRed)
                                    obj2.SetLineColor(ROOT.kBlue)
                                    obj2.SetMarkerColor(ROOT.kBlue)
                                    obj1.Draw("HIST SAME")
                                    obj2.Draw("HIST SAME")
                                    
                                    pad2.cd()
                                    # pad2.SetGridy()
                                    diff = obj1.Clone(f"diff_{obj_name}")
                                    diff.Add(obj2,-1)
                                    diff.SetLineColor(ROOT.kBlack)
                                    diff.SetMarkerColor(ROOT.kBlack)
                                    diff.SetMinimum(-2)
                                    diff.SetMaximum(2)
                                    # diff.GetXaxis().SetLabelSize(0.1)
                                    # diff.GetYaxis().SetLabelSize(0.1)
                                    diff.Draw("PE")
                                    
                                    line = ROOT.TLine(obj1.GetXaxis().GetXmin(), 1, obj1.GetXaxis().GetXmax(), 1)
                                    line.SetLineStyle(2)
                                    line.Draw("SAME")
                                    
                                    pad2.Modified()
                                    pad2.Update()

                                    canvas.Modified()
                                    canvas.Update()
                                    
                                    canvas.Write()
                                    pdf_canvas.cd()
                                    canvas.Print(pdf_output)
                        elif "V%s" % harmonic in pt_key.GetName():
                            pt_name = pt_key.GetName()
                            print(pt_name)
                            pt_class = pt_key.GetClassName()
                            print(pt_class)
                            v2_hist1 = default1.Get(pt_name)
                            v2_hist2 = default2.Get(pt_name)

                            if not v2_hist1 or not v2_hist2:
                                print(f"Skipping {pt_name} due to missing histograms.")
                                continue

                            if v2_hist1.Integral() == 0 or v2_hist2.Integral() == 0:
                                print(f"Skipping {pt_name} due to zero integral.")
                                continue

                            max_val = max(v2_hist1.GetMaximum(), v2_hist2.GetMaximum()) * 1.1

                            canvas = ROOT.TCanvas("v2_canvas", "v2_comparison", 800, 800)
                            pad1 = ROOT.TPad("pad1", "pad1", 0, 0.4, 1, 1)
                            pad2 = ROOT.TPad("pad2", "pad2", 0, 0, 1, 0.4)
                            pad1.SetBottomMargin(0.0)
                            pad1.SetLeftMargin(0.15)
                            pad2.SetTopMargin(0.0)
                            pad2.SetBottomMargin(0.3)
                            pad2.SetLeftMargin(0.15)
                            pad1.Draw()
                            pad2.Draw()

                            pad1.cd()
                            frame = pad1.DrawFrame(v2_hist1.GetXaxis().GetXmin(), 0, v2_hist1.GetXaxis().GetXmax(), max_val)
                            frame.SetTitle(f"Comparison of {pt_name}")
                            frame.GetXaxis().SetTitle(v2_hist1.GetXaxis().GetTitle())
                            frame.GetYaxis().SetTitle(v2_hist1.GetYaxis().GetTitle())
                            v2_hist1.SetLineColor(ROOT.kRed)
                            v2_hist1.SetMarkerColor(ROOT.kRed)
                            v2_hist2.SetLineColor(ROOT.kBlue)
                            v2_hist2.SetMarkerColor(ROOT.kBlue)
                            v2_hist1.Draw("PE")
                            v2_hist2.Draw("SAME")

                            info_panel = ROOT.TPaveText(0.40, 0.15, 0.60, 0.0, "NDC")
                            info_panel.SetBorderSize(0)
                            info_panel.SetFillStyle(0)
                            info_panel.SetTextAlign(22)
                            info_panel.SetTextFont(42)
                            info_panel.AddText("Centrality: %s" % centrality[i])
                            i=i+1

                            v2_hist1.GetListOfFunctions().Add(info_panel)

                            legend = ROOT.TLegend(0.40, 0.12, 0.60, 0.20, "", "brNDC")
                            legend.SetBorderSize(0)
                            legend.AddEntry(v2_hist1, r"v_{%s} with EP method" % harmonic, "LPE")
                            legend.AddEntry(v2_hist2, r"v_{%s} with SP method" % harmonic, "LPE")
                            legend.Draw()

                            pad2.cd()
                            # pad2.SetGridy()
                            diff = v2_hist1.Clone(f"diff_{pt_name}")
                            diff.GetListOfFunctions().RemoveLast()
                            diff.Add(v2_hist2,-1)
                            diff.SetLineColor(ROOT.kBlack)
                            diff.SetMarkerColor(ROOT.kBlack)
                            diff_min = abs(diff.GetMinimum()) + diff.GetBinError(diff.GetMinimumBin())
                            diff_max = abs(diff.GetMaximum()) + diff.GetBinError(diff.GetMaximumBin())
                            diff.SetMinimum(-1.1*abs(diff_min))
                            diff.SetMaximum(1.1*abs(diff_max))
                            diff.GetYaxis().SetTitle(r"v_{%s}^{EP} - v_{%s}^{SP}" % (harmonic, harmonic))
                            diff.Draw("PE")

                            line = ROOT.TLine(v2_hist1.GetXaxis().GetXmin(), 0, v2_hist1.GetXaxis().GetXmax(), 0)
                            line.SetLineStyle(2)
                            line.Draw("SAME")

                            pad2.Modified()
                            pad2.Update()

                            canvas.Modified()
                            canvas.Update()
                            
                            # Cambia la directory corrente a "default_dir" prima di salvare
                            default_dir.cd()
                            canvas.Write()
                            pdf_canvas.cd()
                            canvas.Print(pdf_output)
                                    
    pdf_canvas.Print(f"{pdf_output}]")
    output_file.Close()
    file1.Close()
    file2.Close()

def main():
    if len(sys.argv) < 4:
        print("Usage: python combine_root_files.py file1.root file2.root output.root")
        sys.exit(1)
    
    combine_root_files(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5] if len(sys.argv) > 5 else True)

if __name__ == "__main__":
    main()
