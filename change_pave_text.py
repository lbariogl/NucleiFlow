import ROOT
import sys
import os

# Apri il file ROOT
harmonic = sys.argv[2]
file_in = ROOT.TFile.Open(sys.argv[1], "READ")
file_out = ROOT.TFile("output.root", "RECREATE")

for cent_key in file_in.GetListOfKeys():
    if cent_key.GetClassName() == "TDirectoryFile":
        cent_name = cent_key.GetName()

        cent_dir = file_out.mkdir(cent_name)
        cent_dir.cd()

        cent = file_in.Get(cent_name)
            
        for default_key in cent.GetListOfKeys():
            if default_key.GetClassName() == "TDirectoryFile":
                default_name = default_key.GetName()

                default_dir = cent_dir.mkdir(default_name)
                default_dir.cd()

                default = cent.Get(default_name)

                for pt_key in default.GetListOfKeys():
                    if pt_key.GetClassName() == "TDirectoryFile":
                        pt_name = pt_key.GetName()

                        pt_dir = default_dir.mkdir(pt_name)
                        pt_dir.cd()

                        pt = default.Get(pt_name)

                        for obj_key in pt.GetListOfKeys():
                            
                            c = ROOT.TCanvas("c", "c", 1400, 900)
                            out_dir = f"png/{cent_name}/default/{pt_name}"
                            if not os.path.exists(out_dir):
                                os.makedirs(out_dir)

                            if "hV" in obj_key.GetName():
                                histo = pt.Get(obj_key.GetName())
                                
                                for obj in histo.GetListOfFunctions():
                                    if obj.InheritsFrom("TPaveText"):
                                        obj.GetListOfLines().RemoveLast()
                                        obj.AddText("v_%s = %.4f \pm %.4f" % (harmonic, histo.GetMean(), histo.GetMeanError()))
                                histo.Draw()
                                c.Modified()
                                c.Update()
                                c.SaveAs(f"{out_dir}/{obj_key.GetName()}.png")
                                histo.Write()
                            else:
                                histo = pt.Get(obj_key.GetName())
                                if histo.InheritsFrom("TH1"):
                                    histo.Draw()
                                    c.Modified()
                                    c.Update()
                                    c.SaveAs(f"{out_dir}/{obj_key.GetName()}.png")
                                else:
                                    histo.SaveAs(f"{out_dir}/{obj_key.GetName()}.png")
                                histo.Write()
                    else:

                        pt_name = pt_key.GetName()
                        pt = default.Get(pt_name)

                        c = ROOT.TCanvas("c", "c", 1400, 900)
                        out_dir = f"png/{cent_name}/default"
                        if not os.path.exists(out_dir):
                            os.makedirs(out_dir)

                        default_dir.cd()
                        if pt.InheritsFrom("TH1"):
                            pt.Draw()

                            c.Modified()
                            c.Update()
                            c.SaveAs(f"{out_dir}/{pt_name}.png")
                        else:
                            pt.SaveAs(f"{out_dir}/{pt_name}.png")

                        pt.Write()
    else:
        cent_name = cent_key.GetName()
        cent = file_in.Get(cent_name)

        c = ROOT.TCanvas("c", "c", 1400, 900)
        out_dir = f"png"
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        file_out.cd()
        if cent.InheritsFrom("TH1"):
            cent.Draw()

            c.Modified()
            c.Update()
            c.SaveAs(f"{out_dir}/{cent_name}.png")
        else:
            cent.SaveAs(f"{out_dir}/{cent_name}.png")
        
        cent.Write()

file_out.Close()
file_in.Close()