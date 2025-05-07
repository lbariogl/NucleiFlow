from itertools import combinations
import yaml
import argparse
import numpy as np
import ROOT

import sys

sys.path.append("utils")
import utils as utils

parser = argparse.ArgumentParser(description="Configure the parameters of the script.")
parser.add_argument(
    "--config",
    dest="config_file",
    help="path to the YAML file with configuration.",
    default="",
)
args = parser.parse_args()
if args.config_file == "":
    print("** No config file provided. Exiting. **")
    exit()

config_file = open(args.config_file, "r")
config = yaml.full_load(config_file)

output_dir = config["output_dir_name"]

input_file_AR_name = config["input_file_AR_name"]
input_file = ROOT.TFile(input_file_AR_name)

ref_names = ["FT0C", "FT0A", "TPCl", "TPCr", "TPC"]

use_EP_tables = True

hSP_dict_EP = {}
hProfile_dict_EP = {}
cSpProfile_dict_EP = {}
hResolution_dict_EP = {}

hSP_dict_SP = {}
hProfile_dict_SP = {}
cSpProfile_dict_SP = {}
hResolution_dict_SP = {}

cent_bins = np.array(
    [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0], dtype=np.float64
)
n_cent_bins = len(cent_bins) - 1


def getResolution(
    histo_resolution, det1_det2_name, det2_det3_name, det1_det3_name, hProfile_dict
):
    histo_resolution.GetYaxis().SetRangeUser(0.0, 1.5)
    n_bins = histo_resolution.GetNbinsX()
    for ibin in range(1, n_bins + 1):
        val1_2 = hProfile_dict[det1_det2_name].GetBinContent(ibin)
        val2_3 = hProfile_dict[det2_det3_name].GetBinContent(ibin)
        val1_3 = hProfile_dict[det1_det3_name].GetBinContent(ibin)

        err1_2 = hProfile_dict[det1_det2_name].GetBinError(ibin)
        err2_3 = hProfile_dict[det2_det3_name].GetBinError(ibin)
        err1_3 = hProfile_dict[det1_det3_name].GetBinError(ibin)

        if val2_3 > 0:
            val = val1_2 * val1_3 / val2_3
        else:
            val = -999

        if val2_3 > 0 and val1_2 > 0 and val2_3 > 0:
            err = (
                np.sqrt(
                    err1_2 * err1_2 / val1_2 / val1_2
                    + err2_3 * err2_3 / val2_3 / val2_3
                    + err1_3 * err1_3 / val1_3 / val1_3
                )
                * val
            )
        else:
            err = 0
        if val > 0:
            histo_resolution.SetBinContent(ibin, np.sqrt(val))
            histo_resolution.SetBinError(ibin, np.sqrt(err))


def doAllPlots(
    det_name1,
    det_name2,
    det_name3,
    ep_table=True,
    use_SP=False,
    hSP_dict=hSP_dict_EP,
    hProfile_dict=hProfile_dict_EP,
    cSpProfile_dict=cSpProfile_dict_EP,
    hResolution_dict=hResolution_dict_EP,
):

    if ep_table:
        print("** Using EP tables **")
        directory = "flow_ep"
        input_suffix = "EP"
    else:
        print("** Using Qvector tables **")
        directory = "flow_qvec"
        input_suffix = "Qvec"

    if use_SP:
        print("** Using scalar product **")
        histo_name = "hScalarProduct"
        output_suffix = "SP"
    else:
        print("** Using event plane **")
        histo_name = "hNormalisedScalarProduct"
        output_suffix = "EP"

    # first combination
    det1_det2_name = f"{det_name1}_{det_name2}"
    if det1_det2_name not in hSP_dict.keys():
        hSP_dict[det1_det2_name] = input_file.Get(
            f"flow-qc/{directory}/{histo_name}_{det1_det2_name}_{input_suffix}"
        )
        hSP_dict[det1_det2_name].RebinX(10)

        hProfile_dict[det1_det2_name] = hSP_dict[det1_det2_name].ProfileX(
            f"hProfile_{det1_det2_name}_{output_suffix}"
        )
        utils.setHistStyle(hProfile_dict[det1_det2_name], ROOT.kRed)

        cSpProfile_dict[det1_det2_name] = ROOT.TCanvas(
            f"cSpProfile_{det1_det2_name}_{output_suffix}",
            f"cSpProfile_{det1_det2_name}_{output_suffix}",
            800,
            600,
        )
        hSP_dict[det1_det2_name].Draw("colz")
        hProfile_dict[det1_det2_name].Draw("pe same")

    # second combination
    det2_det3_name = f"{det_name2}_{det_name3}"
    if det2_det3_name not in hSP_dict.keys():
        hSP_dict[det2_det3_name] = input_file.Get(
            f"flow-qc/{directory}/{histo_name}_{det2_det3_name}_{input_suffix}"
        )
        hSP_dict[det2_det3_name].RebinX(10)

        hProfile_dict[det2_det3_name] = hSP_dict[det2_det3_name].ProfileX(
            f"hProfile_{det2_det3_name}_{output_suffix}"
        )
        utils.setHistStyle(hProfile_dict[det2_det3_name], ROOT.kRed)

        cSpProfile_dict[det2_det3_name] = ROOT.TCanvas(
            f"cSpProfile_{det2_det3_name}_{output_suffix}",
            f"cSpProfile_{det2_det3_name}_{output_suffix}",
            800,
            600,
        )
        hSP_dict[det2_det3_name].Draw("colz")
        hProfile_dict[det2_det3_name].Draw("pe same")

    # third combination
    det1_det3_name = f"{det_name1}_{det_name3}"
    if det1_det3_name not in hSP_dict.keys():
        hSP_dict[det1_det3_name] = input_file.Get(
            f"flow-qc/{directory}/{histo_name}_{det1_det3_name}_{input_suffix}"
        )
        hSP_dict[det1_det3_name].RebinX(10)

        hProfile_dict[det1_det3_name] = hSP_dict[det1_det3_name].ProfileX(
            f"hProfile_{det1_det3_name}_{output_suffix}"
        )
        utils.setHistStyle(hProfile_dict[det1_det3_name], ROOT.kRed)

        cSpProfile_dict[det1_det3_name] = ROOT.TCanvas(
            f"cSpProfile_{det1_det3_name}_{output_suffix}",
            f"cSpProfile_{det1_det3_name}_{output_suffix}",
            800,
            600,
        )
        hSP_dict[det1_det3_name].Draw("colz")
        hProfile_dict[det1_det3_name].Draw("pe same")

    # evaluate resolutions
    n_bins = hSP_dict[det1_det2_name].GetNbinsX()
    cent_axis_limits = [
        hSP_dict[det1_det2_name].GetXaxis().GetBinLowEdge(1),
        hSP_dict[det1_det2_name].GetXaxis().GetBinUpEdge(n_bins),
    ]
    cent_axis_title = hSP_dict[det1_det2_name].GetXaxis().GetTitle()

    resolution_name1 = f"{det_name1}_{det_name2}_{det_name3}"
    resolution_title1 = r"R_{2} " + f"({det_name1}#; {det_name2}, {det_name3})"

    print(f"Creating resolution {resolution_name1}")

    hResolution_dict[resolution_name1] = ROOT.TH1F(
        f"hResolution_{resolution_name1}_{output_suffix}",
        f";{cent_axis_title};{resolution_title1}",
        n_bins,
        cent_axis_limits[0],
        cent_axis_limits[1],
    )
    getResolution(
        hResolution_dict[resolution_name1],
        det1_det2_name,
        det2_det3_name,
        det1_det3_name,
        hProfile_dict,
    )
    utils.setHistStyle(hResolution_dict[resolution_name1], ROOT.kRed)

    resolution_name2 = f"{det_name2}_{det_name3}_{det_name1}"
    resolution_title2 = r"R_{2} " + f"({det_name2}#; {det_name1}, {det_name3})"

    print(f"Creating resolution {resolution_name2}")

    hResolution_dict[resolution_name2] = ROOT.TH1F(
        f"hResolution_{resolution_name2}_{output_suffix}",
        f";{cent_axis_title};{resolution_title2}",
        n_bins,
        cent_axis_limits[0],
        cent_axis_limits[1],
    )
    getResolution(
        hResolution_dict[resolution_name2],
        det2_det3_name,
        det1_det3_name,
        det1_det2_name,
        hProfile_dict,
    )
    utils.setHistStyle(hResolution_dict[resolution_name2], ROOT.kRed)

    resolution_name3 = f"{det_name3}_{det_name1}_{det_name2}"
    resolution_title3 = r"R_{2} " + f"({det_name3}#; {det_name1}, {det_name2})"

    print(f"Creating resolution {resolution_name3}")

    hResolution_dict[resolution_name3] = ROOT.TH1F(
        f"hResolution_{resolution_name3}_{output_suffix}",
        f";{cent_axis_title};{resolution_title3}",
        n_bins,
        cent_axis_limits[0],
        cent_axis_limits[1],
    )
    getResolution(
        hResolution_dict[resolution_name3],
        det1_det3_name,
        det1_det2_name,
        det2_det3_name,
        hProfile_dict,
    )
    utils.setHistStyle(hResolution_dict[resolution_name3], ROOT.kRed)


# rebin histograms


det_combos = list(combinations(ref_names, 3))

for det_combo in det_combos:
    doAllPlots(
        det_combo[0],
        det_combo[1],
        det_combo[2],
        use_EP_tables,
        use_SP=False,
        hSP_dict=hSP_dict_EP,
        hProfile_dict=hProfile_dict_EP,
        cSpProfile_dict=cSpProfile_dict_EP,
        hResolution_dict=hResolution_dict_EP,
    )
    doAllPlots(
        det_combo[0],
        det_combo[1],
        det_combo[2],
        use_EP_tables,
        use_SP=True,
        hSP_dict=hSP_dict_SP,
        hProfile_dict=hProfile_dict_SP,
        cSpProfile_dict=cSpProfile_dict_SP,
        hResolution_dict=hResolution_dict_SP,
    )

if use_EP_tables:
    output_file_name = output_dir + "resolution_EP.root"
else:
    output_file_name = output_dir + "resolution_Qvec.root"

print(f"Writing output to {output_file_name}")

output_file = ROOT.TFile(output_file_name, "recreate")
SP_dir_EP = output_file.mkdir("SP_normalised")
SP_dir_SP = output_file.mkdir("SP")
Resolution_EP_dir = output_file.mkdir("Resolution_EP")
Resolution_SP_dir = output_file.mkdir("Resolution_SP")

SP_dir_EP.cd()
for key in hSP_dict_EP.keys():
    hSP_dict_EP[key].Write()
    hProfile_dict_EP[key].Write()
    cSpProfile_dict_EP[key].Write()

SP_dir_SP.cd()
for key in hSP_dict_SP.keys():
    hSP_dict_SP[key].Write()
    hProfile_dict_SP[key].Write()
    cSpProfile_dict_SP[key].Write()

Resolution_EP_dir.cd()
for res in hResolution_dict_EP.values():
    res.Write()

Resolution_SP_dir.cd()
for res in hResolution_dict_SP.values():
    res.Write()

# check sourav
input_file_sourav = ROOT.TFile("../sourav_check/AnalysisResults.root")

print("Sourav check")
hSP_dict_Sourav = {}
hProfile_dict_Sourav = {}
cSpProfile_dict_Sourav = {}
hResolution_dict_Sourav = {}

hSP_FT0C_TPC = input_file_sourav.Get("phipbpb/ResSPFT0CTPC").Project3D("zx")
hSP_FT0C_TPC.SetName("hScalarProduct_FT0C_TPC_Sourav")
hSP_FT0C_TPC.SetTitle(r";centrality (%); #vec{Q}_{2}^{FT0C} #upoint #vec{Q}_{2}^{TPC}")
# hSP_FT0C_TPC.GetYaxis().SetRangeUser(-2.0, 2.0)
hSP_dict_Sourav["FT0C_TPC"] = hSP_FT0C_TPC
hProfile_dict_Sourav["FT0C_TPC"] = hSP_FT0C_TPC.ProfileX("hProfile_FT0C_TPC_Sourav")
utils.setHistStyle(hProfile_dict_Sourav["FT0C_TPC"], ROOT.kRed)
cSpProfile_dict_Sourav["FT0C_TPC"] = ROOT.TCanvas(
    "cSpProfile_FT0C_TPC_Sourav",
    "cSpProfile_FT0C_TPC_Sourav",
    800,
    600,
)
hSP_FT0C_TPC.Draw("colz")
hProfile_dict_Sourav["FT0C_TPC"].Draw("pe same")

hSP_FT0C_FT0A = input_file_sourav.Get("phipbpb/ResSPFT0CFT0A").Project3D("zx")
hSP_FT0C_FT0A.SetName("hScalarProduct_FT0C_FT0A_Sourav")
hSP_FT0C_FT0A.SetTitle(
    r";centrality (%); #vec{Q}_{2}^{FT0C} #upoint #vec{Q}_{2}^{FT0A}"
)
# hSP_FT0C_FT0A.GetYaxis().SetRangeUser(-2.0, 2.0)
hSP_dict_Sourav["FT0C_FT0A"] = hSP_FT0C_FT0A
hProfile_dict_Sourav["FT0C_FT0A"] = hSP_FT0C_FT0A.ProfileX("hProfile_FT0C_FT0A_Sourav")
utils.setHistStyle(hProfile_dict_Sourav["FT0C_FT0A"], ROOT.kRed)
cSpProfile_dict_Sourav["FT0C_FT0A"] = ROOT.TCanvas(
    "cSpProfile_FT0C_FT0A_Sourav",
    "cSpProfile_FT0C_FT0A_Sourav",
    800,
    600,
)
hSP_FT0C_FT0A.Draw("colz")
hProfile_dict_Sourav["FT0C_FT0A"].Draw("pe same")


hSP_FT0A_TPC = input_file_sourav.Get("phipbpb/ResSPFT0ATPC").Project3D("zx")
hSP_FT0A_TPC.SetName("hScalarProduct_FT0A_TPC_Sourav")
hSP_FT0A_TPC.SetTitle(r";centrality (%); #vec{Q}_{2}^{FT0A} #upoint #vec{Q}_{2}^{TPC}")
# hSP_FT0A_TPC.GetYaxis().SetRangeUser(-2.0, 2.0)
hSP_dict_Sourav["FT0A_TPC"] = hSP_FT0A_TPC
hProfile_dict_Sourav["FT0A_TPC"] = hSP_FT0A_TPC.ProfileX("hProfile_FT0A_TPC_Sourav")
utils.setHistStyle(hProfile_dict_Sourav["FT0A_TPC"], ROOT.kRed)
cSpProfile_dict_Sourav["FT0A_TPC"] = ROOT.TCanvas(
    "cSpProfile_FT0A_TPC_Sourav",
    "cSpProfile_FT0A_TPC_Sourav",
    800,
    600,
)
hSP_FT0A_TPC.Draw("colz")
hProfile_dict_Sourav["FT0A_TPC"].Draw("pe same")

n_bins = hSP_dict_Sourav["FT0C_TPC"].GetNbinsX()
cent_axis_limits = [
    hSP_dict_Sourav["FT0C_TPC"].GetXaxis().GetBinLowEdge(1),
    hSP_dict_Sourav["FT0C_TPC"].GetXaxis().GetBinUpEdge(n_bins),
]
cent_axis_title = hSP_dict_Sourav["FT0C_TPC"].GetXaxis().GetTitle()

resolution_name1 = "FT0C_FT0A_TPC"
resolution_title1 = r"R_{2} (FT0C#; FT0A, TPC)"

print(f"Creating resolution {resolution_name1}")

hResolution_dict_Sourav[resolution_name1] = ROOT.TH1F(
    f"hResolution_{resolution_name1}_Sourav",
    f";{cent_axis_title};{resolution_title1}",
    n_bins,
    cent_axis_limits[0],
    cent_axis_limits[1],
)
getResolution(
    hResolution_dict_Sourav[resolution_name1],
    "FT0C_TPC",
    "FT0A_TPC",
    "FT0C_FT0A",
    hProfile_dict_Sourav,
)
utils.setHistStyle(hResolution_dict_Sourav[resolution_name1], ROOT.kRed)

Sourav_dir = output_file.mkdir("Sourav")
Sourav_dir.cd()
for key in hSP_dict_Sourav.keys():
    hSP_dict_Sourav[key].Write()
    hProfile_dict_Sourav[key].Write()
    cSpProfile_dict_Sourav[key].Write()
for res in hResolution_dict_Sourav.values():
    res.Write()
