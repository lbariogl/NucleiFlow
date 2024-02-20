import utils as utils
import ROOT

import sys
sys.path.append('utils')


input_file_name = './run_task/AnalysisResults_small.root'
input_dir_name = 'nucleiFlow'
output_file_name = 'output.root'

cent_detector_label = 'FT0C'
cent_limits = [30, 50]

input_file = ROOT.TFile(input_file_name)
output_file = ROOT.TFile(output_file_name, 'RECREATE')

hCentFV0A = input_file.Get(f'{input_dir_name}/hCentFV0A')
hCentFV0A.SetDirectory(0)
hCentFT0A = input_file.Get(f'{input_dir_name}/hCentFT0A')
hCentFT0A.SetDirectory(0)
hCentFT0C = input_file.Get(f'{input_dir_name}/hCentFT0C')
hCentFT0C.SetDirectory(0)
hCentFT0M = input_file.Get(f'{input_dir_name}/hCentFT0M')
hCentFT0M.SetDirectory(0)

hSpFT0AvsNsigmaHe3VsPtvsCent = input_file.Get(
    f'{input_dir_name}/hSpFT0AvsNsigmaHe3VsPtvsCent')
hSpFT0CvsNsigmaHe3VsPtvsCent = input_file.Get(
    f'{input_dir_name}/hSpFT0CvsNsigmaHe3VsPtvsCent')
hSpFV0AvsNsigmaHe3VsPtvsCent = input_file.Get(
    f'{input_dir_name}/hSpFV0AvsNsigmaHe3VsPtvsCent')

pt_axis = hSpFT0AvsNsigmaHe3VsPtvsCent.GetAxis(2)
n_pt_bins = pt_axis.GetNbins()

v2_ft0c_val = []
v2_ft0c_err = []
v2_ft0a_val = []
v2_ft0a_err = []
v2_fv0a_val = []
v2_fv0a_err = []

for i_pt in range(1, n_pt_bins + 1):

    pt_dir = output_file.mkdir(f'pt_{i_pt}')
    pt_dir.cd()

    cent_low = cent_limits[0]
    cent_up = cent_limits[1]

    pt_low = pt_axis.GetBinLowEdge(i_pt)
    pt_up = pt_axis.GetBinUpEdge(i_pt)

    utils.getCompleteCanvas(hSpFT0CvsNsigmaHe3VsPtvsCent, cent_low, cent_up, i_pt, i_pt,
                            pt_dir, v2_ft0c_val, v2_ft0c_err, qvec_detector_label='FT0C', suffix=i_pt)

    utils.getCompleteCanvas(hSpFT0AvsNsigmaHe3VsPtvsCent, cent_low, cent_up, i_pt, i_pt,
                            pt_dir, v2_ft0a_val, v2_ft0a_err, qvec_detector_label='FT0A', suffix=i_pt)

    utils.getCompleteCanvas(hSpFV0AvsNsigmaHe3VsPtvsCent, cent_low, cent_up, i_pt, i_pt,
                            pt_dir, v2_fv0a_val, v2_fv0a_err, qvec_detector_label='FV0A', suffix=i_pt)
