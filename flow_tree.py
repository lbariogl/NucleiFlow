import ROOT
from hipe4ml.tree_handler import TreeHandler
import pandas as pd

import sys
sys.path.append('utils')
import utils as utils
from flow import FlowMaker


input_file_name = '~/test_new_format/AnalysisResults_trees.root'
output_file_name = 'output.root'
outuput_file = ROOT.TFile(output_file_name, 'recreate')

cent_detector_label = 'FT0C'
cent_limits = [30, 50]

nuclei_hdl = TreeHandler(input_file_name, 'O2nucleitable', folder_name='DF')
nuclei_df = nuclei_hdl._full_data_frame

nucleiflow_hdl = TreeHandler(input_file_name, 'O2nucleitableflow', folder_name='DF')
nucleiflow_df = nucleiflow_hdl._full_data_frame

complete_df = pd.concat([nuclei_df, nucleiflow_df], axis=1, join='inner')
utils.redifineColumns(complete_df)

# select centrality
complete_df.query(f'fCentFT0C > {cent_limits[0]} and fCentFT0C < {cent_limits[1]}')

# other selections
selections = 'fSign < 0 and abs(fEta) < 0.8 and abs(fDCAxy) < 0.1 and fAvgItsClusSize > 4.5 and fTrackedAsHe == True'
complete_df.query(selections, inplace=True)

flow_maker = FlowMaker()
flow_maker.data_df = complete_df
flow_maker.pt_bins = [1, 1.5, 2., 2.5, 3., 3.5]
flow_maker.output_dir = outuput_file

flow_maker.make_flow()
flow_maker.dump_to_output_dir()

