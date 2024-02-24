import ROOT
from hipe4ml.tree_handler import TreeHandler
import pandas as pd

import sys
sys.path.append('utils')
import utils as utils


input_file_name = '~/debug/AnalysisResults_trees.root'
input_dir_name = 'nucleiFlow'
output_file_name = 'output.root'

cent_detector_label = 'FT0C'
cent_limits = [30, 50]

nuclei_hdl = TreeHandler(input_file_name, 'O2nucleitable', folder_name='DF')
nuclei_df = nuclei_hdl._full_data_frame

nucleiflow_hdl = TreeHandler(input_file_name, 'O2nucleitableflow', folder_name='DF')
nucleiflow_df = nucleiflow_hdl._full_data_frame

complete_df = pd.concat([nuclei_df, nucleiflow_df], axis=1, join='inner')






