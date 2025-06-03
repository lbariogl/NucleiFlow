import ROOT
import pandas as pd
import yaml
import argparse
import os

import sys

from datetime import datetime

sys.path.append("utils")
import utils as utils

parser = argparse.ArgumentParser(description="Create a DataFrame from ROOT files.")
parser.add_argument(
    "--config",
    dest='config_file',
    help='path to the YAML configuration file',
    default="",
)
args = parser.parse_args()
if args.config_file == "":
    print("** No config file provided. Exiting. **")
    exit()

config_file = open(args.config_file, "r")
config = yaml.full_load(config_file)

input_file_name = config["input_file_name"]
df_output_name = config["df_file_name"]
nuclei_tree_name = config["nuclei_tree_name"]
ep_tree_name = config["ep_tree_name"]

useSP = config["useSP"]


print("** Using scalar product **" if useSP else "** Using event plane **")

mandatory_selections = config["mandatory_selections"]

# harmonic parameter
if "harmonic" in config .keys():
    harmonic = config["harmonic"]
else:
    harmonic = 2
    print("** No 'harmonic' key found in config file. Using default value: 2 **")

p_train = config["p_train"]



nuclei_df = utils.get_df_from_tree(input_file_name, nuclei_tree_name)

nucleiflow_df = utils.get_df_from_tree(input_file_name, ep_tree_name)

complete_df = pd.concat([nuclei_df, nucleiflow_df], axis=1, join="inner")

utils.redefineColumns(
    complete_df,
    mass=utils.mass_helion,
    parameters=p_train,
    useSP=useSP,
    harmonic=harmonic,
)

complete_df.query(mandatory_selections, inplace=True)

complete_df.to_csv(df_output_name, index=False)