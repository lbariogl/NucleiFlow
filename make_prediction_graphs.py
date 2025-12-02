import ROOT
import pandas as pd
import numpy as np

import sys

import glob
import os

sys.path.append("utils")
import utils as utils

# Define centrality classes and colors
centrality_classes = [
    (0, 10, ROOT.kCyan + 2),
    (10, 20, ROOT.kGreen),
    (20, 30, ROOT.kOrange),
    (30, 40, ROOT.kRed + 1),
    (40, 60, ROOT.kViolet + 1),
]

input_dir = "../theoretical_models/cluster_thermal_v2/"

# Define particle types and their file paths
models = {
    "Coal": "{input_dir}IP_v2_{low}-{high}_he3.csv",
    "Therm": "{input_dir}IP_v2_{low}-{high}_he3_thermal.csv",
    # "Proton": "../theoretical_models/predictions_proton/v2_pt_proton_5360GeV_PbPb_{low}{high}.txt",
}

# Initialize a list to store TGraphErrors objects
graphs = []

# Loop over particle types and centrality classes
for model, file_template in models.items():
    for low, high, color in centrality_classes:
        input_file_name = file_template.format(input_dir=input_dir, low=low, high=high)

        # Read data
        df = pd.read_csv(input_file_name, index_col=False, sep=",")
        if model == "Coal":
            df = df.iloc[:-2]

        pt = df["pt"].tolist()
        v2 = df["v2"].tolist()
        v2_err = df["std"].tolist()

        n_points = len(v2_err)

        # Convert to numpy arrays
        pt_array = np.array(pt, dtype=np.float64)
        v2_array = np.array(v2, dtype=np.float64)
        v2_err_array = np.array(v2_err, dtype=np.float64)
        zeros_array = np.zeros(n_points, dtype=np.float64)

        # Create TGraphErrors
        graph = ROOT.TGraphErrors(
            n_points, pt_array, v2_array, zeros_array, v2_err_array
        )
        graph.SetName(f"gPred{model}_{low}_{high}")
        utils.setHistStyle(graph, color)

        # Add graph to the list
        graphs.append(graph)

# Save graphs to file
output_file = ROOT.TFile(
    "../theoretical_models/Predictions_october2025.root", "recreate"
)
for graph in graphs:
    graph.SetTitle(";#it{p}_{T} (GeV/#it{c}); v_{2}")
    graph.Write()


# # --- NEW CODE: Create TGraphErrors for all pt*.csv files in ALICE_he3_5360 ---

# pt_csv_files = glob.glob(os.path.join(input_dir, "pt*.csv"))

# for csv_path in pt_csv_files:
#     # Get the base name without .csv
#     base = os.path.basename(csv_path)
#     name = base.replace(".csv", "")
#     hist_name = f"hPredPhi_{name}"

#     # Read the csv using pandas (header expected: 3 columns: x, y, yerr)
#     df = pd.read_csv(csv_path, index_col=False, sep=",")
#     x = np.array(df.iloc[:, 0], dtype=np.float64)
#     y = np.array(df.iloc[:, 1], dtype=np.float64)
#     yerr = np.array(df.iloc[:, 2], dtype=np.float64)
#     n = len(x)

#     # Assume x are bin centers, estimate bin edges
#     if n > 1:
#         bin_widths = np.diff(x)
#         avg_width = np.mean(bin_widths)
#         bin_edges = np.concatenate(
#             ([x[0] - avg_width / 2], x[:-1] + bin_widths / 2, [x[-1] + avg_width / 2])
#         )
#     else:
#         bin_edges = np.array([x[0] - 0.5, x[0] + 0.5])

#     hist = ROOT.TH1D(hist_name, ";#phi (rad); counts", n, bin_edges)
#     for i in range(n):
#         hist.SetBinContent(i + 1, y[i])
#         hist.SetBinError(i + 1, yerr[i])
#     output_file.cd()
#     hist.Write()
