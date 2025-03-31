import ROOT
import pandas as pd
import numpy as np

import sys

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

# Define particle types and their file paths
particle_types = {
    "Helium3": "../theoretical_models/predictions_helium/v2_pt_Helium3_5360GeV_PbPb_{low}{high}.txt",
    "Proton": "../theoretical_models/predictions_proton/v2_pt_proton_5360GeV_PbPb_{low}{high}.txt",
}

# Initialize a list to store TGraphErrors objects
graphs = []

# Loop over particle types and centrality classes
for particle, file_template in particle_types.items():
    for low, high, color in centrality_classes:
        input_file_name = file_template.format(low=low, high=high)

        # Read data
        df = pd.read_csv(input_file_name, index_col=False, sep=" ")

        pt = df["pt"].tolist()
        v2 = df["v2"].tolist()
        v2_err = df["v2_err"].tolist()

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
        graph.SetName(f"gPred{particle}_{low}_{high}")
        utils.setHistStyle(graph, color)

        # Add graph to the list
        graphs.append(graph)

# Save graphs to file
output_file = ROOT.TFile("../theoretical_models/Predictions.root", "recreate")
for graph in graphs:
    graph.Write()
output_file.Close()
