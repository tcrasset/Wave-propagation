import pandas as pd
import numpy as np
import os

directory = "/home/tom/Documents/Uliege/Master2/HPC/Project2/Results"
files = [filename for filename in os.listdir(directory) if filename.startswith("statistics")]


all_dataframes = []
for f in files:

    part_df = pd.read_csv(os.path.join(directory,f), header=0)

    if not part_df.empty:
        all_dataframes.append(part_df)

combined_df = pd.concat(all_dataframes)

