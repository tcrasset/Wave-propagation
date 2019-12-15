import pandas as pd

filename = "/home/tom/Documents/Uliege/Master2/HPC/Project2/Results/project2_out_nbproc_2-nbthreads_1.txt"
df = pd.read_csv(filename, header=0)

print(df)