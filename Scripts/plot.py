
import os
import pandas as pd
import numpy as np
import os
from matplotlib import pyplot as plt
from matplotlib import lines
from matplotlib import colors as mcolors

# %% [markdown]
# # 1.1 Explicit strong scaling

directory = "/home/tom/Documents/Uliege/Master2/HPC/Project2/Report/stats/strongscaling/explicit"
files = [filename for filename in os.listdir(directory) if filename.startswith("statistics_strongscaling_explicit_10000_it")]

all_dataframes = []
for f in files:
    part_df = pd.read_csv(os.path.join(directory,f), header=0)
    if not part_df.empty:
        all_dataframes.append(part_df[0:1])

combined_df = pd.concat(all_dataframes)
combined_df.sort_values(by=['Number of processes', 'Number of threads'])
groupby_df = combined_df.groupby(['Number of processes','Number of threads']).median().reset_index()
groupby_df['Number of processes']

# %% [markdown]
# ## 1.1.1 Plot average time as a function of process number number
# 

fig, (ax1, ax2) = plt.subplots(1,2, figsize=(12,5))

meanTimePerNbThreads = []
theoreticalMeanTime = []
uniqueProcesses = groupby_df['Number of processes'].unique()
for index,i in enumerate(uniqueProcesses):
    condition = (groupby_df['Number of processes'] == i) & (groupby_df['Number of threads'] == 1)
    meanTimePerNbThreads.append(groupby_df[condition]['Time per process'].values)
    
meanTime = [x for y, x in sorted(zip(uniqueProcesses, meanTimePerNbThreads))]
nbProcess = sorted(uniqueProcesses)
ax1.plot(nbProcess, meanTimePerNbThreads, marker='x',markersize=10)

for i in nbProcess:
    theoreticalMeanTime.append(meanTime[0]/i)

ax1.plot(nbProcess, theoreticalMeanTime, marker='x',markersize=10)
ax1.set(ylabel="Runtime (s)",xlabel="Number of processes")
ax1.legend(["In practice","In theory"], loc='best')

# Plot speedup factor as a function of process number
speedUp = []
theoreticalSpeedUp = []
for i in range(0,len(meanTime)):
    speedUp.append(meanTime[0]/meanTime[i])
    theoreticalSpeedUp.append(meanTime[0]/meanTime[i])

ax2.plot(nbProcess, speedUp, marker='x',markersize=10)
ax2.plot(nbProcess, nbProcess, marker='x',markersize=10)
ax2.set(ylabel="Speedup factor", xlabel="Number of processes")
ax2.legend(["Practical speedup", "Theoretical speedup"], loc  ='best')

fig.savefig("Explicit_runtime_Speedup_processes.svg")

# %% [markdown]
# ## 1.1.2 Plot average time as a function of thread number
# 

fig, (ax1, ax2) = plt.subplots(1,2, figsize=(12,5))

meanTimePerNbThreads = []
theoreticalMeanTime = []
uniqueThreads = groupby_df['Number of threads'].unique()
for index,i in enumerate(uniqueThreads):
    condition = (groupby_df['Number of threads'] == i) & (groupby_df['Number of processes'] == 1)
    meanTimePerNbThreads.append(groupby_df[condition]['Time per process'].values)

    
meanTime = [x for y, x in sorted(zip(groupby_df['Number of threads'].unique(), meanTimePerNbThreads))]
nbThreads = sorted(groupby_df['Number of threads'].unique())
ax1.plot(nbThreads, meanTime, marker='x',markersize=10)

for i in nbThreads:
    theoreticalMeanTime.append(meanTime[0]/i)
ax1.plot(nbThreads, theoreticalMeanTime, marker='x',markersize=10)
ax1.set(ylabel="Runtime (s)",xlabel="Number of threads")
ax1.legend(["In practice","In theory"],loc='best')

# Plot speedup factor as a function of thread number
speedUp = []
for i in range(0,len(meanTime)):
    speedUp.append(meanTime[0]/meanTime[i])
ax2.plot(nbThreads, speedUp, marker='x',markersize=10)
ax2.plot(nbThreads, nbThreads, marker='x',markersize=10)
ax2.set(ylabel="Speedup factor", xlabel="Number of threads")
ax2.legend(["Practical speedup", "Theoretical speedup"], loc  ='best')

fig.savefig("Explicit_runtime_Speedup_threads.svg")

# %% [markdown]
# # 1.2 Explicit weak scaling
# 

directory = "/home/tom/Documents/Uliege/Master2/HPC/Project2/Report/stats/weakscaling/explicit/process"

files = [filename for filename in os.listdir(directory) if filename.startswith("statistics_explicit_weakscaling")]

all_dataframes = []
for f in files:
    part_df = pd.read_csv(os.path.join(directory,f), header=0)
    if not part_df.empty:
        all_dataframes.append(part_df[0:1])

combined_df = pd.concat(all_dataframes)
combined_df = combined_df.dropna()
combined_df.sort_values(by=['Number of processes', 'Number of threads'], inplace=True)

# %%
groupby_df = combined_df.groupby(['Number of processes','Number of threads']).median().reset_index()
groupby_df.sort_values(by=['Number of processes', 'Number of threads'])

# %% [markdown]
# ## 1.2.1 Plot average time as a function of process number
# 

# %%
# Plot average time as a function of process number
fig, (ax1, ax2) = plt.subplots(1,2, figsize=(12,5))

meanTimePerNbThreads = []
theoreticalMeanTime = []
uniqueProcesses = groupby_df['Number of processes'].unique()
for index,i in enumerate(uniqueProcesses):
    condition = (groupby_df['Number of processes'] == i) & (groupby_df['Number of threads'] == 1)
    meanTimePerNbThreads.append(groupby_df[condition]['Time per process'].values)
    
meanTime = [x for y, x in sorted(zip(uniqueProcesses, meanTimePerNbThreads))]
nbProcess = sorted(uniqueProcesses)
ax1.plot(nbProcess, meanTime, marker='x',markersize=10)

theoreticalMeanTime = [meanTime[0]] * len(nbProcess)

ax1.plot(nbProcess, theoreticalMeanTime, marker='x',markersize=10)
ax1.set(ylabel="Average runtime (s)",xlabel="Number of processes")
ax1.legend(["In practice","In theory"],loc='best')

# Plot speedup factor as a function of process number
speedUp = []
theoreticalSpeedUp = []
for i in range(0,len(meanTime)):
    speedUp.append(meanTime[0]/meanTime[i])
    theoreticalSpeedUp.append(meanTime[0]/meanTime[i])
ax2.plot(nbProcess, speedUp, marker='x',markersize=10)
ax2.plot(nbProcess, [1]*len(nbProcess), marker='x',markersize=10)
ax2.set(ylabel="Average speedup factor", xlabel="Number of processes")
ax2.set_ylim([0.8,1.2])
ax2.legend(["Practical speedup", "Theoretical speedup"], loc  ='best')

fig.savefig("Weakscaling_processes_runtime_speedup.svg")




directory = "/home/tom/Documents/Uliege/Master2/HPC/Project2/Report/stats/weakscaling/explicit/threads"

files = [filename for filename in os.listdir(directory) if filename.startswith("statistics_explicit_weakscaling")]

all_dataframes = []
for f in files:
    part_df = pd.read_csv(os.path.join(directory,f), header=0)
    if not part_df.empty:
        all_dataframes.append(part_df[0:1])

combined_df = pd.concat(all_dataframes)
combined_df = combined_df.dropna()
combined_df.sort_values(by=['Number of processes', 'Number of threads'], inplace=True)

groupby_df = combined_df.groupby(['Number of processes','Number of threads']).median().reset_index()

# %% [markdown]
# ## 1.2.2 Plot average time/speedup as a function of thread number
# 

# Plot average time as a function of thread number
fig, (ax1, ax2) = plt.subplots(1,2, figsize=(12,5))

meanTimePerNbThreads = []
theoreticalMeanTime = []
uniqueThreads = groupby_df['Number of threads'].unique()
for index,i in enumerate(uniqueThreads):
    condition = (groupby_df['Number of threads'] == i) & (groupby_df['Number of processes'] == 1)
    meanTimePerNbThreads.append(groupby_df[condition]['Time per process'].values)

meanTime = [x for y, x in sorted(zip(groupby_df['Number of threads'].unique(), meanTimePerNbThreads))]
nbThreads = sorted(groupby_df['Number of threads'].unique())

ax1.plot(nbThreads, meanTime, marker='x',markersize=10)

theoreticalMeanTime = [meanTime[0]] * len(nbThreads)

ax1.plot(nbThreads, theoreticalMeanTime, marker='x',markersize=10)
ax1.set(ylabel="Average runtime (s)",xlabel="Number of threads")
ax1.legend(["In practice","In theory"],loc='best')

# Plot speedup factor as a function of process number
speedUp = []
theoreticalSpeedUp = []
for i in range(0,len(meanTime)):
    speedUp.append(meanTime[0]/meanTime[i])
    theoreticalSpeedUp.append(meanTime[0]/meanTime[i])
ax2.plot(nbThreads, speedUp, marker='x',markersize=10)
ax2.plot(nbThreads, [1]*len(nbThreads), marker='x',markersize=10)
ax2.set(ylabel="Average speedup factor", xlabel="Number of threads")
ax2.legend(["Practical speedup", "Theoretical speedup"], loc  ='best')

fig.savefig("Weakscaling_threads_runtime_speedup.svg")

# %% [markdown]
# # 2.1 Implicit strong scaling

directory = "/home/tom/Documents/Uliege/Master2/HPC/Project2/Report/stats/strongscaling/implicit/process"
files = [filename for filename in os.listdir(directory) if filename.startswith("statistics_strongscaling_implicit_correct")]

all_dataframes = []
for f in files:
    part_df = pd.read_csv(os.path.join(directory,f), header=0)

    if not part_df.empty:
        all_dataframes.append(part_df[0:1])

combined_df = pd.concat(all_dataframes)
combined_df.sort_values(by=['Number of processes', 'Number of threads'])
groupby_df = combined_df.groupby(['Number of processes','Number of threads']).median().reset_index()

# %% [markdown]
# ## 2.1.1 Plot average time as a function of process number number
# 

fig, (ax1, ax2) = plt.subplots(1,2, figsize=(12,5))

meanTimePerNbThreads = []
theoreticalMeanTime = []
uniqueProcesses = groupby_df['Number of processes'].unique()
for index,i in enumerate(uniqueProcesses):
    condition = (groupby_df['Number of processes'] == i) & (groupby_df['Number of threads'] == 1)
    meanTimePerNbThreads.append(groupby_df[condition]['Time per process'].values)
    
meanTime = [x for y, x in sorted(zip(uniqueProcesses, meanTimePerNbThreads))]
nbProcess = sorted(uniqueProcesses)
ax1.plot(nbProcess, meanTimePerNbThreads, marker='x',markersize=10)

for i in nbProcess:
    theoreticalMeanTime.append(meanTime[0]/i)

ax1.plot(nbProcess, theoreticalMeanTime, marker='x',markersize=10)
ax1.set(ylabel="Runtime (s)",xlabel="Number of processes")
ax1.legend(["In practice","In theory"],loc='best')

# Plot speedup factor as a function of process number
speedUp = []
theoreticalSpeedUp = []
for i in range(0,len(meanTime)):
    speedUp.append(meanTime[0]/meanTime[i])
    theoreticalSpeedUp.append(meanTime[0]/meanTime[i])
ax2.plot(nbProcess, speedUp, marker='x',markersize=10)
ax2.plot(nbProcess, nbProcess, marker='x',markersize=10)
ax2.set(ylabel="Speedup factor", xlabel="Number of processes")
ax2.legend(["Practical speedup", "Theoretical speedup"], loc  ='best')

fig.savefig("Implicit_strongscaling_runtime_speedup_processes.svg")

# %% [markdown]
# ## 2.1.2 Plot average time as a function of thread number
# 

directory = "/home/tom/Documents/Uliege/Master2/HPC/Project2/Report/stats/strongscaling/implicit/threads"
files = [filename for filename in os.listdir(directory) if filename.startswith("statistics_strongscaling_implicit_correct")]

all_dataframes = []
for f in files:
    part_df = pd.read_csv(os.path.join(directory,f), header=0)

    if not part_df.empty:
        all_dataframes.append(part_df[0:1])

combined_df = pd.concat(all_dataframes)
combined_df.sort_values(by=['Number of processes', 'Number of threads'])
groupby_df = combined_df.groupby(['Number of processes','Number of threads']).median().reset_index()

# %%
# Plot average time as a function of thread number
fig, (ax1, ax2) = plt.subplots(1,2, figsize=(12,5))

meanTimePerNbThreads = []
theoreticalMeanTime = []
uniqueThreads = groupby_df['Number of threads'].unique()
for index,i in enumerate(uniqueThreads):
    condition = (groupby_df['Number of threads'] == i) & (groupby_df['Number of processes'] == 1)
    meanTimePerNbThreads.append(groupby_df[condition]['Time per process'].values)

    
meanTime = [x for y, x in sorted(zip(groupby_df['Number of threads'].unique(), meanTimePerNbThreads))]
nbThreads = sorted(groupby_df['Number of threads'].unique())
ax1.plot(nbThreads, meanTime, marker='x',markersize=10)

for i in nbThreads:
    theoreticalMeanTime.append(meanTime[0]/i)

ax1.plot(nbThreads, theoreticalMeanTime, marker='x',markersize=10)
ax1.set(ylabel="Runtime (s)",xlabel="Number of threads")
ax1.legend(["In practice","In theory"],loc='best')

# Plot speedup factor as a function of thread number
speedUp = []
for i in range(0,len(meanTime)):
    speedUp.append(meanTime[0]/meanTime[i])
ax2.plot(nbThreads, speedUp, marker='x',markersize=10)
ax2.plot(nbThreads, nbThreads, marker='x',markersize=10)
ax2.set(ylabel="Speedup factor", xlabel="Number of threads")
ax2.legend(["Practical speedup", "Theoretical speedup"], loc  ='best')

fig.savefig("Implicit_strongscaling_runtime_speedup_threads.svg")

# %% [markdown]
# # 2.2 Implicit weak scaling

directory = "/home/tom/Documents/Uliege/Master2/HPC/Project2/Report/stats/weakscaling/implicit"
files = [filename for filename in os.listdir(directory) if filename.startswith("statistics_implicit_weakscaling_correct") and filename.find('23') != -1]

all_dataframes = []
for f in files:
    part_df = pd.read_csv(os.path.join(directory,f), header=0)

    if not part_df.empty:
        all_dataframes.append(part_df[0:1])

combined_df = pd.concat(all_dataframes)
combined_df.sort_values(by=['Number of processes', 'Number of threads'])
groupby_df = combined_df.groupby(['Number of processes','Number of threads']).median().reset_index()

# %% [markdown]
# ## 2.2.1 Plot average time as a function of process number

fig, (ax1, ax2) = plt.subplots(1,2, figsize=(12,5))

meanTimePerNbThreads = []
theoreticalMeanTime = []
uniqueProcesses = groupby_df['Number of processes'].unique()
for index,i in enumerate(uniqueProcesses):
    condition = (groupby_df['Number of processes'] == i) & (groupby_df['Number of threads'] == 1)
    meanTimePerNbThreads.append(groupby_df[condition]['Time per process'].values)

    
meanTime = [x for y, x in sorted(zip(uniqueProcesses, meanTimePerNbThreads))]
nbProcess = sorted(uniqueProcesses)
theoreticalMeanTime = [meanTime[0]] * len(nbProcess)
worstTime = nbProcess * meanTime[0]
ax1.plot(nbProcess, meanTime, marker='x',markersize=10)
ax1.plot(nbProcess, theoreticalMeanTime, marker='x',markersize=10)
ax1.plot(nbProcess, worstTime, marker='x',markersize=10)

ax1.set(ylabel="Average runtime (s)",xlabel="Number of processes")
ax1.legend(["In practice","In theory", "Worst case"],loc='best')

# Plot speedup factor as a function of process number
speedUp = []
theoreticalSpeedUp = []
for i in range(0,len(meanTime)):
    speedUp.append(meanTime[0]/meanTime[i])
    theoreticalSpeedUp.append(meanTime[0]/meanTime[i])
worstSpeedup = list(map(lambda x : 1/x, nbProcess))
ax2.plot(nbProcess, speedUp, marker='x',markersize=10)
ax2.plot(nbProcess, [1]*len(nbProcess), marker='x',markersize=10)
ax2.plot(nbProcess, worstSpeedup, marker='x',markersize=10)
ax2.set(ylabel="Average speedup factor", xlabel="Number of processes")
# ax2.set_ylim([0.8,1.2])
ax2.legend(["Practical speedup", "Theoretical speedup", "No speedup"], loc  ='best')

fig.savefig("Implicit_weakscaling_processes_runtime_speedup.svg")

# %% [markdown]
# ## 2.2.2 Plot average time as a function of thread number
# 
fig, (ax1, ax2) = plt.subplots(1,2, figsize=(12,5))

meanTimePerNbThreads = []
theoreticalMeanTime = []
uniqueThreads = groupby_df['Number of threads'].unique()
for index,i in enumerate(uniqueThreads):
    condition = (groupby_df['Number of threads'] == i) & (groupby_df['Number of processes'] == 1)
    meanTimePerNbThreads.append(groupby_df[condition]['Time per process'].values)
    
meanTime = [x for y, x in sorted(zip(uniqueThreads, meanTimePerNbThreads))]
nbThreads = sorted(groupby_df['Number of threads'].unique())
worstTime = nbThreads * meanTime[0]
theoreticalMeanTime = [meanTime[0]] * len(nbThreads)

ax1.plot(nbThreads, meanTime, marker='x',markersize=10)
ax1.plot(nbThreads, theoreticalMeanTime, marker='x',markersize=10)
ax1.plot(nbThreads, worstTime, marker='x',markersize=10)
ax1.set(ylabel="Average runtime (s)",xlabel="Number of threads")
ax1.legend(["In practice","In theory", "Worst case"],loc='best')

# Plot speedup factor as a function of process number
speedUp = []
theoreticalSpeedUp = []
for i in range(0,len(meanTime)):
    speedUp.append(meanTime[0]/meanTime[i])
    theoreticalSpeedUp.append(meanTime[0]/meanTime[i])

worstSpeedup = list(map(lambda x : 1/x, nbThreads))
ax2.plot(nbThreads, speedUp, marker='x',markersize=10)
ax2.plot(nbThreads, [1]*len(nbThreads), marker='x',markersize=10)
ax2.plot(nbProcess, worstSpeedup, marker='x',markersize=10)

ax2.set(ylabel="Average speedup factor", xlabel="Number of threads")
ax2.legend(["Practical speedup", "Theoretical speedup", "No speedup"], loc  ='best')

fig.savefig("Implicit_weakscaling_threads_runtime_speedup.svg")
