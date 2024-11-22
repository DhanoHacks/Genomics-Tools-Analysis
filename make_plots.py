import matplotlib.pyplot as plt
import numpy as np

# Data for plotting (timing results)
tools = ['plyranges', 'pyranges', 'gffutils-pybedtools', 'sql-pybedtools']

# Human and Mouse data
human_load_times = [28.811, 64.6853404045105, 449.48835945129395, 49.05035901069641]
mouse_load_times = [20.365, 44.9543559551239, 319.6335127353668, 36.47839665412903]

# Aggregate query times (mean) in seconds
human_agg1 = [15.091, 5.12, 44, 1.4]
mouse_agg1 = [14.152, 4.08, 34.2, 1.4]

human_agg2 = [15.164, 5.77, 43.8, 1.59]
mouse_agg2 = [14.072, 4.81, 34.7, 1.51]

human_agg3 = [0.085, 0.953, 5.52, 0.531]
mouse_agg3 = [0.08, 0.915, 5.77, 0.542]

# Interval query times (mean) in seconds
human_int1 = [0.174, 1.76, 60, 6.63]
mouse_int1 = [0.124, 3.23, 65, 6.46]

human_int2 = [0.024, 0.399, 89, 7.41]
mouse_int2 = [0.011, 0.519, 64, 7.28]

human_int3 = [0.405, 4.38, 88, 7.28]
mouse_int3 = [0.299, 3.43, 64, 6.83]

# Plot 1: Load times
plt.figure(figsize=(6, 4))
plt.plot(tools, human_load_times, marker='o', label='Human', color='blue')
plt.plot(tools, mouse_load_times, marker='o', label='Mouse', color='red')
plt.title('Time Taken to Load Data')
plt.ylabel('Time (s)')
plt.xlabel('Tools')
plt.legend()
plt.savefig("plot1.png",dpi = 300, bbox_inches="tight")

# Plot 2: Aggregate queries
plt.figure(figsize=(6, 4))
plt.plot(tools, human_agg1, marker='o', label='Human - Query 1', color='blue')
plt.plot(tools, mouse_agg1, marker='o', label='Mouse - Query 1', color='red')
plt.plot(tools, human_agg2, marker='o', label='Human - Query 2', color='green')
plt.plot(tools, mouse_agg2, marker='o', label='Mouse - Query 2', color='orange')
plt.plot(tools, human_agg3, marker='o', label='Human - Query 3', color='purple')
plt.plot(tools, mouse_agg3, marker='o', label='Mouse - Query 3', color='brown')
plt.title('Time Taken for Aggregate Queries')
plt.ylabel('Time (s)')
plt.xlabel('Tools')
plt.legend()
plt.savefig("plot2.png",dpi = 300, bbox_inches="tight")

# Plot 3: Interval queries
plt.figure(figsize=(6, 4))
plt.plot(tools, human_int1, marker='o', label='Human - Query 1', color='blue')
plt.plot(tools, mouse_int1, marker='o', label='Mouse - Query 1', color='red')
plt.plot(tools, human_int2, marker='o', label='Human - Query 2', color='green')
plt.plot(tools, mouse_int2, marker='o', label='Mouse - Query 2', color='orange')
plt.plot(tools, human_int3, marker='o', label='Human - Query 3', color='purple')
plt.plot(tools, mouse_int3, marker='o', label='Mouse - Query 3', color='brown')
plt.title('Time Taken for Interval Queries')
plt.ylabel('Time (s)')
plt.xlabel('Tools')
plt.legend()

# Show all plots
plt.savefig("plot3.png",dpi = 300, bbox_inches="tight")
