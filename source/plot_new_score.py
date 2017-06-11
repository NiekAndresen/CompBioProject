import numpy as np
import matplotlib.pyplot as plt

output_fname = "/home/niek/Computational Biology/CompBioProject/results/new_score_precision_everything"
fig_fname = "/home/niek/HSA_data/results/new_score_precision_everything.png"

x = []
y_ns = []
y_av = []
with open(output_fname, 'r') as f:
    for line in f:
        if line.startswith('x'):
            continue
        l = line.strip().split()
        x += [int(float(l[0]))]
        y_ns += [l[1]]
        y_av += [l[2]]

plt.figure()
plt.plot(x[:1000],y_ns[:1000], label="new score")
plt.plot(x[:1000],y_av[:1000], label="average score")
plt.xlabel("number of top scores used")
plt.ylabel("precision")
plt.legend()
plt.savefig(fig_fname)
plt.show()
        
