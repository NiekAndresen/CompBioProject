import numpy as np
import matplotlib.pyplot as plt

freqs = np.array([])
dists = np.array([])
with open("/home/niek/HSA_data/results/foundFrequency", 'r') as f:
    for line in f:
        arr = line.split()
        freqs = np.append(freqs, [int(arr[0][:-1])])
        dists = np.append(dists, [float(arr[1])])

pos = dists<=20.
posfreqs = freqs[pos]
posdists = dists[pos]
negfreqs = freqs[~pos]
negdists = dists[~pos]
#print(len(freqs), len(posfreqs), len(negfreqs)) #13133 2498 10635 for ex_1_2

plt.figure()
ax = plt.gca()
ax.set_ylim(0,.2)
weights1 = np.ones_like(posfreqs)/float(len(posfreqs))
weights2 = np.ones_like(negfreqs)/float(len(negfreqs))
n, bins, patches1 = plt.hist(posfreqs, 50, normed=0, facecolor='green', alpha=0.5, weights=weights1, label='"true crosslinks" (<20A)')
n, bins, patches2 = plt.hist(negfreqs, 50, normed=0, facecolor='red', alpha=0.5, weights=weights2, label='not crosslinks')
plt.legend()
plt.show()

