import numpy as np

input_fname = "/home/niek/HSA_data/results/tmp"

#get number of lines
nofLines = -1 #-1 for the header line
with open(input_fname, 'r') as f:
    for i,l in enumerate(f):
        if len(l) > 1:
            nofLines += 1

ex_names = []
correlations = np.zeros(nofLines)
top1 = np.zeros(nofLines)
top2 = np.zeros(nofLines)
top3 = np.zeros(nofLines)
top5 = np.zeros(nofLines)
top10 = np.zeros(nofLines)
top20 = np.zeros(nofLines)
topall = np.zeros(nofLines)

with open(input_fname, 'r') as f:
    for i,line in enumerate(f):
        if 'top' in line or len(line) < 2:
            continue
        arr = line.strip().split()
        ex_names.append(arr[0])
        correlations[i], top1[i], top2[i], top3[i], top5[i], top10[i], top20[i], topall[i] = arr[1:]

print("number of experiments:", len(ex_names))

print("correlation:", correlations.mean())
print("top1:", top1.mean())
print("top2:", top2.mean())
print("top3:", top3.mean())
print("top5:", top5.mean())
print("top10:", top10.mean())
print("top20:", top20.mean())
print("topall:", topall.mean())
