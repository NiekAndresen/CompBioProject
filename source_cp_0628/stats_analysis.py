import numpy as np

input_fname = "/home/niek/Computational Biology/CompBioProject/results/everything_statistics"

#get number of lines
nofLines = 0
with open(input_fname, 'r') as f:
    for i,l in enumerate(f):
        if 'top' in l or len(l) < 2:
            continue
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

next = 0
with open(input_fname, 'r') as f:
    for line in f:
        if 'top' in line or len(line) < 2:
            continue
        arr = line.strip().split()
        ex_names.append(arr[0])
        correlations[next], top1[next], top2[next], top3[next], top5[next], top10[next], top20[next], topall[next] = arr[1:]
        next += 1

print("number of experiments:", len(ex_names))

x = [1,2,3,5,10,20]
y = [top1.mean(), top2.mean(), top3.mean(), top5.mean(), top10.mean(), top20.mean()]
import matplotlib.pyplot as plt
plt.figure()
plt.plot(x,y,'x')
plt.show()
print("correlation:", correlations.mean())
print("top1:", top1.mean())
print("top2:", top2.mean())
print("top3:", top3.mean())
print("top5:", top5.mean())
print("top10:", top10.mean())
print("top20:", top20.mean())
print("topall:", topall.mean())
