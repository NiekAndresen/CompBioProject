import matplotlib.pyplot as plt

output_fname = "/home/niek/Computational Biology/CompBioProject/results/matches_found_top1"
fig_fname = "/home/niek/HSA_data/results/true_matches_by_number_of_experiments_top1.png"

x = []
y = []
with open(output_fname, 'r') as result:
    for line in result:
        if line.startswith('number'):
            continue
        numbers = line.strip().split()
        x += [numbers[0]]
        y += [numbers[1]]

plt.figure()
plt.plot(x,y, 'x')
plt.xlabel("number of experiments used\n")
plt.ylabel("number of distinct correct matches found\n")
plt.savefig(fig_fname)
