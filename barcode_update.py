import cPickle as pic
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
from collections import Counter
from collections import defaultdict
import pandas as pd
import numpy as np
import hclust as hc
import scipy.sparse as sp
import os
from sys import argv

barcode_dict = argv[1]

data = pic.load(open(barcode_dict, "rb"))
translate_data = pic.load(open("translate.pkl", "rb"))
aminoNum = pic.load(open("aminotonumber.pkl", "rb"))

# # # # # # # # # # # # # # # # # # #
# Get counts

# Keep count of each residue position/mutant type
res2aa = defaultdict(list)
freq = {}

for key in data.keys():
	value = data[key][0]
	residue, codonMut = value.split("_")
	complement = codonMut.replace("T", "U")
	aminoAcid = translate_data[complement]

	# Save dict of each residue position and the mutations there
	res2aa[int(residue)].append(aminoAcid)

	if aminoAcid == "STOP":
		aminoAcid = "*"

	new_val = "_".join([residue, aminoAcid])
	if new_val not in freq.keys():
		freq[new_val] = 1
	else:
		freq[new_val] += 1

res_num = []
aa_name = []
freq_list = []

for key in freq:	# for the list of barcodes
	r, aa = key.split("_")		# get residue position & amino acid

	# Revert back to old "STOP" identifier
	if aa == "*":
		aa = "STOP"

	res_num.append(int(r))
	aa_name.append(aa)
	freq_list.append(freq[key])

# # # # # # # # # # # # # # # # # # #
# Create heatmap centered on average count
df = pd.DataFrame({'Residue Position': res_num, 'Amino Acid': aa_name, 'Frequency': freq_list})
df = df.pivot("Amino Acid", "Residue Position", "Frequency")

#plt.figure(figsize=(20, 8))
#ax = sns.heatmap(df, center = np.mean(freq_list))

#plt.setp(ax.get_xticklabels(), rotation=90)
#plt.setp(ax.get_yticklabels(), rotation=0)
#plt.title("Frequency of Barcodes (Centered at Average)")
#plt.show()
#plt.savefig("heatmap.png")

# # # # # # # # # # # # # # # # # # # # # #
# Cluster heatmap
ndf = df.fillna(0)

plt.figure(figsize=(20, 8))
cg = sns.clustermap(ndf, xticklabels = 1, center = np.mean(freq_list))

# Rotate tick labels
plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
plt.setp(cg.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)

plt.title("Frequency of Barcodes \nClustered by Euclidean Dist\n(color center = average # barcodes)")
#plt.show()
cg.savefig("cluster_heatmap.png")

# # # # # # # # # # # # # # # # # # # # # #
# Sanity Check
all_aas = [i for i in aminoNum]
missing_pms = {}

if os.path.exists('SanityCheck.txt'):
    os.remove('SanityCheck.txt')

f = open('SanityCheck.txt', 'a')

for k in sorted(res2aa.iterkeys()):

	# Iterate through the list of aa mutations and count
	counts = {}
	aa_list = res2aa[k]

	for aa in aa_list:
		if aa in counts:
			counts[aa] +=  1
		else:
			counts[aa] =   1

	# Confirm that the residue has entries for every amino acid
	if not set(counts) == set(all_aas):
		f.write("Residue %r is missing the following amino acid mutations: %s\n" %(k, ', '.join(set(all_aas) - set(counts))))
	else:
		f.write("Residue %r has mutations to all possible amino acids.\n" %k)

	# Create a tuple to pair each residue to its list of counts
	missing_pms[k] = len(set(all_aas) - set(counts))

# # # # # # # # # # # # # # # # # # #
# Bar chart showing number of missing mutations
d = missing_pms
fig, ax = plt.subplots()
plt.bar(range(len(d)), d.values(), align='center')
plt.xticks(range(len(d)), d.keys(), rotation = 90)
plt.yticks(xrange(0,25,5))
plt.xlabel("Residue Number")
plt.title("Number of Missing AA-Substitutions by Residue Number")
plt.savefig("missing_count.png")


#NOTES:
# (barcode: aa residue mutated _ resulting in codon)
# 150 aa protein
#count bar codes => counting aa mutations at residue n =>
	#if mutation x found more, then this implies mutation x is more liveable
#codon bias
