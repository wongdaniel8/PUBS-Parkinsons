import cPickle as pic
import matplotlib.pyplot as plt
import numpy as np

import cPickle as pic 
data = pic.load(open("allele_dic.pkl", "rb"))
translate_data = pic.load(open("translate.pkl", "rb"))
aminoNum = pic.load(open("aminotonumber.pkl", "rb"))
# print(data) 
# print(translate_data)
ogDict = {}
for key in data.keys():
	value = data[key][0]
	residue, codonMut = value.split("_")
	complement = codonMut.replace("T", "U")
	aminoAcid = translate_data[complement]
	num = aminoNum[aminoAcid]
	ogDict[key] = [residue, codonMut, aminoAcid, num]


x = []
y = []
for i in ogDict:
	a = ogDict[i]
	x += [int(a[0])]
	y += [a[3]]


heatmap, xedges, yedges = np.histogram2d(x, y, bins=50)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
# plt.colorbar(heatmap) 
plt.title("Frequency of Mutations at each Residue Number")
plt.xlabel("residue number")
plt.ylabel("amino acid")
plt.clf()
plt.imshow(heatmap.T, extent=extent, origin='lower')
plt.show()





#NOTES:
# (barcode: aa residue mutated _ resulting in codon)
# 150 aa protein
#count bar codes => counting aa mutations at residue n => 
	#if mutation x found more, then this implies mutation x is more liveable
#codon bias








