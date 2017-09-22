import cPickle as pic
import matplotlib.pyplot as plt
import numpy as np

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
# print(ogDict)


x = [1,1,1,2,2,2]
y = [1,1,1,2,2,2]
heatmap, xedges, yedges = np.histogram2d(x, y)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
plt.clf()
plt.imshow(heatmap.T, extent=extent, origin='lower')
plt.show()


# plt.imshow(a, cmap='hot', interpolation='nearest')
# plt.show()
#NOTES:
# (barcode: aa residue mutated _ resulting in codon)
# 150 aa protein
#count bar codes => counting aa mutations at residue n => 
	#if mutation x found more, then this implies mutation x is more liveable
#codon bias