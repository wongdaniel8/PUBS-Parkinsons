import cPickle as pic
barcodes = ['GATGCTTATGTACGTAGA','AGGAGCAACTCCAACGGT']

#unpack data structures
data = pic.load(open("Allele_dic.pkl", "rb"))
trans = pic.load(open("translate.pkl","rb"))

#for each barcodes in a large list of barcodes received from experiment
barcounts = {}
for m in barcodes:
	if m not in barcounts:
		barcounts[m] = 1
	else:
		barcounts[m] 

for i in barcodes:
	mut = data[i]
	mut = mut[0]
	nummut = mut.split('_')

	#convert to RNA
	RNAnummut = ''
	for j in nummut[1]:
		if j == 'T':
			j = 'U'
		RNAnummut += j

	#add to list: now AA number, DNA codon, AA converted
	nummut += trans[RNAnummut]
	print('The amino acid number ' + str(mut[0]) + 
		' was changed to the amino acid with letter code ' +
		str(nummut[2]) + '.')



=======
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
print(ogDict)



#NOTES:
# (barcode: aa residue mutated _ resulting in codon)
# 150 aa protein
#count bar codes => counting aa mutations at residue n => 
	#if mutation x found more, then this implies mutation x is more liveable
#codon bias
>>>>>>> 26675a6cf5c9c6635a04d3bdfad994d4551b893a
