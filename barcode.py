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



