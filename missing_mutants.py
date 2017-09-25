import cPickle as pic
from collections import defaultdict
from collections import OrderedDict

bar2cod = pic.load(open("allele_dic.pkl", "rb"))
cod2aa = pic.load(open("translate.pkl", "rb"))
aa2num = pic.load(open("aminotonumber.pkl", "rb"))

## Create dictionary to match residue positions to all mutations there
res2aa = defaultdict(list)

for k,v in bar2cod.iteritems():
    res, codon = v[0].split("_")
    aa = cod2aa[codon.replace("T","U")]
    res2aa[res].append(aa)

# Find counts for missing point mutations for each residue
all_aas = [i for i in aa2num]
missing_pms = {}

for k, aa_list in sorted(res2aa.iteritems()):

    # Iterate through the list of aa mutations and count
    counts = {}
    for aa in aa_list:
        if aa in counts:
            counts[aa] +=  1
        else:
            counts[aa] =   1

    # Confirm that the residue has entries for every amino acid
    if not set(counts) == set(all_aas):
        print "Residue %r is missing the following amino acid mutations: %s" %(k, ', '.join(set(all_aas) - set(counts)))
        missing_pms[k] = len(set(all_aas) - set(counts))
    else:
        print "Residue %r has mutations to all possible amino acids." %k

# Bar chart of residue position and number of missing point mutation barcodes
import matplotlib.pyplot as plt

d = OrderedDict(sorted(missing_pms.items()))

plt.bar(range(len(d)), d.values(), align='center')
plt.xticks(range(len(d)), d.keys())
plt.xlabel("Residue Number")
plt.title("Number of Missing AA-Substitutions by Residue Number")

plt.show()
