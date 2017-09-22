import cPickle as pic
from collections import defaultdict


bar2cod = pic.load(open("allele_dic.pkl", "rb"))
cod2aa = pic.load(open("translate.pkl", "rb"))
aa2num = pic.load(open("aminotonumber.pkl", "rb"))

# take a look at data structures
#print bar2cod
#print cod2aa
#print aa2num

# bar2cod format:{'ACCAAAACCCGGAGAACA': ['50_CGG'], ...}
# cod2aa format:{'ACC': 'T', 'GUC': 'V', 'ACA': 'T'...}
# aa2num format: {'A': 9, 'C': 8...}

# For each barcode in the barcode-to-codon conversion dictionary, pair
# it to the proper amino acid, residue, and aa "number" ID
bar2aa = {}
bar2num = {}
res2aa = defaultdict(list)

for k,v in bar2cod.iteritems():

    # Separate residue number and codon
    res, codon = v[0].split("_")

    # Save dictionary with barcode to amino acid conversion, substituting T -> U
    aa = cod2aa[codon.replace("T","U")]
    bar2aa[k] = aa

    # Save dictionary with barcode to amino acid "number" conversion
    bar2num[k] = aa2num[aa]

    # Save dict of each residue position and the mutations there
    res2aa[res].append(aa)

# Print out key-value pairs
#for k, v in bar2aa.items():
#    print(k, v)
#for k, v in bar2num.items():
#    print(k, v)
#for k, v in res2aa.items():
#    print (k,v)

#################################################################
# Count the number of each amino acid substitution at each residue position

# Initialize list that will hold residue/amino acid count tuples
count_list = []

# Create list of the keys of the amino acid dictionary, our list of all possible aas
all_aas = [i for i in aa2num]
#print all_aas

# For each residue position
for k, aa_list in res2aa.iteritems():

    # Iterate through the list of aa mutations and count
    counts = {}
    for aa in aa_list:
        if aa in counts:
            counts[aa] +=  1
        else:
            counts[aa] =   1

    # Confirm that the residue has entries for every amino acid
    #print set(counts)
    #print set(aa_list)
    if not set(counts) == set(all_aas):
        print "Residue %r is missing the following amino acid mutations: %r" %(k, set(all_aas) - set(counts))
    else:
        print "Residue %r has mutations to all possible amino acids." %k

    # Create a tuple to pair each residue to its list of counts
    res_count = (k, counts)
    count_list.append(res_count)

# Print the dictionary of residue counts
print count_list
