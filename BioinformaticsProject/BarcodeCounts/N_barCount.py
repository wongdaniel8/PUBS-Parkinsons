import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# Input file: R2 to view barcodes. If "l" argument used, filename is SampleR2.
if sys.argv[1] == "l":
    fname = 'R2'
if sys.argv[1] == "s":
    fname = 'Undetermined_S0_L001_R2_001.fastq'

# Initialize dict of 26 instances of 0
freq = {}
for i in range(0,27):
    freq[i] = 0

with open(fname, "rU") as f:
    count = -1
    for line in f:
        if count == -1 or count %4 !=0:
            count +=1
            continue
        line = line.rstrip("\n")
        if count %4 == 0:
            nCount = line.count('N')
            freq[nCount] += 1
            count += 1

if sys.argv[1] == "l":
    filename = 'N_barCount_Sample.txt'
if sys.argv[1] == "s":
    filename = 'N_barCount_full.txt'

if os.path.exists(filename):
    os.remove(filename)

with open(filename, 'a') as wf:
    for k in freq:
        wf.write('%r: %r\n' %(k, freq[k]))

# Histogram of list of N counts - only renders locally
if sys.argv[1] == "l":
    plt.bar(freq.keys(), freq.values())
    plt.title('Barcode Frequency of N-counts')
    plt.show()
