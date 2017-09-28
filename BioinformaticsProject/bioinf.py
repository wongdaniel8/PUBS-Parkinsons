from Bio import SeqIO
from Bio.Seq import Seq
import sys

R1List = []
R2List = []
R3List = []

location = sys.argv[1]
wt = "ATGGATGTATTCATGAAAGGACTTTCAAAGGCCAAGGAGGGAGTTGTGGCTGCTGCTGAGAAAACCAAACAGGGTGTGGCAGAAGCAGCAGGAAAGACAAAAGAGGGTGTTCTCTATGTAGGCTCCAAAACCAAGGAGGGAGTGGTGCATGGTGTGGCAACAGTGGCTGAGAAGACCAAAG\
AGCAAGTGACAAATGTTGGAGGAGCAGTGGTGACGGGTGTGACAGCAGTAGCCCAGAAGACAGTGGAGGGAGCAGGGAGCATTGCAGCAGCCACTGGCTTTGTCAAAAAGGACCAGTTGGGCAAGAATGAAGAAGGAGCCCCACAGGAAGGAATTCTGGAAGATATGCCTGTGGATCCTGA\
CAATGAGGCTTATGAAATGCCTTCTGAGGAAGGGTATCAAGACTACGAACCTGAAGCC"
if location == "l": 
	files = ["R1", "R2", "R3"] #local
else:
	files = ["Undetermined_S0_L001_R1_001.fastq", "Undetermined_S0_L001_R2_001.fastq", "Undetermined_S0_L001_R3_001.fastq"] #global file names on server for R1, R2, R3

for file in files:
	if location == "l":
		if file == "R1":
			l = R1List
		if file == "R2":
			l = R2List
		if file == "R3":
			l = R3List
	else:
		if file == "Undetermined_S0_L001_R1_001.fastq":
			l = R1List
		if file == "Undetermined_S0_L001_R2_001.fastq":
			l = R2List
		if file == "Undetermined_S0_L001_R3_001.fastq":
			l = R3List

	with open(file, "rU") as f:
	    count = -1
	    for line in f:
	    	if count == -1 or count % 4 != 0:
	    		count += 1
	    		continue
	    	line = line.rstrip("\n") 
	    	if count % 4 == 0:
	    		l.append(line)
	   		count += 1

def merge(r1, r3):
	"""
	Returns a string of the same length as the WT sequence, resolving any conflicts between forward read r1 and backward read r3.
	Returned string should be free of any Ns.
	"""
	seq = Seq(r3)
	r3rc = seq.reverse_complement()
	returnSeq = "" #sequence to build and return

	#append r1 to returnSeq
	for i in range(0, len(r1)):
		if r1[i] == "N":
			returnSeq += wt[i]
		else:
			returnSeq += r1[i]

	r1overlap = r1[220:300] # a substring of r1 representing the overlapping region of r1 with the other read r3rc
	
	#build overlap region, if both are identical and non-N we take that, 
	#if the other read has a valid nucleotide, we take the valid one
	#if both are N, we input from the wt sequence
	for i in range(0, len(r1overlap)):
		if r1overlap[i] == r3rc[i] != "N":
			returnSeq += r1[i]
		elif r1overlap[i] == "N" and r3rc[i] != "N":
			returnSeq += r3rc[i]
		elif r1overlap[i] != "N" and r3rc[i] == "N":
			returnSeq += r1[i]
		else:
			returnSeq += wt[220 + i]

	#append r3rc to rteturnSeq
	for i in range(0, len(r3rc)):
		if r3rc[i] == "N":
			returnSeq += wt[220 + i]
		else:
			returnSeq += r3rc[i]

	print(returnSeq)
	return returnSeq


def findMutation(seq):
	"""
	Returns the index of the mutation found in seq based on heuristics
	"""
	return 


def findAminoAcid(seq, mutationIndex):
	"""
	returns the amino acid that corresponds to the mutation at index mutationIndex found in the sequence seq
	"""
	return 



#COMMANDS TO RUN THE SCRIPT===================================

merge (R1List[1], R3List[1])
# print(R1List[0])
# print(R3List[0])
# print(R3List)
# print(len(R1List))
# print(len(R2List))

# print(len(wt))
# print(len(R1List[0]))
# print(len(R2List[0]))
# print(len(R3List[0]))

# for i in range(0, len(R1List)):
	# print(R1List[i])
	# print(R3List[i])

#NOTES==========================================================
##indices relative to wt
# R1 = [0,300]
# R3 = [220, 420]

