from Bio import SeqIO
from Bio.Seq import Seq
import sys
import numpy as np
import pickle


R1List = []
R2List = []
R3List = []
r1Score = []
r2Score = []
r3Score = []
wt = "ATGGATGTATTCATGAAAGGACTTTCAAAGGCCAAGGAGGGAGTTGTGGCTGCTGCTGAGAAAACCAAACAGGGTGTGGCAGAAGCAGCAGGAAAGACAAAAGAGGGTGTTCTCTATGTAGGCTCCAAAACCAAGGAGGGAGTGGTGCATGGTGTGGCAACAGTGGCTGAGAAGACCAAAG\
AGCAAGTGACAAATGTTGGAGGAGCAGTGGTGACGGGTGTGACAGCAGTAGCCCAGAAGACAGTGGAGGGAGCAGGGAGCATTGCAGCAGCCACTGGCTTTGTCAAAAAGGACCAGTTGGGCAAGAATGAAGAAGGAGCCCCACAGGAAGGAATTCTGGAAGATATGCCTGTGGATCCTGA\
CAATGAGGCTTATGAAATGCCTTCTGAGGAAGGGTATCAAGACTACGAACCTGAAGCC"

##load from pickle files
R1List = pickle.load( open( "R1List.p", "rb" ) )
R2List = pickle.load( open( "R2List.p", "rb" ) )
R3List = pickle.load( open( "R3List.p", "rb" ) )
# r1Score = pickle.load( open( "r1Score.p", "rb" ) )
# r2Score = pickle.load( open( "r2Score.p", "rb" ) )
# r3Score = pickle.load( open( "r3Score.p", "rb" ) )

location = sys.argv[1]
if location == "l": 
	files = ["R1", "R2", "R3"] #local
else:
	files = ["Undetermined_S0_L001_R1_001.fastq", "Undetermined_S0_L001_R2_001.fastq", "Undetermined_S0_L001_R3_001.fastq"] #global file names on server for R1, R2, R3

# for file in files:
# 	if location == "l":
# 		if file == "R1":
# 			l = R1List
# 		if file == "R2":
# 			l = R2List
# 		if file == "R3":
# 			l = R3List
# 	else:
# 		if file == "Undetermined_S0_L001_R1_001.fastq":
# 			l = R1List
# 		if file == "Undetermined_S0_L001_R2_001.fastq":
# 			l = R2List
# 		if file == "Undetermined_S0_L001_R3_001.fastq":
# 			l = R3List
# 	with open(file, "rU") as f:
# 	    count = -1
# 	    for line in f:
# 	    	if count == -1 or count % 4 != 0:
# 	    		count += 1
# 	    		continue
# 	    	line = line.rstrip("\n") 
# 	    	if count % 4 == 0:
# 	    		l.append(line)
# 	   		count += 1

for file in files:
	if location == "l":
		if file == "R1":
			l1 = r1Score
		if file == "R2":
			l1 = r2Score
		if file =="R3":
			l1 = r3Score
	else:
		if file == "Undetermined_S0_L001_R1_001.fastq":
			l1 = r1Score
		if file == "Undetermined_S0_L001_R2_001.fastq":
			l1 = r2Score
		if file == "Undetermined_S0_L001_R3_001.fastq":
			l1 = r3Score


	for record in SeqIO.parse(file, "fastq"):
		scores = record.letter_annotations["phred_quality"]
		median = 0#np.median(scores)
		l1.append(scores)
		# if location == "l":
		# 	if file == "R1":
		# 		r1Score.append((scores, median))
		# 	if file == "R2":
		# 		r2Score.append((scores, median))
		# 	if file == "R3":
		# 		r3Score.append((scores, median))
		# else:
		# 	if file == "Undetermined_S0_L001_R1_001.fastq":
		# 		r1Score.append((scores, median))
		# 	if file == "Undetermined_S0_L001_R2_001.fastq":
		# 		r2Score.append((scores, median))
		# 	if file == "Undetermined_S0_L001_R3_001.fastq":
		# 		r3Score.append((scores, median))




##dump into pickle files
# pickle.dump( R1List, open( "R1List.p", "wb" ) )
# pickle.dump( R2List, open( "R2List.p", "wb" ) )
# pickle.dump( R3List, open( "R3List.p", "wb" ) )
pickle.dump( r1Score, open( "r1Score.p", "wb" ) )
pickle.dump( r2Score, open( "r2Score.p", "wb" ) )
pickle.dump( r3Score, open( "r3Score.p", "wb" ) )




# print(len(R1List))
print(len(r1Score))
# print(len(R1List[0]))
# print(len(r3Score[0]))



# print(R1List)
# __N____NN____NNN_________NNNNNNNNNNNN
#                 NNNNNNNN__________NNNN______N____N______N________

#Ns at end mostly
#

def merge(r1, r3):
	"""
	Returns a string of the same length as the WT sequence, resolving any conflicts between forward read r1 and backward read r3.
	Returned string should be free of any Ns.
	"""
	##indices relative to wt
		# R1 = [0,299]
		# R3 = [219, 419] #inclusive bounds
	seq = Seq(r3)
	r3rc = seq.reverse_complement()
	returnSeq = "" #sequence to build and return

	r1mismatches = [i for i in range(len(r1)) if r1[i] != wt[i]]
	r3rcmismatches = [i for i in range(len(r3rc)) if r3rc[i] != wt[219 + i]]

	print(r1mismatches)
	print("r3rc last 50", r3rc[-50:])
	print("wt last 50: ", wt[-50:])
	print(r3rcmismatches)

	#append r1 to returnSeq
	for i in range(0, len(r1)):
		if r1[i] == "N":
			returnSeq += wt[i]
		else:
			returnSeq += r1[i]

	r1overlap = r1[219:300] #299 inclusive # a substring of r1 representing the overlapping region of r1 with the other read r3rc
	
	#build overlap region, if both are identical and non-N we take that, 
	#if the other read has a valid nucleotide, we take the valid one
	#if both are N, we input from the wt sequence
	for i in range(0, len(r1overlap)):
		if (r1overlap[i] == r3rc[i]) and r1overlap[i] != "N":
			returnSeq += r1overlap[i]
		elif r1overlap[i] == "N" and r3rc[i] != "N":
			returnSeq += r3rc[i]
		elif r1overlap[i] != "N" and r3rc[i] == "N":
			returnSeq += r1overlap[i]
		else:
			returnSeq += wt[219 + i] 

	#append r3rc to rteturnSeq
	for i in range(0, len(r3rc)):
		if r3rc[i] == "N":
			returnSeq += wt[219 + i]
		else:
			returnSeq += r3rc[i]

	for i in range(0, len(returnSeq)):
		if returnSeq[i] == "N":
			print(i)

	mismatches = [i for i in range(len(wt)) if returnSeq[i] != wt[i]]
	print("msi", mismatches)

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

# print(len(wt))
# print(len(R1List[0]))
# print(len(R2List[0]))
# print(len(R3List[0]))

# for i in range(0, len(R1List)):
	# print(R1List[i])
	# print(R3List[i])

#NOTES==========================================================
##indices relative to wt
	# R1 = [0,299]
	# R3 = [219, 419]


# count Ns in R1 N1
# count Ns in R2 N2







