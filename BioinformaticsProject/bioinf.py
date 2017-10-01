from Bio import SeqIO
from Bio.Seq import Seq
import sys
import numpy as np
import cPickle as pickle
import matplotlib.pyplot as plt

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
translate = pickle.load(open("translate.pkl", "rb"))
aminoNum = pickle.load(open("aminotonumber.pkl", "rb"))
R1List = pickle.load( open( "R1List.p", "rb" ) )
R2List = pickle.load( open( "R2List.p", "rb" ) )
R3List = pickle.load( open( "R3List.p", "rb" ) )
r1Score = pickle.load( open( "r1Score.p", "rb" ) )
r2Score = pickle.load( open( "r2Score.p", "rb" ) )
r3Score = pickle.load( open( "r3Score.p", "rb" ) )

location = sys.argv[1]
if location == "l": 
	files = ["R1", "R2", "R3"] #local
else:
	files = ["Undetermined_S0_L001_R1_001.fastq", "Undetermined_S0_L001_R2_001.fastq", "Undetermined_S0_L001_R3_001.fastq"] #global file names on server for R1, R2, R3

##parse out R(n)List containing sequence reads
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


##parse out Qscores and instantiate r(n)Score lists
# Qscore = dict((chr(i),i-33) for i in range(33,74)) 
# for file in files:
# 	if location == "l":
# 		if file == "R1":
# 			l1 = r1Score
# 		if file == "R2":
# 			l1 = r2Score
# 		if file =="R3":
# 			l1 = r3Score
# 	else:
# 		if file == "Undetermined_S0_L001_R1_001.fastq":
# 			l1 = r1Score
# 		if file == "Undetermined_S0_L001_R2_001.fastq":
# 			l1 = r2Score
# 		if file == "Undetermined_S0_L001_R3_001.fastq":
# 			l1 = r3Score
# 	with open(file, "rU") as f:
# 	    count = -3
# 	    for line in f:
# 	    	if count == -1 or count % 4 != 0:
# 	    		count += 1
# 	    		continue
# 	    	if count % 4 == 0:
# 		    	line = line.rstrip("\n") 
# 		    	scores = []
# 		    	for char in line:
# 		    		scores.append(Qscore[char])
# 	    		l1.append(scores)
# 	   		count += 1


##dump into pickle files
# pickle.dump( R1List, open( "R1List.p", "wb" ) )
# pickle.dump( R2List, open( "R2List.p", "wb" ) )
# pickle.dump( R3List, open( "R3List.p", "wb" ) )
# pickle.dump( r1Score, open( "r1Score.p", "wb" ) )
# pickle.dump( r2Score, open( "r2Score.p", "wb" ) )
# pickle.dump( r3Score, open( "r3Score.p", "wb" ) )


# print(r3Score)
# print(len(r1Score[0]))
# print(len(R1List[0]))
# print(r3Score)
# print(len(R1List[0]))
# print(len(r3Score[0]))



def scoreHistogram(read):
	if read == "r1":
		scoreMatrix = np.asarray(r1Score)
	if read == "r3":
		scoreMatrix = np.asarray(r3Score)
	means = np.mean(scoreMatrix, axis=0)
	print(means)
	plt.bar(range(len(means)), means)
	plt.show()
# scoreHistogram("r1")
# scoreHistogram("r3")



def merge(r1, r3, index):
	"""
	Retursn an index of the predicted mutation along with the amino acid correlate, returns -1 if not possible
	"""
	##indices relative to wt
		# R1 = [0,299]
		# R3 = [219, 419] #inclusive bounds
	seq = Seq(r3)
	r3rc = seq.reverse_complement()
	returnSeq = "" #sequence to build and return

	
	# print(r1mismatches)
	# print("r3rc last 50", r3rc[-50:])
	# print("wt last 50: ", wt[-50:])
	# print(r3rcmismatches)

	# #easy case of no N's in first half of r1 or last half of wt? (check this should be half or some proportion of read1 and r3rc, not necessarily 420 / 2 = 210)
	# if "N" not in r1[0:210] and "N" not in r3rc[-210:]:
	# 	return dummyImpute()

	r1mismatches = [i for i in range(len(r1)) if r1[i] != wt[i] and r1[i] != "N"]
	r3rcmismatches = [i for i in range(len(r3rc)) if r3rc[i] != wt[219 + i] and r3rc[i] != "N"]


	#eliminate mismatches that have an N in reading frame:
	newr1 = []
	for mis in r1mismatches:
		if mis < len(r1) - 3 and mis > 1: #won't work for edge indices 
			if mis % 3 == 0:
				if r1[mis + 1] != "N" and r1[mis + 2] != "N":
					newr1.append(mis)
			if mis % 3 == 1:
				if r1[mis - 1] != "N" and r1[mis + 1] != "N":
					newr1.append(mis)
			if mis % 3 == 2:
				if r1[mis - 1] != "N" and r1[mis - 2] != "N":
					newr1.append(mis)
	r1mismatches = newr1
	newr3rc = []
	for mis in r3rcmismatches: #check frame for r3rc!!!!!!!!!!!!!! also wobble bp
		if mis < len(r3rc) - 3 and mis > 1: #won't work for edge indices 
			if mis % 3 == 0:
				if r3rc[mis + 1] != "N" and r3rc[mis + 2] != "N":
					newr3rc.append(mis)
			if mis % 3 == 1:
				if r3rc[mis - 1] != "N" and r3rc[mis + 1] != "N":
					newr3rc.append(mis)
			if mis % 3 == 2:
				if r3rc[mis - 1] != "N" and r3rc[mis - 2] != "N":
					newr3rc.append(mis)
	r3rcmismatches = newr3rc

	print("AAA",r1mismatches, r3rcmismatches)
	
	#get most confident mutant index  #bug is here
	confidences = []
	for mis in r1mismatches:
		confidences.append(r1Score[index][mis])
	maxR1Index = np.argmax(confidences)
	maxR1Value = confidences[maxR1Index]
	mutIndexR1 = r1mismatches[maxR1Index]
	# repeat for r3rc
	confidences = []
	for mis in r3rcmismatches:
		confidences.append(r3Score[index][mis])
	maxR3Index = np.argmax(confidences)
	maxR3Value = confidences[maxR3Index]
	mutIndexR3 = r3rcmismatches[maxR3Index]



	if maxR1Value > maxR3Value:
		if maxR1Index % 3 == 0:
			s1 = r1[mutIndexR1] + r1[mutIndexR1 + 1] + r1[mutIndexR1 + 2]
			s1 = s1.replace("T", "U")
			# print(mutIndexR1, translate[s1])
			return (mutIndexR1, s1, translate[s1])
		if maxR1Index % 3 == 1:
			s1 = r1[mutIndexR1 - 1] + r1[mutIndexR1] + r1[mutIndexR1 + 1]
			s1 = s1.replace("T", "U")
			return (mutIndexR1, s1, translate[s1])
		if maxR1Index % 3 == 2:
			s1 = r1[mutIndexR1 - 2] + r1[mutIndexR1 - 1] + r1[mutIndexR1]
			s1 = s1.replace("T", "U")
			return (mutIndexR1, s1, translate[s1])
	
	else: #check frame for r3rc!!!!!!!!!!!!!!
		if mutIndexR3 % 3 == 0:
			s1 = r3rc[mutIndexR3] + r3rc[mutIndexR3 + 1] + r3rc[mutIndexR3 + 2]
			s1 = s1.replace("T", "U")
			return (mutIndexR3, s1, translate[s1])
		if mutIndexR3 % 3 == 1:
			s1 = r3rc[mutIndexR3 - 1] + r3rc[mutIndexR3] + r3rc[mutIndexR3 + 1]
			s1 = s1.replace("T", "U")
			return (mutIndexR3, s1, translate[s1])
		if mutIndexR3 % 3 == 2:
			s1 = r3rc[mutIndexR3 - 2] + r3rc[mutIndexR3 - 1] + r3rc[mutIndexR3]
			s1 = s1.replace("T", "U")
			return (mutIndexR3, s1, translate[s1])
		
	if len(r1mismatches) == 0 and len(r3rcmismatches) == 0: #no apparent mutation outside of the N-cases
		print("thrown out")
		return -1

	#def dummyImpute():
		##append r1 to returnSeq
		# for i in range(0, len(r1)):
		# 	if r1[i] == "N":
		# 		returnSeq += wt[i]
		# 	else:
		# 		returnSeq += r1[i]

		# r1overlap = r1[219:300] #299 inclusive # a substring of r1 representing the overlapping region of r1 with the other read r3rc
		
		# #build overlap region, if both are identical and non-N we take that, 
		# #if the other read has a valid nucleotide, we take the valid one
		# #if both are N, we input from the wt sequence
		# for i in range(0, len(r1overlap)):
		# 	if (r1overlap[i] == r3rc[i]) and r1overlap[i] != "N":
		# 		returnSeq += r1overlap[i]
		# 	elif r1overlap[i] == "N" and r3rc[i] != "N":
		# 		returnSeq += r3rc[i]
		# 	elif r1overlap[i] != "N" and r3rc[i] == "N":
		# 		returnSeq += r1overlap[i]
		# 	else:
		# 		returnSeq += wt[219 + i] 

		# #append r3rc to returnSeq
		# for i in range(0, len(r3rc)):
		# 	if r3rc[i] == "N":
		# 		returnSeq += wt[219 + i]
		# 	else:
		# 		returnSeq += r3rc[i]

		# for i in range(0, len(returnSeq)):
		# 	if returnSeq[i] == "N":
		# 		print(i)
		# return returnSeq


	# mismatches = [i for i in range(len(wt)) if returnSeq[i] != wt[i]]
	# print("msi", len(mismatches), mismatches)

	# return returnSeq


#COMMANDS TO RUN THE SCRIPT===================================

# merge (R1List[1], R3List[1])
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
barcodes = {}
for i in range(0, len(R1List)):
	if R2List[i].count("N") > 15:
		continue
	merged = merge(R1List[i], R3List[i], i)
	if merged == -1:
		continue
	else: #residue, codonMut, aminoAcid, num]
		barcodes[R2List[i]] = [merged[0], merged[1], merged[2], aminoNum[merged[2]]]
print(barcodes) 

pickle.dump( barcodes, open( "barcodes.p", "wb" ) )
	




# for record in SeqIO.parse(file, "fastq"):
	# 	scores = record.letter_annotations["phred_quality"]
	# 	median = 0#np.median(scores)
	# 	l1.append(scores)
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






