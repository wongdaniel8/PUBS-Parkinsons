from Bio import SeqIO
R1List = []
R2List = []
R3List = []
files = ["R1", "R2", "R3"]
for file in files:
	if file == "R1":
		l = R1List
	if file == "R2":
		l = R2List
	if file == "R3":
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

# print(R1List)
# print(R2List)
# print(R3List)

# print(len(R1List))
# print(len(R2List))
# print(len(R3List))

# for i in range(0, len(R1List)):
	# print(R1List[i])
	# print(R3List[i])
