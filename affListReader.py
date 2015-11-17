file = open("./affParam.txt", "r")
header1 = file.readline()
lines = file.readlines()
file.close()

element = "O1-"

elementList = []
a1V = []
b1V = []
a2V = []
b2V = []
a3V = []
b3V = []
a4V = []
b4V = []
cV = []

for line in lines:
    columns = line.split()
    if columns[0] == element:
        a1 = float(columns[1])
        b1 = float(columns[2])
        a2 = float(columns[3])
        b2 = float(columns[4])
        a3 = float(columns[5])
        b3 = float(columns[6])
        a4 = float(columns[7])
        b4 = float(columns[8])
        c = float(columns[9])
    
print(a1, b1, a2, b2, a3, b3, a4, b4, c)    
    
    
	# elementList.append(columns[0])
	# a1V.append(float(columns[1]))
	# b1V.append(float(columns[2]))
	# a2V.append(float(columns[3]))
	# b2V.append(float(columns[4]))
	# a3V.append(float(columns[5]))
	# b3V.append(float(columns[6]))
	# a4V.append(float(columns[7]))
	# b4V.append(float(columns[8]))
	# cV.append(float(columns[9]))

# for idx, elem in enumerate(elementList):
	# if elem == element:
		# indexElem = idx
	

# print(indexElem, a1V[indexElem], b1V[indexElem], a2V[indexElem], b2V[indexElem], a3V[indexElem], b3V[indexElem], a4V[indexElem], b4V[indexElem], cV[indexElem])