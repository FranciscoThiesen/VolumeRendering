import os
sum = 199
for x in range(1,100):
	creu = "mv " + "imagem" + str(int(100+x-1))+ ".pgm " + " imagem" + str(int(100+x-1 + sum)) + ".pgm"
	sum = sum - 2
	#print(creu)
	os.system(creu)

