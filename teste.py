import os

fileName = ['imagem1e-1.pgm', 'imagem1e-2.pgm', 'imagem1e-3.pgm', 'imagem1e-4.pgm','imagem1e-5.pgm']
tol      = ['1e-1', '1e-2', '1e-3', '1e-4', '1e-5']



for x in range(5):
	os.system("./fuck " + fileName[x] + " " + tol[x]) 	



for x in range(1, 100):
	passo = x*0.5
	nome = "imagem" + str(passo) + ".pgm" 
	commando = "./fuck2 " + nome + " " + str(passo)
	print(commando)
	os.system(commando)


