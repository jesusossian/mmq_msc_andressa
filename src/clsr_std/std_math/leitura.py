import numpy as np

def leitura_instance(file_name):

	arq = open(file_name)

	N = int(arq.readline())

	FR_ = [float(arq.readline())]*N
	FP_ = [float(arq.readline())]*N
	
	HR_ = [float(arq.readline())]*N
	HP_ = [float(arq.readline())]*N

	D_ = [int(i) for i in arq.readline().split()]
	R_ = [int(i) for i in arq.readline().split()]
	C = float(arq.readline().rstrip('\n'))

	#PP = [0]*N
	PP = np.zeros(N)
	#PR  = [0]*N
	PR = np.zeros(N)
	
	FP = np.zeros(N)
	FR = np.zeros(N)
	
	HP = np.zeros(N)
	HR = np.zeros(N)

	D = np.zeros(N)
	R = np.zeros(N)

	for  i in range(N):
		FP[i] = FP_[i]

	for  i in range(N):
		FR[i] = FR_[i]

	for  i in range(N):
		HP[i] = HP_[i]

	for  i in range(N):
		HR[i] = HR_[i]

	for  i in range(N):
		D[i] = D_[i]
	
	for  i in range(N):
		R[i] = R_[i]

	return N, PP, PR, FP, FR, HP, HR, D, R, C
