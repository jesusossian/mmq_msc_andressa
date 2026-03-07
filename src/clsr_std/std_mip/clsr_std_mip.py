from pathlib import Path
import os
import leitura as ler
import optimization as opt
import numpy as np
import pandas as pd
import sys
from datetime import datetime, date

file_name = sys.argv[1]

# path
result_path   = Path('result/')
instance_path = Path('../../../data/c1sifa')

# global variables
objv = 0
objb = 0
nnode = 0
rtime = 0
gap = 0

def main():

	N, PP, PR, FP, FR, HR, HP, D, R, C = ler.leitura_instance(os.path.join(instance_path,file_name))

	SD = (np.zeros((N,N))).tolist()
	SR = (np.zeros((N,N))).tolist()
	for  i in range(N):
		SD[i][i] = D[i]
		SR[i][i] = R[i]
		for j in range(i+1, N):
			SD[i][j] = SD[i][j-1] + D[j]
			SR[i][j] = SR[i][j-1] + R[j]
	
	objv, objb, gap, rtime, nnode, status = opt.clsr_std(N,PP,PR,FP,FR,HR,HP,D,R,SD,SR,C)
	
	arquivo = open(os.path.join(result_path,'clsr_std_mip_first_feasible.txt'),'a')
	arquivo.write(file_name+';'+str(round(objv,2))+';'+str(round(objb,2))+ \
					';'+str(round(gap,2))+';'+str(round(rtime,2))+';'+str(round(nnode,2))+ \
					';'+str(round(status,2))+
					'\n')
	arquivo.close()

if __name__== "__main__" :
	main()
