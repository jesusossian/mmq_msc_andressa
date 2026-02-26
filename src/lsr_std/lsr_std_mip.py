from pathlib import Path
import os
import leitura as ler
import optimization as opt
import relax_fix as rf
import fix_optimize as fop
import gera_particoes as gera
import numpy as np
import pandas as pd
import sys
from datetime import datetime, date
import time
import itertools

# parameters
file_name = sys.argv[1]

# paths
result_path = Path('result/')
solution_path = Path('solution/')
instance_path = Path('../../data/sifaleras')

def main():
    
    # read instances
	N, PP, PR, FP, FR, HP, HR, D, R = ler.leitura_instance(os.path.join(instance_path,file_name))
	
	SD = (np.zeros((N,N))).tolist()
	SR = (np.zeros((N,N))).tolist()
	for  i in range(N):
		SD[i][i] = D[i]
		SR[i][i] = R[i]
		for j in range(i+1, N):
			SD[i][j] = SD[i][j-1] + D[j]
			SR[i][j] = SR[i][j-1] + R[j]
		
	objval, objbound, mgap, rtime, ncount, tmp = opt.lsr_std_mip(N, PP, PR, FP, FR, HP, HR, D, R, SD, SR)
	#arquivo = open(os.path.join(result_path,'lsr_std_mip.txt'),'a')
	#arquivo.write(file_name+';'
	# 	+str(round(objval,2))+';'
	# 	+str(round(objbound,2))+';'
	# 	+str(round(mgap,2))+';'
	# 	+str(round(rtime,2))+';'
	# 	+str(round(ncount,2))+';'
	# 	+str(round(tmp,2))
	# 	+'\n')
	#arquivo.close()

if __name__== "__main__" :
	main()
