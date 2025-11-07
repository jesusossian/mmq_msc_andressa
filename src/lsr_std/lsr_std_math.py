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
	
	yp_val = np.zeros(N)
	yr_val = np.zeros(N)

	SD = (np.zeros((N,N))).tolist()
	SR = (np.zeros((N,N))).tolist()
	for  i in range(N):
		SD[i][i] = D[i]
		SR[i][i] = R[i]
		for j in range(i+1, N):
			SD[i][j] = SD[i][j-1] + D[j]
			SR[i][j] = SR[i][j-1] + R[j]
	
	# parameters rf
	tam_partp_rf = 5
	num_fixp_rf = 2
	tam_partr_rf = 6
	num_fixr_rf = 3
	
	# parameters fo
	tam_partp_fo = 7
	num_fixp_fo = 2
	tam_partr_fo = 8
	num_fixr_fo = 3

	subsetp_rf = gera.gera_particoes(N, tam_partp_rf, num_fixp_rf)
	subsetr_rf = gera.gera_particoes(N, tam_partr_rf, num_fixr_rf)
	subsetp_fo = gera.gera_particoes(N, tam_partp_fo, num_fixp_fo)
	subsetr_fo = gera.gera_particoes(N, tam_partr_fo, num_fixr_fo)
	
	start_time = time.time()
	for conjp, conjr in itertools.zip_longest(subsetp_rf,subsetr_rf): #zip
		rf_objval, yp_val, yr_val = rf.relax_fix(conjp, conjr, yp_val, yr_val, N, PP, PR, FP, FR, HP, HR, D, R, SD, SR)
	rf_rtime = time.time() - start_time
	
	start_time = time.time()
	for conjp, conjr in itertools.zip_longest(subsetp_fo,subsetr_fo):
		fop_objval, yp_val, yr_val = fop.fix_and_optimize(conjp, conjr, yp_val, yr_val, N, PP, PR, FP, FR, HP, HR, D, R, SD, SR)
	fop_rtime = time.time() - start_time
	
	arquivo = open(os.path.join(result_path,'lsr_std_math.txt'),'a')
	arquivo.write(file_name+';'
		+str(round(rf_objval,2))+';'
		+str(round(rf_rtime,2))+';'
		+str(round(fop_objval,2))+';'
		+str(round(fop_rtime,2))
		+'\n')
	arquivo.close()

if __name__== "__main__" :
	main()
