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

file_name = sys.argv[1]

# paths
result_path = Path('result/')
solution_path = Path('solution/')
instance_path = Path('../../../data/csifa')

def main():

    # read instances
	N, PP, PR, FP, FR, HP, HR, D, R, C = ler.leitura_instance(os.path.join(instance_path,file_name))

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
	tam_partp_rf = 5 # tamanho particoes
	num_fixp_rf = 2 # numero de fixadas
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
		rf_objval, yp_val, yr_val = rf.relax_fix(conjp, conjr, yp_val, yr_val, N, PP, PR, FP, FR, HP, HR, D, R, SD, SR, C)
	rf_rtime = time.time() - start_time
	
	start_time = time.time()
	for conjp, conjr in itertools.zip_longest(subsetp_fo,subsetr_fo):
		fo_objval, yp_val, yr_val = fop.fix_and_optimize(conjp, conjr, yp_val, yr_val, N, PP, PR, FP, FR, HP, HR, D, R, SD, SR, C)
	fo_rtime = time.time() - start_time
	
	arquivo = open(os.path.join(result_path,'clsr_mc_rf_fo.txt'),'a')
	arquivo.write(file_name+';'
		+str(round(rf_objval,2))+';'
		+str(round(rf_rtime,2))+';'
		+str(round(fo_objval,2))+';'
		+str(round(fo_rtime,2))
		+'\n')
	arquivo.close()

	objval, objbound, mgap, rtime, ncount, tmp = opt.clsr_mc(N, PP, PR, FP, FR, HP, HR, D, R, SD, SR, C, yp_val, yr_val)

	arquivo = open(os.path.join(result_path,'clsr_mc_rf_fo_mip.txt'),'a')
	arquivo.write(file_name+';'
	 	+str(round(objval,2))+';'
	 	+str(round(objbound,2))+';'
	 	+str(round(mgap,2))+';'
	 	+str(round(rtime,2))+';'
	 	+str(round(ncount,2))+';'
	 	+str(round(tmp,2))
	 	+'\n')
	arquivo.close()
	
	objval, objbound, mgap, rtime, ncount, tmp = opt.clsr_mc_mip(N, PP, PR, FP, FR, HP, HR, D, R, SD, SR, C)
	
	arquivo = open(os.path.join(result_path,'clsr_mc_mip.txt'),'a')
	arquivo.write(file_name+';'
	 	+str(round(objval,2))+';'
	 	+str(round(objbound,2))+';'
	 	+str(round(mgap,2))+';'
	 	+str(round(rtime,2))+';'
	 	+str(round(ncount,2))+';'
	 	+str(round(tmp,2))
	 	+'\n')
	arquivo.close()

	objval, rtime = opt.clsr_mc_lp(N, PP, PR, FP, FR, HP, HR, D, R, SD, SR, C)
	
	arquivo = open(os.path.join(result_path,'clsr_mc_lp.txt'),'a')
	arquivo.write(file_name+';'
	 	+str(round(objval,2))+';'
	 	+str(round(rtime,2))+';'
	 	+'\n')
	arquivo.close()

if __name__== "__main__" :
	main()
