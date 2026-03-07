def gera_particoes(N, tam_part, num_part_fix):
	
	tam_janela = tam_part - num_part_fix
	subset = []
	 
	for i in range(0, N, tam_janela):
		if i + tam_part > N:
			subset.append([k for k in range(i,N)])
		else:
			subset.append([k for k in range(i,i+tam_part)])

	return subset