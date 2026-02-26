#!/bin/bash
#52 periodos, 108 instancias

for id in {97..108} #$(seq 1)
do
	python3 lsr_std_mip.py 52_${id}.txt >> report/out_std_mip_52_${id}.txt
done
