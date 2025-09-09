#!/bin/bash
#52 periodos, 108 instancias

for id in {1..1} #$(seq 1)
do
	python3 clsr_math_std_mip.py c52_${id}.txt >> report/out_std_c52_${id}.txt
done
