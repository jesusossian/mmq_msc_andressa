#!/bin/bash
#52 periodos, 108 instancias

for id in {1..2} #$(seq 1)
do
	python3 lsr_std_math.py 52_${id}.txt >> report/out_std_52_${id}.txt
done
