#!/bin/bash
#52 periodos, 108 instancias

<<<<<<< HEAD
for id in {1..1} #$(seq 1)
=======
for id in {1..2} #$(seq 1)
>>>>>>> 2c7fcbaa55dac7d8380557019eee6fe94ea83a5d
do
	python3 lsr_std_lp.py 52_${id}.txt >> report/out_std_52_${id}.txt
done
