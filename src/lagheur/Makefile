# Copyright 2020, Gurobi Optimization, LLC

PLATFORM = linux64
INC      = /opt/gurobi1201/linux64/include/
CLIB     = -L/opt/gurobi1201/linux64/lib -lgurobi120
CPPLIB   = -L/opt/gurobi1201/linux64/lib -lgurobi_c++ $(CLIB)



all: modelo

run: run_modelo

modelo:main.cpp
	g++ -m64 -g -std=c++11 main.cpp -o out -I$(INC) $(CPPLIB) -lpthread -lm

run_modelo: modelo
	modelo

# ------------------------------------------------------------

clean :
	/bin/rm -rf *.o *~
	/bin/rm -rf out
	/bin/rm -rf *.mps *.ord *.sos *.lp *.sav *.net *.msg *.log *.clp

# ------------------------------------------------------------

# Local Variables:
# mode: makefile
# End:
