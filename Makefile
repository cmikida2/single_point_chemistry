#----------------------------------------------------------------------------
#          Makefile for PyJac Single-Point Example Programs
#----------------------------------------------------------------------------


COMPILER = g++

CANTERA_INC = /home/cmikida2/Work/Fall_20/cantera_install/include
CANTERA_LIB = /home/cmikida2/Work/Fall_20/cantera_install/lib

OPTS = -I$(CANTERA_INC)
HDRS = $(CANTERA_INC)

all: test_rk4 test_leap_imex test_leap_mr_imex test_trap

# RK4 with/without Cantera
test_rk4: test_rk4.o species_rates.hpp jacobian.hpp
	$(COMPILER) -L$(CANTERA_LIB) -o test_rk4 test_rk4.o jacobian.cpp species_rates.cpp chem_utils.cpp -pthread -lm -lcantera -llapack -lblas

test_rk4.o: test_rk4.cpp species_rates.cpp jacobian.cpp chem_utils.cpp $(HDRS)
	$(COMPILER) $(OPTS) -g -c test_rk4.cpp

# Leap IMEX with/without Cantera
test_leap_imex: test_leap_imex.o species_rates.hpp jacobian.hpp LeapIMEXMethod.H
	$(COMPILER) -L$(CANTERA_LIB) -o test_leap_imex test_leap_imex.o jacobian.cpp species_rates.cpp chem_utils.cpp -pthread -lm -lcantera -llapack -lblas

test_leap_imex.o: test_leap_imex.cpp species_rates.cpp jacobian.cpp chem_utils.cpp LeapIMEXMethod.H $(HDRS)
	$(COMPILER) $(OPTS) -g -c test_leap_imex.cpp

# Leap IMEX MR with/without Cantera
test_leap_mr_imex: test_leap_mr_imex.o species_rates.hpp jacobian.hpp LeapMRIMEXMethod.H
	$(COMPILER) -L$(CANTERA_LIB) -o test_leap_mr_imex test_leap_mr_imex.o jacobian.cpp species_rates.cpp chem_utils.cpp -pthread -lm -lcantera -llapack -lblas

test_leap_mr_imex.o: test_leap_mr_imex.cpp species_rates.cpp jacobian.cpp chem_utils.cpp LeapMRIMEXMethod.H $(HDRS)
	$(COMPILER) $(OPTS) -g -c test_leap_mr_imex.cpp

test_trap: test_trap.cpp jacobian.cpp species_rates.cpp chem_utils.cpp
	$(COMPILER) -g -o test_trap test_trap.cpp jacobian.cpp species_rates.cpp chem_utils.cpp -llapack -lblas 
	
