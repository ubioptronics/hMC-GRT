#DEBUGFLAGS=-ggdb -pg 
#DEBUGFLAGS=-ggdb -pg -a -Q 

FLAGS=-ggdb -Wall
LIBS=-lm
CC=x86_64-w64-mingw32-gcc $(FLAGS)

objchecker = trmc_2fls_tc_checker.o trmc_2fls_tc.o 
objcore = trmc_2fls_core.o trmc_2fls_tc.o 

trmc_2fls_tc_checker : $(objchecker)
	$(CC) -o trmc_2fls_tc_checker $(objchecker) $(LIBS)
trmc_2fls_tc_checker.o: trmc_2fls_globals.h trmc_2fls_base.h
trmc_2fls_tc.o: trmc_2fls_opts.h trmc_2fls_base.h

time: $(objcore)
	$(CC) -o trmc_mlmf $(objcore) $(LIBS)
trmc_2fls_core.o: trmc_2fls_globals.h trmc_2fls_base.h
trmc_2fls_tc.o: trmc_2fls_opts.h trmc_2fls_base.h


.PHONY: clean
clean :
	rm -f $(objcore) $(objchecker)

