include ../Makefile.lumi.in

OBJS=tests_mod.o run_tests.o

all: run_tests

run_tests: $(OBJS)
	$(FC) $(FOPTS) $(INCS) -I ../src -o run_tests $(OBJS) -L../src -leztrans $(LIBS)

%.o: %.F90
	$(FC) $(FOPTS) $(INCS) -I ../src -c $?

