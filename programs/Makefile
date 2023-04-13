include ../Makefile.lumi.in

OBJS=driver.o

all: driver

driver: $(OBJS)
	$(FC) $(FOPTS) $(INCS) -I ../src -o driver $? -L../src -leztrans $(LIBS)

driver.o: ../src/libeztrans.a driver.F90
	$(FC) $(FOPTS) $(INCS) -I ../src -c driver.F90

