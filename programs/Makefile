include ../Makefile.in

OBJS=driver.o

all: driver

driver: $(OBJS)
	$(FC) $(FOPTS) $(INCS) -I ../src -o driver $? -L../src -leztrans $(LIBS)

%.o: %.F90
	$(FC) $(FOPTS) $(INCS) -I ../src -c $?

