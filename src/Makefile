include ../Makefile.lumi.in

OBJS=aux_mod.o fftw3_mod.o config_mod.o trgtol_mod.o trltog_mod.o trltom_mod.o trmtol_mod.o trmtos_mod.o trstom_mod.o edir_trans_mod.o einv_trans_mod.o eztrans_setup.o eztrans_end.o parkind1_dummy.o set99b.o fft992.o

all: libeztrans.a

libeztrans.a: $(OBJS)
	ar -rc libeztrans.a *.o
	
aux_mod.mod: aux_mod.o
fftw3_mod.mod: fftw3_mod.o
config_mod.mod: aux_mod.mod fftw3_mod.mod config_mod.o
trgtol_mod.mod: config_mod.mod trgtol_mod.o
trltog_mod.mod: config_mod.mod trltog_mod.o
trltom_mod.mod: config_mod.mod trltom_mod.o
trmtol_mod.mod: config_mod.mod trmtol_mod.o
trmtos_mod.mod: config_mod.mod trmtos_mod.o
trstom_mod.mod: config_mod.mod trstom_mod.o
edir_trans_mod.mod: config_mod.mod trgtol_mod.mod trltom_mod.mod edir_trans_mod.o
%.o: %.F90
	$(FC) $(FOPTS) $(INCS) -c $?

