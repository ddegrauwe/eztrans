all: libeztrans programs test

libeztrans:
	make -C src/

programs: libeztrans
	make -C programs/

test: libeztrans
	make -C test/

clean:
	rm -f src/*.mod src/*.o programs/driver.o programs/driver test/run_tests.o test/run_tests
