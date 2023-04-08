all: libeztrans programs test

libeztrans:
	make -C src/

programs: src/libeztrans.a
	make -C programs/

test: src/libeztrans.a
	make -C test/

clean:
	rm -f src/*.mod src/*.o programs/driver.o programs/driver test/run_tests.o test/run_tests
