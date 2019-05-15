
.PHONY: test


main: main.cpp logic_gate.h
	mpicxx -std=c++98 -o main main.cpp -fopenmp


tobin: tobin.cpp
	g++ -o tobin tobin.cpp

totxt: totxt.cpp
	g++ -o totxt totxt.cpp

compare: compare.cpp
	g++ -o compare compare.cpp

clean : 
	rm -f *.dat *.txt  *.svg print generate main random compare *.map totxt a.out tobin hadamar N_hadamar Not Rw C_not C_rw *.bin

test: main totxt compare
	mpirun -n 4 ./main 3 1 1 ./test/in.bin out.bin
	./compare out.bin ./test/out.bin
	./totxt out.bin out.txt