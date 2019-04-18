
.PHONY: test


main: main.cpp logic_gate.h
	mpicxx -std=c++98 -o main main.cpp -fopenmp


tobin: tobin.cpp
	g++ -o tobin tobin.cpp

totxt: totxt.cpp
	g++ -o totxt totxt.cpp

hadamar: hadamar.cpp logic_gate.h
	mpicxx -std=c++98 -o hadamar hadamar.cpp -fopenmp

N_hadamar: N_hadamar.cpp logic_gate.h
	mpicxx -std=c++98 -o N_hadamar N_hadamar.cpp -fopenmp

Not: Not.cpp logic_gate.h
	mpicxx -std=c++98 -o Not Not.cpp -fopenmp

C_not: C_not.cpp logic_gate.h
	mpicxx -std=c++98 -o C_not C_not.cpp -fopenmp

Rw: Rw.cpp logic_gate.h
	mpicxx -std=c++98 -o Rw Rw.cpp -fopenmp


C_rw: C_rw.cpp logic_gate.h
	mpicxx -std=c++98 -o C_rw C_rw.cpp -fopenmp

compare: compare.cpp
	g++ -o compare compare.cpp

testhadamar: hadamar compare
	mpirun -n 1 ./hadamar 6 1 1 2 1 hadamar.bin
	./compare test/hadamar.bin hadamar.bin
	mpirun -n 2 ./hadamar 6 1 1 2 1 hadamar.bin
	./compare test/hadamar.bin hadamar.bin
	mpirun -n 4 ./hadamar 6 1 1 2 1 hadamar.bin
	./compare test/hadamar.bin hadamar.bin
	mpirun -n 4 ./hadamar 6 1 1 2 4 hadamar.bin
	./compare test/hadamar.bin hadamar.bin
	mpirun -n 4 ./hadamar 6 1 1 2 8 hadamar.bin
	./compare test/hadamar.bin hadamar.bin
	mpirun -n 8 ./hadamar 6 1 1 2 1 hadamar.bin
	./compare test/hadamar.bin hadamar.bin
	rm hadamar.bin

testN_hadamar: N_hadamar compare
	mpirun -n 1 ./N_hadamar 6 1 1 2 1 N_hadamar.bin
	./compare test/N_hadamar.bin N_hadamar.bin
	mpirun -n 2 ./N_hadamar 6 1 1 2 1 N_hadamar.bin
	./compare test/N_hadamar.bin N_hadamar.bin
	mpirun -n 4 ./N_hadamar 6 1 1 2 1 N_hadamar.bin
	./compare test/N_hadamar.bin N_hadamar.bin
	mpirun -n 4 ./N_hadamar 6 1 1 2 4 N_hadamar.bin
	./compare test/N_hadamar.bin N_hadamar.bin
	mpirun -n 4 ./N_hadamar 6 1 1 2 8 N_hadamar.bin
	./compare test/N_hadamar.bin N_hadamar.bin
	mpirun -n 8 ./N_hadamar 6 1 1 2 1 N_hadamar.bin
	./compare test/N_hadamar.bin N_hadamar.bin
	rm N_hadamar.bin

testNot: Not compare
	mpirun -n 1 ./Not 6 1 1 2 1 Not.bin
	./compare test/Not.bin Not.bin
	mpirun -n 2 ./Not 6 1 1 2 1 Not.bin
	./compare test/Not.bin Not.bin
	mpirun -n 4 ./Not 6 1 1 2 1 Not.bin
	./compare test/Not.bin Not.bin
	mpirun -n 4 ./Not 6 1 1 2 4 Not.bin
	./compare test/Not.bin Not.bin
	mpirun -n 4 ./Not 6 1 1 2 8 Not.bin
	./compare test/Not.bin Not.bin
	mpirun -n 8 ./Not 6 1 1 2 1 Not.bin
	./compare test/Not.bin Not.bin
	rm Not.bin

testRw: Rw compare
	mpirun -n 1 ./Rw 6 1 1 2 1 Rw.bin 3.1415
	./compare test/Rw.bin Rw.bin
	mpirun -n 2 ./Rw 6 1 1 2 1 Rw.bin 3.1415
	./compare test/Rw.bin Rw.bin
	mpirun -n 4 ./Rw 6 1 1 2 1 Rw.bin 3.1415
	./compare test/Rw.bin Rw.bin
	mpirun -n 4 ./Rw 6 1 1 2 4 Rw.bin 3.1415
	./compare test/Rw.bin Rw.bin
	mpirun -n 4 ./Rw 6 1 1 2 8 Rw.bin 3.1415
	./compare test/Rw.bin Rw.bin
	mpirun -n 8 ./Rw 6 1 1 2 1 Rw.bin 3.1415
	./compare test/Rw.bin Rw.bin
	rm Rw.bin

testC_not: C_not compare
	mpirun -n 1 ./C_not 6 6 1  2 1 C_not.bin
	./compare test/C_not.bin C_not.bin
	mpirun -n 2 ./C_not 6 6 1  2 1 C_not.bin
	./compare test/C_not.bin C_not.bin
	mpirun -n 4 ./C_not 6 6 1  2 1 C_not.bin
	./compare test/C_not.bin C_not.bin
	mpirun -n 4 ./C_not 6 6 1  2 4 C_not.bin
	./compare test/C_not.bin C_not.bin
	mpirun -n 4 ./C_not 6 6 1 2 8 C_not.bin
	./compare test/C_not.bin C_not.bin
	mpirun -n 8 ./C_not 6 6 1  2 1 C_not.bin
	./compare test/C_not.bin C_not.bin
	mpirun -n 64 ./C_not 6 6 1  2 1 C_not.bin
	./compare test/C_not.bin C_not.bin
	rm C_not.bin

testC_rw: C_rw compare
	mpirun -n 1 ./C_rw 6 1 6 2 1 C_rw.bin 3.1415
	./compare test/C_rw.bin C_rw.bin
	mpirun -n 2 ./C_rw 6 1 6 2 1 C_rw.bin 3.1415
	./compare test/C_rw.bin C_rw.bin
	mpirun -n 4 ./C_rw 6 1 6 2 1 C_rw.bin 3.1415
	./compare test/C_rw.bin C_rw.bin
	mpirun -n 4 ./C_rw 6 1 6 2 4 C_rw.bin 3.1415
	./compare test/C_rw.bin C_rw.bin
	mpirun -n 4 ./C_rw 6 1 6 2 8 C_rw.bin 3.1415
	./compare test/C_rw.bin C_rw.bin
	mpirun -n 8 ./C_rw 6 1 6 2 1 C_rw.bin 3.1415
	./compare test/C_rw.bin C_rw.bin
	mpirun -n 64 ./C_rw 6 1 6 2 1 C_rw.bin 3.1415
	./compare test/C_rw.bin C_rw.bin
	rm C_rw.bin

test: testhadamar testN_hadamar testNot testRw testC_not testC_rw


lint: *.cpp
	cpplint logic_gate.h

clean : 
	rm -f *.dat *.txt  *.svg print generate main random compare *.map totxt a.out *.bin tobin hadamar N_hadamar Not Rw C_not C_rw

