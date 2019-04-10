main: main.cpp
	mpicxx -std=c++98 -o main main.cpp -fopenmp


tobin: tobin.cpp
	g++ -o tobin tobin.cpp

totxt: totxt.cpp
	g++ -o totxt totxt.cpp


clean : 
	rm -f *.dat *.txt  *.svg print generate main random compare *.map totxt a.out *.bin tobin

test: main totxt
	mpirun -n 1 ./main 2 2 1 2 2
	./totxt out.bin out.txt
