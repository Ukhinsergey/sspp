main: main.cpp
	mpicxx -std=c++98 -O3 -o main main.cpp
generate: generate.cpp
	g++ -O3 -o generate generate.cpp
print: print.cpp
	g++ -O3 -o print print.cpp
random: random.cpp
	g++ -O3 -o random random.cpp
compare: compare.cpp
	g++ -O3 -o compare compare.cpp


genmatr: generate
	./generate d 10 10 A.dat
	./generate d 10 10 B.dat
printmatr: print
	./print A.dat A.txt
	./print B.dat B.txt

test: main compare 
	mpirun -np 8 main "./../test/A.dat" "./../test/B.dat" "./../test/C0.dat" 
	mpirun -np 27 main "./../test/A.dat" "./../test/B.dat" "./../test/C1.dat" 
	./compare "./../test/C0.dat" "./../test/C.dat"
	./compare "./../test/C1.dat" "./../test/C.dat"
clean : 
	rm -f *.dat *.txt  *.svg print generate main random *.map