main: main.o 
	g++ -o main main.o
generate: generate.o
	g++ -o generate generate.o
compare: compare.o
	g++ -o compare compare.o
print: print.o
	g++ -o print print.o


generate.o: generate.cpp
	g++ -c -o generate.o generate.cpp
main.o: main.cpp
	g++ -c -o main.o main.cpp
compare.o: compare.o
	g++ -c -o compare.o compare.cpp
print.o: print.cpp
	g++ -c -o print.o print.cpp

test: compare main
	./main "./test/A.dat" "./test/B.dat" "./test/C0.dat" 0
	./main "./test/A.dat" "./test/B.dat" "./test/C1.dat" 1
	./main "./test/A.dat" "./test/B.dat" "./test/C2.dat" 2
	./main "./test/A.dat" "./test/B.dat" "./test/C3.dat" 3
	./main "./test/A.dat" "./test/B.dat" "./test/C4.dat" 4
	./main "./test/A.dat" "./test/B.dat" "./test/C5.dat" 5
	./compare "./test/C0.dat" "./test/C1.dat"
	./compare "./test/C0.dat" "./test/C2.dat"
	./compare "./test/C0.dat" "./test/C3.dat"
	./compare "./test/C0.dat" "./test/C4.dat"
	./compare "./test/C0.dat" "./test/C5.dat"
	rm ./test/C*.dat


clean : 
	rm *.dat *.txt *.o