.PHONY: test plot


main: main.o 
	g++ -o main -I/usr/local/include main.o /usr/local/lib/libpapi.a
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

compare.o: compare.cpp
	g++ -c -o compare.o compare.cpp
print.o: print.cpp
	g++ -c -o print.o print.cpp

test: compare main 
	./main "./test/A.dat" "./test/B.dat" "./test/C0.dat" 0
	./main "./test/A.dat" "./test/B.dat" "./test/C1.dat" 1
	./main "./test/A.dat" "./test/B.dat" "./test/C2.dat" 2
	./compare "./test/C0.dat" "./test/C.dat"
	./compare "./test/C1.dat" "./test/C.dat"
	./compare "./test/C2.dat" "./test/C.dat"
	rm -f ./test/C0.dat ./test/C1.dat ./test/C2.dat


report: generate main
	./generate f 1000 1000 A.dat
	./generate f 1000 1000 B.dat
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2

plot :
	gnuplot plotset
	



clean : 
	rm -f *.dat *.txt *.o *.svg  print compare generate main a.out