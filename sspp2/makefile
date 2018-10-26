.PHONY: test plot


main: main.cpp 
	g++ -std=c++11 -o main -I/usr/local/include main.cpp -lpapi
generate: generate.cpp
	g++ -std=c++11 -o generate generate.cpp
compare: compare.cpp
	g++ -std=c++11 -o compare compare.cpp
print: print.cpp
	g++ -std=c++11 -o print print.cpp
main1: main1.cpp
	g++ -std=c++11 -o main1 -I/usr/local/include main1.cpp -lpapi

	
test: compare main 
	./main "./test/A.dat" "./test/B.dat" "./test/C0.dat" 0 0
	./main "./test/A.dat" "./test/B.dat" "./test/C1.dat" 1 0
	./main "./test/A.dat" "./test/B.dat" "./test/C2.dat" 2 0
	./main "./test/A.dat" "./test/B.dat" "./test/C0.dat" 0 1
	./main "./test/A.dat" "./test/B.dat" "./test/C1.dat" 1 1
	./main "./test/A.dat" "./test/B.dat" "./test/C2.dat" 2 1
	./compare "./test/C0.dat" "./test/C.dat"
	./compare "./test/C1.dat" "./test/C.dat"
	./compare "./test/C2.dat" "./test/C.dat"
	rm -f ./test/C0.dat ./test/C1.dat ./test/C2.dat


report: generate main
	for number in 1000 2000 3000 4000 5000; do \
	./generate f $$number $$number A.dat ;\
	./generate f $$number $$number B.dat ;\
	./main A.dat B.dat C.dat 0 ;\
	./main A.dat B.dat C.dat 1 ;\
	./main A.dat B.dat C.dat 2 ;\
	done

easdsad: 
	./generate f 1000 1000 A.dat
	./generate f 1000 1000 B.dat
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./generate f 2000 2000 A.dat
	./generate f 2000 2000 B.dat
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./generate f 3000 2000 A.dat
	./generate f 2000 2000 B.dat
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./generate f 2000 2000 A.dat
	./generate f 2000 2000 B.dat
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./generate f 2000 2000 A.dat
	./generate f 2000 2000 B.dat
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
plot :
	gnuplot plotset
	



clean : 
	rm -f *.dat *.txt *.o *.svg  print compare generate main a.out
