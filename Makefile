all: a.out

a.out: nlo_main.o nlo.o nlo_functions.o
	gcc -L/usr/local/lib nlo_main.o nlo.o nlo_functions.o -lgsl -lgslcblas -lm
	
nlo_main.o: nlo_main.c
	gcc -Wall -I/usr/local/include -c nlo_main.c

nlo.o: nlo.c
	gcc -O0 -Wall -I/usr/local/include -c nlo.c

nlo_functions.o: nlo_functions.c
	gcc -O2 -Wall -I/usr/local/include -c nlo_functions.c
	
clean:
	rm -rf *.o