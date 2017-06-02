CFLAGS=-c -Wall -I/usr/local/include
LDFLAGS=-L/usr/local/lib 
LIBS=-lgsl -lgslcblas -lm 

all: a.out

a.out: main.o nlo2.o B1.o B2.o B3.o T.o R.o IR.o C.o nlo_functions.o
	gcc $(LDFLAGS) main.o nlo2.o B1.o B2.o B3.o T.o R.o IR.o C.o nlo_functions.o  $(LIBS)
	

main.o:
	gcc $(CFLAGS) main.c


nlo2.o:
	gcc $(CFLAGS) nlo2.c


B1.o:
	gcc $(CFLAGS) B1.c


B2.o:
	gcc $(CFLAGS) B2.c


B3.o:
	gcc $(CFLAGS) B3.c
	

T.o:
	gcc $(CFLAGS) T.c
	

R.o:
	gcc $(CFLAGS) R.c
	

IR.o:
	gcc $(CFLAGS) IR.c
	

C.o:
	gcc $(CFLAGS) C.c
	
nlo_functions.o:
	gcc $(CFLAGS) nlo_functions.c