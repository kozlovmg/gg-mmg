CFLAGS=-c -Wall  -I/usr/local/include
LDFLAGS=-L/usr/local/lib 
LIBS=-lgsl -lgslcblas -lm 

all: a.out

a.out: main.o nlo2.o B1.o B2.o B3.o T.o R.o IR.o Csum.o Cpart1.o Cpart2.o Cpart3.o Cpart4.o Cpart5.o Cpart6.o nlo_functions.o
	gcc $(LDFLAGS) main.o nlo2.o B1.o B2.o B3.o T.o R.o IR.o Csum.o Cpart1.o Cpart2.o Cpart3.o Cpart4.o Cpart5.o Cpart6.o nlo_functions.o  $(LIBS)
	

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

Cpart1.o:
	gcc $(CFLAGS) Cpart1.c
	
Cpart2.o:
	gcc $(CFLAGS) Cpart2.c
	
Cpart3.o:
	gcc $(CFLAGS) Cpart3.c
	
Cpart4.o:
	gcc $(CFLAGS) Cpart4.c
	
Cpart5.o:
	gcc $(CFLAGS) Cpart5.c
	
Cpart6.o:
	gcc $(CFLAGS) Cpart6.c
	
Csum.o:
	gcc $(CFLAGS) Csum.c
	
nlo_functions.o:
	gcc $(CFLAGS) nlo_functions.c