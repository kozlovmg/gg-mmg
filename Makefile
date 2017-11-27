CFLAGS=-c -Wall  -I/usr/local/include
LDFLAGS=-L/usr/local/lib 
LIBS=-lgsl -lgslcblas -lm 

all: a.out 


main.o: main.c
	gcc $(CFLAGS) main.c
	
born.o: born.c
	gcc $(CFLAGS) born.c

a.out: main.o nlo2.o B1.o B2.o B3.o T.o R.o IR.o Csum.o Cpart1.o Cpart2.o Cpart3.o Cpart4.o Cpart5.o Cpart6.o Cpart7.o nlo_functions.o born.o
	gcc $(LDFLAGS) main.o nlo2.o B1.o B2.o B3.o T.o R.o IR.o Csum.o Cpart1.o Cpart2.o Cpart3.o Cpart4.o Cpart5.o Cpart6.o Cpart7.o nlo_functions.o born.o $(LIBS)
	


nlo2.o: nlo2.c
	gcc $(CFLAGS) nlo2.c


B1.o: B1.c
	gcc $(CFLAGS) B1.c


B2.o: B2.c
	gcc $(CFLAGS) B2.c


B3.o: B3.c
	gcc $(CFLAGS) B3.c
	

T.o: T.c
	gcc $(CFLAGS) T.c
	

R.o: R.c
	gcc $(CFLAGS) R.c
	

IR.o: IR.c
	gcc $(CFLAGS) IR.c

Cpart1.o: Cpart1.c
	gcc $(CFLAGS) Cpart1.c
	
Cpart2.o: Cpart2.c
	gcc $(CFLAGS) Cpart2.c
	
Cpart3.o: Cpart3.c
	gcc $(CFLAGS) Cpart3.c
	
Cpart4.o: Cpart4.c
	gcc $(CFLAGS) Cpart4.c
	
Cpart5.o: Cpart5.c
	gcc $(CFLAGS) Cpart5.c
	
Cpart6.o: Cpart6.c
	gcc $(CFLAGS) Cpart6.c
	
Cpart7.o: Cpart7.c
	gcc $(CFLAGS) Cpart7.c
	
Csum.o: Csum.c
	gcc $(CFLAGS) Csum.c
	
nlo_functions.o: nlo_functions.c
	gcc $(CFLAGS) nlo_functions.c