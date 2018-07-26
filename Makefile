
#CC=icc
#CFLAGS=-c -Wall  -I/usr/local/include  -xHost -O2 -fp-model source

CC=gcc
CFLAGS=-c -Wall  -I/usr/local/include -march=native


LDFLAGS=-L. -L/usr/local/lib 
LIBS= -lgsl -lgslcblas -lm -lquadmath 
MYLIBS=-lbubl -ltri -lbox -lpent  -llbl -lnlofunc


SOURCES=born.c   ir.c  main.c  nlo.c 









OBJECTS=$(SOURCES:.c=.o)

EXECUTABLE=nlo


	
all: libbubl libtri libbox libpent liblbl nlofunc $(SOURCES) $(EXECUTABLE) 

libbubl:
	$(MAKE) -C bubl

libtri:
	$(MAKE) -C tri


libbox:
	$(MAKE) -C box
	
libpent:
	$(MAKE) -C pent
	
liblbl:
	$(MAKE) -C lbl
	
nlofunc:
	$(MAKE) -C func
	
$(EXECUTABLE): $(OBJECTS) libbox.a liblbl.a libpent.a libnlofunc.a libtri.a libbubl.a
	$(CC) $(LDFLAGS)  $(OBJECTS) -o $@  $(MYLIBS) $(LIBS)  

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm *.o
	rm ./bubl/*.o
	rm ./tri/*.o
	rm ./box/*.o
	rm ./pent/*.o
	rm ./lbl/*.o
	rm ./func/*.o
	rm *.a