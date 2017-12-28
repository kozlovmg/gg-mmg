CFLAGS=-c -Wall  -I/usr/local/include
LDFLAGS=-L/usr/local/lib 
LIBS=-lgsl -lgslcblas -lm 

SOURCES=main.c  nlo_functions.c  nlo.c  born.c\
B1.c  B3.c  B2.c  \
T.c  IR.c  R.c  \
Cpart1.c Cpart2.c Cpart3.c Cpart4.c Cpart5.c Cpart6.c Csum.c \
lbl-1_1.c  lbl-1_2.c  lbl-1_3.c \
lbl-2_1.c  lbl-2_2.c  lbl-2_3.c \
lbl-3_1.c  lbl-3_2.c  lbl-3_3.c \
lbl-4_1.c  lbl-4_2.c  lbl-4_3.c \
lbl-5_1.c  lbl-5_2.c  lbl-5_3.c \
lbl-6_1.c  lbl-6_2.c  lbl-6_3.c
         

OBJECTS=$(SOURCES:.c=.o)

EXECUTABLE=nlo

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	gcc $(LDFLAGS) $(OBJECTS) -o $@  $(LIBS)

.c.o:
	gcc $(CFLAGS) $< -o $@