#  For IBM-RS6000
#CC     = xlC -+
#CFLAGS = -g 

#  For SUN - g++
CC     = g++ 
CFLAGS = -g -DGNU

INCL   = -I../include -I../template
LIBDIR  = -L../

.c.o:
	$(CC) $(CFLAGS) $(INCL) -c $<

all: totex makedep

totex: totex.o
	$(CC) $(CCFLAGS) $(INCL) -o totex   totex.o   $(LIBDIR) -ldo -lm
triang: triang.o
	$(CC) $(CCFLAGS) $(INCL) -o triang  triang.o   $(LIBDIR) -ldo -lm
makedep: makedep.o 
	$(CC) $(CCFLAGS) $(INCL) -o makedep makedep.o $(LIBDIR) -ldo -lm

triang.o:	triang.c
totex.o:        totex.c 
makedep.o: 	makedep.c
