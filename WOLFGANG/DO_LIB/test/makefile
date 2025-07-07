
HERE   = $(HOME)/lib


#  For IBM-RS6000

CC     = xlC -+
CFLAGS = -O3 -DXLC

#  For LINUX

CC       = g++
CFLAGS   = -g -DXLC -DGNU
COMPFLAG = -x c++

#  For SUN - g++

#CC     = g++
#CFLAGS = -g -DGNU

INCL   = -I$(HERE)/include -I$(HERE)/template
LIB    = -L$(HERE) -ldo

########################################################


all: test_tnumarray test_tarray test_tarray2d testvec testmat

ttest   : ttest.o 
	$(CC) $(CFLAGS) -o ttest ttest.o  $(LIB)
test_tarray   : test_tarray.o 
	$(CC) $(CFLAGS) -o test_tarray   test_tarray.o  $(LIB)
testmat       : testmat.o 
	$(CC) $(CFLAGS) -o testmat   testmat.o   $(LIB)
test_svd       : test_svd.o 
	$(CC) $(CFLAGS) -o test_svd  test_svd.o   $(LIB)
testvec       : testvec.o 
	$(CC) $(CFLAGS) -o testvec   testvec.o  $(LIB) -lblas
test_tarray2d : test_tarray2d.o 
	$(CC) $(CFLAGS) -o test_tarray2d test_tarray2d.o  $(LIB)
test_tnumarray : test_tnumarray.o 
	$(CC) $(CFLAGS) -o test_tnumarray test_tnumarray.o   $(LIB)

matprint : matprint.o
	$(CC) $(CFLAGS) -o matprint matprint.o $(LIB)


depend:
	makedepend $(INCL) *.c

.c.o:
	$(CC) $(CFLAGS) $(COMPFLAG) $(INCL) -c $<
    

clean:  
	rm *.o test_array test_numarray testvec test_array2d *~ *.bak

# DO NOT DELETE THIS LINE -- make depend depends on it.
