.SUFFIXES: .SUFFIXES .tex

TEX = tarray.tex tnumarray.tex numer.tex stringc.tex index.tex intro.tex \
      inverse1.tex inverse2.tex determinant1.tex determinant2.tex \
      diagonalize1.tex diagonalize2.tex  blas.tex

libdoc.dvi: $(TEX) libdoc.tex
	latex libdoc.tex

blas.tex  : ../include/blas.h
	../util/totex ../include/blas.h
	mv ../include/blas.tex .
tarray.tex: ../template/tarray.h
	../util/totex ../template/tarray.h
	mv ../template/tarray.tex .
inverse2.tex: ../inverse.c
	../util/totex ../inverse.c
	mv ../inverse.tex inverse2.tex
inverse1.tex: ../include/inverse.h
	../util/totex ../include/inverse.h
	mv ../include/inverse.tex inverse1.tex
determinant2.tex: ../determinant.c
	../util/totex ../determinant.c
	mv ../determinant.tex determinant2.tex
determinant1.tex: ../include/determinant.h
	../util/totex ../include/determinant.h
	mv ../include/determinant.tex determinant1.tex
diagonalize2.tex: ../diagonalize.c
	../util/totex ../diagonalize.c
	mv ../diagonalize.tex diagonalize2.tex
diagonalize1.tex: ../include/diagonalize.h
	../util/totex ../include/diagonalize.h
	mv ../include/diagonalize.tex diagonalize1.tex
tnumarray.tex: ../template/tnumarray.h
	../util/totex ../template/tnumarray.h
	mv ../template/tnumarray.tex .
numer.tex: ../include/numer.h
	../util/totex ../include/numer.h
	mv ../include/numer.tex .
stringc.tex: ../stringc.c ../include/stringc.h
	../util/totex ../include/stringc.h
	mv ../include/stringc.tex .
	../util/totex ../stringc.c
	cat ../stringc.tex >> stringc.tex
	rm ../stringc.tex
index.tex: libdoc.idx
	makeindex libdoc.idx
	mv libdoc.ind index.tex
