.SUFFIXES: .SUFFIXES .tex .dvi .c,v .h,v
%.c : RCS/%.c,v
	-co $*.c

%.h : RCS/%.h,v
	-co $*.h
#
# update all RCS files if neccesary
#
all: *.c *.h
	echo " Updating Files in Generic "

