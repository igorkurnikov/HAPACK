#############################################################################
# Makefile for building: libambfmt.a
# Generated by qmake (2.01a) (Qt 4.2.1) on: Thu Dec 21 22:03:11 2006
# Project:  ambfmt.pro
# Template: lib
# Command: /usr/bin/qmake -unix -o Makefile ambfmt.pro
#############################################################################

####### Compiler, tools and options

include ../../../config.h
ifeq ($(CXX),g++)
DEFINES       = 
CFLAGS        = -O2 -Wall -W -D_REENTRANT $(DEFINES)
CXXFLAGS      = -O2 -Wall -W -D_REENTRANT $(DEFINES)
endif
INCPATH       = -I.. -I../../freelib -I.
CFLAGS       += $(INCPATH)
CXXFLAGS     += $(INCPATH)
LEXFLAGS      = 
YACCFLAGS     = -d
AR            = ar cqs
RANLIB        = 
QMAKE         = /usr/bin/qmake
TAR           = tar -cf
COMPRESS      = gzip -9f
COPY          = cp -f
COPY_FILE     = $(COPY)
COPY_DIR      = $(COPY) -r
INSTALL_FILE  = install -m 644 -p
INSTALL_DIR   = $(COPY_DIR)
INSTALL_PROGRAM = install -m 755 -p
DEL_FILE      = rm -f
SYMLINK       = ln -sf
DEL_DIR       = rmdir
MOVE          = mv -f
CHK_DIR_EXISTS= test -d
MKDIR         = mkdir -p

####### Output directory

OBJECTS_DIR   = ./

####### Files

SOURCES       = ambfrc-amoeba.cpp \
		amboff.cpp \
		chain.cpp \
		energee.cpp \
		exclude.cpp \
		extent.cpp \
		les.cpp \
		parm.cpp \
		parmset.cpp \
		prep.cpp \
		prmtop-amoeba.cpp \
		prmtop.cpp
OBJECTS       = ambfrc-amoeba.o \
		amboff.o \
		chain.o \
		energee.o \
		exclude.o \
		extent.o \
		les.o \
		parm.o \
		parmset.o \
		prep.o \
		prmtop-amoeba.o \
		prmtop.o
DESTDIR       = 
TARGET        = libambfmt.a

first: all
####### Build rules

all: Makefile $(TARGET) 

staticlib: $(TARGET)

$(TARGET):  $(OBJECTS) $(OBJCOMP) 


dist: 
	@$(CHK_DIR_EXISTS) .tmp/ambfmt1.0.0 || $(MKDIR) .tmp/ambfmt1.0.0 
	$(COPY_FILE) --parents $(SOURCES) $(DIST) .tmp/ambfmt1.0.0/ && $(COPY_FILE) --parents energy.hpp exclude.hpp extent.hpp les.hpp parm.hpp parmset.hpp prmtop.hpp .tmp/ambfmt1.0.0/ && $(COPY_FILE) --parents ambfrc-amber.cpp ambfrc-amoeba.cpp amboff.cpp chain.cpp exclude.cpp extent.cpp les.cpp parm.cpp parmset.cpp prmtop-amoeba.cpp prmtop.cpp .tmp/ambfmt1.0.0/ && (cd `dirname .tmp/ambfmt1.0.0` && $(TAR) ambfmt1.0.0.tar ambfmt1.0.0 && $(COMPRESS) ambfmt1.0.0.tar) && $(MOVE) `dirname .tmp/ambfmt1.0.0`/ambfmt1.0.0.tar.gz . && $(DEL_FILE) -r .tmp/ambfmt1.0.0


yaccclean:
lexclean:
clean:compiler_clean 
	-$(DEL_FILE) $(OBJECTS)
	-$(DEL_FILE) *~ core *.core


####### Sub-libraries

distclean: clean
	-$(DEL_FILE) $(TARGET) 
	-$(DEL_FILE) Makefile


mocclean: compiler_moc_header_clean compiler_moc_source_clean

mocables: compiler_moc_header_make_all compiler_moc_source_make_all

compiler_moc_header_make_all:
compiler_moc_header_clean:
compiler_rcc_make_all:
compiler_rcc_clean:
compiler_image_collection_make_all: qmake_image_collection.cpp
compiler_image_collection_clean:
	-$(DEL_FILE) qmake_image_collection.cpp
compiler_moc_source_make_all:
compiler_moc_source_clean:
compiler_uic_make_all:
compiler_uic_clean:
compiler_clean: compiler_moc_header_clean compiler_rcc_clean compiler_image_collection_clean compiler_moc_source_clean compiler_uic_clean 

####### Compile

include depend

####### Install

install:   FORCE

uninstall:   FORCE

FORCE:

