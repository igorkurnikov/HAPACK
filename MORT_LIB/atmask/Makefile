#############################################################################
# Makefile for building: libatmask.a
# Generated by qmake (2.01a) (Qt 4.2.1) on: Thu Dec 21 22:03:19 2006
# Project:  atmask.pro
# Template: lib
# Command: /usr/bin/qmake -unix -o Makefile atmask.pro
#############################################################################

####### Compiler, tools and options

include ../../../config.h
ifeq ($(CXX),g++)
DEFINES       = 
CFLAGS        = -O2 -Wall -W -D_REENTRANT $(DEFINES)
CXXFLAGS      = -O2 -Wall -W -D_REENTRANT $(DEFINES)
endif
INCPATH       = -I. -I.. -I../../freelib
CFLAGS       += $(INCPATH)
CXXFLAGS     += $(INCPATH)
LEXFLAGS      = 
YACCFLAGS     = -d
INCPATH       = -I. -I.. -I../../freelib
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

SOURCES       = atom.cpp \
		condition.cpp \
		dist.cpp \
		logic.cpp \
		mask.cpp \
		paren.cpp \
		resd.cpp 
OBJECTS       = atom.o \
		condition.o \
		dist.o \
		logic.o \
		mask.o \
		paren.o \
		resd.o
DESTDIR       = 
TARGET        = libatmask.a

first: all
####### Build rules

all: Makefile $(TARGET) 

staticlib: $(TARGET)

$(TARGET):  $(OBJECTS) $(OBJCOMP) 


dist: 
	@$(CHK_DIR_EXISTS) .tmp/atmask1.0.0 || $(MKDIR) .tmp/atmask1.0.0 
	$(COPY_FILE) --parents $(SOURCES) $(DIST) .tmp/atmask1.0.0/ && $(COPY_FILE) --parents condition.hpp logic.hpp mask.hpp .tmp/atmask1.0.0/ && $(COPY_FILE) --parents atom.cpp condition.cpp dist.cpp logic.cpp mask.cpp paren.cpp resd.cpp .tmp/atmask1.0.0/ && (cd `dirname .tmp/atmask1.0.0` && $(TAR) atmask1.0.0.tar atmask1.0.0 && $(COMPRESS) atmask1.0.0.tar) && $(MOVE) `dirname .tmp/atmask1.0.0`/atmask1.0.0.tar.gz . && $(DEL_FILE) -r .tmp/atmask1.0.0


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

