# this makefile is for LINUX machines only 
#  ILUC removed.. 
OBJ	=  fgmr.o iluk.o ilut.o arms2.o vbiluk.o vbilut.o auxill.o

LIB  = ./LIB/PQ.o ./LIB/piluNEW.o ./LIB/indsetC.o \
       ./LIB/sets.o ./LIB/tools.o ./LIB/systimer.o ./LIB/misc.o  \
       ./LIB/MatOps.o ./LIB/ilutpC.o ./LIB/setblks.o ./LIB/svdInvC.o
    
AR = ar -rcv

#
FC      =  f77
FCFLAGS =  -c -g -Wall
CC      =  gcc
CCFLAGS =  -c -g -DLINUX -Wall -O3
LD      =  f77 
LDFLAGS = 
#
# clear list of default suffixes, and declare default suffixes
.SUFFIXES:
.SUFFIXES: .f .c .o
# default rule to make .o files from .f files
.f.o  : ;       $(FC) $(FCFLAGS) $*.f -o $*.o
.c.o  : ;       $(CC) $(CCFLAGS) $*.c -o $*.o
#

lib libitsol.a: $(OBJ) $(LIB)
	$(AR) libitsol.a  $(OBJ) $(LIB)
##       ranlib libitsol.a
	     
#
clean :
	rm -f ${OBJ} *.o *.ex *~ core *.cache ${LIB} LIB/*~ LIB/*.cache OUT/*

cleanall :
	rm -f ${OBJ} *.o *.ex *.a *.cache *~ core ${LIB} \
        LIB/*~ LIB/*.cache OUT/* TESTS_COO/*.o TESTS_COO/*.ex \
        TESTS_HB/*.o TESTS_HB/*.ex 
