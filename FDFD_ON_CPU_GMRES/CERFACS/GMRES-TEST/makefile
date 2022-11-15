# Mac OS-X
FORTRAN = mpif77
CC = cc
OPTS    = -u
BLAS   = /sw/lib/libf77blas.a /sw/lib/libatlas.a /sw/lib/libf2c.a /usr/lib/libcc_dynamic.a blas_rot.o
LAPACK = -llapack
OTHER  = -lmpi
# user's object files
#
sGMRES = sPackgmres.o
dGMRES = dPackgmres.o
cGMRES = cPackgmres.o
zGMRES = zPackgmres.o
#
.SUFFIXES: .f .o
.SUFFIXES: .F .f
#
sarith: sDriver sParDriver
darith: dDriver dParDriver
carith: cDriver cParDriver
zarith: zDriver zParDriver
all: sDriver dDriver cDriver zDriver sParDriver dParDriver cParDriver zParDriver
#
#
sDriver: sDriver.o $(sGMRES) 
	$(FORTRAN) $(LOPTS) -o $@ $@.o $(sGMRES) $(BLAS) $(OTHER)
#
dDriver: dDriver.o $(dGMRES) 
	$(FORTRAN) $(LOPTS) -o $@ $@.o $(dGMRES) $(BLAS) $(OTHER)
#
cDriver: cDriver.o $(cGMRES) 
	$(FORTRAN) $(LOPTS) -o $@ $@.o $(cGMRES) $(BLAS) $(OTHER)
#
zDriver: zDriver.o $(zGMRES) 
	$(FORTRAN) $(LOPTS) -o $@ $@.o $(zGMRES) $(BLAS) $(OTHER)
#
sParDriver: sParDriver.o $(sGMRES) 
	$(FORTRAN) $(LOPTS) -o $@ $@.o $(sGMRES) $(BLAS) $(OTHER)
#
dParDriver: dParDriver.o $(dGMRES) 
	$(FORTRAN) $(LOPTS) -o $@ $@.o $(dGMRES) $(BLAS) $(OTHER)
#
cParDriver: cParDriver.o $(cGMRES) 
	$(FORTRAN) $(LOPTS) -o $@ $@.o $(cGMRES) $(BLAS) $(OTHER)
#
zParDriver: zParDriver.o $(zGMRES) 
	$(FORTRAN) $(LOPTS) -o $@ $@.o $(zGMRES) $(BLAS) $(OTHER)
#
#
.f.o :
	$(FORTRAN) $(FOPTS) -c $*.f
clean :
	/bin/rm -f *.o sDriver dDriver cDriver zDriver sParDriver dParDriver cParDriver zParDriver
