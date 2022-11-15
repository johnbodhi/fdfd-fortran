*************************************************************************
**                 DRIVER EXAMPLE FOR THE GMRes CODE
*************************************************************************
      program validation
*
*
      integer lda, ldstrt, lwork
      parameter (lda = 20, ldstrt = 6)
      parameter (lwork = ldstrt**2 + ldstrt*(lda+5) + 5*lda + 1)
*
      integer i, j, n, m
      integer revcom, colx, coly, colz, nbscal
      integer irc(5), icntl(8), info(3)
*
      integer matvec, precondLeft, precondRight, dotProd
      parameter (matvec=1, precondLeft=2, precondRight=3, dotProd=4)
*
      integer nout
*
      complex  a(lda,lda+2), work(lwork), aux(lda+2)
      real  cntl(5), rinfo(2)
*
      complex ZERO, ONE
      parameter (ZERO = (0.0e0, 0.0e0), ONE = (1.0e0, 0.0e0))
*
* variables required by the parallel implementation
      include 'mpif.h'
      integer type, token, status(MPI_STATUS_SIZE)
      integer nproc, infompi, comm, me, nloc, iconf(2)
      integer nbcol, istart, jstart
*
*
* MPI initialization
*
      call MPI_INIT(infompi)
      comm = MPI_COMM_WORLD
      call MPI_COMM_SIZE(comm,nproc,infompi)
      call MPI_COMM_RANK(comm,me,infompi)
*
***************************************************************
** Generate the test matrix a and set the right-hand side
** in positions (n+1) to 2n of the array work.
** The right-hand side is chosen such that the exact solution
** is the vector of all ones.
***************************************************************
* The solution of the tridiagonal system is performed in parallel
* The matrix is decomposed by block of rows so that the matrix-vector
* can be easily performed.
* Each processor is in charged of a block of rows and stored the
* corresponding entries of the initial guess and the right hand-sides.
* We give below an example of the data distribution for a system of
* dimension 8 in complex arithmetic solved on 2 processors.
*
*                  A                             x  =  b
*
*
*       |  4  -2                                |   x1    b1 
*  P0   | -1+i 4  -2+i                          |   x2    b2  
*       |     -1+i   4  -2+i                    |   x3    b3      
*       |         -1+i   4  -2+i                |   x4    b4          
*  --------------------------------------------------------                 
*       |             -1+i   4  -2+i            |   x5    b5              
*       |                 -1+i   4  -2+i        |   x6    b6                
*  P1   |                     -1+i   4  -2+i    |   x7    b7               
*       |                         -1+i   4  -2+i|   x8    b8                          
*
*
* For the sake of simplicity each processor will have the same number of
* row denoted nloc. Consequently the size of the linear systems will be
* nloc times the number of processors
*
      if (me.eq.0) then
        write(*,*) '***********************************************'
        write(*,*) 'This code is an example of use of GMRES'
        write(*,*) 'Results are written in output files'
        write(*,*)  'fort.31','  ','sol_cTest_Par'
        write(*,*) '***********************************************'
        write(*,*)
        write(*,*) 'Local matrix size < ', lda
        read(*,*) nloc
        if (nloc.gt.lda) then
          write(*,*) 'You are asking for a too large matrix'
          goto 100
        endif
        write(*,*) ' Global matrix size ',nloc*nproc
*********************************
** Choose the restart parameter
*********************************
*
        write(*,*) 'Restart  <', ldstrt
        read(*,*) m
        iconf(1) = nloc
        iconf(2) = m
      endif
*
      call MPI_BCAST(iconf,2,MPI_INTEGER,0,comm,infompi)
      nloc = iconf(1)
      m    = iconf(2)
      n    = nloc*nproc
*
* Initialize the local matrix (nloc x (nloc+2) ) matrix
* only part of it might be used by the different processor
* depending on its rank
       do j = 1,nloc+2
          do i = 1,nloc
            a(i,j) = ZERO
          enddo
       enddo
*
       do i = 1,nloc
          a(i,i)   = (-1.0e0, 1.0e0)
          a(i,i+1) = (4.0e0,  0.0e0)
          a(i,i+2) = (-2.0e0, 1.0e0)
       enddo
*
* Intialise the column index of the first column 
* of the submatrix that will be involved in the mat-vec depending
* on the processor rank.
* Similarly initialize the index of the first entry of the local
* vector in the vector to be involved in the mat-vec.
*
*       jstart = 1 for all the processors but the first
*         ||       that does not have predecessor
*         \/
*                                           x1
*       |-1+i   4  -2+i             |       x2  <-- istart=2 for all the processors
*   A = |    -1+i   4  -2+i         |   x = x3               but the first
*       |        -1+i   4  -2+i     |       x4
*       |            -1+i   4  -2+i |       x5    
*                                           x6
*
       if (me.eq.0) then
        jstart = 2
        istart = 1
       else
        jstart = 1
        istart = 2
       endif
       nbCol = nloc+2
       if (me.eq.0) then
         nbCol = nbCol - 1
       endif
       if (me.eq.(nproc-1)) then
         nbCol = nbCol - 1
       endif
*
** Initialise the right hand side
      do j = 1,nloc+2
        aux(j) = ONE
      enddo
      call cgemv('N',nloc,nbCol,ONE,A(1,jstart),lda,aux,1,
     &            ZERO,work(nloc+1),1)
      do j = 1,nloc
        work(j) = ONE/2.0
      enddo
*
*
*
*******************************************************
** Initialize the control parameters to default value
*******************************************************
*
      call init_cgmres(icntl,cntl)
*
*************************
*c Tune some parameters
*************************
*
* Save the convergence history standard output
      if (me.eq.0) then
        icntl(3) = 31
      else
        icntl(1) = 0
        icntl(2) = 0
        icntl(3) = 0
      endif
* Maximum number of iterations
      icntl(7) = 100 
*
* preconditioner location
      icntl(4) = 1
* orthogonalization scheme
      icntl(5)=0
* initial guess
      icntl(6) = 0
* residual calculation strategy at restart
      icntl(8) = 1
*
*****************************************
** Reverse communication implementation
*****************************************
*
10     call drive_cgmres(n,nloc,m,lwork,work,
     &         irc,icntl,cntl,info,rinfo)
       revcom = irc(1)
       colx   = irc(2)
       coly   = irc(3)
       colz   = irc(4)
       nbscal = irc(5)
*
       if (revcom.eq.matvec) then
* perform the matrix vector product
         call ccopy(nloc,work(colx),1,aux(istart),1) 
* Send the entry of aux required to perform the parallel tridiagonal matrix-vector
* product
         if (me.ne.(nproc-1)) then
*  send the last entry of y local vector to the next processor that
*   needs this entry to perform its part of the matrix-vector product
            call MPI_SEND(aux(istart+nloc-1),1,MPI_COMPLEX,
     &                me+1,2,comm,infompi)
         endif
         if (me.ne.0) then
*  send the first entry of y local vector to the previous processor that
*   needs this entry to perform its part of the matrix-vector product
            call MPI_SEND(aux(istart),1,MPI_COMPLEX,me-1,
     &                3,comm,infompi)
         endif
* Receive the entry of aux required to perform the parallel tridiagonal matrix-vector
* product
         if (me.ne.(nproc-1)) then
*  receive the last entry of the vector to be involved in the mat-vec.
* this entry is computed by the next processor
            call MPI_RECV(aux(istart+nloc),1,MPI_COMPLEX,
     &                me+1,3,comm,status,infompi)
         endif
         if (me.ne.0) then
*  receive the first entry of the vector to be involved in the mat-vec.
* this entry is computed by the previous processor
            call MPI_RECV(aux(1),1,MPI_COMPLEX,me-1,
     &                2,comm,status,infompi)
         endif
*  Compute the local matrix-vector product
         call cgemv('N',nloc,nbCol,ONE,A(1,jstart),lda,aux,1,
     &            ZERO,work(colz),1)
*
         goto 10
*
       else if (revcom.eq.precondLeft) then
* perform the left preconditioning
*         work(colz) <-- M^{-1} * work(colx)
         call ccopy(nloc,work(colx),1,work(colz),1)
         goto 10
*
       else if (revcom.eq.precondRight) then
* perform the right preconditioning
         call ccopy(nloc,work(colx),1,work(colz),1)
         goto 10
*
       else if (revcom.eq.dotProd) then
*      perform the scalar product
*      work(colz) <-- work(colx) work(coly)
*
         call cgemv('C',nloc,nbscal,ONE,work(colx),nloc,
     &               work(coly),1,ZERO,aux,1)
         call MPI_ALLREDUCE(aux,work(colz),nbscal,
     &          MPI_COMPLEX, MPI_SUM,comm,infompi)
         goto 10
       endif
*
*******************************
* dump the solution on a file
*******************************
*
      if (me.eq.0) then 
        nout = 11
        open(nout,FILE='sol_cTest_Par',STATUS='unknown')
        if (icntl(5).eq.0) then
          write(nout,*) 'Orthogonalisation : MGS'
        elseif (icntl(5).eq.1) then
          write(nout,*) 'Orthogonalisation : IMGS'
        elseif (icntl(5).eq.2) then
          write(nout,*) 'Orthogonalisation : CGS'
        elseif (icntl(5).eq.3) then
          write(nout,*) 'Orthogonalisation : ICGS'
        endif
        write(nout,*) 'Restart : ', m
        write(nout,*) 'info(1) = ',info(1),'  info(2) = ',info(2)
        write(nout,*) 'rinfo(1) = ',rinfo(1),'  rinfo(2) = ',rinfo(2)
        write(nout,*) 'Optimal local workspace = ', info(3)
        write(nout,*) ' ********************************************* '
        write(nout,*)
        write(nout,*) (work(i),i=1,min(10,nloc))
        close(nout)
      endif
*
*
100    continue
*
      call MPI_FINALIZE(infompi)
*
      stop
      end
