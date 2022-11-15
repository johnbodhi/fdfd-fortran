        !COMPILER-GENERATED INTERFACE MODULE: Fri Nov 12 22:15:05 2021
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZBCG2__genmod
          INTERFACE 
            SUBROUTINE ZBCG2(PRINT_RESID,L,N,X,NONZERO_X,RHS,MATVEC,    &
     &PRECOND,TOLER,MXMATVEC,WORK,INFO)
              INTEGER(KIND=4), INTENT(IN) :: N
              INTEGER(KIND=4), INTENT(IN) :: L
              LOGICAL(KIND=4), INTENT(IN) :: PRINT_RESID
              COMPLEX(KIND=8), INTENT(INOUT) :: X(N)
              LOGICAL(KIND=4), INTENT(IN) :: NONZERO_X
              COMPLEX(KIND=8), INTENT(IN) :: RHS(N)
              EXTERNAL MATVEC
              EXTERNAL PRECOND
              REAL(KIND=8), INTENT(INOUT) :: TOLER
              INTEGER(KIND=4), INTENT(INOUT) :: MXMATVEC
              COMPLEX(KIND=8), INTENT(OUT) :: WORK(N,3+2*(L+1))
              INTEGER(KIND=4), INTENT(OUT) :: INFO
            END SUBROUTINE ZBCG2
          END INTERFACE 
        END MODULE ZBCG2__genmod
