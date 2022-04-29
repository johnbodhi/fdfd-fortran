        !COMPILER-GENERATED INTERFACE MODULE: Sun Mar 13 22:33:51 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZDOT_BCG__genmod
          INTERFACE 
            FUNCTION ZDOT_BCG(N,ZX,ZY)
              INTEGER(KIND=4), INTENT(IN) :: N
              COMPLEX(KIND=8), INTENT(IN) :: ZX(N)
              COMPLEX(KIND=8), INTENT(IN) :: ZY(N)
              COMPLEX(KIND=8) :: ZDOT_BCG
            END FUNCTION ZDOT_BCG
          END INTERFACE 
        END MODULE ZDOT_BCG__genmod
