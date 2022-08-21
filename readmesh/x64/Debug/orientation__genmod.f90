        !COMPILER-GENERATED INTERFACE MODULE: Sat Apr 30 04:30:14 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ORIENTATION__genmod
          INTERFACE 
            SUBROUTINE ORIENTATION(P,Q,R,NODES,N,ORIENT)
              INTEGER(KIND=4), INTENT(IN) :: N
              INTEGER(KIND=4), INTENT(IN) :: P
              INTEGER(KIND=4), INTENT(IN) :: Q
              INTEGER(KIND=4), INTENT(IN) :: R
              REAL(KIND=8), INTENT(IN) :: NODES(N,3)
              INTEGER(KIND=4), INTENT(OUT) :: ORIENT
            END SUBROUTINE ORIENTATION
          END INTERFACE 
        END MODULE ORIENTATION__genmod
