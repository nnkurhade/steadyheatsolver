        !COMPILER-GENERATED INTERFACE MODULE: Sat Apr 30 04:30:12 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LEFT_INDEX__genmod
          INTERFACE 
            SUBROUTINE LEFT_INDEX(NODES,N,MINN)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: NODES(N,3)
              INTEGER(KIND=4), INTENT(OUT) :: MINN
            END SUBROUTINE LEFT_INDEX
          END INTERFACE 
        END MODULE LEFT_INDEX__genmod
