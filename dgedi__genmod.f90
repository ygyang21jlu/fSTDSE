        !COMPILER-GENERATED INTERFACE MODULE: Thu Mar 30 09:34:45 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DGEDI__genmod
          INTERFACE 
            SUBROUTINE DGEDI(A,LDA,N,IPVT,DET,WORK,JOB)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: LDA
              REAL(KIND=8) :: A(LDA,N)
              INTEGER(KIND=4) :: IPVT(N)
              REAL(KIND=8) :: DET(2)
              REAL(KIND=8) :: WORK(N)
              INTEGER(KIND=4) :: JOB
            END SUBROUTINE DGEDI
          END INTERFACE 
        END MODULE DGEDI__genmod
