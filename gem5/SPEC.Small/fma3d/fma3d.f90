      PROGRAM FMA_MAIN
!!
!! Copyright (c) by KEY Associates; 13-DEC-1993 21:41:48.14
!! Copyright (c) by KEY Associates;  1-APR-1997 18:15:29.00
!!
!! MAXIF = maximum number of root and included input files.
!!
      INTEGER :: MAXIF = 100

      CALL FMA_3D ( MAXIF )

      PRINT *, ' FMA-3D> Normal Termination.'
      STOP

      END
