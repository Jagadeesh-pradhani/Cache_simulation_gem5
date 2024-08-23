      SUBROUTINE IVSORT (IV,IP,N,IOUNIT)
!!
!! Module: Integer_Vector_SORT
!!
!! Purpose: Partition sorting algorithm to put the integer array IV in
!! numerical order, smallest to largest. The integer array IP is permuted
!! exactly the same as IV. If IP is initialized to IP(1:N) =(1,2,3,...,N),
!! then the pair IV,IP is returned as a reverse map of what was input.
!!
!! Reference: COLLECTED ALGORITHMS OF THE ACM - 63,64,65
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, PARAMETER :: MAXSGS = 32
!!
      INTEGER                 &
     &          IV(1:N),      & ! In/Out Array to be sorted
     &          IP(1:N),      & ! In/Out Reordered the same as IV
     &          N,            & ! Input  Length of the arrays IV and IP
     &          IHIGH(MAXSGS),& ! Local  Segment pointers
     &          ILOW(MAXSGS), & ! Local  Segment pointers
     &          IOUNIT          ! Printer output logical unit number
!!
!! Initialize.
!!
      NSEGS = 1
      IL = 1
      IH = N
!!
!! If no elements in this segment, do nothing.
!!
   10 IF (IL .GE. IH) GO TO 80
!!
!! Choose ISEP (Separation entry):
!!  Make IV(IL) <= IV((IL+IH)/2) <= IV(IH) by interchange.
!!  Set ISEP = IV((IL+IH)/2)
!!
   20 ISEPX = (IH + IL) / 2
      ISEP  = IV(ISEPX)
!!
!! IXL is lower segment index (current).
!!
      IXL = IL
!!
!! Make IV(IL) <= IV(ISEPX)
!!
      IF (IV(IL) .LE. ISEP) GO TO 30
      IV(ISEPX) = IV(IL)
      IV(IL) = ISEP
      ISEP = IV(ISEPX)
      IT = IP(ISEPX)
      IP(ISEPX) = IP(IL)
      IP(IL) = IT
!!
!! IXH is highest segment index (current).
!!
   30 IXH = IH
!!
!! Make IV(IH) >= IV(ISEPX)
!!
      IF (IV(IH) .GE. ISEP) GO TO 50
      IV(ISEPX) = IV(IH)
      IV(IH) = ISEP
      ISEP = IV(ISEPX)
      IT = IP(ISEPX)
      IP(ISEPX) = IP(IH)
      IP(IH) = IT
!!
!! Make IV(IL) <= IV(ISEPX)
!!
      IF (IV(IL) .LE. ISEP) GO TO 50
      IV(ISEPX) = IV(IL)
      IV(IL) = ISEP
      ISEP = IV(ISEPX)
      IT = IP(ISEPX)
      IP(ISEPX) = IP(IL)
      IP(IL) = IT
      GO TO 50
!!
!! Exchange low part entry which is greater than separator with high
!! part entry which is less than or equal to the separator value.
!!
   40 IT = IV(IXH)
      IV(IXH) = IV(IXL)
      IV(IXL) = IT
      IT = IP(IXH)
      IP(IXH) = IP(IXL)
      IP(IXL) = IT
!!
!! Move down upper segment as far as we can.
!!
   50 IXH = IXH - 1
      IF (IV(IXH) .GT. ISEP) GO TO 50
!!
!! Move up lower segment as far as we can.
!!
   60 IXL = IXL + 1
      IF (IV(IXL) .LT. ISEP) GO TO 60
!!
!! Nothing to do if both segments have at most one entry in common.
!!
      IF (IXL .LE. IXH) GO TO 40
!!
!! If both segments overlap, then they are separated. In this case
!! continue with shorter segment, storing the longer.
!!
      IF (IXH - IL .LE. IH - IXL) GO TO 70
!!
!! Lower segment longer, continue with upper after saving lower.
!!
      IF (NSEGS .GT. MAXSGS) THEN
        WRITE (IOUNIT,'(A/16X,A/16X,A,I3)')                                    &
     &    '0*** ERROR *** (Module IVSORT)',                                    &
     &    'Number of segments required for vector sort exceeds '//             &
     &    'pointer storage.','Increase MAXSGS, currently: ',MAXSGS
        PRINT *, ' Pointer Storage Overflow In IVSORT.'
        STOP
      ENDIF
!!
      ILOW(NSEGS) = IL
      IHIGH(NSEGS) = IXH
      IL = IXL
      NSEGS = NSEGS + 1
      GO TO 90
!!
!! Upper segment longer, continue with lower after saving upper.
!!
   70 IF (NSEGS .GT. MAXSGS) THEN
        WRITE (IOUNIT,'(A/16X,A/16X,A,I3)')                                    &
     &    '0*** ERROR *** (Module IVSORT)',                                    &
     &    'Number of segments required for vector sort exceeds '//             &
     &    'pointer storage.','Increase MAXSGS, currently: ',MAXSGS
        PRINT *, ' Pointer Storage Overflow In IVSORT.'
        STOP
      ENDIF
!!
      ILOW(NSEGS) = IXL
      IHIGH(NSEGS) = IH
      IH = IXH
      NSEGS = NSEGS + 1
      GO TO 90
!!
!! Get another segment for processing if there are any more.
!!
   80 NSEGS = NSEGS - 1
      IF (NSEGS .EQ. 0) RETURN
      IL = ILOW(NSEGS)
      IH = IHIGH(NSEGS)
!!
!! Continue to segment as long as length is greater than 11.
!!
   90 IF (IH - IL .GE. 11) GO TO 20
      IF (IL .EQ. 1) GO TO 10
      GO TO 110
!!
!! Sort elements within segment by interchange of adjacent pairs.
!!
  100 IL = IL + 1
  110 IF (IL .EQ. IH) GO TO 80
      ISEP = IV(IL + 1)
      IF (IV(IL) .LE. ISEP) GO TO 100
      IT = IP(IL + 1)
      IXL = IL
  120 IV(IXL + 1) = IV(IXL)
      IP(IXL + 1) = IP(IXL)
      IXL = IXL - 1
      IF (ISEP .LT. IV(IXL)) GO TO 120
      IV(IXL + 1) = ISEP
      IP(IXL + 1) = IT
      GO TO 100
!!
      END
!!_
      INTEGER FUNCTION INTFEXT (IEXT,JINT,ILEN,LEXT)
!!
!! Module: INTernal_From_EXTernal       Coded by: S W Key, 4-APR-1989 20:10:38
!!
!! Purpose: Returns the internal index given the external index LEXT.
!!
!! Using a binary search, the location of LEXT in the array IEXT(1:ILEN) is
!! found. The internal number for LEXT is then given by JINT(1:ILEN). The
!! search assumes IEXT(1) is the smallest value and IEXT(ILEN) is the largest.
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                 &
     &          IEXT(1:ILEN), & ! Input, Sorted vector of external ID's.
     &          JINT(1:ILEN), & ! Input, Internal pointer for external ID's.
     &          ILEN,         & ! Input, Lenght of IEXT and JINT
     &          LEXT            ! Input, Specific external ID for which
                                !        internal pointer is sought.
!!
!! Return a value of zero in the event that LEXT is not found in IEXT(1:ILEN).
!!
      INTFEXT = 0
!!
!! Initialize beginning IB and end IE of search range.
!!
      IB = 1
      IE = ILEN
!!
      IM = 0
 100  IL = IM
      IM = (IB + IE) / 2
!!
!! Test for repeated access to the same location (IM .EQ. IL) to avoid
!! an infinite loop. An infinite loop will occur if LEXT is not found
!! in the array IEXT(1:ILEN).
!!
      DO WHILE (IM .NE. IL)
        ID = LEXT - IEXT(IM)
        IF (ID .LT. 0) THEN
          IE = IM - 1
          GO TO 100
        ENDIF
        IF (ID .GT. 0) THEN
          IB = IM + 1
          GO TO 100
        ENDIF
        INTFEXT = JINT(IM)
!!
!! Force an exit from the DO WHILE when LEXT is found in IVEC(1:ILEN,1).
!!
        IL = IM
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE RVSORT (RV,IP,N,IOUNIT)
!!
!! Module: Real_Vector_SORT
!!
!! Purpose: Partition sorting algorithm to put the real array RV in
!! numerical order, smallest to largest. The integer array IP is permuted
!! exactly the same as RV. If IP is initialized to IP(1:N) =(1,2,3,...,N),
!! then the pair RV,IP is returned as a reverse map of what was input.
!!
!! Reference: COLLECTED ALGORITHMS OF THE ACM - 63,64,65
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, PARAMETER :: MAXSGS = 32
!!
      REAL(KIND(0D0))          &
     &          RV(1:N)          ! In/Out Array to be sorted
      INTEGER                  &
     &          IP(1:N),       & ! In/Out Reordered the same as RV
     &          N,             & ! Input  Length of the arrays RV and IP
     &          IHIGH(MAXSGS), & ! Local  Segment pointers
     &          ILOW(MAXSGS),  & ! Local  Segment pointers
     &          IOUNIT           ! Printer output logical unit number
!!
!! Initialize.
!!
      NSEGS = 1
      IL = 1
      IH = N
!!
!! If no elements in this segment, do nothing.
!!
   10 IF (IL .GE. IH) GO TO 80
!!
!! Choose SEP (Separation entry):
!!  Make RV(IL) <= RV((IL+IH)/2) <= RV(IH) by interchange.
!!  Set SEP = RV((IL+IH)/2)
!!
   20 ISEP = (IH + IL) / 2
      SEP  = RV(ISEP)
!!
!! IXL is lower segment index (current).
!!
      IXL = IL
!!
!! Make RV(IL) <= RV(ISEP)
!!
      IF (RV(IL) .LE. SEP) GO TO 30
      RV(ISEP) = RV(IL)
      RV(IL) = SEP
      SEP = RV(ISEP)
      IT = IP(ISEP)
      IP(ISEP) = IP(IL)
      IP(IL) = IT
!!
!! IXH is highest segment index (current).
!!
   30 IXH = IH
!!
!! Make RV(IH) >= RV(ISEP)
!!
      IF (RV(IH) .GE. SEP) GO TO 50
      RV(ISEP) = RV(IH)
      RV(IH) = SEP
      SEP = RV(ISEP)
      IT = IP(ISEP)
      IP(ISEP) = IP(IH)
      IP(IH) = IT
!!
!! Make RV(IL) <= RV(ISEP)
!!
      IF (RV(IL) .LE. SEP) GO TO 50
      RV(ISEP) = RV(IL)
      RV(IL) = SEP
      SEP = RV(ISEP)
      IT = IP(ISEP)
      IP(ISEP) = IP(IL)
      IP(IL) = IT
      GO TO 50
!!
!! Exchange low part entry which is greater than separator with high
!! part entry which is less than or equal to the separator value.
!!
   40 RT = RV(IXH)
      RV(IXH) = RV(IXL)
      RV(IXL) = RT
      IT = IP(IXH)
      IP(IXH) = IP(IXL)
      IP(IXL) = IT
!!
!! Move down upper segment as far as we can.
!!
   50 IXH = IXH - 1
      IF (RV(IXH) .GT. SEP) GO TO 50
!!
!! Move up lower segment as far as we can.
!!
   60 IXL = IXL + 1
      IF (RV(IXL) .LT. SEP) GO TO 60
!!
!! Nothing to do if both segments have at most one entry in common.
!!
      IF (IXL .LE. IXH) GO TO 40
!!
!! If both segments overlap, then they are separated. In this case
!! continue with shorter segment, storing the longer.
!!
      IF (IXH - IL .LE. IH - IXL) GO TO 70
!!
!! Lower segment longer, continue with upper after saving lower.
!!
      IF (NSEGS .GT. MAXSGS) THEN
        WRITE (IOUNIT,'(A/16X,A/16X,A,I3)')                                    &
     &    '0*** ERROR *** (Module RVSORT)',                                    &
     &    'Number of segments required for vector sort exceeds '//             &
     &    'pointer storage.','Increase MAXSGS, currently: ',MAXSGS
        PRINT *, ' Pointer Storage Overflow In RVSORT.'
        STOP
      ENDIF
!!
      ILOW(NSEGS) = IL
      IHIGH(NSEGS) = IXH
      IL = IXL
      NSEGS = NSEGS + 1
      GO TO 90
!!
!! Upper segment longer, continue with lower after saving upper.
!!
   70 IF (NSEGS .GT. MAXSGS) THEN
        WRITE (IOUNIT,'(A/16X,A/16X,A,I3)')                                    &
     &    '0*** ERROR *** (Module RVSORT)',                                    &
     &    'Number of segments required for vector sort exceeds '//             &
     &    'pointer storage.','Increase MAXSGS, currently: ',MAXSGS
        PRINT *, ' Pointer Storage Overflow In RVSORT.'
        STOP
      ENDIF
!!
      ILOW(NSEGS) = IXL
      IHIGH(NSEGS) = IH
      IH = IXH
      NSEGS = NSEGS + 1
      GO TO 90
!!
!! Get another segment for processing if there are any more.
!!
   80 NSEGS = NSEGS - 1
      IF (NSEGS .EQ. 0) RETURN
      IL = ILOW(NSEGS)
      IH = IHIGH(NSEGS)
!!
!! Continue to segment as long as length is greater than 11.
!!
   90 IF (IH - IL .GE. 11) GO TO 20
      IF (IL .EQ. 1) GO TO 10
      GO TO 110
!!
!! Sort elements within segment by interchange of adjacent pairs.
!!
  100 IL = IL + 1
  110 IF (IL .EQ. IH) GO TO 80
      SEP = RV(IL + 1)
      IF (RV(IL) .LE. SEP) GO TO 100
      IT = IP(IL + 1)
      IXL = IL
  120 RV(IXL + 1) = RV(IXL)
      IP(IXL + 1) = IP(IXL)
      IXL = IXL - 1
      IF (SEP .LT. RV(IXL)) GO TO 120
      RV(IXL + 1) = SEP
      IP(IXL + 1) = IT
      GO TO 100
!!
      END
!!_
      INTEGER FUNCTION INTFVAL (RVEC,IVEC,ILEN,VALUE)
!!
!! Module: Index_For_Value              Coded by: S W Key, 4-APR-1989 20:10:38
!!
!! Purpose: Returns the index given the real VALUE.
!!
!! Using a binary search, the location of VALUE in the array RVEC(1:ILEN) is
!! found. The index number for VALUE is then given by IVEC(*). The search
!! assumes RVEC(1) is the smallest value and RVEC(ILEN) is the largest. The
!! search is begun by interpolating linearly between RVEC(1) and RVEC(ILEN).
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      REAL(KIND(0D0))         &
     &          RVEC(1:ILEN), & ! Input, Sorted vector of real values
     &          VALUE           ! Input, Specific search value
      INTEGER                 &
     &          IVEC(1:ILEN), & ! Input, index values of unsorted RVEC(1:ILEN)
     &          ILEN            ! Input, Lenght of RVEC
!!
!! Return a value of zero in the event that VALUE is not found in RVEC(1:ILEN).
!!
      INTFVAL = 0
!!
!! Initialize beginning IB and end IE of search range.
!!
      IB = 1
      IE = ILEN
!!
      IM = 0
 100  IL = IM
      IM = (IB + IE) / 2
!!
!! Test for repeated access to the same location (IM .EQ. IL) to avoid
!! an infinite loop. An infinite loop will occur if VALUE is not found
!! in the array RVEC(1:ILEN).
!!
      DO WHILE (IM .NE. IL)
        Difference = VALUE - RVEC(IM)
        IF (Difference .LT. 0.0) THEN
          IE = IM - 1
          GO TO 100
        ENDIF
        IF (Difference .GT. 0.0) THEN
          IB = IM + 1
          GO TO 100
        ENDIF
        INTFVAL = IVEC(IM)
!!
!! Force an exit from the DO WHILE when VALUE is found in RVEC(1:ILEN).
!!
        IL = IM
      ENDDO
!!
      RETURN
      END
