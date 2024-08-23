      SUBROUTINE GETIRV (NVALS,LCRD,LPRT,EOF,IERR,TEXT,IRINDEX,ILCOUNT)
!!
!! Module: GET_Input_Record_Values   Coded by: S W Key, 25-MAR-1989 15:26:24
!!     Copyright (c) by KEY Associates;  9-DEC-1997 18:57:09.00
!!     Copyright (c) by KEY Associates;  1-JAN-1999 11:42:38:00
!!
!! Purpose: Read next input record, identify and encode all values present.
!!
!! Free-Field Record Reader. This module expects to find and process
!! RCRD(1:80) = <value><delimiter>... It will read the next input record
!! and any continuation lines, ignoring any preceding or interposed
!! comment or blank (empty) lines. The grammar rules are as follows:
!!
!! 1. At the start ALL values are initialized to null.
!!
!! 2. Each input line is 80 characters or less in length.
!!
!! 3. Values are separated with delimiters. A delimiter is any of the
!!    following:
!!      a. a single tab surrounded by any number of spaces,
!!      b. a single comma surrounded by any number of spaces,
!!      c. a single equal surrounded by any number of spaces,
!!      d. a single asterisk surrounded by any number of spaces,
!!      e. any number of consecutive spaces not bounded by
!!         commas, tabs, equals, or asterisks,
!!      f. End-Of-Record character ($).
!!
!! 4. Blank lines and comment lines are ignored. Leading blanks before
!!    the first value occurring in a line are ignored.
!!
!! 5. All values must be separated by a delimiter. Repeated delimiters
!!    with the exception of blanks produce a null value for the parameter
!!    in that position.
!!
!! 6. An asterisk "*" is SIMULTANEOUSLY a delimiter and a continuation
!!    character meaning additional values will be found on the following
!!    line.
!!
!! 7. A dollar sign "$" denotes a comment. All characters to the right of
!!    the dollar sign are ignored.
!!      a. If a dollar sign is the first non-space character, the line is a
!!         comment line and processing proceeds to the next line for data.
!!      b. If a dollar sign is not the first non-space character, the data
!!         expected from the input record is assumed to be complete.
!!
!! 8. Any sequence of character constants, signed integer constants, and
!!    floating point constants is acceptable. Only the first 32 characters
!!    in a character constant are reported. If a floating point constant
!!    is supplied for an integer variable, NINT is used to perform a mode
!!    conversion.
!!
      USE value_
!!
!! Arguments.
      INTEGER,       INTENT(INOUT) :: NVALS   ! Num values expected/num returned
      INTEGER,       INTENT(IN)    :: LCRD    ! Input file logical unit number
      INTEGER,       INTENT(IN)    :: LPRT    ! Print file logical unit number
      LOGICAL,       INTENT(OUT)   :: EOF     ! Flag for End-Of-File
      INTEGER,       INTENT(INOUT) :: IERR    ! Error counter
      CHARACTER(80), INTENT(OUT)   :: TEXT    ! First line of input record
      INTEGER,       INTENT(OUT)   :: IRINDEX ! Input line number for keyword
      INTEGER,       INTENT(INOUT) :: ILCOUNT ! Input line counter
!!
!! Local variables.
      CHARACTER( 1), SAVE :: NUL = "N"  ! Null value indicator
      CHARACTER( 1), SAVE :: UNK = "U"  ! Unknown value indicator
      CHARACTER( 1), SAVE :: DOL = "$"  ! Comment/blank line value
      CHARACTER( 1), SAVE :: SPA = " "  ! Space value
      CHARACTER(81), SAVE :: DTS        ! Dots for error message
      CHARACTER(81), SAVE :: PTS        ! Up-arrows for error message
      CHARACTER(81)       :: UMK        ! Error message construction buffer
      CHARACTER(81)       :: RCRD       ! Input record buffer
      LOGICAL             :: EOR        ! Flag for end of record
      LOGICAL             :: CONT       ! Flag for continuation line
      LOGICAL             :: SPFND      ! Flag for space
      LOGICAL             :: CLINE      ! Flag for comment line
      LOGICAL             :: CHECK      ! Flag for unrecognized input entries
!!
!! External functions.
      CHARACTER(32)       :: C_VALUE    ! Returns character value
      EXTERNAL               C_VALUE

      LOGICAL, SAVE :: FIRST = .TRUE.
!!
      IF (FIRST) THEN
        DO N = 1,81
          DTS(N:N) = '.'
          PTS(N:N) = '^'
        ENDDO
        FIRST = .FALSE.
      ENDIF
!!
!! Initialize all values to "null."
!!
      DO N = 1,NVALS
        VALUE(N)%VTYP = NUL
      ENDDO
!!
!! Read input lines until a non-comment or a non-blank line is found.
!!
      CLINE = .TRUE.
      DO WHILE (CLINE)
        EOF = .TRUE.
        ILCOUNT = ILCOUNT + 1
        READ (LCRD,'(A)',ERR=100) RCRD
        EOF = .FALSE.
 100    IF (EOF) RETURN
        RCRD(81:81) = '$'
!!
!! Find first non-space character (ignore leading spaces).
!!
        I = 0
        SPFND = .TRUE.
        DO WHILE (SPFND)
          I = I + 1
          SPFND = (RCRD(I:I) .EQ. SPA)
        ENDDO
        CLINE = (RCRD(I:I) .EQ. DOL)
      ENDDO
!!
!! Save first line of the input record (no continuation lines are saved)
!! in case a character string (title or descriptive text) is expected.
!!
      TEXT = RCRD(1:80)
!!
!! Report line number in input file where the input record starts.
!!
      IRINDEX = ILCOUNT
!!
!! Extract the requested number of values from RCRD, VALUE(1:NVALS).
!!
      I = 1
      N = 1
      ILEN = 81
      EOR  = .FALSE.
      CONT = .FALSE.
      DO WHILE (.NOT.EOR .AND. N .LE. NVALS)
         CALL RCRDRD (RCRD(I:81),ILEN,N,INXT,CONT,EOR)
!!
!! If a continuation delimiter was found, read the next line.
!! If the next line is a comment or blank line, read another
!! line.
!!
         IF (CONT) THEN
            CLINE = .TRUE.
            DO WHILE (CLINE)
              EOF = .TRUE.
              ILCOUNT = ILCOUNT + 1
              READ (LCRD,'(A)',ERR=110) RCRD
              EOF = .FALSE.
 110          IF (EOF) RETURN
              RCRD(81:81) = '$'
              I = 0
              SPFND = .TRUE.
              DO WHILE (SPFND)
                I = I + 1
                SPFND = (RCRD(I:I) .EQ. SPA)
              ENDDO
              CLINE = (RCRD(I:I) .EQ. DOL)
            ENDDO
            I = 1
            ILEN = 81
            CONT = .FALSE.
         ELSE
            I = I + INXT
            ILEN = ILEN - INXT
         ENDIF
!!
         N = N + 1
      ENDDO
!!
!! Report actual number of entries found in input record.
!!
      NVALS = N - 1
!!
!! If any values of an unknown type were found, write a message and
!! increment the error counter for each occurrence. Skip input records
!! known to be text lines.
!!
      CHECK =                                &
     &  (                                    &
     &  INDEX(C_VALUE(1),'TITLE')   .EQ. 0   &
     &  .AND.                                &
     &  INDEX(C_VALUE(1),'QAREC')   .EQ. 0   &
     &  .AND.                                &
     &  INDEX(C_VALUE(1),'INCLUDE') .EQ. 0   &
     &  )
      IF (CHECK) THEN
        DO N = 1,NVALS
          IF (VALUE(N)%VTYP .EQ. UNK) THEN
            J = VALUE(N)%LOC
            K = INDEX(VALUE(N)%CVAL(1:32),SPA)
            IF (K .EQ. 0) THEN
              K = J + 31
            ELSE
              K = J + K - 1
            ENDIF
            UMK = DTS
            UMK(J:K) = PTS(J:K)
            WRITE (LPRT,'(/A/(18X,A/))')                       &
     &        '*** WARNING *** Message From Module: GETIRV',   &
     &        'Unknown Value Type Detected In Input Record:',  &
     &        '['//RCRD(1:81)//']',                            &
     &        '['//UMK//']'
            IERR = IERR + 1
          ENDIF
        ENDDO
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE RCRDRD (RCRD,ILEN,NVAL,INXT,CONT,EOR)
!!
!! Module: ReCoRD_ReaD          Coded by: S W Key, 11-MAR-1989 15:09:41
!!     Copyright (c) by KEY Associates;  1-JAN-1999 11:42:38:00
!!
!! Called by: GETIRV
!!
!! Purpose: Read the next entry in the record, identify it, and encode it.
!!
!! This module expects to find and process RCRD(1:) = <value><delimiter>...
!! where <value>
!!      = a character constant
!!      = a signed integer constant, or
!!      = a real constant,
!! and <delimiter>
!!      = a single comma w/ any number of spaces,
!!      = a single tab w/ any number of spaces,
!!      = a single equal w/ any number of spaces,
!!      = a single asterisk w/ any number of spaces,
!!      = any number of consecutive spaces, or
!!      = an End-Of-Record character ($).
!!
!! If <value> is missing, it is identified as "null" (VALUE(NVAL)%VTYP = 'N').
!! This module will leave INXT pointing to the first position of the next
!! <value><delimiter> pair in RCRD.
!!
      USE value_
!!
!! Arguments.
      CHARACTER(*), INTENT(IN)  :: RCRD  ! Character string
      INTEGER,      INTENT(IN)  :: ILEN  ! Length of character string RCRD
      INTEGER,      INTENT(IN)  :: NVAL  ! Index for value being processed.
      INTEGER,      INTENT(OUT) :: INXT  ! Off-set to next entry
      LOGICAL,      INTENT(OUT) :: CONT  ! Flag for continuation
      LOGICAL,      INTENT(OUT) :: EOR   ! Flag for End-Of-Record
!!
!! Local variables.
      CHARACTER(1), SAVE :: NUL = "N"  ! Value flag, null
      CHARACTER(1), SAVE :: CHR = "C"  ! Value flag, character
      CHARACTER(1), SAVE :: INT = "I"  ! Value flag, integer
      CHARACTER(1), SAVE :: REA = "R"  ! Value flag, real
      CHARACTER(1), SAVE :: DBL = "D"  ! Value flag, double precision
      CHARACTER(1), SAVE :: UNK = "U"  ! Value flag, unknown
      CHARACTER(1), SAVE :: TAB = " "  ! Tab
      CHARACTER(1), SAVE :: COM = ","  ! Comma
      CHARACTER(1), SAVE :: EQU = "="  ! Equal
      CHARACTER(1), SAVE :: SPA = " "  ! Space
      CHARACTER(1), SAVE :: AST = "*"  ! Asterisk
      CHARACTER(1), SAVE :: DOL = "$"  ! Dollar sign
      LOGICAL            :: FOUND      ! Flag for next delimiter
      LOGICAL            :: ERFND      ! Flag for End-Of-Record
      LOGICAL            :: CNFND      ! Flag for continue delimiter *
      LOGICAL            :: DLFND      ! Flag for space,comma,tab,*,$,=
      LOGICAL            :: TBFND      ! Flag for TaB_FouND
      LOGICAL            :: CMFND      ! Flag for CoMma_FouND
      LOGICAL            :: EQFND      ! Flag for EQual_FouND
      LOGICAL            :: SPFND      ! Flag for SPace_FouND
      LOGICAL            :: NTASP      ! Flag for NoT_A_SPace
!!
!! External functions.
      LOGICAL            :: INFND      ! Ftn for identifying an integer value
      LOGICAL            :: REFND      ! Ftn for identifying a real value
      LOGICAL            :: DPFND      ! Ftn for identifying a double precision
      LOGICAL            :: CHFND      ! Ftn for identifying a character string
      EXTERNAL INFND, REFND, DPFND, CHFND
!!
      TAB = ACHAR(9)
!!
!! Step over leading spaces before <value>.
!!
      I = 0
      FOUND = .FALSE.
      DO WHILE (.NOT.FOUND .AND. I .LT. ILEN)
        I = I + 1
        FOUND = (RCRD(I:I) .NE. SPA)
      ENDDO
!!
!! Look for the next delimiter.
!!
      J = I
      I = I - 1
      FOUND = .FALSE.
      DO WHILE (.NOT.FOUND .AND. I .LT. ILEN)
        K = I
        I = I + 1
        FOUND = DLFND (RCRD(I:I))
      ENDDO
!!
!! Examine RCRD(J:K) for variable type and encode numeric values.
!!
      LJK = K - J + 1
      IF (LJK .EQ. 0) THEN
        VALUE(NVAL)%VTYP = NUL
      ELSE
        VALUE(NVAL)%VTYP = UNK
        VALUE(NVAL)%LOC  = J + 81 - ILEN
!!
!! In any event save input record entry as character string.
!!
        VALUE(NVAL)%CVAL = RCRD(J:J+MIN(31,LJK-1))
!!
!! Test for valid/recognized entry types.
!!
        IF (INFND(RCRD(J:K),LJK)) THEN
          READ (RCRD(J:K),*,ERR=100) VALUE(NVAL)%IVAL
          VALUE(NVAL)%VTYP = INT
 100      CONTINUE
        ELSE IF (REFND(RCRD(J:K),LJK)) THEN
          READ (RCRD(J:K),*,ERR=200) VALUE(NVAL)%RVAL
          VALUE(NVAL)%VTYP = REA
 200      CONTINUE
        ELSE IF (DPFND(RCRD(J:K),LJK)) THEN
          READ (RCRD(J:K),*,ERR=300) VALUE(NVAL)%DVAL
          VALUE(NVAL)%VTYP = DBL
 300      CONTINUE
        ELSE IF (CHFND(RCRD(J:K),LJK)) THEN
          VALUE(NVAL)%VTYP = CHR
        ENDIF
      ENDIF
!!
!! Locate the next delimiter. Look for one of the following:
!!      1. An End-Of-Record character ($)
!!      2. One or more consecutive spaces
!!      3. A comma with or without leading spaces
!!      4. A tab with or without leading spaces
!!      5. An equal with or without leading spaces
!!      6. An asterisk with or without leading spaces
!!
      I = K
      FOUND = .FALSE.
      DO WHILE (.NOT.FOUND .AND. I .LT. ILEN)
        I = I + 1
        ERFND = (RCRD(I:I) .EQ. DOL)
        CNFND = (RCRD(I:I) .EQ. AST)
        CMFND = (RCRD(I:I) .EQ. COM)
        TBFND = (RCRD(I:I) .EQ. TAB)
        EQFND = (RCRD(I:I) .EQ. EQU)
        NTASP = (RCRD(I:I) .NE. SPA)
        FOUND = ERFND .OR. CNFND .OR. CMFND .OR. TBFND .OR. NTASP .OR. EQFND
      ENDDO
!!
      SPFND = NTASP .AND. .NOT.(CMFND .OR. TBFND .OR. EQFND)
      IF (SPFND) THEN
        INXT = I - 1
      ELSE
        INXT = I
      ENDIF
      CONT = CNFND
      EOR  = ERFND
!!
      RETURN
      END
!!_
      LOGICAL FUNCTION DLFND ( CRX )
!!
!! Copyright (c) by KEY Associates;  1-JAN-1999 11:42:38:00
!!
!! Called by: RCRDRD
!!
!! Purpose: Identify character CRX as a delimiter.
!!
!! Description: If CRX equals a delimiter (Tab,Space,$,*,Comma,=), return
!! a value of .TRUE. otherwise return .FALSE.
!!
!! Arguments.
      CHARACTER(1), INTENT(IN) :: CRX       ! Character to be tested
!!
!! Local variables.
      CHARACTER(6), SAVE :: DEL = "  $*,="  ! Delimiter set
!!
      DEL(1:1) = ACHAR(9)
!!
      DLFND = (INDEX(DEL,CRX) .NE. 0)
!!
      RETURN
      END
!!_
      LOGICAL FUNCTION INFND ( CRX,ILEN )
!!
!! Copyright (c) by KEY Associates;  1-JAN-1999 11:42:38:00
!!
!! Called by: RCRDRD
!!
!! Purpose: Identify character string CRX as an integer.
!!
!! Description: If CRX equals an integer, return a value of .TRUE.
!! otherwise return .FALSE.
!!
!! Arguments.
      CHARACTER(*), INTENT(IN) :: CRX   ! String to be tested
      INTEGER,      INTENT(IN) :: ILEN  ! String length
!!
!! Local variables.
      CHARACTER(10), SAVE :: NTEGER = "1234567890"
      CHARACTER( 2), SAVE :: OTHER  = "+-"
!!
!! Do the simplest thing first: look for all integers in CRX
!!
      INFND = (VERIFY(CRX(1:ILEN),NTEGER) .EQ. 0)
!!
!! If something other than an integer was found, examine in more 
!! detail. We know CRX contains something other than integers. 
!! The only other acceptable possiblity is a signed integer.
!!
      IF (.NOT.INFND) THEN
!!
!! If ILEN is greater than one, check for a leading sign.
!!
        IF (ILEN .GT. 1) THEN
          INFND = (INDEX(OTHER,CRX(1:1)) .NE. 0) .AND.  &
     &            (VERIFY(CRX(2:ILEN),NTEGER) .EQ. 0)
        ENDIF
      ENDIF
!!
      RETURN
      END
!!_
      LOGICAL FUNCTION REFND ( CRX,ILEN )
!!
!! Copyright (c) by KEY Associates;  1-JAN-1999 11:42:38:00
!!
!! Called by: RCRDRD
!!
!! Purpose: Identify character string CRX as a real constant.
!!
!! Description: If CRX equals a real, return a value of .TRUE.
!! otherwise return .FALSE.
!!
!! Assumptions: CRX(1:ILEN) has already been tested as an integer.
!!
!! Arguments.
      CHARACTER(*), INTENT(IN) :: CRX   ! String to be tested
      INTEGER,      INTENT(IN) :: ILEN  ! String length
!!
!! Local variables.
      CHARACTER(10), SAVE :: NTEGER = "1234567890"
      LOGICAL             :: VALIDC  ! Valid character set
      LOGICAL             :: SIMPLE  !
      LOGICAL             :: ONEINT  ! One integer character
      LOGICAL             :: ONEDPT  ! One decimal point character
      LOGICAL             :: SIGNED  ! A + or - sign is present
      LOGICAL             :: ONEEXP  ! One exponent character
      LOGICAL             :: EXPINT  ! Exponent is an integer
      LOGICAL             :: INTBFE  ! Integer before the exponent
      LOGICAL, PARAMETER  :: BACK = .TRUE.
!!
!! External functions.
      LOGICAL             :: INFND
      EXTERNAL INFND
!!
!! Assure we have only valid characters for a real constant.
!!
      IF (ILEN .EQ. 1) THEN
        VALIDC = (VERIFY(CRX(1:ILEN),NTEGER) .EQ. 0)
      ELSE
        VALIDC = (VERIFY(CRX(1:ILEN),NTEGER//".Ee+-") .EQ. 0)
      ENDIF
      IF (.NOT.VALIDC) THEN
        REFND = .FALSE.
        RETURN
      ENDIF
!!
!! Run a simple validation test -- one decimal point.
!!
      ONEDPT = (INDEX(CRX(1:ILEN),".") .EQ. INDEX(CRX(1:ILEN),".",BACK))
      IF (.NOT.ONEDPT) THEN
        REFND = .FALSE.
        RETURN
      ENDIF
!!
!! Assume the simplest case first: integers with a decimal point.
!!
      SIMPLE = (VERIFY(CRX(1:ILEN),NTEGER//".") .EQ. 0)
      IF (SIMPLE) THEN
        REFND = .TRUE.
        RETURN
      ENDIF
!!
!! Test for next simplest case, a leading +/- sign.
!!
      SIGNED = (INDEX("+-",CRX(1:1)) .NE. 0)
      SIMPLE = (VERIFY(CRX(2:ILEN),NTEGER//".") .EQ. 0)
      ONEINT = (SCAN(CRX(2:ILEN),NTEGER) .NE. 0)
      IF (SIGNED .AND. SIMPLE .AND. ONEINT) THEN
        REFND = .TRUE.
        RETURN
      ENDIF
!!
!! Test for a single exponent character, for example, 1E6
!!
      CALL CUPPER(CRX(1:ILEN))
      ONEEXP = (INDEX(CRX(1:ILEN),"E") .EQ. INDEX(CRX(1:ILEN),"E",BACK))
      IF (.NOT.ONEEXP) THEN
        REFND = .FALSE.
        RETURN
      ENDIF
!!
!! Divide CRX(1:ILEN) at the exponent.
!!
      KEXP = INDEX(CRX(1:ILEN),"E")
      IF (KEXP.EQ.ILEN .OR. KEXP.EQ.1) THEN
        REFND = .FALSE.
        RETURN
      ENDIF
!!
!! Test for exponent being at most a signed integer.
!!
      EXPINT = INFND(CRX(KEXP+1:ILEN),ILEN-KEXP)
      IF (.NOT.EXPINT) THEN
        REFND = .FALSE.
        RETURN
      ENDIF
!!
!! Assume the simplest case first: a non-decimal number before E,
!! for example, +1E-6 or -1E+4 or 2E-7
!!
      SIMPLE = INFND(CRX(1:KEXP-1),KEXP-1)
      IF (SIMPLE) THEN
        REFND = .TRUE.
        RETURN
      ENDIF
!!
!! Last case; there must be a decimal point present before E, for
!! example, -.1E6. First insure there is at least an integer before E.
!!
      INTBFE = (SCAN(CRX(1:KEXP-1),NTEGER) .NE. 0)
      IF (.NOT.INTBFE) THEN
        REFND = .FALSE.
        RETURN
      ENDIF
!!
!! The only thing left is to insure is that any sign present is in
!! the first location.
!!
      IMIN = INDEX(CRX(1:KEXP-1),"-")
      IPOS = INDEX(CRX(1:KEXP-1),"+")
      IF (IMIN.EQ.0 .AND. IPOS.EQ.0) THEN
        REFND = .TRUE.
      ELSE IF (IMIN.EQ.1 .AND. IPOS.EQ.0) THEN
        REFND = .TRUE.
      ELSE IF (IMIN.EQ.0 .AND. IPOS.EQ.1) THEN
        REFND = .TRUE.
      ELSE
        REFND = .FALSE.
      ENDIF
!!
      RETURN
      END
!!_
      LOGICAL FUNCTION DPFND ( CRX,ILEN )
!!
!! Copyright (c) by KEY Associates;  1-JAN-1999 11:42:38:00
!!
!! Module: Double_Precision_FouND
!!
!! Called by: RCRDRD
!!
!! Purpose: Identify character string CRX as a double precision.
!!
!! Description: If CRX equals a double precision, return a value 
!! of .TRUE. otherwise return .FALSE.
!!
!! Assumptions: CRX(1:ILEN) has already been tested as an integer.
!!
!! Arguments.
      CHARACTER(*), INTENT(IN) :: CRX   ! String to be tested
      INTEGER,      INTENT(IN) :: ILEN  ! String length
!!
!! Local variables.
      CHARACTER(10), SAVE :: NTEGER = "1234567890"
      LOGICAL             :: VALIDC  ! Valid character set
      LOGICAL             :: ONEEXP  ! One exponent character
      LOGICAL             :: EXPINT  ! Exponent is integer
      LOGICAL             :: SIMPLE
      LOGICAL             :: ONEDPT  ! One decimal point
      LOGICAL             :: INTBFD  ! Integer before D
      LOGICAL, PARAMETER  :: BACK = .TRUE.
      LOGICAL             :: INFND
      EXTERNAL INFND
!!
!! Assure we have only valid characters for a double precision constant.
!!
      VALIDC = (VERIFY(CRX(1:ILEN),NTEGER//".Dd+-") .EQ. 0)
      IF (.NOT.VALIDC) THEN
        DPFND = .FALSE.
        RETURN
      ENDIF
!!
!! To have a double precision value, a D must be present. 
!!
      CALL CUPPER(CRX(1:ILEN))
      KEXP = INDEX(CRX(1:ILEN),"D")
      IF (KEXP .EQ. 0) THEN
        DPFND = .FALSE.
        RETURN
      ENDIF
!!
!! Test for a single exponent character, for example, 1D6
!!
      ONEEXP = (INDEX(CRX(1:ILEN),"D") .EQ. INDEX(CRX(1:ILEN),"D",BACK))
      IF (.NOT.ONEEXP) THEN
        DPFND = .FALSE.
        RETURN
      ENDIF
!!
!! Test for one decimal point character, for example, 1.D6
!!
      ONEDPT = (INDEX(CRX(1:ILEN),".") .EQ. INDEX(CRX(1:ILEN),".",BACK))
      IF (.NOT.ONEDPT) THEN
        DPFND = .FALSE.
        RETURN
      ENDIF
!!
!! Divide CRX(1:ILEN) at the exponent D; it must be preceeded by and 
!! followed by numeric values, for example, 1D6
!!
      IF (KEXP.EQ.1 .OR. KEXP.EQ.ILEN) THEN
        DPFND = .FALSE.
        RETURN
      ENDIF
!!
!! Test for an exponent being at most a signed integer.
!!
      EXPINT = INFND(CRX(KEXP+1:ILEN),ILEN-KEXP)
      IF (.NOT.EXPINT) THEN
        DPFND = .FALSE.
        RETURN
      ENDIF
!!
!! Assume the simplest case first: at most an integer before D,
!! for example, +1D+2 or 1D-2
!!
      SIMPLE = INFND(CRX(1:KEXP-1),KEXP-1)
      IF (SIMPLE) THEN
        DPFND = .TRUE.
        RETURN
      ENDIF
!!
!! Last case; there must be a decimal point present before D, for
!! example, -.1D6. First insure there is at least an integer before D.
!!
      INTBFD = (SCAN(CRX(1:KEXP-1),NTEGER) .NE. 0)
      IF (.NOT.INTBFD) THEN
        DPFND = .FALSE.
        RETURN
      ENDIF
!!
!! The only thing left is to insure is that any sign present is in
!! the first location.
!!
      IMIN = INDEX(CRX(1:KEXP-1),"-")
      IPOS = INDEX(CRX(1:KEXP-1),"+")
      IF (IMIN.EQ.0 .AND. IPOS.EQ.0) THEN
        DPFND = .TRUE.
      ELSE IF (IMIN.EQ.1 .AND. IPOS.EQ.0) THEN
        DPFND = .TRUE.
      ELSE IF (IMIN.EQ.0 .AND. IPOS.EQ.1) THEN
        DPFND = .TRUE.
      ELSE
        DPFND = .FALSE.
      ENDIF
!!
      RETURN
      END
!!_
      LOGICAL FUNCTION CHFND ( CRX,ILEN )
!!
!! Module: CHaracter_FouND      Coded by: S W Key, 11-MAR-1989 15:11:53
!!     Copyright (c) by KEY Associates;  1-JAN-1999 11:42:38:00
!!
!! Called by: RCRDRD
!!
!! Purpose: Identify character string CRX as a keyword or variable name.
!!
!! Description: If CRX equals a character, return a value of .TRUE.
!! otherwise return .FALSE.
!!
!! Arguments.
      CHARACTER(*), INTENT(IN) :: CRX   ! String to be tested
      INTEGER,      INTENT(IN) :: ILEN  ! String length
!!
!! Local variables.
      CHARACTER(10), SAVE :: NTEGER = "1234567890"
      CHARACTER(27), SAVE :: UCALPH = "ABCDEFGHIJKLMNOPQRSTUVWXYZ_"
      CHARACTER(27), SAVE :: LCALPH = "abcdefghijklmnopqrstuvwxyz_"
!!
!! Look for acceptable leading character: alphabet character.
!!
      CHFND = (VERIFY(CRX(1:1),UCALPH//LCALPH) .EQ. 0)
!!
!! Verify that the remaining characters are valid.
!!
      IF (CHFND) THEN
        CHFND = (VERIFY(CRX(1:ILEN),NTEGER//UCALPH//LCALPH) .EQ. 0)
      ENDIF
!!
      RETURN
      END
!!_
      REAL(KIND(0E0)) FUNCTION R_VALUE ( NVAL )
!!
      USE shared_common_data
      USE value_
!!
      R_VALUE = 0.0
      IF (VALUE(NVAL)%VTYP .EQ. 'N') THEN
        R_VALUE = 0.0
      ELSE IF (VALUE(NVAL)%VTYP .EQ. 'R') THEN
        R_VALUE = VALUE(NVAL)%RVAL
      ELSE IF (VALUE(NVAL)%VTYP .EQ. 'D') THEN
        R_VALUE = REAL (VALUE(NVAL)%DVAL)
      ELSE IF (VALUE(NVAL)%VTYP .EQ. 'I') THEN
        R_VALUE = REAL (VALUE(NVAL)%IVAL)
      ELSE
!!
!! Unexpected value type found.
!!
        WRITE (MSG1,'(I8)') NVAL
        CALL USER_MESSAGE                                         &
     &    (                                                       &
     &    MSGL//'FATAL'//                                         &
     &    MSGL//'R_VALUE.001.00'//                                &
     &    MSGL//'VALUE Key Word (Entry 1): '//VALUE(1)%CVAL//     &
     &    MSGL//'Unexpected C_Value Found: '//VALUE(NVAL)%CVAL//  &
     &    MSGL//'VALUE Entry Number, NVAL: '//MSG1//              &
     &    MSGL//'VALUE Entry Type,   VTYP: '//VALUE(NVAL)%VTYP    &
     &    )
      ENDIF
!!
      RETURN
      END
!!_
      REAL(KIND(0D0)) FUNCTION D_VALUE ( NVAL )
!!
!! Copyright (c) by KEY Associates;  1-JAN-1999 11:42:38:00
!!
      USE shared_common_data
      USE value_
!!
!! Arguments.
      INTEGER, INTENT(IN) :: NVAL
!!
      D_VALUE = 0.0D+00
      IF (VALUE(NVAL)%VTYP .EQ. 'N') THEN
        D_VALUE = 0.0D+00
      ELSE IF (VALUE(NVAL)%VTYP .EQ. 'D') THEN
        D_VALUE = VALUE(NVAL)%DVAL
      ELSE IF (VALUE(NVAL)%VTYP .EQ. 'R') THEN
        D_VALUE = DBLE (VALUE(NVAL)%RVAL)
      ELSE IF (VALUE(NVAL)%VTYP .EQ. 'I') THEN
        D_VALUE = DBLE (VALUE(NVAL)%IVAL)
      ELSE
!!
!! Unexpected value type found.
!!
        WRITE (MSG1,'(I8)') NVAL
        CALL USER_MESSAGE                                         &
     &    (                                                       &
     &    MSGL//'FATAL'//                                         &
     &    MSGL//'D_VALUE.001.00'//                                &
     &    MSGL//'VALUE Key Word (Entry 1): '//VALUE(1)%CVAL//     &
     &    MSGL//'Unexpected D_Value Found: '//VALUE(NVAL)%CVAL//  &
     &    MSGL//'VALUE Entry Number, NVAL: '//MSG1//              &
     &    MSGL//'VALUE Entry Type,   VTYP: '//VALUE(NVAL)%VTYP    &
     &    )
      ENDIF
!!
      RETURN
      END
!!_
      INTEGER FUNCTION I_VALUE ( NVAL )
!!
!! Copyright (c) by KEY Associates;  1-JAN-1999 11:42:38:00
!!
      USE value_
!!
!! Arguments.
      INTEGER, INTENT(IN) :: NVAL
!!
!! Local variables.
      LOGICAL :: L_VALUE
      EXTERNAL   L_VALUE
!!
      I_VALUE = 0
      IF (VALUE(NVAL)%VTYP .EQ. 'N') THEN
        I_VALUE = 0
      ELSE IF (VALUE(NVAL)%VTYP .EQ. 'I') THEN
        I_VALUE = VALUE(NVAL)%IVAL
      ELSE IF (VALUE(NVAL)%VTYP .EQ. 'R') THEN
        I_VALUE = NINT (VALUE(NVAL)%RVAL)
      ELSE IF (VALUE(NVAL)%VTYP .EQ. 'D') THEN
        I_VALUE = NINT (VALUE(NVAL)%DVAL)
      ELSE IF (L_VALUE(NVAL)) THEN
        I_VALUE = 1
      ELSE
        I_VALUE = 0
      ENDIF
!!
      RETURN
      END
!!_
      CHARACTER*(*) FUNCTION C_VALUE ( NVAL )
!!
!! Copyright (c) by KEY Associates;  1-JAN-1999 11:42:38:00
!!
      USE shared_common_data
      USE value_
!!
!! Arguments.
      INTEGER, INTENT(IN) :: NVAL
!!
      IF (VALUE(NVAL)%VTYP .EQ. 'N') THEN
        C_VALUE = '<null!>'
      ELSE IF (VALUE(NVAL)%VTYP .EQ. 'C') THEN
        CALL CUPPER (VALUE(NVAL)%CVAL)
        C_VALUE = VALUE(NVAL)%CVAL
      ELSE
!!
!! Unknown value type found.
!!
        WRITE (MSG1,'(I8)') NVAL
        CALL USER_MESSAGE                                          &
     &    (                                                        &
     &    MSGL//'FATAL'//                                          &
     &    MSGL//'C_VALUE.001.00'//                                 &
     &    MSGL//'VALUE Key Word (Entry 1): '//VALUE(1)%CVAL//      &
     &    MSGL//'Unexpected C_Value Found: '//VALUE(NVAL)%CVAL//   &
     &    MSGL//'VALUE Entry Number, NVAL: '//MSG1//               &
     &    MSGL//'VALUE Entry Type,   VTYP: '//VALUE(NVAL)%VTYP     &
     &    )
      ENDIF
!!
      RETURN
      END
!!_
      LOGICAL FUNCTION L_VALUE ( NVAL )
!!
!! Copyright (c) by KEY Associates;  1-JAN-1999 11:42:38:00
!!
      USE shared_common_data
      USE value_
!!
!! Arguments.
      INTEGER, INTENT(IN) :: NVAL
!!
      L_VALUE = .FALSE.
      IF (VALUE(NVAL)%VTYP .EQ. 'N') THEN
        L_VALUE = .FALSE.
      ELSE IF (VALUE(NVAL)%VTYP .EQ. 'D') THEN
        L_VALUE = (VALUE(NVAL)%DVAL .NE. 0.0D+00)
      ELSE IF (VALUE(NVAL)%VTYP .EQ. 'R') THEN
        L_VALUE = (VALUE(NVAL)%RVAL .NE. 0.0E+00)
      ELSE IF (VALUE(NVAL)%VTYP .EQ. 'I') THEN
        L_VALUE = (VALUE(NVAL)%IVAL .NE. 0)
      ELSE IF (VALUE(NVAL)%VTYP.EQ.'C'.OR.VALUE(NVAL)%VTYP.EQ.'U') THEN
        CALL CUPPER (VALUE(NVAL)%CVAL)
        L_VALUE =                                      &
     &          (                                      &
     &          (INDEX(VALUE(NVAL)%CVAL,'T')  .NE. 0)  &
     &          .OR.                                   &
     &          (INDEX(VALUE(NVAL)%CVAL,'Y')  .NE. 0)  &
     &          .OR.                                   &
     &          (INDEX(VALUE(NVAL)%CVAL,'ON') .NE. 0)  &
     &          )
      ELSE
!!
!! Unexpected value type found.
!!
        WRITE (MSG1,'(I8)') NVAL
        CALL USER_MESSAGE                                          &
     &    (                                                        &
     &    MSGL//'FATAL'//                                          &
     &    MSGL//'L_VALUE.001.00'//                                 &
     &    MSGL//'VALUE Key Word (Entry 1): '//VALUE(1)%CVAL//      &
     &    MSGL//'Unexpected L_Value Found: '//VALUE(NVAL)%CVAL//   &
     &    MSGL//'VALUE Entry Number, NVAL: '//MSG1//               &
     &    MSGL//'VALUE Entry Type,   VTYP: '//VALUE(NVAL)%VTYP     &
     &    )
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE GETIRC (NVALS,LCRD,LPRT,EOF,IERR,TEXT,IRINDEX,ILCOUNT)
!!
!! Module: GET_Input_Record_Count  Coded by: S W Key, 5-APR-1989 21:22:11
!!     Copyright (c) by KEY Associates;  1-JAN-1999 11:42:38:00
!!
!! Purpose: Read next input record, get key word, and count values present.
!!
!! Free-Field Record Reader. This module expects to find and process
!! RCRD(1:80) = <value><delimiter>... It will read the next input record
!! and any continuation lines, ignoring any preceding or interposed
!! comment or blank (empty) lines. The grammar rules are as follows:
!!
!! 1. At the start ALL values are initialized to null.
!!
!! 2. Each input line is 80 characters or less in length.
!!
!! 3. Values are separated with delimiters. A delimiter is any of the
!!    following:
!!      a. a single tab surrounded by any number of spaces,
!!      b. a single comma surrounded by any number of spaces,
!!      c. a single equal surrounded by any number of spaces,
!!      d. a single asterisk surrounded by any number of spaces,
!!      e. any number of consecutive spaces not bounded by
!!         commas, tabs, equals, or asterisks,
!!      f. End-Of-Record character ($).
!!
!! 4. Blank lines and comment lines are ignored. Leading blanks before
!!    the first value occurring in a line are ignored.
!!
!! 5. All values must be separated by a delimiter. Repeated delimiters
!!    with the exception of blanks produce a null value for the parameter
!!    in that position.
!!
!! 6. An asterisk "*" is SIMULTANEOUSLY a delimiter and a continuation
!!    character meaning additional values will be found on the following
!!    line.
!!
!! 7. A dollar sign "$" denotes a comment. All characters to the right of
!!    the dollar sign are ignored.
!!      a. If a dollar sign is the first non-space character, the line is a
!!         comment line and processing proceeds to the next line for data.
!!      b. If a dollar sign is not the first non-space character, the data
!!         expected from the input record is assumed to be complete.
!!
!! 8. Any sequence of character constants, signed integer constants, and
!!    floating point constants is acceptable. Only the first 32 characters
!!    in a character constant are reported. If a floating point constant
!!    is supplied for an integer variable, NINT is used to perform a mode
!!    conversion.
!!
      USE value_
!!
!! Arguments.
      INTEGER,       INTENT(INOUT) :: NVALS   ! Num values expected/num returned
      INTEGER,       INTENT(IN)    :: LCRD    ! Input file logical unit number
      INTEGER,       INTENT(IN)    :: LPRT    ! Print file logical unit number
      LOGICAL,       INTENT(OUT)   :: EOF     ! Flag for End-Of-File
      INTEGER,       INTENT(INOUT) :: IERR    ! Error counter
      CHARACTER(80), INTENT(OUT)   :: TEXT    ! First line of input record
      INTEGER,       INTENT(OUT)   :: IRINDEX ! Input line number for keyword
      INTEGER,       INTENT(INOUT) :: ILCOUNT ! Input line counter
!!
!! Local variables.
      CHARACTER( 1), SAVE :: NUL = "N"  ! Null value indicator
      CHARACTER( 1), SAVE :: UNK = "U"  ! Unknown value indicator
      CHARACTER( 1), SAVE :: DOL = "$"  ! Comment/blank line value
      CHARACTER( 1), SAVE :: SPA = " "  ! Space value
      CHARACTER(81), SAVE :: DTS        ! Dots for error message
      CHARACTER(81), SAVE :: PTS        ! Up-arrows for error message
      CHARACTER(81)       :: UMK        ! Error message construction buffer
      CHARACTER(81)       :: RCRD       ! Input record buffer
      LOGICAL             :: EOR        ! Flag for end of record
      LOGICAL             :: CONT       ! Flag for continuation line
      LOGICAL             :: SPFND      ! Flag for space
      LOGICAL             :: CLINE      ! Flag for comment line
      LOGICAL             :: CHECK      ! Flag for unrecognized input entries
      LOGICAL             :: FRSTV      ! Flag for reading keyword, 1st value
!!
!! External functions.
      CHARACTER(32)       :: C_VALUE    ! Returns character value
      EXTERNAL               C_VALUE

      LOGICAL, SAVE :: FIRST = .TRUE.
!!
      IF (FIRST) THEN
        DO N = 1,81
          DTS(N:N) = '.'
          PTS(N:N) = '^'
        ENDDO
        FIRST = .FALSE.
      ENDIF
!!
!! Initialize key word value to "null."
!!
      VALUE(1)%VTYP = NUL
!!
!! Initialize count of values in record to zero.
!!
      NVALS = 0
!!
!! Read input lines until a non-comment or a non-blank line is found.
!!
      CLINE = .TRUE.
      DO WHILE (CLINE)
        EOF = .TRUE.
        ILCOUNT = ILCOUNT + 1
        READ (LCRD,'(A)',ERR=100) RCRD
        EOF = .FALSE.
 100    IF (EOF) RETURN
        RCRD(81:81) = '$'
        I = 0
        SPFND = .TRUE.
        DO WHILE (SPFND)
          I = I + 1
          SPFND = (RCRD(I:I) .EQ. SPA)
        ENDDO
        CLINE = (RCRD(I:I) .EQ. DOL)
      ENDDO
!!
!! Save first line of the input record (no continuation lines) in case a
!! character string (title or descriptive text) is expected.
!!
      TEXT = RCRD(1:80)
!!
!! Report line number in input file where the input record starts.
!!
      IRINDEX = ILCOUNT
!!
!! Extract values from RCRD -- key word plus a count.
!!
      I = 1
      N = 1
      ILEN  = 81
      EOR   = .FALSE.
      CONT  = .FALSE.
      FRSTV = .TRUE.
      DO WHILE (.NOT.EOR)
         IF (FRSTV) THEN
            CALL RCRDRD (RCRD(I:81),ILEN,1,INXT,CONT,EOR)
            FRSTV = .FALSE.
         ELSE
            CALL SKPTNV (RCRD(I:81),ILEN,INXT,CONT,EOR)
         ENDIF
!!
!! If a continuation delimiter was found, read the next line.
!! If the next line is a comment or blank line, read another
!! line.
!!
         IF (CONT) THEN
            CLINE = .TRUE.
            DO WHILE (CLINE)
              EOF = .TRUE.
              ILCOUNT = ILCOUNT + 1
              READ (LCRD,'(A)',ERR=200) RCRD
              EOF = .FALSE.
 200          IF (EOF) RETURN
              RCRD(81:81) = '$'
              I = 0
              SPFND = .TRUE.
              DO WHILE (SPFND)
                I = I + 1
                SPFND = (RCRD(I:I) .EQ. SPA)
              ENDDO
              CLINE = (RCRD(I:I) .EQ. DOL)
            ENDDO
            I = 1
            ILEN = 81
            CONT = .FALSE.
         ELSE
            I = I + INXT
            ILEN = ILEN - INXT
         ENDIF
!!
         N = N + 1
      ENDDO
      NVALS = N - 1
!!
!! If non-character key word was found (an unknown type), write a message and
!! increment the error counter.
!!
      IF (VALUE(1)%VTYP .EQ. UNK) THEN
        J = VALUE(1)%LOC
        K = INDEX(VALUE(1)%CVAL(1:32),SPA)
        IF (K .EQ. 0) THEN
          K = J + 31
        ELSE
          K = J + K - 1
        ENDIF
        UMK = DTS
        UMK(J:K) = PTS(J:K)
        WRITE (LPRT,'(/A/(18X,A/))')                       &
     &    '*** WARNING *** Message From Module: GETIRC',   &
     &    'Unknown Value Type Detected In Input Record:',  &
     &    '['//RCRD(1:81)//']',                            &
     &    '['//UMK//']'
        IERR = IERR + 1
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE SKPTNV (RCRD,ILEN,INXT,CONT,EOR)
!!
!! Module: SKiP_To_Next_Value    Coded by: S W Key, 5-APR-1989 21:24:36
!!     Copyright (c) by KEY Associates;  1-JAN-1999 11:42:38:00
!!
!! Called by: GETIRC
!!
!! Purpose: Find the next entry in the record and skip over it.
!!
!! This module expects to find and process RCRD(1:) = <value><delimiter>...
!! where <value>
!!      = a character constant
!!      = an integer constant, or
!!      = a real constant,
!! and <delimiter>
!!      = a single comma w/ any number of spaces,
!!      = a single tab w/ any number of spaces,
!!      = a single equal w/ any number of spaces,
!!      = a single asterisk w/ any number of spaces,
!!      = any number of consecutive spaces, or
!!      = an End-Of-Record character ($).
!!
!! This module will leave INXT pointing to the first position of the next
!! <value><delimiter> pair in RCRD.
!!
!! Arguments.
      CHARACTER(*), INTENT(IN)  :: RCRD  ! Character string
      INTEGER,      INTENT(IN)  :: ILEN  ! Length of character string RCRD
      INTEGER,      INTENT(OUT) :: INXT  ! Off-set to next entry
      LOGICAL,      INTENT(OUT) :: CONT  ! Flag for continuation
      LOGICAL,      INTENT(OUT) :: EOR   ! Flag for End-Of-Record
!!
!! Local variables.
      CHARACTER(1), SAVE :: NUL = "N"  ! Value flag, null     
      CHARACTER(1), SAVE :: CHR = "C"  ! Value flag, character
      CHARACTER(1), SAVE :: INT = "I"  ! Value flag, integer  
      CHARACTER(1), SAVE :: REA = "R"  ! Value flag, real     
      CHARACTER(1), SAVE :: UNK = "U"  ! Value flag, unknown  
      CHARACTER(1), SAVE :: TAB = " "  ! Tab                  
      CHARACTER(1), SAVE :: COM = ","  ! Comma                
      CHARACTER(1), SAVE :: EQU = "="  ! Equal                
      CHARACTER(1), SAVE :: SPA = " "  ! Space                
      CHARACTER(1), SAVE :: AST = "*"  ! Asterisk             
      CHARACTER(1), SAVE :: DOL = "$"  ! Dollar sign
      LOGICAL            :: FOUND      ! Flag for next delimiter 
      LOGICAL            :: ERFND      ! Flag for End-Of-Record
      LOGICAL            :: CNFND      ! Flag for continue delimiter * 
      LOGICAL            :: DLFND      ! Flag for space,comma,tab,*,$,=
      LOGICAL            :: TBFND      ! Flag for TaB_FouND    
      LOGICAL            :: CMFND      ! Flag for CoMma_FouND  
      LOGICAL            :: EQFND      ! Flag for EQual_FouND  
      LOGICAL            :: SPFND      ! Flag for SPace_FouND  
      LOGICAL            :: NTASP      ! Flag for NoT_A_SPace  
!!
      TAB = ACHAR(9)
!!
!! Step over leading spaces before <value>.
!!
      I = 0
      FOUND = .FALSE.
      DO WHILE (.NOT.FOUND .AND. I .LE. ILEN)
        I = I + 1
        FOUND = (RCRD(I:I) .NE. SPA)
      ENDDO
!!
!! Look for the next delimiter.
!!
      J = I
      I = I - 1
      FOUND = .FALSE.
      DO WHILE (.NOT.FOUND .AND. I .LE. ILEN)
        K = I
        I = I + 1
        FOUND = DLFND (RCRD(I:I))
      ENDDO
!!
!! Skip over RCRD(J:K). Locate the next delimiter. Look for one of the
!! following:
!!      1. An End-Of-Record character ($)
!!      2. One or more consecutive spaces
!!      3. A comma with or without leading spaces
!!      4. A tab with or without leading spaces
!!      5. An equal with or without leading spaces
!!      6. An asterisk with or without leading spaces
!!
      I = K
      FOUND = .FALSE.
      DO WHILE (.NOT.FOUND .AND. I .LE. ILEN)
        I = I + 1
        ERFND = (RCRD(I:I) .EQ. DOL)
        CNFND = (RCRD(I:I) .EQ. AST)
        CMFND = (RCRD(I:I) .EQ. COM)
        TBFND = (RCRD(I:I) .EQ. TAB)
        EQFND = (RCRD(I:I) .EQ. EQU)
        NTASP = (RCRD(I:I) .NE. SPA)
        FOUND = ERFND .OR. CNFND .OR. CMFND .OR. TBFND .OR. NTASP .OR. EQFND
      ENDDO
!!
      SPFND = NTASP .AND. .NOT.(CMFND .OR. TBFND .OR. EQFND)
      IF (SPFND) THEN
       INXT = I - 1
      ELSE
       INXT = I
      ENDIF
      CONT = CNFND
      EOR  = ERFND
!!
      RETURN
      END
