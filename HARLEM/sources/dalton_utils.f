C
C Axxiliary Utilities to communicate to Dalton 
C
C  /* Deck opendx */
      SUBROUTINE OPENDX (LUDX,NAME,NELEM,NREC,STATUS,LRDX,NBDX,OLDDX)
C
C 15-Jun-1985 hjaaj
C
C Revisions :  9-Dec-1987 hjaaj (Alliant version)
C
C Purpose:
C   Open files for direct access through WRITDX and READDX routines.
C   The ....DX routines enables direct access, even when the number
C   of elements per record (the logical record length) is greater
C   than the maximum physical record length.
C   >>> THIS IS MACHINE DEPENDENT <<<
C
C Input:
C  LUDX     file unit number
C  NELEM    number of integer words per logical record
C  NREC     number of logical records
C  STATUS   file status: 'OLD', 'NEW', or 'UNKNOWN'
C
C Output:
C  LRDX     physical record length (in integers)
C  NBDX     number of physical records per logical record
C  OLDDX    logical, true if old LUDX file was opened
C
C
      CHARACTER*(*) NAME, STATUS
      LOGICAL OLDDX
      PARAMETER (LUCMD = 5, LUPRI = 6)

      NBDX   = 1
      LRDX   = NELEM
C#ifdef SYS_CRAY
C
C     CRAY has 8 byte integers.
C
C      LRECL  = 8*LRDX
C#elif SYS_T3D
C	LRECL  = 8*LRDX
C#elif SYS_IRIX
C	LRECL  = LRDX
C#else
	LRECL  = 4*LRDX
C#endif
C
C     ILDX   = NREC*NBDX * LRECL for "initialsize" in bytes
C
      IF (STATUS .EQ. 'NEW') GO TO 300
      IF (STATUS .NE. 'OLD' .AND. STATUS .NE. 'UNKNOWN') GO TO 9000
C
C     OPEN OLD FILE
C
         OPEN(LUDX,FILE=NAME,STATUS='OLD',FORM='UNFORMATTED',ERR=300,
     *        ACCESS='DIRECT',RECL=LRECL)
         OLDDX = .TRUE.
      GO TO 600
C
  300 CONTINUE
      IF (STATUS .EQ. 'OLD') GO TO 9100
C
C     OPEN NEW FILE
C
         OPEN(LUDX,FILE=NAME,STATUS='NEW',FORM='UNFORMATTED',
     *        ACCESS='DIRECT',RECL=LRECL)
         OLDDX = .FALSE.
  600  CONTINUE
      RETURN
C
C error branches
C
 9000 CONTINUE
      WRITE (LUPRI,'(//A,A/A,I5)')
     *   ' *** ERROR (OPENDX) INVALID STATUS KEYWORD: ',STATUS,
     *   '                    FILE NUMBER =',LUDX
C     CALL QTRACE(LUPRI)
C      CALL QUIT('*** ERROR (OPENDX) INVALID STATUS KEYWORD')
C
 9100 CONTINUE
      WRITE (LUPRI,'(//A/A,I5/A)')
     *   ' *** ERROR (OPENDX) OLD FILE NOT FOUND',
     *   '                    FILE NUMBER =',LUDX,
     *   ' --- or wrong record length on old file.'
C      CALL QTRACE(LUPRI)
C      CALL QUIT('*** ERROR (OPENDX) FILE NOT FOUND')
C
C end of OPENDX
C
      END
C  /* Deck finddx */
      LOGICAL FUNCTION FINDDX(LU,LRDX,I,LEN,IVEC)
C
C 27-Jun-1985 Hans Jorgen Aa. Jensen
C
C For direct access find record,
C when LEN may be greater than maximum record length.
C
      INTEGER IVEC(LEN)
      IF (LEN .LE. LRDX) THEN
         READ (LU, REC=I, IOSTAT=IOS) IVEC
         IF (IOS .NE. 0) GO TO 900
      ELSE
         NBUF = (LEN-1)/LRDX + 1
         IREC = 1 + NBUF*(I-1)
         JADD = 0
         DO 100 IBUF = 1,NBUF-1
            READ (LU, REC=IREC, IOSTAT=IOS) (IVEC(JADD+J), J = 1,LRDX)
            IF (IOS .NE. 0) GO TO 900
            IREC = IREC + 1
            JADD = JADD + LRDX
  100    CONTINUE
         READ (LU, REC=IREC, IOSTAT=IOS) (IVEC(J), J = JADD+1,LEN)
         IF (IOS .NE. 0) GO TO 900
      END IF
      FINDDX = .TRUE.
      RETURN
C
  900 CONTINUE
      FINDDX = .FALSE.
      RETURN
      END
C  /* Deck readdx */
      SUBROUTINE READDX(LU,LRDX,I,LEN,IVEC)
C
C 30-Apr-1985 Hans Jorgen Aa. Jensen
C
C For direct access read
C when LEN may be greater than maximum record length.
C
      INTEGER IVEC(LEN)
      IF (LEN .LE. LRDX) THEN
         READ (LU, REC = I) IVEC
      ELSE
         NBUF = (LEN-1)/LRDX + 1
         IREC = 1 + NBUF*(I-1)
         JADD = 0
         DO 100 IBUF = 1,NBUF-1
            READ (LU, REC = IREC) (IVEC(JADD+J), J = 1,LRDX)
            IREC = IREC + 1
            JADD = JADD + LRDX
  100    CONTINUE
         READ (LU, REC = IREC) (IVEC(J), J = JADD+1,LEN)
      END IF
      RETURN
      END

C  /* Deck writdx */
      SUBROUTINE WRITDX(LU,LRDX,I,LEN,IVEC)
C
C 30-Apr-1985 Hans Jorgen Aa. Jensen
C
C For direct access write
C when LEN may be greater than maximum record length.
C
      INTEGER IVEC(LEN)
      IF (LEN .LE. LRDX) THEN
         WRITE (LU, REC = I) IVEC
      ELSE
         NBUF = (LEN-1)/LRDX + 1
         IREC = 1 + NBUF*(I-1)
         JADD = 0
         DO 100 IBUF = 1,NBUF-1
            WRITE (LU, REC = IREC) (IVEC(JADD+J), J = 1,LRDX)
            IREC = IREC + 1
            JADD = JADD + LRDX
  100    CONTINUE
         WRITE (LU, REC = IREC) (IVEC(J), J = JADD+1,LEN)
      END IF
      RETURN
      END


