*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     File  f06.f
*     A subset of the NAG F06 Chapter with some modifications.
*     The routines perform the same function as the NAG F06 routines.
*
*                         Level 0  F06  Scalar routines
*                         -------  ---- ---------------
*     f06aaz+         f06baf/drot3g+  f06bcf/dcsg+    f06blf/ddiv+
*     f06bmf/dnorm+
*
*                         Level 1  BLAS Vector routines
*                         -------  ---- ---------------
*     daxpy           dcopy           ddot            dnrm2    
*     dscal           dswap           idamax          drot
*
*                         Level 1  F06  Vector routines
*                         -------  ---  ---------------
*     f06dbf/iload    f06fbf/dload    dddiv           f06fcf/ddscl 
*     f06dff/icopy    f06fjf/dssq+    f06fkf          f06flf/dcond
*     f06klf/idrank+  f06fqf          f06frf/dgrfg+
*
*                         Level 2  BLAS Matrix-vector routines
*                         -------  ---  ----------------------
*     dgemv           dger            dsymv           dsyr     
*     dtrmv           dtrsv
*
*                         Level 2  F06  Matrix routines
*                         -------  ---  ---------------
*     f06qff          f06qgf          f06qhf          f06qjf   
*     f06qkf          f06qnf          f06qrf          f06qsf   
*     f06qtf          f06qvf          f06qwf          f06qxf   
*     f06qzf
*
*    +Differs from the Nag F06 version.
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

*UPTODATE F06AAZTEXT
      SUBROUTINE F06AAZ( SRNAME, INFO )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     Modified by PEG 9/25/88.
C     .. Scalar Arguments ..
      INTEGER            INFO
      CHARACTER*13       SRNAME
C     ..
C
C  Purpose
C  =======
C
C  F06AAZ  is an error handler for the Level 2 BLAS routines.
C
C  It is called by the Level 2 BLAS routines if an input parameter is
C  invalid.
C
C  Parameters
C  ==========
C
C  SRNAME - CHARACTER*13.
C           On entry, SRNAME specifies the name of the routine which
C           called F06AAZ.
C
C  INFO   - INTEGER.
C           On entry, INFO specifies the position of the invalid
C           parameter in the parameter-list of the calling routine.
C
C
C  Auxiliary routine for Level 2 Blas.
C
C  Written on 20-July-1986.
C
C     .. Local Scalars ..
      CHARACTER*80       REC (1)
C     ..
C     .. Executable Statements ..
C      WRITE (REC (1),99999) SRNAME, INFO
C
      RETURN
C
99999 FORMAT ( ' ** On entry to ', A13, ' parameter number ', I2,
     $         ' had an illegal value' )
C
C     End of F06AAZ.
C
      END
** END OF F06AAZTEXT
      SUBROUTINE F06BAF( X, Y, CS, SN )

      DOUBLE PRECISION   X, Y, CS, SN

      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/

C
C  Note: f06baf/drot3g is different from the Nag routine f06baf.
C
C  DROT3G  generates a plane rotation that reduces the vector (X, Y) to
C  the vector (A, 0),  where A is defined as follows...
C
C     If both X and Y are negligibly small, or
C     if Y is negligible relative to Y,
C     then  A = X,  and the identity rotation is returned.
C
C     If X is negligible relative to Y,
C     then  A = Y,  and the swap rotation is returned.
C
C     Otherwise,  A = sign(X) * sqrt( X**2 + Y**2 ).
C
C  In all cases,  X and Y are overwritten by A and 0,  and CS will lie
C  in the closed interval (0, 1).  Also,  the absolute value of CS and
C  SN (if nonzero) will be no less than the machine precision,  EPS.
C
C  DROT3G  guards against overflow and underflow.
C  It is assumed that  FLMIN .lt. EPS**2  (i.e.  RTMIN .lt. EPS).
C
C  Systems Optimization Laboratory, Stanford University.
C  Original version dated January 1982.
C  F77 version dated 28-June-1986.
C  This version of DROT3G dated 28-June-1986.
C
      DOUBLE PRECISION   A, B, EPS, ONE, RTMIN, ZERO
      LOGICAL            FIRST
      INTRINSIC          ABS, MAX, SQRT
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0 )

      SAVE               FIRST , EPS   , RTMIN
      DATA               FIRST / .TRUE. /

      IF( FIRST )THEN
         FIRST = .FALSE.
         EPS    = WMACH(3)
         RTMIN  = WMACH(6)
      END IF

      IF (Y .EQ. ZERO) THEN

         CS = ONE
         SN = ZERO

      ELSE IF (X .EQ. ZERO) THEN

         CS = ZERO
         SN = ONE
         X  = Y

      ELSE

         A      = ABS(X)
         B      = ABS(Y)
         IF (MAX(A,B) .LE. RTMIN) THEN
            CS = ONE
            SN = ZERO
         ELSE
            IF (A .GE. B) THEN
               IF (B .LE. EPS*A) THEN
                  CS = ONE
                  SN = ZERO
                  GO TO 900
               ELSE
                  A  = A * SQRT( ONE + (B/A)**2 )
               END IF
            ELSE
               IF (A .LE. EPS*B) THEN
                  CS = ZERO
                  SN = ONE
                  X  = Y
                  GO TO 900
               ELSE
                  A  = B * SQRT( ONE + (A/B)**2 )
               END IF
            END IF
            IF (X .LT. ZERO) A = - A
            CS = X/A
            SN = Y/A
            X  = A
         END IF
      END IF

  900 Y  = ZERO

      RETURN

*     End of  F06BAF (DROT3G).

      END
*UPTODATE F06BCFTEXT
      SUBROUTINE F06BCF( T, C, S )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 13 REVISED. IER-602 (MAR 1988).
C     .. Scalar Arguments ..
      DOUBLE PRECISION   C, S, T
C     ..
C
C  F06BCF returns values c and s such that
C
C     c = cos( theta ),   s = sin( theta )
C
C  for a given value of
C
C     t = tan( theta ).
C
C  c is always non-negative and s has the same sign as t, so that
C
C     c = 1.0/sqrt( 1.0 + t**2 ),   s = t/sqrt( 1.0 + t**2 ).
C
C  Nag Fortran 77 O( 1 ) basic linear algebra routine.
C
C  -- Written on 28-February-1986.
C     Sven Hammarling, Nag Central Office.
C  -- Modified 19-August-1987.
C     Sven Hammarling and Jeremy Du Croz, Nag Central Office.
C        No longer sets s to zero when t is less than eps.
C  -- Modified 24-July-1991.
C     Philip E. Gill, UCSD.
C        Modified to call mchpar instead of x02ajf
C
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
C     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER        ( ONE = 1.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   ABST, EPS, REPS, RRTEPS, RTEPS
      LOGICAL            FIRST
C     .. External Functions ..
C+    DOUBLE PRECISION   X02AJF
C+    EXTERNAL           X02AJF
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT
C     .. Save statement ..
      SAVE               FIRST, EPS, REPS, RTEPS, RRTEPS
C     .. Data statements ..
      DATA               FIRST/ .TRUE. /
C     ..
C     .. Executable Statements ..
      IF( FIRST )THEN
         FIRST  = .FALSE.
         EPS    = WMACH(3)
C+       EPS    = X02AJF( )
         REPS   =  1/EPS
         RTEPS  =  SQRT( EPS )
         RRTEPS =  1/RTEPS
      END IF
C
      ABST = ABS( T )
      IF( ABST.LT.RTEPS )THEN
         C = ONE
         S = T
      ELSE IF( ABST.GT.RRTEPS )THEN
         C = 1/ABST
         S = SIGN( ONE, T )
      ELSE
         C = 1/SQRT( 1 + ABST**2 )
         S = C*T 
      END IF
C
      RETURN
C
C     End of F06BCF. ( SCSG )
C
      END    
** END OF F06BCFTEXT
*UPTODATE F06BLFTEXT
      DOUBLE PRECISION FUNCTION F06BLF( A, B, FAIL )
      DOUBLE PRECISION          DDIV
      ENTRY                     DDIV  ( A, B, FAIL )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     Modified by PEG 9/25/88.
C     .. Scalar Arguments ..
      DOUBLE PRECISION                  A, B
      LOGICAL                           FAIL
C     ..
C
C  F06BLF returns the value div given by
C
C     div = ( a/b                 if a/b does not overflow,
C           (
C           ( 0.0                 if a .eq. 0.0,
C           (
C           ( sign( a/b )*flmax   if a .ne. 0.0  and a/b would overflow,
C
C  where  flmax  is a large value, via the function name. In addition if
C  a/b would overflow then  fail is returned as true, otherwise  fail is
C  returned as false.
C
C  Note that when  a and b  are both zero, fail is returned as true, but
C  div  is returned as  0.0. In all other cases of overflow  div is such
C  that  abs( div ) = flmax.
C
C  When  b = 0  then  sign( a/b )  is taken as  sign( a ).
C
C  Nag Fortran 77 O( 1 ) basic linear algebra routine.
C
C  -- Written on 26-October-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION      ONE         , ZERO
      PARAMETER           ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION      ABSB, DIV, FLMAX, FLMIN
      LOGICAL               FIRST
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
C     .. Intrinsic Functions ..
      INTRINSIC             ABS, SIGN
C     .. Save statement ..
      SAVE                  FIRST, FLMIN, FLMAX
C     .. Data statements ..
      DATA                  FIRST/ .TRUE. /
C     ..
C     .. Executable Statements ..
      IF( A.EQ.ZERO )THEN
         DIV = ZERO
         IF( B.EQ.ZERO )THEN
            FAIL = .TRUE.
         ELSE
            FAIL = .FALSE.
         END IF
      ELSE
C
         IF( FIRST )THEN
            FIRST  = .FALSE.
            FLMIN  = WMACH( 5 )
            FLMAX  = WMACH( 7 )
         END IF
C
         IF( B.EQ.ZERO )THEN
            DIV  =  SIGN( FLMAX, A )
            FAIL = .TRUE.
         ELSE
            ABSB = ABS( B )
            IF( ABSB.GE.ONE )THEN
               FAIL = .FALSE.
               IF( ABS( A ).GE.ABSB*FLMIN )THEN
                  DIV = A/B
               ELSE
                  DIV = ZERO
               END IF
            ELSE
               IF( ABS( A ).LE.ABSB*FLMAX )THEN
                  FAIL = .FALSE.
                  DIV  =  A/B
               ELSE
                  FAIL = .TRUE.
                  DIV  = FLMAX
                  IF( ( ( A.LT.ZERO ).AND.( B.GT.ZERO ) ).OR.
     $                ( ( A.GT.ZERO ).AND.( B.LT.ZERO ) )     )
     $               DIV = -DIV
               END IF
            END IF
         END IF
      END IF
C
      F06BLF = DIV
      RETURN
C
C     End of F06BLF. ( DDIV )
C
      END
** END OF F06BLFTEXT
*UPTODATE F06BMFTEXT
      DOUBLE PRECISION FUNCTION F06BMF( SCALE, SSQ )
      DOUBLE PRECISION          DNORM
      ENTRY                     DNORM ( SCALE, SSQ )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     Modified by PEG 9/25/88.
C     .. Scalar Arguments ..
      DOUBLE PRECISION                  SCALE, SSQ
C     ..
C
C  F06BMF returns the value norm given by
C
C     norm = ( scale*sqrt( ssq ), scale*sqrt( ssq ) .lt. flmax
C            (
C            ( flmax,             scale*sqrt( ssq ) .ge. flmax
C
C  via the function name.
C
C
C  Nag Fortran 77 O( 1 ) basic linear algebra routine.
C
C  -- Written on 22-October-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      DOUBLE PRECISION      FLMAX, NORM, SQT
      LOGICAL               FIRST
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
C     .. Intrinsic Functions ..
      INTRINSIC             SQRT
C     .. Save statement ..
      SAVE                  FIRST, FLMAX
C     .. Data statements ..
      DATA                  FIRST/ .TRUE. /
C     ..
C     .. Executable Statements ..
      IF( FIRST )THEN
         FIRST = .FALSE.
         FLMAX = WMACH( 7 )
      END IF
C
      SQT = SQRT( SSQ )
      IF( SCALE.LT.FLMAX/SQT )THEN
         NORM = SCALE*SQT
      ELSE
         NORM = FLMAX
      END IF
C
      F06BMF = NORM
      RETURN
C
C     End of F06BMF. ( DNORM )
C
      END
** END OF F06BMFTEXT
*UPTODATE F06ECFTEXT
      SUBROUTINE F06ECF( N, ALPHA, X, INCX, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     Modified by PEG 9/25/88.
C     .. Entry Points ..
      ENTRY      DAXPY ( N, ALPHA, X, INCX, Y, INCY )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * ), Y( * )
C     ..
C
C  F06ECF performs the operation
C
C     y := alpha*x + y
C
C
C  Nag Fortran 77 version of the Blas routine DAXPY.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 3-September-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      INTEGER            I, IX, IY
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ALPHA.NE.ZERO )THEN
            IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
               DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
                  Y( IX ) = ALPHA*X( IX ) + Y( IX )
   10          CONTINUE
            ELSE
               IF( INCY.GE.0 )THEN
                  IY = 1
               ELSE
                  IY = 1 - ( N - 1 )*INCY
               END IF
               IF( INCX.GT.0 )THEN
                  DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
                     Y( IY ) = ALPHA*X( IX ) + Y( IY )
                     IY      = IY            + INCY
   20             CONTINUE
               ELSE
                  IX = 1 - ( N - 1 )*INCX
                  DO 30, I = 1, N
                     Y( IY ) = ALPHA*X( IX ) + Y( IY )
                     IX      = IX            + INCX
                     IY      = IY            + INCY
   30             CONTINUE
               END IF
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06ECF. ( DAXPY )
C
      END
** END OF F06ECFTEXT
*UPTODATE F06EFFTEXT
      SUBROUTINE F06EFF( N, X, INCX, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      ENTRY      DCOPY ( N, X, INCX, Y, INCY )
C     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * ), Y( * )
C     ..
C
C  F06EFF performs the operation
C
C     y := x
C
C
C  Nag Fortran 77 version of the Blas routine DCOPY.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 26-November-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      INTEGER            I, IX, IY
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ( INCX.EQ.INCY ).AND.( INCY.GT.0 ) )THEN
            DO 10, IY = 1, 1 + ( N - 1 )*INCY, INCY
               Y( IY ) = X( IY )
   10       CONTINUE
         ELSE
            IF( INCX.GE.0 )THEN
               IX = 1
            ELSE
               IX = 1 - ( N - 1 )*INCX
            END IF
            IF( INCY.GT.0 )THEN
               DO 20, IY = 1, 1 + ( N - 1 )*INCY, INCY
                  Y( IY ) = X( IX )
                  IX      = IX      + INCX
   20          CONTINUE
            ELSE
               IY = 1 - ( N - 1 )*INCY
               DO 30, I = 1, N
                  Y( IY ) = X( IX )
                  IY      = IY      + INCY
                  IX      = IX      + INCX
   30          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06EFF. ( DCOPY )
C
      END
** END OF F06EFFTEXT
*UPTODATE F06EAFTEXT
      DOUBLE PRECISION FUNCTION F06EAF( N, X, INCX, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      DOUBLE PRECISION          DDOT
      ENTRY                     DDOT  ( N, X, INCX, Y, INCY )
C     .. Scalar Arguments ..
      INTEGER                           INCX, INCY, N
C     .. Array Arguments ..
      DOUBLE PRECISION                  X( * ), Y( * )
C     ..
C
C  F06EAF returns the value
C
C     F06EAF = x'y
C
C
C  Nag Fortran 77 version of the Blas routine DDOT.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 21-September-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION      ZERO
      PARAMETER           ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION      SUM
      INTEGER               I, IX, IY
C     ..
C     .. Executable Statements ..
      SUM = ZERO
      IF( N.GT.0 )THEN
         IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
            DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
               SUM = SUM + X( IX )*Y( IX )
   10       CONTINUE
         ELSE
            IF( INCY.GE.0 )THEN
               IY = 1
            ELSE
               IY = 1 - ( N - 1 )*INCY
            END IF
            IF( INCX.GT.0 )THEN
               DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
                  SUM = SUM + X( IX )*Y( IY )
                  IY  = IY  + INCY
   20          CONTINUE
            ELSE
               IX = 1 - ( N - 1 )*INCX
               DO 30, I = 1, N
                  SUM = SUM + X( IX )*Y( IY )
                  IX  = IX  + INCX
                  IY  = IY  + INCY
   30          CONTINUE
            END IF
         END IF
      END IF
C
      F06EAF = SUM
      RETURN
C
C     End of F06EAF. ( DDOT )
C
      END
** END OF F06EAFTEXT
*UPTODATE F06EJFTEXT
      DOUBLE PRECISION FUNCTION F06EJF( N, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      DOUBLE PRECISION          DNRM2
      ENTRY                     DNRM2 ( N, X, INCX )
C     .. Scalar Arguments ..
      INTEGER                           INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION                  X( * )
C     ..
C
C  F06EJF returns the euclidean norm of a vector via the function
C  name, so that
C
C     F06EJF := sqrt( x'*x )
C
C
C  Nag Fortran 77 version of the Blas routine DNRM2.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 25-October-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION      ONE         , ZERO
      PARAMETER           ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION      NORM, SCALE, SSQ
C     .. External Functions ..
      DOUBLE PRECISION      F06BMF
      EXTERNAL              F06BMF
C     .. External Subroutines ..
      EXTERNAL              F06FJF
C     .. Intrinsic Functions ..
      INTRINSIC             ABS
C     ..
C     .. Executable Statements ..
      IF( N.LT.1 )THEN
         NORM  = ZERO
      ELSE IF( N.EQ.1 )THEN
         NORM  = ABS( X( 1 ) )
      ELSE
         SCALE = ZERO
         SSQ   = ONE
         CALL F06FJF( N, X, INCX, SCALE, SSQ )
         NORM  = F06BMF( SCALE, SSQ )
      END IF
C
      F06EJF = NORM
      RETURN
C
C     End of F06EJF. ( DNRM2 )
C
      END
** END OF F06EJFTEXT
*UPTODATE F06EDFTEXT
      SUBROUTINE F06EDF( N, ALPHA, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      ENTRY      DSCAL ( N, ALPHA, X, INCX )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
C     ..
C
C  F06EDF performs the operation
C
C     x := alpha*x
C
C
C  Nag Fortran 77 version of the Blas routine DSCAL.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 26-November-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      INTEGER            IX
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ALPHA.EQ.ZERO )THEN
            DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = ZERO
   10       CONTINUE
         ELSE IF( ALPHA.EQ.( -ONE ) )THEN
            DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = -X( IX )
   20       CONTINUE
         ELSE IF( ALPHA.NE.ONE )THEN
            DO 30, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = ALPHA*X( IX )
   30       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06EDF. ( DSCAL )
C
      END
** END OF F06EDFTEXT
*UPTODATE F06EGFTEXT
      SUBROUTINE F06EGF( N, X, INCX, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      ENTRY      DSWAP ( N, X, INCX, Y, INCY )
C     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * ), Y( * )
C     ..
C
C  F06EGF performs the operations
C
C     temp := x,   x := y,   y := temp.
C
C
C  Nag Fortran 77 version of the Blas routine DSWAP.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 26-November-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, IX, IY
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ( INCX.EQ.INCY ).AND.( INCY.GT.0 ) )THEN
            DO 10, IY = 1, 1 + ( N - 1 )*INCY, INCY
               TEMP    = X( IY )
               X( IY ) = Y( IY )
               Y( IY ) = TEMP
   10       CONTINUE
         ELSE
            IF( INCX.GE.0 )THEN
               IX = 1
            ELSE
               IX = 1 - ( N - 1 )*INCX
            END IF
            IF( INCY.GT.0 )THEN
               DO 20, IY = 1, 1 + ( N - 1 )*INCY, INCY
                  TEMP    = X( IX )
                  X( IX ) = Y( IY )
                  Y( IY ) = TEMP
                  IX      = IX      + INCX
   20          CONTINUE
            ELSE
               IY = 1 - ( N - 1 )*INCY
               DO 30, I = 1, N
                  TEMP    = X( IX )
                  X( IX ) = Y( IY )
                  Y( IY ) = TEMP
                  IY      = IY      + INCY
                  IX      = IX      + INCX
   30          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06EGF. ( DSWAP )
C
      END
** END OF F06EGFTEXT
*UPTODATE F06JLFTEXT
      INTEGER FUNCTION F06JLF( N, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      INTEGER          IDAMAX
      ENTRY            IDAMAX( N, X, INCX )
C     .. Scalar Arguments ..
      INTEGER                  INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION         X( * )
C     ..
C
C  F06JLF returns the smallest value of i such that
C
C     abs( x( i ) ) = max( abs( x( j ) ) )
C                      j
C
C
C  Nag Fortran 77 version of the Blas routine IDAMAX.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 31-May-1983.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      DOUBLE PRECISION         XMAX
      INTEGER                  I, IMAX, IX
C     .. Intrinsic Functions ..
      INTRINSIC                ABS
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IMAX = 1
         IF( N.GT.1 )THEN
            XMAX = ABS( X( 1 ) )
            IX   = 1
            DO 10, I = 2, N
               IX = IX + INCX
               IF( XMAX.LT.ABS( X( IX ) ) )THEN
                  XMAX = ABS( X( IX ) )
                  IMAX = I
               END IF
   10       CONTINUE
         END IF
      ELSE
         IMAX = 0
      END IF
C
      F06JLF = IMAX
      RETURN
C
C     End of F06JLF. ( IDAMAX )
C
      END
** END OF F06JLFTEXT
*UPTODATE F06EPFTEXT
      SUBROUTINE F06EPF( N, X, INCX, Y, INCY, C, S )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      ENTRY      DROT  ( N, X, INCX, Y, INCY, C, S )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   C, S
      INTEGER            INCX, INCY, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * ), Y( * )
C     ..
C
C  F06EPF performs the plane rotation
C
C     ( x  y ) = ( x  y )*( c  -s ).
C                         ( s   c )
C
C
C  Nag Fortran 77 version of the Blas routine DROT.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 23-January-1984.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP1
      INTEGER            I, IX, IY
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ( S.NE.ZERO ).OR.( C.NE.ONE ) )THEN
            IF( ( C.EQ.ZERO ).AND.( S.EQ.ONE ) )THEN
               IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
                  DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
                     TEMP1   = -X( IX )
                     X( IX ) =  Y( IX )
                     Y( IX ) =  TEMP1
   10             CONTINUE
               ELSE
                  IF( INCY.GE.0 )THEN
                     IY = 1
                  ELSE
                     IY = 1 - ( N - 1 )*INCY
                  END IF
                  IF( INCX.GT.0 )THEN
                     DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
                        TEMP1   = -X( IX )
                        X( IX ) =  Y( IY )
                        Y( IY ) =  TEMP1
                        IY      =  IY       + INCY
   20                CONTINUE
                  ELSE
                     IX = 1 - ( N - 1 )*INCX
                     DO 30, I = 1, N
                        TEMP1   = -X( IX )
                        X( IX ) =  Y( IY )
                        Y( IY ) =  TEMP1
                        IX      =  IX      + INCX
                        IY      =  IY      + INCY
   30                CONTINUE
                  END IF
               END IF
            ELSE IF( ( C.EQ.ZERO ).AND.( S.EQ.( -ONE ) ) )THEN
               IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
                  DO 40, IX = 1, 1 + ( N - 1 )*INCX, INCX
                     TEMP1   =  X( IX )
                     X( IX ) = -Y( IX )
                     Y( IX ) =  TEMP1
   40             CONTINUE
               ELSE
                  IF( INCY.GE.0 )THEN
                     IY = 1
                  ELSE
                     IY = 1 - ( N - 1 )*INCY
                  END IF
                  IF( INCX.GT.0 )THEN
                     DO 50, IX = 1, 1 + ( N - 1 )*INCX, INCX
                        TEMP1   =  X( IX )
                        X( IX ) = -Y( IY )
                        Y( IY ) =  TEMP1
                        IY      =  IY       + INCY
   50                CONTINUE
                  ELSE
                     IX = 1 - ( N - 1 )*INCX
                     DO 60, I = 1, N
                        TEMP1   =  X( IX )
                        X( IX ) = -Y( IY )
                        Y( IY ) =  TEMP1
                        IX      =  IX      + INCX
                        IY      =  IY      + INCY
   60                CONTINUE
                  END IF
               END IF
            ELSE
               IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
                  DO 70, IX = 1, 1 + ( N - 1 )*INCX, INCX
                     TEMP1   = X( IX )
                     X( IX ) = S*Y( IX ) + C*TEMP1
                     Y( IX ) = C*Y( IX ) - S*TEMP1
   70             CONTINUE
               ELSE
                  IF( INCY.GE.0 )THEN
                     IY = 1
                  ELSE
                     IY = 1 - ( N - 1 )*INCY
                  END IF
                  IF( INCX.GT.0 )THEN
                     DO 80, IX = 1, 1 + ( N - 1 )*INCX, INCX
                        TEMP1   = X( IX )
                        X( IX ) = S*Y( IY ) + C*TEMP1
                        Y( IY ) = C*Y( IY ) - S*TEMP1
                        IY      = IY        + INCY
   80                CONTINUE
                  ELSE
                     IX = 1 - ( N - 1 )*INCX
                     DO 90, I = 1, N
                        TEMP1   = X( IX )
                        X( IX ) = S*Y( IY ) + C*TEMP1
                        Y( IY ) = C*Y( IY ) - S*TEMP1
                        IX      = IX        + INCX
                        IY      = IY        + INCY
   90                CONTINUE
                  END IF
               END IF
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06EPF. ( DROT )
C
      END
** END OF F06EPFTEXT
*UPTODATE F06DBFTEXT
      SUBROUTINE F06DBF( N, CONST, X, INCX )
      ENTRY      ILOAD ( N, CONST, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      INTEGER            CONST, INCX, N
C     .. Array Arguments ..
      INTEGER            X( * )
C     ..
C                      
C  F06DBF performs the operation
C
C     x = const*e,   e' = ( 1  1 ... 1 ).
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 18-February-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      INTEGER            IX
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( CONST.NE.0 )THEN
            DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = CONST
   10       CONTINUE
         ELSE
            DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = 0
   20       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06DBF. ( ILOAD )
C
      END
** END OF F06DBFTEXT
*UPTODATE F06FBFTEXT
      SUBROUTINE F06FBF( N, CONST, X, INCX )
      ENTRY      DLOAD ( N, CONST, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION   CONST
      INTEGER            INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
C     ..
C
C  F06FBF performs the operation
C
C     x = const*e,   e' = ( 1  1 ... 1 ).
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 22-September-1983.
C     Sven Hammarling, Nag Central Office.
C
C                      
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      INTEGER            IX
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( CONST.NE.ZERO )THEN
            DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = CONST
   10       CONTINUE
         ELSE
            DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = ZERO
   20       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06FBF. ( DLOAD )
C
      END
** END OF F06FBFTEXT
*UPTODATE F06FCFTEXT
      SUBROUTINE F06FCF( N, D, INCD, X, INCX )
      ENTRY      DDSCL ( N, D, INCD, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      INTEGER            INCD, INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), X( * )
C     ..
C
C  F06FCF performs the operation
C
C     x := diag( d )*x
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 22-September-1983.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      INTEGER            I, ID, IX
C     .. External Subroutines ..
      EXTERNAL           DSCAL
C     .. Intrinsic Functions ..
      INTRINSIC          ABS
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ( INCD.EQ.0 ).AND.( INCX.NE.0 ) )THEN
            CALL DSCAL( N, D( 1 ), X, ABS( INCX ) )
         ELSE IF( ( INCD.EQ.INCX ).AND.( INCD.GT.0 ) )THEN
            DO 10, ID = 1, 1 + ( N - 1 )*INCD, INCD
               X( ID ) = D( ID )*X( ID )
   10       CONTINUE
         ELSE
            IF( INCX.GE.0 )THEN
               IX = 1
            ELSE
               IX = 1 - ( N - 1 )*INCX
            END IF
            IF( INCD.GT.0 )THEN
               DO 20, ID = 1, 1 + ( N - 1 )*INCD, INCD
                  X( IX ) = D( ID )*X( IX )
                  IX      = IX              + INCX
   20          CONTINUE
            ELSE
               ID = 1 - ( N - 1 )*INCD
               DO 30, I = 1, N
                  X( IX ) = D( ID )*X( IX )
                  ID      = ID              + INCD
                  IX      = IX              + INCX
   30          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06FCF. ( DDSCL )
C            
      END
** END OF F06FCFTEXT
      subroutine dddiv ( n, d, incd, x, incx )

      implicit           double precision (a-h,o-z)
      double precision   d(*), x(*)

*     dddiv  performs the diagonal scaling  x  =  x / d.

      integer            i, id, ix
      external           dscal
      intrinsic          abs
      parameter        ( one = 1.0d+0 )

      if (n .gt. 0) then
         if (incd .eq. 0  .and.  incx .ne. 0) then
            call dscal ( n, one/d(1), x, abs(incx) )
         else if (incd .eq. incx  .and.  incd .gt. 0) then
            do 10 id = 1, 1 + (n - 1)*incd, incd
               x(id) = x(id) / d(id)
   10       continue
         else
            if (incx .ge. 0) then
               ix = 1
            else
               ix = 1 - (n - 1)*incx
            end if
            if (incd .gt. 0) then
               do 20 id = 1, 1 + (n - 1)*incd, incd
                  x(ix) = x(ix) / d(id)
                  ix    = ix   + incx
   20          continue
            else
               id = 1 - (n - 1)*incd
               do 30  i = 1, n
                  x(ix) = x(ix) / d(id)
                  id    = id + incd
                  ix    = ix + incx
   30          continue
            end if
         end if
      end if

*     end of dddiv
      end
*UPTODATE F06DFFTEXT
      SUBROUTINE F06DFF( N, X, INCX, Y, INCY )
      ENTRY      ICOPY ( N, X, INCX, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
C     .. Array Arguments ..
      INTEGER            X( * ), Y( * )
C     ..
C
C  F06DFF performs the operation
C
C     y := x
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 10-February-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      INTEGER            I, IX, IY
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ( INCX.EQ.INCY ).AND.( INCY.GT.0 ) )THEN
            DO 10, IY = 1, 1 + ( N - 1 )*INCY, INCY
               Y( IY ) = X( IY )
   10       CONTINUE
         ELSE
            IF( INCX.GE.0 )THEN
               IX = 1
            ELSE
               IX = 1 - ( N - 1 )*INCX
            END IF
            IF( INCY.GT.0 )THEN
               DO 20, IY = 1, 1 + ( N - 1 )*INCY, INCY
                  Y( IY ) = X( IX )
                  IX      = IX      + INCX
   20          CONTINUE
            ELSE
               IY = 1 - ( N - 1 )*INCY
               DO 30, I = 1, N
                  Y( IY ) = X( IX )
                  IY      = IY      + INCY
                  IX      = IX      + INCX
   30          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06DFF. ( ICOPY )
C
      END
** END OF F06DFFTEXT
*UPTODATE F06FJFTEXT
      SUBROUTINE F06FJF( N, X, INCX, SCALE, SUMSQ )
      ENTRY      DSSQ  ( N, X, INCX, SCALE, SUMSQ )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION   SCALE, SUMSQ
      INTEGER            INCX, N
C     .. Array Arguments ..                   
      DOUBLE PRECISION   X( * )
C     ..
C
C  F06FJF returns the values scl and smsq such that
C
C     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
C
C  where x( i ) = X( 1 + ( i - 1 )*INCX ). The value of sumsq is assumed
C  to be at least unity and the value of smsq will then satisfy
C
C     1.0 .le. smsq .le. ( sumsq + n ) .
C
C  scale is assumed to be non-negative and scl returns the value
C
C     scl = max( scale, abs( x( i ) ) ) .
C
C  scale and sumsq must be supplied in SCALE and SUMSQ respectively.
C  scl and smsq are overwritten on SCALE and SUMSQ respectively.
C
C  The routine makes only one pass through the vector X.
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 22-October-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   ABSXI
      INTEGER            IX
C     .. Intrinsic Functions ..
      INTRINSIC          ABS
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
            IF( X( IX ).NE.ZERO )THEN
               ABSXI = ABS( X( IX ) )
               IF( SCALE.LT.ABSXI )THEN
                  SUMSQ = 1     + SUMSQ*( SCALE/ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SUMSQ = SUMSQ +       ( ABSXI/SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
      RETURN
C
C     End of F06FJF. ( DSSQ )
C
      END
** END OF F06FJFTEXT
*UPTODATE F06FKFTEXT
      DOUBLE PRECISION FUNCTION F06FKF( N, W, INCW, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      INTEGER                           INCW, INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION                  W( * ), X( * )
C     ..
C
C  F06FKF returns the weighted euclidean norm of a vector via the
C  function name, so that
C
C     F06FKF := sqrt( x'*W*x ),   where   W = diag( w ).
C
C  The elements of w are assumed to be non-negative.
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 25-June-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION      ONE         , ZERO
      PARAMETER           ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION      ABSYI, NORM, SCALE, SSQ
      INTEGER               I, IW, IX
C     .. External Functions ..
      DOUBLE PRECISION      F06BMF
      EXTERNAL              F06BMF
C     .. Intrinsic Functions ..
      INTRINSIC             ABS, SQRT
C     ..
C     .. Executable Statements ..
      IF( N.LT.1 )THEN
         NORM  = ZERO
      ELSE IF( N.EQ.1 )THEN
         NORM  = SQRT( W( 1 ) )*ABS( X( 1 ) )
      ELSE
         IF( INCW.GT.0 )THEN
            IW = 1
         ELSE
            IW = 1 - ( N - 1 )*INCW
         END IF
         IF( INCX.GT.0 )THEN
            IX = 1
         ELSE
            IX = 1 - ( N - 1 )*INCX
         END IF
         SCALE = ZERO
         SSQ   = ONE
         DO 10, I = 1, N
            IF( ( W( IW ).NE.ZERO ).AND.( X( IX ).NE.ZERO ) )THEN
               ABSYI = SQRT( W( IW ) )*ABS( X( IX ) )
               IF( SCALE.LT.ABSYI )THEN
                  SSQ   = 1     + SSQ*( SCALE/ABSYI )**2
                  SCALE = ABSYI
               ELSE
                  SSQ   = SSQ   +     ( ABSYI/SCALE )**2
               END IF
            END IF
            IW = IW + INCW
            IX = IX + INCX
   10    CONTINUE
         NORM = F06BMF( SCALE, SSQ )
      END IF
C
      F06FKF = NORM
      RETURN
C
C     End of F06FKF.
C
      END
** END OF F06FKFTEXT
*UPTODATE F06FLFTEXT
      SUBROUTINE F06FLF( N, X, INCX, XMAX, XMIN )
      ENTRY      DCOND ( N, X, INCX, XMAX, XMIN )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION   XMAX, XMIN
      INTEGER            INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
C     ..
C
C  F06FLF returns the values xmax and xmin given by
C
C     xmax = max( abs( x( i ) ) ),   xmin = min( abs( x( i ) ) ).
C             i                              i
C
C  If n is less than unity then xmax and xmin are returned as zero.
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 27-February-1986.
C     Sven Hammarling, Nag Central Office.
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      INTEGER            IX
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
C     ..
C     .. Executable Statements ..
      IF( N.LT.1 )THEN
         XMAX = ZERO
         XMIN = ZERO
      ELSE
         XMAX = ABS( X( 1 ) )
         XMIN = XMAX
         DO 10 IX = 1 + INCX, 1 + ( N - 1 )*INCX, INCX
            XMAX = MAX( XMAX, ABS( X( IX ) ) )
            XMIN = MIN( XMIN, ABS( X( IX ) ) )
   10    CONTINUE
      END IF
C
      RETURN
C
C     End of F06FLF. ( DCOND )
C
      END
** END OF F06FLFTEXT
*UPTODATE F06KLFTEXT
      INTEGER FUNCTION F06KLF( N, X, INCX, TOL )
      INTEGER          IDRANK                    
      ENTRY            IDRANK( N, X, INCX, TOL )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     Modified by PEG 9/25/88.
C     .. Scalar Arguments ..
      DOUBLE PRECISION         TOL
      INTEGER                  INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION         X( * )
C     ..
C
C  F06KLF finds the first element of the n element vector x for which
C
C     abs( x( k ) ).le.tol*max( abs( x( 1 ) ), ..., abs( x( k - 1 ) ) )
C
C  and returns the value ( k - 1 ) in the function name F06KLF. If no
C  such k exists then F06KLF is returned as n.
C
C  If tol is supplied as less than zero then the value epsmch, where
C  epsmch is the relative machine precision, is used in place of tol.
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 27-February-1986.
C     Sven Hammarling, Nag Central Office.
C                        
C     .. Parameters ..
      DOUBLE PRECISION         ZERO
      PARAMETER              ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION         TL, XMAX
      INTEGER                  IX, K
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
C     .. Intrinsic Functions ..
      INTRINSIC                ABS, MAX
C     ..
C     .. Executable Statements ..
      K = 0
      IF( N.GE.1 )THEN
         IX = 1
         IF( TOL.LT.ZERO )THEN
            TL = WMACH(3)
         ELSE
            TL = TOL
         END IF
         XMAX = ABS( X( IX ) )
C
C+       WHILE( K.LT.N )LOOP
   10    IF   ( K.LT.N )THEN
            IF( ABS( X( IX ) ).LE.TL*XMAX )
     $         GO TO 20
            XMAX = MAX( XMAX, ABS( X( IX ) ) )
            K    = K  + 1
            IX   = IX + INCX
            GO TO 10
         END IF
C+       END WHILE
C
      END IF
C
   20 F06KLF = K
      RETURN
C
C     End of F06KLF. ( IDRANK )
C
      END
** END OF F06KLFTEXT
*UPTODATE F06FQFTEXT
      SUBROUTINE F06FQF( PIVOT, DIRECT, N, ALPHA, X, INCX, C, S )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, N
      CHARACTER*1        DIRECT, PIVOT
C     .. Array Arguments ..
      DOUBLE PRECISION   C( * ), S( * ), X( * )
C     ..
C
C  F06FQF generates the parameters of an orthogonal matrix P such that
C
C     when   PIVOT = 'F' or 'f'   and   DIRECT = 'F' or 'f'
C     or     PIVOT = 'V' or 'v'   and   DIRECT = 'B' or 'b'
C
C        P*( alpha ) = ( beta ),
C          (   x   )   (   0  )
C
C     when   PIVOT = 'F' or 'f'   and   DIRECT = 'B' or 'b'
C     or     PIVOT = 'V' or 'v'   and   DIRECT = 'F' or 'f'
C
C        P*(   x   ) = (   0  ),
C          ( alpha ) = ( beta )
C
C  where alpha is a scalar and x is an n element vector.
C
C  When  PIVOT = 'F' or 'f'  ( fixed pivot )
C  and  DIRECT = 'F' or 'f'  ( forward sequence ) then
C        
C     P is given as the sequence of plane rotation matrices
C
C        P = P( n )*P( n - 1 )*...*P( 1 )
C
C     where P( k ) is a plane rotation matrix for the ( 1, k + 1 ) plane
C     designed to annihilate the kth element of x.
C
C  When  PIVOT = 'V' or 'v'  ( variable pivot )
C  and  DIRECT = 'B' or 'b'  ( backward sequence ) then
C
C     P is given as the sequence of plane rotation matrices
C
C        P = P( 1 )*P( 2 )*...*P( n )
C
C     where P( k ) is a plane rotation matrix for the ( k, k + 1 ) plane
C     designed to annihilate the kth element of x.
C
C  When  PIVOT = 'F' or 'f'  ( fixed pivot )
C  and  DIRECT = 'B' or 'b'  ( backward sequence ) then
C
C     P is given as the sequence of plane rotation matrices
C
C        P = P( 1 )*P( 2 )*...*P( n )
C
C     where P( k ) is a plane rotation matrix for the ( k, n + 1 ) plane
C     designed to annihilate the kth element of x.
C
C  When  PIVOT = 'V' or 'v'  ( variable pivot )
C  and  DIRECT = 'F' or 'f'  ( forward sequence ) then
C
C     P is given as the sequence of plane rotation matrices
C
C        P = P( n )*P( n - 1 )*...*P( 1 )
C
C     where P( k ) is a plane rotation matrix for the ( k, k + 1 ) plane
C     designed to annihilate the kth element of x.
C
C  The routine returns the cosine, c( k ), and sine, s( k ) that define
C  the matrix P( k ), such that the two by two rotation part of P( k ),
C  R( k ), has the form
C
C     R( k ) = (  c( k )  s( k ) ).
C              ( -s( k )  c( k ) )
C
C  On entry, ALPHA must contain  the scalar alpha and on exit, ALPHA is
C  overwritten by beta. The cosines and sines are returned in the arrays
C  C and S and the vector x is overwritten by the tangents of the plane
C  rotations ( t( k ) = s( k )/c( k ) ).
C
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 19-April-1985.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      INTEGER            I, IX
C     .. External Subroutines ..
      EXTERNAL           F06BAF
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
            IX = 1 + ( N - 1 )*INCX
            IF( ( PIVOT.EQ.'V' ).OR.( PIVOT.EQ.'v' ) )THEN
               DO 10, I = N, 2, -1
                  CALL F06BAF( X( IX - INCX ), X( IX ), C( I ), S( I ) )
                  IX = IX - INCX
   10          CONTINUE
               CALL F06BAF( ALPHA, X( IX ), C( 1 ), S( 1 ) )
            ELSE IF( ( PIVOT.EQ.'F' ).OR.( PIVOT.EQ.'f' ) )THEN
C
C              Here we choose c and s so that
C
C                 ( alpha ) := (  c  s )*( alpha  )
C                 (   0   )    ( -s  c ) ( x( i ) )
C
C              which is equivalent to
C
C                 (   0   ) := ( c  -s )*( x( i ) )
C                 ( alpha )    ( s   c ) ( alpha  )
C
C              and so we need to return  s( i ) = -s  in order to make
C              R( i ) look like
C
C                 R( i ) = (  c( i )  s( i ) ).
C                          ( -s( i )  c( i ) )
C
               DO 20, I = N, 1, -1
                  CALL F06BAF( ALPHA, X( IX ), C( I ), S( I ) )
                  S( I )  = -S( I )
                  X( IX ) = -X( IX )
                  IX      =  IX      - INCX
   20          CONTINUE
            END IF
         ELSE IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
            IX = 1
            IF( ( PIVOT.EQ.'V' ).OR.( PIVOT.EQ.'v' ) )THEN
C
C              Here we choose c and s so that
C
C                 ( x( i + 1 ) ) := (  c  s )*( x( i + 1 ) )
C                 (    0       )    ( -s  c ) ( x( i )     )
C
C              which is equivalent to
C
C                 (    0       ) := ( c  -s )*( x( i )     )
C                 ( x( i + 1 ) )    ( s   c ) ( x( i + 1 ) )
C
C              and so we need to return  s( i ) = -s  in order to make
C              R( i ) look like
C
C                 R( i ) = (  c( i )  s( i ) ).
C                          ( -s( i )  c( i ) )
C
               DO 30, I = 1, N - 1
                  CALL F06BAF( X( IX + INCX ), X( IX ), C( I ), S( I ) )
                  S( I )  = -S( I )
                  X( IX ) = -X( IX )
                  IX      =  IX      + INCX
   30          CONTINUE
               CALL F06BAF( ALPHA, X( IX ), C( N ), S( N ) )
               S( N )  = -S( N )
               X( IX ) = -X( IX )
            ELSE IF( ( PIVOT.EQ.'F' ).OR.( PIVOT.EQ.'f' ) )THEN
               DO 40, I = 1, N
                  CALL F06BAF( ALPHA, X( IX ), C( I ), S( I ) )
                  IX = IX + INCX
   40          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06FQF. ( SSROTG )
C
      END
** END OF F06FQFTEXT
*UPTODATE F06FRFTEXT
      SUBROUTINE F06FRF( N, ALPHA, X, INCX, TOL, ZETA )
      ENTRY      DGRFG ( N, ALPHA, X, INCX, TOL, ZETA )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     Modified by PEG 9/25/88.
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, TOL, ZETA
      INTEGER            INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
C     ..
C
C  F06FRF generates details of a generalized Householder reflection such
C  that
C
C     P*( alpha ) = ( beta ),   P'*P = I.
C       (   x   )   (   0  )
C
C  P is given in the form
C
C     P = I - ( zeta )*( zeta  z' ),
C             (   z  )
C
C  where z is an n element vector and zeta is a scalar that satisfies
C
C     1.0 .le. zeta .le. sqrt( 2.0 ).
C
C  zeta is returned in ZETA unless x is such that
C
C     max( abs( x( i ) ) ) .le. max( eps*abs( alpha ), tol )
C
C  where eps is the relative machine precision and tol is the user
C  supplied value TOL, in which case ZETA is returned as 0.0 and P can
C  be taken to be the unit matrix.
C
C  beta is overwritten on alpha and z is overwritten on x.
C  the routine may be called with  n = 0  and advantage is taken of the
C  case where  n = 1.
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 30-August-1984.
C     Sven Hammarling, Nag Central Office.
C     This version dated 28-September-1984.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   BETA, EPS, SCALE, SSQ
      LOGICAL            FIRST
      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/
C     .. External Subroutines ..
      EXTERNAL           F06FJF, DSCAL
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SIGN, SQRT
C     .. Save statement ..
      SAVE               EPS, FIRST
C     .. Data statements ..
      DATA               FIRST/ .TRUE. /
C     ..
C     .. Executable Statements ..
      IF( N.LT.1 )THEN
         ZETA = ZERO
      ELSE IF( ( N.EQ.1 ).AND.( X( 1 ).EQ.ZERO ) )THEN
         ZETA = ZERO
      ELSE
C
         IF( FIRST )THEN
            FIRST = .FALSE.
            EPS   =  WMACH(3)
         END IF
C
C        Treat case where P is a 2 by 2 matrix specially.
C
         IF( N.EQ.1 )THEN
C
C           Deal with cases where  ALPHA = zero  and
C           abs( X( 1 ) ) .le. max( EPS*abs( ALPHA ), TOL )  first.
C
            IF( ALPHA.EQ.ZERO )THEN
               ZETA   =  ONE
               ALPHA  =  ABS ( X( 1 ) )
               X( 1 ) = -SIGN( ONE, X( 1 ) )
            ELSE IF( ABS( X( 1 ) ).LE.MAX( EPS*ABS( ALPHA ), TOL ) )THEN
               ZETA   =  ZERO
            ELSE
               IF( ABS( ALPHA ).GE.ABS( X( 1 ) ) )THEN
                  BETA = ABS( ALPHA ) *SQRT( 1 + ( X( 1 )/ALPHA )**2 )
               ELSE
                  BETA = ABS( X( 1 ) )*SQRT( 1 + ( ALPHA/X( 1 ) )**2 )
               END IF
               ZETA = SQRT( ( ABS( ALPHA ) + BETA )/BETA )
               IF( ALPHA.GE.ZERO )
     $            BETA = -BETA
               X( 1 ) = -X( 1 )/( ZETA*BETA )
               ALPHA  = BETA
            END IF
         ELSE
C
C           Now P is larger than 2 by 2.
C
            SSQ   = ONE
            SCALE = ZERO
            CALL F06FJF( N, X, INCX, SCALE, SSQ )
C
C           Treat cases where  SCALE = zero,
C           SCALE .le. max( EPS*abs( ALPHA ), TOL )  and
C           ALPHA = zero  specially.
C           Note that  SCALE = max( abs( X( i ) ) ).
C
            IF( ( SCALE.EQ.ZERO ).OR.
     $          ( SCALE.LE.MAX( EPS*ABS( ALPHA ), TOL ) ) )THEN
               ZETA  = ZERO
            ELSE IF( ALPHA.EQ.ZERO )THEN
               ZETA  = ONE
               ALPHA = SCALE*SQRT( SSQ )
               CALL DSCAL( N, -1/ALPHA, X, INCX )
            ELSE
               IF( SCALE.LT.ABS( ALPHA ) )THEN
                  BETA = ABS( ALPHA )*SQRT( 1 + SSQ*( SCALE/ALPHA )**2 )
               ELSE
                  BETA = SCALE       *SQRT( SSQ +   ( ALPHA/SCALE )**2 )
               END IF
               ZETA = SQRT( ( BETA + ABS( ALPHA ) )/BETA )
               IF( ALPHA.GT.ZERO )
     $            BETA = -BETA
               CALL DSCAL( N, -1/( ZETA*BETA ), X, INCX )
               ALPHA = BETA
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06FRF. ( DGRFG )
C
      END
** END OF F06FRFTEXT
*UPTODATE F06PAFTEXT
      SUBROUTINE F06PAF( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     .. Entry Points ..
      ENTRY      DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
C     ..
C
C  Purpose
C  =======
C
C  DGEMV  performs one of the matrix-vector operations
C
C     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
C
C  where alpha and beta are scalars, x and y are vectors and A is an
C  m by n matrix.
C
C  Parameters
C  ==========
C
C  TRANS  - CHARACTER*1.
C           On entry, TRANS specifies the operation to be performed as
C           follows:
C
C              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
C
C              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
C
C              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C           On entry, M specifies the number of rows of the matrix A.
C           M must be at least zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the number of columns of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - DOUBLE PRECISION.
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry, the leading m by n part of the array A must
C           contain the matrix of coefficients.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, m ).
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of DIMENSION at least
C           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
C           and at least
C           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
C           Before entry, the incremented array X must contain the
C           vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  BETA   - DOUBLE PRECISION.
C           On entry, BETA specifies the scalar beta. When BETA is
C           supplied as zero then Y need not be set on input.
C           Unchanged on exit.
C
C  Y      - DOUBLE PRECISION array of DIMENSION at least
C           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
C           and at least
C           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
C           Before entry with BETA non-zero, the incremented array Y
C           must contain the vector y. On exit, Y is overwritten by the
C           updated vector y.
C
C  INCY   - INTEGER.
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.(TRANS.EQ.'N' .OR. TRANS.EQ.'n').AND.
     $         .NOT.(TRANS.EQ.'T' .OR. TRANS.EQ.'t').AND.
     $         .NOT.(TRANS.EQ.'C' .OR. TRANS.EQ.'c')      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06PAF/DGEMV ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
C
C     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
C     up the start points in  X  and  Y.
C
      IF( (TRANS.EQ.'N' .OR. TRANS.EQ.'n') )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through A.
C
C     First form  y := beta*y.
C
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( (TRANS.EQ.'N' .OR. TRANS.EQ.'n') )THEN
C
C        Form  y := alpha*A*x + y.
C
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
C
C        Form  y := alpha*A'*x + y.
C
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO 100, J = 1, N
               TEMP = ZERO
               DO 90, I = 1, M
                  TEMP = TEMP + A( I, J )*X( I )
   90          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  100       CONTINUE
         ELSE
            DO 120, J = 1, N
               TEMP = ZERO
               IX   = KX
               DO 110, I = 1, M
                  TEMP = TEMP + A( I, J )*X( IX )
                  IX   = IX   + INCX
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06PAF (DGEMV ).
C
      END
** END OF F06PAFTEXT
*UPTODATE F06PMFTEXT
      SUBROUTINE F06PMF( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     .. Entry Points ..
      ENTRY      DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, LDA, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
C     ..
C
C  Purpose
C  =======
C
C  DGER   performs the rank 1 operation
C
C     A := alpha*x*y' + A,
C
C  where alpha is a scalar, x is an m element vector, y is an n element
C  vector and A is an m by n matrix.
C
C  Parameters
C  ==========
C
C  M      - INTEGER.
C           On entry, M specifies the number of rows of the matrix A.
C           M must be at least zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the number of columns of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - DOUBLE PRECISION.
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( m - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the m
C           element vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  Y      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCY ) ).
C           Before entry, the incremented array Y must contain the n
C           element vector y.
C           Unchanged on exit.
C
C  INCY   - INTEGER.
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry, the leading m by n part of the array A must
C           contain the matrix of coefficients. On exit, A is
C           overwritten by the updated matrix.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, m ).
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JY, KX
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06PMF/DGER  ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through A.
C
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
C
      RETURN
C
C     End of F06PMF (DGER  ).
C
      END
** END OF F06PMFTEXT
*UPTODATE F06PCFTEXT
      SUBROUTINE F06PCF( UPLO, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     .. Entry Points ..
      ENTRY      DSYMV ( UPLO, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, N
      CHARACTER*1        UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
C     ..
C
C  Purpose
C  =======
C
C  DSYMV  performs the matrix-vector  operation
C
C     y := alpha*A*x + beta*y,
C
C  where alpha and beta are scalars, x and y are n element vectors and
C  A is an n by n symmetric matrix.
C
C  Parameters
C  ==========
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the upper or lower
C           triangular part of the array A is to be referenced as
C           follows:
C
C              UPLO = 'U' or 'u'   Only the upper triangular part of A
C                                  is to be referenced.
C
C              UPLO = 'L' or 'l'   Only the lower triangular part of A
C                                  is to be referenced.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the order of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - DOUBLE PRECISION.
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry with  UPLO = 'U' or 'u', the leading n by n
C           upper triangular part of the array A must contain the upper
C           triangular part of the symmetric matrix and the strictly
C           lower triangular part of A is not referenced.
C           Before entry with UPLO = 'L' or 'l', the leading n by n
C           lower triangular part of the array A must contain the lower
C           triangular part of the symmetric matrix and the strictly
C           upper triangular part of A is not referenced.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, n ).
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the n
C           element vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  BETA   - DOUBLE PRECISION.
C           On entry, BETA specifies the scalar beta. When BETA is
C           supplied as zero then Y need not be set on input.
C           Unchanged on exit.
C
C  Y      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCY ) ).
C           Before entry, the incremented array Y must contain the n
C           element vector y. On exit, Y is overwritten by the updated
C           vector y.
C
C  INCY   - INTEGER.
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP1, TEMP2
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.(UPLO.EQ.'U' .OR. UPLO.EQ.'u').AND.
     $         .NOT.(UPLO.EQ.'L' .OR. UPLO.EQ.'l')      )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 5
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 7
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 10
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06PCF/DSYMV ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( N.EQ.0 ).OR.( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
C
C     Set up the start points in  X  and  Y.
C
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( N - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( N - 1 )*INCY
      END IF
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through the triangular part
C     of A.
C
C     First form  y := beta*y.
C
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, N
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, N
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, N
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, N
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
C
C        Form  y  when A is stored in upper triangle.
C
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 60, J = 1, N
               TEMP1 = ALPHA*X( J )
               TEMP2 = ZERO
               DO 50, I = 1, J - 1
                  Y( I ) = Y( I ) + TEMP1*A( I, J )
                  TEMP2  = TEMP2  + A( I, J )*X( I )
   50          CONTINUE
               Y( J ) = Y( J ) + TEMP1*A( J, J ) + ALPHA*TEMP2
   60       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 80, J = 1, N
               TEMP1 = ALPHA*X( JX )
               TEMP2 = ZERO
               IX    = KX
               IY    = KY
               DO 70, I = 1, J - 1
                  Y( IY ) = Y( IY ) + TEMP1*A( I, J )
                  TEMP2   = TEMP2   + A( I, J )*X( IX )
                  IX      = IX      + INCX
                  IY      = IY      + INCY
   70          CONTINUE
               Y( JY ) = Y( JY ) + TEMP1*A( J, J ) + ALPHA*TEMP2
               JX      = JX      + INCX
               JY      = JY      + INCY
   80       CONTINUE
         END IF
      ELSE
C
C        Form  y  when A is stored in lower triangle.
C
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 100, J = 1, N
               TEMP1  = ALPHA*X( J )
               TEMP2  = ZERO
               Y( J ) = Y( J )       + TEMP1*A( J, J )
               DO 90, I = J + 1, N
                  Y( I ) = Y( I ) + TEMP1*A( I, J )
                  TEMP2  = TEMP2  + A( I, J )*X( I )
   90          CONTINUE
               Y( J ) = Y( J ) + ALPHA*TEMP2
  100       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 120, J = 1, N
               TEMP1   = ALPHA*X( JX )
               TEMP2   = ZERO
               Y( JY ) = Y( JY )       + TEMP1*A( J, J )
               IX      = JX
               IY      = JY
               DO 110, I = J + 1, N
                  IX      = IX      + INCX
                  IY      = IY      + INCY
                  Y( IY ) = Y( IY ) + TEMP1*A( I, J )
                  TEMP2   = TEMP2   + A( I, J )*X( IX )
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP2
               JX      = JX      + INCX
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06PCF (DSYMV ).
C
      END
** END OF F06PCFTEXT
*UPTODATE F06PPFTEXT
      SUBROUTINE F06PPF( UPLO, N, ALPHA, X, INCX, A, LDA )
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     .. Entry Points ..
      ENTRY      DSYR  ( UPLO, N, ALPHA, X, INCX, A, LDA )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, LDA, N
      CHARACTER*1        UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * )
C     ..
C
C  Purpose
C  =======
C
C  DSYR   performs the symmetric rank 1 operation
C
C     A := alpha*x*x' + A,
C
C  where alpha is a real scalar, x is an n element vector and A is an
C  n by n symmetric matrix.
C
C  Parameters
C  ==========
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the upper or lower
C           triangular part of the array A is to be referenced as
C           follows:
C
C              UPLO = 'U' or 'u'   Only the upper triangular part of A
C                                  is to be referenced.
C
C              UPLO = 'L' or 'l'   Only the lower triangular part of A
C                                  is to be referenced.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the order of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - DOUBLE PRECISION.
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the n
C           element vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry with  UPLO = 'U' or 'u', the leading n by n
C           upper triangular part of the array A must contain the upper
C           triangular part of the symmetric matrix and the strictly
C           lower triangular part of A is not referenced. On exit, the
C           upper triangular part of the array A is overwritten by the
C           upper triangular part of the updated matrix.
C           Before entry with UPLO = 'L' or 'l', the leading n by n
C           lower triangular part of the array A must contain the lower
C           triangular part of the symmetric matrix and the strictly
C           upper triangular part of A is not referenced. On exit, the
C           lower triangular part of the array A is overwritten by the
C           lower triangular part of the updated matrix.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, n ).
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JX, KX
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.(UPLO.EQ.'U' .OR. UPLO.EQ.'u').AND.
     $         .NOT.(UPLO.EQ.'L' .OR. UPLO.EQ.'l')      )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 7
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06PPF/DSYR  ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
C
C     Set the start point in X if the increment is not unity.
C
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through the triangular part
C     of A.
C
      IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
C
C        Form  A  when A is stored in upper triangle.
C
         IF( INCX.EQ.1 )THEN
            DO 20, J = 1, N
               IF( X( J ).NE.ZERO )THEN
                  TEMP = ALPHA*X( J )
                  DO 10, I = 1, J
                     A( I, J ) = A( I, J ) + X( I )*TEMP
   10             CONTINUE
               END IF
   20       CONTINUE
         ELSE
            JX = KX
            DO 40, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IX   = KX
                  DO 30, I = 1, J
                     A( I, J ) = A( I, J ) + X( IX )*TEMP
                     IX        = IX        + INCX
   30             CONTINUE
               END IF
               JX = JX + INCX
   40       CONTINUE
         END IF
      ELSE
C
C        Form  A  when A is stored in lower triangle.
C
         IF( INCX.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( J ).NE.ZERO )THEN
                  TEMP = ALPHA*X( J )
                  DO 50, I = J, N
                     A( I, J ) = A( I, J ) + X( I )*TEMP
   50             CONTINUE
               END IF
   60       CONTINUE
         ELSE
            JX = KX
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IX   = JX
                  DO 70, I = J, N
                     A( I, J ) = A( I, J ) + X( IX )*TEMP
                     IX        = IX        + INCX
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06PPF (DSYR  ).
C
      END
** END OF F06PPFTEXT
*UPTODATE F06PFFTEXT
      SUBROUTINE F06PFF( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     .. Entry Points ..
      ENTRY      DTRMV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
C     .. Scalar Arguments ..
      INTEGER            INCX, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * )
C     ..
C
C  Purpose
C  =======
C
C  DTRMV  performs one of the matrix-vector operations
C
C     x := A*x,   or   x := A'*x,
C
C  where x is n element vector and A is an n by n unit, or non-unit,
C  upper or lower triangular matrix.
C
C  Parameters
C  ==========
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the matrix is an upper or
C           lower triangular matrix as follows:
C
C              UPLO = 'U' or 'u'   A is an upper triangular matrix.
C
C              UPLO = 'L' or 'l'   A is a lower triangular matrix.
C
C           Unchanged on exit.
C
C  TRANS  - CHARACTER*1.
C           On entry, TRANS specifies the operation to be performed as
C           follows:
C
C              TRANS = 'N' or 'n'   x := A*x.
C
C              TRANS = 'T' or 't'   x := A'*x.
C
C              TRANS = 'C' or 'c'   x := A'*x.
C
C           Unchanged on exit.
C
C  DIAG   - CHARACTER*1.
C           On entry, DIAG specifies whether or not A is unit
C           triangular as follows:
C
C              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
C
C              DIAG = 'N' or 'n'   A is not assumed to be unit
C                                  triangular.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the order of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry with  UPLO = 'U' or 'u', the leading n by n
C           upper triangular part of the array A must contain the upper
C           triangular matrix and the strictly lower triangular part of
C           A is not referenced.
C           Before entry with UPLO = 'L' or 'l', the leading n by n
C           lower triangular part of the array A must contain the lower
C           triangular matrix and the strictly upper triangular part of
C           A is not referenced.
C           Note that when  DIAG = 'U' or 'u', the diagonal elements of
C           A are not referenced either, but are assumed to be unity.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, n ).
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the n
C           element vector x. On exit, X is overwritten with the
C           tranformed vector x.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JX, KX
      LOGICAL            NOUNIT
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.(UPLO .EQ.'U' .OR. UPLO .EQ.'u').AND.
     $         .NOT.(UPLO .EQ.'L' .OR. UPLO .EQ.'l')      )THEN
         INFO = 1
      ELSE IF( .NOT.(TRANS.EQ.'N' .OR. TRANS.EQ.'n').AND.
     $         .NOT.(TRANS.EQ.'T' .OR. TRANS.EQ.'t').AND.
     $         .NOT.(TRANS.EQ.'C' .OR. TRANS.EQ.'c')      )THEN
         INFO = 2
      ELSE IF( .NOT.(DIAG .EQ.'U' .OR. DIAG .EQ.'u').AND.
     $         .NOT.(DIAG .EQ.'N' .OR. DIAG .EQ.'n')      )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06PFF/DTRMV ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 )
     $   RETURN
C
      NOUNIT = (DIAG.EQ.'N' .OR. DIAG.EQ.'n')
C
C     Set up the start point in X if the increment is not unity. This
C     will be  ( N - 1 )*INCX  too small for descending loops.
C
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through A.
C
      IF( (TRANS.EQ.'N' .OR. TRANS.EQ.'n') )THEN
C
C        Form  x := A*x.
C
         IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
            IF( INCX.EQ.1 )THEN
               DO 20, J = 1, N
                  IF( X( J ).NE.ZERO )THEN
                     TEMP = X( J )
                     DO 10, I = 1, J - 1
                        X( I ) = X( I ) + TEMP*A( I, J )
   10                CONTINUE
                     IF( NOUNIT )
     $                  X( J ) = X( J )*A( J, J )
                  END IF
   20          CONTINUE
            ELSE
               JX = KX
               DO 40, J = 1, N
                  IF( X( JX ).NE.ZERO )THEN
                     TEMP = X( JX )
                     IX   = KX
                     DO 30, I = 1, J - 1
                        X( IX ) = X( IX ) + TEMP*A( I, J )
                        IX      = IX      + INCX
   30                CONTINUE
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )*A( J, J )
                  END IF
                  JX = JX + INCX
   40          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 60, J = N, 1, -1
                  IF( X( J ).NE.ZERO )THEN
                     TEMP = X( J )
                     DO 50, I = N, J + 1, -1
                        X( I ) = X( I ) + TEMP*A( I, J )
   50                CONTINUE
                     IF( NOUNIT )
     $                  X( J ) = X( J )*A( J, J )
                  END IF
   60          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 80, J = N, 1, -1
                  IF( X( JX ).NE.ZERO )THEN
                     TEMP = X( JX )
                     IX   = KX
                     DO 70, I = N, J + 1, -1
                        X( IX ) = X( IX ) + TEMP*A( I, J )
                        IX      = IX      - INCX
   70                CONTINUE
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )*A( J, J )
                  END IF
                  JX = JX - INCX
   80          CONTINUE
            END IF
         END IF
      ELSE
C
C        Form  x := A'*x.
C
         IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
            IF( INCX.EQ.1 )THEN
               DO 100, J = N, 1, -1
                  TEMP = X( J )
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 90, I = J - 1, 1, -1
                     TEMP = TEMP + A( I, J )*X( I )
   90             CONTINUE
                  X( J ) = TEMP
  100          CONTINUE
            ELSE
               JX = KX + ( N - 1 )*INCX
               DO 120, J = N, 1, -1
                  TEMP = X( JX )
                  IX   = JX
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 110, I = J - 1, 1, -1
                     IX   = IX   - INCX
                     TEMP = TEMP + A( I, J )*X( IX )
  110             CONTINUE
                  X( JX ) = TEMP
                  JX      = JX   - INCX
  120          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 140, J = 1, N
                  TEMP = X( J )
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 130, I = J + 1, N
                     TEMP = TEMP + A( I, J )*X( I )
  130             CONTINUE
                  X( J ) = TEMP
  140          CONTINUE
            ELSE
               JX = KX
               DO 160, J = 1, N
                  TEMP = X( JX )
                  IX   = JX
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 150, I = J + 1, N
                     IX   = IX   + INCX
                     TEMP = TEMP + A( I, J )*X( IX )
  150             CONTINUE
                  X( JX ) = TEMP
                  JX      = JX   + INCX
  160          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06PFF (DTRMV ).
C
      END
** END OF F06PFFTEXT
*UPTODATE F06PJFTEXT
      SUBROUTINE F06PJF( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     .. Entry Points ..
      ENTRY      DTRSV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
C     .. Scalar Arguments ..
      INTEGER            INCX, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * )
C     ..
C
C  Purpose
C  =======
C
C  DTRSV  solves one of the systems of equations
C
C     A*x = b,   or   A'*x = b,
C
C  where b and x are n element vectors and A is an n by n unit, or
C  non-unit, upper or lower triangular matrix.
C
C  No test for singularity or near-singularity is included in this
C  routine. Such tests must be performed before calling this routine.
C
C  Parameters
C  ==========
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the matrix is an upper or
C           lower triangular matrix as follows:
C
C              UPLO = 'U' or 'u'   A is an upper triangular matrix.
C
C              UPLO = 'L' or 'l'   A is a lower triangular matrix.
C
C           Unchanged on exit.
C
C  TRANS  - CHARACTER*1.
C           On entry, TRANS specifies the equations to be solved as
C           follows:
C
C              TRANS = 'N' or 'n'   A*x = b.
C
C              TRANS = 'T' or 't'   A'*x = b.
C
C              TRANS = 'C' or 'c'   A'*x = b.
C
C           Unchanged on exit.
C
C  DIAG   - CHARACTER*1.
C           On entry, DIAG specifies whether or not A is unit
C           triangular as follows:
C
C              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
C
C              DIAG = 'N' or 'n'   A is not assumed to be unit
C                                  triangular.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the order of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry with  UPLO = 'U' or 'u', the leading n by n
C           upper triangular part of the array A must contain the upper
C           triangular matrix and the strictly lower triangular part of
C           A is not referenced.
C           Before entry with UPLO = 'L' or 'l', the leading n by n
C           lower triangular part of the array A must contain the lower
C           triangular matrix and the strictly upper triangular part of
C           A is not referenced.
C           Note that when  DIAG = 'U' or 'u', the diagonal elements of
C           A are not referenced either, but are assumed to be unity.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, n ).
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the n
C           element right-hand side vector b. On exit, X is overwritten
C           with the solution vector x.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JX, KX
      LOGICAL            NOUNIT
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.(UPLO .EQ.'U' .OR. UPLO .EQ.'u').AND.
     $         .NOT.(UPLO .EQ.'L' .OR. UPLO .EQ.'l')      )THEN
         INFO = 1
      ELSE IF( .NOT.(TRANS.EQ.'N' .OR. TRANS.EQ.'n').AND.
     $         .NOT.(TRANS.EQ.'T' .OR. TRANS.EQ.'t').AND.
     $         .NOT.(TRANS.EQ.'C' .OR. TRANS.EQ.'c')      )THEN
         INFO = 2
      ELSE IF( .NOT.(DIAG .EQ.'U' .OR. DIAG .EQ.'u').AND.
     $         .NOT.(DIAG .EQ.'N' .OR. DIAG .EQ.'n')      )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06PJF/DTRSV ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 )
     $   RETURN
C
      NOUNIT = (DIAG.EQ.'N' .OR. DIAG.EQ.'n')
C
C     Set up the start point in X if the increment is not unity. This
C     will be  ( N - 1 )*INCX  too small for descending loops.
C
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through A.
C
      IF( (TRANS.EQ.'N' .OR. TRANS.EQ.'n') )THEN
C
C        Form  x := inv( A )*x.
C
         IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
            IF( INCX.EQ.1 )THEN
               DO 20, J = N, 1, -1
                  IF( X( J ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( J ) = X( J )/A( J, J )
                     TEMP = X( J )
                     DO 10, I = J - 1, 1, -1
                        X( I ) = X( I ) - TEMP*A( I, J )
   10                CONTINUE
                  END IF
   20          CONTINUE
            ELSE
               JX = KX + ( N - 1 )*INCX
               DO 40, J = N, 1, -1
                  IF( X( JX ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )/A( J, J )
                     TEMP = X( JX )
                     IX   = JX
                     DO 30, I = J - 1, 1, -1
                        IX      = IX      - INCX
                        X( IX ) = X( IX ) - TEMP*A( I, J )
   30                CONTINUE
                  END IF
                  JX = JX - INCX
   40          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 60, J = 1, N
                  IF( X( J ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( J ) = X( J )/A( J, J )
                     TEMP = X( J )
                     DO 50, I = J + 1, N
                        X( I ) = X( I ) - TEMP*A( I, J )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE
               JX = KX
               DO 80, J = 1, N
                  IF( X( JX ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )/A( J, J )
                     TEMP = X( JX )
                     IX   = JX
                     DO 70, I = J + 1, N
                        IX      = IX      + INCX
                        X( IX ) = X( IX ) - TEMP*A( I, J )
   70                CONTINUE
                  END IF
                  JX = JX + INCX
   80          CONTINUE
            END IF
         END IF
      ELSE
C
C        Form  x := inv( A' )*x.
C
         IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
            IF( INCX.EQ.1 )THEN
               DO 100, J = 1, N
                  TEMP = X( J )
                  DO 90, I = 1, J - 1
                     TEMP = TEMP - A( I, J )*X( I )
   90             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( J, J )
                  X( J ) = TEMP
  100          CONTINUE
            ELSE
               JX = KX
               DO 120, J = 1, N
                  TEMP = X( JX )
                  IX   = KX
                  DO 110, I = 1, J - 1
                     TEMP = TEMP - A( I, J )*X( IX )
                     IX   = IX   + INCX
  110             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( J, J )
                  X( JX ) = TEMP
                  JX      = JX   + INCX
  120          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 140, J = N, 1, -1
                  TEMP = X( J )
                  DO 130, I = N, J + 1, -1
                     TEMP = TEMP - A( I, J )*X( I )
  130             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( J, J )
                  X( J ) = TEMP
  140          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 160, J = N, 1, -1
                  TEMP = X( JX )
                  IX   = KX
                  DO 150, I = N, J + 1, -1
                     TEMP = TEMP - A( I, J )*X( IX )
                     IX   = IX   - INCX
  150             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( J, J )
                  X( JX ) = TEMP
                  JX      = JX   - INCX
  160          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06PJF (DTRSV ).
C
      END
** END OF F06PJFTEXT
*UPTODATE F06QFFTEXT
      SUBROUTINE F06QFF( MATRIX, M, N, A, LDA, B, LDB )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      CHARACTER*1        MATRIX
      INTEGER            M, N, LDA, LDB
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
C     ..
C
C  F06QFF  copies  the  m by n  matrix  A  into  the  m by n  matrix  B.
C
C  If   MATRIX = 'G' or 'g'   then  A  and  B  are  regarded as  general
C                             matrices,
C  if   MATRIX = 'U' or 'u'   then  A  and  B  are  regarded  as   upper
C                             triangular,  and only  elements  for which
C                             i.le.j  are referenced,
C  if   MATRIX = 'L' or 'l'   then  A  and  B  are  regarded  as   lower
C                             triangular,  and only  elements  for which
C                             i.ge.j  are referenced.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 21-November-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      INTEGER            I, J
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( ( MATRIX.EQ.'G' ).OR.( MATRIX.EQ.'g' ) )THEN
         DO 20 J = 1, N
            DO 10 I = 1, M
               B( I, J ) = A( I, J )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( ( MATRIX.EQ.'U' ).OR.( MATRIX.EQ.'u' ) )THEN
         DO 40 J = 1, N
            DO 30 I = 1, MIN( M, J )
               B( I, J ) = A( I, J )
   30       CONTINUE
   40    CONTINUE
      ELSE IF( ( MATRIX.EQ.'L' ).OR.( MATRIX.EQ.'l' ) )THEN
         DO 60 J = 1, MIN( M, N )
            DO 50 I = J, M
               B( I, J ) = A( I, J )
   50       CONTINUE
   60    CONTINUE
      END IF
C
      RETURN
C
C     End of F06QFF. ( SMCOPY )
C
      END
** END OF F06QFFTEXT
*UPTODATE F06QGFTEXT
      DOUBLE PRECISION FUNCTION F06QGF( NORM, MATRIX, M, N, A, LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER                           LDA, M, N
      CHARACTER*1                       MATRIX, NORM
C     .. Array Arguments ..
      DOUBLE PRECISION                  A( LDA, * )
C     ..
C
C  Purpose
C  =======
C
C  F06QGF  returns the value of the one norm,  or the Frobenius norm, or
C  the element of  largest absolute value of a real matrix A.  A  may be
C  rectangular,  or square,  or triangular,  or symmetric.
C
C  Description
C  ===========
C
C  F06QGF returns the value
C
C     F06QGF = ( max( abs( a( i, j ) ) ) , NORM = 'M' or 'm'
C              (
C              ( norm1( A )  ,             NORM = '1', 'O' or 'o'
C              (
C              ( normF( A ) ,              NORM = 'F', 'f', 'E' or 'e'
C
C  where norm1 denotes the one norm of a matrix (maximum column sum) and
C  normF denotes the  Frobenius norm of a matrix  (square root of sum of
C  squares).  Note that  max( abs( a( i, j ) ) )  is not a  matrix norm.
C
C  The type of matrix for which  F06QGF is returned is determined by the
C  parameter MATRIX.
C
C  If   MATRIX = 'G' or 'g'   then  A  is regarded as  a general matrix,
C  If   MATRIX = 'U' or 'u'   then  A  is regarded as  upper triangular,
C  If   MATRIX = 'L' or 'l'   then  A  is regarded as  lower triangular,
C  If   MATRIX = 'S' or 's'   then  A  is regarded as symmetric and only
C             or 'H' or 'h'   the  upper triangular part of the array  A
C                             is referenced,
C  If   MATRIX = 'Y' or 'y'   then  A  is regarded as symmetric and only
C             or 'E' or 'e'   the  lower triangular part of the array  A
C                             is referenced.
C
C  Parameters
C  ==========
C
C  NORM  -  CHARACTER*1.
C
C           On entry,  NORM specifies the value to be returned in F06QGF
C           as described above.
C
C           Unchanged on exit.
C
C  MATRIX - CHARACTER*1.
C
C           On entry,  MATRIX  specifies the type of matrix and,  in the
C           case of a  symmetric matrix,  the part of the array in which
C           the matrix is stored as described above.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C
C           On entry,  M  specifies the number of rows of the matrix  A.
C           M  must be at least  zero and when the  matrix is  symmetric
C           then  M must be equal to  N. When  M = 0  then F06QGF is set
C           to zero and an immediate return is effected.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the number of columns of the matrix A.
C           N  must be at least zero. When  N = 0  then F06QGF is set to
C           zero and an immediate return is effected.
C
C           Unchanged on exit.
C
C  A      - REAL array of DIMENSION ( LDA, n ).
C
C           Before entry,  A  must contain the  m by n  matrix for which
C           F06QGF is required.
C
C           If  MATRIX = 'U' or 'u' or 'S' or 's' or 'H' or 'h' then the
C           strictly lower triangular part of A is not referenced.
C
C           If  MATRIX = 'L' or 'l' or 'Y' or 'y' or 'E' or 'e' then the
C           strictly upper triangular part of A is not referenced.
C
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in  the  calling  (sub)  program.  LDA  must be at least  M.
C
C           Unchanged on exit.
C
C  Further comments
C  ================
C
C  If A is part of a matrix B partitioned as
C
C     B = ( B1  B2 ) ,
C         ( B3  A  )
C
C  where  B1 is an l by k matrix  ( l.ge.0, k.ge.0 ),  then this routine
C  may be called with the parameter  A as  b( l + 1, k + 1 ) and  LDA as
C  the first dimension of  B  as declared in the calling  (sub) program.
C
C  This routine  can be inefficient on  paged machines when the one norm
C  is required, the matrix is symmetric and N is large.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 13-January-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION                  ONE, ZERO
      PARAMETER                         ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION                  SCALE, SUM, VALUE
      INTEGER                           I, J
C     .. External Functions ..
      DOUBLE PRECISION                  F06BMF
      EXTERNAL                          F06BMF
C     .. External Subroutines ..
      EXTERNAL                          F06FJF
C     .. Intrinsic Functions ..
      INTRINSIC                         ABS, MAX, MIN
C     ..
C     .. Executable Statements ..
      IF( MIN( M, N ).EQ.0 )THEN
         VALUE = ZERO
      ELSE IF( ( NORM.EQ.'M' ).OR.( NORM.EQ.'m' ) )THEN
C
C        Find  max( abs( a( i, j ) ) ).
C
         VALUE = ZERO
         IF( ( MATRIX.EQ.'G' ).OR.( MATRIX.EQ.'g' ) )THEN
            DO 20 J = 1, N
               DO 10 I = 1, M
                  VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   10          CONTINUE
   20       CONTINUE
         ELSE IF( ( MATRIX.EQ.'U' ).OR.( MATRIX.EQ.'u' ).OR.
     $            ( MATRIX.EQ.'S' ).OR.( MATRIX.EQ.'s' ).OR.
     $            ( MATRIX.EQ.'H' ).OR.( MATRIX.EQ.'h' ) )THEN
            DO 40 J = 1, N
               DO 30 I = 1, MIN( M, J )
                  VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   30          CONTINUE
   40       CONTINUE
         ELSE IF( ( MATRIX.EQ.'L' ).OR.( MATRIX.EQ.'l' ).OR.
     $            ( MATRIX.EQ.'Y' ).OR.( MATRIX.EQ.'y' ).OR.
     $            ( MATRIX.EQ.'E' ).OR.( MATRIX.EQ.'e' ) )THEN
            DO 60 J = 1, MIN( M, N )
               DO 50 I = J, M
                  VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   50          CONTINUE
   60       CONTINUE
         END IF
      ELSE IF( ( NORM.EQ.'1' ).OR.( NORM.EQ.'O' ).OR.
     $         ( NORM.EQ.'o' ) )THEN
C
C        Find  norm1( A ).
C
         VALUE = ZERO
         IF( ( MATRIX.EQ.'G' ).OR.( MATRIX.EQ.'g' ) )THEN
            DO 80 J = 1, N
               SUM = ZERO
               DO 70 I = 1, M
                  SUM = SUM + ABS( A( I, J ) )
   70          CONTINUE
               VALUE = MAX( VALUE, SUM )
   80       CONTINUE
         ELSE IF( ( MATRIX.EQ.'U' ).OR.( MATRIX.EQ.'u' ) )THEN
            DO 100 J = 1, N
               SUM = ZERO
               DO 90 I = 1, MIN( M, J )
                  SUM = SUM + ABS( A( I, J ) )
   90          CONTINUE
               VALUE = MAX( VALUE, SUM )
  100       CONTINUE
         ELSE IF( ( MATRIX.EQ.'L' ).OR.( MATRIX.EQ.'l' ) )THEN
            DO 120 J = 1, MIN( M, N )
               SUM = ZERO
               DO 110 I = J, M
                  SUM = SUM + ABS( A( I, J ) )
  110          CONTINUE
               VALUE = MAX( VALUE, SUM )
  120       CONTINUE
         ELSE IF( ( MATRIX.EQ.'S' ).OR.( MATRIX.EQ.'s' ).OR.
     $            ( MATRIX.EQ.'H' ).OR.( MATRIX.EQ.'h' ) )THEN
            DO 150 J = 1, N
               SUM = ZERO
               DO 130 I = 1, J
                  SUM = SUM + ABS( A( I, J ) )
  130          CONTINUE
               DO 140 I = J + 1, N
                  SUM = SUM + ABS( A( J, I ) )
  140          CONTINUE
               VALUE = MAX( VALUE, SUM )
  150       CONTINUE
         ELSE IF( ( MATRIX.EQ.'Y' ).OR.( MATRIX.EQ.'y' ).OR.
     $            ( MATRIX.EQ.'E' ).OR.( MATRIX.EQ.'e' ) )THEN
            DO 180 J = 1, N
               SUM = ZERO
               DO 160 I = 1, J - 1
                  SUM = SUM + ABS( A( J, I ) )
  160          CONTINUE
               DO 170 I = J, N
                  SUM = SUM + ABS( A( I, J ) )
  170          CONTINUE
               VALUE = MAX( VALUE, SUM )
  180       CONTINUE
         END IF
      ELSE IF( ( NORM.EQ.'F' ).OR.( NORM.EQ.'f' ).OR.( NORM.EQ.'E' ).OR.
     $         ( NORM.EQ.'e' ) )THEN
C
C        Find  normF( A ).
C
         SCALE = ZERO
         SUM = ONE
         IF( ( MATRIX.EQ.'G' ).OR.( MATRIX.EQ.'g' ) )THEN
            DO 190 J = 1, N
               CALL F06FJF( M, A( 1, J ), 1, SCALE, SUM )
  190       CONTINUE
         ELSE IF( ( MATRIX.EQ.'U' ).OR.( MATRIX.EQ.'u' ) )THEN
            DO 200 J = 1, N
               CALL F06FJF( MIN( M, J ), A( 1, J ), 1, SCALE, SUM )
  200       CONTINUE
         ELSE IF( ( MATRIX.EQ.'L' ).OR.( MATRIX.EQ.'l' ) )THEN
            DO 210 J = 1, MIN( M, N )
               CALL F06FJF( M - J + 1, A( J, J ), 1, SCALE, SUM )
  210       CONTINUE
         ELSE IF( ( MATRIX.EQ.'S' ).OR.( MATRIX.EQ.'s' ).OR.
     $            ( MATRIX.EQ.'H' ).OR.( MATRIX.EQ.'h' ).OR.
     $            ( MATRIX.EQ.'Y' ).OR.( MATRIX.EQ.'y' ).OR.
     $            ( MATRIX.EQ.'E' ).OR.( MATRIX.EQ.'e' ) )THEN
            IF( ( MATRIX.EQ.'S' ).OR.( MATRIX.EQ.'s' ).OR.
     $          ( MATRIX.EQ.'H' ).OR.( MATRIX.EQ.'h' ) )THEN
               DO 220 J = 2, N
                  CALL F06FJF( J - 1, A( 1, J ), 1, SCALE, SUM )
  220          CONTINUE
            ELSE
               DO 230 J = 1, N - 1
                  CALL F06FJF( N - J, A( J + 1, J ), 1, SCALE, SUM )
  230          CONTINUE
            END IF
            SUM = 2*SUM
            CALL F06FJF( N, A( 1, 1 ), LDA + 1, SCALE, SUM )
         END IF
         VALUE = F06BMF( SCALE, SUM )
      END IF
C
      F06QGF = VALUE
      RETURN
C
C     End of F06QGF. ( SMNRM )
C
      END
** END OF F06QGFTEXT
*UPTODATE F06QHFTEXT
      SUBROUTINE F06QHF( MATRIX, M, N, CONST, DIAG, A, LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      CHARACTER*1        MATRIX
      DOUBLE PRECISION   CONST, DIAG
      INTEGER            LDA, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
C     ..
C
C  F06QHF forms the m by n matrix A given by
C
C     a( i, j ) = (  diag  i.eq.j,
C                 (
C                 ( const  i.ne.j.
C
C  If   MATRIX = 'G' or 'g'   then  A  is regarded  as a general matrix,
C  if   MATRIX = 'U' or 'u'   then  A  is regarded  as upper triangular,
C                             and only  elements  for which  i.le.j  are
C                             referenced,
C  if   MATRIX = 'L' or 'l'   then  A  is regarded  as lower triangular,
C                             and only  elements  for which  i.ge.j  are
C                             referenced.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 21-November-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      INTEGER            I, J
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( ( MATRIX.EQ.'G' ).OR.( MATRIX.EQ.'g' ) )THEN
         DO 20 J = 1, N
            DO 10 I = 1, M
               A( I, J ) = CONST
   10       CONTINUE
   20    CONTINUE
         IF( CONST.NE.DIAG )THEN
            DO 30 I = 1, MIN( M, N )
               A( I, I ) = DIAG
   30       CONTINUE
         END IF
      ELSE IF( ( MATRIX.EQ.'U' ).OR.( MATRIX.EQ.'u' ) )THEN
         DO 50 J = 1, N
            DO 40 I = 1, MIN( M, J )
               A( I, J ) = CONST
   40       CONTINUE
   50    CONTINUE
         IF( CONST.NE.DIAG )THEN
            DO 60 I = 1, MIN( M, N )
               A( I, I ) = DIAG
   60       CONTINUE
         END IF
      ELSE IF( ( MATRIX.EQ.'L' ).OR.( MATRIX.EQ.'l' ) )THEN
         DO 80 J = 1, MIN( M, N )
            DO 70 I = J, M
               A( I, J ) = CONST
   70       CONTINUE
   80    CONTINUE
         IF( CONST.NE.DIAG )THEN
            DO 90 I = 1, MIN( M, N )
               A( I, I ) = DIAG
   90       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06QHF. ( SMLOAD )
C
      END
** END OF F06QHFTEXT
*UPTODATE F06QJFTEXT
      SUBROUTINE F06QJF( SIDE, TRANS, N, PERM, K, B, LDB )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K, LDB, N
      CHARACTER*1        SIDE, TRANS
C     .. Array Arguments ..
      DOUBLE PRECISION   B( LDB, * )
      INTEGER            PERM( * )
C     ..
C
C  Purpose
C  =======
C
C  F06QJF performs one of the transformations
C
C     B := P'*B   or   B := P*B,   where B is an m by k matrix,
C
C  or
C
C     B := B*P'   or   B := B*P,   where B is a k by m matrix,
C
C  P being an m by m permutation matrix of the form
C
C     P = P( 1, index( 1 ) )*P( 2, index( 2 ) )*...*P( n, index( n ) ),
C
C  where  P( i, index( i ) ) is the permutation matrix that interchanges
C  items i and index( i ). That is P( i, index( i ) ) is the unit matrix
C  with rows and columns  i and index( i )  interchanged.  Of course, if
C  index( i ) = i  then  P( i, index( i ) ) = I.
C
C  This routine  is intended for use in  conjunction with  Nag auxiliary
C  routines that  perform  interchange  operations,  such  as  pivoting.
C
C  Parameters
C  ==========
C
C  SIDE   - CHARACTER*1.
C  TRANS
C           On entry,  SIDE  ( Left-hand side, or Right-hand side )  and
C           TRANS  ( Transpose, or No transpose )  specify the operation
C           to be performed as follows.
C
C           SIDE = 'L' or 'l'   and   TRANS = 'T' or 't'
C
C              Perform the operation   B := P'*B.
C
C           SIDE = 'L' or 'l'   and   TRANS = 'N' or 'n'
C
C              Perform the operation   B := P*B.
C
C           SIDE = 'R' or 'r'   and   TRANS = 'T' or 't'
C
C              Perform the operation   B := B*P'.
C
C           SIDE = 'R' or 'r'   and   TRANS = 'N' or 'n'
C
C              Perform the operation   B := B*P.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C
C           On entry,  N must specify the value of n. N must be at least
C           zero.  When  N = 0  then an  immediate  return is  effected.
C
C           Unchanged on exit.
C
C  PERM   - INTEGER array of DIMENSION at least ( n ).
C
C           Before  entry,  PERM  must  contain the  n  indices  for the
C           permutation matrices. index( i ) must satisfy
C
C              1 .le. index( i ) .le. m.
C
C           It is usual for index( i ) to be at least i, but this is not
C           necessary for this routine.
C
C           Unchanged on exit.
C
C  K      - INTEGER.
C
C           On entry with  SIDE = 'L' or 'l',  K must specify the number
C           of columns of B and on entry with  SIDE = 'R' or 'r', K must
C           specify the  number of rows of B.  K must be at least  zero.
C           When  K = 0  then an immediate return is effected.
C
C           Unchanged on exit.
C
C  B      - REAL  array of  DIMENSION  ( LDB, ncolb ),  where  ncolb = k
C           when  SIDE = 'L' or 'l'  and  ncolb = m  when  SIDE = 'R' or
C           'r'.
C
C           Before entry  with  SIDE = 'L' or 'l',  the  leading  m by K
C           part  of  the  array   B  must  contain  the  matrix  to  be
C           transformed  and  before entry with  SIDE = 'R' or 'r',  the
C           leading  K by m part of the array  B must contain the matrix
C           to  be  transformed.  On  exit,  B  is  overwritten  by  the
C           transformed matrix.
C
C  LDB    - INTEGER.
C
C           On entry,  LDB  must specify  the  leading dimension  of the
C           array  B  as declared  in the  calling  (sub) program.  When
C           SIDE = 'L' or 'l'   then  LDB  must  be  at  least  m,  when
C           SIDE = 'R' or 'r'   then  LDB  must  be  at  least  k.
C           Unchanged on exit.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 13-January-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, J, L
      LOGICAL            LEFT, NULL, RIGHT, TRNSP
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( MIN( N, K ).EQ.0 )
     $   RETURN
      LEFT = ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' )
      RIGHT = ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' )
      NULL = ( TRANS.EQ.'N' ).OR.( TRANS.EQ.'n' )
      TRNSP = ( TRANS.EQ.'T' ).OR.( TRANS.EQ.'t' )
      IF( LEFT )THEN
         IF( TRNSP )THEN
            DO 20 I = 1, N
               IF( PERM( I ).NE.I )THEN
                  L = PERM( I )
                  DO 10 J = 1, K
                     TEMP = B( I, J )
                     B( I, J ) = B( L, J )
                     B( L, J ) = TEMP
   10             CONTINUE
               END IF
   20       CONTINUE
         ELSE IF( NULL )THEN
            DO 40 I = N, 1, -1
               IF( PERM( I ).NE.I )THEN
                  L = PERM( I )
                  DO 30 J = 1, K
                     TEMP = B( L, J )
                     B( L, J ) = B( I, J )
                     B( I, J ) = TEMP
   30             CONTINUE
               END IF
   40       CONTINUE
         END IF
      ELSE IF( RIGHT )THEN
         IF( TRNSP )THEN
            DO 60 J = N, 1, -1
               IF( PERM( J ).NE.J )THEN
                  L = PERM( J )
                  DO 50 I = 1, K
                     TEMP = B( I, J )
                     B( I, J ) = B( I, L )
                     B( I, L ) = TEMP
   50             CONTINUE
               END IF
   60       CONTINUE
         ELSE IF( NULL )THEN
            DO 80 J = 1, N
               IF( PERM( J ).NE.J )THEN
                  L = PERM( J )
                  DO 70 I = 1, K
                     TEMP = B( I, L )
                     B( I, L ) = B( I, J )
                     B( I, J ) = TEMP
   70             CONTINUE
               END IF
   80       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06QJF. ( SGEAP )
C
      END
** END OF F06QJFTEXT
*UPTODATE F06QKFTEXT
      SUBROUTINE F06QKF( SIDE, TRANS, N, PERM, K, B, LDB )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K, LDB, N
      CHARACTER*1        SIDE, TRANS
C     .. Array Arguments ..
      DOUBLE PRECISION   PERM( * ), B( LDB, * )
C     ..
C
C  Purpose
C  =======
C
C  F06QKF performs one of the transformations
C
C     B := P'*B   or   B := P*B,   where B is an m by k matrix,
C
C  or
C
C     B := B*P'   or   B := B*P,   where B is a k by m matrix,
C
C  P being an m by m permutation matrix of the form
C
C     P = P( 1, index( 1 ) )*P( 2, index( 2 ) )*...*P( n, index( n ) ),
C
C  where  P( i, index( i ) ) is the permutation matrix that interchanges
C  items i and index( i ). That is P( i, index( i ) ) is the unit matrix
C  with rows and columns  i and  index( i )  interchanged. Of course, if
C  index( i ) = i  then  P( i, index( i ) ) = I.
C
C  This  routine is intended  for use in conjunction with  Nag auxiliary
C  routines  that  perform  interchange  operations,  such  as  sorting.
C
C  Parameters
C  ==========
C
C  SIDE   - CHARACTER*1.
C  TRANS
C           On entry,  SIDE  ( Left-hand side, or Right-hand side )  and
C           TRANS  ( Transpose, or No transpose )  specify the operation
C           to be performed as follows.
C
C           SIDE = 'L' or 'l'   and   TRANS = 'T' or 't'
C
C              Perform the operation   B := P'*B.
C
C           SIDE = 'L' or 'l'   and   TRANS = 'N' or 'n'
C
C              Perform the operation   B := P*B.
C
C           SIDE = 'R' or 'r'   and   TRANS = 'T' or 't'
C
C              Perform the operation   B := B*P'.
C
C           SIDE = 'R' or 'r'   and   TRANS = 'N' or 'n'
C
C              Perform the operation   B := B*P.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C
C           On entry, N must specify the value of n.  N must be at least
C           zero.  When  N = 0  then an  immediate  return  is effected.
C
C           Unchanged on exit.
C
C  PERM   - REAL             array of DIMENSION at least ( n ).
C
C           Before  entry,  PERM  must  contain  the  n indices  for the
C           permutation matrices. index( i ) must satisfy
C
C              1 .le. index( i ) .le. m.
C
C           It is usual for index( i ) to be at least i, but this is not
C           necessary for this routine. It is assumed that the statement
C           INDEX = PERM( I )  returns the correct integer in  INDEX, so
C           that,  if necessary,  PERM( I )  should contain a real value
C           slightly larger than  INDEX.
C
C           Unchanged on exit.
C
C  K      - INTEGER.
C
C           On entry with  SIDE = 'L' or 'l',  K must specify the number
C           of columns of B and on entry with  SIDE = 'R' or 'r', K must
C           specify the number of rows of  B.  K must be at least  zero.
C           When  K = 0  then an immediate return is effected.
C
C           Unchanged on exit.
C
C  B      - REAL  array  of  DIMENSION ( LDB, ncolb ),  where  ncolb = k
C           when  SIDE = 'L' or 'l'  and  ncolb = m  when  SIDE = 'R' or
C           'r'.
C
C           Before entry  with  SIDE = 'L' or 'l',  the  leading  m by K
C           part  of  the  array   B  must  contain  the  matrix  to  be
C           transformed  and before  entry with  SIDE = 'R' or 'r',  the
C           leading  K by m part of the array  B must contain the matrix
C           to  be  transformed.  On exit,   B  is  overwritten  by  the
C           transformed matrix.
C
C  LDB    - INTEGER.
C
C           On entry,  LDB  must specify  the  leading dimension  of the
C           array  B  as declared  in the  calling  (sub) program.  When
C           SIDE = 'L' or 'l'   then  LDB  must  be  at  least  m,  when
C           SIDE = 'R' or 'r'   then  LDB  must  be  at  least  k.
C           Unchanged on exit.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 11-August-1987.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      LOGICAL            LEFT, NULL, RIGHT, TRNSP
      INTEGER            I, J, L
      DOUBLE PRECISION   TEMP
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( MIN( N, K ).EQ.0 )
     $   RETURN
      LEFT = ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' )
      RIGHT = ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' )
      NULL = ( TRANS.EQ.'N' ).OR.( TRANS.EQ.'n' )
      TRNSP = ( TRANS.EQ.'T' ).OR.( TRANS.EQ.'t' )
      IF( LEFT )THEN
         IF( TRNSP )THEN
            DO 20 I = 1, N
               L = PERM( I )
               IF( L.NE.I )THEN
                  DO 10 J = 1, K
                     TEMP = B( I, J )
                     B( I, J ) = B( L, J )
                     B( L, J ) = TEMP
   10             CONTINUE
               END IF
   20       CONTINUE
         ELSE IF( NULL )THEN
            DO 40 I = N, 1, -1
               L = PERM( I )
               IF( L.NE.I )THEN
                  DO 30 J = 1, K
                     TEMP = B( L, J )
                     B( L, J ) = B( I, J )
                     B( I, J ) = TEMP
   30             CONTINUE
               END IF
   40       CONTINUE
         END IF
      ELSE IF( RIGHT )THEN
         IF( TRNSP )THEN
            DO 60 J = N, 1, -1
               L = PERM( J )
               IF( L.NE.J )THEN
                  DO 50 I = 1, K
                     TEMP = B( I, J )
                     B( I, J ) = B( I, L )
                     B( I, L ) = TEMP
   50             CONTINUE
               END IF
   60       CONTINUE
         ELSE IF( NULL )THEN
            DO 80 J = 1, N
               L = PERM( J )
               IF( L.NE.J )THEN
                  DO 70 I = 1, K
                     TEMP = B( I, L )
                     B( I, L ) = B( I, J )
                     B( I, J ) = TEMP
   70             CONTINUE
               END IF
   80       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06QKF. ( SGEAPR )
C
      END
** END OF F06QKFTEXT
*UPTODATE F06QNFTEXT
      SUBROUTINE F06QNF( SIDE, N, K1, K2, S, A, LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K1, K2, LDA, N
      CHARACTER*1        SIDE
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), S( * )
C     ..
C
C  F06QNF applies a  sequence  of  pairwise interchanges to either  the
C  left,  or the right,  of the  n by n  upper triangular matrix  U,  to
C  transform U to an  upper Hessenberg matrix. The interchanges are
C  applied in planes k1 up to k2.
C
C  The upper Hessenberg matrix, H, is formed as
C
C     H = P*U,    when   SIDE = 'L' or 'l',  (  Left-hand side )
C
C  where P is a permutation matrix of the form
C
C     P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 )
C
C  and is formed as
C
C     H = U*P',   when   SIDE = 'R' or 'r',  ( Right-hand side )
C
C  where P is a permutation matrix of the form
C
C     P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C  P( k ) being a pairwise interchange for the  ( k, k + 1 ) plane.
C  The  two by two
C  interchange part of P( k ), R( k ), is assumed to have the form
C
C     R( k ) = ( 0  1 ).
C              ( 1  0 )
C
C  The matrix  U must be supplied in the n by n leading upper triangular
C  part of the array  A, and this is overwritten by the upper triangular
C  part of  H.
C
C  The  sub-diagonal elements of  H, h( k + 1, k ),  are returned in the
C  elements s( k ),  k = k1, k1 + 1, ..., k2 - 1.
C
C  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
C  greater than n then an immediate return is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 16-May-1988.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   AIJ, TEMP
      INTEGER            I, J
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( ( MIN( N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.
     $   ( K2.GT.N ) )RETURN
      IF( ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' ) )THEN
C
C        Apply the permutations to columns n back to k1.
C
         DO 20 J = N, K1, -1
            IF( J.GE.K2 )THEN
               AIJ = A( K2, J )
            ELSE
C
C              Form  the  additional sub-diagonal element  h( j + 1, j )
C              and store it in s( j ).
C
               AIJ    = ZERO
               S( J ) = A( J, J )
            END IF
            DO 10 I = MIN( K2, J ) - 1, K1, -1
               TEMP          = A( I, J )
               A( I + 1, J ) = TEMP
               AIJ           = AIJ
   10       CONTINUE
            A( K1, J ) = AIJ
   20    CONTINUE
      ELSE IF( ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' ) )THEN
C
C        Apply  the  plane interchanges to  columns  k1  up to
C        ( k2 - 1 ) and  form   the   additional  sub-diagonal
C        elements,   storing  h( j + 1, j ) in s( j ).
C
         DO 40 J = K1, K2 - 1
            DO 30 I = 1, J
               TEMP = A( I, J + 1 )
               A( I, J + 1 ) = A( I, J )
               A( I, J )     = TEMP
   30       CONTINUE
            S( J )            = A( J + 1, J + 1 )
            A( J + 1, J + 1 ) = ZERO
   40    CONTINUE
      END IF
C
      RETURN
C
C     End of F06QNF. ( SUTSRH )
C
      END
** END OF F06QNFTEXT
*UPTODATE F06QRFTEXT
      SUBROUTINE F06QRF( SIDE, N, K1, K2, C, S, A, LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K1, K2, LDA, N
      CHARACTER*1        SIDE
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
C     ..
C
C  F06QRF restores an upper Hessenberg matrix H to upper triangular form
C  by  applying a sequence of  plane rotations  from either the left, or
C  the right.  The matrix  H  is assumed to have  non-zero  sub-diagonal
C  elements  in  positions  h( k + 1, k ),  k = k1, k1 + 1, ..., k2 - 1,
C  only  and  h( k + 1, k )  must  be  supplied  in  s( k ).
C
C  H is restored to the upper triangular matrix R either as
C
C     R = P*H,   when   SIDE = 'L' or 'l'  (  Left-hand side )
C
C  where P is an orthogonal matrix of the form
C
C     P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C  or as
C
C     R = H*P',  when   SIDE = 'R' or 'r'  ( Right-hand side )
C
C  where P is an orthogonal matrix of the form
C
C     P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C  in both cases  P( k )  being a  plane rotation  for the  ( k, k + 1 )
C  plane.  The cosine and sine that define P( k ) are returned in c( k )
C  and  s( k )  respectively.  The two by two  rotation part of  P( k ),
C  Q( k ), is of the form
C
C     Q( k ) = (  c( k )  s( k ) ).
C              ( -s( k )  c( k ) )
C
C  The upper triangular part of the matrix  H  must be supplied in the n
C  by n  leading upper triangular part of  A, and this is overwritten by
C  the upper triangular matrix R.
C
C  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
C  greater than n then an immediate return is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 13-January-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   AIJ, CTEMP, STEMP, SUBH, TEMP
      INTEGER            I, J
C     .. External Subroutines ..
      EXTERNAL           F06BAF
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( ( MIN( N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.
     $   ( K2.GT.N ) )RETURN
      IF( ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' ) )THEN
C
C        Restore   H  to  upper  triangular  form  by  annihilating  the
C        sub-diagonal elements of H.  The jth rotation is chosen so that
C
C           ( h( j, j ) ) := (  c  s )*( h( j, j )     ).
C           (     0     )    ( -s  c ) ( h( j + 1, j ) )
C
C        Apply the rotations in columns k1 up to n.
C
         DO 20 J = K1, N
            AIJ = A( K1, J )
            DO 10 I = K1, MIN( J, K2 ) - 1
               TEMP = A( I + 1, J )
               A( I, J ) = S( I )*TEMP + C( I )*AIJ
               AIJ = C( I )*TEMP - S( I )*AIJ
   10       CONTINUE
            IF( J.LT.K2 )THEN
C
C              Set up the rotation.
C
               SUBH = S( J )
               CALL F06BAF( AIJ, SUBH, C( J ), S( J ) )
               A( J, J ) = AIJ
            ELSE
               A( K2, J ) = AIJ
            END IF
   20    CONTINUE
      ELSE IF( ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' ) )THEN
C
C        Restore   H  to  upper  triangular  form  by  annihilating  the
C        sub-diagonal elements of H.  The jth rotation is chosen so that
C
C           ( h( j + 1, j + 1 ) ) := (  c  s )*( h( j + 1, j + 1 ) ),
C           (         0         )    ( -s  c ) ( h( j + 1, j )     )
C
C        which can be expressed as
C
C           ( 0  h( j + 1, j + 1 ) ) :=
C
C               ( h( j + 1, j )  h( j + 1, j + 1 ) )*(  c  s ).
C                                                    ( -s  c )
C
C        Thus we return  c( j ) = c  and  s( j ) = -s  to make the plane
C        rotation matrix look like
C
C           Q( j ) = (  c( j )  s( j ) ).
C                    ( -s( j )  c( j ) )
C
         DO 40 J = K2 - 1, K1, -1
            SUBH = S( J )
            CALL F06BAF( A( J + 1, J + 1 ), SUBH, CTEMP, STEMP )
            STEMP = -STEMP
            S( J ) = STEMP
            C( J ) = CTEMP
            IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
               DO 30 I = J, 1, -1
                  TEMP = A( I, J + 1 )
                  A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                  A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
   30          CONTINUE
            END IF
   40    CONTINUE
      END IF
C
      RETURN
C
C     End of F06QRF. ( SUHQR )
C
      END
** END OF F06QRFTEXT
*UPTODATE F06QSFTEXT
      SUBROUTINE F06QSF( SIDE, N, K1, K2, C, S, A, LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K1, K2, LDA, N
      CHARACTER*1        SIDE
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
C     ..
C
C  F06QSF restores an upper spiked matrix  H to upper triangular form by
C  applying a sequence of plane rotations, in planes  k1 up to k2,  from
C  either the left, or the right.
C
C  The matrix  H is assumed to have non-zero elements only in the spiked
C  positions, h( k2, k ) for a row spike and h( k + 1, k1 ) for a column
C  spike, k = k1, k1 + 1, ..., k2 - 1, and these must be supplied in the
C  elements s( k ).
C
C  When  SIDE = 'L' or 'l'  (  Left-hand side )
C
C     H  is  assumed  to have a  row spike  and is restored to the upper
C     triangular matrix  R as
C
C        R = P*H,
C
C     where P is an orthogonal matrix of the form
C
C        P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C     P( k )  being a  plane rotation  matrix for the  ( k, k2 )  plane.
C
C  When  SIDE = 'R' or 'r'  ( Right-hand side )
C
C     H  is assumed to have a  column spike and is restored to the upper
C     triangular matrix R as
C
C        R = H*P',
C
C     where P is an orthogonal matrix of the form
C
C        P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C     P( k ) being a plane rotation matrix for the  ( k1, k + 1 ) plane.
C
C  The  two by two  rotation  part of  P( k ),  Q( k ),  is of  the form
C
C     Q( k ) = (  c( k )  s( k ) )
C              ( -s( k )  c( k ) )
C
C  and  c( k ) and s( k ) are returned in the kth elements of the arrays
C  C and S respectively.
C
C  The upper triangular part of the matrix  H must be supplied in the  n
C  by n  leading upper triangular part of  A, and this is overwritten by
C  the upper triangular matrix R.
C
C  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
C  greater than n then an immediate return is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 13-January-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   AIJ, CTEMP, SPIKE, STEMP, TEMP
      INTEGER            I, J
C     .. External Subroutines ..
      EXTERNAL           F06BAF
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( ( MIN( N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.
     $   ( K2.GT.N ) )RETURN
      IF( ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' ) )THEN
C
C        Restore H to upper triangular form by annihilating the elements
C        in  the  spike  of  H.  The  jth rotation  is  chosen  so  that
C
C        ( h( j, j ) ) := (  c  s )*( h( j , j ) ).
C        (     0     )    ( -s  c ) ( h( k2, j ) )
C
C        Apply the rotations in columns k1 up to ( k2 - 1 ).
C
         DO 20 J = K1, K2 - 1
            SPIKE = S( J )
            DO 10 I = K1, J - 1
               AIJ = A( I, J )
               A( I, J ) = S( I )*SPIKE + C( I )*AIJ
               SPIKE = C( I )*SPIKE - S( I )*AIJ
   10       CONTINUE
C
C           Set up the rotation.
C
            CALL F06BAF( A( J, J ), SPIKE, C( J ), S( J ) )
   20    CONTINUE
C
C        Apply the rotations to columns k2 up to n.
C
         DO 40 J = K2, N
            TEMP = A( K2, J )
            DO 30 I = K1, K2 - 1
               AIJ = A( I, J )
               A( I, J ) = S( I )*TEMP + C( I )*AIJ
               TEMP = C( I )*TEMP - S( I )*AIJ
   30       CONTINUE
            A( K2, J ) = TEMP
   40    CONTINUE
      ELSE IF( ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' ) )THEN
C
C        Restore H to upper triangular form by annihilating the spike of
C        H. The jth rotation is chosen so that
C
C           ( h( j, j ) ) := (  c  s )*( h( j, j )  ),
C           (     0     )    ( -s  c ) ( h( j, k1 ) )
C
C        which can be expressed as
C
C           ( 0  h( j, j ) ) := ( h( j, k1 )  h( j, j ) )*(  c  s ).
C                                                         ( -s  c )
C
C        Thus we return  c( j ) = c  and  s( j ) = -s  to make the plane
C        rotation matrix look like
C
C           Q( j ) = (  c( j )  s( j ) ).
C                    ( -s( j )  c( j ) )
C
         DO 70 J = K2, K1 + 1, -1
            CALL F06BAF( A( J, J ), S( J - 1 ), CTEMP, STEMP )
            STEMP = -STEMP
            S( J - 1 ) = STEMP
            C( J - 1 ) = CTEMP
            IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
               DO 50 I = J - 1, K1 + 1, -1
                  SPIKE = S( I - 1 )
                  S( I - 1 ) = STEMP*A( I, J ) + CTEMP*SPIKE
                  A( I, J ) = CTEMP*A( I, J ) - STEMP*SPIKE
   50          CONTINUE
               DO 60 I = K1, 1, -1
                  TEMP = A( I, K1 )
                  A( I, K1 ) = STEMP*A( I, J ) + CTEMP*TEMP
                  A( I, J ) = CTEMP*A( I, J ) - STEMP*TEMP
   60          CONTINUE
            END IF
   70    CONTINUE
      END IF
C
      RETURN
C
C     End of F06QSF. ( SUSQR )
C
      END
** END OF F06QSFTEXT
*UPTODATE F06QTFTEXT
      SUBROUTINE F06QTF( SIDE, N, K1, K2, C, S, A, LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K1, K2, LDA, N
      CHARACTER*1        SIDE
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
C     ..
C
C  F06QTF performs the transformation
C
C     R := P*U*Q'  when  SIDE = 'L' or 'l'  (  Left-hand side )
C
C     R := Q*U*P'  when  SIDE = 'R' or 'r'  ( Right-hand side ),
C
C  where  U and R  are  n by n  upper  triangular  matrices,   P  is  an
C  orthogonal matrix,  consisting of a given sequence of plane rotations
C  to be  applied  in  planes  k1 to k2,  and  Q  is  a  unitary  matrix
C  consisting of a sequence of plane rotations, applied in planes  k1 to
C  k2,  chosen to make  R  upper triangular.
C
C  When  SIDE = 'L' or 'l'  then  P  is  given  as a  sequence of  plane
C  rotation matrices
C
C     P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C  where  P( k ) is a plane rotation matrix for the  ( k, k + 1 ) plane.
C  In this case the matrix Q is given as
C
C     Q = Q( k2 - 1 )*...*Q( k1 + 1 )*Q( k1 ),
C
C  where  Q( k ) is a plane rotation matrix for the  ( k, k + 1 ) plane.
C
C  When  SIDE = 'R' or 'r'  then  P  is  given  as a  sequence of  plane
C  rotation matrices
C
C     P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C  where  P( k ) is a plane rotation matrix for the  ( k, k + 1 ) plane.
C  In this case the matrix Q is given as
C
C     Q = Q( k1 )*Q( k1 + 1 )*...*Q( k2 - 1 ),
C
C  where  Q( k ) is a plane rotation matrix for the  ( k, k + 1 ) plane.
C
C  The  upper  triangular  matrix  U  must  be  supplied  in the  n by n
C  leading upper triangular part of  A,  and this  is overwritten by the
C  upper triangular matrix  R.  The cosine  and  sine  that  define  the
C  plane rotation matrix  P( k )  must be supplied in  c( k ) and s( k )
C  respectively,  and  the two by two rotation part of  P( k ),  T( k ),
C  is assumed to be of the form
C
C     T( k ) = (  c( k )  s( k ) ).
C              ( -s( k )  c( k ) )
C
C  The cosine  and  sine that define  Q( k )  are overwritten on  c( k )
C  and  s( k )  respectively and the two by two rotation part of  Q( k )
C  will have the form of  T( k )  above.
C
C  If  n or k1  are less  than  unity, or  k1  is not  less than  k2, or
C  k2  is greater than  n  then an immediate return is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 26-November-1987.
C     Sven Hammarling and Mick Pont, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   AIJ, CTEMP, FILL, STEMP, TEMP
      INTEGER            I, I1, J
C     .. External Subroutines ..
      EXTERNAL           F06BAF
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( ( MIN( N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.
     $   ( K2.GT.N ) )RETURN
      IF( ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' ) )THEN
C
C        Apply the left-hand transformations,  column by column,  to the
C        triangular part of  U,  but not to  anywhere  that would  cause
C        fill.
C
         DO 20 J = K1 + 1, N
C
C           Apply  P( k1 ) ... P( j - 1 )  to column j.
C
            AIJ = A( K1, J )
            DO 10 I = K1, MIN( J - 1, K2 - 1 )
               A( I, J ) = S( I )*A( I + 1, J ) + C( I )*AIJ
               AIJ = C( I )*A( I + 1, J ) - S( I )*AIJ
   10       CONTINUE
            A( I, J ) = AIJ
   20    CONTINUE
C
C           Now apply each  left-hand tranformation  to form the fill-in
C           elements and apply a  right-hand transformation to eliminate
C           the fill-in element.
C
         DO 40 J = K1, K2 - 1
C
C           Apply  P( j )  to the jth diagonal element  and the  fill-in
C           position.
C
            FILL = -S( J )*A( J, J )
            A( J, J ) = C( J )*A( J, J )
C
C           Now  set up  the rotation  Q( j )  to eliminate the  fill-in
C           element,  and  apply  Q( j )  to  the  jth  and  ( j + 1 )th
C           columns.
C
            CALL F06BAF( A( J + 1, J + 1 ), FILL, CTEMP, STEMP )
            C( J ) = CTEMP
            S( J ) = -STEMP
            IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
               STEMP = -STEMP
               DO 30 I = 1, J
                  TEMP = A( I, J + 1 )
                  A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                  A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
   30          CONTINUE
            END IF
   40    CONTINUE
      ELSE IF( ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' ) )THEN
C
C        We intermingle the  left and right hand transformations so that
C        at the kth step we form
C
C           A := Q( k )*A*P( k )'.
C
C        First  apply  the  transformations  in  columns  k2 back to k1.
C
         DO 60 J = K2 - 1, K1, -1
C
C           First apply  P( j ).
C
            IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
               CTEMP = C( J )
               STEMP = S( J )
               DO 50 I = 1, J
                  TEMP = A( I, J + 1 )
                  A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                  A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
   50          CONTINUE
C
C              Next form the fill-in element  a( j + 1, j )  by applying
C              P( j ).
C
               FILL = S( J )*A( J + 1, J + 1 )
               A( J + 1, J + 1 ) = C( J )*A( J + 1, J + 1 )
C
C              Now set up the rotation  Q( j )  to eliminate the fill-in
C              element.
C
               CALL F06BAF( A( J, J ), FILL, C( J ), S( J ) )
            END IF
   60    CONTINUE
C
C        Finally  apply  Q( k2 - 1 ) ... Q( k1 )  to columns  n  back to
C        ( k1 + 1 ).
C
         DO 80 J = N, K1 + 1, -1
            I1 = MIN( K2, J )
            AIJ = A( I1, J )
            DO 70 I = I1 - 1, K1, -1
               TEMP = A( I, J )
               A( I + 1, J ) = C( I )*AIJ - S( I )*TEMP
               AIJ = S( I )*AIJ + C( I )*TEMP
   70       CONTINUE
            A( K1, J ) = AIJ
   80    CONTINUE
      END IF
      RETURN
C
C     End of F06QTF. ( SUTSQR )
C
      END
** END OF F06QTFTEXT
*UPTODATE F06QVFTEXT
      SUBROUTINE F06QVF( SIDE, N, K1, K2, C, S, A, LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K1, K2, LDA, N
      CHARACTER*1        SIDE
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
C     ..
C
C  F06QVF applies a  given sequence  of  plane rotations  to either  the
C  left,  or the right,  of the  n by n  upper triangular matrix  U,  to
C  transform U to an  upper Hessenberg matrix. The rotations are applied
C  in planes k1 up to k2.
C
C  The upper Hessenberg matrix, H, is formed as
C
C     H = P*U,    when   SIDE = 'L' or 'l',  (  Left-hand side )
C
C  where P is an orthogonal matrix of the form
C
C     P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 )
C
C  and is formed as
C
C     H = U*P',   when   SIDE = 'R' or 'r',  ( Right-hand side )
C
C  where P is an orthogonal matrix of the form
C
C     P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C  P( k ) being a plane rotation matrix for the  ( k, k + 1 ) plane. The
C  cosine and sine that define P( k ), k = k1, k1 + 1, ..., k2 - 1, must
C  be  supplied  in  c( k )  and  s( k )  respectively.  The  two by two
C  rotation part of P( k ), R( k ), is assumed to have the form
C
C     R( k ) = (  c( k )  s( k ) ).
C              ( -s( k )  c( k ) )
C
C  The matrix  U must be supplied in the n by n leading upper triangular
C  part of the array  A, and this is overwritten by the upper triangular
C  part of  H.
C
C  The  sub-diagonal elements of  H, h( k + 1, k ),  are returned in the
C  elements s( k ),  k = k1, k1 + 1, ..., k2 - 1.
C
C  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
C  greater than n then an immediate return is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 13-January-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   AIJ, CTEMP, STEMP, TEMP
      INTEGER            I, J
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( ( MIN( N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.
     $   ( K2.GT.N ) )RETURN
      IF( ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' ) )THEN
C
C        Apply the plane rotations to columns n back to k1.
C
         DO 20 J = N, K1, -1
            IF( J.GE.K2 )THEN
               AIJ = A( K2, J )
            ELSE
C
C              Form  the  additional sub-diagonal element  h( j + 1, j )
C              and store it in s( j ).
C
               AIJ = C( J )*A( J, J )
               S( J ) = -S( J )*A( J, J )
            END IF
            DO 10 I = MIN( K2, J ) - 1, K1, -1
               TEMP = A( I, J )
               A( I + 1, J ) = C( I )*AIJ - S( I )*TEMP
               AIJ = S( I )*AIJ + C( I )*TEMP
   10       CONTINUE
            A( K1, J ) = AIJ
   20    CONTINUE
      ELSE IF( ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' ) )THEN
C
C        Apply  the  plane rotations  to  columns  k1  up to  ( k2 - 1 )
C        and  form   the   additional  sub-diagonal  elements,   storing
C        h( j + 1, j ) in s( j ).
C
         DO 40 J = K1, K2 - 1
            IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
               STEMP = S( J )
               CTEMP = C( J )
               DO 30 I = 1, J
                  TEMP = A( I, J + 1 )
                  A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                  A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
   30          CONTINUE
               S( J ) = STEMP*A( J + 1, J + 1 )
               A( J + 1, J + 1 ) = CTEMP*A( J + 1, J + 1 )
            END IF
   40    CONTINUE
      END IF
C
      RETURN
C
C     End of F06QVF. ( SUTSRH )
C
      END
** END OF F06QVFTEXT
*UPTODATE F06QWFTEXT
      SUBROUTINE F06QWF( SIDE, N, K1, K2, C, S, A, LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K1, K2, LDA, N
      CHARACTER*1        SIDE
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
C     ..
C
C  F06QWF applies a  given sequence  of  plane rotations  to either  the
C  left,  or the right,  of the  n by n  upper triangular  matrix  U  to
C  transform  U  to an upper spiked matrix. The rotations are applied in
C  planes k1 up to k2.
C
C  The upper spiked matrix, H, is formed as
C
C     H = P*U,   when   SIDE = 'L' or 'l',  ( Left-hand side )
C
C  where P is an orthogonal matrix of the form
C
C     P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C  P( k ) being a plane rotation matrix for the ( k, k2 ) plane, and is
C  formed as
C
C     H = U*P',   when   SIDE = 'R' or 'r',  ( Right-hand side )
C
C  where P is an orthogonal matrix of the form
C
C     P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C  P( k )  being a  plane rotation matrix for the  ( k1, k + 1 )  plane.
C
C  The cosine and sine that define  P( k ), k = k1, k1 + 1, ..., k2 - 1,
C  must be  supplied  in  c( k ) and s( k ) respectively. The two by two
C  rotation part of P( k ), R( k ), is assumed to have the form
C
C     R( k ) = (  c( k )  s( k ) ).
C              ( -s( k )  c( k ) )
C
C  The matrix  U must be supplied in the n by n leading upper triangular
C  part of the array  A, and this is overwritten by the upper triangular
C  part of H.
C
C  When  SIDE = 'L' or 'l'  then a  row spike  is  generated  in  H  and
C  when  SIDE = 'R' or 'r'  then a  column spike is generated. For a row
C  spike the elements  h( k2, k )  and for a  column spike  the elements
C  h( k + 1, k1 ), k = k1, k1 + 1, ..., k2 - 1, are returned in  s( k ).
C
C  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
C  greater than n then an immediate return is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 13-January-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   AIJ, CTEMP, SPIKE, STEMP, TEMP
      INTEGER            I, J
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( ( MIN( N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.
     $   ( K2.GT.N ) )RETURN
      IF( ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' ) )THEN
C
C        Apply the plane rotations to columns n back to k2.
C
         DO 20 J = N, K2, -1
            TEMP = A( K2, J )
            DO 10 I = K2 - 1, K1, -1
               AIJ = A( I, J )
               A( I, J ) = S( I )*TEMP + C( I )*AIJ
               TEMP = C( I )*TEMP - S( I )*AIJ
   10       CONTINUE
            A( K2, J ) = TEMP
   20    CONTINUE
C
C        Form  the spike  and apply the rotations in columns  ( k2 - 1 )
C        back to k1.
C
         DO 40 J = K2 - 1, K1, -1
            SPIKE = -S( J )*A( J, J )
            A( J, J ) = C( J )*A( J, J )
            DO 30 I = J - 1, K1, -1
               AIJ = A( I, J )
               A( I, J ) = S( I )*SPIKE + C( I )*AIJ
               SPIKE = C( I )*SPIKE - S( I )*AIJ
   30       CONTINUE
            S( J ) = SPIKE
   40    CONTINUE
      ELSE IF( ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' ) )THEN
C
C        Apply the  plane rotations to columns  ( k1 + 1 ) up to k2  and
C        form the spike.
C
         DO 70 J = K1 + 1, K2
            CTEMP = C( J - 1 )
            STEMP = S( J - 1 )
            IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
               DO 50 I = 1, K1
                  TEMP = A( I, K1 )
                  A( I, K1 ) = STEMP*A( I, J ) + CTEMP*TEMP
                  A( I, J ) = CTEMP*A( I, J ) - STEMP*TEMP
   50          CONTINUE
               DO 60 I = K1 + 1, J - 1
                  SPIKE = S( I - 1 )
                  S( I - 1 ) = STEMP*A( I, J ) + CTEMP*SPIKE
                  A( I, J ) = CTEMP*A( I, J ) - STEMP*SPIKE
   60          CONTINUE
               S( J - 1 ) = STEMP*A( J, J )
               A( J, J ) = CTEMP*A( J, J )
            END IF
   70    CONTINUE
      END IF
C
      RETURN
C
C     End of F06QWF. ( SUTSRS )
C
      END
** END OF F06QWFTEXT
*UPTODATE F06QXFTEXT
      SUBROUTINE F06QXF( SIDE, PIVOT, DIRECT, M, N, K1, K2, C, S, A,
     $                   LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K1, K2, LDA, M, N
      CHARACTER*1        DIRECT, PIVOT, SIDE
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
C     ..
C
C  F06QXF  performs the transformation
C
C     A := P*A,   when   SIDE = 'L' or 'l'  (  Left-hand side )
C
C     A := A*P',  when   SIDE = 'R' or 'r'  ( Right-hand side )
C
C  where A is an m by n matrix and P is an orthogonal matrix, consisting
C  of a  sequence  of  plane  rotations,  applied  in  planes  k1 to k2,
C  determined by the parameters PIVOT and DIRECT as follows:
C
C     When  PIVOT  = 'V' or 'v'  ( Variable pivot )
C     and   DIRECT = 'F' or 'f'  ( Forward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C        where  P( k )  is a plane rotation matrix for the  ( k, k + 1 )
C        plane.
C
C     When  PIVOT  = 'V' or 'v'  ( Variable pivot )
C     and   DIRECT = 'B' or 'b'  ( Backward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C        where  P( k )  is a plane rotation matrix for the  ( k, k + 1 )
C        plane.
C
C     When  PIVOT  = 'T' or 't'  ( Top pivot )
C     and   DIRECT = 'F' or 'f'  ( Forward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k2 - 1 )*P( k2 - 2 )*...*P( k1 ),
C
C        where  P( k )  is a plane rotation matrix for the ( k1, k + 1 )
C        plane.
C
C     When  PIVOT  = 'T' or 't'  ( Top pivot )
C     and   DIRECT = 'B' or 'b'  ( Backward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C        where  P( k )  is a plane rotation matrix for the ( k1, k + 1 )
C        plane.
C
C     When  PIVOT  = 'B' or 'b'  ( Bottom pivot )
C     and   DIRECT = 'F' or 'f'  ( Forward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k2 - 1 )*P( k2 - 2 )*...*P( k1 ),
C
C        where  P( k )  is a  plane rotation  matrix  for the  ( k, k2 )
C        plane.
C
C     When  PIVOT  = 'B' or 'b'  ( Bottom pivot )
C     and   DIRECT = 'B' or 'b'  ( Backward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C        where  P( k )  is a  plane rotation  matrix  for the  ( k, k2 )
C        plane.
C
C  c( k ) and s( k )  must contain the  cosine and sine  that define the
C  matrix  P( k ).  The  two by two  plane rotation  part of the  matrix
C  P( k ), R( k ), is assumed to be of the form
C
C     R( k ) = (  c( k )  s( k ) ).
C              ( -s( k )  c( k ) )
C
C  If m, n or k1 are less than unity,  or k2 is not greater than k1,  or
C  SIDE = 'L' or 'l'  and  k2  is greater than  m, or  SIDE = 'R' or 'r'
C  and  k2  is greater than  n,  then an  immediate return  is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 20-November-1986.
C     Sven Hammarling and Mick Pont, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   AIJ, CTEMP, STEMP, TEMP
      INTEGER            I, J
      LOGICAL            LEFT, RIGHT
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      LEFT = ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' )
      RIGHT = ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' )
      IF( ( MIN( M, N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.
     $    ( ( LEFT ).AND.( K2.GT.M ) ).OR.
     $    ( ( RIGHT ).AND.( K2.GT.N ) ) )RETURN
      IF( LEFT )THEN
         IF( ( PIVOT.EQ.'V' ).OR.( PIVOT.EQ.'v' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 20 J = 1, N
                  AIJ = A( K1, J )
                  DO 10 I = K1, K2 - 1
                     TEMP = A( I + 1, J )
                     A( I, J ) = S( I )*TEMP + C( I )*AIJ
                     AIJ = C( I )*TEMP - S( I )*AIJ
   10             CONTINUE
                  A( K2, J ) = AIJ
   20          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 40 J = 1, N
                  AIJ = A( K2, J )
                  DO 30 I = K2 - 1, K1, -1
                     TEMP = A( I, J )
                     A( I + 1, J ) = C( I )*AIJ - S( I )*TEMP
                     AIJ = S( I )*AIJ + C( I )*TEMP
   30             CONTINUE
                  A( K1, J ) = AIJ
   40          CONTINUE
            END IF
         ELSE IF( ( PIVOT.EQ.'T' ).OR.( PIVOT.EQ.'t' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 60 J = 1, N
                  TEMP = A( K1, J )
                  DO 50 I = K1, K2 - 1
                     AIJ = A( I + 1, J )
                     A( I + 1, J ) = C( I )*AIJ - S( I )*TEMP
                     TEMP = S( I )*AIJ + C( I )*TEMP
   50             CONTINUE
                  A( K1, J ) = TEMP
   60          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 80 J = 1, N
                  TEMP = A( K1, J )
                  DO 70 I = K2 - 1, K1, -1
                     AIJ = A( I + 1, J )
                     A( I + 1, J ) = C( I )*AIJ - S( I )*TEMP
                     TEMP = S( I )*AIJ + C( I )*TEMP
   70             CONTINUE
                  A( K1, J ) = TEMP
   80          CONTINUE
            END IF
         ELSE IF( ( PIVOT.EQ.'B' ).OR.( PIVOT.EQ.'b' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 100 J = 1, N
                  TEMP = A( K2, J )
                  DO 90 I = K1, K2 - 1
                     AIJ = A( I, J )
                     A( I, J ) = S( I )*TEMP + C( I )*AIJ
                     TEMP = C( I )*TEMP - S( I )*AIJ
   90             CONTINUE
                  A( K2, J ) = TEMP
  100          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 120 J = 1, N
                  TEMP = A( K2, J )
                  DO 110 I = K2 - 1, K1, -1
                     AIJ = A( I, J )
                     A( I, J ) = S( I )*TEMP + C( I )*AIJ
                     TEMP = C( I )*TEMP - S( I )*AIJ
  110             CONTINUE
                  A( K2, J ) = TEMP
  120          CONTINUE
            END IF
         END IF
      ELSE IF( RIGHT )THEN
         IF( ( PIVOT.EQ.'V' ).OR.( PIVOT.EQ.'v' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 140 J = K1, K2 - 1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     DO 130 I = 1, M
                        TEMP = A( I, J + 1 )
                        A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  130                CONTINUE
                  END IF
  140          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 160 J = K2 - 1, K1, -1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     DO 150 I = M, 1, -1
                        TEMP = A( I, J + 1 )
                        A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  150                CONTINUE
                  END IF
  160          CONTINUE
            END IF
         ELSE IF( ( PIVOT.EQ.'T' ).OR.( PIVOT.EQ.'t' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 180 J = K1 + 1, K2
                  CTEMP = C( J - 1 )
                  STEMP = S( J - 1 )
                  IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
                     DO 170 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, K1 )
                        A( I, K1 ) = STEMP*TEMP + CTEMP*A( I, K1 )
  170                CONTINUE
                  END IF
  180          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 200 J = K2, K1 + 1, -1
                  CTEMP = C( J - 1 )
                  STEMP = S( J - 1 )
                  IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
                     DO 190 I = M, 1, -1
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, K1 )
                        A( I, K1 ) = STEMP*TEMP + CTEMP*A( I, K1 )
  190                CONTINUE
                  END IF
  200          CONTINUE
            END IF
         ELSE IF( ( PIVOT.EQ.'B' ).OR.( PIVOT.EQ.'b' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 220 J = K1, K2 - 1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     DO 210 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, K2 ) + CTEMP*TEMP
                        A( I, K2 ) = CTEMP*A( I, K2 ) - STEMP*TEMP
  210                CONTINUE
                  END IF
  220          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 240 J = K2 - 1, K1, -1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     DO 230 I = M, 1, -1
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, K2 ) + CTEMP*TEMP
                        A( I, K2 ) = CTEMP*A( I, K2 ) - STEMP*TEMP
  230                CONTINUE
                  END IF
  240          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06QXF. ( SGESRC )
C
      END
** END OF F06QXFTEXT
      SUBROUTINE F06QZF( HESS, N, K1, K2, C, S, A, LDA )
*     .. Scalar Arguments ..
      CHARACTER*1        HESS
      INTEGER            K1, K2, LDA, N
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
*     ..
*
*  F06QZF  either applies a  given sequence  of  plane rotations  to the
*  right of the n by n reverse lower triangular matrix T, to transform T
*  to a  reverse lower Hessenberg matrix  H, or restores a reverse lower
*  Hessenberg matrix H to reverse lower triangular form T, by applying a
*  sequence of plane rotations from the right.
*
*  The rotations are applied  in planes k1 up to k2.
*
*  When   HESS = 'C' or 'c',   ( Create ),  then   the   reverse   lower
*  Hessenberg matrix, H, is formed as
*
*     H = T*P',
*
*  where P is an orthogonal matrix of the form
*
*     P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
*
*  P( k ) being a plane rotation matrix for the  ( k, k + 1 ) plane. The
*  cosine and sine that define P( k ), k = k1, k1 + 1, ..., k2 - 1, must
*  be  supplied  in  c( k )  and  s( k )  respectively.  The  two by two
*  rotation part of P( k ), R( k ), is assumed to have the form
*
*     R( k ) = (  c( k )  s( k ) ).
*              ( -s( k )  c( k ) )
*
*  The matrix  T must be supplied in the n by n reverse lower triangular
*  part  of the array  A,  and this is overwritten by the  reverse lower
*  triangular part of  H.
*
*  The super-diagonal elements of  H, h( n - k, k ), are returned in the
*  elements s( k ),  k = k1, k1 + 1, ..., k2 - 1.
*
*  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
*  greater than n then an immediate return is effected.
*
*  When   HESS = 'R' or 'r',   ( Remove ),  then   the   reverse   lower
*  Hessenberg matrix  H  is  assumed  to  have  non-zero  super-diagonal
*  elements  in  positions  h( n - k, k ),  k = k1, k1 + 1, ..., k2 - 1,
*  only and  h( n - k, k ) must be supplied in  s( k ). H is restored to
*  the reverse lower triangular matrix T as
*
*     T = H*P',
*
*  where P is an orthogonal matrix of the form
*
*     P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
*
*  P( k ) being a plane rotation for the  ( k, k + 1 ) plane. The cosine
*  and  sine  that  define  P( k )  are  returned  in  c( k ) and s( k )
*  respectively.  The  two by two  rotation part of  P( k ),  R( k ), is
*  of the form
*
*     R( k ) = (  c( k )  s( k ) ).
*              ( -s( k )  c( k ) )
*
*  The reverse lower triangular part of the matrix H must be supplied in
*  the  n by n  reverse  lower  triangular  part  of  A,   and  this  is
*  overwritten by the reverse triangular matrix T.
*
*  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
*  greater than n then an immediate return is effected.
*
*  When   n = 7, k1 = 2 and k2 = 5   then  T  and  H  are  of  the  form
*
*     T = ( 0  0  0  0  0  0  X ),   H = ( 0  0  0  0  0  0  X ).
*         ( 0  0  0  0  0  X  X )        ( 0  0  0  0  X  X  X )
*         ( 0  0  0  0  X  X  X )        ( 0  0  0  X  X  X  X )
*         ( 0  0  0  X  X  X  X )        ( 0  0  X  X  X  X  X )
*         ( 0  0  X  X  X  X  X )        ( 0  X  X  X  X  X  X )
*         ( 0  X  X  X  X  X  X )        ( 0  X  X  X  X  X  X )
*         ( X  X  X  X  X  X  X )        ( X  X  X  X  X  X  X )
*
*
*  This routine  is  principally intended  for use  with the  non-linear
*  optimization routines such as E04UCF, in order to help vectorization.
*  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
*
*  -- Written on 10-May-1988.
*     Sven Hammarling, Nag Central Office.
*
*
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     .. External Subroutines ..
      EXTERNAL           F06BAF
*     .. Local Scalars ..
      DOUBLE PRECISION   CTEMP, STEMP, SUPH, TEMP
      INTEGER            I, J
*     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Executable Statements ..
      IF( ( MIN( N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.( K2.GT.N ) )
     $   RETURN
      IF( ( HESS.EQ.'C' ).OR.( HESS.EQ.'c' ) )THEN
*
*        Apply  the  plane rotations  to  columns  k1  up to  ( k2 - 1 )
*        and  form   the  additional  super-diagonal  elements,  storing
*        h( n - j, j ) in s( j ).
*
         DO 20, J = K1, K2 - 1
            IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
               STEMP             = S( J )
               CTEMP             = C( J )
               S( J )            = STEMP*A( N - J, J + 1 )
               A( N - J, J + 1 ) = CTEMP*A( N - J, J + 1 )
               DO 10, I = N - J + 1, N
                  TEMP          = A( I, J + 1 )
                  A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                  A( I, J )     = STEMP*TEMP + CTEMP*A( I, J )
   10          CONTINUE
            END IF
   20    CONTINUE
      ELSE IF( ( HESS.EQ.'R' ).OR.( HESS.EQ.'r' ) )THEN
*
*        Restore  H to reverse lower triangular form by annihilating the
*        super-diagonal elements of  H.  The  jth rotation  is chosen so
*        that
*
*          ( h( n - j, n - j ) ) := (  c  s )*( h( n - j, n - j     ) ),
*          (         0         )    ( -s  c ) ( h( n - j, n - j - 1 ) )
*
*        which can be expressed as
*
*           ( 0  h( n - j, n - j ) ) :=
*
*               ( h( n - j, n - j - 1 )  h( n - j, n - j ) )*(  c  s ).
*                                                            ( -s  c )
*
*        Thus we return  c( j ) = c  and  s( j ) = -s  to make the plane
*        rotation matrix look like
*
*           R( j ) = (  c( j )  s( j ) ).
*                    ( -s( j )  c( j ) )
*
         DO 40, J = K2 - 1, K1, -1
            SUPH   =  S( J )
            CALL F06BAF( A( N - J, J + 1 ), SUPH, CTEMP, STEMP )
            STEMP  = -STEMP
            S( J ) =  STEMP
            C( J ) =  CTEMP
            IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
               DO 30, I = N - J + 1, N
                  TEMP          = A( I, J + 1 )
                  A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                  A( I, J )     = STEMP*TEMP + CTEMP*A( I, J )
   30          CONTINUE
            END IF
   40    CONTINUE
      END IF
*
      RETURN
*
*     End of F06QZF.
*
      END
** END OF F06QZFTEXT
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     file  lssubs.f
*
*     lssol    lsadd    lsadds   lsbnds   lschol   lscore   lscrsh
*     lsdel    lsdflt   lsfeas   lsfile   lsgetp   lsgset   lskey
*     lsloc    lsmove   lsmuls   lsoptn   lsprt    lssetx
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lssol ( mm, n,
     $                   nclin, ldA, ldR,
     $                   A, bl, bu, cvec,
     $                   istate, kx, x, R, b,
     $                   inform, iter, obj, clamda,
     $                   iw, leniw, w, lenw )

      implicit           double precision(a-h,o-z)
      integer            leniw, lenw
      integer            istate(n+nclin), kx(n)
      integer            iw(leniw)
      double precision   bl(n+nclin), bu(n+nclin), A(ldA,*)
      double precision   clamda(n+nclin), cvec(*)
      double precision   R(ldR,*), x(n), b(*)
      double precision   w(lenw)

*     ==================================================================
*     lssol  solves problems of the form
*
*           Minimize               f(x)
*              x
*                                 (  x )
*           subject to    bl  .le.(    ).le.  bu,
*                                 ( Ax )
*
*     where  '  denotes the transpose of a column vector,  x  denotes
*     the n-vector of parameters and  f(x) is one of the following
*     functions...
*
*    FP = None                        (find a feasible point).
*    LP = c'x
*    QP1=       1/2 x'Rx               R  n times n, symmetric pos. def.
*    QP2= c'x + 1/2 x'Rx               .  .   ..        ..       ..  ..
*    QP3=       1/2 x'R'Rx             R  m times n, upper triangular.
*    QP4= c'x + 1/2 x'R'Rx             .  .   ..  .   ..      ...
*    LS1=       1/2 (b - Rx)'(b - Rx)  R  m times n, rectangular.
*    LS2= c'x + 1/2 (b - Rx)'(b - Rx)  .  .   ..  .     ...
*    LS3=       1/2 (b - Rx)'(b - Rx)  R  m times n, upper triangular.
*    LS4= c'x + 1/2 (b - Rx)'(b - Rx)  .  .   ..  .   ..      ...
*
*     The matrix  R  is entered as the two-dimensional array  R  (of
*     row dimension  ldR).  If  ldR = 0,  R  is not accessed.
*
*     The vector  c  is entered in the one-dimensional array  cvec.
*
*     nclin  is the number of general linear constraints (rows of  A).
*     (nclin may be zero.)
*
*     The first  n  components of  bl  and   bu  are lower and upper
*      bounds on the variables.  The next  nclin  components are
*     lower and upper bounds on the general linear constraints.
*
*     The matrix  A  of coefficients in the general linear constraints
*     is entered as the two-dimensional array  A  (of dimension
*     ldA by n).  If nclin = 0, A is not accessed.
*
*     The vector  x  must contain an initial estimate of the solution,
*     and will contain the computed solution on output.
*
*
*     Complete documentation for LSSOL is contained in Report SOL 86-1,
*     Users Guide for LSSOL (Version 1.0), by P.E. Gill,
*     S.J. Hammarling, W. Murray, M.A. Saunders and M.H. Wright,
*     Department of Operations Research, Stanford University, Stanford, 
*     California 94305.
*
*     Systems Optimization Laboratory, Stanford University.
*     Version 1.00 Dated  30-Jan-1986.
*     Version 1.01 Dated  30-Jun-1986.   Level-2 BLAS added
*     Version 1.02 Dated  13-May-1988.   Level-2 matrix routines added.
*     Version 1.03 Dated  19-Jun-1989.   Some obscure bugs fixed.      
*     Version 1.04 Dated  26-Aug-1991.   nrank bug fixed.      
*     Version 1.05 Dated  20-Sep-1992.   Output modified. 
*                         20-Oct-1992    Summary file included.
*                         
*     Copyright  1984/1993  Stanford University.
*
*     This material may be reproduced by or for the U.S. Government 
*     pursuant to the copyright license under DAR Clause 7-104.9(a)
*     (1979 Mar).
*
*     This material is based upon work partially supported by the 
*     National Science Foundation under Grants MCS-7926009 and 
*     ECS-8312142; the Department of Energy Contract AM03-76SF00326,
*     PA No. DE-AT03-76ER72018; the Army Research Office Contract
*     DAA29-84-K-0156; and the Office of Naval Research Grant
*     N00014-75-C-0267.
*     ==================================================================
      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2
      common    /sol3cm/ lennam, ldT, ncolT, ldQ
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9
      common    /sol5cm/ Asize, dTmax, dTmin

      parameter         (lenls = 20)
      common    /sol1ls/ locls(lenls)

      logical            lsdbg
      parameter         (ldbg = 5)
      common    /lsdebg/ ilsdbg(ldbg), lsdbg
*     +Include lsparm+++++++++++++++++++++++++++++++++++++++++++++++++++
      parameter         (mxparm = 30)
      integer            iprmls(mxparm), ipsvls
      double precision   rprmls(mxparm), rpsvls

      common    /lspar1/ ipsvls(mxparm),
     $                   idbgls, iPrnt , iSumry, itmax1, itmax2, lcrash,
     $	                 ldbgls, lprob , msgls , nn    , nnclin, nprob , 
     $                   ipadls(18)

      common    /lspar2/ rpsvls(mxparm),
     $                   bigbnd, bigdx , bndlow, bndupp, tolact, tolfea,
     $                   tolrnk, rpadls(23)

      equivalence       (iprmls(1), idbgls), (rprmls(1), bigbnd)

      save      /lspar1/, /lspar2/
*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      equivalence   (msgls , msglvl), (idbgls, idbg), (ldbgls, msgdbg)

      intrinsic          max, min

*     Local variables.

      logical            cold  , factrz, linObj, named , rowerr,
     $                   unitQ , vertex
      character*2        prbtyp
      character*8        names(1)
      parameter         (zero   =0.0d+0, point3 =3.3d-1, point8 =0.8d+0)
      parameter         (point9 =0.9d+0, one    =1.0d+0, hundrd =1.0d+2)

      character*40       title
      data               title
     $                 / 'SOL/LSSOL  ---  Version 1.05-2  April 93' /
*                         123456789|123456789|123456789|123456789|
*
*     Set the machine-dependent constants.

      call lsoptn('Problem Type Qp2')
      call mchpar()

      epsmch = wmach( 3)
      rteps  = wmach( 4)
      nout   = wmach(11)

      epspt3 = epsmch**point3
      epspt5 = rteps
      epspt8 = epsmch**point8
      epspt9 = epsmch**point9

      named  = .false.

      inform = 0
      iter   = 0

C-->  condmx = max( one/epspt5, hundrd )
C-->  condmx = max( one/epspt3, hundrd )
      condmx = max( one/epspt5, hundrd )

      nctotl = n + nclin

*     Set the default values of the parameters.

      call lsdflt( mm, n, nclin, title )

*     Set all parameters determined by the problem type.

      if      (lprob .eq. 1 ) then
         prbtyp    = 'FP'
         m      = 0
         linObj = .false.
         factrz = .true.
      else if (lprob .eq. 2 ) then
         prbtyp    = 'LP'
         m      = 0
         linObj = .true.
         factrz = .true.
      else if (lprob .eq. 3 ) then
         prbtyp    = 'QP'
         m      = mm
         linObj = .false.
         factrz = .true.
      else if (lprob .eq. 4 ) then
         prbtyp    = 'QP'
         m      = mm
         linObj = .true.
         factrz = .true.
      else if (lprob .eq. 5 ) then
         prbtyp    = 'QP'
         m      = mm
         linObj = .false.
         factrz = .false.
      else if (lprob .eq. 6 ) then
         prbtyp    = 'QP'
         m      = mm
         linObj = .true.
         factrz = .false.
      else if (lprob .eq. 7 ) then
         prbtyp    = 'LS'
         m      = mm
         linObj = .false.
         factrz = .true.
      else if (lprob .eq. 8 ) then
         prbtyp    = 'LS'
         m      = mm
         linObj = .true.
         factrz = .true.
      else if (lprob .eq. 9 ) then
         prbtyp    = 'LS'
         m      = mm
         linObj = .false.
         factrz = .false.
      else if (lprob .eq. 10) then
         prbtyp    = 'LS'
         m      = mm
         linObj = .true.
         factrz = .false.
      end if

*     Assign the dimensions of arrays in the parameter list of lscore.
*     Economies of storage are possible if the minimum number of active
*     constraints and the minimum number of fixed variables are known in
*     advance.  The expert user should alter minact and minfxd
*     accordingly.
*     If a linear program is being solved and the matrix of general
*     constraints is fat,  i.e.,  nclin .lt. n,  a non-zero value is
*     known for minfxd.  Note that in this case, vertex must be
*     set  .true..

      minact = 0
      minfxd = 0

      vertex = .false.
      if (      (prbtyp .eq. 'LP'  .or.  prbtyp .eq. 'FP')
     $    .and.  nclin  .lt. n   ) then
         minfxd = n - nclin - 1
         vertex = .true.
      end if

      mxfree = n - minfxd
      maxact = max( 1, min( n, nclin ) )
      maxnZ  = n - ( minfxd + minact )

      if (nclin .eq. 0) then
         ldQ    = 1
         ldT    = 1
         ncolt  = 1
         vertex = .false.
      else
         ldQ    = max( 1, mxfree )
         ldT    = max( maxnZ, maxact )
         ncolt  = mxfree
      end if

*     Allocate certain arrays that are not done in lsloc.

      litotl = 0

      lAx    = 1
      lwtotl = lAx + nclin  - 1

*     Allocate remaining work arrays.

      call lsloc ( lprob, n, nclin, litotl, lwtotl )

      cold  = lcrash .eq. 0

*     Check input parameters and storage limits.

      ncnln  = 0
      lennam = 1

      call cmchk ( nerror, msglvl, lcrash, (.not.factrz),
     $             leniw, lenw, litotl, lwtotl,
     $             n, nclin, ncnln,
     $             istate, kx, named, names,
     $             bigbnd, bl, bu, x )

      if (nerror .gt. 0) then
         inform = 6
         go to 800
      end if

      lkactv = locls( 1)

      lanorm = locls( 2)
      lpx    = locls( 4)
      lres   = locls( 5)
      lres0  = locls( 6)
      lgq    = locls( 8)
      lcq    = locls( 9)
      lrlam  = locls(10)
      lT     = locls(11)
      lQ     = locls(12)
      lwtinf = locls(13)
      lwrk   = locls(14)
      lfeatl = locls(15)

      if (tolfea .gt. zero)
     $   call dload ( n+nclin, (tolfea), w(lfeatl), 1 )

      ianrmj = lanorm
      do 200, j = 1, nclin
         w(ianrmj) = dnrm2 ( n, A(j,1), ldA )
         ianrmj    = ianrmj + 1
  200 continue
      if (nclin .gt. 0)
     $   call dcond ( nclin, w(lanorm), 1, Asize, amin )

      call dcond ( nctotl, w(lfeatl), 1, feamax, feamin )
      call dcopy ( nctotl, w(lfeatl), 1, w(lwtinf), 1 )
      call dscal ( nctotl, (one/feamin), w(lwtinf), 1 )

      ssq1   = zero

      if (factrz) then
*        ===============================================================
*        Factorize R using QR or Cholesky.  kx must be initialized.
*        ===============================================================
         do 210, i = 1, n
            kx(i) = i
  210    continue

         if      (prbtyp .eq. 'LP'  .or.  prbtyp .eq. 'FP') then
            nrank = 0
         else if (prbtyp .eq. 'QP') then
*           ------------------------------------------------------------
*           Compute the Cholesky factorization of R.  The Hessian is
*           m by m and resides in the upper left-hand corner of R.
*           ------------------------------------------------------------
            do 220, j = m+1, n
               call dload ( m, (zero), R(1,j), 1 )
  220       continue

            call lschol( ldR, m, nrank, tolrnk, kx, R, info )

            if (nrank .gt. 0)
     $         call dload ( nrank, (zero), w(lres0), 1 )
         else if (prbtyp .eq. 'LS') then
*           ------------------------------------------------------------
*           Compute the orthogonal factorization PRQ = ( U ),  where P
*                                                      ( 0 )
*           is an orthogonal matrix and Q is a permutation matrix.
*           Overwrite R with the upper-triangle U.  The orthogonal
*           matrix P is applied to the residual and discarded.  The
*           permutation is stored in the array KX.  Once U has been
*           computed we need only work with vectors of length N within
*           lscore.  However, it is necessary to store the sum of
*           squares of the terms  b(nrank+1),...,b(m),  where B = Pr.
*           ------------------------------------------------------------
            call dgeqrp( 'Column iterchanges', m, n, R, ldR,
     $                   w(lwrk), iw(lkactv), w(lgq), info )

            lj  = lkactv
            do 230, j = 1, n
               jmax = iw(lj)
               if (jmax .gt. j) then
                  jsave    = kx(jmax)
                  kx(jmax) = kx(j)
                  kx(j)    = jsave
               end if
               lj = lj + 1
  230       continue

            call dgeapq( 'Transpose', 'Separate', m, min( n,m-1 ), 
     $                   R, ldR, w(lwrk), 1, b, m, w(lgq), info )

            rownrm = dnrm2 ( n, R(1,1), ldR )
            if (          rownrm  .le.        tolrnk
     $          .or.  abs(R(1,1)) .le. rownrm*tolrnk) then
               nrank = 0
            else
               nrank = idrank( min(n, m), R, ldR+1, tolrnk )
            end if

            if (m .gt. nrank) ssq1 = dnrm2 ( m-nrank, b(nrank+1), 1 )

            if (nrank .gt. 0)
     $         call dcopy ( nrank, b, 1, w(lres0), 1 )
         end if
      else
*        ===============================================================
*        R is input as an upper-triangular matrix with m rows.
*        ===============================================================
         nrank = m
         if (nrank .gt. 0) then
            if      (prbtyp .eq. 'QP') then
               call dload ( nrank, (zero), w(lres0), 1 )
            else if (prbtyp .eq. 'LS') then
               call dcopy ( nrank, b, 1, w(lres0), 1 )
            end if
         end if
      end if

      if (       msglvl .gt. 0     .and.  nrank  .lt. n
     $    .and.  prbtyp .ne. 'LP'  .and.  prbtyp .ne. 'FP') then
*         if (iPrint .gt. 0) write(iPrint, 9000) nrank
      end if
*     ------------------------------------------------------------------
*     Find an initial working set.
*     ------------------------------------------------------------------
      call lscrsh( cold, vertex,
     $             nclin, nctotl, nactiv, nartif,
     $             nfree, n, ldA,
     $             istate, iw(lkactv),
     $             bigbnd, tolact,
     $             A, w(lAx), bl, bu, x, w(lgq), w(lwrk) )

*     ------------------------------------------------------------------
*     Compute the TQ factorization of the constraints while keeping R in
*     upper-triangular form.  Transformations associated with Q are
*     applied to cq.  Transformations associated with P are applied to
*     res0.  If some simple bounds are in the working set,  kx is
*     re-ordered so that the free variables come first.
*     ------------------------------------------------------------------
*     First, add the bounds. To save a bit of work, cq is not loaded
*     until after kx has been re-ordered.

      ngq   = 0
      nres  = 0
      if (nrank .gt. 0) nres = 1
      unitQ = .true.

      call lsbnds( unitQ,
     $             inform, nZ, nfree, nrank, nres, ngq,
     $             n, ldQ, ldA, ldR, ldT,
     $             istate, kx, condmx,
     $             A, R, w(lT), w(lres0), w(lcq), w(lQ),
     $             w(lwrk), w(lpx), w(lrlam) )

      if (linObj) then

*        Install the transformed linear term in cq.
*        cmqmul applies the permutations in kx to cvec.

         ngq = 1
         call dcopy ( n, cvec, 1, w(lcq), 1 )
         call cmqmul( 6, n, nZ, nfree, ldQ, unitQ,
     $                kx, w(lcq), w(lQ), w(lwrk) )
      end if

      if (nactiv .gt. 0) then
         nact1  = nactiv
         nactiv = 0

         call lsadds( unitQ, vertex,
     $                inform, 1, nact1, nactiv, nartif, nZ, nfree,
     $                nrank, nrejtd, nres, ngq,
     $                n, ldQ, ldA, ldR, ldT,
     $                istate, iw(lkactv), kx, condmx,
     $                A, R, w(lT), w(lres0), w(lcq), w(lQ),
     $                w(lwrk), w(lpx), w(lrlam) )
      end if

*     ------------------------------------------------------------------
*     Move the initial  x  onto the constraints in the working set.
*     Compute the transformed residual vector  Pr = Pb - RQ'x.
*     ------------------------------------------------------------------
      call lssetx( linObj, rowerr, unitQ,
     $             nclin, nactiv, nfree, nrank, nZ,
     $             n, nctotl, ldQ, ldA, ldR, ldT,
     $             istate, iw(lkactv), kx,
     $             jmax, errmax, ctx, xnorm,
     $             A, w(lAx), bl, bu, w(lcq), w(lres), w(lres0),
     $             w(lfeatl), R, w(lT), x, w(lQ), w(lpx), w(lwrk) )

      if (rowerr) then
*         if (iPrint .gt. 0) write(iPrint, 9010)
         inform = 3
         numinf = 1
         suminf = errmax
         go to 800
      end if

      jinf = 0

      call lscore( prbtyp, named, names, linObj, unitQ,
     $             inform, iter, jinf, nclin, nctotl,
     $             nactiv, nfree, nrank, nZ, nZr,
     $             n, ldA, ldR,
     $             istate, iw(lkactv), kx,
     $             ctx, obj, ssq1,
     $             suminf, numinf, xnorm,
     $             bl, bu, A, clamda, w(lAx),
     $             w(lfeatl), R, x, w )

      obj    = obj    + ctx
      if (prbtyp .eq. 'LS'  .and.  nrank .gt. 0)
     $   call dcopy ( nrank, w(lres), 1, b, 1 )

*     ==================================================================
*     Print messages if required.
*     ==================================================================
  800 if (msglvl .gt.   0) then
         if (inform .eq.   0) then
            if (prbtyp .eq. 'FP') then
*               if (iPrint .gt. 0) write(iPrint, 2001)
*               if (iSumm  .gt. 0) write(iSumm , 2001)
            else
*               if (iPrint .gt. 0) write(iPrint, 2002) prbtyp
*               if (iSumm  .gt. 0) write(iSumm , 2002) prbtyp
            end if
         end if
         if (iPrint .gt. 0) then
*            if (inform .eq.   1) write(iPrint, 2010) prbtyp
*            if (inform .eq.   2) write(iPrint, 2020) prbtyp
*            if (inform .eq.   3) write(iPrint, 2030)
*            if (inform .eq.   4) write(iPrint, 2040)
*            if (inform .eq.   5) write(iPrint, 2050)
*            if (inform .eq.   6) write(iPrint, 2060) nerror
         end if

*         if (iSumm  .gt. 0) then
*            if (inform .eq.   1) write(iSumm , 2010) prbtyp
*            if (inform .eq.   2) write(iSumm , 2020) prbtyp
*            if (inform .eq.   3) write(iSumm , 2030)
*            if (inform .eq.   4) write(iSumm , 2040)
*            if (inform .eq.   5) write(iSumm , 2050)
*            if (inform .eq.   6) write(iSumm , 2060) nerror
*         end if
            
         if (inform .lt.   6) then
            if      (numinf .eq. 0) then
                if (prbtyp .ne. 'FP') then
*                   if (iPrint .gt. 0) write(iPrint, 3000) prbtyp, obj
*                   if (iSumm  .gt. 0) write(iSumm , 3000) prbtyp, obj
                end if
            else if (inform .eq. 3) then
*               if (iPrint .gt. 0) write(iPrint, 3010) suminf
*               if (iSumm  .gt. 0) write(iSumm , 3010) suminf
            else
*
*               if (iPrint .gt. 0) write(iPrint, 3020) suminf
*               if (iSumm  .gt. 0) write(iSumm , 3020) suminf
            end if
            if (numinf .gt. 0) obj = suminf
         end if
      end if

*     Recover the optional parameters set by the user.

      call icopy ( mxparm, ipsvls, 1, iprmls, 1 )
      call dcopy ( mxparm, rpsvls, 1, rprmls, 1 )

      return

 2001 format(/ ' Exit LSSOL - Feasible point found.     ')
 2002 format(/ ' Exit LSSOL - Optimal ', A2, ' solution.')
 2010 format(/ ' Exit LSSOL - Weak ',    A2, ' solution.')
 2020 format(/ ' Exit LSSOL - ', A2,         ' solution is unbounded.' )
 2030 format(/ ' Exit LSSOL - Cannot satisfy the linear constraints. ' )
 2040 format(/ ' Exit LSSOL - Too many iterations.')
 2050 format(/ ' Exit LSSOL - Too many iterations without changing X.' )
 2060 format(/ ' Exit LSSOL - ', I10, ' errors found in the input',
     $         ' parameters.  Problem abandoned.'         )
 3000 format(/ ' Final ', A2, ' objective value =', G16.7 )
 3010 format(/ ' Minimum sum of infeasibilities =', G16.7 )
 3020 format(/ ' Final sum of infeasibilities =',   G16.7 )

 9000 format(/ ' Rank of the objective function data matrix = ', I5 )
 9010 format(  ' XXX  Cannot satisfy the working-set constraints to',
     $         ' the accuracy requested.')

*     end of lssol
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lsadd ( unitQ,
     $                   inform, ifix, iadd, jadd,
     $                   nactiv, nZ, nfree, nrank, nres, ngq,
     $                   n, ldA, ldQ, ldR, ldT,
     $                   kx, condmx,
     $                   A, R, T, res, gqm, Q,
     $                   w, c, s )

      implicit           double precision(a-h,o-z)
      logical            unitQ
      integer            kx(n)
      double precision   A(ldA,*), R(ldR,*), T(ldT,*),
     $                   res(n,*), gqm(n,*), Q(ldQ,*)
      double precision   w(n), c(n), s(n)

*     ==================================================================
*     lsadd   updates the factorization,  A(free) * (Z Y) = (0 T),  when
*     a constraint is added to the working set.  If  nrank .gt. 0, the
*     factorization  ( R ) = PCQ  is also updated,  where  C  is the
*                    ( 0 )
*     least squares matrix,  R  is upper-triangular,  and  P  is an
*     orthogonal matrix.  The matrices  C  and  P  are not stored.
*
*     There are three separate cases to consider (although each case
*     shares code with another)...
*
*     (1) A free variable becomes fixed on one of its bounds when there
*         are already some general constraints in the working set.
*
*     (2) A free variable becomes fixed on one of its bounds when there
*         are only bound constraints in the working set.
*
*     (3) A general constraint (corresponding to row  iadd  of  A) is
*         added to the working set.
*
*     In cases (1) and (2), we assume that  kx(ifix) = jadd.
*     In all cases,  jadd  is the index of the constraint being added.
*
*     If there are no general constraints in the working set,  the
*     matrix  Q = (Z Y)  is the identity and will not be touched.
*
*     If  nres .gt. 0,  the row transformations are applied to the rows
*     of the  (n by nres)  matrix  res.
*     If  ngq .gt. 0,  the column transformations are applied to the
*     columns of the  (ngq by n)  matrix  gqm'.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 31-October--1984.
*     Level-2 matrix routines added 25-Apr-1988.
*     This version of  lsadd  dated 14-Sep-92.
*     ==================================================================
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9
      common    /sol5cm/ Asize, dTmax, dTmin

      logical            lsdbg
      parameter         (ldbg = 5)
      common    /lsdebg/ ilsdbg(ldbg), lsdbg

      logical            bound , overfl
      external           ddiv  , dnrm2
      intrinsic          max   , min
      parameter         (zero = 0.0d+0, one = 1.0d+0)

*     if the condition estimator of the updated factors is greater than
*     condbd,  a warning message is printed.

      condbd = one / epspt9

      overfl = .false.
      bound  = jadd .le. n

      if (bound) then
*        ===============================================================
*        A simple bound has entered the working set.  iadd  is not used.
*        ===============================================================
*         if (lsdbg  .and.  ilsdbg(1) .gt. 0)
*     $      write(iPrint, 1010) nactiv, nZ, nfree, ifix, jadd, unitQ
         nanew = nactiv

         if (unitQ) then

*           Q  is not stored, but kx defines an ordering of the columns
*           of the identity matrix that implicitly define  Q.
*           Define the sequence of pairwise interchanges p that moves
*           the newly-fixed variable to position nfree.
*           Reorder kx accordingly.

            do 100, i = 1, nfree-1
               if (i .ge. ifix) then
                  w (i) = i + 1
                  kx(i) = kx(i+1)
               else
                  w(i) = i
               end if
  100       continue
         else
*           ------------------------------------------------------------
*           Q  is stored explicitly.
*           ------------------------------------------------------------
*           Set  w = the  (ifix)-th  row of  Q.
*           Move the  (nfree)-th  row of  Q  to position  ifix.

            call dcopy ( nfree, Q(ifix,1), ldQ, w, 1 )
            if (ifix .lt. nfree) then
               call dcopy ( nfree, Q(nfree,1), ldQ, Q(ifix,1), ldQ )
               kx(ifix) = kx(nfree)
            end if
         end if
         kx(nfree) = jadd
      else
*        ===============================================================
*        A general constraint has entered the working set.
*        ifix  is not used.
*        ===============================================================
*         if (lsdbg  .and.  ilsdbg(1) .gt. 0)
*     $      write(iPrint, 1020) nactiv, nZ, nfree, iadd, jadd, unitQ

         nanew  = nactiv + 1

*        Transform the incoming row of  A  by  Q'.  Use c as workspace.

         call dcopy ( n, A(iadd,1), ldA, w, 1 )
         call cmqmul( 8, n, nZ, nfree, ldQ, unitQ, kx, w, Q, c )

*        Check that the incoming row is not dependent upon those
*        already in the working set.

         dTnew  = dnrm2 ( nZ, w, 1 )
         if (nactiv .eq. 0) then

*           This is the only general constraint in the working set.

            cond   = ddiv  ( Asize, dTnew, overfl )
            tdTmax = dTnew
            tdTmin = dTnew
         else

*           There are already some general constraints in the working
*           set. Update the estimate of the condition number.

            tdTmax = max( dTnew, dTmax )
            tdTmin = min( dTnew, dTmin )
            cond   = ddiv  ( tdTmax, tdTmin, overfl )
         end if

         if (cond .gt. condmx  .or.  overfl) go to 900

         if (unitQ) then

*           First general constraint added.  Set  Q = I.

            call f06qhf( 'General', nfree, nfree, zero, one, Q, ldQ )
            unitQ  = .false.
         end if
      end if

      if (bound) then
         npiv  = nfree
      else
         npiv  = nZ
      end if

      nT = min( nrank, npiv )

      if (unitQ) then
*        ---------------------------------------------------------------
*        Q (i.e., Q) is not stored explicitly.
*        Apply the sequence of pairwise interchanges P that moves the
*        newly-fixed variable to position nfree.
*        ---------------------------------------------------------------
         if (ngq .gt. 0)
     $      call f06qkf( 'Left', 'Transpose', nfree-1, w, ngq, gqm, n )
            
         if (nrank .gt. 0) then

*           Apply the pairwise interchanges to the triangular part of R.
*           The subdiagonal elements generated by this process are
*           stored in  s(1), s(2), ..., s(nt-1).

            call f06qnf( 'Right', n, ifix, nT, s, R, ldR )

            if (nt .lt. npiv) then

*              R is upper trapezoidal.  Apply the interchanges in
*              columns  nT  thru  npiv.

               do 200, i = ifix, nT-1
                  w(i) = i
  200          continue

               call f06qkf( 'Right', 'Normal', nfree-1, w, nT, R, ldR )
            end if
            
*           Eliminate the subdiagonal elements of R with a left-hand
*           sweep of rotations P2 in planes (1,2), (2,3), ...,(nt-1,nt).
*           Apply P2 to res.

            call f06qrf( 'Left ', n, ifix, nT, c, s, R, ldR )
            if (nres .gt. 0) 
     $         call f06qxf( 'Left', 'Variable', 'Forwards', nT, nres,
     $                      ifix, nT, c, s, res, n )
         end if
      else
*        ---------------------------------------------------------------
*        Full matrix Q.  Define a sweep of plane rotations P such that
*                           Pw = beta*e(npiv).
*        The rotations are applied in the planes (1,2), (2,3), ...,
*        (npiv-1,npiv).  The rotations must be applied to Q, R, T
*        and GQM'.
*        ---------------------------------------------------------------
         call f06fqf( 'Varble', 'Forwrds', npiv-1, w(npiv), w, 1, c, s )

         if (bound  .and.  nactiv .gt. 0) then

            call dcopy ( nactiv, s(nZ), 1, w(nZ), 1 )
         
            s(       nZ  ) = s(nZ)*T(nactiv,nZ+1)
            T(nactiv,nZ+1) = c(nZ)*T(nactiv,nZ+1)
         
            call f06qzf( 'Create', nactiv, 1, nactiv, c(nZ+1), s(nZ+1),
     $                   T(1,nZ+1), ldT )
            call dcopy ( nactiv, s(nZ), 1, T(nactiv,nZ), ldT-1 )
         
            call dcopy ( nactiv, w(nZ), 1, s(nZ), 1 )
         end if

         if (ngq .gt. 0)
     $      call f06qxf( 'Left ', 'Variable', 'Forwards', npiv , ngq,
     $                   1, npiv, c, s, gqm, n )
         call f06qxf( 'Right', 'Variable', 'Forwards', nfree, nfree,
     $                1, npiv, c, s, Q, ldQ )

         if (nrank .gt. 0) then

*           Apply the rotations to the triangular part of R.
*           The subdiagonal elements generated by this process are
*           stored in  s(1),  s(2), ..., s(nt-1).

            nT = min( nrank, npiv )
            call f06qvf( 'Right', n, 1, nT, c, s, R, ldR )

            if (nt .lt. npiv) then

*              R is upper trapezoidal.  Pretend R is (nt x n) and
*              apply the rotations in columns  nT  thru  npiv.

               call f06qxf( 'Right', 'Variable', 'Forwards', nT, n,
     $                      nT, npiv, c, s, R, ldR )
            end if

*           Eliminate the subdiagonal elements of R with a left-hand
*           sweep of rotations P2 in planes (1,2), (2,3), ...,(nt-1,nt).
*           Apply P2 to res.

            call f06qrf( 'Left ', n, 1, nT, c, s, R, ldR )
            if (nres .gt. 0)
     $         call f06qxf( 'Left', 'Variable', 'Forwards', nT, nres,
     $                      1, nT, c, s, res, n )
         end if

         if (bound) then

*           The last row and column of Q has been transformed to plus
*           or minus the unit vector e(nfree).  We can reconstitute the
*           columns of GQM and R corresponding to the new fixed variable.

            if (w(nfree) .lt. zero) then
               nf = min( nrank, nfree )
               if (nf  .gt. 0) call dscal ( nf , -one,   R(1,nfree), 1 )
               if (ngq .gt. 0) call dscal ( ngq, -one, gqm(nfree,1), n )
            end if

*           ------------------------------------------------------------
*           The diagonals of T have been altered.  Recompute the
*           largest and smallest values.
*           ------------------------------------------------------------
            if (nactiv .gt. 0) then
               call dcond( nactiv, T(nactiv,nZ), ldT-1, tdTmax, tdTmin )
               cond   = ddiv  ( tdTmax, tdTmin, overfl )
            end if
         else
*           ------------------------------------------------------------
*           General constraint.  Install the new row of T.
*           ------------------------------------------------------------
            call dcopy ( nanew, w(nZ), 1, T(nanew,nZ), ldT )
         end if
      end if

*     ==================================================================
*     Prepare to exit.  Check the magnitude of the condition estimator.
*     ==================================================================
  900 if (nanew .gt. 0) then
         if (cond .lt. condmx  .and.  .not. overfl) then

*           The factorization has been successfully updated.

            inform = 0
            dTmax  = tdTmax
            dTmin  = tdTmin
            if (cond .ge. condbd) then
*               if (iPrint .gt. 0) write(iPrint, 2000) jadd
            end if
         else

*           The proposed working set appears to be linearly dependent.

            inform = 1
            if (lsdbg  .and.  ilsdbg(1) .gt. 0) then
*               write( iPrint, 3000 )
               if (bound) then
*                  write(iPrint, 3010) Asize, dTmax, dTmin
               else
                  if (nactiv .gt. 0) then
*                     write(iPrint, 3020) Asize, dTmax, dTmin, dTnew
                  else
*                     write(iPrint, 3030) Asize, dTnew
                  end if
               end if
            end if
         end if
      end if

      return

 1010 format(/ ' //lsadd //  Simple bound added.'
     $       / ' //lsadd //  nactiv    nZ nfree  ifix  jadd unitQ'
     $       / ' //lsadd //  ', 5i6, l6 )
 1020 format(/ ' //lsadd //  General constraint added.           '
     $       / ' //lsadd //  nactiv    nZ nfree  iadd  jadd unitQ'
     $       / ' //lsadd //  ', 5i6, l6 )
 2000 format(/ ' XXX  Serious ill-conditioning in the working set',
     $         ' after adding constraint ',  i5
     $       / ' XXX  Overflow may occur in subsequent iterations.'//)
 3000 format(/ ' //lsadd //  Dependent constraint rejected.' )
 3010 format(/ ' //lsadd //     Asize     dTmax     dTmin        '
     $       / ' //lsadd //', 1p, 3e10.2 )
 3020 format(/ ' //lsadd //     Asize     dTmax     dTmin     dTnew'
     $       / ' //lsadd //', 1p, 4e10.2 )
 3030 format(/ ' //lsadd //     Asize     dTnew'
     $       / ' //lsadd //', 1p, 2e10.2 )

*     end of lsadd
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lsadds( unitQ, vertex,
     $                   inform, k1, k2, nactiv, nartif, nZ, nfree,
     $                   nrank, nrejtd, nres, ngq,
     $                   n, ldQ, ldA, ldR, ldT,
     $                   istate, kactiv, kx, condmx,
     $                   A, R, T, res, gqm, Q,
     $                   w, c, s )

      implicit           double precision(a-h,o-z)
      logical            unitQ, vertex
      integer            istate(*), kactiv(n), kx(n)
      double precision   condmx
      double precision   A(ldA,*), R(ldR,*),
     $                   T(ldT,*), res(n,*), gqm(n,*), Q(ldQ,*)
      double precision   w(n), c(n), s(n)

*     ==================================================================
*     lsadds  includes general constraints k1 thru k2 as new rows of
*     the TQ factorization stored in T, Q.  If nrank is nonZero, the
*     changes in Q are reflected in nrank by n triangular factor R such
*     that
*                         C  =  P ( R ) Q,
*                                 ( 0 )
*     where  P  is orthogonal.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written  October-31-1984.
*     This version of lsadds dated  16-May-1988.
*     ==================================================================

      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/
      common    /sol5cm/ Asize, dTmax, dTmin

      external           dnrm2
      parameter        ( zero = 0.0d+0, one = 1.0d+0 )

      rtmax  = wmach(8)

*     Estimate the condition number of the constraints that are not
*     to be refactorized.

      if (nactiv .eq. 0) then
         dTmax = zero
         dTmin = one
      else
         call dcond ( nactiv, T(nactiv,nZ+1), ldT-1, dTmax, dTmin )
      end if

      do 200, k = k1, k2
         iadd = kactiv(k)
         jadd = n + iadd
         if (nactiv .lt. nfree) then

            call lsadd ( unitQ,
     $                   inform, ifix, iadd, jadd,
     $                   nactiv, nZ, nfree, nrank, nres, ngq,
     $                   n, ldA, ldQ, ldR, ldT,
     $                   kx, condmx,
     $                   A, R, T, res, gqm, Q,
     $                   w, c, s )

            if (inform .eq. 0) then
               nactiv = nactiv + 1
               nZ     = nZ     - 1
            else
               istate(jadd) =   0
               kactiv(k)    = - kactiv(k)
            end if
         end if
  200 continue

      if (nactiv .lt. k2) then

*        Some of the constraints were classed as dependent and not
*        included in the factorization.  Re-order the part of kactiv
*        that holds the indices of the general constraints in the
*        working set.  Move accepted indices to the front and shift
*        rejected indices (with negative values) to the end.

         l      = k1 - 1
         do 300, k = k1, k2
            i         = kactiv(k)
            if (i .ge. 0) then
               l      = l + 1
               if (l .ne. k) then
                  iswap     = kactiv(l)
                  kactiv(l) = i
                  kactiv(k) = iswap
               end if
            end if
  300    continue

*        If a vertex is required, add some temporary bounds.
*        We must accept the resulting condition number of the working
*        set.

         if (vertex) then
            cndmax = rtmax
            nZadd  = nZ
            do 320, iartif = 1, nZadd
               if (unitQ) then
                  ifix = nfree
                  jadd = kx(ifix)
               else
                  rowmax = zero
                  do 310, i = 1, nfree
                     rnorm = dnrm2 ( nZ, Q(i,1), ldQ )
                     if (rowmax .lt. rnorm) then
                        rowmax = rnorm
                        ifix   = i
                     end if
  310             continue
                  jadd = kx(ifix)

                  call lsadd ( unitQ,
     $                         inform, ifix, iadd, jadd,
     $                         nactiv, nZ, nfree, nrank, nres, ngq,
     $                         n, ldA, ldQ, ldR, ldT,
     $                         kx, cndmax,
     $                         A, R, T, res, gqm, Q,
     $                         w, c, s  )
               end if
               nfree  = nfree  - 1
               nZ     = nZ     - 1
               nartif = nartif + 1
               istate(jadd) = 4
  320       continue
         end if
      end if

      nrejtd = k2 - nactiv

*     end of lsadds
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lsbnds( unitQ,
     $                   inform, nZ, nfree, nrank, nres, ngq,
     $                   n, ldQ, ldA, ldR, ldT,
     $                   istate, kx, condmx,
     $                   A, R, T, res, gqm, Q,
     $                   w, c, s )

      implicit           double precision(a-h,o-z)
      logical            unitQ
      integer            istate(*), kx(n)
      double precision   condmx
      double precision   A(ldA,*), R(ldR,*),
     $                   T(ldT,*), res(n,*), gqm(n,*), Q(ldQ,*)
      double precision   w(n), c(n), s(n)

*     ==================================================================
*     lsbnds updates the factor R as kx is reordered to reflect the
*     status of the bound constraints given by istate.  kx is reordered
*     so that the fixed variables come last.  One of two alternative
*     methods are used to reorder kx. One method needs fewer accesses 
*     to kx, the other gives a matrix  Rz  with more rows and columns.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written  30-December-1985.
*     This version of lsbnds dated 13-May-88.
*     ==================================================================

      nfixed = n - nfree

      if (nrank .lt. n  .and.  nrank .gt. 0) then
*        ---------------------------------------------------------------
*        R is specified but singular.  Try and keep the dimension of Rz
*        as large as possible.
*        ---------------------------------------------------------------
         nactv = 0
         nfree = n
         nZ    = n

         j     = n
*+       while (j .gt. 0  .and.  n-nfree .lt. nfixed) do
  100    if    (j .gt. 0  .and.  n-nfree .lt. nfixed) then
            if (istate(j) .gt. 0) then
               jadd = j
               do 110, ifix = nfree, 1, -1
                  if (kx(ifix) .eq. jadd) go to 120
  110          continue

*              Add bound jadd.

  120          call lsadd ( unitQ,
     $                      inform, ifix, iadd, jadd,
     $                      nactv, nZ, nfree, nrank, nres, ngq,
     $                      n, ldA, ldQ, ldR, ldT,
     $                      kx, condmx,
     $                      A, R, T, res, gqm, Q,
     $                      w, c, s )

               nfree = nfree - 1
               nZ    = nZ    - 1
            end if
            j = j - 1
            go to 100
*+       end while
         end if
      else     
*        ---------------------------------------------------------------
*        R is of full rank,  or is not specified.
*        ---------------------------------------------------------------
         if (nfixed .gt. 0) then

*           Order kx so that the free variables come first.

            lstart = nfree + 1
            do 250, k = 1, nfree
               j = kx(k)
               if (istate(j) .gt. 0) then
                  do 220, l = lstart, n
                     j2 = kx(l)
                     if (istate(j2) .eq. 0) go to 230
  220             continue

  230             kx(k)  = j2
                  kx(l)  = j
                  lstart = l + 1

                  if (nrank .gt. 0)
     $               call cmrswp( n, nres, nrank, ldR, k, l,
     $                            R, res, c, s )
               end if
  250       continue

         end if
         nZ = nfree
      end if

*     end of lsbnds
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lschol( ldh, n, nrank, tolrnk, kx, h, inform )

      implicit           double precision (a-h,o-z)
      integer            kx(*)
      double precision   h(ldh,*)

*     ==================================================================
*     LSCHOL  forms the Cholesky factorization of the positive
*     semi-definite matrix H such that
*                   PHP'  =  R'R
*     where  P  is a permutation matrix and  R  is upper triangular.
*     The permutation P is chosen to maximize the diagonal of R at each
*     stage.  Only the diagonal and super-diagonal elements of H are
*     used.
*
*     Output:
*
*         inform = 0   the factorization was computed successfully,
*                      with the Cholesky factor written in the upper
*                      triangular part of H and P stored in kx.
*                  1   the matrix H was indefinite.
*
*     Original version of lschol dated  2-February-1981.
*     Level 2 Blas added 29-June-1986.
*     This version of lschol dated  26-Jun-1989. 
*     ==================================================================
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2
      intrinsic          abs   , max   , sqrt
      external           idamax
      parameter        ( zero = 0.0d+0, one = 1.0d+0 )

      inform = 0
      nrank  = 0

*     Main loop for computing rows of  R.

      do 200, j = 1, n

*        Find maximum available diagonal.    

         kmax = j - 1 + idamax( n-j+1, h(j,j), ldh+1 )
         dmax = h(kmax,kmax)

         if (dmax .le. tolrnk*abs(h(1,1))) go to 300

*        Perform a symmetric interchange if necessary.

         if (kmax .ne. j) then
            k        = kx(kmax)
            kx(kmax) = kx(j)
            kx(j)    = k

            call dswap ( kmax-j, h(j+1,kmax), 1, h(j,j+1 ), ldh )
            call dswap ( j     , h(1  ,j   ), 1, h(1,kmax), 1   )
            call dswap ( n-kmax+1, h(kmax,kmax), ldh,
     $                             h(j,kmax)   , ldh )

         end if

*        Set the diagonal of  R.

         d      = sqrt( dmax )
         h(j,j) = d
         nrank  = nrank + 1

         if (j .lt. n) then

*           Set the super-diagonal elements of this row of R and update
*           the elements of the block that is yet to be factorized.
                                                          
            call dscal ( n-j,   (one/d), h(j  ,j+1), ldh )
            call dsyr  ( 'u', n-j, -one, h(j  ,j+1), ldh,
     $                                   h(j+1,j+1), ldh )
         end if

  200 continue
*     ------------------------------------------------------------------
*     Check for the semi-definite case.
*     ------------------------------------------------------------------
  300 if (nrank .lt. n) then

*        Find the largest element in the unfactorized block.

         supmax = zero
         do 310, i = j, n-1
            k      = i + idamax( n-i, h(i,i+1), ldh )
            supmax = max( supmax, abs(h(i,k)) )
  310    continue

         if (supmax .gt. tolrnk*abs(h(1,1))) then
*            if (iPrint .gt. 0) write(iPrint, 1000) dmax, supmax
            inform = 1
         end if
      end if

      return

 1000 format(' XXX  Hessian appears to be indefinite.'
     $      /' XXX  Maximum diagonal and off-diagonal ignored',
     $             ' in the Cholesky factorization:', 1p, 2e22.14 )

*     end of lschol
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lscore( prbtyp, named, names, linObj, unitQ,
     $                   inform, iter, jinf, nclin, nctotl,
     $                   nactiv, nfree, nrank, nZ, nZr,
     $                   n, ldA, ldR,
     $                   istate, kactiv, kx,
     $                   ctx, ssq, ssq1, suminf, numinf, xnorm,
     $                   bl, bu, A, clamda, Ax,
     $                   featol, R, x, w )

      implicit           double precision(a-h,o-z)
      character*2        prbtyp
      character*8        names(*)
      integer            istate(nctotl), kactiv(n), kx(n)
      double precision   bl(nctotl), bu(nctotl), A(ldA,*),
     $                   clamda(nctotl), Ax(*),
     $                   featol(nctotl), R(ldR,*), x(n)
      double precision   w(*)
      logical            named, linObj, unitQ

*     ==================================================================
*     lscore  is a subroutine for linearly constrained linear-least
*     squares.  On entry, it is assumed that an initial working set of
*     linear constraints and bounds is available.
*     The arrays istate, kactiv and kx will have been set accordingly
*     and the arrays t and Q will contain the TQ factorization of
*     the matrix whose rows are the gradients of the active linear
*     constraints with the columns corresponding to the active bounds
*     removed.  the TQ factorization of the resulting (nactiv by nfree)
*     matrix is  A(free)*Q = (0 T),  where Q is (nfree by nfree) and t
*     is reverse-triangular.
*
*     Values of istate(j) for the linear constraints.......
*
*     istate(j)
*     ---------
*          0    constraint j is not in the working set.
*          1    constraint j is in the working set at its lower bound.
*          2    constraint j is in the working set at its upper bound.
*          3    constraint j is in the working set as an equality.
*
*     Constraint j may be violated by as much as featol(j).
*
*     Systems Optimization Laboratory, Stanford University.
*     This version of  lscore  dated  14-Sep-92.
*
*     Copyright  1984/1993  Stanford University.
*
*     This material may be reproduced by or for the U.S. Government 
*     pursuant to the copyright license under DAR Clause 7-104.9(a)
*     (1979 Mar).
*
*     This material is based upon work partially supported by the 
*     National Science Foundation under Grants MCS-7926009 and 
*     ECS-8312142; the Department of Energy Contract AM03-76SF00326,
*     PA No. DE-AT03-76ER72018; the Army Research Office Contract
*     DAA29-84-K-0156; and the Office of Naval Research Grant
*     N00014-75-C-0267.
*     ==================================================================
      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2
      common    /sol3cm/ lennam, ldT   , ncolt , ldQ
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9
      common    /sol5cm/ Asize, dTmax  , dTmin

      integer            locls
      parameter         (lenls = 20)
      common    /sol1ls/ locls(lenls)

      logical            cmdbg, lsdbg
      parameter         (ldbg = 5)
      common    /lsdebg/ ilsdbg(ldbg), lsdbg
      common    /cmdebg/ icmdbg(ldbg), cmdbg
*     +Include lsparm+++++++++++++++++++++++++++++++++++++++++++++++++++
      parameter         (mxparm = 30)
      integer            iprmls(mxparm), ipsvls
      double precision   rprmls(mxparm), rpsvls

      common    /lspar1/ ipsvls(mxparm),
     $                   idbgls, iPrnt , iSumry, itmax1, itmax2, lcrash,
     $	                 ldbgls, lprob , msgls , nn    , nnclin, nprob , 
     $                   ipadls(18)

      common    /lspar2/ rpsvls(mxparm),
     $                   bigbnd, bigdx , bndlow, bndupp, tolact, tolfea,
     $                   tolrnk, rpadls(23)

      equivalence       (iprmls(1), idbgls), (rprmls(1), bigbnd)

      save      /lspar1/, /lspar2/
*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      equivalence   (msgls , msglvl), (idbgls, idbg), (ldbgls, msgdbg)

      external           ddiv  , dnrm2
      intrinsic          abs   , max   , sqrt
      logical            convrg, cyclin, error , firstv, hitcon,
     $                   hitlow, needfg, overfl, prnt  , rowerr
      logical            singlr, stall , statpt, unbndd, uncon , unitgZ,
     $                   weak
      parameter        ( zero   =0.0d+0, half   =0.5d+0, one   =1.0d+0 )
      parameter        ( mrefn  =1     , mstall =50                    )

*     Specify the machine-dependent parameters.

      flmax  = wmach(7)

      lanorm = locls( 2)
      lAp    = locls( 3)
      lpx    = locls( 4)
      lres   = locls( 5)
      lres0  = locls( 6)
      lhZ    = locls( 7)
      lgq    = locls( 8)
      lcq    = locls( 9)
      lrlam  = locls(10)
      lT     = locls(11)
      lQ     = locls(12)
      lwtinf = locls(13)
      lwrk   = locls(14)

*     Set up the adresses of the contiguous arrays  ( res0, res )
*     and  ( gq, cq ).

      nres   = 0
      if (nrank .gt. 0) nres = 2
      ngq    = 1
      if (linObj) ngq = 2

*     Initialize.

      irefn  =   0
      iter   =   0
      itmax  = itmax1

      jadd   =   0
      jdel   =   0
      nphase =   1
      nstall =   0
      numinf = - 1
      nZr    =   0

      alfa   = zero
      condmx = flmax
      dRzmax = one
      dRzmin = one
      ssq    = zero

      cyclin = .false.
      error  = .false.
      firstv = .false.
      prnt   = .true.
      needfg = .true.
      stall  = .true.
      uncon  = .false.
      unbndd = .false.

*     If debug output is required,  print nothing until iteration IDBG.

      msgsvd = msglvl
      if (idbg .gt. 0  .and.  idbg .le. itmax) then
         msglvl = 0
      end if

*======================== start of the main loop =======================
*
*      cyclin = false
*      unbndd = false
*      error  = false
*      k      = 0
*
*      repeat
*            repeat
*                  compute Z'g,  print details of this iteration
*                  stat pt = (Z'g .eq. 0)
*                  if (not stat pt) then
*                     error =  k .ge. itmax
*                     if (not error) then
*                        compute p, alfa
*                        error = unbndd  or  cyclin
*                        if (not error) then
*                           k = k + 1
*                           x = x + alfa p
*                           if (feasible) update Z'g
*                           if necessary, add a constraint
*                        end if
*                     end if
*                  end if
*            until  stat pt  or  error
*
*            compute lam1, lam2, smllst
*            optmul =  smllst .gt. 0
*            if ( not (optmul .or. error) ) then
*                  delete an artificial or regular constraint
*            end if
*      until optmul  or  error
*
*=======================================================================

*     repeat
*        repeat
  100       if (needfg) then
               if (nrank .gt. 0) then
                  resnrm = dnrm2 ( nrank, w(lres), 1 )
                  ssq    = half*(ssq1**2 + resnrm**2 )
               end if

               if (numinf .ne. 0) then

*                 Compute the transformed gradient of either the sum of
*                 of infeasibilities or the objective.  Initialize
*                 singlr and unitgZ.

                  call lsgset( prbtyp, linObj, singlr, unitgZ, unitQ,
     $                         n, nclin, nfree,
     $                         ldA, ldQ, ldR, nrank, nZ, nZr,
     $                         istate, kx,
     $                         bigbnd, tolrnk, numinf, suminf,
     $                         bl, bu, A, w(lres), featol,
     $                         w(lgq), w(lcq), R, x, w(lwtinf),
     $                         w(lQ), w(lwrk) )

                  if (prbtyp .ne. 'FP'  .and.  numinf .eq. 0
     $                                  .and.  nphase .eq. 1) then
                     itmax  = iter + itmax2
                     nphase = 2
                  end if
               end if
            end if

            gznorm = zero
            if (nZ  .gt. 0 ) gznorm = dnrm2 ( nZ, w(lgq), 1 )

            if (nZr .eq. nZ) then
               gZrnrm = gznorm
            else
               gZrnrm = zero
               if (nZr .gt. 0) gZrnrm = dnrm2 ( nZr, w(lgq), 1 )
            end if

            gfnorm = gznorm
            if (nfree .gt. 0  .and.  nactiv .gt. 0)
     $         gfnorm = dnrm2 ( nfree, w(lgq), 1 )

*           ------------------------------------------------------------
*           Print the details of this iteration.
*           ------------------------------------------------------------
*           Define small quantities that reflect the size of x, R and
*           the constraints in the working set.  If feasible,  estimate
*           the rank and condition number of Rz.
*           Note that nZr .le. nrank + 1.

            if (nZr .eq. 0) then
               singlr = .false.
            else
               if (numinf .gt. 0  .or.  nZr .gt. nrank) then
                  absrzz = zero
                  singlr = .true.
               else
                  call dcond ( nZr, R, ldR+1, dRzmax, dRzmin )
                  absrzz = abs( R(nZr,nZr) )
                  rownrm = dnrm2 ( n, R(1,1), ldR )
                  singlr =       absrzz      .le. dRzmax*tolrnk 
     $                     .or.  rownrm      .le.        tolrnk 
     $                     .or.  abs(R(1,1)) .le. rownrm*tolrnk
               end if

*               if (lsdbg  .and.  ilsdbg(1) .gt. 0)
*     $            write(iPrint, 9100) singlr, absrzz, dRzmax, dRzmin
            end if

            condRz = ddiv  ( dRzmax, dRzmin, overfl )
            condT  = one
            if (nactiv .gt. 0)
     $         condT  = ddiv  ( dTmax , dTmin , overfl )

            if (prnt) then
               call lsprt ( prbtyp, isdel, iter, jadd, jdel,
     $                      msglvl, nactiv, nfree, n, nclin,
     $                      nrank, ldR, ldT, nZ, nZr, istate,
     $                      alfa, condRz, condT, gfnorm, gZrnrm,
     $                      numinf, suminf, ctx, ssq,
     $                      Ax, R, w(lT), x, w(lwrk) )
               jdel  = 0
               jadd  = 0
               alfa  = zero
            end if

            if (numinf .gt. 0) then
               dinky  = zero
            else
               objsiz = one  + abs( ssq + ctx )
               wssize = zero
               if (nactiv .gt. 0) wssize = dTmax
               dinky  = epspt8 * max( wssize, objsiz, gfnorm )
               if (uncon) then
                  unitgZ = gZrnrm .le. dinky
               end if
            end if

*            if (lsdbg  .and.  ilsdbg(1) .gt. 0)
*     $         write(iPrint, 9000) unitgZ, irefn, gZrnrm, dinky

*           If the projected gradient  Z'g  is small and Rz is of full
*           rank, X is a minimum on the working set.  An additional
*           refinement step is allowed to take care of an inaccurate
*           value of dinky.

            statpt = .not. singlr  .and.  gZrnrm .le. dinky
     $                             .or.   irefn  .gt. mrefn

            if (.not. statpt) then
*              ---------------------------------------------------------
*              Compute a search direction.
*              ---------------------------------------------------------
               prnt  = .true.

               error = iter .ge. itmax
               if (.not. error) then

                  irefn = irefn + 1
                  iter  = iter  + 1

                  if (iter .eq. idbg  .and.  iPrint .gt. 0) then
                     lsdbg  = .true.
                     cmdbg  =  lsdbg
                     msglvl =  msgsvd
                  end if

                  call lsgetp( linObj, singlr, unitgZ, unitQ,
     $                         n, nclin, nfree,
     $                         ldA, ldQ, ldR, nrank, numinf, nZr,
     $                         kx, ctp, pnorm,
     $                         A, w(lAp), w(lres), w(lhZ), w(lpx),
     $                         w(lgq), w(lcq), R, w(lQ), w(lwrk) )

*                 ------------------------------------------------------
*                 Find the constraint we bump into along p.
*                 Update x and Ax if the step alfa is nonZero.
*                 ------------------------------------------------------
*                 alfhit is initialized to bigalf.  If it remains
*                 that way after the call to cmalf, it will be
*                 regarded as infinite.

                  bigalf = ddiv  ( bigdx, pnorm, overfl )

                  call cmalf ( firstv, hitlow,
     $                         istate, inform, jadd, n, nctotl, numinf,
     $                         alfhit, palfa, atphit, 
     $                         bigalf, bigbnd, pnorm,
     $                         w(lAnorm), w(lAp), Ax, bl, bu, 
     $                         featol, w(lpx), x )

*                 If  Rz  is nonsingular,  alfa = 1.0  will be the
*                 step to the least-squares minimizer on the
*                 current subspace. If the unit step does not violate
*                 the nearest constraint by more than featol,  the
*                 constraint is not added to the working set.

                  hitcon = singlr  .or.  palfa  .le. one
                  uncon  = .not. hitcon

                  if (hitcon) then
                     alfa = alfhit
                  else
                     jadd   = 0
                     alfa   = one
                  end if

*                 Check for an unbounded solution or negligible step.

                  unbndd =  alfa .ge. bigalf
                  stall  = abs( alfa*pnorm ) .le. epspt9*xnorm
                  if (stall) then
                     nstall = nstall + 1
                     cyclin = nstall .gt. mstall
                  else
                     nstall = 0
                  end if

                  error = unbndd  .or.  cyclin
                  if (.not.  error) then
*                    ---------------------------------------------------
*                    Set x = x + alfa*p.  Update Ax, gq, res and ctx.
*                    ---------------------------------------------------
                     if (alfa .ne. zero)
     $                  call lsmove( hitcon, hitlow, linObj, unitgZ,
     $                               nclin, nrank, nZr,
     $                               n, ldR, jadd, numinf,
     $                               alfa, ctp, ctx, xnorm,
     $                               w(lAp), Ax, bl, bu, w(lgq),
     $                               w(lhZ), w(lpx), w(lres),
     $                               R, x, w(lwrk) )

                     if (hitcon) then
*                       ------------------------------------------------
*                       Add a constraint to the working set.
*                       Update the TQ factors of the working set.
*                       Use p as temporary work space.
*                       ------------------------------------------------
*                       Update  istate.

                        if (bl(jadd) .eq. bu(jadd)) then
                           istate(jadd) = 3
                        else if (hitlow) then
                           istate(jadd) = 1
                        else
                           istate(jadd) = 2
                        end if
                        iadd = jadd - n
                        if (jadd .le. n) then

                           do 510, ifix = 1, nfree
                              if (kx(ifix) .eq. jadd) go to 520
  510                      continue
  520                   end if

                        call lsadd ( unitQ,
     $                               inform, ifix, iadd, jadd,
     $                               nactiv, nZ, nfree, nrank, nres,ngq,
     $                               n, ldA, ldQ, ldR, ldT,
     $                               kx, condmx,
     $                               A, R, w(lT), w(lres),w(lgq),w(lQ),
     $                               w(lwrk), w(lrlam), w(lpx) )

                        nZr    = nZr - 1
                        nZ     = nZ  - 1

                        if (jadd .le. n) then

*                          A simple bound has been added.

                           nfree  = nfree  - 1
                        else

*                          A general constraint has been added.

                           nactiv = nactiv + 1
                           kactiv(nactiv) = iadd
                        end if
                        irefn  = 0
                     end if
*                    ---------------------------------------------------
*                    Check the feasibility of constraints with non-
*                    negative istate values.  If some violations have
*                    occurred.  Refine the current x and set inform so
*                    that feasibility is checked in lsgset.
*                    ---------------------------------------------------
                     call lsfeas( n, nclin, istate,
     $                            bigbnd, cnorm, err1, jmax1, nviol,
     $                            Ax, bl, bu, featol, x, w(lwrk) )

                     if (err1 .gt. featol(jmax1)) then
                        call lssetx( linObj, rowerr, unitQ,
     $                               nclin, nactiv, nfree, nrank, nZ,
     $                               n, nctotl, ldQ, ldA, ldR, ldT,
     $                               istate, kactiv, kx,
     $                               jmax1, err2, ctx, xnorm,
     $                               A, Ax, bl, bu, w(lcq),
     $                               w(lres), w(lres0), featol, R,
     $                               w(lT), x, w(lQ), w(lpx), w(lwrk) )

*                        if (lsdbg  .and.  ilsdbg(1) .gt. 0)
*     $                     write(iPrint, 2100) err1, err2
                        if (rowerr) then
*                           if (iPrint .gt. 0) write(iPrint, 2200)
                        end if
                        uncon  =   .false.
                        irefn  =   0
                        numinf = - 1
                     end if
                     needfg = alfa .ne. zero
                  end if
               end if
            end if
*        until      statpt  .or.  error
         if (.not. (statpt  .or.  error) ) go to 100

*        ===============================================================
*        Try and find the index jdel of a constraint to drop from
*        the working set.
*        ===============================================================
         jdel   = 0

         if (numinf .eq. 0  .and.  prbtyp .eq. 'FP') then
            if (n .gt. nZ)
     $         call dload ( n-nZ, (zero), w(lrlam), 1 )
            jtiny  = 0
            jsmlst = 0
            jbigst = 0
         else

            call lsmuls( prbtyp,
     $                   msglvl, n, nactiv, nfree,
     $                   ldA, ldT, numinf, nZ, nZr,
     $                   istate, kactiv, kx, dinky,
     $                   jsmlst, ksmlst, jinf, jtiny,
     $                   jbigst, kbigst, trulam,
     $                   A, w(lanorm), w(lgq), w(lrlam),
     $                   w(lT), w(lwtinf) )
         end if

         if (.not. error) then
            if (     jsmlst .gt. 0) then

*              LSMULS found a regular constraint with multiplier less
*              than (-dinky).

               jdel   = jsmlst
               kdel   = ksmlst
               isdel  = istate(jdel)
               istate(jdel) = 0

            else if (jsmlst .lt. 0) then

               jdel   = jsmlst

            else if (numinf .gt. 0  .and.  jbigst .gt. 0) then

*              No feasible point exists for the constraints but the
*              sum of the constraint violations may be reduced by
*              moving off constraints with multipliers greater than 1.

               jdel   = jbigst
               kdel   = kbigst
               isdel  = istate(jdel)
               if (trulam .le. zero) is = - 1
               if (trulam .gt. zero) is = - 2
               istate(jdel) = is
               firstv = .true.
               numinf = numinf + 1
            end if

            if      (jdel .ne. 0  .and.  singlr) then

*              Cannot delete a constraint when Rz is singular.
*              Probably a weak minimum.

               jdel = 0
            else if (jdel .ne. 0               ) then

*              Constraint jdel has been deleted.
*              Update the matrix factorizations.

               call lsdel ( unitQ,
     $                      n, nactiv, nfree, nres, ngq, nZ, nZr,
     $                      ldA, ldQ, ldR, ldT, nrank,
     $                      jdel, kdel, kactiv, kx,
     $                      A, w(lres), R, w(lT), w(lgq), w(lQ),
     $                      w(lwrk), w(lpx) )
            end if
         end if

         irefn  =  0
         convrg =  jdel .eq. 0
         prnt   = .false.
         uncon  = .false.
         needfg = .false.

*     until       convrg  .or.  error
      if (.not.  (convrg  .or.  error)) go to 100

*  .........................End of main loop............................
      weak = jtiny .gt. 0  .or.  singlr

      if (error) then
         if (unbndd) then
            inform = 2
            if (numinf .gt. 0) inform = 3
         else if (iter .ge. itmax) then
            inform = 4
         else if (cyclin) then
            inform = 5
         end if
      else if (convrg) then
         inform = 0
         if (numinf .gt. 0) then
            inform = 3
         else if (prbtyp .ne. 'FP'  .and.  weak) then
            inform = 1
         end if
      end if

*     ------------------------------------------------------------------
*     Set   clamda.  Print the full solution.
*     ------------------------------------------------------------------
      msglvl = msgsvd
*      if (msglvl .gt. 0  .and.  iPrint .gt. 0)
*     $     write(iPrint, 2000) prbtyp, iter, inform

      call cmprt ( msglvl, nfree, ldA,
     $             n, nclin, nctotl, bigbnd,
     $             named, names,
     $             nactiv, istate, kactiv, kx,
     $             A, bl, bu, x, clamda, w(lrlam), x )

      return

 2000 format(/ ' Exit from ', a2, ' problem after ', i4, ' iterations.',
     $         '  inform =', i3 )
 2100 format(  ' XXX  Iterative refinement.  Maximum errors before and',
     $         ' after refinement are ',  1p, 2e14.2 )
 2200 format(  ' XXX  Warning.  Cannot satisfy the constraints to the',
     $         ' accuracy requested.')
 9000 format(/ ' //lscore//  unitgZ irefn     gZrnrm      dinky'
     $       / ' //lscore//  ', l6, i6, 1p, 2e11.2 )
 9100 format(/ ' //lscore//  singlr    abs(Rzz)      dRzmax      dRzmin'
     $       / ' //lscore//  ', l6,     1p, 3e12.4 )

*     end of lscore
      end                         

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lscrsh( cold, vertex,
     $                   nclin, nctotl, nactiv, nartif,
     $                   nfree, n, ldA,
     $                   istate, kactiv,
     $                   bigbnd, tolact,
     $                   A, Ax, bl, bu, x, wx, work )

      implicit           double precision(a-h,o-z)
      logical            cold, vertex
      integer            istate(nctotl), kactiv(n)
      double precision   A(ldA,*), Ax(*), bl(nctotl), bu(nctotl),
     $                   x(n), wx(n), work(n)

*     ==================================================================
*     lscrsh  computes the quantities  istate (optionally), kactiv,
*     nactiv, nZ and nfree  associated with the working set at x.
*     The computation depends upon the value of the input parameter
*     cold,  as follows...
*
*     cold = true.  An initial working set will be selected. First,
*                   nearly-satisfied or violated bounds are added.
*                   Next,  general linear constraints are added that
*                   have small residuals.
*
*     cold = false. The quantities kactiv, nactiv, nZ and nfree are
*                   computed from istate,  specified by the user.
*
*     Values of istate(j)....
*
*        - 2         - 1         0           1          2         3
*     a'x lt bl   a'x gt bu   a'x free   a'x = bl   a'x = bu   bl = bu
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 31-October-1984.
*     This version of lscrsh dated 14-May-1992.
*     ==================================================================
      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2

      logical            lsdbg
      parameter         (ldbg = 5)
      common    /lsdebg/ ilsdbg(ldbg), lsdbg

      external           ddot
      intrinsic          abs, min
      parameter        ( zero = 0.0d+0, one = 1.0d+0 )

      flmax  =   wmach(7)
      biglow = - bigbnd
      bigupp =   bigbnd

*     ------------------------------------------------------------------
*     Move the variables inside their bounds.
*     ------------------------------------------------------------------
      do 10, j = 1, n
         b1    = bl(j)
         b2    = bu(j)

         if (b1 .gt. biglow) then
            if (x(j) .lt. b1) x(j) = b1
         end if
      
         if (b2 .lt. bigupp) then
            if (x(j) .gt. b2) x(j) = b2
         end if
   10 continue

      call dcopy ( n, x, 1, wx, 1 )

      if (lsdbg) then
*         if (ilsdbg(1) .gt. 0)
*     $      write(iPrint, 1000) cold, nclin, nctotl
*         if (ilsdbg(2) .gt. 0)
*     $      write(iPrint, 1100) (wx(j), j = 1, n)
      end if

      nfixed = 0
      nactiv = 0
      nartif = 0

*     If a cold start is being made, initialize  istate.
*     If  bl(j) = bu(j),  set  istate(j)=3  for all variables and linear
*     constraints.

      if (cold) then
         do 100, j = 1, nctotl
            istate(j) = 0
            if (bl(j) .eq. bu(j)) istate(j) = 3
  100    continue
      else
         do 110, j = 1, nctotl
            if (istate(j) .gt. 3  .or.  istate(j) .lt. 0) istate(j) = 0
  110    continue
      end if

*     Initialize nfixed, nfree and kactiv.
*     Ensure that the number of bounds and general constraints in the
*     working set does not exceed n.

      do 200, j = 1, nctotl
         if (nfixed + nactiv .eq. n) istate(j) = 0
         if (istate(j) .gt. 0) then
            if (j .le. n) then
               nfixed = nfixed + 1
               if (istate(j) .eq. 1) wx(j) = bl(j)
               if (istate(j) .ge. 2) wx(j) = bu(j)
            else
               nactiv = nactiv + 1
               kactiv(nactiv) = j - n
            end if
         end if
  200 continue

*     ------------------------------------------------------------------
*     If a cold start is required,  attempt to add as many
*     constraints as possible to the working set.
*     ------------------------------------------------------------------
      if (cold) then

*        See if any bounds are violated or nearly satisfied.
*        If so,  add these bounds to the working set and set the
*        variables exactly on their bounds.

         j = n
*+       while (j .ge. 1  .and.  nfixed + nactiv .lt. n) do
  300    if    (j .ge. 1  .and.  nfixed + nactiv .lt. n) then
            if (istate(j) .eq. 0) then
               b1     = bl(j)
               b2     = bu(j)
               is     = 0
               if (b1 .gt. biglow) then
                  if (wx(j) - b1 .le. (one + abs( b1 ))*tolact) is = 1
               end if
               if (b2 .lt. bigupp) then
                  if (b2 - wx(j) .le. (one + abs( b2 ))*tolact) is = 2
               end if
               if (is .gt. 0) then
                  istate(j) = is
                  if (is .eq. 1) wx(j) = b1
                  if (is .eq. 2) wx(j) = b2
                  nfixed = nfixed + 1
               end if
            end if
            j = j - 1
            go to 300
*+       end while
         end if

*        ---------------------------------------------------------------
*        The following loop finds the linear constraint (if any) with
*        smallest residual less than or equal to tolact  and adds it
*        to the working set.  This is repeated until the working set
*        is complete or all the remaining residuals are too large.
*        ---------------------------------------------------------------
*        First, compute the residuals for all the constraints not in the
*        working set.

         if (nclin .gt. 0  .and.  nactiv+nfixed .lt. n) then
            do 410, i = 1, nclin
               if (istate(n+i) .le. 0)
     $         Ax(i) = ddot  (n, A(i,1), ldA, wx, 1 )
  410       continue

            is     = 1
            toobig = tolact + tolact

*+          while (is .gt. 0  .and.  nfixed + nactiv .lt. n) do
  500       if    (is .gt. 0  .and.  nfixed + nactiv .lt. n) then
               is     = 0
               resmin = tolact

               do 520, i = 1, nclin
                  j      = n + i
                  if (istate(j) .eq. 0) then
                     b1     = bl(j)
                     b2     = bu(j)
                     resl   = toobig
                     resu   = toobig
                     if (b1 .gt. biglow)
     $                  resl  = abs( Ax(i) - b1 ) / (one + abs( b1 ))
                     if (b2 .lt. bigupp)
     $                  resu  = abs( Ax(i) - b2 ) / (one + abs( b2 ))
                     residl   = min( resl, resu )
                     if(residl .lt. resmin) then
                        resmin = residl
                        imin   = i
                        is     = 1
                        if (resl .gt. resu) is = 2
                     end if
                  end if
  520          continue

               if (is .gt. 0) then
                  nactiv = nactiv + 1
                  kactiv(nactiv) = imin
                  j         = n + imin
                  istate(j) = is
               end if
               go to 500
*+          end while
            end if
         end if
      end if
            
      if (vertex  .and.  nactiv+nfixed .lt. n) then
*        ---------------------------------------------------------------
*        Find an initial vertex by temporarily fixing some variables.
*        ---------------------------------------------------------------
*        Compute lengths of columns of selected linear constraints
*        (just the ones corresponding to variables eligible to be
*        temporarily fixed).        

         do 630, j = 1, n
            if (istate(j) .eq. 0) then
               colsiz = zero
               do 620, k = 1, nclin
                  if (istate(n+k) .gt. 0)
     $            colsiz = colsiz + abs( A(k,j) )
  620          continue
               work(j) = colsiz
            end if
  630    continue
         
*        Find the  nartif  smallest such columns.
*        This is an expensive loop.  Later we can replace it by a
*        4-pass process (say), accepting the first col that is within
*        t  of  colmin, where  t = 0.0, 0.001, 0.01, 0.1 (say).
*        (This comment written in 1980).
         
*+       while (nfixed + nactiv .lt. n) do
  640    if    (nfixed + nactiv .lt. n) then
            colmin = flmax
            do 650, j = 1, n
               if (istate(j) .eq. 0) then
                  if (nclin .eq. 0) go to 660
                  colsiz = work(j)
                  if (colmin .gt. colsiz) then
                     colmin = colsiz
                     jmin   = j
                  end if
               end if
  650       continue
            j         = jmin
  660       istate(j) = 4
            nartif    = nartif + 1
            nfixed    = nfixed + 1
            go to 640
*+       end while
         end if
      end if
      
      nfree = n - nfixed

*      if (lsdbg) then
*         if (ilsdbg(1) .gt. 0)
*     $       write(iPrint, 1300) nfixed, nactiv, nartif
*         if (ilsdbg(2) .gt. 0)
*     $       write(iPrint, 1200) (wx(j), j = 1, n)
*      end if

      return

 1000 format(/ ' //lscrsh// cold nclin nctotl'
     $       / ' //lscrsh// ', l4, i6, i7 )
 1100 format(/ ' //lscrsh// Variables before crash... '/ (5g12.3))
 1200 format(/ ' //lscrsh// Variables after  crash... '/ (5g12.3))
 1300 format(/ ' //lscrsh// Working set selected ...             '
     $       / ' //lscrsh// nfixed nactiv nartif      '
     $       / ' //lscrsh// ', i6, 2i7 )

*     end of lscrsh
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lsdel ( unitQ,
     $                   n, nactiv, nfree, nres, ngq, nZ, nZr,
     $                   ldA, ldQ, ldR, ldT, nrank,
     $                   jdel, kdel, kactiv, kx,
     $                   A, res, R, T, gq, Q,
     $                   c, s )

      implicit           double precision(a-h,o-z)
      logical            unitQ
      integer            kactiv(n), kx(n)
      double precision   A(ldA,*), res(n,*), R(ldR,*), T(ldT,*),
     $                   gq(n,*), Q(ldQ,*)
      double precision   c(n), s(n)

*     ==================================================================
*     lsdel   updates the least-squares factor R and the factorization
*     A(free) (Z Y) = (0 T) when a regular, temporary or artificial
*     constraint is deleted from the working set.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 31-October-1984.
*     Level-2 matrix routines added 25-Apr-1988.
*     This version of lsdel dated  10-Sep-92.
*     ==================================================================
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2
      common    /sol5cm/ Asize, dTmax, dTmin

      logical            lsdbg
      parameter         (ldbg = 5)
      common    /lsdebg/ ilsdbg(ldbg), lsdbg
      
      intrinsic          max   , min
      external           idamax
      parameter        ( zero = 0.0d+0, one = 1.0d+0 )

      if (jdel .gt. 0) then
*        ---------------------------------------------------------------
*        Regular constraint or temporary bound deleted.
*        ---------------------------------------------------------------

         if (jdel .le. n) then

*           Case 1.  A simple bound has been deleted.
*           =======  Columns nfree+1 and ir of r must be swapped.

            ir     = nZ    + kdel
*            if (lsdbg  .and.  ilsdbg(1) .gt. 0)
*     $         write(iPrint, 1100) nactiv, nZ, nfree, ir, jdel, unitQ

            itdel  = 1
            nfree  = nfree + 1

            if (nfree .lt. ir) then
               kx(ir)    = kx(nfree)
               kx(nfree) = jdel
               if (nrank .gt. 0)
     $            call cmrswp( n, nres, nrank, ldR, nfree, ir,
     $                         R, res, c, s )
               call dswap ( ngq, gq(nfree,1), n, gq(ir,1), n )
            end if

            if (.not. unitQ) then

*              Copy the incoming column of  A(free)  into the end of T.

               do 130, ka = 1, nactiv
                  i = kactiv(ka)
                  T(ka,nfree) = A(i,jdel)
  130          continue

*              Expand Q by adding a unit row and column.

               if (nfree .gt. 1) then
                  call dload ( nfree-1, zero, Q(nfree,1), ldQ )
                  call dload ( nfree-1, zero, Q(1,nfree), 1  )
               end if
               Q(nfree,nfree) = one
            end if
         else        

*           Case 2.  A general constraint has been deleted.
*           =======

*            if (lsdbg  .and.  ilsdbg(1) .gt. 0)
*     $         write(iPrint, 1200) nactiv, nZ, nfree, kdel, jdel, unitQ

            itdel  = kdel
            nactiv = nactiv - 1

*           Delete row  kdel  of T and move up the ones below it.
*           T becomes reverse lower Hessenberg.

            do 220, i = kdel, nactiv
               kactiv(i) = kactiv(i+1)
               ld        = nfree - i
               call dcopy ( i+1, T(i+1,ld), ldT, T(i,ld), ldT )
  220       continue
         end if

         nZ    = nZ     + 1

         if (nactiv .eq. 0) then
            dTmax = one
            dTmin = one
         else
*           ------------------------------------------------------------
*           Restore the nactiv by (nactiv+1) reverse-Hessenberg matrix 
*           T  to reverse-triangular form.  The last  nactiv  super-
*           diagonal elements are removed using a backward sweep of
*           plane rotations.  The rotation for the singleton in the 
*           first column is generated separately.
*           ------------------------------------------------------------
            nsup   = nactiv - itdel + 1
                                              
            if (nsup .gt. 0) then
               npiv   = nfree  - itdel + 1

               if (nsup .gt. 1) then
                  call dcopy ( nsup-1, T(nactiv-1,nZ+1), ldT-1, 
     $                         s(nZ+1), 1)
                  call f06qzf( 'Remove', nactiv, 1, nsup, 
     $                         c(nZ+1), s(nZ+1), T(1,nZ+1), ldT )
               end if

               call f06baf( T(nactiv,nZ+1), T(nactiv,nZ), cs, sn )
               T(nactiv,nZ) = zero
               s(nZ)   = - sn
               c(nZ)   =   cs
         
               call f06qxf( 'Right', 'Variable', 'Backwards', 
     $                      nfree, nfree, nZ, npiv, c, s, Q, ldQ )
               call f06qxf( 'Left ', 'Variable', 'Backwards', 
     $                      npiv , ngq  , nZ, npiv, c, s, gq, n    )
            
               nT = min( nrank, npiv )
               
               if (nT .lt. npiv  .and.  nT .gt. 0) then
               
*                 R is upper trapezoidal, pretend R is (nt x n) and 
*                 apply the rotations in columns  max(nt,nZ)  thru npiv.
               
                  call f06qxf( 'Right', 'Variable', 'Backwards', nT, n,
     $                         max(nt,nZ), npiv, c, s, R, ldR )
               end if
               
*              Apply the column transformations to the triangular part 
*              of  R.  The arrays  c  and  s  containing the column
*              rotations are overwritten by the row rotations that
*              restore  R  to upper-triangular form.
               
               if (nZ .lt. nT) then
                  call f06qtf( 'Right', nT, nZ, nT, c, s, R, ldR )
               end if
               
*              Apply the row rotations to the remaining rows of R.
               
               if (n .gt. nT)
     $            call f06qxf( 'Left', 'Variable', 'Backwards', 
     $                         nT, n-nT, nZ, nT, c, s, R(1,nT+1), ldR )
               
               if (nres .gt. 0)
     $            call f06qxf( 'Left', 'Variable', 'Backwards', 
     $                         nT, nres, nZ, nT, c, s, res, n )
            end if
            
            call dcond ( nactiv, T(nactiv,nZ+1), ldT-1, dTmax, dTmin )
         end if
      end if

      nZr1 = nZr + 1

      if (nZ .gt. nZr) then
         if (jdel .gt. 0) then
            jart =   nZr1 - 1 + idamax( nZ-nZr1+1, gq(nZr1,1), 1 )
         else
            jart = - jdel
         end if

*         if (lsdbg  .and.  ilsdbg(1) .gt. 0)
*     $      write( iPrint, 1000 ) nZ, nZr1, jart

         if (jart .gt. nZr1) then

*           Swap columns NZR1 and JART of R.

            if (unitQ) then
               k        = kx(nZr1)
               kx(nZr1) = kx(jart)
               kx(jart) = k
            else
               call dswap ( nfree, Q(1,nZr1), 1, Q(1,jart), 1 )
            end if

            call dswap ( ngq, gq(nZr1,1), n, gq(jart,1), n )
            if (nrank .gt. 0)
     $         call cmrswp( n, nres, nrank, ldR, nZr1, jart,
     $                      R, res, c, s )
         end if
      end if

      nZr = nZr1

      return

 1000 format(/ ' //lsdel //  Artificial constraint deleted.      '
     $       / ' //lsdel //      nZ   nZr   jart                 '
     $       / ' //lsdel //  ', 3i6 )
 1100 format(/ ' //lsdel //  Simple bound deleted.               '
     $       / ' //lsdel //  nactiv    nZ nfree    ir  jdel unitQ'
     $       / ' //lsdel //  ', 5i6, l6 )
 1200 format(/ ' //lsdel //  General constraint deleted.         '
     $       / ' //lsdel //  nactiv    nZ nfree  kdel  jdel unitQ'
     $       / ' //lsdel //  ', 5i6, l6 )

*     end of lsdel
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lsdflt( m, n, nclin, title )

      implicit           double precision(a-h,o-z)

      character*(*)      title

*     ==================================================================
*     lsdflt  loads the default values of the parameters not set by 
*     the user.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original Fortran 77 version written 17-September-1985.
*     This version of lsdflt dated  21-Mar-93.
*     ==================================================================
      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9

      logical            cmdbg, lsdbg
      parameter         (ldbg = 5)
      common    /lsdebg/ ilsdbg(ldbg), lsdbg
      common    /cmdebg/ icmdbg(ldbg), cmdbg

      logical            newopt
      common    /sol3ls/ newopt
      save      /sol3ls/

*     +Include lsparm+++++++++++++++++++++++++++++++++++++++++++++++++++
      parameter         (mxparm = 30)
      integer            iprmls(mxparm), ipsvls
      double precision   rprmls(mxparm), rpsvls

      common    /lspar1/ ipsvls(mxparm),
     $                   idbgls, iPrnt , iSumry, itmax1, itmax2, lcrash,
     $	                 ldbgls, lprob , msgls , nn    , nnclin, nprob , 
     $                   ipadls(18)

      common    /lspar2/ rpsvls(mxparm),
     $                   bigbnd, bigdx , bndlow, bndupp, tolact, tolfea,
     $                   tolrnk, rpadls(23)

      equivalence       (iprmls(1), idbgls), (rprmls(1), bigbnd)

      save      /lspar1/, /lspar2/
*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      equivalence   (msgls , msglvl), (idbgls, idbg), (ldbgls, msgdbg)

      character*4        icrsh(0:2)
      character*3        lstype(1:10)
      character*16       key
      parameter         (zero   =  0.0d+0, ten    = 10.0d+0)
      parameter         (hundrd =100.0d+0)
      parameter         (rdummy = -11111., idummy = -11111)
      parameter         (gigant = 1.0d+20*.99999          )
      parameter         (wrktol = 1.0d-2                  )
      data               icrsh(0), icrsh(1), icrsh(2)
     $                 /'Cold'   ,'Warm'   ,'Hot '   /
      data               lstype(1), lstype(2)
     $                 /' FP'     ,' LP'     /
      data               lstype(3), lstype(4), lstype(5), lstype(6)
     $                 /'QP1'     ,'QP2'     ,'QP3'     ,'QP4'     /
      data               lstype(7), lstype(8), lstype(9), lstype(10)
     $                 /'LS1'     ,'LS2'     ,'LS3'     ,'LS4'     /

      epsmch = wmach( 3)

*     Make a dummy call to lskey to ensure that the defaults are set.

      call lskey ( nout, '*', key )
      newopt = .true.

*     Save the optional parameters set by the user.  The values in
*     rprmls and iprmls may be changed to their default values.

      call icopy ( mxparm, iprmls, 1, ipsvls, 1 )
      call dcopy ( mxparm, rprmls, 1, rpsvls, 1 )

      if (       iPrnt  .lt. 0      )  iPrnt   = nout
      if (       iSumry .lt. 0      )  iSumry  = 6
      if (       iSumry .eq. iPrnt  )  iSumry  = 0
                                       iPrint  = iPrnt
                                       iSumm   = iSumry
      if (       lprob  .lt. 0      )  lprob   = 7
      if (       lcrash .lt. 0
     $    .or.   lcrash .gt. 2      )  lcrash  = 0
      if (       itmax1 .lt. 0      )  itmax1  = max(50, 5*(n+nclin))
      if (       itmax2 .lt. 0      )  itmax2  = max(50, 5*(n+nclin))
      if (       msglvl .eq. idummy )  msglvl  = 10
      if (       idbg   .lt. 0
     $    .or.   idbg   .gt. itmax1 + itmax2
     $                              )  idbg    = 0
      if (       msgdbg .lt. 0      )  msgdbg  = 0
      if (       msgdbg .eq. 0      )  idbg    = itmax1 + itmax2 + 1
      if (       tolact .lt. zero   )  tolact  = wrktol
      if (       tolfea .eq. rdummy
     $    .or.  (tolfea .ge. zero
     $    .and.  tolfea .lt. epsmch))  tolfea  = epspt5
      if (       tolrnk .le. zero
     $    .and. (lprob  .eq. 5  .or.
     $           lprob  .eq. 7  .or.
     $           lprob  .eq. 9)     )  tolrnk  = hundrd*epsmch
      if (       tolrnk .le. zero   )  tolrnk  =    ten*epspt5
      if (       bigbnd .le. zero   )  bigbnd  = gigant
      if (       bigdx  .le. zero   )  bigdx   = max(gigant, bigbnd)

      lsdbg = idbg .eq. 0  .and.  iPrint .gt. 0
      cmdbg = lsdbg
      k     = 1
      msg   = msgdbg
      do 200, i = 1, ldbg
         ilsdbg(i) = mod( msg/k, 10 )
         icmdbg(i) = ilsdbg(i)
         k = k*10
  200 continue

      if (msglvl .gt. 0) then

*        Print the title.

         lenT = len( title )
         if (lenT .gt. 0) then
            nspace = (81 - lenT)/2 + 1
*            if (iPrint .gt. 0) then
*               write(iPrint, '(///// (80a1) )')
*     $            (' ', j=1, nspace), (title(j:j), j=1,lenT)
*               write(iPrint, '(80a1 //)')
*     $            (' ', j=1, nspace), ('='       , j=1,lenT)
*            end if

*            if (iSumm .gt. 0) then
*               write(iSumm , '(///// (80a1) )')
*     $            (' ', j=1, nspace), (title(j:j), j=1,lenT)
*               write(iSumm , '(80a1 //)')
*     $            (' ', j=1, nspace), ('='       , j=1,lenT)
*            end if
         end if

*         if (iPrint .gt. 0) then
*            write(iPrint, 2000)
*            write(iPrint, 2100) lstype(lprob),
*     $                          nclin , tolfea, icrsh(lcrash),
*     $                          n     , bigbnd, tolact,
*     $                          m     , bigdx , tolrnk
*            write(iPrint, 2200) msglvl, iPrnt , itmax1, 
*     $                          epsmch, iSumry, itmax2
*         end if
      end if

      return

 2000 format(
     $//' Parameters'
     $/ ' ----------' )
 2100 format(
     $/ ' Problem type...........', 7x, a3
     $/ ' Linear constraints.....',     i10,   6x,
     $  ' Feasibility tolerance..', 1p, e10.2, 6x,
     $  1x, a4, ' start.............'
     $/ ' Variables..............',     i10,   6x,
     $  ' Infinite bound size....', 1p, e10.2, 6x,
     $  ' Crash tolerance........',     e10.2
     $/ ' Objective matrix rows..',     i10,   6x,
     $  ' Infinite step size.....', 1p, e10.2, 6x,
     $  ' Rank tolerance.........',     e10.2 )
 2200 format(
     $/ ' Print level............',     i10,   6x,
     $  ' Print file.............',     i10,   6x,
     $  ' Feasibility phase itns.',     i10
     $/ ' eps (machine precision)', 1p, e10.2, 6x,
     $  ' Summary file...........',     i10,   6x,
     $  ' Optimality  phase itns.',     i10 )

*     end of lsdflt
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lsfeas( n, nclin, istate,
     $                   bigbnd, cvnorm, errmax, jmax, nviol,
     $                   Ax, bl, bu, featol, x, work )

      implicit           double precision(a-h,o-z)
      integer            istate(n+nclin)
      double precision   Ax(*), bl(n+nclin), bu(n+nclin)
      double precision   featol(n+nclin), x(n)
      double precision   work(n+nclin)

*     ==================================================================
*     lsfeas  computes the following...
*     (1)  The number of constraints that are violated by more
*          than  featol  and the 2-norm of the constraint violations.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version      April    1984.
*     This version of  lsfeas  dated  17-October-1985.
*     ==================================================================
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2

      logical            lsdbg
      parameter         (ldbg = 5)
      common    /lsdebg/ ilsdbg(ldbg), lsdbg

      external           idamax, dnrm2
      intrinsic          abs
      parameter        ( zero = 0.0d+0 )

      biglow = - bigbnd
      bigupp =   bigbnd

*     ==================================================================
*     Compute nviol,  the number of constraints violated by more than
*     featol,  and cvnorm,  the 2-norm of the constraint violations and
*     residuals of the constraints in the working set.
*     ==================================================================
      nviol  = 0

      do 200, j = 1, n+nclin
         feasj  = featol(j)
         is     = istate(j)
         res    = zero

         if (is .ge. 0  .and.  is .lt. 4) then
            if (j .le. n) then
               con =  x(j)
            else
               i   = j - n
               con = Ax(i)
            end if

            tolj   = feasj

*           Check for constraint violations.

            if (bl(j) .gt. biglow) then
               res    = bl(j) - con
               if (res .gt.   feasj ) nviol = nviol + 1
               if (res .gt.    tolj ) go to 190
            end if

            if (bu(j) .lt. bigupp) then
               res    = bu(j) - con
               if (res .lt. (-feasj)) nviol = nviol + 1
               if (res .lt.  (-tolj)) go to 190
            end if

*           This constraint is satisfied,  but count the residual as a
*           violation if the constraint is in the working set.

            if (is .le. 0) res = zero
            if (is .eq. 1) res = bl(j) - con
            if (is .ge. 2) res = bu(j) - con
            if (abs( res ) .gt. feasj) nviol = nviol + 1
         end if
  190    work(j) = res
  200 continue

      jmax   = idamax( n+nclin, work, 1 )
      errmax = abs ( work(jmax) )

*      if (lsdbg  .and.  ilsdbg(1) .gt. 0)
*     $   write(iPrint, 1000) errmax, jmax

      cvnorm  = dnrm2 ( n+nclin, work, 1 )

      return

 1000 format(/ ' //lsfeas//  The maximum violation is ', 1pe14.2,
     $                     ' in constraint', i5 )

*     end of lsfeas
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lsfile( ioptns, inform )
      integer            ioptns, inform

*     ==================================================================
*     lsfile  reads the options file from unit  ioptns  and loads the
*     options into the relevant elements of  iprmls  and  rprmls.
*
*     If  ioptns .lt. 0  or  ioptns .gt. 99  then no file is read,
*     otherwise the file associated with unit  ioptns  is read.
*
*     Output:
*
*         inform = 0  if a complete  options  file was found
*                     (starting with  begin  and ending with  end);
*                  1  if  ioptns .lt. 0  or  ioptns .gt. 99;
*                  2  if  begin  was found, but end-of-file
*                     occurred before  end  was found;
*                  3  if end-of-file occurred before  begin  or
*                     endrun  were found;
*                  4  if  endrun  was found before  begin.
*     ==================================================================
      logical             newopt
      common     /sol3ls/ newopt
      save       /sol3ls/

      double precision    wmach(15)
      common     /solmch/ wmach
      save       /solmch/

      external            mchpar, lskey
      logical             first
      save                first , nout
      data                first /.true./

*     If first time in, set nout.
*     newopt is true first time into lsfile or lsoptn
*     and just after a call to lssol.

      if (first) then
         first  = .false.
         newopt = .true.
         call mchpar()
         nout = wmach(11)
      end if

      call opfile( ioptns, nout, inform, lskey )


*     end of lsfile
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lsgetp( linObj, singlr, unitgZ, unitQ,
     $                   n, nclin, nfree,
     $                   ldA, ldQ, ldR, nrank, numinf, nZr,
     $                   kx, ctp, pnorm,
     $                   A, Ap, res, hZ, p,
     $                   gq, cq, R, Q, work )

      implicit           double precision(a-h,o-z)
      logical            linObj, singlr, unitgZ, unitQ
      integer            kx(n)
      double precision   A(ldA,*), Ap(*), res(*), hZ(*), p(n),
     $                   gq(n), cq(*), R(ldR,*), Q(ldQ,*)
      double precision   work(n)

*     ==================================================================
*     lsgetp  computes the following quantities for  lscore.
*     (1) The vector  (hZ) = (Rz)(pz).
*         If X is not yet feasible,  the product is computed directly.
*         If  Rz is singular,  hZ  is zero.  Otherwise  hZ  satisfies
*         the equations
*                        Rz'hZ = -gz,
*         where  g  is the total gradient.  If there is no linear term
*         in the objective,  hZ  is set to  dz  directly.
*     (2) The search direction P (and its 2-norm).  The vector P is
*         defined as  Z*(pz), where  (pz)  depends upon whether or
*         not x is feasible and the nonsingularity of  (Rz).
*         If  numinf .gt. 0,  (pz)  is the steepest-descent direction.
*         Otherwise,  x  is the solution of the  nZr*nZr  triangular
*         system   (Rz)*(pz) = (hZ).
*     (3) The vector Ap,  where A is the matrix of linear constraints.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 31-October-1984.
*     Level 2 Blas added 11-June-1986.
*     This version of lsgetp dated  23-Oct-92.
*     ==================================================================
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2

      logical            lsdbg
      parameter         (ldbg = 5)
      common    /lsdebg/ ilsdbg(ldbg), lsdbg

      external           ddot  , dnrm2
      intrinsic          min
      parameter        ( zero = 0.0d+0, one  = 1.0d+0 )

      if (singlr) then
*        ---------------------------------------------------------------
*        The triangular factor for the current objective function is
*        singular,  i.e., the objective is linear along the last column
*        of Zr.  This can only occur when unitgZ is true.
*        ---------------------------------------------------------------
         if (nZr .gt. 1) then
            call dcopy ( nZr-1, R(1,nZr), 1, p, 1 )
            call dtrsv ( 'U', 'N', 'N', nZr-1, R, ldR, p, 1 )
         end if
         p(nZr) = - one

         gtp = ddot  ( nZr, gq, 1, p, 1 )
         if (gtp .gt. zero) call dscal ( nZr, (-one), p, 1 )

         if (nZr .le. nrank) then
            if (numinf .eq. 0) then
               if (unitgZ) then
                  hZ(nZr) = R(nZr,nZr)*p(nZr)
               else
                  call dload ( nZr, (zero), hZ, 1 )
               end if
            else
               hZ(1)   = R(1,1)*p(1)
            end if
         end if
      else
*        ---------------------------------------------------------------
*        The objective is quadratic in the space spanned by Zr.
*        ---------------------------------------------------------------
         if (linObj) then
            if (unitgZ) then
               if (nZr .gt. 1)
     $            call dload ( nZr-1, (zero), hZ, 1 )
               hZ(nZr) = - gq(nZr)/R(nZr,nZr)
            else
               call dcopy ( nZr, gq  , 1, hZ, 1 )
               call dscal ( nZr, (-one), hZ, 1 )
               call dtrsv ( 'U', 'T', 'N', nZr, R, ldR, hZ, 1 )
            end if
         else
            call dcopy ( nZr, res, 1, hZ, 1 )
         end if

*        Solve  Rz*pz = hZ.

         call dcopy ( nZr, hZ, 1, p, 1 )
         call dtrsv ( 'U', 'N', 'N', nZr, R, ldR, p, 1 )
      end if

*     Compute  p = Zr*pz  and its norm.

      if (linObj)
     $   ctp = ddot  ( nZr, cq, 1, p, 1 )
      pnorm  = dnrm2 ( nZr, p, 1 )

      call cmqmul( 1, n, nZr, nfree, ldQ, unitQ, kx, p, Q, work )

*      if (lsdbg  .and.  ilsdbg(2) .gt. 0)
*     $   write(iPrint, 1000) (p(j), j = 1, n)

*     Compute  Ap.

      if (nclin .gt. 0) then
         call dgemv ( 'No transpose', nclin, n, one, A, ldA,
     $                p, 1, zero, Ap, 1 )

*         if (lsdbg  .and.  ilsdbg(2) .gt. 0)
*     $      write(iPrint, 1100) (ap(i), i = 1, nclin)
      end if

      return

 1000 format(/ ' //lsgetp//   p ... ' / (1p, 5e15.5))
 1100 format(/ ' //lsgetp//  Ap ... ' / (1p, 5e15.5))

*     end of lsgetp
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lsgset( prbtyp, linObj, singlr, unitgZ, unitQ,
     $                   n, nclin, nfree,
     $                   ldA, ldQ, ldR, nrank, nZ, nZr,
     $                   istate, kx,
     $                   bigbnd, tolrnk, numinf, suminf,
     $                   bl, bu, A, res, featol,
     $                   gq, cq, R, x, wtinf, Q, wrk )

      implicit           double precision(a-h,o-z)
      character*2        prbtyp
      logical            linObj, singlr, unitgZ, unitQ
      integer            istate(*), kx(n)
      double precision   bl(*), bu(*), A(ldA,*),
     $                   res(*), featol(*)
      double precision   gq(n), cq(*), R(ldR,*), x(n), wtinf(*),
     $                   Q(ldQ,*)
      double precision   wrk(n)

*     ==================================================================
*     lsgset  finds the number and weighted sum of infeasibilities for
*     the bounds and linear constraints.   An appropriate transformed
*     gradient vector is returned in  GQ.
*
*     Positive values of  ISTATE(j)  will not be altered.  These mean
*     the following...
*
*               1             2           3
*           a'x = bl      a'x = bu     bl = bu
*
*     Other values of  ISTATE(j)  will be reset as follows...
*           a'x lt bl     a'x gt bu     a'x free
*              - 2           - 1           0
*
*     If  x  is feasible,  LSGSET computes the vector Q(free)'g(free),
*     where  g  is the gradient of the the sum of squares plus the
*     linear term.  The matrix Q is of the form
*                    ( Q(free)  0       ),
*                    (   0      I(fixed))
*     where  Q(free)  is the orthogonal factor of  A(free)  and  A  is
*     the matrix of constraints in the working set.  The transformed
*     gradients are stored in GQ.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 31-October-1984.
*     Level 2 Blas added 11-June-1986.
*     This version of lsgset dated 14-Sep-92.
*     ==================================================================
      external           ddot  , idrank
      intrinsic          abs   , max   , min
      parameter        ( zero = 0.0d+0, one = 1.0d+0 )
                     
      bigupp =   bigbnd
      biglow = - bigbnd

      numinf =   0
      suminf =   zero
      call dload ( n, zero, gq, 1 )

      do 200, j = 1, n+nclin
         if (istate(j) .le. 0) then
            feasj  = featol(j)
            if (j .le. n) then
               ctx = x(j)
            else
               k   = j - n
               ctx = ddot  ( n, A(k,1), ldA, x, 1 )
            end if
            istate(j) = 0

*           See if the lower bound is violated.

            if (bl(j) .gt. biglow) then
               s = bl(j) - ctx
               if (s     .gt. feasj ) then
                  istate(j) = - 2
                  weight    = - wtinf(j)
                  go to 160
               end if
            end if

*           See if the upper bound is violated.

            if (bu(j) .ge. bigupp) go to 200
            s = ctx - bu(j)
            if (s     .le. feasj ) go to 200
            istate(j) = - 1
            weight    =   wtinf(j)

*           Add the infeasibility.

  160       numinf = numinf + 1
            suminf = suminf + abs( weight ) * s
            if (j .le. n) then
               gq(j) = weight
            else
               call daxpy ( n, weight, A(k,1), ldA, gq, 1 )
            end if
         end if
  200 continue

*     ------------------------------------------------------------------
*     Install  gq,  the transformed gradient.
*     ------------------------------------------------------------------
      singlr = .false.
      unitgZ = .true.

      if (numinf .gt. 0) then
         call cmqmul( 6, n, nZ, nfree, ldQ, unitQ, kx, gq, Q, wrk )
      else if (numinf .eq. 0  .and.  prbtyp .eq. 'FP') then
         call dload ( n, zero, gq, 1 )
      else

*        Ready for the Optimality Phase.
*        Set nZr so that Rz is nonsingular.

         if (nrank .eq. 0) then
            if (linObj) then
               call dcopy ( n, cq, 1, gq, 1 )
            else
               call dload ( n, zero, gq, 1 )
            end if
            nZr    = 0
         else

*           Compute  gq = - R' * (transformed residual)

            call dcopy ( nrank, res, 1, gq, 1 )
            call dscal ( nrank, (-one), gq, 1 )
            call dtrmv ( 'U', 'T', 'N', nrank, R, ldR, gq, 1 )
            if (nrank .lt. n)
     $         call dgemv( 'T', nrank, n-nrank, -one,R(1,nrank+1),ldR,
     $                      res, 1, zero, gq(nrank+1), 1 )

            if (linObj) call daxpy ( n, one, cq, 1, gq, 1 )
            unitgZ = .false.

            rownrm = dnrm2 ( n, R(1,1), ldR )
            if (         rownrm  .le.        tolrnk
     $          .or. abs(R(1,1)) .le. rownrm*tolrnk) then
               nZr = 0
            else
               nZr = idrank( min(nrank, nZ), R, ldR+1, tolrnk )
            end if
         end if
      end if

*     end of lsgset
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lskey ( nout, buffer, key )

      implicit           double precision(a-h,o-z)
      character*(*)      buffer

*     ==================================================================
*     lskey   decodes the option contained in  buffer  in order to set
*     a parameter value in the relevant element of  iprmls  or  rprmls.
*
*
*     Input:
*
*     nout   a unit number for printing error messages.
*            if nout = 0 no error messages are printed.
*
*     Output:
*
*     key    The first keyword contained in BUFFER.
*
*
*     lskey  calls opnumb and the subprograms
*                 lookup, scannr, tokens, upcase
*     (now called oplook, opscan, optokn, opuppr)
*     supplied by Informatics General, Inc., Palo Alto, California.
*
*     Systems Optimization Laboratory, Stanford University.
*     This version of  lskey dated 19-Oct-92.
*     ==================================================================
*     +Include lsparm+++++++++++++++++++++++++++++++++++++++++++++++++++
      parameter         (mxparm = 30)
      integer            iprmls(mxparm), ipsvls
      double precision   rprmls(mxparm), rpsvls

      common    /lspar1/ ipsvls(mxparm),
     $                   idbgls, iPrnt , iSumry, itmax1, itmax2, lcrash,
     $	                 ldbgls, lprob , msgls , nn    , nnclin, nprob , 
     $                   ipadls(18)

      common    /lspar2/ rpsvls(mxparm),
     $                   bigbnd, bigdx , bndlow, bndupp, tolact, tolfea,
     $                   tolrnk, rpadls(23)

      equivalence       (iprmls(1), idbgls), (rprmls(1), bigbnd)

      save      /lspar1/, /lspar2/
*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      external           opnumb
      logical            first , more  , number, opnumb, sorted
      save               first

      parameter         (     maxkey = 28,  maxtie = 12,   maxtok = 10,
     $                        maxtyp = 16)
      character*16       keys(maxkey), ties(maxtie), token(maxtok),
     $                   type(maxtyp)
      character*16       key, key2, key3, value

      parameter         (idummy = -11111,  rdummy = -11111.0,
     $                   sorted = .true.,  zero   =  0.0     )

      data                first
     $                  /.true./
      data   keys
     $ / 'BEGIN           ',
     $   'COLD            ', 'CONSTRAINTS     ', 'CRASH           ',
     $   'DEBUG           ', 'DEFAULTS        ', 'END             ',
     $   'FEASIBILITY     ', 'HOT             ', 'INFINITE        ',
     $   'IPRMLS          ', 'ITERATIONS      ', 'ITERS:ITERATIONS',
     $   'ITNS :ITERATIONS', 'LINEAR          ', 'LIST            ',
     $   'LOWER           ', 'NOLIST          ', 'OPTIMALITY      ',
     $   'PRINT           ', 'PROBLEM         ', 'RANK            ',
     $   'RPRMLS          ', 'START           ', 'SUMMARY         ',
     $   'UPPER           ', 'VARIABLES       ', 'WARM            '/

      data   ties
     $ / 'BOUND           ', 'CONSTRAINTS     ', 'FILE            ',
     $   'LEVEL           ',
     $   'NO              ', 'NO.      :NUMBER', 'NUMBER          ',
     $   'PHASE           ', 'STEP            ',
     $   'TOLERANCE       ', 'TYPE            ', 'YES             '/

      data   type
     $ / 'FP              ',
     $   'LEAST       :LS1', 'LINEAR       :LP', 'LP              ',
     $   'LS          :LS1', 'LS1             ', 'LS2             ',
     $   'LS3             ', 'LS4             ', 'LSQ         :LS1',
     $   'QP          :QP2', 'QP1             ', 'QP2             ',
     $   'QP3             ', 'QP4             ', 'QUADRATIC   :QP2'/
*-----------------------------------------------------------------------

      if (first) then
         first  = .false.
         do 10, i = 1, mxparm
            iprmls(i) = idummy
            rprmls(i) = rdummy
   10    continue
      end if

*     Eliminate comments and empty lines.
*     A '*' appearing anywhere in buffer terminates the string.

      i      = index( buffer, '*' )
      if (i .eq. 0) then
         lenbuf = len( buffer )
      else
         lenbuf = i - 1
      end if
      if (lenbuf .le. 0) then
         key = '*'
         go to 900
      end if

*     ------------------------------------------------------------------
*     Extract up to maxtok tokens from the record.
*     ntoken returns how many were actually found.
*     key, key2, key3 are the first tokens if any, otherwise blank.
*     ------------------------------------------------------------------
      ntoken = maxtok
      call optokn( buffer(1:lenbuf), ntoken, token )
      key    = token(1)
      key2   = token(2)
      key3   = token(3)

*     Certain keywords require no action.

      if (key .eq. ' '     .or.  key .eq. 'BEGIN' ) go to 900
      if (key .eq. 'LIST'  .or.  key .eq. 'NOLIST') go to 900
      if (key .eq. 'END'                          ) go to 900

*     Most keywords will have an associated integer or real value,
*     so look for it no matter what the keyword.

      i      = 1
      number = .false.

   50 if (i .lt. ntoken  .and.  .not. number) then
         i      = i + 1
         value  = token(i)
         number = opnumb( value )
         go to 50
      end if

      if (number) then
*         read (value, '(bn, e16.0)') rvalue
      else
         rvalue = zero
      end if

*     Convert the keywords to their most fundamental form
*     (upper case, no abbreviations).
*     sorted says whether the dictionaries are in alphabetic order.
*     LOCi   says where the keywords are in the dictionaries.
*     LOCi = 0 signals that the keyword wasn't there.

      call oplook( maxkey, keys, sorted, key , loc1 )
      call oplook( maxtie, ties, sorted, key2, loc2 )

*     ------------------------------------------------------------------
*     Decide what to do about each keyword.
*     The second keyword (if any) might be needed to break ties.
*     Some seemingly redundant testing of more is used
*     to avoid compiler limits on the number of consecutive else ifs.
*     ------------------------------------------------------------------
      more   = .true.
      if (more) then
         more   = .false.
         if (key .eq. 'COLD        ') then
            lcrash = 0
         else if (key .eq. 'CONSTRAINTS ') then
            nnclin = rvalue
         else if (key .eq. 'CRASH       ') then
            tolact = rvalue
         else if (key .eq. 'DEBUG       ') then
            ldbgls = rvalue
         else if (key .eq. 'DEFAULTS    ') then
            do 20, i = 1, mxparm
               iprmls(i) = idummy
               rprmls(i) = rdummy
   20       continue
         else if (key .eq. 'FEASIBILITY ') then
              if (key2.eq. 'PHASE       ') itmax1 = rvalue
              if (key2.eq. 'TOLERANCE   ') tolfea = rvalue
              if (loc2.eq.  0            ) then
*                 if (nout .gt. 0)          write(nout, 2320) key2
              end if
         else
            more   = .true.
         end if
      end if

      if (more) then
         more   = .false.
         if (key .eq. 'HOT         ') then
            lcrash = 2
         else if (key .eq. 'INFINITE    ') then
              if (key2.eq. 'BOUND       ') bigbnd = rvalue * 0.99999
              if (key2.eq. 'STEP        ') bigdx  = rvalue
              if (loc2.eq.  0            ) then
*                 if (nout .gt. 0)          write(nout, 2320) key2
              end if
         else if (key .eq. 'IPRMLS      ') then
*           Allow things like  iprmls 21 = 100  to set iprmls(21) = 100
            ivalue = rvalue
            if (ivalue .ge. 1  .and. ivalue .le. mxparm) then
*               read (key3, '(bn, i16)') iprmls(ivalue)
            else
*               if (nout .gt. 0) write(nout, 2400) ivalue
            end if
         else if (key .eq. 'ITERATIONS  ') then
            itmax2 = rvalue
         else if (key .eq. 'LINEAR      ') then
            nnclin = rvalue
         else if (key .eq. 'LOWER       ') then
            bndlow = rvalue
         else
            more   = .true.
         end if
      end if

      if (more) then
         more   = .false.
         if      (key .eq. 'OPTIMALITY  ') then
            itmax2 = rvalue
         else if (key .eq. 'PROBLEM     ') then
            if      (key2 .eq. 'NUMBER') then
               nprob  = rvalue
            else if (key2 .eq. 'TYPE  ') then

*              Recognize     Problem type = LP     etc.

               call oplook( maxtyp, type, sorted, key3, loc3 )
               if (key3 .eq. 'FP' ) lprob = 1
               if (key3 .eq. 'LP' ) lprob = 2
               if (key3 .eq. 'QP1') lprob = 3
               if (key3 .eq. 'QP2') lprob = 4
               if (key3 .eq. 'QP3') lprob = 5
               if (key3 .eq. 'QP4') lprob = 6
               if (key3 .eq. 'LS1') lprob = 7
               if (key3 .eq. 'LS2') lprob = 8
               if (key3 .eq. 'LS3') lprob = 9
               if (key3 .eq. 'LS4') lprob = 10
               if (loc3 .eq.  0   ) then
*                  if (nout .gt. 0)  write(nout, 2330) key3
               end if
            else
*               if (nout .gt. 0) write(nout, 2320) key2
            end if
         else
            more   = .true.
         end if
      end if

      if (more) then
         more   = .false.
         if      (key .eq. 'PRINT       ') then
              if (key2.eq. 'FILE        ') iPrnt  = rvalue
              if (key2.eq. 'LEVEL       ') msgls  = rvalue
              if (loc2.eq.  0            ) then
*                 if (nout .gt. 0)          write(nout, 2320) key2
              end if
         else if (key .eq. 'RANK        ') then
            tolrnk = rvalue
         else if (key .eq. 'RPRMLS      ') then
*           Allow things like  rprmls 21 = 2  to set rprmls(21) = 2.0
            ivalue = rvalue
            if (ivalue .ge. 1  .and. ivalue .le. mxparm) then
*               read (key3, '(bn, e16.0)') rprmls(ivalue)
            else
*               if (nout .gt. 0) write(nout, 2400) ivalue
            end if
         else if (key .eq. 'START       ') then
            idbgls = rvalue
         else if (key .eq. 'SUMMARY     ') then
            iSumry = rvalue
         else if (key .eq. 'UPPER       ') then
            bndupp = rvalue
         else if (key .eq. 'VARIABLES   ') then
            nn     = rvalue
         else if (key .eq. 'WARM        ') then
            lcrash = 1
         else
*            if (nout .gt. 0) write(nout, 2300) key
         end if
      end if

  900 return

 2300 format(' XXX  Keyword not recognized:         ', a)
 2320 format(' XXX  Second keyword not recognized:  ', a)
 2330 format(' XXX  Third  keyword not recognized:  ', a)
 2400 format(' XXX  The parm subscript is out of range:', i10)

*     end of lskey
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lsloc ( lprob, n, nclin, litotl, lwtotl )

      implicit           double precision(a-h,o-z)

*     ==================================================================
*     lsloc   allocates the addresses of the work arrays for  lscore.
*
*     Note that the arrays  ( gq, cq )  and  ( res, res0, hz )  lie in
*     contiguous areas of workspace.
*     res, res0 and hZ are not needed for LP.
*     CQ is defined when the objective has an explicit linear term.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written  29-October-1984.
*     This version of lsloc dated 16-February-1986.
*     ==================================================================
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2
      common    /sol3cm/ lennam, ldT   , ncolt, ldQ

      parameter        ( lenls = 20 )
      common    /sol1ls/ locls(lenls)

      logical            lsdbg
      parameter        ( ldbg = 5 )
      common    /lsdebg/ ilsdbg(ldbg), lsdbg

      miniw     = litotl + 1
      minw      = lwtotl + 1


*     Assign array lengths that depend upon the problem dimensions.

      if (nclin .eq. 0) then
         lenT  = 0
         lenQ = 0
      else
         lenT  = ldT *ncolt
         lenQ = ldQ*ldQ
      end if

      lencq  = 0
      if (lprob .eq. 2*(lprob/2)) lencq  = n
      lenres = 0
      if (lprob .gt. 2          ) lenres = n

      lkactv    = miniw
      miniw     = lkactv + n

      lanorm    = minw
      lAp       = lanorm + nclin
      lpx       = lAp    + nclin
      lgq       = lpx    + n
      lcq       = lgq    + n
      lres      = lcq    + lencq
      lres0     = lres   + lenres
      lhZ       = lres0  + lenres
      lrlam     = lhZ    + lenres
      lT        = lrlam  + n
      lQ        = lT     + lenT
      lwtinf    = lQ     + lenQ
      lwrk      = lwtinf + n  + nclin
      lfeatl    = lwrk   + n  + nclin
      minw      = lfeatl + n  + nclin

      locls( 1) = lkactv
      locls( 2) = lanorm
      locls( 3) = lAp
      locls( 4) = lpx
      locls( 5) = lres
      locls( 6) = lres0
      locls( 7) = lhZ
      locls( 8) = lgq
      locls( 9) = lcq
      locls(10) = lrlam
      locls(11) = lT
      locls(12) = lQ
      locls(13) = lwtinf
      locls(14) = lwrk
      locls(15) = lfeatl

      litotl    = miniw - 1
      lwtotl    = minw  - 1

*     end of lsloc
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lsmove( hitcon, hitlow, linObj, unitgZ,
     $                   nclin, nrank, nZr,
     $                   n, ldR, jadd, numinf,
     $                   alfa, ctp, ctx, xnorm,
     $                   Ap, Ax, bl, bu, gq, hZ, p, res,
     $                   R, x, work )

      implicit           double precision (a-h,o-z)
      logical            hitcon, hitlow, linObj, unitgZ
      double precision   Ap(*), Ax(*), bl(*), bu(*), gq(*), hZ(*),
     $                   p(n), res(*), R(ldR,*), x(n)
      double precision   work(*)

*     ==================================================================
*     lsmove  changes x to x + alfa*p and updates ctx, Ax, res and gq
*     accordingly.
*
*     If a bound was added to the working set,  move x exactly on to it,
*     except when a negative step was taken (cmalf may have had to move
*     to some other closer constraint.)
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 27-December-1985.
*     Level 2 BLAS added 11-June-1986.
*     This version of lsmove dated 14-Sep-92.
*     ==================================================================
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2

      logical            lsdbg
      parameter         (ldbg = 5)
      common    /lsdebg/ ilsdbg(ldbg), lsdbg

      external           dnrm2
      intrinsic          abs   , min
      parameter        ( zero  = 0.0d+0, one = 1.0d+0 )

      call daxpy ( n, alfa, p, 1, x, 1 )
      if (linObj) ctx = ctx + alfa*ctp

      if (hitcon  .and.  jadd .le. n) then
         bnd = bu(jadd)
         if (hitlow) bnd = bl(jadd)
         if (alfa .ge. zero) x(jadd) = bnd
      end if
      xnorm  = dnrm2 ( n, x, 1 )

      if (nclin .gt. 0)
     $   call daxpy ( nclin, alfa, Ap, 1, Ax, 1 )

      if (nZr .le. nrank) then
         if (unitgZ) then
            res(nZr) = res(nZr) - alfa*hZ(nZr)
         else
            call daxpy ( nZr, (-alfa), hZ, 1, res, 1  )
         end if

         if (numinf .eq. 0) then

*           Update the transformed gradient GQ so that
*           gq = gq + alfa*R'( HZ ).
*                            ( 0  )

            if (unitgZ) then
               call daxpy ( n-nZr+1, alfa*hZ(nZr), R(nZr,nZr), ldR,
     $                                             gq(nZr)   , 1      )
            else
               call dcopy ( nZr, hZ, 1, work, 1 )
               call dtrmv ( 'U', 'T', 'N', nZr, R, ldR, work, 1 )
               if (nZr .lt. n)
     $            call dgemv ( 'T', nZr, n-nZr, one, R(1,nZr+1), ldR,
     $                         hZ, 1, zero, work(nZr+1), 1 )
               call daxpy ( n, alfa, work, 1, gq, 1 )
            end if
         end if
      end if

*     end of lsmove
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lsmuls( prbtyp,
     $                   msglvl, n, nactiv, nfree,
     $                   ldA, ldT, numinf, nZ, nZr,
     $                   istate, kactiv, kx, dinky,
     $                   jsmlst, ksmlst, jinf, jtiny,
     $                   jbigst, kbigst, trulam,
     $                   A, anorms, gq, rlamda, T, wtinf )

      implicit           double precision(a-h,o-z)
      character*2        prbtyp
      integer            istate(*), kactiv(n), kx(n)
      double precision   A(ldA,*), anorms(*),
     $                   gq(n), rlamda(n), T(ldT,*), wtinf(*)

*     ==================================================================
*     lsmuls  first computes the Lagrange multiplier estimates for the
*     given working set.  It then determines the values and indices of
*     certain significant multipliers.  In this process, the multipliers
*     for inequalities at their upper bounds are adjusted so that a
*     negative multiplier for an inequality constraint indicates non-
*     optimality.  All adjusted multipliers are scaled by the 2-norm
*     of the associated constraint row.  In the following, the term
*     minimum refers to the ordering of numbers on the real line,  and
*     not to their magnitude.
*
*     jsmlst  is the index of the minimum of the set of adjusted
*             multipliers with values less than  - dinky.  A negative
*             jsmlst defines the index in Q'g of the artificial
*             constraint to be deleted.
*     ksmlst  marks the position of general constraint jsmlst in kactiv.
*
*     jbigst  is the index of the largest of the set of adjusted
*             multipliers with values greater than (1 + dinky).
*     kbigst  marks its position in kactiv.
*
*     On exit,  elements 1 thru nactiv of rlamda contain the unadjusted
*     multipliers for the general constraints.  Elements nactiv onwards
*     of rlamda contain the unadjusted multipliers for the bounds.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 31-October-1984.
*     This version of lsmuls dated  14-Sep-92.
*     ==================================================================
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2

      logical            lsdbg
      parameter         (ldbg = 5)
      common    /lsdebg/ ilsdbg(ldbg), lsdbg

      intrinsic          abs, min
      parameter        ( one    =1.0d+0 )

      nfixed =   n - nfree

      jsmlst =   0
      ksmlst =   0
      smllst = - dinky

      tinylm =   dinky
      jtiny  =   0

      jbigst =   0
      kbigst =   0
      biggst =   one + dinky

      if (nZr .lt. nZ) then
*        ---------------------------------------------------------------
*        Compute jsmlst for the artificial constraints.
*        ---------------------------------------------------------------
         do 100, j = nZr+1, nZ
            rlam = - abs( gq(j) )
            if (rlam .lt. smllst) then
               smllst =   rlam
               jsmlst = - j
            else if (rlam .lt. tinylm) then
               tinylm =   rlam
               jtiny  =   j
            end if
  100    continue

         if (msglvl .ge. 20) then
*            if (iPrint .gt. 0) write(iPrint, 1000) (gq(k), k=nZr+1,nZ)
         end if
      end if

*     ------------------------------------------------------------------
*     Compute jsmlst for regular constraints and temporary bounds.
*     ------------------------------------------------------------------
*     First, compute the Lagrange multipliers for the general
*     constraints in the working set, by solving  T'*lamda = Y'g.

      if (n .gt. nZ)
     $   call dcopy ( n-nZ, gq(nZ+1), 1, rlamda, 1 )
      if (nactiv .gt. 0)
     $   call cmtsol( 2, ldT, nactiv, T(1,nZ+1), rlamda )

*     -----------------------------------------------------------------
*     Now set elements nactiv, nactiv+1,... of  rlamda  equal to
*     the multipliers for the bound constraints.
*     -----------------------------------------------------------------
      do 190, l = 1, nfixed
         j     = kx(nfree+l)
         blam  = rlamda(nactiv+l)
         do 170, k = 1, nactiv
            i    = kactiv(k)
            blam = blam - A(i,j)*rlamda(k)
  170    continue
         rlamda(nactiv+l) = blam
  190 continue

*     -----------------------------------------------------------------
*     Find jsmlst and ksmlst.
*     -----------------------------------------------------------------
      do 330, k = 1, n - nZ
         if (k .gt. nactiv) then
            j = kx(nZ+k)
         else
            j = kactiv(k) + n
         end if

         is   = istate(j)

         i    = j - n
         if (j .le. n) anormj = one
         if (j .gt. n) anormj = anorms(i)

         rlam = rlamda(k)

*        Change the sign of the estimate if the constraint is in
*        the working set at its upper bound.

         if (is .eq. 2) rlam =      - rlam
         if (is .eq. 3) rlam =   abs( rlam )
         if (is .eq. 4) rlam = - abs( rlam )

         if (is .ne. 3) then
            scdlam = rlam * anormj
            if      (scdlam .lt. smllst) then
               smllst = scdlam
               jsmlst = j
               ksmlst = k
            else if (scdlam .lt. tinylm) then
               tinylm = scdlam
               jtiny  = j
            end if
         end if

         if (numinf .gt. 0  .and.  j .gt. jinf) then
            scdlam = rlam/wtinf(j)
            if (scdlam .gt. biggst) then
               biggst = scdlam
               trulam = rlamda(k)
               jbigst = j
               kbigst = k
            end if
         end if
  330 continue

*     -----------------------------------------------------------------
*     If required, print the multipliers.
*     -----------------------------------------------------------------
*      if (msglvl .ge. 20  .and.  iPrint .gt. 0) then
*         if (nfixed .gt. 0)
*     $      write(iPrint, 1100) prbtyp, (kx(nfree+k),
*     $                         rlamda(nactiv+k), k=1,nfixed)
*         if (nactiv .gt. 0)
*     $      write(iPrint, 1200) prbtyp, (kactiv(k),
*     $                         rlamda(k), k=1,nactiv)
*      end if

*      if (lsdbg  .and.  ilsdbg(1) .gt. 0) then
*         write(iPrint, 9000) jsmlst, smllst, ksmlst
*         write(iPrint, 9100) jbigst, biggst, kbigst
*         write(iPrint, 9200) jtiny , tinylm
*      end if

      return

 1000 format(/ ' Multipliers for the artificial constraints        '
     $       / 4(5x, 1pe11.2))
 1100 format(/ ' Multipliers for the ', a2, ' bound  constraints   '
     $       / 4(i5, 1pe11.2))
 1200 format(/ ' Multipliers for the ', a2, ' linear constraints   '
     $       / 4(i5, 1pe11.2))
 9000 format(/ ' //lsmuls//  jsmlst     smllst     ksmlst (scaled) '
     $       / ' //lsmuls//  ', i6, 1pe11.2, 5x, i6 )
 9100 format(  ' //lsmuls//  jbigst     biggst     kbigst (scaled) '
     $       / ' //lsmuls//  ', i6, 1pe11.2, 5x, i6 )
 9200 format(  ' //lsmuls//   jtiny     tinylm                     '
     $       / ' //lsmuls//  ', i6, 1pe11.2)

*     end of lsmuls
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lsoptn( string )
      character*(*)      string

*     ==================================================================
*     lsoptn  loads the option supplied in  string  into the relevant
*     element of  iprmls  or  rprmls.
*     ==================================================================

      logical             newopt
      common     /sol3ls/ newopt
      save       /sol3ls/

      double precision    wmach(15)
      common     /solmch/ wmach
      save       /solmch/

      external            mchpar
      character*16        key
      character*72        buffer
      logical             first , prnt
      save                first , nout  , prnt
      data                first /.true./

*     If first time in, set  nout.
*     newopt  is true first time into  lsfile  or  lsoptn
*     and just after a call to  lssol.
*     prnt    is set to true whenever  newopt  is true.

      if (first) then
         first  = .false.
         newopt = .true.
         call mchpar()
         nout   =  wmach(11)
      end if
      buffer = string

*     Call  lskey   to decode the option and set the parameter value.
*     If newopt is true, reset prnt and test specially for nolist.

      if (newopt) then
         newopt = .false.
         prnt   = .true.
         call lskey ( nout, buffer, key )

         if (key .eq. 'NOLIST') then
            prnt   = .false.
         else
*            write(nout, '(// a / a /)')
*     $         ' Calls to Option Routine',
*     $         ' -----------------------'
*            write(nout, '( 6x, a )') buffer
         end if
      else
         if (prnt) then
            iPrint = nout
*            write(nout, '( 6x, a )') buffer
         else
            iPrint = 0
         end if

         call lskey ( iPrint, buffer, key )
         if (key .eq.   'list') prnt = .true.
         if (key .eq. 'nolist') prnt = .false.
      end if

*     end of lsoptn
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lsprt ( prbtyp, isdel, iter, jadd, jdel,
     $                   msglvl, nactiv, nfree, n, nclin,
     $                   nrank, ldR, ldT, nZ, nZr, istate,
     $                   alfa, condRz, condT, gfnorm, gZrnrm,
     $                   numinf, suminf, ctx, ssq,
     $                   Ax, R, T, x, work )

      implicit           double precision(a-h,o-z)
      character*2        prbtyp
      integer            istate(*)
      double precision   Ax(*), R(ldR,*), T(ldT,*), x(n)
      double precision   work(n)

*     ==================================================================
*     lsprt  prints various levels of output for  lscore.
*
*           Msg        Cumulative result
*           ---        -----------------
*
*        le   0        no output.
*
*        eq   1        nothing now (but full output later).
*
*        eq   5        one terse line of output.
*
*        ge  10        same as 5 (but full output later).
*
*        ge  20        constraint status,  x  and  Ax.
*
*        ge  30        diagonals of  T  and  R.
*
*
*     Debug printing is performed depending on the logical variable 
*     lsdbg, which is set true when  idbg  major iterations have been
*     performed.  At this point,  printing is done according to a string
*     of binary digits of the form  SVT  (stored in the integer array
*     ilsdbg).
*
*     S  set 'on' gives information from the maximum step routine cmalf.
*     V  set 'on' gives various vectors in  lscore  and its auxiliaries.
*     T  set 'on' gives a trace of which routine was called and an
*                 indication of the progress of the run.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 31-October-1984.
*     This version of lsprt dated 07-Jan-93.
*     ==================================================================
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2

      logical            lsdbg
      parameter         (ldbg = 5)
      common    /lsdebg/ ilsdbg(ldbg), lsdbg

      logical            first , linObj, newSet, prtHdr
      character*2        ladd, ldel
      character*2        lstate(0:5)
      parameter         (mLine1 = 40, mLine2 = 5)
      data               lstate(0), lstate(1), lstate(2)
     $                  /'  '     , 'L '     , 'U '     /
      data               lstate(3), lstate(4), lstate(5)
     $                  /'E '     , 'F '     , 'A '     /

      if (msglvl .ge. 15) then
*         if (iPrint .gt. 0) write(iPrint, 1000) prbtyp, iter
      end if

      if (msglvl .ge. 5) then

         first  = iter  .eq. 0
         linObj = nrank .eq. 0

         Itn    = mod( iter, 1000  )
         ndf    = mod( nZr , 10000 )

         nArt = nZ - nZr

         if      (jdel .gt. 0) then
            kdel =   isdel
         else if (jdel .lt. 0) then
            jdel =   nArt + 1
            kdel =   5
         else
            kdel =   0
         end if

         if (jadd .gt. 0) then
            kadd = istate(jadd)
         else
            kadd = 0
         end if

         ldel   = lstate(kdel)
         ladd   = lstate(kadd)

         if (numinf .gt. 0) then
            obj    = suminf
         else
            obj    = ssq + ctx
         end if

*        ---------------------------------------------------------------
*        If necessary, print a header. 
*        Print a single line of information.
*        ---------------------------------------------------------------
         if (iPrint .gt. 0) then
*           ------------------------------
*           Terse line for the Print file.
*           ------------------------------
            newSet = lines1 .ge. mLine1
            prtHdr = msglvl .ge. 15  .or.  first 
     $                               .or.  newSet

            if (prtHdr) then
               if (linObj) then 
*                  write(iPrint, 1200)
               else
*                  write(iPrint, 1300)
               end if
               lines1 = 0
            end if

            if (linObj) then
*               write(iPrint, 1700) Itn, jdel, ldel, jadd, ladd,
*     $                             alfa, numinf, obj, gZrnrm, ndf, nArt,
*     $                             n-nfree, nactiv, gfnorm, condT
            else
*               write(iPrint, 1700) Itn, jdel, ldel, jadd, ladd,
*     $                             alfa, numinf, obj, gZrnrm, ndf, nArt,
*     $                             n-nfree, nactiv, gfnorm, condT,
*     $                             CondRz
            end if
            lines1 = lines1 + 1
         end if

         if (iSumm .gt. 0) then
*           --------------------------------
*           Terse line for the Summary file.
*           --------------------------------
            newSet = lines2 .ge. mLine2
            prtHdr =                      first 
     $                              .or.  newSet

            if (prtHdr) then
*               write(iSumm , 1100)
               lines2 = 0
            end if
*            write(iSumm , 1700) Itn, jdel, ldel, jadd, ladd,
*     $                          alfa, numinf, obj, gZrnrm, ndf, nArt
            lines2 = lines2 + 1
         end if

         if (msglvl .ge. 20  .and.  iPrint .gt. 0) then
*            write(iPrint, 2000) prbtyp
*            write(iPrint, 2100) (x(j) , istate(j)  ,  j=1,n)
*            if (nclin .gt. 0)
*     $      write(iPrint, 2200) (Ax(k), istate(n+k), k=1,nclin )

            if (msglvl .ge. 30) then
*              ---------------------------------------------------------
*              Print the diagonals of  T  and  R.
*              ---------------------------------------------------------
               if (nactiv .gt. 0) then
                  call dcopy ( nactiv, T(nactiv,nZ+1), ldT-1, work,1 )
*                  write(iPrint, 3000) prbtyp, (work(j), j=1,nactiv)
               end if
*               if (nrank  .gt. 0)
*     $            write(iPrint, 3100) prbtyp, (R(j,j) , j=1,nrank )
            end if
*            write(iPrint, 5000)
         end if
      end if

      return

 1000 format(/// ' ', a2, ' iteration', i5
     $         / ' =================' )
 1100 format(// ' Itn Jdel  Jadd     Step Ninf  Sinf/Objective',
     $          ' Norm gZ   Zr  Art' )
 1200 format(// ' Itn Jdel  Jadd     Step Ninf  Sinf/Objective',
     $          ' Norm gZ   Zr  Art ',
     $          ' Bnd  Lin Norm gf  Cond T' )
 1300 format(// ' Itn Jdel  Jadd     Step Ninf  Sinf/Objective',
     $          ' Norm gZ   Zr  Art ',
     $          ' Bnd  Lin Norm gf  Cond T Cond Rz' )
 1700 format(    i4, i5, a1, i5, a1, 1p, e8.1, i5, e16.8, 
     $           e8.1, 2i5,
     $           2i5, e8.1, 2e8.0 )
 2000 format(/ ' Values and status of the ', a2, ' constraints'
     $       / ' ---------------------------------------' )
 2100 format(/ ' Variables...'                 / (1x, 5(1p, e15.6, i5)))
 2200 format(/ ' General linear constraints...'/ (1x, 5(1p, e15.6, i5)))
 3000 format(/ ' Diagonals of ' , a2,' working set factor T'
     $       /   (1p, 5e15.6))
 3100 format(/ ' Diagonals of ' , a2, ' triangle R         '
     $       /   (1p, 5e15.6))
 5000 format(/// ' ---------------------------------------------------',
     $           '--------------------------------------------' )

*     end of lsprt
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lssetx( linObj, rowerr, unitQ,
     $                   nclin, nactiv, nfree, nrank, nZ,
     $                   n, nctotl, ldQ, ldA, ldR, ldT,
     $                   istate, kactiv, kx,
     $                   jmax, errmax, ctx, xnorm,
     $                   A, Ax, bl, bu, cq, res, res0, featol,
     $                   R, T, x, Q, p, work )

      implicit           double precision (a-h,o-z)
      logical            linObj, rowerr, unitQ
      integer            istate(nctotl), kactiv(n), kx(n)
      double precision   A(ldA,*), Ax(*), bl(nctotl), bu(nctotl),
     $                   cq(*), res(*), res0(*), featol(nctotl), p(n),
     $                   R(ldR,*), T(ldT,*), Q(ldQ,*), x(n)
      double precision   work(nctotl)

*     ==================================================================
*     lssetx  computes the point on a working set that is closest to the
*     input vector  x  (in the least-squares sense).  The norm of  x, 
*     the transformed residual vector  Pr - RQ'x,  and the constraint
*     values  Ax  are also initialized.
*
*     If the computed point gives a row error of more than the
*     feasibility tolerance, an extra step of iterative refinement is
*     used.  If  x  is still infeasible,  rowerr is set to true.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 31-October-1984.
*     This version of lssetx dated 29-December-1985.
*     ==================================================================
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2

      logical            lsdbg
      parameter         (ldbg = 5)
      common    /lsdebg/ ilsdbg(ldbg), lsdbg

      external           idamax, ddot
      intrinsic          abs, min
      parameter        ( ntry  = 2 )
      parameter        ( zero  = 0.0d+0, one = 1.0d+0 )

*     ------------------------------------------------------------------
*     Move  x  onto the simple bounds in the working set.
*     ------------------------------------------------------------------
      do 100, k = nfree+1, n
          j   = kx(k)
          is  = istate(j)
          bnd = bl(j)
          if (is .ge. 2) bnd  = bu(j)
          if (is .ne. 4) x(j) = bnd
  100 continue

*     ------------------------------------------------------------------
*     Move  x  onto the general constraints in the working set.
*     We shall make  ntry  tries at getting acceptable row errors.
*     ------------------------------------------------------------------
      ktry   = 1
      jmax   = 1
      errmax = zero

*     repeat
  200    if (nactiv .gt. 0) then

*           Set  work = residuals for constraints in the working set.
*           Solve for p, the smallest correction to x that gives a point
*           on the constraints in the working set.  Define  p = Y*(py),
*           where  py  solves the triangular system  T*(py) = residuals.

            do 220, i = 1, nactiv
               k   = kactiv(i)
               j   = n + k
               bnd = bl(j)
               if (istate(j) .eq. 2) bnd = bu(j)
               work(i) = bnd - ddot  ( n, A(k,1), ldA, x, 1 )
  220       continue

            call cmtsol( 1, ldT, nactiv, T(1,nZ+1), work )
            call dload ( n, zero, p, 1 )
            call dcopy ( nactiv, work, 1, p(nZ+1), 1 )

            call cmqmul( 2, n, nZ, nfree, ldQ, unitQ, kx, p, Q, work )
            call daxpy ( n, one, p, 1, x, 1 )
         end if

*        ---------------------------------------------------------------
*        Compute the 2-norm of  x.
*        Initialize  Ax  for all the general constraints.
*        ---------------------------------------------------------------
         xnorm  = dnrm2 ( n, x, 1 )
         if (nclin .gt. 0)
     $      call dgemv ( 'N', nclin, n, one, A, ldA,
     $                   x, 1, zero, Ax, 1 )

*        ---------------------------------------------------------------
*        Check the row residuals.
*        ---------------------------------------------------------------
         if (nactiv .gt. 0) then
            do 300, k = 1, nactiv
               i   = kactiv(k)
               j   = n + i
               is  = istate(j)
               if (is .eq. 1) work(k) = bl(j) - Ax(i)
               if (is .ge. 2) work(k) = bu(j) - Ax(i)
  300       continue

            jmax   = idamax( nactiv, work, 1 )
            errmax = abs( work(jmax) )
         end if

         ktry = ktry + 1
*     until    (errmax .le. featol(jmax) .or. ktry .gt. ntry
      if (.not.(errmax .le. featol(jmax) .or. ktry .gt. ntry)) go to 200

      rowerr = errmax .gt. featol(jmax)

*     ==================================================================
*     Compute the linear objective value  c'x  and the transformed
*     residual  Pr  -  RQ'x = res0  -  RQ'x.
*     ==================================================================
      if (nrank .gt. 0  .or.  linObj) then
         call dcopy ( n, x, 1, p, 1 )
         call cmqmul( 6, n, nZ, nfree, ldQ, unitQ, kx, p, Q, work )
      end if

      ctx = zero
      if (linObj)
     $   ctx = ddot  ( n, cq, 1, p, 1 )

      if (nrank .gt. 0) then

         call dtrmv ( 'U', 'N', 'N', nrank, R, ldR, p, 1 )
         if (nrank .lt. n)
     $      call dgemv ( 'N', nrank, n-nrank, one, R(1,nrank+1), ldR,
     $                   p(nrank+1), 1, one, p, 1 )

         call dcopy ( nrank,         res0, 1, res, 1 )
         call daxpy ( nrank, (-one), p   , 1, res, 1 )

      end if

*      if (lsdbg  .and.  ilsdbg(2) .gt. 0)
*     $   write(iPrint, 2200) (x(j), j = 1, n)

      return

 2200 format(/ ' //lssetx// Variables after refinement ... '/ (5g12.3))

*     end of lssetx
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     File  OPSUBS FORTRAN
*
*     OPFILE   OPLOOK   OPNUMB   OPSCAN   OPTOKN   OPUPPR
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE OPFILE( IOPTNS, NOUT, INFORM, OPKEY )
      INTEGER            IOPTNS, NOUT, INFORM
      EXTERNAL           OPKEY

************************************************************************
*     OPFILE  reads the options file from unit  IOPTNS  and loads the
*     options into the relevant elements of the integer and real
*     parameter arrays.
*
*     Systems Optimization Laboratory, Stanford University.
*     This version dated December 18, 1985.
************************************************************************
      LOGICAL             PRNT
      CHARACTER*16        KEY   , TOKEN(1)
      CHARACTER*72        BUFFER, OLDBUF

      PRNT   = .TRUE.

*     Return if the unit number is out of range.

      IF (IOPTNS .LT. 0  .OR.  IOPTNS .GT. 99) THEN
         INFORM = 1
         RETURN
      END IF

*     ------------------------------------------------------------------
*     Look for  BEGIN, ENDRUN  or  SKIP.
*     ------------------------------------------------------------------
      NREAD  = 0
   50    CONTINUE
*        READ (IOPTNS, '(A)', END = 930) BUFFER
         NREAD = NREAD + 1
         NKEY  = 1
         CALL OPTOKN( BUFFER, NKEY, TOKEN )
         KEY   = TOKEN(1)
         IF (KEY .EQ. 'ENDRUN') GO TO 940
*         IF (KEY .NE. 'BEGIN' ) THEN
*            IF (NREAD .EQ. 1  .AND.  KEY .NE. 'SKIP') THEN
*               WRITE (NOUT, 2000) IOPTNS, BUFFER
*            END IF
*           GO TO 50
*         END IF

*     ------------------------------------------------------------------
*     BEGIN found.
*     This is taken to be the first line of an OPTIONS file.
*     Read the second line to see if it is NOLIST.
*     ------------------------------------------------------------------
      OLDBUF = BUFFER
*      READ (IOPTNS, '(A)', END = 920) BUFFER

      CALL OPKEY ( NOUT, BUFFER, KEY )

      IF (KEY .EQ. 'NOLIST') THEN
         PRNT   = .FALSE.
      END IF

*      IF (PRNT) THEN
*         WRITE (NOUT, '(// A / A /)')
*     $      ' OPTIONS file',
*     $      ' ------------'
*         WRITE (NOUT, '(6X, A )') OLDBUF, BUFFER
*      END IF

*     ------------------------------------------------------------------
*     Read the rest of the file.
*     ------------------------------------------------------------------
*+    while (key .ne. 'end') loop
  100 IF    (KEY .NE. 'END') THEN
*         READ (IOPTNS, '(A)', END = 920) BUFFER
*        IF (PRNT)
*     $      WRITE (NOUT, '( 6X, A )') BUFFER

         CALL OPKEY ( NOUT, BUFFER, KEY )

         IF (KEY .EQ.   'LIST') PRNT = .TRUE.
         IF (KEY .EQ. 'NOLIST') PRNT = .FALSE.
         GO TO 100
      END IF
*+    end while

      INFORM =  0
      RETURN

*  920 WRITE (NOUT, 2200) IOPTNS
      INFORM = 2
      RETURN

*  930 WRITE (NOUT, 2300) IOPTNS
      INFORM = 3
      RETURN
940   CONTINUE
*   940 WRITE (NOUT, '(// 6X, A)') BUFFER
      INFORM = 4
      RETURN

 2000 FORMAT(
     $ //' XXX  Error while looking for an OPTIONS file on unit', I7
     $ / ' XXX  The file should start with BEGIN, SKIP or ENDRUN'
     $ / ' XXX  but the first record found was the following:'
     $ //' ---->', A
     $ //' XXX  Continuing to look for OPTIONS file...')
 2200 FORMAT(//' XXX  End-of-file encountered while processing',
     $         ' an OPTIONS file on unit', I6)
 2300 FORMAT(//' XXX  End-of-file encountered while looking for',
     $         ' an OPTIONS file on unit', I6)

*     End of  OPFILE.

      END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE OPLOOK (NDICT, DICTRY, ALPHA, KEY, ENTRY)
C
C
C Description and usage:
C
C       Performs dictionary lookups.  A pointer is returned if a
C    match is found between the input key and the corresponding
C    initial characters of one of the elements of the dictionary.
C    If a "synonym" has been provided for an entry, the search is
C    continued until a match to a primary dictionary entry is found.
C    Cases of no match, or multiple matches, are also provided for.
C
C     Dictionary entries must be left-justified, and may be alphabetized
C    for faster searches.  Secondary entries, if any, are composed of
C    two words separated by one or more characters such as blank, tab,
C    comma, colon, or equal sign which are treated as non-significant
C    by OPSCAN.  The first entry of each such pair serves as a synonym
C    for the second, more fundamental keyword.
C
C       The ordered search stops after the section of the dictionary
C    having the same first letters as the key has been checked, or
C    after a specified number of entries have been examined.  A special
C    dictionary entry, the vertical bar '|', will also terminate the
C    search.  This will speed things up if an appropriate dictionary
C    length parameter cannot be determined.  Both types of search are
C    sequential.  See "Notes" below for some suggestions if efficiency
C    is an issue.
C
C
C Parameters:
C
C    Name    Dimension  Type  I/O/S  Description
C    NDICT               I    I      Number of dictionary entries to be
C                                    examined.
C    DICTRY  NDICT       C    I      Array of dictionary entries,
C                                    left-justified in their fields.
C                                    May be alphabetized for efficiency,
C                                    in which case ALPHA should be
C                                    .TRUE.  Entries with synonyms are
C                                    of the form
C                                    'ENTRY : SYNONYM', where 'SYNONYM'
C                                    is a more fundamental entry in the
C                                    same dictionary.  NOTE: Don't build
C                                    "circular" dictionaries!
C    ALPHA               L    I      Indicates whether the dictionary
C                                    is in alphabetical order, in which
C                                    case the search can be terminated
C                                    sooner.
C    KEY                 C    I/O    String to be compared against the
C                                    dictionary.  Abbreviations are OK
C                                    if they correspond to a unique
C                                    entry in the dictionary.  KEY is
C                                    replaced on termination by its most
C                                    fundamental equivalent dictionary
C                                    entry (uppercase, left-justified)
C                                    if a match was found.
C    ENTRY               I      O    Dictionary pointer.  If > 0, it
C                                    indicates which entry matched KEY.
C                                    In case of trouble, a negative
C                                    value means that a UNIQUE match
C                                    was not found - the absolute value
C                                    of ENTRY points to the second
C                                    dictionary entry that matched KEY.
C                                    Zero means that NO match could be
C                                    found.  ENTRY always refers to the
C                                    last search performed -
C                                    in searching a chain of synonyms,
C                                    a non-positive value will be
C                                    returned if there is any break,
C                                    even if the original input key
C                                    was found.
C
C
C External references:
C
C    Name    Description
C    OPSCAN  Finds first and last significant characters.
C
C
C Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C               Appears to satisfy the ANSI Fortran 77 standard.
C
C
C Notes:
C
C    (1)  IMPLICIT NONE is non-standard.  (Has been commented out.)
C
C    (2)  We have assumed that the dictionary is not too big.  If
C         many searches are to be done or if the dictionary has more
C         than a dozen or so entries, it may be advantageous to build
C         an index array of pointers to the beginning of the section
C         of the dictionary containing each letter, then pass in the
C         portion of the dictionary beginning with DICTRY (INDEX).
C         (This won't generally work for dictionaries with synonyms.)
C         For very large problems, a completely different approach may
C         be advisable, e.g. a binary search for ordered dictionaries.
C
C    (3)  OPLOOK is case sensitive.  In most applications it will be
C         necessary to use an uppercase dictionary, and to convert the
C         input key to uppercase before calling OPLOOK.  Companion
C         routines OPTOKN and PAIRS, available from the author, already
C         take care of this.
C
C    (4)  The key need not be left-justified.  Any leading (or
C         trailing) characters which are "non-significant" to OPSCAN
C         will be ignored.  These include blanks, horizontal tabs,
C         commas, colons, and equal signs.  See OPSCAN for details.
C
C    (5)  The ASCII collating sequence for character data is assumed.
C         (N.B. This means the numerals precede the alphabet, unlike
C         common practice!)  This should not cause trouble on EBCDIC
C         machines if DICTRY just contains alphabetic keywords.
C         Otherwise it may be necessary to use the FORTRAN lexical
C         library routines to force use of the ASCII sequence.
C
C    (6)  Parameter NUMSIG sets a limit on the length of significant
C         dictionary entries.  Special applications may require that
C         this be increased.  (It is 16 in the present version.)
C
C    (7)  No protection against "circular" dictionaries is provided:
C         don't claim that A is B, and that B is A.  All synonym chains
C         must terminate!  Other potential errors not checked for
C         include duplicate or mis-ordered entries.
C
C    (8)  The handling of ambiguities introduces some ambiguity:
C
C            ALPHA = .TRUE.  A potential problem, when one entry
C                            looks like an abbreviation for another
C                            (eg. does 'A' match 'A' or 'AB'?) was
C                            resolved by dropping out of the search
C                            immediately when an "exact" match is found.
C
C            ALPHA = .FALSE. The programmer must ensure that the above
C                            situation does not arise: each dictionary
C                            entry must be recognizable, at least when
C                            specified to full length.  Otherwise, the
C                            result of a search will depend on the
C                            order of entries.
C
C
C Author:  Robert Kennelly, Informatics General Corporation.
C
C
C Development history:
C
C    24 Feb. 1984  RAK/DAS  Initial design and coding.
C    25 Feb. 1984    RAK    Combined the two searches by suitable
C                           choice of terminator FLAG.
C    28 Feb. 1984    RAK    Optional synonyms in dictionary, no
C                           longer update KEY.
C    29 Mar. 1984    RAK    Put back replacement of KEY by its
C                           corresponding entry.
C    21 June 1984    RAK    Corrected bug in error handling for cases
C                           where no match was found.
C    23 Apr. 1985    RAK    Introduced test for exact matches, which
C                           permits use of dictionary entries which
C                           would appear to be ambiguous (for ordered
C                           case).  Return -I to point to the entry
C                           which appeared ambiguous (had been -1).
C                           Repaired loop termination - had to use
C                           equal length strings or risk quitting too
C                           soon when one entry is an abbreviation
C                           for another.  Eliminated HIT, reduced
C                           NUMSIG to 16.
C    15 Nov. 1985    MAS    Loop 20 now tests .LT. FLAG, not .LE. FLAG.
C                           If ALPHA is false, FLAG is now '|', not '{'.
C    26 Jan. 1986    PEG    Declaration of FLAG and TARGET modified to
C                           conform to ANSI-77 standard.
C-----------------------------------------------------------------------


C     Variable declarations.
C     ----------------------

*     IMPLICIT NONE

C     Parameters.

      INTEGER
     $   NUMSIG
      CHARACTER
     $   BLANK, VBAR
      PARAMETER
     $   (BLANK = ' ', VBAR = '|', NUMSIG = 16)

C     Variables.

      LOGICAL
     $   ALPHA
      INTEGER
     $   ENTRY, FIRST, I, LAST, LENGTH, MARK, NDICT
*     CHARACTER
*    $   DICTRY (NDICT) * (*), FLAG * (NUMSIG),
*    $   KEY * (*), TARGET * (NUMSIG)
      CHARACTER
     $   DICTRY (NDICT) * (*), FLAG * 16,
     $   KEY * (*), TARGET * 16

C     Procedures.

      EXTERNAL
     $   OPSCAN


C     Executable statements.
C     ----------------------

      ENTRY = 0

C     Isolate the significant portion of the input key (if any).

      FIRST = 1
      LAST  = MIN( LEN(KEY), NUMSIG )
      CALL OPSCAN (KEY, FIRST, LAST, MARK)

      IF (MARK .GT. 0) THEN
         TARGET = KEY (FIRST:MARK)

C        Look up TARGET in the dictionary.

   10    CONTINUE
            LENGTH = MARK - FIRST + 1

C           Select search strategy by cunning choice of termination test
C           flag.  The vertical bar is just about last in both the
C           ASCII and EBCDIC collating sequences.

            IF (ALPHA) THEN
               FLAG = TARGET
            ELSE
               FLAG = VBAR
            END IF


C           Perform search.
C           ---------------

            I = 0
   20       CONTINUE
               I = I + 1
               IF (TARGET (1:LENGTH) .EQ. DICTRY (I) (1:LENGTH)) THEN
                  IF (ENTRY .EQ. 0) THEN

C                    First "hit" - must still guard against ambiguities
C                    by searching until we've gone beyond the key
C                    (ordered dictionary) or until the end-of-dictionary
C                    mark is reached (exhaustive search).

                     ENTRY = I

C                    Special handling if match is exact - terminate
C                    search.  We thus avoid confusion if one dictionary
C                    entry looks like an abbreviation of another.
C                    This fix won't generally work for un-ordered
C                    dictionaries!

                     FIRST = 1
                     LAST = NUMSIG
                     CALL OPSCAN (DICTRY (ENTRY), FIRST, LAST, MARK)
                     IF (MARK .EQ. LENGTH) I = NDICT
                  ELSE


C                    Oops - two hits!  Abnormal termination.
C                    ---------------------------------------

                     ENTRY = -I
                     RETURN
                  END IF
               END IF

C           Check whether we've gone past the appropriate section of the
C           dictionary.  The test on the index provides insurance and an
C           optional means for limiting the extent of the search.

            IF (DICTRY (I) (1:LENGTH) .LT. FLAG  .AND.  I .LT. NDICT)
     $         GO TO 20


C           Check for a synonym.
C           --------------------

            IF (ENTRY .GT. 0) THEN

C              Look for a second entry "behind" the first entry.  FIRST
C              and MARK were determined above when the hit was detected.

               FIRST = MARK + 2
               CALL OPSCAN (DICTRY (ENTRY), FIRST, LAST, MARK)
               IF (MARK .GT. 0) THEN

C                 Re-set target and dictionary pointer, then repeat the
C                 search for the synonym instead of the original key.

                  TARGET = DICTRY (ENTRY) (FIRST:MARK)
                  ENTRY = 0
                  GO TO 10

               END IF
            END IF

      END IF
      IF (ENTRY .GT. 0) KEY = DICTRY (ENTRY)


C     Normal termination.
C     -------------------

      RETURN

C     End of OPLOOK
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      FUNCTION OPNUMB( STRING )

      LOGICAL          OPNUMB
      CHARACTER*(*)    STRING

************************************************************************
*     Description and usage:
*
*        A simple(-minded) test for numeric data is implemented by
*        searching an input string for legitimate characters:
*                digits 0 to 9, D, E, -, + and .
*        Insurance is provided by requiring that a numeric string
*        have at least one digit, at most one D, E or .
*        and at most two -s or +s.  Note that a few ambiguities remain:
*
*           (a)  A string might have the form of numeric data but be
*                intended as text.  No general test can hope to detect
*                such cases.
*
*           (b)  There is no check for correctness of the data format.
*                For example a meaningless string such as 'E1.+2-'
*                will be accepted as numeric.
*
*        Despite these weaknesses, the method should work in the
*        majority of cases.
*
*
*     Parameters:
*
*        Name    Dimension  Type  I/O/S  Description
*        OPNUMB              L      O    Set .TRUE. if STRING appears
*                                        to be numerical data.
*        STRING              C    I      Input data to be tested.
*
*
*     Environment:  ANSI FORTRAN 77.
*
*
*     Notes:
*
*        (1)  It is assumed that STRING is a token extracted by
*             OPTOKN, which will have converted any lower-case
*             characters to upper-case.
*
*        (2)  OPTOKN pads STRING with blanks, so that a genuine
*             number is of the form  '1234        '.
*             Hence, the scan of STRING stops at the first blank.
*
*        (3)  COMPLEX data with parentheses will not look numeric.
*
*
*     Systems Optimization Laboratory, Stanford University.
*     12 Nov  1985    Initial design and coding, starting from the
*                     routine ALPHA from Informatics General, Inc.
************************************************************************

      LOGICAL         NUMBER
      INTEGER         J, LENGTH, NDIGIT, NEXP, NMINUS, NPLUS, NPOINT
      CHARACTER*1     ATOM

      NDIGIT = 0
      NEXP   = 0
      NMINUS = 0
      NPLUS  = 0
      NPOINT = 0
      NUMBER = .TRUE.
      LENGTH = LEN (STRING)
      J      = 0

   10    J    = J + 1
         ATOM = STRING (J:J)
         IF      (ATOM .GE. '0'  .AND.  ATOM .LE. '9') THEN
            NDIGIT = NDIGIT + 1
         ELSE IF (ATOM .EQ. 'D'  .OR.   ATOM .EQ. 'E') THEN
            NEXP   = NEXP   + 1
         ELSE IF (ATOM .EQ. '-') THEN
            NMINUS = NMINUS + 1
         ELSE IF (ATOM .EQ. '+') THEN
            NPLUS  = NPLUS  + 1
         ELSE IF (ATOM .EQ. '.') THEN
            NPOINT = NPOINT + 1
         ELSE IF (ATOM .EQ. ' ') THEN
            J      = LENGTH
         ELSE
            NUMBER = .FALSE.
         END IF

         IF (NUMBER  .AND.  J .LT. LENGTH) GO TO 10

      OPNUMB = NUMBER
     $         .AND.  NDIGIT .GE. 1
     $         .AND.  NEXP   .LE. 1
     $         .AND.  NMINUS .LE. 2
     $         .AND.  NPLUS  .LE. 2
     $         .AND.  NPOINT .LE. 1

      RETURN

*     End of OPNUMB
      END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE OPSCAN (STRING, FIRST, LAST, MARK)
C
C
C Description and usage:
C
C       Looks for non-blank fields ("tokens") in a string, where the
C    fields are of arbitrary length, separated by blanks, tabs, commas,
C    colons, or equal signs.  The position of the end of the 1st token
C    is also returned, so this routine may be conveniently used within
C    a loop to process an entire line of text.
C
C       The procedure examines a substring, STRING (FIRST : LAST), which
C    may of course be the entire string (in which case just call OPSCAN
C    with FIRST <= 1 and LAST >= LEN (STRING) ).  The indices returned
C    are relative to STRING itself, not the substring.
C
C
C Parameters:
C
C    Name    Dimension  Type  I/O/S  Description
C    STRING              C    I      Text string containing data to be
C                                    scanned.
C    FIRST               I    I/O    Index of beginning of substring.
C                                    If <= 1, the search begins with 1.
C                                    Output is index of beginning of
C                                    first non-blank field, or 0 if no
C                                    token was found.
C    LAST                I    I/O    Index of end of substring.
C                                    If >= LEN (STRING), the search
C                                    begins with LEN (STRING).  Output
C                                    is index of end of last non-blank
C                                    field, or 0 if no token was found.
C    MARK                I      O    Points to end of first non-blank
C                                    field in the specified substring.
C                                    Set to 0 if no token was found.
C
C
C Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C               ANSI Fortran 77, except for the tab character HT.
C
C Notes:
C
C    (1)  IMPLICIT NONE is non-standard.  Constant HT (Tab) is defined
C         in a non-standard way:  the CHAR function is not permitted
C         in a PARAMETER declaration (OK on VAX, though).  For Absoft
C         FORTRAN 77 on 68000 machines, use HT = 9.  In other cases, it
C         may be best to declare HT as a variable and assign
C         HT = CHAR(9) on ASCII machines, or CHAR(5) for EBCDIC.
C
C    (2)  The pseudo-recursive structure was chosen for fun.  It is
C         equivalent to three DO loops with embedded GO TOs in sequence.
C
C    (3)  The variety of separators recognized limits the usefulness of
C         this routine somewhat.  The intent is to facilitate handling
C         such tokens as keywords or numerical values.  In other
C         applications, it may be necessary for ALL printing characters
C         to be significant.  A simple modification to statement
C         function SOLID will do the trick.
C
C
C Author:  Robert Kennelly, Informatics General Corporation.
C
C
C Development history:
C
C    29 Dec. 1984    RAK    Initial design and coding, (very) loosely
C                           based on SCAN_STRING by Ralph Carmichael.
C    25 Feb. 1984    RAK    Added ':' and '=' to list of separators.
C    16 Apr. 1985    RAK    Defined SOLID in terms of variable DUMMY
C                           (previous re-use of STRING was ambiguous).
C
C-----------------------------------------------------------------------


C     Variable declarations.
C     ----------------------

*     IMPLICIT NONE

C     Parameters.

      CHARACTER
     $   BLANK, EQUAL, COLON, COMMA, HT
      PARAMETER
     $   (BLANK = ' ', EQUAL = '=', COLON = ':', COMMA = ',')

C     Variables.

      LOGICAL
     $   SOLID
      INTEGER
     $   BEGIN, END, FIRST, LAST, LENGTH, MARK
      CHARACTER
     $   DUMMY, STRING * (*)

C     Statement functions.

      SOLID (DUMMY) = (DUMMY .NE. BLANK) .AND.
     $                (DUMMY .NE. COLON) .AND.
     $                (DUMMY .NE. COMMA) .AND.
     $                (DUMMY .NE. EQUAL) .AND.
     $                (DUMMY .NE. HT)


C     Executable statements.
C     ----------------------

****  HT     = CHAR(9) for ASCII machines, CHAR(5) for EBCDIC.
      HT     = CHAR(9)
      MARK   = 0
      LENGTH = LEN (STRING)
      BEGIN  = MAX (FIRST, 1)
      END    = MIN (LENGTH, LAST)

C     Find the first significant character ...

      DO 30 FIRST = BEGIN, END, +1
         IF (SOLID (STRING (FIRST : FIRST))) THEN

C           ... then the end of the first token ...

            DO 20 MARK = FIRST, END - 1, +1
               IF (.NOT.SOLID (STRING (MARK + 1 : MARK + 1))) THEN

C                 ... and finally the last significant character.

                  DO 10 LAST = END, MARK, -1
                     IF (SOLID (STRING (LAST : LAST))) THEN
                        RETURN
                     END IF
   10             CONTINUE

C                 Everything past the first token was a separator.

                  LAST = LAST + 1
                  RETURN
               END IF
   20       CONTINUE

C           There was nothing past the first token.

            LAST = MARK
            RETURN
         END IF
   30 CONTINUE

C     Whoops - the entire substring STRING (BEGIN : END) was composed of
C     separators !

      FIRST = 0
      MARK = 0
      LAST = 0
      RETURN

C     End of OPSCAN
      END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE OPTOKN (STRING, NUMBER, LIST)
C
C
C Description and usage:
C
C       An aid to parsing input data.  The individual "tokens" in a
C    character string are isolated, converted to uppercase, and stored
C    in an array.  Here, a token is a group of significant, contiguous
C    characters.  The following are NON-significant, and hence may
C    serve as separators:  blanks, horizontal tabs, commas, colons,
C    and equal signs.  See OPSCAN for details.  Processing continues
C    until the requested number of tokens have been found or the end
C    of the input string is reached.
C
C
C Parameters:
C
C    Name    Dimension  Type  I/O/S  Description
C    STRING              C    I      Input string to be analyzed.
C    NUMBER              I    I/O    Number of tokens requested (input)
C                                    and found (output).
C    LIST    NUMBER      C      O    Array of tokens, changed to upper
C                                    case.
C
C
C External references:
C
C    Name    Description
C    OPSCAN  Finds positions of first and last significant characters.
C    OPUPPR  Converts a string to uppercase.
C
C
C Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C               Appears to satisfy the ANSI Fortran 77 standard.
C
C
C Notes:
C
C    (1)  IMPLICIT NONE is non-standard.  (Has been commented out.)
C
C
C Author:  Robert Kennelly, Informatics General Corporation.
C
C
C Development history:
C
C    16 Jan. 1984    RAK    Initial design and coding.
C    16 Mar. 1984    RAK    Revised header to reflect full list of
C                           separators, repaired faulty WHILE clause
C                           in "10" loop.
C    18 Sep. 1984    RAK    Change elements of LIST to uppercase one
C                           at a time, leaving STRING unchanged.
C
C-----------------------------------------------------------------------


C     Variable declarations.
C     ----------------------

*     IMPLICIT NONE

C     Parameters.

      CHARACTER
     $   BLANK
      PARAMETER
     $   (BLANK = ' ')

C     Variables.

      INTEGER
     $   COUNT, FIRST, I, LAST, MARK, NUMBER
      CHARACTER
     $   STRING * (*), LIST (NUMBER) * (*)

C     Procedures.

      EXTERNAL
     $   OPUPPR, OPSCAN


C     Executable statements.
C     ----------------------

C     WHILE there are tokens to find, loop UNTIL enough have been found.

      FIRST = 1
      LAST = LEN (STRING)

      COUNT = 0
   10 CONTINUE

C        Get delimiting indices of next token, if any.

         CALL OPSCAN (STRING, FIRST, LAST, MARK)
         IF (LAST .GT. 0) THEN
            COUNT = COUNT + 1

C           Pass token to output string array, then change case.

            LIST (COUNT) = STRING (FIRST : MARK)
            CALL OPUPPR (LIST (COUNT))
            FIRST = MARK + 2
            IF (COUNT .LT. NUMBER) GO TO 10

         END IF


C     Fill the rest of LIST with blanks and set NUMBER for output.

      DO 20 I = COUNT + 1, NUMBER
         LIST (I) = BLANK
   20 CONTINUE

      NUMBER = COUNT


C     Termination.
C     ------------

      RETURN

C     End of OPTOKN
      END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE OPUPPR(STRING)
C
C ACRONYM:  UPper CASE
C
C PURPOSE:  This subroutine changes all lower case letters in the
C           character string to upper case.
C
C METHOD:   Each character in STRING is treated in turn.  The intrinsic
C           function INDEX effectively allows a table lookup, with
C           the local strings LOW and UPP acting as two tables.
C           This method avoids the use of CHAR and ICHAR, which appear
C           be different on ASCII and EBCDIC machines.
C
C ARGUMENTS
C    ARG       DIM     TYPE I/O/S DESCRIPTION
C  STRING       *       C   I/O   Character string possibly containing
C                                 some lower-case letters on input;
C                                 strictly upper-case letters on output
C                                 with no change to any non-alphabetic
C                                 characters.
C
C EXTERNAL REFERENCES:
C  LEN    - Returns the declared length of a CHARACTER variable.
C  INDEX  - Returns the position of second string within first.
C
C ENVIRONMENT:  ANSI FORTRAN 77
C
C DEVELOPMENT HISTORY:
C     DATE  INITIALS  DESCRIPTION
C   06/28/83   CLH    Initial design.
C   01/03/84   RAK    Eliminated NCHAR input.
C   06/14/84   RAK    Used integer PARAMETERs in comparison.
C   04/21/85   RAK    Eliminated DO/END DO in favor of standard code.
C   09/10/85   MAS    Eliminated CHAR,ICHAR in favor of LOW, UPP, INDEX.
C
C AUTHOR: Charles Hooper, Informatics General, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      CHARACTER      STRING * (*)
      INTEGER        I, J
      CHARACTER      C*1, LOW*26, UPP*26
      DATA           LOW /'abcdefghijklmnopqrstuvwxyz'/,
     $               UPP /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

      DO 10 J = 1, LEN(STRING)
         C    = STRING(J:J)
         IF (C .GE. 'a'  .AND.  C .LE. 'z') THEN
            I           = INDEX( LOW, C )
            IF (I .GT. 0) STRING(J:J) = UPP(I:I)
         END IF
   10 CONTINUE
      RETURN
     
*     End of OPUPPR

      END 
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     File  cmsubs.f
*
*     cmalf1   cmalf    cmchk    cmperm   cmprt   +cmqmul  +cmr1md
*    +cmrswp   cmtsol
*   (+) Also in cmoptsubs.f
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine cmalf1( firstv, negstp, bigalf, bigbnd, pnorm,
     $                   jadd1 , jadd2 , palfa1, palfa2,
     $                   istate, n, nctotl,
     $                   Anorm, Ap, Ax, bl, bu, featol, p, x )

      implicit           double precision(a-h,o-z)
      logical            firstv, negstp
      integer            istate(nctotl)
      double precision   Anorm(*), Ap(*), Ax(*)
      double precision   bl(nctotl), bu(nctotl), featol(nctotl),
     $                   p(n), x(n)

*     ==================================================================
*     cmalf1  finds steps palfa1, palfa2 such that
*        x + palfa1*p  reaches a linear constraint that is currently not
*                      in the working set but is satisfied.
*        x + palfa2*p  reaches a linear constraint that is currently not
*                      in the working set but is violated.
*     The constraints are perturbed by an amount featol, so that palfa1
*     is slightly larger than it should be,  and palfa2 is slightly
*     smaller than it should be.  This gives some leeway later when the
*     exact steps are computed by cmalf.
*
*     Constraints in the working set are ignored  (istate(j) .ge. 1).
*
*     If negstp is true, the search direction will be taken to be  - p.
*
*
*     Values of istate(j)....
*
*        - 2         - 1         0           1          2         3
*     a'x lt bl   a'x gt bu   a'x free   a'x = bl   a'x = bu   bl = bu
*
*     The values  -2  and  -1  do not occur once a feasible point has
*     been found.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original Fortran 66 version written  May 1980.
*     This version of cmalf1 dated 27-Oct-92.
*     ==================================================================
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9

      logical            cmdbg
      integer            lcmdbg
      parameter         (lcmdbg = 5)
      common    /cmdebg/ icmdbg(lcmdbg), cmdbg

      logical            lastv
      intrinsic          abs
      parameter        ( zero = 0.0d+0, one = 1.0d+0 )

*      if (cmdbg  .and.  icmdbg(3) .gt. 0) write (iPrint, 1100)
      lastv  = .not. firstv
      jadd1  = 0
      jadd2  = 0
      palfa1 = bigalf

      palfa2 = zero
      if (firstv) palfa2 = bigalf

      do 200 j = 1, nctotl
         js = istate(j)
         if (js .le. 0) then
            if (j .le. n) then
               atx    = x(j)
               atp    = p(j)
               rownrm = one
            else
               i      = j - n
               atx    = Ax(i)
               atp    = Ap(i)
               rownrm = one  +  Anorm(i)
            end if
            if (negstp) atp = - atp

            if ( abs( atp ) .le. epspt9*rownrm*pnorm) then

*              This constraint appears to be constant along P.  It is
*              not used to compute the step.  Give the residual a value
*              that can be spotted in the debug output.

               res = - one
            else if (atp .le. zero  .and.  js .ne. -2) then
*              ---------------------------------------------------------
*              a'x  is decreasing and the lower bound is not violated.
*              ---------------------------------------------------------
*              First test for smaller PALFA1.

               absatp = - atp
               if (bl(j) .gt. (-bigbnd)) then
                  res    = atx - bl(j) + featol(j)
                  if (bigalf*absatp .gt. abs( res )) then
                     if (palfa1*absatp .gt. res)  then
                        palfa1 = res / absatp
                        jadd1  = j
                     end if
                  end if
               end if

               if (js .eq. -1) then

*                 The upper bound is violated.  Test for either larger
*                 or smaller PALFA2, depending on the value of FIRSTV.

                  res    = atx - bu(j) - featol(j)
                  if (bigalf*absatp .gt. abs( res )) then
                     if (firstv  .and.  palfa2*absatp .gt. res  .or.
     $                    lastv  .and.  palfa2*absatp .lt. res) then
                        palfa2 = res / absatp
                        jadd2  = j
                     end if
                  end if
               end if
            else if (atp .gt. zero  .and.  js .ne. -1) then
*              ---------------------------------------------------------
*              a'x  is increasing and the upper bound is not violated.
*              ---------------------------------------------------------
*              Test for smaller PALFA1.

               if (bu(j) .lt. bigbnd) then
                  res = bu(j) - atx + featol(j)
                  if (bigalf*atp .gt. abs( res )) then
                     if (palfa1*atp .gt. res) then
                        palfa1 = res / atp
                        jadd1  = j
                     end if
                  end if
               end if

               if (js .eq. -2) then

*                 The lower bound is violated.  Test for a new PALFA2.

                  res  = bl(j) - atx - featol(j)
                  if (bigalf*atp .gt. abs( res )) then
                     if (firstv  .and.  palfa2*atp .gt. res  .or.
     $                    lastv  .and.  palfa2*atp .lt. res) then
                        palfa2 = res / atp
                        jadd2  = j
                     end if
                  end if
               end if
            end if

*            if (cmdbg  .and.  icmdbg(3) .gt. 0)
*     $         write (iPrint, 1200) j, js, featol(j), res,
*     $                            atp, jadd1, palfa1, jadd2, palfa2
         end if
  200 continue

      return

 1100 format(/ '    j  js         featol        res             Ap',
     $         '     jadd1       palfa1     jadd2       palfa2' /)
 1200 format(i5, i4, 3g15.5, 2(i6, g17.7))

*     end of cmalf1
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine cmalf ( firstv, hitlow, 
     $                   istate, inform, jadd, n, nctotl, numinf,
     $                   alfa, palfa, atphit, 
     $                   bigalf, bigbnd, pnorm,
     $                   Anorm, Ap, Ax, bl, bu, 
     $                   featol, p, x )

      implicit           double precision(a-h,o-z)
      integer            istate(nctotl)
      double precision   Anorm(*), Ap(*), Ax(*),
     $                   bl(nctotl), bu(nctotl), featol(nctotl),
     $                   p(n), x(n)
      logical            firstv, hitlow

*     ==================================================================
*     cmalf   finds a step alfa such that the point x + alfa*p reaches
*     one of the linear constraints (including bounds).  Two possible
*     steps are defined as follows...
*
*     alfa1   is the maximum step that can be taken without violating
*             one of the linear constraints that is currently satisfied.
*     alfa2   reaches a linear constraint that is currently violated.
*             Usually this will be the furthest such constraint along p,
*             but if firstv = .true. it will be the first one along p.
*             This is used only when the problem has been determined to 
*             be infeasible, and the sum of infeasibilities are being
*             minimized.  (alfa2  is not defined if numinf = 0.)
*
*     alfa will usually be the minimum of alfa1 and alfa2.
*     alfa could be negative (since we allow inactive constraints
*     to be violated by as much as featol).  In such cases, a
*     third possible step is computed, to find the nearest satisfied
*     constraint (perturbed by featol) along the direction  - p.
*     alfa  will be reset to this step if it is shorter.  This is the
*     only case for which the final step  alfa  does not move x exactly
*     onto a constraint (the one denoted by jadd).
*
*     Constraints in the working set are ignored  (istate(j) ge 1).
*
*     jadd    denotes which linear constraint is reached.
*
*     hitlow  indicates whether it is the lower or upper bound that
*             has restricted alfa.
*
*     Values of istate(j)....
*
*        - 2         - 1         0           1          2         3
*     a'x lt bl   a'x gt bu   a'x free   a'x = bl   a'x = bu   bl = bu
*
*     The values -2 and -1 do not occur once a feasible point has been
*     found.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original Fortran 66 version written  May 1980.
*     This version of  cmalf  dated  28-Oct-92.
*     ==================================================================
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9

      logical            cmdbg
      integer            lcmdbg
      parameter         (lcmdbg = 5)
      common    /cmdebg/ icmdbg(lcmdbg), cmdbg

      logical            hlow1, hlow2, lastv, negstp, step2
      intrinsic          abs, min
      parameter        ( zero = 0.0d+0, one = 1.0d+0 )

      inform = 0

*     ------------------------------------------------------------------
*     First pass -- find steps to perturbed constraints, so that
*     palfa1 will be slightly larger than the true step, and
*     palfa2 will be slightly smaller than it should be.
*     In degenerate cases, this strategy gives us some freedom in the
*     second pass.  The general idea follows that described by P.M.J.
*     Harris, p.21 of Mathematical Programming 5, 1 (1973), 1--28.
*     ------------------------------------------------------------------

      negstp = .false.
      call cmalf1( firstv, negstp, bigalf, bigbnd, pnorm,
     $             jadd1, jadd2, palfa1, palfa2,
     $             istate, n, nctotl,
     $             Anorm, Ap, Ax, bl, bu, featol, p, x )

      jsave1 = jadd1
      jsave2 = jadd2

*     ------------------------------------------------------------------
*     Second pass -- recompute step-lengths without perturbation.
*     Amongst constraints that are less than the perturbed steps,
*     choose the one (of each type) that makes the largest angle
*     with the search direction.
*     ------------------------------------------------------------------
*      if (cmdbg  .and.  icmdbg(3) .gt. 0) write (iPrint, 1000)
      alfa1  = bigalf
      alfa2  = zero
      if (firstv) alfa2 = bigalf

      apmax1 = zero
      apmax2 = zero
      atp1   = zero
      atp2   = zero
      hlow1  = .false.
      hlow2  = .false.
      lastv  = .not. firstv

      do 400 j = 1, nctotl
         js = istate(j)
         if (js .le. 0) then
            if (j  .le. n)  then
               atx    = x(j)
               atp    = p(j)
               rownrm = one
            else
               i      = j - n
               atx    = Ax(i)
               atp    = Ap(i)
               rownrm = Anorm(i) + one
            end if

            if ( abs( atp ) .le. epspt9*rownrm*pnorm) then

*              This constraint appears to be constant along p.  It is
*              not used to compute the step.  Give the residual a value
*              that can be spotted in the debug output.

               res = - one
            else if (atp .le. zero  .and.  js .ne. -2) then
*              ---------------------------------------------------------
*              a'x  is decreasing.
*              ---------------------------------------------------------
*              The lower bound is satisfied.  Test for smaller alfa1.

               absatp = - atp
               if (bl(j) .gt. (-bigbnd)) then
                  res    = atx - bl(j)
                  if (palfa1*absatp .ge. res  .or.  j .eq. jsave1) then
                     if (apmax1*rownrm*pnorm .lt. absatp) then
                        apmax1 = absatp / (rownrm*pnorm)
                        alfa1  = res / absatp
                        jadd1  = j
                        atp1   = atp
                        hlow1  = .true.
                     end if
                  end if
               end if

               if (js. eq. -1)  then

*                 The upper bound is violated.  Test for either a bigger
*                 or smaller alfa2,  depending on the value of firstv.

                  res    = atx - bu(j)
                  if (     (firstv  .and.  palfa2*absatp .ge. res
     $                 .or.  lastv  .and.  palfa2*absatp .le. res)
     $                 .or.  j .eq.  jsave2) then
                     if (apmax2*rownrm*pnorm .lt. absatp) then
                        apmax2 = absatp / (rownrm*pnorm)
                        if      (absatp .ge. one          ) then
                           alfa2 = res / absatp
                        else if (res    .lt. bigalf*absatp) then
                           alfa2 = res / absatp
                        else
                           alfa2 = bigalf
                        end if
                        jadd2  = j
                        atp2   = atp
                        hlow2  = .false.
                     end if
                  end if
               end if
            else if (atp .gt. zero  .and.  js .ne.  -1)  then
*              ---------------------------------------------------------
*              a'x  is increasing and the upper bound is not violated.
*              ---------------------------------------------------------
*              Test for smaller alfa1.

               if (bu(j) .lt. bigbnd) then
                  res = bu(j) - atx
                  if (palfa1*atp .ge. res  .or.  j .eq. jsave1) then
                     if (apmax1*rownrm*pnorm .lt. atp) then
                        apmax1 = atp / (rownrm*pnorm)
                        alfa1  = res / atp
                        jadd1  = j
                        atp1   = atp
                        hlow1  = .false.
                     end if
                  end if
               end if

               if (js .eq. -2)  then

*                 The lower bound is violated.  Test for a new ALFA2.

                  res    = bl(j) - atx
                  if (     (firstv  .and.  palfa2*atp .ge. res
     $                 .or.  lastv  .and.  palfa2*atp .le. res)
     $                 .or.  j .eq.  jsave2) then
                     if (apmax2*rownrm*pnorm .lt. atp) then
                        apmax2 = atp / (rownrm*pnorm)
                        if      (atp .ge. one       ) then
                           alfa2 = res / atp
                        else if (res .lt. bigalf*atp) then
                           alfa2 = res / atp
                        else
                           alfa2 = bigalf
                        end if
                        jadd2  = j
                        atp2   = atp
                        hlow2  = .true.
                     end if
                  end if
               end if
            end if

*            if (cmdbg  .and.  icmdbg(3) .gt. 0)
*     $      write (iPrint, 1200) j, js, featol(j), res, atp, jadd1,
*     $                         alfa1, jadd2, alfa2
         end if
  400 continue

*     ==================================================================
*     Determine alfa, the step to be taken.
*     ==================================================================
*     In the infeasible case, check whether to take the step alfa2
*     rather than alfa1...

      step2 = numinf .gt. 0  .and.  jadd2 .gt. 0

*     We do so if alfa2 is less than alfa1 or (if firstv is false)
*     lies in the range  (alfa1, palfa1)  and has a smaller value of
*     atp.

      step2 = step2 .and. (alfa2 .lt. alfa1   .or.   lastv  .and.
     $                     alfa2 .le. palfa1  .and.  apmax2 .ge. apmax1)

      if (step2) then
         alfa   = alfa2
         palfa  = palfa2
         jadd   = jadd2
         atphit = atp2
         hitlow = hlow2
      else
         alfa   = alfa1
         palfa  = palfa1
         jadd   = jadd1
         atphit = atp1
         hitlow = hlow1

*        If alfa1 is negative, the constraint to be added (jadd)
*        remains unchanged, but alfa may be shortened to the step
*        to the nearest perturbed satisfied constraint along  - p.

         negstp = alfa .lt. zero
         if (negstp) then
            call cmalf1( firstv, negstp, bigalf, bigbnd, pnorm,
     $                   jadd1, jadd2, palfa1, palfa2,
     $                   istate, n, nctotl,
     $                   Anorm, Ap, Ax, bl, bu, featol, p, x )

*            if (cmdbg  .and.  icmdbg(1) .gt. 0)
*     $         write (iPrint, 9000) alfa, palfa1

            alfa = - min( abs( alfa ), palfa1 )
         end if
      end if

*     Test for undefined or infinite step.

      if (jadd .eq. 0) then
         alfa   = bigalf
         palfa  = bigalf
      end if

      if (alfa .ge. bigalf) inform = 3
*      if (cmdbg  .and.  icmdbg(1) .gt. 0  .and.  inform .gt. 0)
*     $   write (iPrint, 9010) jadd, alfa
      return

 1000 format(/ ' cmalf  entered'
     $       / '    j  js         featol        res             Ap',
     $         '     jadd1        alfa1     jadd2        alfa2 '/)
 1200 format( i5, i4, 3g15.5, 2(i6, g17.7) )
 9000 format(/ ' //cmalf //  negative step',
     $       / ' //cmalf //           alfa          palfa'
     $       / ' //cmalf //', 2g15.4 )
 9010 format(/ ' //cmalf //  unbounded step.'
     $       / ' //cmalf //  jadd           alfa'
     $       / ' //cmalf //  ', i4, g15.4 )

*     end of cmalf
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine cmchk ( nerror, msglvl, lcrash, userkx,
     $                   liwork, lwork, litotl, lwtotl,
     $                   n, nclin, ncnln,
     $                   istate, kx, named, names,
     $                   bigbnd, bl, bu, x )

      implicit           double precision(a-h,o-z)
      character*8        names(*)
      logical            named, userkx
      integer            istate(n+nclin+ncnln), kx(n)
      double precision   bl(n+nclin+ncnln), bu(n+nclin+ncnln), x(n)

*     ==================================================================
*     cmchk   checks the data input to various optimizers.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original Fortran 66 version written 10-May-1980.
*     Fortran 77 version written  5-October-1984.
*     This version of cmchk dated  21-Mar-93.       
*     ==================================================================
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2

      logical            ok
      intrinsic          abs
      character*5        id(3)
      data                id(1)   ,  id(2)   ,  id(3)
     $                 / 'varbl'  , 'lncon'  , 'nlcon'   /

      nerror = 0

*     ------------------------------------------------------------------
*     Check that there is enough workspace to solve the problem.
*     ------------------------------------------------------------------
      ok     = litotl .le. liwork  .and.  lwtotl .le. lwork
      if (.not. ok)  then
         nerror = nerror + 1
         if (iPrint .gt. 0) then
*            write (iPrint, 1100) liwork, lwork, litotl, lwtotl
*            write (iPrint, 1110)
         end if
      else if (msglvl .gt. 0)  then
         if (iPrint .gt. 0) then
*            write (iPrint, 1100) liwork, lwork, litotl, lwtotl
         end if
      end if

      if (userkx) then
*        ---------------------------------------------------------------
*        Check for a valid KX.
*        ---------------------------------------------------------------
         ifail = 1
         call cmperm( kx, 1, n, ifail )
         if (ifail .ne. 0) then
*            if (iPrint .gt. 0) write (iPrint, 1300)
            nerror = nerror + 1
         end if
      end if

*     ------------------------------------------------------------------
*     Check the bounds on all variables and constraints.
*     ------------------------------------------------------------------
      do 200, j = 1, n+nclin+ncnln
         b1     = bl(j)
         b2     = bu(j)
         ok     = b1 .lt. b2  .or. 
     $            b1 .eq. b2  .and.  abs(b1) .lt. bigbnd

         if (.not. ok)  then
            nerror = nerror + 1
            if (j .gt. n+nclin)  then
               k  = j - n - nclin
               l  = 3
            else if (j .gt. n)  then
               k  = j - n
               l  = 2
            else
               k = j
               l = 1
            end if
            if (iPrint .gt. 0) then
               if (named) then
                  if (b1 .eq. b2) then
*                     write (iPrint, 1210) names(j), b1, bigbnd
                  else
*                     write (iPrint, 1215) names(j), b1, b2
                  end if
               else 
                  if (b1 .eq. b2) then
*                     write (iPrint, 1200) id(l), k, b1, bigbnd
                  else
*                     write (iPrint, 1205) id(l), k, b1, b2
                  end if
               end if
            end if
         end if
  200 continue

*     ------------------------------------------------------------------
*     If warm start, check  istate.
*     ------------------------------------------------------------------
      if (lcrash .eq. 1) then         
         do 420, j = 1, n+nclin+ncnln
            is     = istate(j)
            ok     = is .ge. (- 2)   .and.   is .le. 4
            if (.not. ok)  then
               nerror = nerror + 1
*               if (iPrint .gt. 0) write (iPrint, 1500) j, is
            end if
  420    continue
      end if

      return

 1100 format(/ ' Workspace provided is     iw(', i8,
     $         '),  w(', i8, ').' /
     $         ' To solve problem we need  iw(', i8,
     $         '),  w(', i8, ').')
 1110 format(/ ' XXX  Not enough workspace to solve problem.')
 1200 format(/ ' XXX  The equal bounds on  ', a5, i3,
     $         '  are infinite.   Bounds =', g16.7,
     $         '  bigbnd =', g16.7)
 1205 format(/ ' XXX  The bounds on  ', a5, i3,
     $         '  are inconsistent.   bl =', g16.7, '   bu =', g16.7)
 1210 format(/ ' XXX  The equal bounds on  ', a8,
     $         '  are infinite.   Bounds =', g16.7,
     $         '  bigbnd =', g16.7)
 1215 format(/ ' XXX  The bounds on  ', a8,
     $         '  are inconsistent.   bl =', g16.7, '   bu =', g16.7)
 1300 format(/ ' XXX  kx has not been supplied as a valid',
     $         '  permutation.' )
 1500 format(/ ' XXX  Component', i5, '  of  istate  is out of',
     $         ' range...', i10)

*     end of  cmchk
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine cmperm( kx, m1, m2, ifail )

      integer            ifail, m1, m2
      integer            kx(m2)

*     ==================================================================
*     CMPERM checks that elements M1 to M2 of KX contain a valid
*     permutation of the integers M1 to M2. The contents of KX are
*     unchanged on exit.
*
*     SOL version of NAG Library routine M01ZBF.
*     Written by N.N.Maclaren, University of Cambridge.
*     This version of CMPERM dated 18-June-1986.
*     ==================================================================
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2

      logical            cmdbg
      integer            lcmdbg
      parameter         (lcmdbg = 5)
      common    /cmdebg/ icmdbg(lcmdbg), cmdbg

      integer            i, ierr, j, k
      intrinsic          abs

*     Check the parameters.

      if (m2 .lt. 1  .or.  m1 .lt. 1  .or.  m1 .gt. m2) then
         ierr = 1
*         if (cmdbg  .and.  icmdbg(3) .gt. 0)
*     $      write (iPrint, fmt=1100) m1, m2
      else
         ierr = 0

*        Check that KX is within range.

         do 20 i = m1, m2
            j = kx(i)
            if ((j .lt. m1) .or. (j .gt. m2)) go to 100
            if (i .ne. j) kx(i) = -j
   20    continue

*        Check that no value is repeated.

         do 60 i = m1, m2
            k = - kx(i)
            if (k .ge. 0) then
               j     = i
   40          kx(j) = k
               j     = k
               k     = - kx(j)
               if (k .gt. 0) go to 40
               if (j .ne. i) go to 120
            end if
   60    continue
      end if

*     Return

   80 if (ierr .ne. 0) then
         ifail = ierr
      else
         ifail = 0
      end if
      return
  100 ierr = 2
*      if (iPrint .gt. 0) write (iPrint, fmt=1200) i, j
      go to 140
  120 ierr = 3
*      if (iPrint .gt. 0) write (iPrint, fmt=1300) j

*     Restore KX.

  140 do 160 i = m1, m2
         kx(i) = abs(kx(i))
  160 continue
      go to 80

 1100 format(/ ' //cmperm//  Illegal parameter values,'
     $       / ' //cmperm//    m1    m1'
     $       / ' //cmperm//', 2i6 )
 1200 format(/ ' XXX  kx(',I6,') contains an out-of-range value =', i16)
 1300 format(/ ' XXX  kx contains a duplicate value =',             i16)

*     end of cmperm
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine cmprt ( msglvl, nfree, ldA,
     $                   n, nclin, nctotl, bigbnd,
     $                   named, names,
     $                   nactiv, istate, kactiv, kx,
     $                   A, bl, bu, c, clamda, rlamda, x )

      implicit           double precision(a-h,o-z)
      character*8        names(*)
      logical            named
      integer            istate(nctotl), kactiv(n), kx(n)
      double precision   A(ldA,*), bl(nctotl), bu(nctotl), c(*),
     $                   clamda(nctotl), rlamda(n), x(n)

*     ==================================================================
*     cmprt   creates the expanded Lagrange multiplier vector clamda.
*     If msglvl .eq 1 or msglvl .ge. 10,  cmprt prints  x,  A*x,
*     c(x),  their bounds, the multipliers, and the residuals (distance
*     to the nearer bound).
*     cmprt is called by lscore and npcore just before exiting.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original Fortran 77 version written  October 1984.
*     This version of  cmprt  dated  10-June-1986.
*     ==================================================================
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2

      logical            cmdbg
      integer            lcmdbg
      parameter         (lcmdbg = 5)
      common    /cmdebg/ icmdbg(lcmdbg), cmdbg

      character*2        ls, lstate(7)
      character*5        id(3), id3
      character*8        id4
      external           ddot
      intrinsic          abs

      parameter        ( zero  = 0.0d+0 )
      data               id(1) / 'varbl' /
      data               id(2) / 'lncon' /
      data               id(3) / 'nlcon' /
      data               lstate(1) / '--' /, lstate(2) / '++' /
      data               lstate(3) / 'fr' /, lstate(4) / 'll' /
      data               lstate(5) / 'ul' /, lstate(6) / 'eq' /
      data               lstate(7) / 'tb' /

      nplin  = n     + nclin
      nz     = nfree - nactiv

*     Expand multipliers for bounds, linear and nonlinear constraints
*     into the  clamda  array.

      call dload ( nctotl, zero, clamda, 1 )
      nfixed = n - nfree
      do 150, k = 1, nactiv+nfixed
         if (k .le. nactiv) j = kactiv(k) + n
         if (k .gt. nactiv) j = kx(nz+k)
         clamda(j) = rlamda(k)
  150 continue

      if (iPrint .eq. 0
     $    .or.  (msglvl .lt. 10  .and.  msglvl .ne. 1)) return

*      write (iPrint, 1100)
      id3 = id(1)

      do 500, j = 1, nctotl
         b1     = bl(j)
         b2     = bu(j)
         wlam   = clamda(j)
         is     = istate(j)
         ls     = lstate(is + 3)
         if (j .le. n) then

*           Section 1 -- the variables  x.
*           ------------------------------
            k      = j
            v      = x(j)

         else if (j .le. nplin) then

*           Section 2 -- the linear constraints  A*x.
*           -----------------------------------------
            if (j .eq. n + 1) then
*               write (iPrint, 1200)
               id3 = id(2)
            end if

            k      = j - n
            v      = ddot  ( n, A(k,1), ldA, x, 1 )
         else

*           Section 3 -- the nonlinear constraints  c(x).
*           ---------------------------------------------

            if (j .eq. nplin + 1) then
*               write (iPrint, 1300)
               id3 = id(3)
            end if

            k      = j - nplin
            v      = c(k)
         end if

*        Print a line for the j-th variable or constraint.
*        -------------------------------------------------
         res    = v - b1
         res2   = b2 - v
         if (abs(res) .gt. abs(res2)) res = res2
         ip     = 1
         if (b1 .le. ( - bigbnd )) ip = 2
         if (b2 .ge.     bigbnd  ) ip = ip + 2
         if (named) then

            id4 = names(j)
            if (ip .eq. 1) then
*               write (iPrint, 2100) id4,    ls, v, b1, b2, wlam, res
            else if (ip .eq. 2) then
*               write (iPrint, 2200) id4,    ls, v,     b2, wlam, res
            else if (ip .eq. 3) then
*               write (iPrint, 2300) id4,    ls, v, b1,     wlam, res
            else
*               write (iPrint, 2400) id4,    ls, v,         wlam, res
           end if

         else

            if (ip .eq. 1) then
*               write (iPrint, 3100) id3, k, ls, v, b1, b2, wlam, res
            else if (ip .eq. 2) then
*               write (iPrint, 3200) id3, k, ls, v,     b2, wlam, res
            else if (ip .eq. 3) then
*               write (iPrint, 3300) id3, k, ls, v, b1,     wlam, res
            else
*               write (iPrint, 3400) id3, k, ls, v,         wlam, res
           end if
         end if
  500 continue
      return

 1100 format(// ' Variable        State', 5x, ' Value',
     $   6x, ' Lower bound', 4x, ' Upper bound',
     $   '  Lagr multiplier', '     Residual' /)
 1200 format(// ' Linear constr   State', 5x, ' Value',
     $   6x, ' Lower bound', 4x, ' Upper bound',
     $   '  Lagr multiplier', '     Residual' /)
 1300 format(// ' Nonlnr constr   State', 5x, ' Value',
     $   6x, ' Lower bound', 4x, ' Upper bound',
     $   '  Lagr multiplier', '     Residual' /)
 2100 format(1x, a8, 10x, a2, 3g16.7, g16.7, g16.4)
 2200 format(1x, a8, 10x, a2, g16.7, 5x, ' None', 6x, g16.7,
     $   g16.7, g16.4)
 2300 format(1x, a8, 10x, a2, 2g16.7, 5x, ' None', 6x, g16.7, g16.4)
 2400 format(1x, a8, 10x, a2,  g16.7, 5x, ' None', 11x, ' None',
     $   6x, g16.7, g16.4)
 3100 format(1x, a5, i3, 10x, a2, 3g16.7, g16.7, g16.4)
 3200 format(1x, a5, i3, 10x, a2,  g16.7,
     $   5x, ' None', 6x, g16.7, g16.7, g16.4)
 3300 format(1x, a5, i3, 10x, a2, 2g16.7, 5x, ' None', 6x,
     $   g16.7, g16.4)
 3400 format(1x, a5, i3, 10x, a2,  g16.7,
     $   5x, ' None', 11x, ' None', 6x, g16.7, g16.4)

*     end of  cmprt
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine cmqmul( mode, n, nz, nfree, ldQ, unitQ,
     $                   kx, v, Q, w )

      implicit           double precision(a-h,o-z)
      logical            unitQ
      integer            kx(n)
      double precision   v(n), Q(ldQ,*), w(n)

*     ==================================================================
*     cmqmul  transforms the vector  v  in various ways using the
*     matrix  Q = ( Z  Y )  defined by the input parameters.
*
*        MODE               result
*        ----               ------
*
*          1                v = Z v
*          2                v = Y v
*          3                v = Q v
*
*     On input,  v  is assumed to be ordered as  ( v(free)  v(fixed) ).
*     On output, v  is a full n-vector.
*
*
*          4                v = Z'v
*          5                v = Y'v
*          6                v = Q'v
*
*     On input,  v  is a full n-vector.
*     On output, v  is ordered as  ( v(free)  v(fixed) ).
*
*          7                v = Y'v
*          8                v = Q'v
*
*     On input,  v  is a full n-vector.
*     On output, v  is as in modes 5 and 6 except that v(fixed) is not
*     set.
*
*     Modes  1, 4, 7 and 8  do not involve  v(fixed).
*     Original F66 version  April 1983.
*     Fortran 77 version written  9-February-1985.
*     Level 2 BLAS added 10-June-1986.
*     This version of cmqmul dated 14-Sep-92.
*     ==================================================================
      parameter        ( zero = 0.0d+0, one = 1.0d+0 )

      nfixed = n - nfree
      j1     = 1
      j2     = nfree
      if (mode .eq. 1  .or.  mode .eq. 4) j2 = nz
      if (mode .eq. 2  .or.  mode .eq. 5  .or.  mode .eq. 7) j1 = nz + 1
      lenv   = j2 - j1 + 1
      if (mode .le. 3) then
*        ===============================================================
*        Mode = 1, 2  or  3.
*        ===============================================================

         if (nfree .gt. 0) call dload ( nfree, zero, w, 1 )

*        Copy  v(fixed)  into the end of  wrk.

         if (mode .ge. 2  .and.  nfixed .gt. 0)
     $      call dcopy ( nfixed, v(nfree+1), 1, w(nfree+1), 1 )

*        Set  W  =  relevant part of  Q * V.

         if (lenv .gt. 0)  then
            if (unitQ) then
               call dcopy ( lenv, v(j1), 1, w(j1), 1 )
            else
               call dgemv ( 'n', nfree, j2-j1+1, one, Q(1,j1), ldQ,
     $                      v(j1), 1, one, w, 1 )
            end if
         end if

*        Expand  w  into  v  as a full n-vector.

         call dload ( n, zero, v, 1 )
         do 220, k = 1, nfree
            j      = kx(k)
            v(j)   = w(k)
  220    continue

*        Copy  w(fixed)  into the appropriate parts of  v.

         if (mode .gt. 1)  then
            do 320, l = 1, nfixed
               j       = kx(nfree+l)
               v(j)    = w(nfree+l)
  320       continue
         end if

      else
*        ===============================================================
*        Mode = 4, 5, 6, 7  or  8.
*        ===============================================================
*        Put the fixed components of  v  into the end of  w.

         if (mode .eq. 5  .or.  mode .eq. 6)  then
            do 420, l = 1, nfixed
               j          = kx(nfree+l)
               w(nfree+l) = v(j)
  420       continue
         end if

*        Put the free  components of  v  into the beginning of  w.

         if (nfree .gt. 0)  then
            do 520, k = 1, nfree
               j      = kx(k)     
               w(k) = v(j)
  520       continue

*           Set  v  =  relevant part of  Q' * w.

            if (lenv .gt. 0)  then
               if (unitQ) then
                  call dcopy ( lenv, w(j1), 1, v(j1), 1 )
               else
                  call dgemv ( 'T', nfree, j2-j1+1, one, Q(1,j1), ldQ,
     $                         w, 1, zero, v(j1), 1 )
               end if
            end if
         end if

*        Copy the fixed components of  w  into the end of  v.

         if (nfixed .gt. 0  .and.  (mode .eq. 5  .or.  mode .eq. 6))
     $      call dcopy ( nfixed, w(nfree+1), 1, v(nfree+1), 1 )
      end if

*     end of  cmqmul
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine cmr1md( n, nu, nrank, ldR, lenv, lenw,
     $                   R, u, v, w, c, s )

      implicit           double precision(a-h,o-z)
      integer            n, nu, nrank, ldR, lenv, lenw
      double precision   R(ldR,*), u(n,*), v(n), w(n)
      double precision   c(n), s(n)

*     ==================================================================
*     cmr1md  modifies the  nrank*n  upper-triangular matrix  R  so that
*     Q*(R + v*w')  is upper triangular,  where  Q  is orthogonal,
*     v  and  w  are vectors, and the modified  R  overwrites the old.
*     Q  is the product of two sweeps of plane rotations (not stored).
*     If required,  the rotations are applied to the nu columns of
*     the matrix  U.
*
*     The matrix v*w' is an (lenv) by (lenw) matrix.
*     The vector v is overwritten.
*                                           
*     Systems Optimization Laboratory, Stanford University.
*     Original version   October  1984.
*     Level-2 matrix routines added 22-Apr-1988.
*     This version of  cmr1md  dated 22-Apr-1988.
*     ==================================================================
      intrinsic          min

      j = min( lenv, nrank )
      if (nrank .gt. 0) then
*        ---------------------------------------------------------------
*        Reduce  v to beta*e( j )  using a backward sweep of rotations
*        in planes (j-1, j), (j-2, j), ..., (1, j).
*        ---------------------------------------------------------------
         call f06fqf( 'Fixed', 'Backwards', j-1, v(j), v, 1, c, s )

*        ---------------------------------------------------------------
*        Apply the sequence of rotations to U.
*        ---------------------------------------------------------------
         if (nu .gt. 0)
     $      call f06qxf( 'Left', 'Bottom', 'Backwards', j, nu,
     $                   1, j, c, s, u, n )

*        ---------------------------------------------------------------
*        Apply the sequence of rotations to R. This generates a spike in
*        the j-th row of R, which is stored in s.
*        ---------------------------------------------------------------
         call f06qwf( 'Left', n, 1, j, c, s, R, ldR )

*        ---------------------------------------------------------------
*        Form  beta*e(j)*w' + R.  This a spiked matrix, with a row
*        spike in row j.
*        ---------------------------------------------------------------
         call daxpy( min( j-1, lenw ), v(j), w   , 1, s     , 1     )
         call daxpy( lenw-j+1        , v(j), w(j), 1, R(j,j), ldR )

*        ---------------------------------------------------------------
*        Eliminate the spike using a forward sweep of rotations in
*        planes (1, j), (2, j), ..., (j-1, j).
*        ---------------------------------------------------------------
         call f06qsf( 'Left', n, 1, j, c, s, R, ldR )

*        ---------------------------------------------------------------
*        Apply the rotations to U.
*        ---------------------------------------------------------------
         if (nu .gt. 0)
     $      call f06qxf( 'Left', 'Bottom', 'Forwards', j, nu,
     $                   1, j, c, s, u, n )
      end if

*     end of  cmr1md
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine cmrswp( n, nU, nrank, ldR, i, j, R, U, c, s )

      implicit           double precision(a-h,o-z)
      double precision   R(ldR,*), U(n,*)
      double precision   c(n), s(n)

*     ==================================================================
*     CMRSWP  interchanges the  i-th  and  j-th  (i .lt. j)  columns of
*     an  nrank x n  upper-trapezoidal matrix  R   and restores the
*     resulting matrix to upper-trapezoidal form using two sweeps of 
*     plane rotations applied on the left.  R is overwritten.  
*
*     If nU .gt. 0,  the rotations are applied to the  nU  columns of
*     the matrix  U.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 31-October-1984.
*     Level-2 matrix routines added 13-May-1988.
*     This version of  CMRSWP  dated  26-Aug-1991.
*     ==================================================================

      parameter        ( zero = 0.0d+0 )
      intrinsic          min

*     Swap the elements of the i-th and j-th columns of R on, or above,
*     the main diagonal.

      call dswap( min(i,nrank), R(1,i), 1, R(1,j), 1 )
      lenRj = min(j, nrank)

      if (lenRj .gt. i) then
*        ---------------------------------------------------------------
*        Reduce elements  R(i+1,j), ..., R(lenRj,j)  to  beta*e(lenRj)  
*        using a backward sweep in planes
*        (lenRj-1,lenRj), (lenRj-2,lenRj), ..., (i+1,lenRj).
*        If required, apply the sequence of rotations to U.
*        ---------------------------------------------------------------
         call f06fqf( 'Fixed', 'Backwards', lenRj-i-1, R(lenRj,j),
     $                R(i+1,j), 1, c(i+1), s(i+1) )

         if (nU .gt. 0)
     $      call f06qxf( 'Left', 'Bottom', 'Backwards', n, nU, 
     $                   i+1, lenRj, c, s, U, n )

*        Put zeros into the j-th column of R in positions corresponding 
*        to the sub-diagonals of the i-th column.

         s(i) = R(lenRj,j)
         call dload ( lenRj-i, zero, R(i+1,j), 1 )

*        Apply the sequence of rotations to R.  This generates a spike 
*        in the (lenRj)-th row of R, which is stored in s.

         call f06qwf( 'Left', n, i+1, lenRj, c, s, R, ldR )

*        Eliminate the spike using a forward sweep in planes
*        (i,lenRj), (i+1,lenRj), ..., (lenRj-1,lenRj).
*        If necessary, apply the sequence of rotations to U.

         call f06qsf( 'Left', n, i, lenRj, c, s, R, ldR )

         if (nU .gt. 0)
     $      call f06qxf( 'Left', 'Bottom', 'Forwards', lenRj, nU,
     $                   i, lenRj, c, s, U, n )
      end if

*     end of cmrswp
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine cmtsol( mode, ldt, n, t, y )

      implicit           double precision(a-h,o-z)
      integer            mode, ldt, n
      double precision   t(ldt,*), y(n)

*     ==================================================================
*     cmtsol  solves equations involving a reverse-triangular matrix  T
*     and a right-hand-side vector  y,  returning the solution in  y.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original Fortran 77 version written February-1985.
*     ==================================================================
      parameter        ( zero = 0.0d+0 )

      n1 = n + 1
      if (mode .eq. 1) then

*        Mode = 1  ---  Solve  T * y(new) = y(old).

         do 100, j = 1, n
            jj     = n1 - j
            yj     = y(j)/t(j,jj)
            y(j)   = yj
            l      = jj - 1
            if (l .gt. 0  .and.  yj .ne. zero)
     $         call daxpy( l, (-yj), t(j+1,jj), 1, y(j+1), 1 )
  100    continue
      else

*        Mode = 2  ---  Solve  T' y(new) = y(old).

         do 500, j = 1, n
            jj     = n1 - j
            yj     = y(j)/t(jj,j)
            y(j)   = yj
            l      = jj - 1
            if (l .gt. 0  .and.  yj .ne. zero)
     $         call daxpy( l, (-yj), t(jj,j+1), ldt, y(j+1), 1 )
  500    continue
      end if

*     Reverse the solution vector.

      if (n .gt. 1) then
         l = n/2
         do 800, j = 1, l
            jj     = n1 - j
            yj     = y(j)
            y(j)   = y(jj)
            y(jj)  = yj
  800    continue
      end if

*     end of cmtsol
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     File  mcsubs.f
*
*     mchpar   mcenvn   mcenv2   mcstor   mcmin    mcclos   mcopen
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine mchpar()

*     ==================================================================
*     MCHPAR  must define certain machine parameters as follows:
*
*     wmach(1)  = NBASE  = base of floating-point arithmetic.
*     wmach(2)  = NDIGIT = no. of base wmach(1) digits of precision.
*     wmach(3)  = EPS    = floating-point precision.
*     wmach(4)  = RTEPS  = sqrt(EPS).
*     wmach(5)  = RMIN   = smallest positive normalized floating-point
*                          number.
*     wmach(6)  = RTRMIN = sqrt(RMIN).
*     wmach(7)  = BIG    = a very large number that can be represented
*                          without overflow. BIG is equal to 1.0/SMALL,
*                          where SMALL is a small number for which the
*                          machine can evaluate  1.0/SMALL  without
*                          causing overflow.
*     wmach(8)  = RTBIG  = sqrt(BIG).
*     wmach(9)  = UNDFLW = 0 if underflow is not fatal, +ve otherwise.
*                          (not implemented in post-1982 versions).
*     wmach(10) = NIN    = standard file number of the input stream.
*     wmach(11) = NOUT   = standard file number of the output stream.
*     ==================================================================
      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/

      logical            first , hdwire
      integer            emin  , nbase , ndigit, nin   , nout

      double precision   base  , eps   , big   , rmin
      double precision   small , undflw
      intrinsic          sqrt
      save               first
      data               first / .true. /


      if (first) then
         first = .false.

*        ---------------------------------------------------------------
*        Machine-dependent code.
*        1. Set UNDFLW, NIN, NOUT, HDWIRE as desired.
*        2. If  HDWIRE = .TRUE.  set the machine constants
*               NBASE, NDIGIT, EPS, RMIN, BIG
*           in-line.  Otherwise, they will be computed by MCENVN.
*        ---------------------------------------------------------------
         undflw = 0
         nin    = 5
C-->     nout   = 6
C-->     nout   = 9
         nout   = 9
C-->     call mcopen ( nin , 'IN ' )
C-->     call mcopen ( nout, 'OUT' )

         hdwire = .false.

         if (hdwire) then

*           IEEE standard floating-point arithmetic.
*           (Rounded arithmetic is mandated).

            nbase  = 2
            ndigit = 52
            base   = nbase
            eps    = base**(- ndigit)
            rmin   = base**(- 126)
            big    = base**(+ 127)
         else
            call mcenvn( nbase, ndigit, eps, emin, rmin )
            small  = rmin*nbase**4
            big    = 1.0d0/small
         end if

         wmach( 1) = nbase
         wmach( 2) = ndigit
         wmach( 3) = eps
         wmach( 4) = sqrt( eps )
         wmach( 5) = rmin
         wmach( 6) = sqrt( rmin )
         wmach( 7) = big
         wmach( 8) = sqrt( big )
         wmach( 9) = undflw
         wmach(10) = nin
         wmach(11) = nout
      end if

      return

*     end of  mchpar.

      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine mcenvn( beta, t, eps, emin, rmin )
C     ------------------------------------------------------------------
C     Based on NAG Mark 1.0 release ENVIRN.
C
C     MCENVN returns the machine parameters given by:
C
C        BETA - INTEGER.
C               The base of the machine.
C
C        T    - INTEGER.
C               The number of ( BETA ) digits in the mantissa.
C
C        EPS  - REAL.
C               The smallest positive number such that
C
C                  fl( 1.0 - EPS ) .lt. 1.0,
C
C               where fl denotes the computed value.
C
C        EMIN - INTEGER.
C               The minimum exponent before (gradual) underflow occurs.
C
C        RMIN - REAL.
C               The smallest normalized number for the machine given by
C               BASE**( EMIN - 1 ), where BASE is the floating point
C               value of BETA.
C
C
C     The computation of EPS, EMIN and RMIN is based on a routine,
C     PARANOIA by W. Kahan of the University of California at Berkeley.
C
C
C     Nag Fortran 77 O( 1 ) basic linear algebra routine (ENVIRN).
C
C     -- Written on 2-June-1987.
C     Sven Hammarling, Mick Pont and Janet Welding, Nag Central Office.
C     Modified by PEG, 7-Aug-1990.
C     ------------------------------------------------------------------
      double precision   eps, rmin
      integer            beta, emin, t

      double precision   a, b, leps, lrmin, one, rbase, small, two, zero
      integer            gnmin, gpmin, i, lbeta, lemin, lt, ngnmin,
     *                   ngpmin, nout
      logical            first, iwarn, lrnd

      double precision   mcstor
      external           mcstor

      external           mcenv2, mcmin

      intrinsic          abs, max, min

      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2

      save               first, iwarn, lbeta, lemin, leps, lrmin, lt

      data               first/.true./, iwarn/.false./

      if (first) then
         first = .false.
         zero  = 0
         one   = 1
         two   = 2

*        LBETA, LT, LEPS, LEMIN and LRMIN are the local values of BETA,
*        T, EPS, EMIN and RMIN.
*
*        Throughout this routine we use the function MCSTOR to ensure
*        that relevant values are stored and not held in registers, or
*        are not affected by optimizers.
*
*        MCENV2 returns the parameters LBETA and LT. ( LRND is not used
*        here. )

         call mcenv2( lbeta, lt, lrnd )

*        Start to find EPS.

         b = lbeta
         if (lrnd) then
            leps = (b**(1-lt))/two
         else
            leps = b**(1-lt)
         end if

*        Computation of EPS complete. Now find EMIN.
*        Let a = + or - 1, and + or - (1 + base**(-3)).
*        Keep dividing a by BETA until (gradual) underflow occurs.
*        This is detected when we cannot recover the previous a.

         rbase = one/lbeta
         small = one
         do 20, i = 1, 3
            small = mcstor( small*rbase, zero )
   20    continue
         a     = mcstor( one, small )
         call mcmin ( ngpmin,  one, lbeta )
         call mcmin ( ngnmin, -one, lbeta )
         call mcmin ( gpmin ,    a, lbeta )
         call mcmin ( gnmin ,   -a, lbeta )

         if (ngpmin .eq. ngnmin  .and.  gpmin .eq. gnmin) then

            if (ngpmin .eq. gpmin) then
               lemin = ngpmin
*              ( Non twos-complement machines, no gradual underflow;
*              eg VAX )
            else if (gpmin-ngpmin .eq. 3) then
               lemin = ngpmin - 1 + lt
*              ( Non twos-complement machines, with gradual underflow;
*              eg IEEE standard followers )
            else
               lemin = min( ngpmin, gpmin )
*              ( A guess; no known machine )
               iwarn = .true.
            end if

         else if (ngpmin .eq. gpmin .and. ngnmin .eq. gnmin) then
            if (abs( ngpmin-ngnmin ) .eq. 1) then
               lemin = max( ngpmin, ngnmin )
*              ( Twos-complement machines, no gradual underflow;
*              eg Cyber 205 )
            else
               lemin = min( ngpmin, ngnmin )
*              ( A guess; no known machine )
               iwarn = .true.
            end if

         else if (abs(ngpmin-ngnmin) .eq. 1 .and. gpmin .eq. gnmin) then
            if (gpmin-min( ngpmin, ngnmin ) .eq. 3) then
               lemin = max( ngpmin, ngnmin ) - 1 + lt
*              ( Twos-complement machines with gradual underflow;
*              no known machine )
            else
               lemin = max( ngpmin, ngnmin )
*              ( A guess; no known machine )
               iwarn = .true.
            end if
         else
            lemin = min( ngpmin, ngnmin, gpmin, gnmin )
*           ( A guess; no known machine )
            iwarn = .true.
         end if
*        **
*        Comment out this IF block if Emin is ok
         if (iwarn) then
            first = .true.
*            write (nout,fmt=99999) lemin
         end if
*        **

*        Finally compute RMIN by successive division by BETA.
*        We could compute RMIN as base**( EMIN - 1 ), but some machines
*        underflow during this computation.

         lrmin = 1
         do 40, i = 1, 1 - lemin
            lrmin = lrmin/lbeta
   40    continue
      end if

      beta  = lbeta
      t     = lt
      eps   = leps
      emin  = lemin
      rmin  = lrmin
      return

*     End of mcenvn (envirn).

99999 format( // ' WARNING. The value Emin may be incorrect:-  Emin = ',
     $           I8 / ' If, after inspection, the value Emin looks',
     $           ' acceptable please comment out ' / ' the IF block',
     $           ' as marked within the code of routine MCENVN,' /
     $           ' otherwise contact UCSD. ' / )
      end
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine mcenv2( beta, t, rnd )
C     ------------------------------------------------------------------
C     Based on NAG Mark 1.0 release ENVRON.
C
C     MCENV2 returns the machine parameters given by:
C
C        BETA - INTEGER.
C               The base of the machine.
C
C        T    - INTEGER.
C               The number of ( BETA ) digits in the mantissa.
C
C        RND  - LOGICAL.
C               Whether proper rounding ( RND = .TRUE. ) or chopping
C               ( RND = .FALSE. ) occurs in addition. This may not be a
C               reliable guide to the way in which the machine perfoms
C               its arithmetic.
C
C     The routine is based on the routine of the same name by Malcolm
C     and incorporates suggestions by Gentleman and Marovich. See
C
C        Malcolm M. A. (1972) Algorithms to reveal properties of
C           floating-point arithmetic. Comms. of the ACM, 15, 949-951.
C
C        Gentleman W. M. and Marovich S. B. (1974) More on algorithms
C           that reveal properties of floating point arithmetic units.
C           Comms. of the ACM, 17, 276-277.
C
C
C     Nag Fortran 77 O( 1 ) basic linear algebra routine (envron).
C
C     -- Written on 26-November-1984.
C     Sven Hammarling and Mick Pont, Nag Central Office.
C     Modified by PEG, 7-Aug-1990.
C     ------------------------------------------------------------------
      integer            beta, t
      logical            rnd

      double precision   a, b, c, c1, c2, f, one, qtr, sava
      integer            lbeta, lt
      logical            first, lrnd

      double precision   mcstor
      external           mcstor

      save               first, lbeta, lrnd, lt

      data               first/.true./


      if (first) then
         first = .false.
         one   = 1

*        LBETA, LT and LRND are the local values of BETA, T and RND.
*
*        Throughout this routine we use the function MCSTOR to ensure
*        that relevant values are stored and not held in registers, or
*        are not affected by optimizers.
*
*        Compute  a = 2.0**m  with the smallest positive integer m such
*        that
*
*           fl( a + 1.0 ) = a.

         a  = 1
         c  = 1

*       +       while( c .eq. one )loop
   20    if (c .eq. one) then
            a  = 2*a
            c  = mcstor( a, one )
            c  = mcstor( c, -a  )
            go to 20
         end if
*       +       end while

*        Now compute  b = 2.0**m  with the smallest positive integer m
*        such that
*
*           fl( a + b ) .gt. a.

         b  = 1
         c  = mcstor( a, b )

*       +       while( c .eq. a )loop
   40    if (c .eq. a) then
            b  = 2*b
            c  = mcstor( a, b )
            go to 40
         end if
*       +       end while

*        Now compute the base. a and b are neighbouring floating point
*        numbers in the interval ( beta**t, beta**( t + 1 ) ) and so
*        their difference is beta. Adding 0.25 to c is to ensure that it
*        is truncated to beta and not ( beta - 1 ).

         qtr   = one/4
         c     = mcstor( c, -a )
         lbeta = c + qtr

*        Now determine whether rounding or chopping occurs, by adding
*        a bit less than beta/2 and a bit more than beta/2 to a.

         b    = lbeta
         f    = mcstor( b/2, -b/100 )
         c1   = mcstor( f  ,  a     )
         f    = mcstor( b/2,  b/100 )
         c2   = mcstor( f  ,  a     )
         sava = a

*        Now find the mantissa, t. It should be the integer part of
*        log to the base beta of a, however it is safer to determine t
*        by powering. So we find t as the smallest positive integer
*        for which
*
*           fl( beta**t + 1.0 ) = 1.0.

         lt = 0
         a  = 1
         c  = 1

*        +       while( c .eq. one )loop
   60    if (c .eq. one) then
            lt  = lt + 1
            a   = a*lbeta
            c   = mcstor( a, one )
            c   = mcstor( c,  -a )
            go to 60
         end if
*        +       end while

         lrnd = c1 .eq. sava  .and.  c2 .ne. sava
      end if

      beta = lbeta
      t    = lt
      rnd  = lrnd

      return

*     End of mcenv2 (envron).

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      double precision function mcstor( a, b )
C     ------------------------------------------------------------------
C     Based on NAG Mark 1.0 release.
C
C     MCSTOR is intended to force A and B to be stored prior to doing the
C     addition of A and B. For use in situations where optimizers might
C     hold one of these in a register.
C
C
C     Nag Fortran 77 O( 1 ) basic linear algebra routine (mcstor).
C
C     -- Written on 28-November-1984.
C     Sven Hammarling, Nag Central Office.
C     ------------------------------------------------------------------
      double precision                a, b

      mcstor = a + b

      return

C     end of mcstor.

      end
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine mcmin ( emin, start, base )
C     ------------------------------------------------------------------
C     Based on NAG Mark 1.0 release.
C
C     Service routine for ENVIRN (mcenv2).
C
C
C     Nag Fortran 77 O( 1 ) basic linear algebra routine (getmin).
C
C     -- Written on 2-June-1987.
C     Mick Pont, Nag Central Office.
C     ------------------------------------------------------------------
      double precision   start
      integer            base, emin

      double precision   a, b1, b2, c1, c2, d1, d2, one, rbase, zero
      integer            i

      double precision   mcstor
      external           mcstor

      a     = start
      one   = 1
      rbase = one/base
      zero  = 0
      emin  = 1
      b1    = mcstor( a*rbase, zero )
      c1    = a
      c2    = a
      d1    = a
      d2    = a


   20 if ((c1 .eq. a)  .and.  (c2 .eq. a)  .and.
     $    (d1 .eq. a)  .and.  (d2 .eq. a)      ) then
         emin = emin - 1
         a    = b1
         b1   = mcstor(  a/base, zero )
         c1   = mcstor( b1*base, zero )
         d1   = zero
         do 40, i = 1, base
            d1 = d1 + b1
   40    continue
         b2   = mcstor(  a*rbase, zero )
         c2   = mcstor( b2/rbase, zero )
         d2   = zero
         do 60, i = 1, base
            d2 = d2 + b2
   60    continue
         go to 20
      end if

      return
      
*     End of mcmin (getmin).

      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine mcclos( lun )

*     ==================================================================
*     mcclos  closes the file with logical unit number lun.
*     ==================================================================

*      close ( lun )

*     end of mcclos
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine mcopen( lun, state )

      integer            lun
      character*3        state

*     ==================================================================
*     mcopen  opens the file with logical unit number lun.
*
*     state   is intended to be input as 'IN ' or 'OUT'.  It may
*     be helpful on some systems.  Only  sequential files are used and
*     it is never necessary to read and write to the same file.
*
*     'IN ' refers to an existing input file that will not be altered.
*     'OUT' means that a new file will be output.  If the file already
*     exists, it might be OK to overwrite it, but on some systems it
*     is better to create a new version of the file.  The choice is
*     open (to coin a phrase).
*
*     15-Nov-91: First version based on Minos 5.4 routine m1open.
*     20-Oct-92: Current version.
*     ==================================================================

*      if ( state .eq. 'IN ' ) then

*        Open an input file (e.g., OPTIONS, LOAD).
*        'OLD' means there will be an error if the file does not exist.
*        Since some systems position existing files at the end
*        (rather than the beginning), a rewind is performed.
      
*         open  ( lun, status='OLD' )
*         rewind( lun, err=900 )

*      else if ( state .eq. 'OUT' ) then

*        Open an output file (e.g. SUMMARY, SOLUTION).
*        If it is OK to overwrite an existing file, we could do the
*        same as for input:

*---     open  ( lun )
*---     rewind( lun, err=900 )

*        On DEC VAX/VMS systems it is better to let a new generation
*        be created when the first write occurs, so we do nothing:

*         open( lun, status='unknown' )
*      end if

*  900 return

*     end of mcopen.
      end
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     File  qrsubs.f
*
*     dgeqr    dgeqrp   dgeap    dgeapq
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE DGEQR ( M, N, A, LDA, ZETA, INFORM )
      INTEGER            M, N, LDA, INFORM
      DOUBLE PRECISION   A( LDA, * ), ZETA( * )
C
C  1. Purpose
C     =======
C
C  DGEQR  reduces the  m by n, m.ge.n, matrix A to upper triangular form
C  by means of orthogonal transformations.
C
C  2. Description
C     ===========
C
C  The m by n matrix A is factorized as
C
C     A = Q*( R )   when   m.gt.n,
C           ( 0 )
C
C     A = Q*R       when   m = n,
C
C  where  Q  is an  m by m  orthogonal matrix and  R  is an n by n upper
C  triangular matrix.
C
C  The  factorization  is  obtained  by  Householder's  method. The  kth
C  transformation matrix, Q( k ), which is used to introduce zeros  into
C  the kth column of A is given in the form
C
C     Q( k ) = ( I     0   ),
C              ( 0  T( k ) )
C
C  where
C
C     T( k ) = I - u( k )*u( k )',   u( k ) = ( zeta( k ) ),
C                                             (    z( k ) )
C
C  zeta( k )  is a scalar and  z( k )  is an  ( m - k )  element vector.
C  zeta( k )  and  z( k ) are chosen to annhilate the elements below the
C  triangular part of  A.
C
C  The vector  u( k )  is returned in the kth element of ZETA and in the
C  kth column of A, such that zeta( k ) is in ZETA( k ) and the elements
C  of z( k ) are in a( k + 1, k ), ..., a( m, k ). The elements of R are
C  returned in the upper triangular part of  A.
C
C  Q is given by
C
C     Q = ( Q( p )*Q( p - 1 )*...*Q( 1 ) )',
C
C  where p = min( n, m - 1 ).
C
C  3. Parameters
C     ==========
C
C  M      - INTEGER.
C
C           On entry, M must specify the number of rows of  A. M must be
C           at least  n.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C
C           On entry, N must specify the number of columns of  A. N must
C           be  at  least zero. When  N = 0  then an immediate return is
C           effected.
C
C           Unchanged on exit.
C
C  A      - 'real' array of DIMENSION ( LDA, n ).
C
C           Before entry, the leading  M by N  part of the array  A must
C           contain the matrix to be factorized.
C
C           On exit, the  N by N upper triangular part of A will contain
C           the  upper  triangular  matrix  R  and the  M by N  strictly
C           lower triangular part of  A  will  contain  details  of  the
C           factorization as described above.
C
C  LDA    - INTEGER.
C
C           On entry, LDA  must  specify  the  leading dimension of  the
C           array  A  as declared in the calling (sub) program. LDA must
C           be at least  m.
C
C           Unchanged on exit.
C
C  ZETA   - 'real' array of DIMENSION at least ( n ).
C
C           On  exit, ZETA( k )  contains the scalar  zeta( k )  for the
C           kth  transformation.  If  T( k ) = I  then   ZETA( k ) = 0.0
C           otherwise  ZETA( k )  contains  zeta( k ) as described above
C           and is always in the range ( 1.0, sqrt( 2.0 ) ).
C
C  INFORM - INTEGER.
C
C           On successful  exit  INFORM  will be zero, otherwise  INFORM
C           will  be set to unity indicating that an input parameter has
C           been  incorrectly  set. See  the  next section  for  further
C           details.
C
C  4. Diagnostic Information
C     ======================
C
C  INFORM = 1
C
C     One or more of the following conditions holds:
C
C        M   .lt. N
C        N   .lt. 0
C        LDA .lt. M
C
C  5. Further information
C     ===================
C
C  Following the use of this routine the operations
C
C     B := Q'*B   and   B := Q*B,
C
C  where  B  is an  m by k  matrix, can  be  performed  by calls to  the
C  auxiliary  linear  algebra routine  DGEAPQ. The  operation  B := Q'*B
C  can be obtained by the call:
C
C     INFORM = 0
C     CALL DGEAPQ( 'Transpose', 'Separate', M, N, A, LDA, ZETA,
C    $             K, B, LDB, WORK, INFORM )
C
C  and  B := Q*B  can be obtained by the call:
C
C     INFORM = 0
C     CALL DGEAPQ( 'No transpose', 'Separate', M, N, A, LDA, ZETA,
C    $             K, B, LDB, WORK, INFORM )
C
C  In  both  cases  WORK  must be a  k  element array  that  is used  as
C  workspace. If  B  is a one-dimensional array (single column) then the
C  parameter  LDB  can be replaced by  M. See routine DGEAPQ for further
C  details.
C
C  Operations involving the matrix  R  are performed by  the
C  Level 2 BLAS  routines  DTRMV  and DTRSV . Note that no test for near
C  singularity of R is incorporated in this routine or in routine  DTRSV
C  and  so it is  strongly recommended that the auxiliary linear algebra
C  routine  DUTCO  be called, prior to solving equations involving R, in
C  order  to determine whether  or not  R  is nearly singular. If  R  is
C  nearly  singular  then  the  auxiliary linear algebra  routine  DUTSV
C  can  be used to  determine  the  singular value decomposition  of  R.
C
C
C  Nag Fortran 77 Auxiliary linear algebra routine.
C
C  -- Written on 13-December-1984.
C     Sven Hammarling, Nag Central Office.
C
      EXTERNAL           DGEMV , DGER  , DGRFG
      INTRINSIC          MIN
      INTEGER            K     , LA
      DOUBLE PRECISION   TEMP
      DOUBLE PRECISION   ONE   ,         ZERO
      PARAMETER        ( ONE   = 1.0D+0, ZERO  = 0.0D+0 )

*     Check the input parameters.

      IF( N.EQ.0 )THEN
         INFORM = 0
         RETURN
      END IF
      IF( ( M.LT.N ).OR.( N.LT.0 ).OR.( LDA.LT.M ) )THEN
         INFORM = 1
         RETURN
      END IF

*     Perform the factorization.

      LA = LDA
      DO 20, K = 1, MIN( M - 1, N )

*        Use a Householder reflection to zero the kth column of A.
*        First set up the reflection.

         CALL DGRFG ( M - K, A( K, K ), A( K + 1, K ), 1, ZERO,
     $                ZETA( K ) )
         IF( ( ZETA( K ).GT.ZERO ).AND.( K.LT.N ) )THEN
            IF( ( K + 1 ).EQ.N )
     $         LA = M - K + 1
            TEMP      = A( K, K )
            A( K, K ) = ZETA( K )

*           We now perform the operation  A := Q( k )*A.

*           Let B denote the bottom ( m - k + 1 ) by ( n - k ) part
*           of A.

*           First form  work = B'*u. ( work is stored in the elements
*           ZETA( k + 1 ), ..., ZETA( n ). )

            CALL DGEMV ( 'Transpose', M - K + 1, N - K,
     $                   ONE, A( K, K + 1 ), LA, A( K, K ), 1,
     $                   ZERO, ZETA( K + 1 ), 1 )

*           Now form  B := B - u*work'.

            CALL DGER  ( M - K + 1, N - K, -ONE, A( K, K ), 1,
     $                   ZETA( K + 1 ), 1, A( K, K + 1 ), LA )

*           Restore beta.

            A( K, K ) = TEMP
         END IF
   20 CONTINUE

*     Store the final zeta when m.eq.n.

      IF( M.EQ.N )
     $   ZETA( N ) = ZERO

      INFORM = 0
      RETURN

*     End of DGEQR .

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE DGEQRP( PIVOT, M, N, A, LDA, ZETA, PERM, WORK, INFORM )
      CHARACTER*1        PIVOT
      INTEGER            M, N, LDA, INFORM
      INTEGER            PERM( * )
      DOUBLE PRECISION   A( LDA, * ), ZETA( * ), WORK( * )

      DOUBLE PRECISION   WMACH
      COMMON    /SOLMCH/ WMACH(15)
      SAVE      /SOLMCH/

C  1. Purpose
C     =======
C
C  DGEQRP reduces the  m by n matrix A to upper triangular form by means
C  of orthogonal transformations and column permutations.
C
C  2. Description
C     ===========
C
C  The m by n matrix A is factorized as
C
C     A = Q*( R )*P'      when   m.gt.n,
C           ( 0 )
C
C     A = Q*R*P'          when   m = n,
C
C     A = Q*( R  X )*P'   when   m.lt.n,
C
C  where  Q  is  an  m by m  orthogonal matrix, R  is a  min( m, n )  by
C  min( m, n )  upper triangular matrix and  P is an  n by n permutation
C  matrix.
C
C  The  factorization  is  obtained  by  Householder's  method. The  kth
C  transformation matrix, Q( k ),  which is used to introduce zeros into
C  the kth column of A is given in the form
C
C     Q( k ) = ( I     0   ),
C              ( 0  T( k ) )
C
C  where
C
C     T( k ) = I - u( k )*u( k )',   u( k ) = ( zeta( k ) ),
C                                             (    z( k ) )
C
C  zeta( k )  is a scalar and  z( k )  is an  ( m - k )  element vector.
C  zeta( k )  and  z( k ) are chosen to annhilate the elements below the
C  triangular part of  A.
C
C  The vector  u( k )  is returned in the kth element of ZETA and in the
C  kth column of A, such that zeta( k ) is in ZETA( k ) and the elements
C  of z( k ) are in a( k + 1, k ), ..., a( m, k ). The elements of R are
C  returned in the upper triangular part of A.
C
C  Q is given by
C
C     Q = ( Q( p )*Q( p - 1 )*...*Q( 1 ) )',
C
C  where p = min( m - 1, n ).
C
C  Two options are available for the column permutations. In either case
C  the column for which the  sub-diagonal elements are to be annihilated
C  at the  kth step is chosen from the remaining ( n - k + 1 )  columns.
C  The  particular column chosen as the pivot column is either that  for
C  which  the  unreduced  part  ( elements k onwards )  has the  largest
C  Euclidean  length, or  is that for  which the ratio of the  Euclidean
C  length  of the  unreduced part  to the  Euclidean length of the whole
C  column is a maximum.
C
C  3. Parameters
C     ==========
C
C  PIVOT  - CHARACTER*1.
C
C           On  entry, PIVOT  specifies  the  pivoting  strategy  to  be
C           performed as follows.
C
C           PIVOT = 'C' or 'c'
C
C              Column  interchanges  are  to be  incorporated  into  the
C              factorization, such that the  column whose unreduced part
C              has  maximum  Euclidean  length  is chosen  as the  pivot
C              column at each step.
C
C           PIVOT = 'S' or 's'
C
C              Scaled  column interchanges  are to be  incorporated into
C              the  factorization, such  that the  column for which  the
C              ratio  of the  Euclidean  length of the unreduced part of
C              the column to the original Euclidean length of the column
C              is a maximum is chosen as the  pivot column at each step.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C
C           On entry, M  must specify the number of rows of A. M must be
C           at  least  zero. When  M = 0  then  an  immediate return  is
C           effected.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C
C           On entry, N  must specify the number of columns of A. N must
C           be  at least zero. When  N = 0  then an immediate return  is
C           effected.
C
C           Unchanged on exit.
C
C  A      - 'real' array of DIMENSION ( LDA, n ).
C
C           Before entry, the leading  M by N  part of the array  A must
C           contain the matrix to be factorized.
C
C           On  exit, the  min( M, N ) by min( M, N )  upper  triangular
C           part of A will contain the upper triangular matrix R and the
C           M by min( M, N )  strictly lower triangular part of  A  will
C           contain details  of the  factorization  as  described above.
C           When m.lt.n then the remaining M by ( N - M ) part of A will
C           contain the matrix X.
C
C  LDA    - INTEGER.
C
C           On  entry, LDA  must  specify  the leading dimension of  the
C           array  A  as declared in the calling (sub) program. LDA must
C           be at least  m.
C
C           Unchanged on exit.
C
C  ZETA   - 'real' array of DIMENSION at least ( n ).
C
C           On exit, ZETA( k )  contains the scalar  zeta  for  the  kth
C           transformation. If T( k ) = I then ZETA( k) = 0.0, otherwise
C           ZETA( k )  contains the scalar  zeta( k ) as described above
C           and  is  always  in  the  range  ( 1.0, sqrt( 2.0 ) ).  When
C           n .gt. m  the  elements  ZETA( m + 1 ),  ZETA( m + 2 ), ...,
C           ZETA( n )  are used as internal workspace.
C
C  PERM   - INTEGER array of DIMENSION at least min( m, n ).
C
C           On exit, PERM  contains details of the permutation matrix P,
C           such  that  PERM( k ) = k  if no  column interchange occured
C           at  the  kth  step  and  PERM( k ) = j, ( k .lt. j .le. n ),
C           if  columns  k and j  were  interchanged at  the  kth  step.
C           Note  that, although  there are  min( m - 1, n )  orthogonal
C           transformations, there are min( m, n ) permutations.
C
C  WORK   - 'real' array of DIMENSION at least ( 2*n ).
C
C           Used as internal workspace.
C
C           On exit, WORK( j ), j = 1, 2, ..., n, contains the Euclidean
C           length  of the  jth  column  of the  permuted  matrix  A*P'.
C
C  INFORM - INTEGER.
C
C           On  successful exit, INFORM  will be zero, otherwise  INFORM
C           will  be set to unity indicating that an input parameter has
C           been  incorrectly supplied. See the next section for further
C           details.
C
C  4. Diagnostic Information
C     ======================
C
C  INFORM = 1
C
C     One or more of the following conditions holds:
C
C        PIVOT .ne. 'C' or 'c' or 'S' or 's'
C        M     .lt. 0
C        N     .lt. 0
C        LDA   .lt. M
C
C  5. Further information
C     ===================
C
C  Following the use of this routine the operations
C
C     B := Q'*B   and   B := Q*B,
C
C  where  B  is an  m by k  matrix, can  be  performed  by calls to  the
C  auxiliary  linear algebra  routine  DGEAPQ. The  operation  B := Q'*B
C  can be obtained by the call:
C
C     INFORM = 0
C     CALL DGEAPQ( 'Transpose', 'Separate', M, N, A, LDA, ZETA,
C    $             K, B, LDB, WORK, INFORM )
C
C  and  B := Q*B  can be obtained by the call:
C
C     INFORM = 0
C     CALL DGEAPQ( 'No transpose', 'Separate', M, N, A, LDA, ZETA,
C    $             K, B, LDB, WORK, INFORM )
C
C  In  both  cases  WORK  must be  a  k  element array  that is used  as
C  workspace. If B is a one-dimensional array ( single column ) then the
C  parameter  LDB  can be replaced by  M. See routine DGEAPQ for further
C  details.
C
C  Also following the use of this routine the operations
C
C     B := P'*B   and   B := P*B,
C
C  where B is an n by k matrix, and the operations
C
C     B := B*P    and   B := B*P',
C
C  where  B is a k by n  matrix, can  be performed by calls to the basic
C  linear  algebra  routine  DGEAP .  The  operation  B := P'*B  can  be
C  obtained by the call:
C
C     CALL DGEAP ( 'Left', 'Transpose', N, MIN( M, N ), PERM,
C    $             K, B, LDB )
C
C  the operation  B := P*B  can be obtained by the call:
C
C     CALL DGEAP ( 'Left', 'No transpose', N, MIN( M, N ), PERM,
C    $             K, B, LDB )
C
C  If  B is a one-dimensional array ( single column ) then the parameter
C  LDB  can be replaced by  N  in the above two calls.
C  The operation  B := B*P  can be obtained by the call:
C
C     CALL DGEAP ( 'Right', 'No transpose', K, MIN( M, N ), PERM,
C    $             M, B, LDB )
C
C  and  B := B*P'  can be obtained by the call:
C
C     CALL DGEAP ( 'Right', 'Transpose', K, MIN( M, N ), PERM,
C    $             M, B, LDB )
C
C  If  B is a one-dimensional array ( single column ) then the parameter
C  LDB  can be replaced by  K  in the above two calls.
C  See routine DGEAP for further details.
C
C  Operations involving  the matrix  R  are performed by  the
C  Level 2 BLAS  routines  DTRSV  and DTRMV.  Note that no test for near
C  singularity of  R is incorporated in this routine or in routine DTRSV
C  and  so it is  strongly recommended that the auxiliary linear algebra
C  routine  DUTCO  be called, prior to solving equations involving R, in
C  order  to determine whether  or not  R  is nearly singular. If  R  is
C  nearly  singular then  the  auxiliary  linear algebra  routine  DUTSV
C  can  be  used  to  determine  the  singular value decomposition of R.
C  Operations  involving  the  matrix  X  can also be  performed  by the
C  Level 2  BLAS  routines.  Matrices  of  the  form   ( R  X )  can  be
C  factorized as
C
C     ( R  X ) = ( T  0 )*S',
C
C  where  T is upper triangular and S is orthogonal, using the auxiliary
C  linear algebra routine  DUTRQ .
C
C
C  Nag Fortran 77 Auxiliary linear algebra routine.
C
C  -- Written on 13-December-1984.
C     Sven Hammarling, Nag Central Office.
C
      EXTERNAL           MCHPAR, DGEMV , DGER  , DGRFG , DNRM2 , DSWAP
      INTRINSIC          ABS   , MAX   , MIN   , SQRT
      INTEGER            J     , JMAX  , K     , LA
      DOUBLE PRECISION   EPS   , MAXNRM, NORM  , DNRM2 , TEMP  , TOL
      DOUBLE PRECISION   LAMDA
      PARAMETER        ( LAMDA = 1.0D-2 )
      DOUBLE PRECISION   ONE   ,         ZERO
      PARAMETER        ( ONE   = 1.0D+0, ZERO  = 0.0D+0 )

*     Check the input parameters.

      IF( MIN( M, N ).EQ.0 )THEN
         INFORM = 0
         RETURN
      END IF
      IF( ( ( PIVOT.NE.'C' ).AND.( PIVOT.NE.'c' ).AND.
     $      ( PIVOT.NE.'S' ).AND.( PIVOT.NE.'s' )      ).OR.
     $    ( M.LT.0 ).OR.( N.LT.0 ).OR.( LDA.LT.M )           )THEN
         INFORM = 1
         RETURN
      END IF

*     Compute eps and the initial column norms.

      CALL MCHPAR()
      EPS = WMACH( 3 )
      DO 10, J = 1, N
         WORK( J )     = DNRM2 ( M, A( 1, J ), 1 )
         WORK( J + N ) = WORK( J )
   10 CONTINUE

*     Perform the factorization. TOL is the tolerance for DGRFG .

      LA = LDA
      DO 50, K = 1, MIN( M, N )

*        Find the pivot column.

         MAXNRM = ZERO
         JMAX   = K
         DO 20, J = K, N
            IF( ( PIVOT.EQ.'C' ).OR.( PIVOT.EQ.'c' ) )THEN
               IF( WORK( J + N  ).GT.MAXNRM )THEN
                  MAXNRM = WORK( J + N )
                  JMAX   = J
               END IF
            ELSE IF( WORK( J ).GT.ZERO )THEN
               IF( ( WORK( J + N )/WORK( J ) ).GT.MAXNRM )THEN
                  MAXNRM = WORK( J + N )/WORK( J )
                  JMAX   = J
               END IF
            END IF
   20    CONTINUE
         PERM( K ) = JMAX
         IF( JMAX.GT.K )THEN
            CALL DSWAP ( M, A( 1, K ), 1, A( 1, JMAX ), 1 )
            TEMP             = WORK( K )
            WORK( K )        = WORK( JMAX )
            WORK( JMAX )     = TEMP
            WORK( JMAX + N ) = WORK( K + N )
            PERM( K )        = JMAX
         END IF
         TOL = EPS*WORK( K )
         IF( K.LT.M )THEN

*           Use a Householder reflection to zero the kth column of A.
*           First set up the reflection.

            CALL DGRFG ( M - K, A( K, K ), A( K + 1, K ), 1, TOL,
     $                   ZETA( K ) )
            IF( K.LT.N )THEN
               IF( ZETA( K ).GT.ZERO )THEN
                  IF( ( K + 1 ).EQ.N )
     $               LA = M - K + 1
                  TEMP      = A( K, K )
                  A( K, K ) = ZETA( K )

*                 We now perform the operation  A := Q( k )*A.

*                 Let B denote the bottom ( m - k + 1 ) by ( n - k )
*                 part of A.

*                 First form  work = B'*u. ( work is stored in the
*                 elements ZETA( k + 1 ), ..., ZETA( n ). )

                  CALL DGEMV ( 'Transpose', M - K + 1, N - K,
     $                         ONE, A( K, K + 1 ), LA, A( K, K ), 1,
     $                         ZERO, ZETA( K + 1 ), 1 )

*                 Now form  B := B - u*work'.

                  CALL DGER  ( M - K + 1, N - K, -ONE, A( K, K ), 1,
     $                         ZETA( K + 1 ), 1, A( K, K + 1 ), LA )

*                 Restore beta.

                  A( K, K ) = TEMP
               END IF

*              Update the unreduced column norms. Use the Linpack
*              criterion for when to recompute the norms, except that
*              we retain the original column lengths throughout and use
*              a smaller lamda.

               DO 40, J = K + 1, N
                  IF( WORK( J + N ).GT.ZERO )THEN
                     TEMP = ABS( A( K, J ) )/WORK( J + N )
                     TEMP = MAX( ( ONE + TEMP )*( ONE - TEMP ), ZERO )
                     NORM = TEMP
                     TEMP = ONE +
     $                      LAMDA*TEMP*( WORK( J + N )/WORK( J ) )**2
                     IF( TEMP.GT.ONE )THEN
                        WORK( J + N ) = WORK( J + N )*SQRT( NORM )
                     ELSE
                        WORK( J + N ) = DNRM2 ( M - K,
     $                                          A( K + 1, J ), 1 )
                     END IF
                  END IF
   40          CONTINUE
            END IF
         END IF
   50 CONTINUE

*     Store the final zeta when m.le.n.

      IF( M.LE.N )
     $   ZETA( M ) = ZERO

      INFORM = 0
      RETURN

*     End of DGEQRP.

      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE DGEAP ( SIDE, TRANS, M, N, PERM, K, B, LDB )
*     .. Scalar Arguments ..
      INTEGER            K, LDB, M, N
      CHARACTER*1        SIDE, TRANS
*     .. Array Arguments ..
      DOUBLE PRECISION   B( LDB, * )
      INTEGER            PERM( * )
*     ..
*
*  Purpose
*  =======
*
*  DGEAP  performs one of the transformations
*
*     B := P'*B   or   B := P*B,   where B is an m by k matrix,
*
*  or
*
*     B := B*P'   or   B := B*P,   where B is a k by m matrix,
*
*  P being an m by m permutation matrix of the form
*
*     P = P( 1, index( 1 ) )*P( 2, index( 2 ) )*...*P( n, index( n ) ),
*
*  where  P( i, index( i ) ) is the permutation matrix that interchanges
*  items i and index( i ). That is P( i, index( i ) ) is the unit matrix
*  with rows and columns  i and index( i )  interchanged.  Of course, if
*  index( i ) = i  then  P( i, index( i ) ) = I.
*
*  This routine  is intended for use in  conjunction with  Nag auxiliary
*  routines that  perform  interchange  operations,  such  as  pivoting.
*
*  Parameters
*  ==========
*
*  SIDE   - CHARACTER*1.
*  TRANS
*           On entry,  SIDE  ( Left-hand side, or Right-hand side )  and
*           TRANS  ( Transpose, or No transpose )  specify the operation
*           to be performed as follows.
*
*           SIDE = 'L' or 'l'   and   TRANS = 'T' or 't'
*
*              Perform the operation   B := P'*B.
*
*           SIDE = 'L' or 'l'   and   TRANS = 'N' or 'n'
*
*              Perform the operation   B := P*B.
*
*           SIDE = 'R' or 'r'   and   TRANS = 'T' or 't'
*
*              Perform the operation   B := B*P'.
*
*           SIDE = 'R' or 'r'   and   TRANS = 'N' or 'n'
*
*              Perform the operation   B := B*P.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*
*           On entry, M must specify the order of the permutation matrix
*           P.  M must be at least zero.  When  M = 0  then an immediate
*           return is effected.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*
*           On entry,  N must specify the value of n. N must be at least
*           zero.  When  N = 0  then an  immediate  return is  effected.
*
*           Unchanged on exit.
*
*  PERM   - INTEGER array of DIMENSION at least ( n ).
*
*           Before  entry,  PERM  must  contain the  n  indices  for the
*           permutation matrices. index( i ) must satisfy
*
*              1 .le. index( i ) .le. m.
*
*           It is usual for index( i ) to be at least i, but this is not
*           necessary for this routine.
*
*           Unchanged on exit.
*
*  K      - INTEGER.
*
*           On entry with  SIDE = 'L' or 'l',  K must specify the number
*           of columns of B and on entry with  SIDE = 'R' or 'r', K must
*           specify the  number of rows of B.  K must be at least  zero.
*           When  K = 0  then an immediate return is effected.
*
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of  DIMENSION  ( LDB, ncolb ),  where
*           ncolb = k   when   SIDE = 'L' or 'l'  and   ncolb = m   when
*           SIDE = 'R' or 'r'.
*
*           Before entry  with  SIDE = 'L' or 'l',  the  leading  M by K
*           part  of  the  array   B  must  contain  the  matrix  to  be
*           transformed  and  before entry with  SIDE = 'R' or 'r',  the
*           leading  K by M part of the array  B must contain the matrix
*           to  be  transformed.  On  exit,  B  is  overwritten  by  the
*           transformed matrix.
*
*  LDB    - INTEGER.
*
*           On entry,  LDB  must specify  the  leading dimension  of the
*           array  B  as declared  in the  calling  (sub) program.  When
*           SIDE = 'L' or 'l'   then  LDB  must  be  at  least  m,  when
*           SIDE = 'R' or 'r'   then  LDB  must  be  at  least  k.
*           Unchanged on exit.
*
*
*  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
*
*  -- Written on 13-January-1986.
*     Sven Hammarling, Nag Central Office.
*
*
*     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, J, L
      LOGICAL            LEFT, NULL, RIGHT, TRNSP
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
      IF( MIN( M, N, K ).EQ.0 )
     $   RETURN
      LEFT  = ( SIDE .EQ.'L' ).OR.( SIDE .EQ.'l' )
      RIGHT = ( SIDE .EQ.'R' ).OR.( SIDE .EQ.'r' )
      NULL  = ( TRANS.EQ.'N' ).OR.( TRANS.EQ.'n' )
      TRNSP = ( TRANS.EQ.'T' ).OR.( TRANS.EQ.'t' )
      IF( LEFT )THEN
         IF( TRNSP )THEN
            DO 20, I = 1, N
               IF( PERM( I ).NE.I )THEN
                  L = PERM( I )
                  DO 10, J = 1, K
                     TEMP      = B( I, J )
                     B( I, J ) = B( L, J )
                     B( L, J ) = TEMP
   10             CONTINUE
               END IF
   20       CONTINUE
         ELSE IF( NULL )THEN
            DO 40, I = N, 1, -1
               IF( PERM( I ).NE.I )THEN
                  L = PERM( I )
                  DO 30, J = 1, K
                     TEMP      = B( L, J )
                     B( L, J ) = B( I, J )
                     B( I, J ) = TEMP
   30             CONTINUE
               END IF
   40       CONTINUE
         END IF
      ELSE IF( RIGHT )THEN
         IF( TRNSP )THEN
            DO 60, J = 1, N
               IF( PERM( J ).NE.J )THEN
                  L = PERM( J )
                  DO 50, I = 1, K
                     TEMP      = B( I, J )
                     B( I, J ) = B( L, J )
                     B( L, J ) = TEMP
   50             CONTINUE
               END IF
   60       CONTINUE
         ELSE IF( NULL )THEN
            DO 80, J = N, 1, -1
               IF( PERM( J ).NE.J )THEN
                  L = PERM( J )
                  DO 70, I = 1, K
                     TEMP      = B( L, J )
                     B( L, J ) = B( I, J )
                     B( I, J ) = TEMP
   70             CONTINUE
               END IF
   80       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DGEAP . ( F06QJF )
*
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE DGEAPQ( TRANS, WHEREZ, M, N, A, LDA, ZETA,
     $                   NCOLB, B, LDB, WORK, INFORM )
      CHARACTER*1        TRANS, WHEREZ
      INTEGER            M, N, LDA, NCOLB, LDB, INFORM
      DOUBLE PRECISION   A( LDA, * ), ZETA( * ), B( LDB, * ), WORK( * )
C
C  1. Purpose
C     =======
C
C  DGEAPQ performs one of the transformations
C
C     B := Q'*B   or   B := Q*B,
C
C  where B is an m by ncolb matrix and Q is an m by m orthogonal matrix,
C  given as the product of  Householder transformation matrices, details
C  of  which are stored in the  m by n ( m.ge.n )  array  A  and, if the
C  parameter  WHEREZ = 'S' or 's', in the array ZETA.
C
C  This  routine is  intended for use following auxiliary linear algebra
C  routines such as  DGEQR , DGEHES and DSLTRI. ( See those routines for
C  example calls. )
C
C  2. Description
C     ===========
C
C  Q is assumed to be given by
C
C     Q = ( Q( p )*Q( p - 1 )*...*Q( 1 ) )',
C
C  Q( k ) being given in the form
C
C     Q( k ) = ( I     0   ),
C              ( 0  T( k ) )
C
C  where
C
C     T( k ) = I - u( k )*u( k )',   u( k ) = ( zeta( k ) ),
C                                             (    z( k ) )
C
C  zeta( k )  is a scalar and  z( k )  is an  ( m - k )  element vector.
C
C  z( k )  must  be  supplied  in  the  kth  column  of  A  in  elements
C  a( k + 1, k ), ..., a( m, k )  and  zeta( k ) must be supplied either
C  in  a( k, k )  or in  zeta( k ), depending upon the parameter WHEREZ.
C
C  To obtain Q explicitly B may be set to I and premultiplied by Q. This
C  is more efficient than obtaining Q'.
C
C  3. Parameters
C     ==========
C
C  TRANS  - CHARACTER*1.
C
C           On entry, TRANS  specifies the operation to be performed  as
C           follows.
C
C           TRANS = ' ' or 'N' or 'n'
C
C              Perform the operation  B := Q*B.
C
C           TRANS = 'T' or 't' or 'C' or 'c'
C
C              Perform the operation  B := Q'*B.
C
C           Unchanged on exit.
C
C  WHEREZ - CHARACTER*1.
C
C           On entry, WHEREZ specifies where the elements of zeta are to
C           be found as follows.
C
C           WHEREZ = 'I' or 'i'
C
C              The elements of zeta are in A.
C
C           WHEREZ = 'S' or 's'
C
C              The elements of zeta are separate from A, in ZETA.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C
C           On entry, M  must specify the number of rows of A. M must be
C           at least n.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C
C           On entry, N  must specify the number of columns of A. N must
C           be  at least zero. When  N = 0  then an immediate return  is
C           effected.
C
C           Unchanged on exit.
C
C  A      - 'real' array of DIMENSION ( LDA, n ).
C
C           Before entry, the leading  M by N  stricly lower  triangular
C           part of the array  A  must contain details of the matrix  Q.
C           In  addition, when  WHEREZ = 'I' or 'i'  then  the  diagonal
C           elements of A must contain the elements of zeta.
C
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C
C           On  entry, LDA  must specify  the leading dimension  of  the
C           array  A  as declared in the calling (sub) program. LDA must
C           be at least m.
C
C           Unchanged on exit.
C
C  ZETA   - 'real' array of DIMENSION at least min( m - 1, n ).
C
C           Before entry with  WHEREZ = 'S' or 's', the array  ZETA must
C           contain the elements of the vector  zeta.
C
C           When  WHEREZ = 'I' or 'i', the array ZETA is not referenced.
C
C           Unchanged on exit.
C
C  NCOLB  - INTEGER.
C
C           On  entry, NCOLB  must specify  the number of columns of  B.
C           NCOLB  must  be  at  least  zero.  When  NCOLB = 0  then  an
C           immediate return is effected.
C
C           Unchanged on exit.
C
C  B      - 'real' array of DIMENSION ( LDB, ncolb ).
C
C           Before entry, the leading  M by NCOLB  part of  the array  B
C           must  contain  the matrix to be  transformed.
C
C           On  exit,  B  is  overwritten  by  the  transformed  matrix.
C
C  LDB    - INTEGER.
C
C           On  entry, LDB  must specify  the  leading dimension of  the
C           array  B as declared in the calling (sub) program. LDB  must
C           be at least m.
C
C           Unchanged on exit.
C
C  WORK   - 'real' array of DIMENSION at least ( ncolb ).
C
C           Used as internal workspace.
C
C  INFORM - INTEGER.
C
C           On  successful exit  INFORM  will be zero, otherwise  INFORM
C           will  be set to unity indicating that an input parameter has
C           been  incorrectly  set. See  the  next  section  for further
C           details.
C
C  4. Diagnostic Information
C     ======================
C
C  INFORM = 1
C
C     One or more of the following conditions holds:
C
C        TRANS  .ne. ' ' or 'N' or 'n' or 'T' or 't' or 'C' or 'c'
C        WHEREZ .ne. 'I' or 'i' or 'S' or 's'
C        M      .lt. N
C        N      .lt. 0
C        LDA    .lt. M
C        NCOLB  .lt. 0
C        LDB    .lt. M
C
C
C  Nag Fortran 77 Auxiliary linear algebra routine.
C
C  -- Written on 15-November-1984.
C     Sven Hammarling, Nag Central Office.
C
      EXTERNAL           DGEMV , DGER
      INTRINSIC          MIN
      INTEGER            J     , K     , KK    , LB
      DOUBLE PRECISION   TEMP
      DOUBLE PRECISION   ONE   ,         ZERO
      PARAMETER        ( ONE   = 1.0D+0, ZERO  = 0.0D+0 )

*     Check the input parameters.

      IF( MIN( N, NCOLB ).EQ.0 )THEN
         INFORM = 0
         RETURN
      END IF
      IF( ( ( TRANS .NE.' ' ).AND.
     $      ( TRANS .NE.'N' ).AND.( TRANS .NE.'n' ).AND.
     $      ( TRANS .NE.'T' ).AND.( TRANS .NE.'t' ).AND.
     $      ( TRANS .NE.'C' ).AND.( TRANS .NE.'c' )      ).OR.
     $    ( ( WHEREZ.NE.'I' ).AND.( WHEREZ.NE.'i' ).AND.
     $      ( WHEREZ.NE.'S' ).AND.( WHEREZ.NE.'s' )      ).OR.
     $    ( M.LT.N ).OR.( N.LT.0 ).OR.( LDA.LT.M ).OR.
     $    ( NCOLB.LT.0 ).OR.( LDB.LT.M )                      )THEN
         INFORM = 1
         RETURN
      END IF

*     Perform the transformation.

      LB = LDB
      DO 20, KK = 1, MIN( M - 1, N )
         IF( ( TRANS.EQ.'T' ).OR.( TRANS.EQ.'t' ).OR.
     $       ( TRANS.EQ.'C' ).OR.( TRANS.EQ.'c' )     )THEN

*           Q'*B = Q( p )*...*Q( 2 )*Q( 1 )*B,     p = min( m - 1, n ).

            K = KK
         ELSE

*           Q*B  = Q( 1 )'*Q( 2 )'*...*Q( p )'*B,  p = min( m - 1, n ).
*           Note that  Q( k )' = Q( k ).

            K = MIN( N, M - 1 ) + 1 - KK
         END IF
         IF( ( WHEREZ.EQ.'S' ).OR.( WHEREZ.EQ.'s' ) )THEN
            TEMP      = A( K, K )
            A( K, K ) = ZETA( K )
         END IF

*        If ZETA( k ) is zero then Q( k ) = I and we can skip the kth
*        transformation.

         IF( A( K, K ).GT.ZERO )THEN
            IF( NCOLB.EQ.1 )
     $         LB = M - K + 1

*           Let C denote the bottom ( m - k + 1 ) by ncolb part of B.

*           First form  work = C'*u.

            DO 10, J = 1, NCOLB
               WORK( J ) = ZERO
   10       CONTINUE
            CALL DGEMV ( 'Transpose', M - K + 1, NCOLB,
     $                   ONE, B( K, 1 ), LB, A( K, K ), 1,
     $                   ZERO, WORK, 1 )

*           Now form  C := C - u*work'.

            CALL DGER  ( M - K + 1, NCOLB, -ONE, A( K, K ), 1,
     $                   WORK, 1, B( K, 1 ), LB )
         END IF

*        Restore the diagonal element of A.

         IF( ( WHEREZ.EQ.'S' ).OR.( WHEREZ.EQ.'s' ) )
     $      A( K, K ) = TEMP
   20 CONTINUE

      INFORM = 0
      RETURN

*     End of DGEAPQ.

      END
