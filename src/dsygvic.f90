!  DSYGVIC
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSYGVIC( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, &
!                           ETOL, K, W, WORK, LDWORK, WORK2, LWORK, &
!                           IWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER, INTENT(IN)           :: JOBZ, UPLO
!       INTEGER, INTENT(IN)             :: N, ITYPE, LDA, LDB
!       INTEGER, INTENT(IN)             :: LDWORK, LWORK
!       DOUBLE PRECISION, INTENT(IN)    :: ETOL
!       INTEGER, INTENT(OUT)            :: K(2), INFO
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION, INTENT(INOUT) :: A(LDA, *), B(LDB, *)
!       DOUBLE PRECISION, INTENT(OUT)   :: W(*)
!       INTEGER                         :: IWORK(N)
!       DOUBLE PRECISION                :: WORK(LDWORK, *), WORK2(*)
!       ..
!
!  Purpose:
!  ========
!
!  DSYGVIC computes the e-finite eigenvalues and the corresponding
!  eigenvectors of the real generalized positive semi-definite
!  eigenproblem, A*x=(lambda)*B*x, where A and B are symmetric and 
!  B is positive semi-definite with respect to the threshold ETOL.
!  definite.
!
!  Arguments:
!  ==========
!
!  ITYPE    (input) INTEGER
!           Specifies the problem type to be solved:
!           = 1: A*x = (lambda)*B*x
!
!  JOBZ     (input) CHARACTER*1
!           = 'V': Compute eigenvalues and eigenvectors.
!
!  UPLO     (input) CHARACTER*1
!           = 'U': Upper triangles of A and B are stored;
!           = 'L': Lower triangles of A and B are stored.
!
!  N        (input) INTEGER
!           The order of the matrices A and B. N > 0.
!
!  A        (input/output) DOUBLE PRECISION array, dimension (LDA, N)
!           On entry, the symmetric matrix A. If UPLO = 'U', the
!           leading N-by-N upper triangular part of A contains the
!           upper triangular part of the matrix A. If UPLO = 'L',
!           the leading N-by-N lower triangular part of A contains
!           the lower triangular part of the matrix A.
!           On exit, if JOBZ = 'V', then if INFO = k > 0, the first
!           K(1) columns of A contains the matrix X of eigenvectors. 
!
!  LDA      (input) INTEGER
!           The leading dimension of the array A. LDA >= max(1,N).
!
!  B        (input/output) DOUBLE PRECISION array, dimension (LDB, N)
!           On entry, the symmetric semi-positive definite matrix B.
!           If UPLO = 'U', the leading N-by-N upper triangular part of B
!           contains the upper triangular part of the matrix B.
!           If UPLO = 'L', the leading N-by-N lower triangular part of B
!           contains the lower triangular part of the matrix B.
!           On exit, B contains the transformation matrix Q_1*R_1*Q_2*Q_3,
!           depending on the exit stage.
!
!  LDB      (input) INTEGER
!           The leading dimension of the array B.  LDB >= max(1,N).
!
!  ETOL     (input) DOUBLE PRECISION
!           The parameter used to drop small eigenvalues of B.
!
!  K        (output) INTEGER, dimension (2)
!           K(1) indicates the number of finite eigenvalues if INFO = 0.
!           K(2) indicates the case number of different branches of 
!           exits. Specifically, if
!           K(1) > 0  : W contains K(1) eigenvalues; K(2) is 1, 2, 3 or 4.
!           K(1) = 0  : regular pencil but no finite eigenvalue; K(2) is
!                       1, 2 or 3.
!           K(1) = -1 : singular pencil; K(2) is 1, 2, ... , 6 or 7
!           
!
!  W        (output) DOUBLE PRECISION array, dimension (N)
!           If K(1) > 0, W stores the K(1) eigenvalues.
!
!  WORK     (workspace) DOUBLE PRECISION array, dimension (N, N)
!
!  LDWORK   (input) INTEGER
!           The leading dimension of work.  LDWORK >= max(1,N).
!        
!  WORK2    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!           On exit, if LWORK = -1 and INFO = 0, WORK2(1) returns the 
!           optimal size of WORK2.
!
!  LWORK    (input) INTEGER
!           The length of the array WORK2.  LWORK >= 3*N+1.
!           For optimal efficiency, LWORK >= 2*N+( N+1 )*NB,
!           where NB is the blocksize.
!           If LWORK = -1, then a workspace query is assumed; the 
!           routine only calculates the optimal size of the WORK2 array,
!           returns this value as the first entry of the WORK2 array, 
!           and no error message related to LWORK is issued by XERBLA.
!
!  IWORK    (workspace) INTEGER array, dimension (N)
!
!  INFO     (output) INTEGER
!                  INFO  = 0    : successful exit.       
!           -17 <= INFO <= -1   : if INFO = -i, the i-th argument had an 
!                                 illegal value.
!
!  Date: April 2015
!
!  =====================================================================
      SUBROUTINE DSYGVIC( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, ETOL, &
                          K, W, WORK, LDWORK, WORK2, LWORK, IWORK, &
                          INFO )
!
!     .. Scalar Arguments ..
      CHARACTER, INTENT(IN)           :: JOBZ, UPLO
      INTEGER, INTENT(IN)             :: N, ITYPE, LDA, LDB
      INTEGER, INTENT(IN)             :: LDWORK, LWORK
      DOUBLE PRECISION, INTENT(IN)    :: ETOL
      INTEGER, INTENT(OUT)            :: K(2), INFO
!       ..
!       .. Array Arguments ..
      DOUBLE PRECISION, INTENT(INOUT) :: A(LDA, *), B(LDB, *)
      DOUBLE PRECISION, INTENT(OUT)   :: W(*)
      INTEGER                         :: IWORK(N)
      DOUBLE PRECISION                :: WORK(LDWORK, *), WORK2(*)
!       ..
!
!  =====================================================================
!
!     .. Parameter ..
      DOUBLE PRECISION, PARAMETER ::  ONE = 1.0d0, ZERO = 0.0d0
!     ..
!     .. Local Scalar ..
      LOGICAL            :: LQUERY, UPPER, WANTZ
      INTEGER            :: NB, LWKOPT, I, N1, N2, N3, N4, N5
      DOUBLE PRECISION   :: MCHEPS, DLAMCH
!     ..
!     .. External Functions ..
      LOGICAL            :: LSAME
      INTEGER            :: ILAENV
      EXTERNAL           :: LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           :: ABSSORT, DGEMM, DGEQP3, DLAPMT, DLAMCH
      EXTERNAL           :: DORMQR,  DSCAL, DSYEV,  DTRMM,  DTRTRI
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
      WANTZ = LSAME( JOBZ, 'V' )
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
!
      INFO = 0
      K(1) = 0
      K(2) = 0
!
      IF( ITYPE.NE.1 ) THEN
         INFO = -1
      ELSE IF( .NOT.( WANTZ ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( UPPER .OR. LSAME( UPLO, 'L' ) ) ) THEN
         INFO = -3
      ELSE IF( N.LE.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( ETOL.LT.0 ) THEN
         INFO = -9
      ELSE IF ( LDWORK.LT.MAX( 1, N ) ) THEN
         INFO = -13
      END IF
!
      IF( INFO.EQ.0 ) THEN
!
         CALL DGEQP3( N, N, A, LDA, IWORK, W, WORK2, -1, INFO )
!
         IF( LWORK.LT.MAX( 1, 3*N+1 ) .AND. .NOT.LQUERY ) THEN
           INFO = -16
         END IF
      END IF
!
      IF( INFO.NE.0 ) THEN
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Get machine precision
!
      MCHEPS = DLAMCH('E')
!
!  =====================================================================
!                            Phase 1.
!  ---------------------------------------------------------------------
!
!     Diagonalize B by eigenvalue decomposition
!
      CALL DSYEV( 'V', UPLO, N, B, LDB, W, WORK2, LWORK, INFO )
      IF ( INFO .NE. 0 ) THEN
         RETURN
      END IF
!
!     Check if B is not semi-positive definite.
!
      IF ( W(1) <= -ETOL*W(N) ) THEN
         INFO = -7
         RETURN
      END IF
!
!     Sort eigenvalues and eigenvectors in descending order
! 
      DO I = 1, N
         IWORK(I) = N-I+1
      END DO
      CALL DLAPMT( .FALSE., 1, N, W, 1, IWORK )
      CALL DLAPMT( .FALSE., N, N, B, LDB, IWORK )
!
!     Find n1 and n2
!
      IF ( W(1) <= MCHEPS ) THEN
        N1 = 0
        N2 = n
      ELSE
         I = 1
         DO WHILE ( I < N+1 .AND. W(I) > ETOL*W(1) )
            I = I+1
         END DO
         N1 = I-1
         N2 = N-N1
      END IF
!
!     Early Exit -- n1 = 0.
!
      IF ( N1 == 0 ) THEN
         CALL DSYEV( 'N', UPLO, N, A, LDB, W, WORK2, LWORK, INFO )
         IF ( INFO .NE. 0 ) THEN
            RETURN
         END IF
         CALL DLAASRT( N, W, IWORK )
         IF ( W(N) < MCHEPS ) THEN
!           Singular pencil
            K(1) = -1
            K(2) = 1
            RETURN
         ELSE
!           Regular pencil but no finite eigenvalue
            K(2) = 1
            RETURN
         END IF
      END IF
!
!     Update A to A0
!
      CALL DSYMM( 'L', UPLO, N, N, ONE, A, LDA, B, LDB, ZERO, WORK, &
                  LDWORK )
      CALL DGEMM( 'T', 'N', N, N, N, ONE, B, LDB, WORK, LDWORK, ZERO, &
                  A, N )
!
!     Reduce B to identity and zero.
!
      W(1:N1) = 1/SQRT(W(1:N1))
      DO I = 1,N1
         CALL DSCAL( N, W(I), B(1:N,I), 1 )
      END DO
!
!     Update A to A1
!
      DO I = 1,N1
         CALL DSCAL( N, W(I), A(1:N,I), 1 )
      END DO
      DO I = 1,N1
         CALL DSCAL( N, W(I), A(I,1:N), 1 )
      END DO
!
!     Early Exit: n2 = 0, B is e-well conditioned.
!                 #N e-finite eigenvalues!       
!
      IF (N2 == 0) THEN
         CALL DSYEV( 'V', UPLO, N, A, LDA, W, WORK2, LWORK, INFO )
         IF (INFO .NE. 0 ) THEN
            RETURN
         END IF
         CALL DGEMM( 'N', 'N', N, N, N, ONE, B, LDB, A, LDA, ZERO, &
                     WORK, LDWORK )
         A(1:N,1:N) = WORK(1:N,1:N)
         K(1) = N
         K(2) = 1
         RETURN
      END IF
!
!  =====================================================================
!                            Phase 2.
!  ---------------------------------------------------------------------
!
!     Diagonalize A1_22 by eigenvalue decomposition
!
      CALL DSYEV( 'V', UPLO, N2, A(N1+1,N1+1), LDA, W, WORK2, LWORK, &
                  INFO )
      IF ( INFO .NE. 0 ) THEN
         RETURN
      END IF
!
!     Sort eigenvalues and eigenvectors in abs-descending order
!  
      CALL DLAASRT( N2, W, IWORK )
      CALL DLAPMT( .TRUE., N2, N2, A(N1+1,N1+1), LDA, IWORK )
!
!     Find n3 and n4
!
      IF ( ABS(W(1)) <= ETOL ) THEN
         N3 = 0
         N4 = N2
      ELSE
         I = 1
         DO WHILE ( I < N2+1 .AND. ABS(W(I)) > ETOL*ABS(W(1)) )
            I = I+1
         ENDDO
         N3 = I-1
         N4 = N2-N3
      END IF
!
!     Early Exit: n3 = 0.
!
      IF ( N3 == 0 ) THEN
!
!        If n1 < n2, singular pencil, exit. 
!
         IF ( N1 < N2 ) THEN
            K(1) = -1
            K(2) = 2
            RETURN
         END IF
!
!        Reveal the rank of A1_12 by QR decomposition with pivoting
!
         IWORK = 0
         CALL DGEQP3( N1, N2, A(1,N1+1), LDA, IWORK, W, WORK2, LWORK, &
                      INFO )
         IF ( INFO .NE. 0 ) THEN
            RETURN
         END IF
!
!        If n1 = n2:
!           if A1_12 is rank deficient, then singular pencil, exit;
!           else regular pencil but no e-finite eigenvalues, exit.
!
         IF ( N1 == N2 ) THEN
            IF( ABS(A(N2,N)) <= MCHEPS ) THEN
               K(1) = -1
               K(2) = 3
               RETURN
            ELSE
               K(2) = 2
               RETURN
            END IF
         END IF
!
!        If n1 > n2:
!           if A2_12 is rank deficient, then singular pencil, exit;
!           else n1-n2 e-finite eigenvalues.
!
         IF( ABS(A(N2,N)) <= MCHEPS ) THEN
            K(1) = -1
            K(2) = 4
            RETURN
         ELSE
!
!           Update A1_11
!
            CALL DORMQR( 'R', 'N', N1, N1, N2, A(1,N1+1), LDA, W, A, &
                         LDA, WORK2, LWORK, INFO )
            IF ( INFO .NE. 0 ) THEN
               RETURN
            END IF
            CALL DORMQR( 'L', 'T', N1, N1, N2, A(1,N1+1), LDA, W, A, &
                         LDA, WORK2, LWORK, INFO )
            IF ( INFO .NE. 0 ) THEN
               RETURN
            END IF
!
!           B <- Q_1 * R_1 * Q_2
!
            CALL DLAPMT( .TRUE., N, N2, B(1,N1+1), LDB, IWORK )
            CALL DORMQR( 'R', 'N', N, N1, N2, A(1,N1+1), LDA, W, B, &
                         LDB, WORK2, LWORK, INFO )
            IF ( INFO .NE. 0 ) THEN
               RETURN
            END IF
! 
!           Compute U_2
!
            CALL DSYEV( 'V', UPLO, N1-N2, A(N2+1,N2+1), LDA, W, WORK2, &
                        LWORK, INFO )
            IF ( INFO .NE. 0 ) THEN
               RETURN
            END IF
            WORK(N2+1:N1,N2+1:N1) = A(N2+1:N1,N2+1:N1)
!
!           Compute U_3 = -A2_13^(-1) * A2_12 * U_2
!              (1) Compute A2_12 * U_2
!              (2) Compute A2_13^(-1)
!              (3) Compute -A2_13^(-1) * A2_12 * U_2
!
            CALL DGEMM( 'N', 'N', N2, N1-N2, N1-N2, ONE, A(1,N2+1), &
                        LDA, WORK(N2+1,N2+1), LDWORK, ZERO, &
                        WORK(N1+1,N2+1), LDWORK )
!
            CALL DTRTRI( 'U', 'N', N2, A(1,N1+1), N, INFO )
            IF ( INFO .NE. 0 ) THEN
               RETURN
            END IF
!           
            CALL DTRMM( 'L', 'U', 'N', 'N', N2, N1-N2, -ONE, &
                        A(1,N1+1), LDA, WORK(N1+1,N2+1), LDWORK )
!
!           Set U_1 = 0 and compute X = Q_1 * R_1 * Q_2 * U
!
            WORK(1:N2, N2+1:N1) = 0
            CALL DGEMM( 'N', 'N', N, N1-N2, N, ONE, B, LDB, &
                        WORK(1,N2+1), LDWORK, ZERO, A, LDA )
!
!           Output #(N1-N2) eigenvalues
!
            K(1) = N1 - N2
            K(2) = 2
            RETURN
         END IF
      END IF             
!
!     Update A:
!        (1) Update A_12 = A_12*Q_22
!        (2) Update A_21 = Q_22'*A_21
!
      CALL DGEMM( 'N', 'N', N1, N2, N2, ONE, A(1,N1+1), LDA, &
                  A(N1+1,N1+1), LDA, ZERO, WORK(1,N1+1), LDWORK )
      A(1:N1,N1+1:N) = WORK(1:N1,N1+1:N)
!     
      CALL DGEMM( 'T', 'N', N2, N1, N2, ONE, A(N1+1,N1+1), LDA, &
                  A(N1+1,1), LDA, ZERO, WORK(N1+1,1), LDWORK )
      A(N1+1:N,1:N1) = WORK(N1+1:N,1:N1)
!
!     B <- Q_1*R_1*Q_2
!
      CALL DGEMM( 'N', 'N', N, N2, N2, ONE, B(1,N1+1), LDB, &
                  A(N1+1,N1+1), LDA, ZERO, WORK(1,N1+1), LDWORK )
      B(:,N1+1:N) = WORK(:,N1+1:N)
!
!     Early Exit: n4 = 0.
!
      IF (N4 == 0) THEN
!
!        Compute A2_11-A2_12/D2*A2_12'
!        (1) Compute D^(-1)*A2_12' (row scaling)
!        (2) Compute A2_12*D^(-1)*A2_12'
!        (3) Compute A2_11-A2_12/D*A2_12'
!
         W(1:N2) = 1/W(1:N2)     
         DO I = 1,N2
            CALL DSCAL( N1, W(I), A(N1+I,1:N1), 1 )
         END DO
!		 
         CALL DGEMM( 'N', 'N', N1, N1, N2, ONE, A(1,N1+1), LDA, &
                     A(N1+1,1), LDA, ZERO, WORK, LDWORK )
!					 
         A(1:N1,1:N1) = A(1:N1,1:N1) - WORK(1:N1,1:N1)
!
!        Compute U_1 and U_2 = -D^(-1) * A2_12' * U_1
!
         CALL DSYEV( 'V', UPLO, N1, A, N, W, WORK2, LWORK, INFO )
         IF ( INFO .NE. 0 ) THEN
            RETURN
         END IF
         CALL DGEMM( 'N', 'N', N2, N1, N1, -ONE, A(N1+1,1), LDA, A, &
                     LDA, ZERO, WORK(N1+1,1), LDWORK )
!
!        Compute X = Q_1 * R_1 * Q_2 * U and output #eigenvalues
!
         WORK(1:N1,1:N1) = A(1:N1,1:N1)
         CALL DGEMM( 'N', 'N', N, N1, N, ONE, B, LDB, WORK, LDWORK, &
                     ZERO, A, LDA )
!
!        Output #N1 eigenvalues
!
         K(1) = N1
         K(2) = 3
         RETURN
      END IF
!
!  =====================================================================
!                            Phase 3.
!  ---------------------------------------------------------------------
!
!     Since n3 ~= 0 and n4 ~=0, we can view A2 as a 3 by 3 block:
!
!                 n1     n3     n4
!          n1 | A2_11  A2_12  A2_13 |
!     A2 = n3 | A2_12'   D2     0   |
!          n4 | A2_13'   0      0   |
!
!     Early exit: n1 < n4, singular pencil, exit.
! 
      IF ( N1 < N4 ) THEN
         K(1) = -1
         K(2) = 5
         RETURN
      END IF
!
!     If n1 >= n4
!        reveal the rank of A2_13 by QR with column pivoting
!
      A(N1+1:N1+N3,N1+1) = W(1:N3)
!
      IWORK = 0
      CALL DGEQP3( N1, N4, A(1,N1+N3+1), LDA, IWORK, W, WORK2, LWORK, &
                   INFO )
      IF ( INFO .NE. 0 ) THEN
         RETURN
      END IF
!
!     If n1 = n4:
!        if R3_13 rank deficient, singular pencil,
!        else regular pencil but no e-finite eigenvalues.
!
      IF( N1 == N4 ) THEN    
         IF ( ABS(A(N4,N)) <= MCHEPS ) THEN
            K(1) = -1
            K(2) = 6
            RETURN
         ELSE
            K(2) = 3
            RETURN
         END IF
      END IF
!
!     n1 > n4:
!        if R3_13 rank deficient, singular pencil,
!        else n1-n4 e-finite eigenvalues.
!
      IF( ABS(A(N4,N)) <=  MCHEPS ) THEN
          K(1) = -1
          K(2) = 7
          RETURN
      END IF
      
!     Update A to A3
!     (1) A2_11 = Q3_11'*A2_11*Q3_11
!     (2) A2_12 = Q3_11'*A2_12
!     (3) A2_21 = A2_21*Q3_11
!
!     A3 is a 4 by 4 block matrix:
!
!                 n1     n5     n3     n4
!          n1 | A3_11  A3_12  A3_13  A3_14 |
!     A3 = n5 | A3_12' A3_22  A3_23    0   |
!          n3 | A3_13' A3_23' D2_33    0   |
!          n4 | A3_14'   0      0      0   |
!     
      CALL DORMQR( 'R', 'N', N1, N1, N4, A(1,N1+N3+1), LDA, W, &
                   A, LDA, WORK2, LWORK, INFO )
      IF ( INFO .NE. 0 ) THEN
         RETURN
      END IF
      CALL DORMQR( 'L', 'T', N1, N1, N4, A(1,N1+N3+1), LDA, W, &
                   A, LDA, WORK2, LWORK, INFO )
      IF ( INFO .NE. 0 ) THEN
         RETURN
      END IF
!     
      CALL DORMQR( 'L', 'T', N1, N3, N4, A(1,N1+N3+1), LDA, W, &
                   A(1,N1+1), LDA, WORK2, LWORK, INFO )
      IF ( INFO .NE. 0 ) THEN
         RETURN
      END IF
!
      CALL DORMQR( 'R', 'N', N3, N1, N4, A(1,N1+N3+1), LDA, W, &
                   A(N1+1,1), LDA, WORK2, LWORK, INFO )
      IF ( INFO .NE. 0 ) THEN
         RETURN
      END IF
!
!     B <- Q_1 * R_1 * Q_2 * P_3 * Q_3
!
      CALL DLAPMT( .TRUE., N, N4, B(1,N1+N3+1), LDB, IWORK )
      CALL DORMQR( 'R', 'N', N, N1, N4, A(1,N1+N3+1), LDA, W, &
                   B, LDB, WORK2, LWORK, INFO )
      IF ( INFO .NE. 0 ) THEN
         RETURN
      END IF
!
!     Compute A3_22-A3_23/D2*A3_23'
!     (1) Compute D^(-1)*A3_23' --> A_32 (row scaling)
!     (2) Compute A3_23/D*A3_23'
!     (3) Compute A3_22-A3_23/D*A3_23'
!
      W(1:N3) = A(N1+1:N1+N3,N1+1)
      DO I = 1, N3
         W(I) = 1 / W(I)
      END DO
      N5 = N1 - N4
      DO I = 1, N3
         CALL DSCAL( N5, W(I), A(N1+I,N4+1:N1), 1 )
      END DO
!     
      CALL DGEMM( 'N', 'N', N5, N5, N3, ONE, A(N4+1,N1+1), LDA, &
                  A(N1+1,N4+1), LDA, ZERO, WORK(N4+1,N4+1), LDWORK )
!     
      A(N4+1:N1,N4+1:N1) = A(N4+1:N1,N4+1:N1) - WORK(N4+1:N1,N4+1:N1)
!
!     Compute U_2
!
      CALL DSYEV( 'V', UPLO, N5, A(N4+1,N4+1), LDA, W, WORK2, LWORK, &
                  INFO )
      IF ( INFO .NE. 0 ) THEN
         RETURN
      END IF
      WORK(N4+1:N1,N4+1:N1) = A(N4+1:N1,N4+1:N1)
!
!     Compute U_3 = -D \ A3_23' * U_2
!
      CALL DGEMM( 'N', 'N', N3, N5, N5, -ONE, A(N1+1,N4+1), LDA, &
                  A(N4+1,N4+1), LDA, ZERO, WORK(N1+1,N4+1), LDWORK )
!
!     Compute U_4 = -A3_14^(-1) * ( A3_12 * U_2 + A3_13 * U_3 )
!     (1) Compute A3_12 * U_2
!     (2) Compute A3_13 * U_3
!     (3) Compute A3_12*U_2 + A3_13*U_3
!     (4) Compute A3_14^(-1)
!     (5) Compute -A3_14^(-1) * ( A3_12 * U_2 + A3_13 * U_3 )
!
      CALL DGEMM( 'N', 'N', N4, N5, N5, ONE, A(1,N4+1), LDA, &
                  WORK(N4+1,N4+1), LDWORK, ZERO, WORK(1,N4+1), LDWORK )
!
      CALL DGEMM( 'N', 'N', N4, N5, N3, ONE, A(1,N1+1), LDA, &
                  WORK(N1+1,N4+1), LDWORK, ZERO, WORK(N1+N3+1,N4+1), &
                  LDWORK )
!     
      A(N1+N3+1:N,N4+1:N1) = WORK(N1+N3+1:N,N4+1:N1) + &
                                 WORK(1:N4,N4+1:N1)
      WORK(N1+N3+1:N,N4+1:N1) = A(N1+N3+1:N,N4+1:N1)
! 
      CALL DTRTRI( 'U', 'N', N4, A(1,N1+N3+1), N, INFO )
      IF ( INFO .NE. 0 ) THEN
         RETURN
      END IF
! 
      CALL DTRMM( 'L', 'U', 'N', 'N', N4, N5, -ONE, A(1,N1+N3+1), &
                  LDA, WORK(N1+N3+1,N4+1), LDWORK )
!
!     Set U_1 = 0 and compute X = Q_1 * R_1 * Q_2 * P_3 * Q_3 * U
!
      WORK(1:N4, N4+1:N1) = 0
      CALL DGEMM( 'N', 'N', N, N5, N, ONE, B, LDB, WORK(1,N4+1), &
                  LDWORK, ZERO, A, LDA )
!
!     Output #n5 eigenvalues
!
      K(1) = N5
      K(2) = 4
!
!     End of DSYGVIC
!	  
      END
