!  DLAASRT
!
!  Purpose:
!  ========
!
!  Sort the array D(1:N) in absolute-value descending order and
!  record the permutation in array P(1:N) s.t. new D(j) = old D(P(j))
!
!  Arguments
!  =========
!
!  N        (input) INTEGER
!           The length of the array D.
!
!  D        (input/output) DOUBLE PRECISION array, dimension (N)
!           On entry, the array to be sorted.
!           On exit, D has been sorted into absolute-value descending 
!           order.
!
!  P        (output) INTEGER array, dimension (N)
!           P contains the permutation vector.
!
      SUBROUTINE DLAASRT( N, D, P )
!
      INTEGER, INTENT(IN) :: N
      DOUBLE PRECISION, INTENT(INOUT) :: D(N)
      INTEGER, INTENT(OUT) :: P(N)
!
      INTEGER :: I, J, TMP2
      DOUBLE PRECISION :: TMP1
      LOGICAL :: SWAPPED
!  
      SWAPPED = .TRUE.
      J = 0
      DO I = 1,N
         P(I) = I
      ENDDO
      DO WHILE (SWAPPED)
        SWAPPED = .FALSE.
        J = J+1
        DO I = 1, N-J
           IF(ABS(D(I))<ABS(D(I+1))) THEN
              TMP1 = D(I)
              TMP2 = P(I)
              D(I) = D(I+1)
              P(I) = P(I+1)
              D(I+1) = TMP1
              P(I+1) = TMP2
              SWAPPED = .TRUE.     
           END IF
        END DO
      END DO
      END


