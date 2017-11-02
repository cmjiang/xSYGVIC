!  DORMM
!
!  Purpose:
!  ========
!
!  Pre- and post-multiples a random orthogonal matrix on A and B
!
      subroutine dormm( n, A, B )
!      
      integer, intent(in)             :: n
      double precision, intent(inout) :: A(n, *), B(n, *)
!
      integer                         :: i, info, iseed(4)
      double precision, allocatable   :: d(:), work(:), Q(:,:), C(:,:)
      double precision, parameter     :: one = 1.0d0, zero = 0.0d0
!
      iseed(1) = 100
      iseed(2) = 200
      iseed(3) = 300
      iseed(4) = 401
      allocate(Q(n,n), C(n,n))
      allocate(d(n))
      allocate(work(2*n))
      do i = 1, n
         d(i) = 1
      end do
      call dlagge( n, n, n-1, n-1, d, Q, n, iseed, work, info )
      call dsymm( 'L', 'U', n, n, one, A, n, Q, n, zero, C, n )
      call dgemm( 'T', 'N', n, n, n, one, Q, n, C, n, zero, A, n )
      call dsymm( 'L', 'U', n, n, one, B, n, Q, n, zero, C, n )
      call dgemm( 'T', 'N', n, n, n, one, Q, n, C, n, zero, B, n )
      !call dlarge( n ,A, n, iseed, work2, info )
      !call dlarge( n, B, n, iseed, work2, info )
      deallocate(Q,C,d,work)
!
      end

