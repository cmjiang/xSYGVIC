!
!  DGVBOR
!
!  Purpose:
!  ========
!
!  Compute ||X'BX-I|| / ||B||*||X||
!
      subroutine dgvbor( n, k, r, B, X, work )
!
      integer, intent(in)           :: n, k
      double precision, intent(in)  :: B(n,n), X(n,k)
      double precision, intent(out) :: r
      double precision, parameter   :: one = 1.0d0, zero = 0.0d0
      double precision              :: work(*), dlange
!		 
      integer          :: i
      double precision :: rb, rx
      double precision, allocatable :: temp1(:,:), temp2(:,:)
!
      allocate(temp1(n,k), temp2(k,k))
      call dgemm( 'N', 'N', n, k, n, one, B, n, X, n, zero, temp1, n )
      call dgemm( 'T', 'N', k, k, n, one, X, n, temp1, n, zero, temp2, &
                  k )
      do i = 1, k
         temp2(i,i) = temp2(i,i) - 1
      enddo
      r = dlange( 'F', k, k, temp2, k, work )
      rb = dlange( 'F', n, n, B, n, work )
      rx = dlange( 'F', n, k, X, n, work )
      r = r / (rb*rx)
      deallocate(temp1, temp2)
      end
