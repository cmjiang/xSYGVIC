!
!  DGVRES
!
!  Purpose:
!  ========
!
!  Compute ||AX-BXW|| / (||A||*||X||*||B||)
!
      subroutine dgvres( n, k, r, A, B, X, w, work )
!
      integer, intent(in)           :: n, k
      double precision, intent(in)  :: A(n,n), B(n,n), X(n,k), w(k)
      double precision, intent(out) :: r
      double precision, parameter   :: one = 1.0d0, zero = 0.0d0
      double precision              :: work(*), dlange
!		 
      integer          :: i
      double precision :: ra, rb, rx, prec
      double precision, allocatable :: temp1(:,:), temp2(:,:)
!
      allocate(temp1(n,k), temp2(n,k))         
      call dsymm( 'L', 'U', n, k, one, A, n, X, n, zero, temp1, n )
      call dsymm( 'L', 'U', n, k, one, B, n, X, n, zero, temp2, n )
      do i = 1, k
         call dscal( n, w(i), temp2(1:n,i), 1 )
      enddo
      temp1 = temp1 - temp2 
      r  = dlange( 'F', n, k, temp1, n, work )
      ra = dlange( 'F', n, n, A, n, work )
      rb = dlange( 'F', n, n, B, n, work )
      rx = dlange( 'F', n, k, X, n, work )
      r = r / (ra*rx*rb)
      deallocate(temp1, temp2)
      end
