!
!  DPRTM10
!
!  Purpose:
!  ========
!
!  Print the first 10 rows and columns of an n by n matrix
!
      subroutine dprtm10( A, n, out_unit, fileA )
      character(len=32), intent(in) :: fileA
      integer, intent(in)           :: n, out_unit
      double precision, intent(in)  :: A(n,*)
!
      integer :: i
!
      write (out_unit,*) fileA
      do i = 1,10
         write (out_unit, '(10f8.3)') A(i,1), A(i,2), A(i,3), A(i,4), &
            A(i,5), A(i,6), A(i,7), A(i,8), A(i,9), A(i,10)
      end do
      write (out_unit,*) '----------------------------------------'
      end
