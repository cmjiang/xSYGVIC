   program gsep_test
!   
      implicit none
!     
      character(len=32) :: ofile, fileA, fileB
      integer           :: i, j, n, n1, n3, keig(2), info, lwork, na, nb
      integer           :: t1, t2, rate
      double precision  :: e, r1, r2, ro1, ro2, optwork2(1)
      integer, allocatable          :: iwork(:)
      double precision, allocatable :: w1(:), w2(:)
      double precision, allocatable :: work1(:,:), work2(:)
      double precision, allocatable :: A(:,:), B(:,:)
      double precision, allocatable :: A1(:,:), B1(:,:), A2(:,:), B2(:,:)
      integer, parameter            :: out_unit=20
!
      if (iargc() .ne. 3) then
         print *, 'Invalid input'
         stop
      end if
      call getarg(1,fileA)
      call getarg(2,fileB)
      call getarg(3,ofile)
!
      open(unit=10, file = fileA, status='old',action='read')
      open(unit=11, file = fileB, status='old',action='read')
      open(unit=out_unit, file = ofile, status='replace',action='write')
!
      read(10,*) na
      read(11,*) nb
!
      if( na == nb ) then
         n = na
      else
         print *, 'Dimension dismatch'
         stop
      end if
!
      allocate(A(n,n),B(n,n))
      allocate(A1(n,n),B1(n,n),A2(n,n),B2(n,n))
      allocate(w1(n), w2(n))
      allocate(iwork(n))
      allocate(work1(n,n))
!
      do i = 1, n
         read(10,*) A(i,:)
         read(11,*) B(i,:)
      end do
!
      close(10)
      close(11)
!
!     Orthogonal transformation. Not necessary.
!
      call dormm( n, A, B )
!
      write (out_unit,*) ' '
      write (out_unit,*) '========================================'
      write (out_unit,*) '==== Testing Dsygvic Against Dsygv ====='
      write (out_unit,*) '========================================'
      call dprtm10( A, n, out_unit, fileA)
      call dprtm10( B, n, out_unit, fileB)
      A1 = A
      B1 = B
      A2 = A
      B2 = B
!
      e = 1e-12
!
      call dsygvic( 1, 'V', 'L', n, A1, n, B1, n, e, keig, w1, work1, n, &
                    optwork2, -1, iwork, info )  
      lwork = optwork2(1)
      allocate(work2(lwork))
      call system_clock(t1)
      call dsygvic( 1, 'V', 'L', n, A1, n, B1, n, e, keig, w1, work1, n, &
                    work2, lwork, iwork, info )
      call system_clock(t2, rate)
      write (out_unit,*) 'dsygvic'
      write (out_unit,*) ' Info:', info
      write (out_unit,*) ' Number of eigvals:', keig(1)
      write (out_unit,*) ' Case number:', keig(2)
      write (out_unit,*) ' Running time:', (t2-t1)/REAL(rate), '(seconds)'
!
      if (keig(1) > 0) then
         call system_clock(t1)
         call dsygv( 1, 'V', 'L', n, A2, n, B2, n ,w2 ,work2, lwork, &
                     info )
         call system_clock(t2, rate)
         call dgvres( n, keig(1), r1, A, B, A1, w1, work2 )
         call dgvres( n, n, r2, A, B, A2, w2, work2 )
         call dgvbor( n, keig(1), ro1, B, A1, work2 )
         call dgvbor( n, n, ro2, B, A2, work2 )
         write (out_unit,*) ' dsygvic residual:', r1
         write (out_unit,*) ' dsygvic b-orthog:', ro1
         write (out_unit,*) '----------------------------------------'
         write (out_unit,*) 'dsygv'
         write (out_unit,*) ' Info:', info
         write (out_unit,*) ' Running time:',  (t2-t1)/REAL(rate), &
                              '(seconds)'
         write (out_unit,*) ' dsygv residual:', r2
         write (out_unit,*) ' dsygv b-orthog:', ro2
         write (out_unit,*) '-----------eig dsygvic------------------'
         write (out_unit,'(e40.30)') w1(1:keig(1)) 
         write (out_unit,*) '-----------eig dsygv--------------------'
         write (out_unit,'(e40.30)') w2(1:n)
      end if
      close (out_unit)
      deallocate( A, B, A1, B1, A2, B2 )
      deallocate( w1, w2, work1, work2, iwork )
   end
