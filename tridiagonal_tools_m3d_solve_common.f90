! subroutine m3d_solve_rr(mf,r,x)
!   real(rk), intent(in)  :: mf(:,:) ! Factorized tridiagonal matrix
!   real(rk), intent(in)  :: r (:)   ! Right-hand size
!   real(rk), intent(out) :: x (:)   ! Solution vector
    !
    integer(ik) :: i, sz
    !
    sz = size(mf,dim=2)
    if (size(mf,dim=1)<3 .or. size(r)/=sz .or. size(x)/=sz) then
      stop 'tridiagonal_tools%m3d_solve_common - bad input sizes'
    end if
    !
    x(1) = mf(1,1)*r(1)
    transform_rhs: do i=2,sz
      x(i) = mf(1,i)*r(i) + mf(2,i)*x(i-1)
    end do transform_rhs
    !
    backsubstitute: do i=sz-1,1,-1
      x(i) = x(i) - mf(3,i)*x(i+1)
    end do backsubstitute
! end subroutine m3d_solve_rr
