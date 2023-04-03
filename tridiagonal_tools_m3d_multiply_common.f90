! subroutine m3d_multiply_rr(m,v,mv)
!   real(rk), intent(in)  :: m(:,:) ! Tri-diagonal matrix
!   real(rk), intent(in)  :: v(:)   ! Vector 
!   real(rk), intent(out) :: mv(:)  ! Vector 
    !
    integer(ik) :: sz
    !
    sz = size(v)
    if (size(m,dim=1)<3 .or. size(m,dim=2)/=sz .or. size(mv)/=sz) then
      stop 'tridiagonal_tools%m3d_multiply_common - bad input sizes'
    end if
    mv(1:sz  ) =              m(1,:    )*v(1:sz  ) 
    mv(2:sz  ) = mv(2:sz  ) + m(2,:sz-1)*v(1:sz-1)
    mv(1:sz-1) = mv(1:sz-1) + m(3,:sz-1)*v(2:sz  )
! end subroutine m3d_multiply_rr
