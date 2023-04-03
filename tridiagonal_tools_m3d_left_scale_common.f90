! subroutine m3d_left_scale_cr(s,m,sm)
!   complex(rk), intent(in)  :: s   (:) ! Diagonal matrix
!   real(rk), intent(in)     :: m (:,:) ! Tridiagonal matrix
!   complex(rk), intent(out) :: sm(:,:) ! Tridiagonal s . m
    !
    integer(ik) :: sz
    !
    sz = size(m,dim=2)
    if (size(m,dim=1)<3 .or. size(sm,dim=1)<3 .or. size(s)/=sz .or. size(m,dim=2)/=sz) then
      stop 'tridiagonal_tools%m3d_left_scale_common - bad input sizes'
    end if
    !
    sm(1,1:sz  ) = s(1:sz  ) * m(1,1:sz  )
    sm(2,1:sz-1) = s(2:sz  ) * m(2,1:sz-1)
    sm(3,1:sz-1) = s(1:sz-1) * m(3,1:sz-1)
! end subroutine m3d_left_scale_cr
