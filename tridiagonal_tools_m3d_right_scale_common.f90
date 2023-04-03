! subroutine m3d_right_scale_rc(m,s,ms)
!   real(rk), intent(in)     :: m (:,:) ! Tridiagonal matrix
!   complex(rk), intent(in)  :: s   (:) ! Diagonal matrix
!   complex(rk), intent(out) :: ms(:,:) ! Tridiagonal m . s
!   !
    integer(ik) :: sz
    !
    sz = size(m,dim=2)
    if (size(m,dim=1)<3 .or. size(ms,dim=1)<3 .or. size(s)/=sz .or. size(m,dim=2)/=sz) then
      stop 'tridiagonal_tools%m3d_right_scale_common - bad input sizes'
    end if
    !
    ms(1,1:sz  ) = m(1,1:sz  ) * s(1:sz  )
    ms(2,1:sz-1) = m(2,1:sz-1) * s(1:sz-1)
    ms(3,1:sz-1) = m(3,1:sz-1) * s(2:sz  )
! end subroutine m3d_right_scale_rc
