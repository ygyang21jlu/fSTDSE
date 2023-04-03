! subroutine m3d_multiply_rr_g(m,v,mv)
!   real(rk), intent(in)  :: m(:,:)  ! Tri-diagonal matrix
!   real(rk), intent(in)  :: v(:,:)  ! Vectors
!   real(rk), intent(out) :: mv(:,:) ! Vectors
    !
    integer(ik) :: sz, iv
    !
    sz = size(v,dim=1)
    if (size(m,dim=1)<3 .or. size(m,dim=2)/=sz .or. size(mv,dim=1)/=sz .or. size(v,dim=2)/=size(mv,dim=2)) then
      stop 'tridiagonal_tools%m3d_multiply_g_common - bad input sizes'
    end if
    scan_vectors: do iv=1,size(v,dim=2)
      mv(1:sz  ,iv) =                 m(1,:    )*v(1:sz  ,iv) 
      mv(2:sz  ,iv) = mv(2:sz  ,iv) + m(2,:sz-1)*v(1:sz-1,iv)
      mv(1:sz-1,iv) = mv(1:sz-1,iv) + m(3,:sz-1)*v(2:sz  ,iv)
    end do scan_vectors
! end subroutine m3d_multiply_rr_g
