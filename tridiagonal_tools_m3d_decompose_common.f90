! subroutine m3d_decompose_r(m,mf,fail)
!   real(rk), intent(in)           :: m (:,:) ! Tridiagonal matrix, assumed to be diagonally dominant
!                                             ! As the result, we do not need to bother with pivoting
!   real(rk), intent(out)          :: mf(:,:) ! Factorized tridiagonal matrix
!   logical, intent(out), optional :: fail    ! Set to .true. if decomposition fails; 
!                                             ! if fail is absent, abort on decomposition failure.
!   !
!   real(rk)    :: denom
    integer(ik) :: i, sz
    !
    sz = size(m,dim=2)
    if (size(m,dim=1)<3 .or. size(mf,dim=1)<3 .or. size(mf,dim=2)/=sz) then
      stop 'tridiagonal_tools%m3d_decompose_common - bad input sizes'
    end if
    if (present(fail)) fail = .false.
    !
    mf(1,1) = 1._rk/m(1,1)
    mf(2,1) = 0._rk
    mf(3,1) = m(3,1)*mf(1,1)
    factor_m3d: do i=2,sz-1
      denom   = m(1,i)-m(2,i-1)*mf(3,i-1)
      if (too_small(abs(denom))) return
      mf(1,i) = 1._rk/denom
      mf(2,i) = -m(2,i-1)*mf(1,i)
      mf(3,i) = m(3,i)*mf(1,i)
    end do factor_m3d
    if (sz<2) return
    denom    = m(1,sz)-m(2,sz-1)*mf(3,sz-1)
    if (too_small(abs(denom))) return
    mf(1,sz) = 1._rk/denom
    mf(2,sz) = -m(2,sz-1)*mf(1,sz)
    mf(3,sz) = 0._rk
    !
    contains
    logical function too_small(x)
      real(kind(m)), intent(in) :: x
      !
      too_small = .false.
      if (x>100*tiny(x)) return
      too_small = .true.
      if (present(fail)) then
        fail = .true.
      else
        write (out,"('Fatal error in m3d_decompose_common: denominator ',g34.16e3,' is too small.')") x
        stop 'tridiagonal_tools%m3d_decompose_common - decomposition failed'
      end if
    end function too_small
! end subroutine m3d_decompose_r
