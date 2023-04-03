!
!  Simple tools for dealing with three-diagonal matrices and systems of equations.
!
!  All tridiagonal matrices used by this module are in the following format:
!   m(1,i) = Main diagonal; element (i,i)
!   m(2,i) = Sub-diagonal; element (i+1,i)
!   m(3,i) = Super-diagonal; element (i,i+1)
!
module tridiagonal_tools
  use accuracy
  use constants
  use timer
  implicit none
  private
  public m3d_multiply, m3d_decompose, m3d_solve
  public m3d_multiply_x, m3d_decompose_x, m3d_solve_x
  public m3d_left_scale, m3d_right_scale, m3d_transpose
  !
  interface m3d_multiply
    module procedure m3d_multiply_rr
    module procedure m3d_multiply_rc
    module procedure m3d_multiply_cc
    module procedure m3d_multiply_rr_g
    module procedure m3d_multiply_rc_g
    module procedure m3d_multiply_cc_g
  end interface m3d_multiply
  interface m3d_multiply_x
    module procedure m3d_multiply_xrr_g
  end interface m3d_multiply_x
  !
  interface m3d_decompose
    module procedure m3d_decompose_r
    module procedure m3d_decompose_c
  end interface m3d_decompose
  interface m3d_decompose_x
    module procedure m3d_decompose_xr
  end interface m3d_decompose_x
  !
  interface m3d_solve
    module procedure m3d_solve_rr
    module procedure m3d_solve_rc
    module procedure m3d_solve_cc
    module procedure m3d_solve_rr_g
    module procedure m3d_solve_rc_g
    module procedure m3d_solve_cc_g
  end interface m3d_solve
  interface m3d_solve_x
    module procedure m3d_solve_xrr_g
  end interface m3d_solve_x
  !
  interface m3d_left_scale
    module procedure m3d_left_scale_cr
  end interface m3d_left_scale
  !
  interface m3d_right_scale
    module procedure m3d_right_scale_rc
  end interface m3d_right_scale
  !
  interface m3d_transpose
    module procedure m3d_transpose_r
  end interface m3d_transpose
  !
  contains
  !
  subroutine m3d_multiply_rr(m,v,mv)
    real(rk), intent(in)  :: m(:,:) ! Tri-diagonal matrix
    real(rk), intent(in)  :: v(:)   ! Vector 
    real(rk), intent(out) :: mv(:)  ! Vector 
    !
    include "tridiagonal_tools_m3d_multiply_common.f90"
  end subroutine m3d_multiply_rr
  !
  subroutine m3d_multiply_rc(m,v,mv)
    real(rk), intent(in)     :: m(:,:) ! Tri-diagonal matrix
    complex(rk), intent(in)  :: v(:)   ! Vector 
    complex(rk), intent(out) :: mv(:)  ! Vector 
    !
    include "tridiagonal_tools_m3d_multiply_common.f90"
  end subroutine m3d_multiply_rc
  !
  subroutine m3d_multiply_cc(m,v,mv)
    complex(rk), intent(in)  :: m(:,:) ! Tri-diagonal matrix
    complex(rk), intent(in)  :: v(:)   ! Vector 
    complex(rk), intent(out) :: mv(:)  ! Vector 
    !
    include "tridiagonal_tools_m3d_multiply_common.f90"
  end subroutine m3d_multiply_cc
  !
  subroutine m3d_multiply_rr_g(m,v,mv)
    real(rk), intent(in)  :: m(:,:)  ! Tri-diagonal matrix
    real(rk), intent(in)  :: v(:,:)  ! Vectors
    real(rk), intent(out) :: mv(:,:) ! Vectors
    !
    include "tridiagonal_tools_m3d_multiply_g_common.f90"
  end subroutine m3d_multiply_rr_g
  !
  subroutine m3d_multiply_xrr_g(m,v,mv)
    real(xk), intent(in)  :: m(:,:)  ! Tri-diagonal matrix
    real(xk), intent(in)  :: v(:,:)  ! Vectors
    real(xk), intent(out) :: mv(:,:) ! Vectors
    !
    include "tridiagonal_tools_m3d_multiply_g_common.f90"
  end subroutine m3d_multiply_xrr_g
  !
  subroutine m3d_multiply_rc_g(m,v,mv)
    real(rk), intent(in)     :: m(:,:)  ! Tri-diagonal matrix
    complex(rk), intent(in)  :: v(:,:)  ! Vectors
    complex(rk), intent(out) :: mv(:,:) ! Vectors
    !
    include "tridiagonal_tools_m3d_multiply_g_common.f90"
  end subroutine m3d_multiply_rc_g
  !
  subroutine m3d_multiply_cc_g(m,v,mv)
    complex(rk), intent(in)  :: m(:,:)  ! Tri-diagonal matrix
    complex(rk), intent(in)  :: v(:,:)  ! Vectors
    complex(rk), intent(out) :: mv(:,:) ! Vectors
    !
    include "tridiagonal_tools_m3d_multiply_g_common.f90"
  end subroutine m3d_multiply_cc_g
  !
  subroutine m3d_decompose_r(m,mf,fail)
    real(rk), intent(in)           :: m (:,:) ! Tridiagonal matrix, assumed to be diagonally dominant
                                              ! As the result, we do not need to bother with pivoting
    real(rk), intent(out)          :: mf(:,:) ! Factorized tridiagonal matrix
    logical, intent(out), optional :: fail    ! Set to .true. if decomposition fails; 
                                              ! if fail is absent, abort on decomposition failure.
    !
    real(rk)    :: denom
    !
    include "tridiagonal_tools_m3d_decompose_common.f90"
  end subroutine m3d_decompose_r
  !
  subroutine m3d_decompose_xr(m,mf,fail)
    real(xk), intent(in)           :: m (:,:) ! Tridiagonal matrix, assumed to be diagonally dominant
                                              ! As the result, we do not need to bother with pivoting
    real(xk), intent(out)          :: mf(:,:) ! Factorized tridiagonal matrix
    logical, intent(out), optional :: fail    ! Set to .true. if decomposition fails; 
                                              ! if fail is absent, abort on decomposition failure.
    !
    real(xk)    :: denom
    !
    include "tridiagonal_tools_m3d_decompose_common.f90"
  end subroutine m3d_decompose_xr
  !
  subroutine m3d_decompose_c(m,mf,fail)
    complex(rk), intent(in)        :: m (:,:) ! Tridiagonal matrix, assumed to be diagonally dominant
                                              ! As the result, we do not need to bother with pivoting
    complex(rk), intent(out)       :: mf(:,:) ! Factorized tridiagonal matrix
    logical, intent(out), optional :: fail    ! Set to .true. if decomposition fails; 
                                              ! if fail is absent, abort on decomposition failure.
    !
    complex(rk) :: denom
    !
    include "tridiagonal_tools_m3d_decompose_common.f90"
  end subroutine m3d_decompose_c
  !
  subroutine m3d_solve_rr(mf,r,x)
    real(rk), intent(in)  :: mf(:,:) ! Factorized tridiagonal matrix
    real(rk), intent(in)  :: r (:)   ! Right-hand size
    real(rk), intent(out) :: x (:)   ! Solution vector
    !
    include "tridiagonal_tools_m3d_solve_common.f90"
  end subroutine m3d_solve_rr
  !
  subroutine m3d_solve_rc(mf,r,x)
    real(rk), intent(in)     :: mf(:,:) ! Factorized tridiagonal matrix
    complex(rk), intent(in)  :: r (:)   ! Right-hand size
    complex(rk), intent(out) :: x (:)   ! Solution vector
    !
    include "tridiagonal_tools_m3d_solve_common.f90"
  end subroutine m3d_solve_rc
  !
  subroutine m3d_solve_cc(mf,r,x)
    complex(rk), intent(in)  :: mf(:,:) ! Factorized tridiagonal matrix
    complex(rk), intent(in)  :: r (:)   ! Right-hand size
    complex(rk), intent(out) :: x (:)   ! Solution vector
    !
    include "tridiagonal_tools_m3d_solve_common.f90"
  end subroutine m3d_solve_cc
  !
  subroutine m3d_solve_rr_g(mf,r,x)
    real(rk), intent(in)  :: mf(:,:) ! Factorized tridiagonal matrix
    real(rk), intent(in)  :: r (:,:) ! Right-hand sizes
    real(rk), intent(out) :: x (:,:) ! Solution vectors
    !
    include "tridiagonal_tools_m3d_solve_g_common.f90"
  end subroutine m3d_solve_rr_g
  !
  subroutine m3d_solve_xrr_g(mf,r,x)
    real(xk), intent(in)  :: mf(:,:) ! Factorized tridiagonal matrix
    real(xk), intent(in)  :: r (:,:) ! Right-hand sizes
    real(xk), intent(out) :: x (:,:) ! Solution vectors
    !
    include "tridiagonal_tools_m3d_solve_g_common.f90"
  end subroutine m3d_solve_xrr_g
  !
  subroutine m3d_solve_rc_g(mf,r,x)
    real(rk), intent(in)     :: mf(:,:) ! Factorized tridiagonal matrix
    complex(rk), intent(in)  :: r (:,:) ! Right-hand sizes
    complex(rk), intent(out) :: x (:,:) ! Solution vectors
    !
    include "tridiagonal_tools_m3d_solve_g_common.f90"
  end subroutine m3d_solve_rc_g
  !
  subroutine m3d_solve_cc_g(mf,r,x)
    complex(rk), intent(in)  :: mf(:,:) ! Factorized tridiagonal matrix
    complex(rk), intent(in)  :: r (:,:) ! Right-hand sizes
    complex(rk), intent(out) :: x (:,:) ! Solution vectors
    !
    include "tridiagonal_tools_m3d_solve_g_common.f90"
  end subroutine m3d_solve_cc_g
  !
  subroutine m3d_left_scale_cr(s,m,sm)
    complex(rk), intent(in)  :: s   (:) ! Diagonal matrix
    real(rk), intent(in)     :: m (:,:) ! Tridiagonal matrix
    complex(rk), intent(out) :: sm(:,:) ! Tridiagonal s . m
    !
    include "tridiagonal_tools_m3d_left_scale_common.f90"
  end subroutine m3d_left_scale_cr
  !
  subroutine m3d_right_scale_rc(m,s,ms)
    real(rk), intent(in)     :: m (:,:) ! Tridiagonal matrix
    complex(rk), intent(in)  :: s   (:) ! Diagonal matrix
    complex(rk), intent(out) :: ms(:,:) ! Tridiagonal m . s
    !
    include "tridiagonal_tools_m3d_right_scale_common.f90"
  end subroutine m3d_right_scale_rc
  !
  subroutine m3d_transpose_r(m,mt)
    real(rk), intent(in)  :: m (:,:) ! Tridiagonal matrix
    real(rk), intent(out) :: mt(:,:) ! Transpose of m
    !
    if (any(ubound(m)/=ubound(mt))) then
      stop 'tridiagonal_tools%m3d_transpose_r - bad input sizes'
    end if
    !
    mt(1,:) = m(1,:)
    mt(2,:) = m(3,:)
    mt(3,:) = m(2,:)
  end subroutine m3d_transpose_r
end module tridiagonal_tools
