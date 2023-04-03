!
!  Data structures for the H.G.M.-style spherical TDSE solver
!  Note that the actual operators are sometimes different from HGM's
!  This data mayb be shared by a number of modules/routines
!
module spherical_data
  use accuracy
  use constants
  use timer
  use tridiagonal_tools
  use potential_tools
  use cap_tools
  implicit none
  private
  public sd_nradial, sd_nspin, sd_lmax, sd_mmin, sd_mmax
  public sd_rgrid, sd_rgrid_dr, sd_rgrid_r0, sd_rgrid_scale, sd_rgrid_file
  public sd_rgrid_zeta, sd_rgrid_npad
  public sd_rgrid_grad_delta
  public sd_rtab, sd_drtab, sd_pottab, sd_captab, sd_dvdrtab
  public sd_pot_nonlocal
  public sd_capped, sd_cap_start
  public sd_d2n_l0, sd_m2n_l0, sd_m2nf_l0, sd_d2n_lx, sd_m2n_lx, sd_m2nf_lx
  public sd_d1n_l0, sd_m1n_l0, sd_m1nf_l0, sd_d1n_lx, sd_m1n_lx, sd_m1nf_lx
  public sd_d2t_l0, sd_m2t_l0, sd_m2tf_l0, sd_d2t_lx, sd_m2t_lx, sd_m2tf_lx
  public sd_d1t_l0, sd_m1t_l0, sd_m1tf_l0, sd_d1t_lx, sd_m1t_lx, sd_m1tf_lx
  public sd_wfn
  public sd_initialize
  public sd_expand_implicit_operator, sd_laser_clm
  !
  interface sd_expand_implicit_operator
    module procedure sd_expand_implicit_operator_r
    module procedure sd_expand_implicit_operator_c
  end interface sd_expand_implicit_operator
  !
  integer, parameter          :: iu_temp        = 23         ! An arbitrary unit number, which can be used here
  !
  !  Common data, shared by all wavefunctions which are handled by our propagator.
  !
  !  Parameters directly controlled by the user
  ! 
  integer(ik), save           :: sd_nradial     = 220_ik     ! Number of radial points in the grid
  integer(ik), save           :: sd_nspin       = 1_ik       ! Number of spin components in the wavefunctions; either 1 or 2
  integer(ik), save           :: sd_lmax        = 10_ik      ! Largest angular momentum in the grid (lowest is always 0)
  integer(ik), save           :: sd_mmin        =  0_ik      ! Smallest angular momentum projection 
  integer(ik), save           :: sd_mmax        =  0_ik      ! Largest angular momentum projection 
                                                             ! There are only two allowed usage cases:
                                                             ! a) sd_mmin = sd_mmax = constant (linear polarization)
                                                             ! b) sd_mmin = -sd_lmax; sd_mmax = sd_lmax (arbitrary field)
  character(len=clen), save   :: sd_rgrid       = 'uniform'  ! Choice of the radial grid. Can be one of the following:
                                                             ! 'uniform'     - Uniformly distributed points; uses sd_rgrid_dr
                                                             ! 'log'         - Logarithmically distributed points; 
                                                             !                 uses sd_rgrid_r0 and sd_rgrid_scale
                                                             ! 'log-uniform' - Logarithmic grid until step reaches sd_rgrid_dr;
                                                             !                 uniform grid afterwards
                                                             ! 'read'        - Arbitrary user-defined points; uses rd_rgrid_file
                                                             ! Note that both 'log' and 'log-uniform' add a small uniform segment
                                                             ! between the origin and the first 'log' point; otherwise we get
                                                             ! highly oscillatory terms in the derivative operators.
                                                             ! WARNING: Not all grids produced by sd_rgrid/='uniform' are
                                                             ! WARNING: of a good quality. Watch out for rapid oscillations
                                                             ! WARNING: in the left wavefunction: this can (and does!) cause
                                                             ! WARNING: problems with the laser coupling operator
  real(rk), save              :: sd_rgrid_zeta  = 1.0_rk     ! Effective nuclear charge at the origin; this value is only
                                                             ! used to establish the boundary condition on the wavefunction
  real(rk), save              :: sd_rgrid_dr    = 0.2_rk     ! Radial grid spacing (sd_rgrid=='uniform')
  real(rk), save              :: sd_rgrid_r0    = 1e-2_rk    ! First point of the logarithmic grid (sd_rgrid=='log')
  real(rk), save              :: sd_rgrid_scale = 1.05_rk    ! Scaling factor of the logarithmic grid (sd_rgrid=='log')
                                                             ! r(i) = sd_rgrid_r0 * sd_rgrid_scale**(i-1)
  integer(ik), save           :: sd_rgrid_npad  = -1_ik      ! Number of extra "padding" uniformly-spaced grid points at the
                                                             ! origin; this only applies to 'log' and 'log-uniform' grids
                                                             ! (-1) is a special value; it will chose number of uniform padding
                                                             ! points such that grid spacing increases monotonically at the inner
                                                             ! boundary.
  character(len=clen), save   :: sd_rgrid_file  = 'grid.dat' ! File containing a list of radial grid point coordinates (in Bohr), 
                                                             ! one point per line. The initial zero point is assumed implicitly,
                                                             ! and should not be included.
                                                             ! (sd_nradial+1) points are expected (see sd_rtab(:))
                                                             ! Coordinates must be in the ascending order (sd_rgrid=='read')
  real(rk), save              :: sd_rgrid_grad_delta &
                                                = 0.01_rk    ! Step-size control for numerical differential of the potential.
                                                             ! The actual step size for point (i) is computed as:
                                                             ! rd_drtab(i)*sd_rgrid_grad_delta
  !
  !  Derived parameters, computed from the values above
  !
  logical, save               :: sd_globinit    = .false.    ! Global initialialization is complete
  real(rk), allocatable, save :: sd_rtab (:)                 ! Radial coordinates, (sd_nradial+1) points, in Bohr.
                                                             ! The grid starts at zero (implicitly) and ends at sd_rtab(sd_nradiali+1). 
  real(rk), allocatable, save :: sd_drtab(:)                 ! Distance between the current and preceeding grid points, in Bohr.
                                                             ! (1:sd_nradial+1)
                                                             ! sd_drtab(i) = sd_rtab(i)-sd_drtab(i-1), where sd_rtab(0) is taken as zero
  real(rk), allocatable, save :: sd_pottab(:,:,:)            ! Per-channel multiplicative potential at grid points
                                                             ! 1: (1:sd_nradial) - Radial grid indices
                                                             ! 2: (1)            - Spin-orbit sub-channel; currently unused
                                                             ! 3: (0:sd_lmax)    - Angular momentum of the channel
  logical, save               :: sd_capped                   ! .True. if sd_captab(:) is not identically zero
  integer(ik), save           :: sd_cap_start                ! First point of the capping potential
  complex(rk), allocatable, save \
                              :: sd_captab(:)                ! Complex absorbing potential at the edge of the grid
  real(rk), allocatable, save :: sd_dvdrtab(:)               ! Radial gradient of the potential, (sd_nradial) points, in
                                                             ! Hartree per Bohr.
  logical, save               :: sd_pot_nonlocal = .false.   ! Will be set to .true. if the potential appears to be non-local
  !
  !  Differential operators, radial coordinate.
  !  All operators here are in the implicit form, namely:
  !   Op(X) = M^{-1} . D . X
  !  where M and D are three-diagonal real matrices (not necessarily symmetric!).
  !  R=0 boundary conditions are different in L=0 and all the remaining channels,
  !  so that we have to use different implicit operators.
  !
  !  The tri-diagonal matrices are stored as follows:
  !    v(1,:) - the diagonal
  !    v(2,:) - the sub-diagonal; the last element is not used
  !    v(3,:) - the super-diagonal; the last element is not used
  !
  real(rk), allocatable, save :: sd_d2n_l0 (:,:)              ! Laplacian, L=0 channel
  real(rk), allocatable, save :: sd_m2n_l0 (:,:)              ! 
  real(rk), allocatable, save :: sd_m2nf_l0(:,:)              ! sd_m2_l0, factored by m3d_decompose,
  real(rk), allocatable, save :: sd_d2n_lx (:,:)              ! Laplacian, L>0 channel
  real(rk), allocatable, save :: sd_m2n_lx (:,:)              ! 
  real(rk), allocatable, save :: sd_m2nf_lx(:,:)              ! sd_m2_lx, factored by m3d_decompose,
  real(rk), allocatable, save :: sd_d1n_l0 (:,:)              ! Radial gradient, L=0 channel
  real(rk), allocatable, save :: sd_m1n_l0 (:,:)              ! 
  real(rk), allocatable, save :: sd_m1nf_l0(:,:)              ! sd_m1_l0, factored by m3d_decompose,
  real(rk), allocatable, save :: sd_d1n_lx (:,:)              ! Radial gradient, L>0 channel
  real(rk), allocatable, save :: sd_m1n_lx (:,:)              ! 
  real(rk), allocatable, save :: sd_m1nf_lx(:,:)              ! sd_m1_lx, factored by m3d_decompose,
  !
  !  Transposes
  !
  real(rk), allocatable, save :: sd_d2t_l0 (:,:)              ! Laplacian, L=0 channel
  real(rk), allocatable, save :: sd_m2t_l0 (:,:)              ! 
  real(rk), allocatable, save :: sd_m2tf_l0(:,:)              ! sd_m2_l0, factored by m3d_decompose,
  real(rk), allocatable, save :: sd_d2t_lx (:,:)              ! Laplacian, L>0 channel
  real(rk), allocatable, save :: sd_m2t_lx (:,:)              ! 
  real(rk), allocatable, save :: sd_m2tf_lx(:,:)              ! sd_m2_lx, factored by m3d_decompose,
  real(rk), allocatable, save :: sd_d1t_l0 (:,:)              ! Radial gradient, L=0 channel
  real(rk), allocatable, save :: sd_m1t_l0 (:,:)              ! 
  real(rk), allocatable, save :: sd_m1tf_l0(:,:)              ! sd_m1_l0, factored by m3d_decompose,
  real(rk), allocatable, save :: sd_d1t_lx (:,:)              ! Radial gradient, L>0 channel
  real(rk), allocatable, save :: sd_m1t_lx (:,:)              ! 
  real(rk), allocatable, save :: sd_m1tf_lx(:,:)              ! sd_m1_lx, factored by m3d_decompose,
  !
  !  1-electron wavefunction; there may be more than one floating around;
  !
  type sd_wfn          
    complex(rk), allocatable :: wfn(:,:,:,:)     ! The indices are:
                                                 ! 1 - [1:sd_nradial]    - Radial grid points
                                                 ! 2 - [1:sd_nspin]      - Spin components
                                                 ! 3 - [0:sd_lmax]       - Angular momentum components
                                                 !                         Note that only the range [abs(m):sd_lmax] is used
                                                 ! 4 - [sd_mmin:sd_mmax] - Projection of the angular momentum on the local Z axis
  end type sd_wfn
  !
  contains
  !
  subroutine sd_initialize
    if (sd_globinit) return
    call TimerStart('Grid initialization')
    call allocate_arrays
    call initialize_radial_grid
    call initialize_radial_gradient
    call initialize_radial_laplacian
    call initialize_channel_potentials
    sd_globinit = .true.
    call TimerStop('Grid initialization')
  end subroutine sd_initialize 
  !
  !  Construct explict operator matrix for an operator given in implicit form:
  !
  !  Conceptually: block = mop . op
  !
  subroutine sd_expand_implicit_operator_r(op,mop,block)
    real(rk), intent(in)  :: op   (:,:) ! Tridiagonal operator matrix
    real(rk), intent(in)  :: mop  (:,:) ! Modifier operator matrix, given as a factored tridiagonal inverse
    real(rk), intent(out) :: block(:,:) ! Explicit operator matrix
    !
    integer(ik) :: sz, icol
    real(rk)    :: rhs(size(block,dim=1))
    real(rk)    :: tmp(size(block,dim=1))
    !
    sz = size(block,dim=1)
    !
    !  Sanity check
    !
    if (size(op,dim=1)<3 .or. size(mop,dim=1)<3 .or. &
        size(op,dim=2)/=sz .or. size(mop,dim=2)/=sz .or. size(block,dim=2)/=sz) then
      stop 'spherical_data%sd_expand_implicit_operator_r - dimensions mismatch'
    end if
    !
    rhs(:) = 0
    scan_columns: do icol=1,sz
      rhs(icol) = 1._rk
      call m3d_multiply(op,rhs,tmp)
      call m3d_solve(mop,tmp,block(:,icol))
      rhs(icol) = 0._rk
    end do scan_columns
  end subroutine sd_expand_implicit_operator_r
  !
  subroutine sd_expand_implicit_operator_c(op,mop,block)
    complex(rk), intent(in)  :: op   (:,:) ! Tridiagonal operator matrix
    complex(rk), intent(in)  :: mop  (:,:) ! Modifier operator matrix, given as a factored tridiagonal inverse
    complex(rk), intent(out) :: block(:,:) ! Explicit operator matrix
    !
    integer(ik) :: sz, icol
    complex(rk) :: rhs(size(block,dim=1))
    complex(rk) :: tmp(size(block,dim=1))
    !
    sz = size(block,dim=1)
    !
    !  Sanity check
    !
    if (size(op,dim=1)<3 .or. size(mop,dim=1)<3 .or. &
        size(op,dim=2)/=sz .or. size(mop,dim=2)/=sz .or. size(block,dim=2)/=sz) then
      stop 'spherical_data%sd_expand_implicit_operator_c - dimensions mismatch'
    end if
    !
    rhs(:) = 0
    scan_columns: do icol=1,sz
      rhs(icol) = 1._rk
      call m3d_multiply(op,rhs,tmp)
      call m3d_solve(mop,tmp,block(:,icol))
      rhs(icol) = 0._rk
    end do scan_columns
  end subroutine sd_expand_implicit_operator_c
  !
  !  Angular coupling coefficients; see propagator_tools.f90 and wavefunction_tools.f90
  !
  function sd_laser_clm(l,m) result(clm)
    integer(ik), intent(in) :: l, m ! Quantum numbers
    real(xk)                :: clm  ! Coupling coefficient
    !
    clm = sqrt(real(l**2 - m**2,kind=xk)/real(4*l**2 - 1,kind=xk))
  end function sd_laser_clm
  !
  !  Internal routines below
  !
  subroutine allocate_arrays
    integer(ik) :: alloc
    !
    allocate (sd_rtab(sd_nradial+1),sd_drtab(sd_nradial+1), &
              sd_pottab(sd_nradial,1,0:sd_lmax), sd_captab(sd_nradial), &
              sd_dvdrtab(sd_nradial), &
              sd_d2n_l0(3,sd_nradial),sd_m2n_l0(3,sd_nradial), sd_m2nf_l0(3,sd_nradial), &
              sd_d2n_lx(3,sd_nradial),sd_m2n_lx(3,sd_nradial), sd_m2nf_lx(3,sd_nradial), &
              sd_d1n_l0(3,sd_nradial),sd_m1n_l0(3,sd_nradial), sd_m1nf_l0(3,sd_nradial), &
              sd_d1n_lx(3,sd_nradial),sd_m1n_lx(3,sd_nradial), sd_m1nf_lx(3,sd_nradial), &
              sd_d2t_l0(3,sd_nradial),sd_m2t_l0(3,sd_nradial), sd_m2tf_l0(3,sd_nradial), &
              sd_d2t_lx(3,sd_nradial),sd_m2t_lx(3,sd_nradial), sd_m2tf_lx(3,sd_nradial), &
              sd_d1t_l0(3,sd_nradial),sd_m1t_l0(3,sd_nradial), sd_m1tf_l0(3,sd_nradial), &
              sd_d1t_lx(3,sd_nradial),sd_m1t_lx(3,sd_nradial), sd_m1tf_lx(3,sd_nradial), &
              stat=alloc)
    if (alloc/=0) then
      write (out,"('spherical_data%allocate_array: memory allocation failed with code ',i0)") alloc
      stop 'spherical_data%allocate_array - memory allocation failed'
    end if
  end subroutine allocate_arrays
  !
  subroutine initialize_radial_grid
    integer(ik) :: ipt, npt_log
    integer     :: ios
    !
    write (out,"()")
    select case (sd_rgrid)
      case default
        write (out,"('spherical_data%initialize_radial_grid: Radial grid ',a,' is not recognized')") trim(sd_rgrid)
        stop 'spherical_data%initialize_radial_grid - bad grid'
      case ('uniform')
        write (out,"('Creating uniform grid with step ',g24.13,' Bohr')") sd_rgrid_dr
        fill_uniform_grid: do ipt=1,sd_nradial+1
          sd_rtab(ipt) = ipt*sd_rgrid_dr
        end do fill_uniform_grid
      case ('log')
        if (sd_rgrid_npad<0) sd_rgrid_npad = max(0,nint(1._rk/(sd_rgrid_scale-1._rk)))
        write (out,"('Creating logarithmic grid starting at ',g24.13,' Bohr; scale factor ',g24.13)") sd_rgrid_r0, sd_rgrid_scale
        write (out,"('Inserting ',i0,'-point uniform padding grid at the origin')") sd_rgrid_npad
        pad_log_grid: do ipt=1,sd_rgrid_npad
          sd_rtab(ipt) = (ipt*sd_rgrid_r0)/(sd_rgrid_npad+1)
        end do pad_log_grid
        fill_log_grid: do ipt=sd_rgrid_npad+1,sd_nradial+1
          sd_rtab(ipt) = sd_rgrid_r0 * sd_rgrid_scale**(ipt-1-sd_rgrid_npad)
        end do fill_log_grid
      case ('log-uniform')
        if (sd_rgrid_npad<0) sd_rgrid_npad = max(0,nint(1._rk/(sd_rgrid_scale-1._rk)))
        npt_log = 1 + nint(log(sd_rgrid_dr/(sd_rgrid_r0*(sd_rgrid_scale-1)))/log(sd_rgrid_scale))
        write (out,"('Creating log-uniform grid')")
        write (out,"('Logarithmic grid starting at ',g24.13,' Bohr; scale factor ',g24.13)") sd_rgrid_r0, sd_rgrid_scale
        write (out,"('Inserting ',i0,'-point uniform padding grid at the origin')") sd_rgrid_npad
        write (out,"('Uniform grid step ',g24.13,' Bohr')") sd_rgrid_dr
        write (out,"('Logarithmic to uniform switch after point ',i0)") npt_log
        !
        if (npt_log<1 .or. npt_log>=sd_nradial-sd_rgrid_npad) then
          write (out,"('Logarithmic to uniform switching point ',i0,' is not between 1 and ',i0)") npt_log, sd_nradial-sd_rgrid_npad
          stop 'spherical_data%initialize_radial_grid - bad log-uniform parameters'
        end if
        !
        pad_lu_grid: do ipt=1,sd_rgrid_npad
          sd_rtab(ipt) = (ipt*sd_rgrid_r0)/(sd_rgrid_npad+1)
        end do pad_lu_grid
        fill_lu_grid_log_part: do ipt=1+sd_rgrid_npad,npt_log+sd_rgrid_npad
          sd_rtab(ipt) = sd_rgrid_r0 * sd_rgrid_scale**(ipt-1-sd_rgrid_npad)
        end do fill_lu_grid_log_part
        fill_lu_grid_uniform_part: do ipt=npt_log+1+sd_rgrid_npad,sd_nradial+1
          sd_rtab(ipt) = sd_rtab(ipt-1) + sd_rgrid_dr
        end do fill_lu_grid_uniform_part
        !
        write (out,"('Logarithmic grid: ',i6,' points between ',g24.13,' and ',g24.13,' Bohr')") &
               npt_log, sd_rtab(1), sd_rtab(npt_log)
        write (out,"('    Uniform grid: ',i6,' points between ',g24.13,' and ',g24.13,' Bohr')") &
               sd_nradial-npt_log, sd_rtab(npt_log+1), sd_rtab(sd_nradial)
      case ('read')
        write (out,"('Reading user-defined from ',a)") trim(sd_rgrid_file)
        open (iu_temp,form='formatted',action='read',position='rewind',status='old',file=trim(sd_rgrid_file))
        fill_user_grid: do ipt=1,sd_nradial+1
          read (iu_temp,*,iostat=ios) sd_rtab(ipt)
          if (ios/=0) then
            write (out,"('Error ',i0,' reading point ',i0,' of the radial grid.')") ios, ipt
            stop 'spherical_data%initialize_radial_grid - read error'
          end if
        end do fill_user_grid
        close (iu_temp)
    end select
    !
    !  Now calculate the distance table from the grid positions
    ! 
    sd_drtab(1) = sd_rtab(1)
    distance_table: do ipt=2,sd_nradial+1
      if (sd_rtab(ipt)<=sd_rtab(ipt-1)) then
        write (out,"('Radial grid point ',i0,' at r= ',g24.13,' is before the preceeding point at r= ',g24.13,'?!')") &
               ipt, sd_rtab(ipt), sd_rtab(ipt-1)
        stop 'spherical_data%initialize_radial_grid - bad grid'
      end if
      sd_drtab(ipt) = sd_rtab(ipt)-sd_rtab(ipt-1)
    end do distance_table
    !
    write (out,"(/'        Number of radial grid points = ',i0)") sd_nradial
    write (out,"( '    First explicit grid point is at r= ',g24.13)") sd_rtab(1)
    write (out,"( '      Last explicit grid point s at r= ',g24.13)") sd_rtab(sd_nradial)
    write (out,"( '                  The world ends at r= ',g24.13)") sd_rtab(sd_nradial+1)
    write (out,"( 'Effective nuclear charge for L=0 grid= ',g24.13)") sd_rgrid_zeta
    write (out,"()")
  end subroutine initialize_radial_grid
  !
  !  For the derivation of the implicit derivative expressions, see notes in:
  !
  !    hgm-operators-radial-derivatives-2014_Oct_31.pdf
  !    hgm-operators-involving-derivatives-2014_Nov_04.pdf
  !
  subroutine initialize_radial_gradient
    integer(ik) :: ipt
    real(rk)    :: cfac
    !
    !  Interior points; general expression.
    !  For sanity, all expressions are computer-generated.
    !
    grad_lx: do ipt=1,sd_nradial
      ! Diagonal
      sd_d1n_lx(1,ipt) = 1/sd_drtab(ipt) - 1/sd_drtab(ipt+1) + (-sd_drtab(ipt) + sd_drtab(ipt+1))/ &
                         (sd_drtab(ipt)**2 + sd_drtab(ipt)*sd_drtab(ipt+1) + sd_drtab(ipt+1)**2)
      sd_m1n_lx(1,ipt) = (sd_drtab(ipt) + sd_drtab(ipt+1))**2/ &
                         (2*(sd_drtab(ipt)**2 + sd_drtab(ipt)*sd_drtab(ipt+1) + sd_drtab(ipt+1)**2))
      if (ipt>=sd_nradial) cycle grad_lx
      ! Sub-diagonal
      sd_d1n_lx(2,ipt) = -((sd_drtab(ipt+2)**2*(2*sd_drtab(ipt+1) + sd_drtab(ipt+2)))/ &
                         (sd_drtab(ipt+1)*(sd_drtab(ipt+1) + sd_drtab(ipt+2))* &
                         (sd_drtab(ipt+1)**2 + sd_drtab(ipt+1)*sd_drtab(ipt+2) + sd_drtab(ipt+2)**2)))
      sd_m1n_lx(2,ipt) = sd_drtab(ipt+2)**2/(2*(sd_drtab(ipt+1)**2 + sd_drtab(ipt+1)*sd_drtab(ipt+2) + sd_drtab(ipt+2)**2))
      ! Super-diagonal
      sd_d1n_lx(3,ipt) = (sd_drtab(ipt)**2*(sd_drtab(ipt) + 2*sd_drtab(ipt+1)))/(sd_drtab(ipt+1)* &
                         (sd_drtab(ipt) + sd_drtab(ipt+1))*(sd_drtab(ipt)**2 + &
                         sd_drtab(ipt)*sd_drtab(ipt+1) + sd_drtab(ipt+1)**2))
      sd_m1n_lx(3,ipt) = sd_drtab(ipt)**2/(2*(sd_drtab(ipt)**2 + sd_drtab(ipt)*sd_drtab(ipt+1) + sd_drtab(ipt+1)**2))
    end do grad_lx
    sd_d1n_l0 = sd_d1n_lx
    sd_m1n_l0 = sd_m1n_lx
    !
    !  Nothing to do for L>0; general case already covers it.
    !
    !
    !  Special case: L=0.
    !
    sd_d1n_l0(1,1) = -(((sd_drtab(1) + sd_drtab(2))**2*(sd_drtab(1) + 2*sd_drtab(2))* &
                     (-6*sd_drtab(1)**2 + sd_drtab(2)**2 + 2*sd_drtab(1)**3*sd_rgrid_zeta + &
                     sd_drtab(1)*sd_drtab(2)*(3 - 2*sd_drtab(2)*sd_rgrid_zeta)))/ &
                     (sd_drtab(1)*sd_drtab(2)*(sd_drtab(1)**2 + sd_drtab(1)*sd_drtab(2) + sd_drtab(2)**2)* &
                     (-6*sd_drtab(1)**2 - 15*sd_drtab(1)*sd_drtab(2) - 8*sd_drtab(2)**2 + &
                     2*sd_drtab(1)*(sd_drtab(1) + sd_drtab(2))*(sd_drtab(1) + 2*sd_drtab(2))*sd_rgrid_zeta)))
    sd_m1n_l0(1,1) = ((sd_drtab(1) + sd_drtab(2))**2*(sd_drtab(1) + 2*sd_drtab(2))* (-3*sd_drtab(1) - sd_drtab(2) + &
                     sd_drtab(1)*(sd_drtab(1) + sd_drtab(2))*sd_rgrid_zeta))/((sd_drtab(1)**2 + sd_drtab(1)*sd_drtab(2) &
                     + sd_drtab(2)**2)* (-6*sd_drtab(1)**2 - 15*sd_drtab(1)*sd_drtab(2) - 8*sd_drtab(2)**2 + &
                     2*sd_drtab(1)*(sd_drtab(1) + sd_drtab(2))*(sd_drtab(1) + 2*sd_drtab(2))*sd_rgrid_zeta))
    sd_m1n_l0(3,1) = (sd_drtab(1)**2*(sd_drtab(1) + 2*sd_drtab(2))*(-3*sd_drtab(1) - 2*sd_drtab(2) + sd_drtab(1)* &
                     (sd_drtab(1) + sd_drtab(2))*sd_rgrid_zeta))/((sd_drtab(1)**2 + sd_drtab(1)*sd_drtab(2) + &
                     sd_drtab(2)**2)* (-6*sd_drtab(1)**2 - 15*sd_drtab(1)*sd_drtab(2) - 8*sd_drtab(2)**2 + 2*sd_drtab(1)* &
                     (sd_drtab(1) + sd_drtab(2))*(sd_drtab(1) + 2*sd_drtab(2))*sd_rgrid_zeta))
    !
    write (out,"('L>=1 channel d1(1,1)= ',g24.13,' m1(1,1)= ',g24.13)") sd_d1n_lx(1,1), sd_m1n_lx(1,1)
    write (out,"(' L=0 channel d1(1,1)= ',g24.13,' m1(1,1)= ',g24.13)") sd_d1n_l0(1,1), sd_m1n_l0(1,1)
    !
    !  Factor m1 matrices; strictly speaking this is not required for propagation,
    !  since we never take the first radial derivative by itself. However, it is
    !  used in the gradient accuracy check.
    !
    call m3d_decompose(sd_m1n_l0,sd_m1nf_l0)
    call m3d_decompose(sd_m1n_lx,sd_m1nf_lx)
    !
    !  We also need transposes and matching inverses
    !
    call m3d_transpose(sd_d1n_l0,sd_d1t_l0)
    call m3d_transpose(sd_d1n_lx,sd_d1t_lx)
    call m3d_transpose(sd_m1n_l0,sd_m1t_l0)
    call m3d_transpose(sd_m1n_lx,sd_m1t_lx)
    !
    call m3d_decompose(sd_m1t_l0,sd_m1tf_l0)
    call m3d_decompose(sd_m1t_lx,sd_m1tf_lx)
  end subroutine initialize_radial_gradient
  !
  subroutine initialize_radial_laplacian
    integer(ik) :: ipt
    !
    !  Begin with the interior points
    !
    lap_lx: do ipt=1,sd_nradial
      ! Diagonal
      sd_d2n_lx(1,ipt) = -2._rk/(sd_drtab(ipt)*sd_drtab(ipt+1))
      sd_m2n_lx(1,ipt) = (1._rk/6._rk)*(3._rk + sd_drtab(ipt)/sd_drtab(ipt+1) + sd_drtab(ipt+1)/sd_drtab(ipt))
      if (ipt>=sd_nradial) cycle lap_lx
      ! Sub-diagonal
      sd_d2n_lx(2,ipt) = 2._rk/(sd_drtab(ipt+1)*(sd_drtab(ipt+1)+sd_drtab(ipt+2)))
      sd_m2n_lx(2,ipt) = (1._rk/6._rk)*(2._rk - sd_drtab(ipt+2)/sd_drtab(ipt+1) &
                                             - sd_drtab(ipt+1)/(sd_drtab(ipt+1)+sd_drtab(ipt+2)))
      ! Super-diagonal
      sd_d2n_lx(3,ipt) = 2._rk/(sd_drtab(ipt+1)*(sd_drtab(ipt)+sd_drtab(ipt+1)))
      sd_m2n_lx(3,ipt) = (1._rk/6._rk)*(1._rk - sd_drtab(ipt)**2/((sd_drtab(ipt)+sd_drtab(ipt+1))*sd_drtab(ipt+1)))
    end do lap_lx
    sd_d2n_l0 = sd_d2n_lx
    sd_m2n_l0 = sd_m2n_lx
    !
    !  The first row is special; first L>0
    !  Expressions are computer-generated to preserve my sanity.
    !
    sd_d2n_lx(1,1) = (-2*(3*sd_drtab(1) - sd_drtab(2))*(sd_drtab(1) + sd_drtab(2))**2)/ &
                     (sd_drtab(1)**3*sd_drtab(2)*(3*sd_drtab(1) + 4*sd_drtab(2)))
    sd_m2n_lx(1,1) = ((sd_drtab(1) + sd_drtab(2))**2*(sd_drtab(1) + 3*sd_drtab(2)))/ &
                     (3*sd_drtab(1)*sd_drtab(2)*(3*sd_drtab(1) + 4*sd_drtab(2)))
    sd_m2n_lx(3,1) = -((sd_drtab(1) - 2*sd_drtab(2))*(sd_drtab(1) + sd_drtab(2)))/ &
                     (3*sd_drtab(2)*(3*sd_drtab(1) + 4*sd_drtab(2)))
    !
    !  Now L=0
    !
    sd_d2n_l0(1,1) = (2*(sd_drtab(1) + sd_drtab(2))*(6*sd_drtab(1) - (3*sd_drtab(1) - sd_drtab(2))* &
                     (sd_drtab(1) + sd_drtab(2))*sd_rgrid_zeta))/ (sd_drtab(1)**2*sd_drtab(2)* &
                     (-6*(sd_drtab(1) + sd_drtab(2)) + sd_drtab(1)*(3*sd_drtab(1) + 4*sd_drtab(2))*sd_rgrid_zeta))
    sd_m2n_l0(1,1) = ((sd_drtab(1) + sd_drtab(2))*(-3*(sd_drtab(1)**2 + 3*sd_drtab(1)*sd_drtab(2) + &
                     sd_drtab(2)**2) + sd_drtab(1)*(sd_drtab(1) + sd_drtab(2))*(sd_drtab(1) + 3*sd_drtab(2))* &
                     sd_rgrid_zeta))/  (3*sd_drtab(1)*sd_drtab(2)*(-6*(sd_drtab(1) + sd_drtab(2)) +  &
                     sd_drtab(1)*(3*sd_drtab(1) + 4*sd_drtab(2))*sd_rgrid_zeta))
    sd_m2n_l0(3,1) = (-3*sd_drtab(2)**2 - sd_drtab(1)**3*sd_rgrid_zeta + sd_drtab(1)**2* &
                     (3 + sd_drtab(2)*sd_rgrid_zeta) + sd_drtab(1)*sd_drtab(2)*(-3 + 2*sd_drtab(2)*sd_rgrid_zeta))/ &
                     (3*sd_drtab(2)*(-6*(sd_drtab(1) + sd_drtab(2)) + sd_drtab(1)*(3*sd_drtab(1) + 4*sd_drtab(2))*sd_rgrid_zeta))
    !
    write (out,"('L>=1 channel d2(1,1)= ',g24.13,' m2(1,1)= ',g24.13)") sd_d2n_lx(1,1), sd_m2n_lx(1,1)
    write (out,"(' L=0 channel d2(1,1)= ',g24.13,' m2(1,1)= ',g24.13)") sd_d2n_l0(1,1), sd_m2n_l0(1,1)
    !
    !  Factor m2 matrices
    !
    call m3d_decompose(sd_m2n_l0,sd_m2nf_l0)
    call m3d_decompose(sd_m2n_lx,sd_m2nf_lx)
    !
    !  We also need transposes and matching inverses
    !
    call m3d_transpose(sd_d2n_l0,sd_d2t_l0)
    call m3d_transpose(sd_d2n_lx,sd_d2t_lx)
    call m3d_transpose(sd_m2n_l0,sd_m2t_l0)
    call m3d_transpose(sd_m2n_lx,sd_m2t_lx)
    !
    call m3d_decompose(sd_m2t_l0,sd_m2tf_l0)
    call m3d_decompose(sd_m2t_lx,sd_m2tf_lx)
  end subroutine initialize_radial_laplacian
  !
  subroutine initialize_channel_potentials
    integer(ik) :: lval, ipt
    real(rk)    :: rdelta, rup, rdown, fup, fdown, fuplp, fdownlp
    !
    !  Potential itself, including the centrifugal term
    !
    scan_channels: do lval=0,sd_lmax
      scan_points: do ipt=1,sd_nradial
        sd_pottab(ipt,1,lval) = pt_evaluate_potential(lval,0,sd_rtab(ipt),centrifugal=.true.)
      end do scan_points
    end do scan_channels
    !
    !  Radial gradient of the potential (omitting the centrifugal term)
    !  We'll test for non-local potentials here: our dipole acceleration 
    !  expression is only valid for multiplicative potentials.
    !
    sd_pot_nonlocal = .false.
    gradient_points: do ipt=1,sd_nradial
      rdelta = sd_rgrid_grad_delta * sd_drtab(ipt)
      rup    = sd_rtab(ipt) + rdelta
      rdown  = sd_rtab(ipt) - rdelta
      fup    = pt_evaluate_potential(0,0,rup,  centrifugal=.false.)
      fdown  = pt_evaluate_potential(0,0,rdown,centrifugal=.false.)
      gradient_ltest: do lval=1,sd_lmax
        fuplp           = pt_evaluate_potential(lval,0,rup,  centrifugal=.false.)
        fdownlp         = pt_evaluate_potential(lval,0,rdown,centrifugal=.false.)
        sd_pot_nonlocal = sd_pot_nonlocal .or. (abs(fuplp  -fup  )>10._rk*spacing(abs(fup  ))) &
                                          .or. (abs(fdownlp-fdown)>10._rk*spacing(abs(fdown))) 
      end do gradient_ltest
      sd_dvdrtab(ipt) = (fup-fdown)/(rup-rdown)
    end do gradient_points
    if (sd_pot_nonlocal) then
      write (out,"(/'WARNING: Atomic potential appears to be L-dependent.'/)") 
    end if
    !
    initialize_cap: do ipt=1,sd_nradial
      sd_captab(ipt) = cap_evaluate_potential(sd_rtab(ipt),sd_rtab(sd_nradial+1))
    end do initialize_cap
    sd_capped    = .false.
    sd_cap_start = sd_nradial + 1
    scan_cap: do ipt=1,sd_nradial
      if (sd_captab(ipt)==(0._rk,0._rk)) cycle scan_cap
      sd_capped    = .true.
      sd_cap_start = ipt
      exit scan_cap
    end do scan_cap
    if (sd_capped) then
      write (out,"(/'Complex absorbing potential starts at grid point ',i0,' r= ',g24.13)") sd_cap_start, sd_rtab(sd_cap_start)
    end if
  end subroutine initialize_channel_potentials
end module spherical_data
