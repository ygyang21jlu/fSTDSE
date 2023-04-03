!
!  Spherical harmonics rotation.
!
module rotation_tools
  use accuracy
  use constants
  use math
  use spherical_data
  use tridiagonal_tools
  use timer
  implicit none
  private
  public rt_max_rot, rt_blocksize, rt_sense_real
  public rt_rotate
  public rt_rotate_bruteforce
  !
  real(rk), save    :: rt_max_rot   = 1e-3_rk  ! Maximum allowed rotation angle, in Radian
                                               ! The actual rotation threshold will be (rt_max_rot/(sd_lmax+1))

  integer(ik), save :: rt_blocksize = 16_ik    ! Max number of R points to process as a block. This is
                                               ! a performance optimization parameter; it does not change the results
  logical, save     :: rt_sense_real= .true.   ! If .true., try to detect real rotational propagator, and use
                                               ! the corresponding special case
  !
  contains
  !
  subroutine rt_prepare_rotation_steps(from,to,nstep,alp)
    real(xk), intent(in)                 :: from(:)   ! Initial (th,ph) in radian, in at least double precision
    real(xk), intent(in)                 :: to  (:)   ! Final   (th,ph)
    integer(ik), intent(out)             :: nstep     ! Number of subdivided steps
    real(xk), allocatable, intent(inout) :: alp (:,:) ! Final small-rotation angles
    !
    integer(ik) :: is, alloc
    real(xk)    :: delta(2)
    real(xk)    :: step_from(2), step_to(2), step_mid(2), step_delta(2)
    !
    if (size(from)/=2 .or. size(to)/=2) then
      stop 'rotation_tools%rt_prepare_rotation_steps - bad argument sizes'
    end if
    !
    if (allocated(alp)) deallocate(alp)
    !
    delta = to - from
    delta = delta - (2._xk*pi_xk*nint(delta/(2._xk*pi_xk)))
    if (all(abs(delta)<=spacing(1._xk))) then
      nstep = 0
      return
    end if
    !
    nstep = ceiling(sd_lmax*maxval(abs(delta),dim=1)/rt_max_rot)
    ! write (out,"('Need ',i0,' steps to achieve angular step of ',g24.13)") nstep, rt_max_rot/sd_lmax
    allocate (alp(3,nstep),stat=alloc)
    if (alloc/=0) then
      stop 'rotation_tools%rt_prepare_rotation_steps - allocation failed'
    end if
    break_steps: do is=1,nstep
      step_from  = from + (delta*(is-1))/nstep
      step_to    = from + (delta*(is-0))/nstep
      step_mid   = 0.5_xk * (step_from + step_to)
      step_delta = step_to - step_from
      !
      alp(1,is) = -step_delta(2) * sin(step_mid(1))
      alp(2,is) = -step_delta(1)
      alp(3,is) =  step_delta(2) * cos(step_mid(1))
    end do break_steps
  end subroutine rt_prepare_rotation_steps
  !
  !  rt_rotate can be applied to either left or right wavefunctions; the only
  !  difference is the sign of the time step and swapping of the upper/lower diagonals
  !  of the coupling matrix
  !
  subroutine rt_rotate(psi,from,to,left)
    type(sd_wfn), intent(inout) :: psi       ! Wavefunction we act upon
    real(xk), intent(in)        :: from(:)   ! Initial (th,ph) in radian
    real(xk), intent(in)        :: to  (:)   ! Final   (th,ph)
    logical, intent(in)         :: left      ! Set to .true. if rotating left wavefunction
    !
    !
    integer(ik)              :: lval, mval, nstep, istep, ispin
    integer(ik)              :: mmin, mmax            ! Smallesr and largest M for the current L
    integer(ik)              :: ir1, irn, irc
    integer(ik)              :: alloc
    complex(rk)              :: r_tmp (sd_mmin:sd_mmax,rt_blocksize)
    complex(rk)              :: r_tmp2(sd_mmin:sd_mmax,rt_blocksize)
    real(xk), allocatable    :: alp(:,:)              ! Table of rotation angles;
                                                      ! first index:  (ax,ay,az)
                                                      ! second index: rotation step
    complex(rk), allocatable :: prop_fwd(:,:,:)       ! Forward part of the propagator
    complex(rk), allocatable :: prop_rev(:,:,:)       ! Inverse part of the propagator
    real(rk), allocatable    :: prop_fwd_r(:,:,:)     ! Forward part of the propagator (special real case)
    real(rk), allocatable    :: prop_rev_r(:,:,:)     ! Inverse part of the propagator (special real case)
    complex(rk), allocatable :: srm(:,:)              ! Dense rotation matrix, from small-angle expansion
    logical                  :: use_dense             ! Use dense rotation branch; it's faster
    real(rk)                 :: eps
    logical                  :: use_real              ! Propagator is real; use faster special case
    !
    call TimerStart('Small-angle rotate')
    !
    call rt_prepare_rotation_steps(from,to,nstep,alp)
    if (nstep==0) then
      ! Early exit; nothing to do
      call TimerStop('Small-angle rotate')
      return
    end if
    !
    !$omp parallel default(none) &
    !$omp& shared(sd_mmin,sd_mmax,sd_lmax,sd_nradial,sd_nspin,nstep,alp) &
    !$omp& shared(rt_blocksize,psi,left,rt_sense_real) &
    !$omp& private(prop_fwd,prop_rev,prop_fwd_r,prop_rev_r,srm,alloc,lval,mmin,mmax) &
    !$omp& private(use_dense,irn,irc,r_tmp,r_tmp2,eps,use_real)
    !
    !  Parallelizing over L allows non-interleaved memory access from each
    !  thread. The alternative would be to parallelize over R, in which point
    !  we'll run a risk of false sharing and cache-line contention between
    !  the threads executing on separate CPUs.
    !
    !  We'll run in the -decreasing- order of L, since this keeps smaller
    !  work units at the end of the loop
    !
    allocate (prop_fwd  (3,sd_mmin:sd_mmax,nstep),prop_rev  (3,sd_mmin:sd_mmax,nstep), &
              prop_fwd_r(3,sd_mmin:sd_mmax,nstep),prop_rev_r(3,sd_mmin:sd_mmax,nstep), &
              srm(sd_mmin:sd_mmax,sd_mmin:sd_mmax),stat=alloc)
    if (alloc/=0) stop 'rotation_tools%rt_rotate - allocation failed'
    !
    ! Note that L=0 does not require rotation: it is an invariant
    !$omp do schedule(dynamic,1)
    process_l: do lval=sd_lmax,1,-1
      ! Range of allowable M values is from -L to L
      mmin = max(-lval,sd_mmin)
      mmax = min( lval,sd_mmax)
      call build_small_angle_propagators(lval,mmin,mmax,alp,prop_fwd(:,mmin:mmax,:),prop_rev(:,mmin:mmax,:))
      !
      ! Would it be faster to collapse series of sparse propagators into a dense matrix?
      !   "sparse" version costs 26*(mmax-mmin+1)*nstep FLOP
      !   "dense" version costs 5*(mmax-mmin+1)**3 FLOP; it is also more regular and vectorizes better.
      !   So, if nstep is much greater than (mmax-mmin+1)**2/5 or so, we should be switching to the dense code
      ! "much" is the measure of the relative efficiency of the code implementing dense matrix multiply and 
      ! tri-diagonal matrices; it is probably OK to assume dense linear algebra is considerably more efficient.
      !
      use_dense = nstep>(mmax-mmin+1)**2/200
      if (use_dense) then
        call expand_small_angle_propagators(prop_fwd(:,mmin:mmax,:),prop_rev(:,mmin:mmax,:),srm(mmin:mmax,mmin:mmax))
      end if
      use_real = .false.
      if (rt_sense_real) then
        eps = spacing(max(maxval(abs(prop_fwd)),maxval(abs(prop_rev))))
        use_real = all(abs(aimag(prop_fwd))<eps) .and. all(abs(aimag(prop_rev))<eps)
        if (use_real) then
          prop_fwd_r = real(prop_fwd,kind=rk)
          prop_rev_r = real(prop_rev,kind=rk)
        end if
      end if
      !
      !  Now comes the expensive part: the actual transformation, which must be done for each R point. 
      !
      process_radial: do ir1=1,sd_nradial,rt_blocksize
        irn = min(sd_nradial,ir1+rt_blocksize-1)
        irc = irn - ir1 + 1
        process_spin: do ispin=1,sd_nspin
          ! Fetch wavefunction to cache
          fetch_radial: do mval=mmin,mmax
            r_tmp(mval,:irc) = psi%wfn(ir1:irn,ispin,lval,mval)
          end do fetch_radial
          ! If we are dealing with the left wavefunction, we need to conjugate it here,
          ! then again after rotation is complete
          if (left) r_tmp(mmin:mmax,:irc) = conjg(r_tmp(mmin:mmax,:irc))
          ! Apply the propagators to each R point in turn
          if (use_dense) then
            r_tmp(mmin:mmax,1:irc) = matmul(srm(mmin:mmax,mmin:mmax),r_tmp(mmin:mmax,1:irc))
          else
            ! process_radial_cache: do irx=1,irc
            !   apply_rotations: do istep=1,nstep
            !     call m3d_multiply(prop_fwd(:,mmin:mmax,istep),r_tmp (mmin:mmax,irx),r_tmp2(mmin:mmax,irx))
            !     call m3d_solve   (prop_rev(:,mmin:mmax,istep),r_tmp2(mmin:mmax,irx),r_tmp (mmin:mmax,irx))
            !   end do apply_rotations
            ! end do process_radial_cache
            apply_rotations: do istep=1,nstep
              if (use_real) then
                call m3d_multiply(prop_fwd_r(:,mmin:mmax,istep),r_tmp (mmin:mmax,1:irc),r_tmp2(mmin:mmax,1:irc))
                call m3d_solve   (prop_rev_r(:,mmin:mmax,istep),r_tmp2(mmin:mmax,1:irc),r_tmp (mmin:mmax,1:irc))
              else
                call m3d_multiply(prop_fwd  (:,mmin:mmax,istep),r_tmp (mmin:mmax,1:irc),r_tmp2(mmin:mmax,1:irc))
                call m3d_solve   (prop_rev  (:,mmin:mmax,istep),r_tmp2(mmin:mmax,1:irc),r_tmp (mmin:mmax,1:irc))
              end if
            end do apply_rotations
          end if
          ! Undo the conjugation for the left wavefunction
          if (left) r_tmp(mmin:mmax,:irc) = conjg(r_tmp(mmin:mmax,:irc))
          ! Store transformed wavefunction back
          store_radial: do mval=mmin,mmax
            psi%wfn(ir1:irn,ispin,lval,mval) = r_tmp(mval,:irc)
          end do store_radial
        end do process_spin
      end do process_radial
    end do process_l
    !$omp end do
    deallocate (prop_fwd,prop_rev,prop_fwd_r,prop_rev_r,srm)
    !$omp end parallel
    deallocate (alp)
    call TimerStop('Small-angle rotate')
  end subroutine rt_rotate
  !
  subroutine build_small_angle_propagators(lval,mmin,mmax,alp,prop_fwd,prop_rev)
    integer(ik), intent(in)  :: lval            ! Angular momentum L
    integer(ik), intent(in)  :: mmin            ! Lower bound of the angular momentum projection
    integer(ik), intent(in)  :: mmax            ! Upper bound ..
    real(xk), intent(in)     :: alp(:,:)        ! Rotation steps
    complex(rk), intent(out) :: prop_fwd(:,:,:) ! Forward halfs of the rotation propagators
    complex(rk), intent(out) :: prop_rev(:,:,:) ! Reverse halfs of the rotation propagators
    !
    integer(ik) :: nstep, istep
    integer(ik) :: mval
    real(rk)    :: ax, ay, az
    real(rk)    :: d0    (mmin:mmax)   ! D0 rotation coefficients: coupling to same L
    real(rk)    :: dp    (mmin:mmax)   ! D+ rotation coefficients: coupling to higher L
    real(rk)    :: dm    (mmin:mmax)   ! D- rotation coefficients: coupling to lower L
    complex(rk) :: tmp (3,mmin:mmax)
    complex(rk) :: tmp2(3,mmin:mmax)
    !
    nstep = size(alp,dim=2)
    !
    !  First, form the L-dependent coupling coefficients
    !
    fill_d_tables: do mval=mmin,mmax
      d0(mval) = mval
      dp(mval) = 0.5_rk*sqrt(real(lval-mval,kind=rk)*real(lval+mval+1,kind=rk))
      dm(mval) = 0.5_rk*sqrt(real(lval+mval,kind=rk)*real(lval-mval+1,kind=rk))
    end do fill_d_tables
    !
    !  Now, form propagator matrices at each time step
    !
    build_propagators: do istep=1,nstep
      ax = -0.5_rk * real(alp(1,istep),kind=rk)
      ay = -0.5_rk * real(alp(2,istep),kind=rk)
      az = -0.5_rk * real(alp(3,istep),kind=rk)
      tmp(1,mmin:mmax  ) = d0(mmin  :mmax  ) * (-az          )
      tmp(2,mmin:mmax-1) = dp(mmin  :mmax-1) * (-ax +(0,1)*ay)
      tmp(3,mmin:mmax-1) = dm(mmin+1:mmax  ) * (-ax -(0,1)*ay)
      ! Construct direct and inverse propagators at this time step
      ! [1 - i H dt/2]
      prop_fwd(:,:,istep) = -(0,1)*tmp(:,:)
      prop_fwd(1,:,istep) = 1._rk + prop_fwd(1,:,istep)
      ! [1 + i H dt/2]^-1
      tmp2    (:,:)       = +(0,1)*tmp(:,:)
      tmp2    (1,:)       = 1._rk + tmp2(1,:)
      call m3d_decompose(tmp2,prop_rev(:,:,istep))
    end do build_propagators
  end subroutine build_small_angle_propagators
  !
  subroutine expand_small_angle_propagators(prop_fwd,prop_rev,rm)
    complex(rk), intent(in)  :: prop_fwd(:,:,:) ! Forward halfs of the rotation propagators
    complex(rk), intent(in)  :: prop_rev(:,:,:) ! Reverse halfs of the rotation propagators
    complex(rk), intent(out) :: rm(:,:)         ! Explicit, dense rotation matrix
    !
    integer(ik) :: nstep, istep
    integer(ik) :: sz, icol
    complex(rk) :: rhs(size(rm,dim=1),size(rm,dim=1))
    !
    if (any(ubound(prop_fwd)/=ubound(prop_rev)) .or. &
        size(rm,dim=1)/=size(rm,dim=2) .or. size(rm,dim=1)/=size(prop_fwd,dim=2)) then
      stop 'rotation_tools%expand_small_angle_propagators - bad argument sizes'
    end if
    !
    nstep = size(prop_fwd,dim=3)
    sz    = size(rm,dim=1)
    !
    rhs(:,:) = 0
    fill_columns: do icol=1,sz
      rhs(icol,icol) = 1._rk
    end do fill_columns
    apply_terms: do istep=1,nstep
      call m3d_multiply(prop_fwd(:,:,istep),rhs,rm)
      call m3d_solve(prop_rev(:,:,istep),rm,rhs)
    end do apply_terms
    rm = rhs
  end subroutine expand_small_angle_propagators
  !
  subroutine rt_finite_rotation_matrix(from,to,mult,rm)
    real(xk), intent(in)     :: from(:)  ! Initial theta, phi in radian
    real(xk), intent(in)     :: to  (:)  ! Final theta, phi in radian
    integer(ik), intent(in)  :: mult     ! Multiplicity, 2*j+1 (half-integer j is OK)
    complex(rk), intent(out) :: rm(:,:)  ! Finite-rotation matrix, transforming spherical
                                         ! harmonics at [from] to spherical harmonics at [to]
    !
    complex(rk) :: rm_from(mult,mult), rm_to(mult,mult)
    integer(ik) :: im1, im2
    !
    !  Conceptually, we need two rotations:
    !  1) From the initial point [from] to the lab system
    !  2) From the lab system to the final point [to]
    !
    !  It looks like there is an error in L&L expression for the Wigner functions,
    !  so that we have to flip the sign of the beta Euler angle to get the correct
    !  rotation matrix. Oh well.
    !
    call MathYJMRotationMatrix(euler_angles=real((/from(2),-from(1),0._xk/),kind=rk),mult=mult,mat=rm_from)
    call MathYJMRotationMatrix(euler_angles=real((/to(2),  -to(1),  0._xk/),kind=rk),mult=mult,mat=rm_to)
    rm = matmul(conjg(rm_to),transpose(rm_from))
    !
    !  MathYJMRotationMatrix assumes harmonics use L&L III phase conventions.
    !  Our propagator uses Mathematica phase conventions, which differ by
    !  an extra factor (-I)**(L+2M). Let's fix up the matrix now.
    !
    !  The factor we need is:  (-I)**(L+2*M1)/(-I)**(L+2*M2) = (-1)**(M1-M2),
    !  so it's sufficient to change the sign of the rotation matrix whenever 
    !  m1 and m2 have different parities.
    !
    fix_m2: do im2=1,mult,2
      fix_m1_even: do im1=2,mult,2
        rm(im1,im2) = -rm(im1,im2)
      end do fix_m1_even
      if (im2+1>mult) cycle fix_m2
      fix_m1_odd: do im1=1,mult,2
        rm(im1,im2+1) = -rm(im1,im2+1)
      end do fix_m1_odd
    end do fix_m2
  end subroutine rt_finite_rotation_matrix
  !
  subroutine rt_rotate_bruteforce(psi,from,to,left)
    type(sd_wfn), intent(inout) :: psi       ! Wavefunction we act upon
    real(xk), intent(in)        :: from(:)   ! Initial (th,ph) in radian
    real(xk), intent(in)        :: to  (:)   ! Final   (th,ph)
    logical, intent(in)         :: left      ! Set to .true. if rotating left wavefunction
    !
    integer(ik) :: lval, mval, ispin
    integer(ik) :: mmin, mmax
    integer(ik) :: ir1, irn, irc, irx
    complex(rk) :: rm(-sd_lmax:sd_lmax,-sd_lmax:sd_lmax)
    complex(rk) :: r_tmp (sd_mmin:sd_mmax,rt_blocksize)
    complex(rk) :: r_tmp2(sd_mmin:sd_mmax,rt_blocksize)
    !
    if (size(from)/=2 .or. size(to)/=2) stop 'rotation_tools%rt_rotate_bruteforce - bad dimensions'
    !
    !  Detect a no-op
    !
    if (all(abs(from-to)<=spacing(1._xk))) return
    !
    call TimerStart('Rotate brute force')
    !$omp parallel default(none) &
    !$omp& shared(sd_lmax,sd_mmin,sd_mmax,sd_nradial,sd_nspin,rt_blocksize) &
    !$omp& shared(from,to,psi,left) &
    !$omp& private(lval,rm,mmin,mmax,ir1,irn,irc,irx,r_tmp,r_tmp2)
    !$omp do
    process_l: do lval=sd_lmax,0,-1
      call rt_finite_rotation_matrix(from,to,2*lval+1,rm(-lval:lval,-lval:lval))
      mmin = max(sd_mmin,-lval)
      mmax = min(sd_mmax, lval)
      process_radial: do ir1=1,sd_nradial,rt_blocksize
        irn  = min(sd_nradial,ir1+rt_blocksize-1)
        irc  = irn - ir1 + 1
        process_spin: do ispin=1,sd_nspin
          ! Fetch wavefunction to cache
          fetch_radial: do mval=mmin,mmax
            r_tmp(mval,:irc) = psi%wfn(ir1:irn,ispin,lval,mval)
          end do fetch_radial
          ! If we are dealing with the left wavefunction, we need to conjugate it here,
          ! then again after rotation is complete
          if (left) r_tmp(mmin:mmax,:irc) = conjg(r_tmp(mmin:mmax,:irc))
          ! Apply the rotation to each R point in turn
          process_radial_cache: do irx=1,irc
            r_tmp2(mmin:mmax,irx) = matmul(rm(mmin:mmax,mmin:mmax),r_tmp(mmin:mmax,irx))
          end do process_radial_cache
          ! Undo the conjugation for the left wavefunction
          if (left) r_tmp2(mmin:mmax,:irc) = conjg(r_tmp2(mmin:mmax,:irc))
          ! Store transformed wavefunction back
          store_radial: do mval=mmin,mmax
            psi%wfn(ir1:irn,ispin,lval,mval) = r_tmp2(mval,:irc)
          end do store_radial
        end do process_spin
      end do process_radial
    end do process_l
    !$omp end do
    !$omp end parallel
    call TimerStop('Rotate brute force')
  end subroutine rt_rotate_bruteforce
  !
end module rotation_tools
