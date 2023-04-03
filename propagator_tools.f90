!
!  Various components of the propagator
!  Names of externally-visible functions in this module follow the following convention:
!
!   pt_fwd_*_n - Evalute [1-I*H*dt]
!   pt_rev_*_n - Evalute [1+I*H*dt]^{-1}
!   pt_fwd_*_t - Evalute [1+I*H^T*dt]
!   pt_rev_*_t - Evalute [1-I*H^T*dt]^{-1}
!
!  where H is a component of the total Hamilonian. Note that all the transpose operators
!  change the sign in front of the Hamiltonian matrix, which is equivalent to time
!  running backwards for the left eigenfunctions. This is not a typo.
!
!  The currenly implemented components are:
!
!   atomic  - Field-free atomic Hamiltonian
!   laser   - Laser coupling Hamiltonian
!
!  Note that in all laser coupling terms, we use vector potential A in units 
!  where the speed of light c is unity.
!
!  Although the non-adiabatic term due to coordinate-system rotation would formally
!  belong here is well, it's memory access pattern and calling conventions are
!  quite different. It is therefore implemented in rotation_tools module instead
!
module propagator_tools
  use accuracy
  use constants
  use timer
  use spherical_data
  use tridiagonal_tools
  use lapack
  use bicg_tools
  implicit none
  private
  public pt_fwd_atomic_n, pt_rev_atomic_n
  public pt_fwd_atomic_t, pt_rev_atomic_t
  public pt_fwd_laser_n, pt_rev_laser_n
  public pt_fwd_laser_t, pt_rev_laser_t
  public pt_mix_solver, pt_chunk, pt_sense_real, pt_eps_dt, pt_force_par_l
  !
  character(len=20), save :: pt_mix_solver  = 'default' ! Solver to use for L=0 Hmix inverse. Can be:
                                                        ! 'default' - will use the fastest solver
                                                        ! 'SM'      - tri-diagonal solver with Sherman-Morrison
                                                        !             correction.
                                                        ! 'bi-CG'   - bi-conjugate gradient solver
                                                        ! 'lapack'  - Lapack zgesv 
                                                        ! There is no sane reason to use 'lapack' except
                                                        ! for debugging the much faster, and just as
                                                        ! accurate bi-CG solver.
  integer, save           :: pt_chunk       = 1         ! Size of dynamic chunk in loop scheduling
  logical, save           :: pt_sense_real  = .true.    ! Use special-case optimizations for real time steps.
                                                        ! Real time step may lead to better optimizations and/or
                                                        ! switch to a different algorithm, depending on the routine
                                                        ! There is no reason to use .false. except for debugging.
  real(rk), save          :: pt_eps_dt      = -1._rk    ! Time-step tolerance for caching atomic propagator.
                                                        ! Negative value will use spacing(1000._xk)
  logical, save           :: pt_force_par_l = .false.   ! Set to .true. to force use of L-parallel routine in
                                                        ! laser propagator. This may be necessary if compiler
                                                        ! generates lousy code for the par_m routine.
                       
  !
  !  Cache of atomic propagators for the iverse half of the operator.
  !  The last index is angular momentum; the first two are decomposed tri-diagonal matrix.
  !  The cache must be updated each time step is changed.
  !
  complex(xk), save              :: cache_atomic_dt_n
  complex(xk), save              :: cache_atomic_dt_t
  complex(rk), save, allocatable :: cache_atomic_leqf_n(:,:,:)
  complex(rk), save, allocatable :: cache_atomic_leqf_t(:,:,:)
  !
  contains
  !
  !  Calculates forward half of the right atomic propagator: [1-I*Hat*dt] Psi
  !    Hat = -(1/(2me)) m2^-1 d2 + v
  !  Hat operator does not couple neither M nor L blocks. It's evaluation
  !  can proceed in parallel, with no synchronization points. It can however
  !  mix spin channels; these will have to be implemented later.
  !
  subroutine pt_fwd_atomic_n(psi,dt)
    type(sd_wfn), intent(inout) :: psi ! Wavefunction we act upon
    complex(xk), intent(in)     :: dt  ! Timestep, in atomic units. Use at least double precision here
                                       ! We allow dt to be complex to enable solving the
                                       ! ground-state problem using imaginary time propagation
    !
    integer(ik) :: lval, mval
    complex(rk) :: tmp(sd_nradial), hpsi(sd_nradial)
    !
    call TimerStart('Atomic fwd n')
    if (sd_nspin/=1) stop 'propagator_tools%pt_fwd_atomic_n - spin-orbit coupling not implemented'
    !
    !$omp parallel default(none) &
    !$omp& shared(sd_mmin,sd_mmax,sd_lmax,sd_d2n_l0,sd_m2nf_l0,sd_d2n_lx,sd_m2nf_lx,dt,psi) &
    !$omp& shared(sd_pottab,sd_capped,sd_captab,sd_cap_start) &
    !$omp& shared(pt_chunk,pt_sense_real) &
    !$omp& private(mval,lval,tmp,hpsi)
    process_m: do mval=sd_mmin,sd_mmax
      !$omp do schedule(dynamic,pt_chunk)
      process_l: do lval=abs(mval),sd_lmax
        !
        !  Evaluate the Laplacian; this one branches according to L
        !
        if (lval==0) then
          call m3d_multiply(sd_d2n_l0,psi%wfn(:,1,lval,mval),tmp)
          call m3d_solve(sd_m2nf_l0,tmp,hpsi)
        else
          call m3d_multiply(sd_d2n_lx,psi%wfn(:,1,lval,mval),tmp)
          call m3d_solve(sd_m2nf_lx,tmp,hpsi)
        end if
        hpsi = (-0.5_rk/electron_mass) * hpsi
        hpsi = hpsi + sd_pottab(:,1,lval)*psi%wfn(:,1,lval,mval)
        if (sd_capped) then
          hpsi(sd_cap_start:) = hpsi(sd_cap_start:) &
                              + sd_captab(sd_cap_start:)*psi%wfn(sd_cap_start:,1,lval,mval)
        end if
        if (pt_sense_real .and. aimag(dt)==0._xk) then
          hpsi = (-(0,1)*real(dt,kind=rk))*hpsi
        else
          hpsi = (-(0,1)*cmplx(dt,kind=rk))*hpsi
        end if
        psi%wfn(:,1,lval,mval) = psi%wfn(:,1,lval,mval) + hpsi
      end do process_l
      !$omp end do nowait
    end do process_m
    !$omp end parallel
    call TimerStop('Atomic fwd n')
  end subroutine pt_fwd_atomic_n
  !
  !  Calculates implicit half of the right atomic propagator: Psi' = [1+I*Hat*dt]^{-1} Psi
  !  The propagator is evaluated by solving a tri-diagonal system of linear equations:
  !
  !   [ m2 (1+I dt v) - (I dt/(2 me)) d2 ] Psi' = m2 Psi
  !
  !  Note that right now we construct the linear system on each time step; as long as
  !  time step remains constant during the simulation, this linear system matrix can
  !  be precomputed and stored, if necessary.
  !
  subroutine update_atomic_cache_n(dt)
    complex(xk), intent(in) :: dt  ! Timestep, in atomic units
    !
    integer(ik) :: lval, alloc
    complex(rk) :: tmp(sd_nradial)
    complex(rk) :: leq(3,sd_nradial)
    real(xk)    :: eps_dt
    !
    if (allocated(cache_atomic_leqf_n)) then
      eps_dt = pt_eps_dt
      if (eps_dt<0) eps_dt = spacing(1000._xk)
      if (abs(cache_atomic_dt_n-dt)<eps_dt) return
      ! write (out,"('N: delta = ',2g30.20,' dt = ',2g30.20)") cache_atomic_dt_n-dt, dt
      if (allocated(cache_atomic_leqf_n)) deallocate(cache_atomic_leqf_n)
    end if
    !
    call TimerStart('Atomic rev n (update)')
    if (sd_nspin/=1) stop 'propagator_tools%update_atomic_cache_n - spin-orbit coupling not implemented'
    cache_atomic_dt_n = dt
    allocate (cache_atomic_leqf_n(3,sd_nradial,0:sd_lmax),stat=alloc)
    if (alloc/=0) then
      stop 'propagator_tools%update_atomic_cache_n - allocation failed'
    end if
    !
    !$omp parallel default(none) &
    !$omp& shared(sd_lmax,dt,sd_m2n_l0,sd_d2n_l0,sd_m2n_lx,sd_d2n_lx) &
    !$omp& shared(sd_pottab,sd_capped,sd_captab,sd_cap_start) &
    !$omp& shared(cache_atomic_leqf_n,pt_sense_real,pt_chunk) &
    !$omp& private(lval,tmp,leq)
    !$omp do schedule(dynamic,pt_chunk)
    process_l: do lval=0,sd_lmax
      !
      !  Construct linear system matrix 
      !
      if (pt_sense_real .and. aimag(dt)==0._rk) then
        tmp = 1 + (0,1)*real(dt,kind=rk)*sd_pottab(:,1,lval)
        if (sd_capped) then
          tmp(sd_cap_start:) = tmp(sd_cap_start:) + (0,1)*real(dt,kind=rk)*sd_captab(sd_cap_start:)
        end if
      else
        tmp = 1 + (0,1)*cmplx(dt,kind=rk)*sd_pottab(:,1,lval)
        if (sd_capped) then
          tmp(sd_cap_start:) = tmp(sd_cap_start:) + (0,1)*cmplx(dt,kind=rk)*sd_captab(sd_cap_start:)
        end if
      end if
      if (lval==0) then
        call m3d_right_scale(sd_m2n_l0,tmp,leq)
        if (pt_sense_real .and. aimag(dt)==0._rk) then
          leq = leq - ((0,1)*real(dt,kind=rk)/(2*electron_mass))*sd_d2n_l0
        else
          leq = leq - ((0,1)*cmplx(dt,kind=rk)/(2*electron_mass))*sd_d2n_l0
        end if
      else
        call m3d_right_scale(sd_m2n_lx,tmp,leq)
        if (pt_sense_real .and. aimag(dt)==0._rk) then
          leq = leq - ((0,1)*real(dt,kind=rk)/(2*electron_mass))*sd_d2n_lx
        else
          leq = leq - ((0,1)*cmplx(dt,kind=rk)/(2*electron_mass))*sd_d2n_lx
        end if
      end if
      call m3d_decompose(leq,cache_atomic_leqf_n(:,:,lval))
    end do process_l
    !$omp end do
    !$omp end parallel
    call TimerStop('Atomic rev n (update)')
  end subroutine update_atomic_cache_n
  !
  subroutine pt_rev_atomic_n(psi,dt)
    type(sd_wfn), intent(inout) :: psi ! Wavefunction we act upon
    complex(xk), intent(in)     :: dt  ! Timestep, in atomic units
                                       ! We allow dt to be complex to enable solving the
                                       ! ground-state problem using imaginary time propagation
    !
    integer(ik) :: lval, mval
    complex(rk) :: rhs(sd_nradial), tmp(sd_nradial)
    !
    call TimerStart('Atomic rev n')
    if (sd_nspin/=1) stop 'propagator_tools%pt_rev_atomic_n - spin-orbit coupling not implemented'
    call update_atomic_cache_n(dt)
    !
    !$omp parallel default(none) &
    !$omp& shared(sd_mmin,sd_mmax,sd_lmax,dt,sd_m2n_l0,sd_m2n_lx,psi) &
    !$omp& shared(pt_chunk,cache_atomic_leqf_n) &
    !$omp& private(mval,lval,tmp,rhs)
    process_m: do mval=sd_mmin,sd_mmax
      !$omp do schedule(dynamic,pt_chunk)
      process_l: do lval=abs(mval),sd_lmax
        !
        !  Rescale the rhs
        !
        if (lval==0) then
          call m3d_multiply(sd_m2n_l0,psi%wfn(:,1,lval,mval),rhs)
        else
          call m3d_multiply(sd_m2n_lx,psi%wfn(:,1,lval,mval),rhs)
        end if
        call m3d_solve(cache_atomic_leqf_n(:,:,lval),rhs,psi%wfn(:,1,lval,mval))
      end do process_l
      !$omp end do nowait
    end do process_m
    !$omp end parallel
    call TimerStop('Atomic rev n')
  end subroutine pt_rev_atomic_n
  !
  !  Calculates forward half of the left atomic propagator: [1-I*Hat^T*dt] Psi
  !    Hat^t = -(1/(2me)) d2^t m2^t^-1 + v
  !
  subroutine pt_fwd_atomic_t(psi,dt)
    type(sd_wfn), intent(inout) :: psi ! Wavefunction we act upon
    complex(xk), intent(in)     :: dt  ! Timestep, in atomic units
                                       ! We allow dt to be complex to enable solving the
                                       ! ground-state problem using imaginary time propagation
    !
    integer(ik) :: lval, mval
    complex(rk) :: tmp(sd_nradial), hpsi(sd_nradial)
    !
    call TimerStart('Atomic fwd t')
    if (sd_nspin/=1) stop 'propagator_tools%pt_fwd_atomic_t - spin-orbit coupling not implemented'
    !
    !$omp parallel default(none) &
    !$omp& shared(sd_mmin,sd_mmax,sd_lmax,sd_d2t_l0,sd_m2tf_l0,sd_d2t_lx,sd_m2tf_lx,dt,psi) &
    !$omp& shared(sd_pottab,sd_capped,sd_captab,sd_cap_start) &
    !$omp& shared(pt_chunk,pt_sense_real) &
    !$omp& private(mval,lval,tmp,hpsi)
    process_m: do mval=sd_mmin,sd_mmax
      !$omp do schedule(dynamic,pt_chunk)
      process_l: do lval=abs(mval),sd_lmax
        !
        !  Evaluate the Laplacian; this one branches according to L
        !
        if (lval==0) then
          call m3d_solve(sd_m2tf_l0,psi%wfn(:,1,lval,mval),tmp)
          call m3d_multiply(sd_d2t_l0,tmp,hpsi)
        else
          call m3d_solve(sd_m2tf_lx,psi%wfn(:,1,lval,mval),tmp)
          call m3d_multiply(sd_d2t_lx,tmp,hpsi)
        end if
        hpsi = (-0.5_rk/electron_mass) * hpsi
        hpsi = hpsi + sd_pottab(:,1,lval)*psi%wfn(:,1,lval,mval)
        if (sd_capped) then
          hpsi(sd_cap_start:) = hpsi(sd_cap_start:) &
                              + conjg(sd_captab(sd_cap_start:))*psi%wfn(sd_cap_start:,1,lval,mval)
        end if
        if (pt_sense_real .and. aimag(dt)==0._rk) then
          hpsi = (+(0,1)*real(dt,kind=rk))*hpsi
        else
          hpsi = (+(0,1)*cmplx(dt,kind=rk))*hpsi
        end if
        psi%wfn(:,1,lval,mval) = psi%wfn(:,1,lval,mval) + hpsi
      end do process_l
      !$omp end do nowait
    end do process_m
    !$omp end parallel
    call TimerStop('Atomic fwd t')
  end subroutine pt_fwd_atomic_t
  !
  !  Calculates implicit half of the left atomic propagator: Psi' = [1+I*Hat^t*dt]^{-1} Psi
  !  The propagator is evaluated by solving a tri-diagonal system of linear equations:
  !
  !   [ (1+I dt v) m2t - (I dt/(2 me)) d2t ] Y = Psi
  !   Psi' = m2t Y
  !
  subroutine update_atomic_cache_t(dt)
    complex(xk), intent(in)     :: dt  ! Timestep, in atomic units
    !
    integer(ik) :: lval, alloc
    complex(rk) :: tmp(sd_nradial)
    complex(rk) :: leq(3,sd_nradial)
    real(xk)    :: eps_dt
    !
    if (allocated(cache_atomic_leqf_t)) then
      eps_dt = pt_eps_dt
      if (eps_dt<0) eps_dt = spacing(1000._xk)
      if (abs(cache_atomic_dt_t-dt)<eps_dt) return
      ! write (out,"('T: delta = ',2g30.20,' dt = ',2g30.20)") cache_atomic_dt_t-dt, dt
      if (allocated(cache_atomic_leqf_t)) deallocate(cache_atomic_leqf_t)
    end if
    !
    call TimerStart('Atomic rev t (update)')
    if (sd_nspin/=1) stop 'propagator_tools%update_atomic_cache_t - spin-orbit coupling not implemented'
    cache_atomic_dt_t = dt
    allocate (cache_atomic_leqf_t(3,sd_nradial,0:sd_lmax),stat=alloc)
    if (alloc/=0) then
      stop 'propagator_tools%update_atomic_cache_t - allocation failed'
    end if
    !
    !$omp parallel default(none) &
    !$omp& shared(sd_lmax,dt,sd_m2t_l0,sd_d2t_l0,sd_m2t_lx,sd_d2t_lx) &
    !$omp& shared(sd_pottab,sd_capped,sd_captab,sd_cap_start) &
    !$omp& shared(pt_chunk,pt_sense_real,cache_atomic_leqf_t) &
    !$omp& private(lval,tmp,leq)
    !$omp do schedule(dynamic,pt_chunk)
    process_l: do lval=0,sd_lmax
      !
      !  Construct linear system matrix
      !
      if (pt_sense_real .and. aimag(dt)==0._rk) then
        tmp = 1 - (0,1)*real(dt,kind=rk)*sd_pottab(:,1,lval)
        if (sd_capped) then
          tmp(sd_cap_start:) = tmp(sd_cap_start:) - (0,1)*real(dt,kind=rk)*conjg(sd_captab(sd_cap_start:))
        end if
      else
        tmp = 1 - (0,1)*cmplx(dt,kind=rk)*sd_pottab(:,1,lval)
        if (sd_capped) then
          tmp(sd_cap_start:) = tmp(sd_cap_start:) - (0,1)*cmplx(dt,kind=rk)*conjg(sd_captab(sd_cap_start:))
        end if
      end if
      if (lval==0) then
        call m3d_left_scale(tmp,sd_m2t_l0,leq)
        if (pt_sense_real .and. aimag(dt)==0._rk) then
          leq = leq + ((0,1)*real(dt,kind=rk)/(2*electron_mass))*sd_d2t_l0
        else
          leq = leq + ((0,1)*cmplx(dt,kind=rk)/(2*electron_mass))*sd_d2t_l0
        end if
      else
        call m3d_left_scale(tmp,sd_m2t_lx,leq)
        if (pt_sense_real .and. aimag(dt)==0._rk) then
          leq = leq + ((0,1)*real(dt,kind=rk)/(2*electron_mass))*sd_d2t_lx
        else
          leq = leq + ((0,1)*cmplx(dt,kind=rk)/(2*electron_mass))*sd_d2t_lx
        end if
      end if
      call m3d_decompose(leq,cache_atomic_leqf_t(:,:,lval))
    end do process_l
    !$omp end do
    !$omp end parallel
    call TimerStop('Atomic rev t (update)')
  end subroutine update_atomic_cache_t
  !
  subroutine pt_rev_atomic_t(psi,dt)
    type(sd_wfn), intent(inout) :: psi ! Wavefunction we act upon
    complex(xk), intent(in)     :: dt  ! Timestep, in atomic units
                                       ! We allow dt to be complex to enable solving the
                                       ! ground-state problem using imaginary time propagation
    !
    integer(ik) :: lval, mval
    complex(rk) :: tmp(sd_nradial)
    !
    call TimerStart('Atomic rev t')
    if (sd_nspin/=1) stop 'propagator_tools%pt_rev_atomic_t - spin-orbit coupling not implemented'
    call update_atomic_cache_t(dt)
    !
    !$omp parallel default(none) &
    !$omp& shared(sd_mmin,sd_mmax,sd_lmax,dt,sd_m2t_l0,sd_m2t_lx,psi) &
    !$omp& shared(pt_chunk,cache_atomic_leqf_t) &
    !$omp& private(mval,lval,tmp)
    process_m: do mval=sd_mmin,sd_mmax
      !$omp do schedule(dynamic,pt_chunk)
      process_l: do lval=abs(mval),sd_lmax
        !
        !  Recall the linear system matrix and rescale the rhs
        !
        call m3d_solve(cache_atomic_leqf_t(:,:,lval),psi%wfn(:,1,lval,mval),tmp)
        if (lval==0) then
          call m3d_multiply(sd_m2t_l0,tmp,psi%wfn(:,1,lval,mval))
        else
          call m3d_multiply(sd_m2t_lx,tmp,psi%wfn(:,1,lval,mval))
        end if
      end do process_l
      !$omp end do nowait
    end do process_m
    !$omp end parallel
    call TimerStop('Atomic rev t')
  end subroutine pt_rev_atomic_t
  !
  !  Calculates forward half of the right laser propagator.
  !  Hlas contains three distinct tri-diagonal components, which are processed
  !  separaterly. The components are: Hlang, Hlmix, and Hla2, detailed below.
  !  Hlang and Hlmix are decomposed in 2x2 subblocks, are are processed in two
  !  waves, offset by 1. Hla2 is purely diagonal, and is computed in an extra
  !  pass. Note that neither of the three laser interaction components couple
  !  different Ms.
  !
  !  Note that the forward and reverse half of the propagator use exactly the
  !  same operator blocks, but apply them in a different order to maintain 
  !  exact unitarity.
  !
  subroutine pt_fwd_laser_n(psi,a,dt)
    type(sd_wfn), intent(inout) :: psi ! Wavefunction we act upon
    real(xk), intent(in)        :: a   ! Vector potential along the Z axis
    complex(xk), intent(in)     :: dt  ! Timestep, in atomic units
    !
    if (sd_mmax-sd_mmin+1>sd_lmax+1 .and. .not. pt_force_par_l) then
      call pt_fwd_laser_n_parallel_m(psi,a,dt)
    else
      call pt_fwd_laser_n_parallel_l(psi,a,dt)
    end if
  end subroutine pt_fwd_laser_n
  !
  subroutine pt_fwd_laser_n_parallel_m(psi,a,dt)
    type(sd_wfn), intent(inout) :: psi ! Wavefunction we act upon
    real(xk), intent(in)        :: a   ! Vector potential along the Z axis
    complex(xk), intent(in)     :: dt  ! Timestep, in atomic units
    !
    integer(ik) :: lval, mval, sval
    ! Intel Fortran serializes on the memory allocator. To get reasonable
    ! parallel performance, we are forced to pull scratch allocation upstream
    ! from the worker routines :(
    complex(rk) :: lang_c1c(sd_nradial), lang_c2c(sd_nradial), lang_c3c(sd_nradial)
    real(rk)    :: lang_c1r(sd_nradial), lang_c2r(sd_nradial), lang_c3r(sd_nradial)
    complex(rk) :: lang_xl(sd_nradial), lang_xp(sd_nradial)
    complex(rk) :: lmix_zp(sd_nradial,2), lmix_zm(sd_nradial,2)
    complex(rk) :: lmix_xp(sd_nradial,2), lmix_xm(sd_nradial,2)
    complex(rk) :: lmix_rhs(sd_nradial)
    complex(rk) :: lmix_leqp_c(3,sd_nradial), lmix_leqpf_c(3,sd_nradial)
    complex(rk) :: lmix_leqm_c(3,sd_nradial), lmix_leqmf_c(3,sd_nradial)
    real(rk)    :: lmix_leqp_r(3,sd_nradial), lmix_leqpf_r(3,sd_nradial)
    real(rk)    :: lmix_leqm_r(3,sd_nradial), lmix_leqmf_r(3,sd_nradial)
    !
    if (a==0._rk) return ! Nothing to do!
    call TimerStart('Laser fwd n (par M)')
    !$omp parallel default(none) &
    !$omp& shared(sd_mmin,sd_mmax,sd_lmax,sd_nspin,a,dt,psi,pt_chunk) &
    !$omp& private(mval,lval,sval) &
    !$omp& private(lang_c1c,lang_c2c,lang_c3c,lang_c1r,lang_c2r,lang_c3r,lang_xl,lang_xp) &
    !$omp& private(lmix_zp,lmix_zm,lmix_xp,lmix_xm,lmix_rhs,lmix_leqp_c,lmix_leqpf_c) &
    !$omp& private(lmix_leqm_c,lmix_leqmf_c,lmix_leqp_r,lmix_leqpf_r,lmix_leqm_r,lmix_leqmf_r)
    !$omp do schedule(dynamic,pt_chunk)
    process_m: do mval=sd_mmin,sd_mmax
      process_first_wave_l: do lval=abs(mval),sd_lmax-1,2
        spin_part_1: do sval=1,sd_nspin
          call op_lang_n(mval,lval,a,dt,psi%wfn(:,sval,lval,mval),psi%wfn(:,sval,lval+1,mval), &
                         lang_c1c,lang_c2c,lang_c3c,lang_c1r,lang_c2r,lang_c3r,lang_xl,lang_xp)
          call op_lmix_n(mval,lval,a,dt,psi%wfn(:,sval,lval,mval),psi%wfn(:,sval,lval+1,mval), &
                         lmix_zp,lmix_zm,lmix_xp,lmix_xm,lmix_rhs,lmix_leqp_c,lmix_leqpf_c, &
                         lmix_leqm_c,lmix_leqmf_c,lmix_leqp_r,lmix_leqpf_r,lmix_leqm_r,lmix_leqmf_r)
        end do spin_part_1
      end do process_first_wave_l
      !
      process_second_wave_l: do lval=abs(mval)+1,sd_lmax-1,2
        spin_part_2: do sval=1,sd_nspin
          call op_lang_n(mval,lval,a,dt,psi%wfn(:,sval,lval,mval),psi%wfn(:,sval,lval+1,mval), &
                         lang_c1c,lang_c2c,lang_c3c,lang_c1r,lang_c2r,lang_c3r,lang_xl,lang_xp)
          call op_lmix_n(mval,lval,a,dt,psi%wfn(:,sval,lval,mval),psi%wfn(:,sval,lval+1,mval), &
                         lmix_zp,lmix_zm,lmix_xp,lmix_xm,lmix_rhs,lmix_leqp_c,lmix_leqpf_c, &
                         lmix_leqm_c,lmix_leqmf_c,lmix_leqp_r,lmix_leqpf_r,lmix_leqm_r,lmix_leqmf_r)
        end do spin_part_2
      end do process_second_wave_l
      !
      process_diagonal_term: do lval=abs(mval),sd_lmax
        spin_part_3: do sval=1,sd_nspin
          call op_la2(a,dt,psi%wfn(:,sval,lval,mval))
        end do spin_part_3
      end do process_diagonal_term
    end do process_m
    !$omp end do nowait
    !$omp end parallel
    call TimerStop('Laser fwd n (par M)')
  end subroutine pt_fwd_laser_n_parallel_m
  !
  subroutine pt_fwd_laser_n_parallel_l(psi,a,dt)
    type(sd_wfn), intent(inout) :: psi ! Wavefunction we act upon
    real(xk), intent(in)        :: a   ! Vector potential along the Z axis
    complex(xk), intent(in)     :: dt  ! Timestep, in atomic units
    !
    integer(ik) :: lval, mval, sval
    ! Intel Fortran serializes on the memory allocator. To get reasonable
    ! parallel performance, we are forced to pull scratch allocation upstream
    ! from the worker routines :(
    complex(rk) :: lang_c1c(sd_nradial), lang_c2c(sd_nradial), lang_c3c(sd_nradial)
    real(rk)    :: lang_c1r(sd_nradial), lang_c2r(sd_nradial), lang_c3r(sd_nradial)
    complex(rk) :: lang_xl(sd_nradial), lang_xp(sd_nradial)
    complex(rk) :: lmix_zp(sd_nradial,2), lmix_zm(sd_nradial,2)
    complex(rk) :: lmix_xp(sd_nradial,2), lmix_xm(sd_nradial,2)
    complex(rk) :: lmix_rhs(sd_nradial)
    complex(rk) :: lmix_leqp_c(3,sd_nradial), lmix_leqpf_c(3,sd_nradial)
    complex(rk) :: lmix_leqm_c(3,sd_nradial), lmix_leqmf_c(3,sd_nradial)
    real(rk)    :: lmix_leqp_r(3,sd_nradial), lmix_leqpf_r(3,sd_nradial)
    real(rk)    :: lmix_leqm_r(3,sd_nradial), lmix_leqmf_r(3,sd_nradial)
    !
    if (a==0._rk) return ! Nothing to do!
    call TimerStart('Laser fwd n (par L)')
    !$omp parallel default(none) &
    !$omp& shared(sd_mmin,sd_mmax,sd_lmax,sd_nspin,a,dt,psi,pt_chunk) &
    !$omp& private(mval,lval,sval) &
    !$omp& private(lang_c1c,lang_c2c,lang_c3c,lang_c1r,lang_c2r,lang_c3r,lang_xl,lang_xp) &
    !$omp& private(lmix_zp,lmix_zm,lmix_xp,lmix_xm,lmix_rhs,lmix_leqp_c,lmix_leqpf_c) &
    !$omp& private(lmix_leqm_c,lmix_leqmf_c,lmix_leqp_r,lmix_leqpf_r,lmix_leqm_r,lmix_leqmf_r)
    process_m: do mval=sd_mmin,sd_mmax
      !$omp do schedule(dynamic,pt_chunk)
      process_first_wave_l: do lval=abs(mval),sd_lmax-1,2
        spin_part_1: do sval=1,sd_nspin
          call op_lang_n(mval,lval,a,dt,psi%wfn(:,sval,lval,mval),psi%wfn(:,sval,lval+1,mval), &
                         lang_c1c,lang_c2c,lang_c3c,lang_c1r,lang_c2r,lang_c3r,lang_xl,lang_xp)
          call op_lmix_n(mval,lval,a,dt,psi%wfn(:,sval,lval,mval),psi%wfn(:,sval,lval+1,mval), &
                         lmix_zp,lmix_zm,lmix_xp,lmix_xm,lmix_rhs,lmix_leqp_c,lmix_leqpf_c, &
                         lmix_leqm_c,lmix_leqmf_c,lmix_leqp_r,lmix_leqpf_r,lmix_leqm_r,lmix_leqmf_r)
        end do spin_part_1
      end do process_first_wave_l
      !$omp end do
      !
      !$omp do schedule(dynamic,pt_chunk)
      process_second_wave_l: do lval=abs(mval)+1,sd_lmax-1,2
        spin_part_2: do sval=1,sd_nspin
          call op_lang_n(mval,lval,a,dt,psi%wfn(:,sval,lval,mval),psi%wfn(:,sval,lval+1,mval), &
                         lang_c1c,lang_c2c,lang_c3c,lang_c1r,lang_c2r,lang_c3r,lang_xl,lang_xp)
          call op_lmix_n(mval,lval,a,dt,psi%wfn(:,sval,lval,mval),psi%wfn(:,sval,lval+1,mval), &
                         lmix_zp,lmix_zm,lmix_xp,lmix_xm,lmix_rhs,lmix_leqp_c,lmix_leqpf_c, &
                         lmix_leqm_c,lmix_leqmf_c,lmix_leqp_r,lmix_leqpf_r,lmix_leqm_r,lmix_leqmf_r)
        end do spin_part_2
      end do process_second_wave_l
      !$omp end do
      !
      !$omp do schedule(dynamic,pt_chunk)
      process_diagonal_term: do lval=abs(mval),sd_lmax
        spin_part_3: do sval=1,sd_nspin
          call op_la2(a,dt,psi%wfn(:,sval,lval,mval))
        end do spin_part_3
      end do process_diagonal_term
      !$omp end do nowait
    end do process_m
    !$omp end parallel
    call TimerStop('Laser fwd n (par L)')
  end subroutine pt_fwd_laser_n_parallel_l
  !
  subroutine pt_rev_laser_n(psi,a,dt)
    type(sd_wfn), intent(inout) :: psi ! Wavefunction we act upon
    real(xk), intent(in)        :: a   ! Vector potential along the Z axis
    complex(xk), intent(in)     :: dt  ! Timestep, in atomic units
    !
    if (sd_mmax-sd_mmin+1>sd_lmax+1 .and. .not. pt_force_par_l) then
      call pt_rev_laser_n_parallel_m(psi,a,dt)
    else
      call pt_rev_laser_n_parallel_l(psi,a,dt)
    end if
  end subroutine pt_rev_laser_n
  !
  subroutine pt_rev_laser_n_parallel_m(psi,a,dt)
    type(sd_wfn), intent(inout) :: psi ! Wavefunction we act upon
    real(xk), intent(in)        :: a   ! Vector potential along the Z axis
    complex(xk), intent(in)     :: dt  ! Timestep, in atomic units
    !
    integer(ik) :: lval, mval, sval
    ! Intel Fortran serializes on the memory allocator. To get reasonable
    ! parallel performance, we are forced to pull scratch allocation upstream
    ! from the worker routines :(
    complex(rk) :: lang_c1c(sd_nradial), lang_c2c(sd_nradial), lang_c3c(sd_nradial)
    real(rk)    :: lang_c1r(sd_nradial), lang_c2r(sd_nradial), lang_c3r(sd_nradial)
    complex(rk) :: lang_xl(sd_nradial), lang_xp(sd_nradial)
    complex(rk) :: lmix_zp(sd_nradial,2), lmix_zm(sd_nradial,2)
    complex(rk) :: lmix_xp(sd_nradial,2), lmix_xm(sd_nradial,2)
    complex(rk) :: lmix_rhs(sd_nradial)
    complex(rk) :: lmix_leqp_c(3,sd_nradial), lmix_leqpf_c(3,sd_nradial)
    complex(rk) :: lmix_leqm_c(3,sd_nradial), lmix_leqmf_c(3,sd_nradial)
    real(rk)    :: lmix_leqp_r(3,sd_nradial), lmix_leqpf_r(3,sd_nradial)
    real(rk)    :: lmix_leqm_r(3,sd_nradial), lmix_leqmf_r(3,sd_nradial)
    !
    if (a==0._rk) return ! Nothing to do!
    call TimerStart('Laser rev n (par M)')
    !$omp parallel default(none) &
    !$omp& shared(sd_mmin,sd_mmax,sd_lmax,sd_nspin,a,dt,psi,pt_chunk) &
    !$omp& private(mval,lval,sval) &
    !$omp& private(lang_c1c,lang_c2c,lang_c3c,lang_c1r,lang_c2r,lang_c3r,lang_xl,lang_xp) &
    !$omp& private(lmix_zp,lmix_zm,lmix_xp,lmix_xm,lmix_rhs,lmix_leqp_c,lmix_leqpf_c) &
    !$omp& private(lmix_leqm_c,lmix_leqmf_c,lmix_leqp_r,lmix_leqpf_r,lmix_leqm_r,lmix_leqmf_r)
    !$omp do schedule(dynamic,pt_chunk)
    process_m: do mval=sd_mmin,sd_mmax
      process_diagonal_term: do lval=abs(mval),sd_lmax
        spin_part_3: do sval=1,sd_nspin
          call op_la2(a,dt,psi%wfn(:,sval,lval,mval))
        end do spin_part_3
      end do process_diagonal_term
      !
      process_second_wave_l: do lval=abs(mval)+1,sd_lmax-1,2
        spin_part_2: do sval=1,sd_nspin
          call op_lmix_n(mval,lval,a,dt,psi%wfn(:,sval,lval,mval),psi%wfn(:,sval,lval+1,mval), &
                         lmix_zp,lmix_zm,lmix_xp,lmix_xm,lmix_rhs,lmix_leqp_c,lmix_leqpf_c, &
                         lmix_leqm_c,lmix_leqmf_c,lmix_leqp_r,lmix_leqpf_r,lmix_leqm_r,lmix_leqmf_r)
          call op_lang_n(mval,lval,a,dt,psi%wfn(:,sval,lval,mval),psi%wfn(:,sval,lval+1,mval), &
                         lang_c1c,lang_c2c,lang_c3c,lang_c1r,lang_c2r,lang_c3r,lang_xl,lang_xp)
        end do spin_part_2
      end do process_second_wave_l
      !
      process_first_wave_l: do lval=abs(mval),sd_lmax-1,2
        spin_part_1: do sval=1,sd_nspin
          call op_lmix_n(mval,lval,a,dt,psi%wfn(:,sval,lval,mval),psi%wfn(:,sval,lval+1,mval), &
                         lmix_zp,lmix_zm,lmix_xp,lmix_xm,lmix_rhs,lmix_leqp_c,lmix_leqpf_c, &
                         lmix_leqm_c,lmix_leqmf_c,lmix_leqp_r,lmix_leqpf_r,lmix_leqm_r,lmix_leqmf_r)
          call op_lang_n(mval,lval,a,dt,psi%wfn(:,sval,lval,mval),psi%wfn(:,sval,lval+1,mval), &
                         lang_c1c,lang_c2c,lang_c3c,lang_c1r,lang_c2r,lang_c3r,lang_xl,lang_xp)
        end do spin_part_1
      end do process_first_wave_l
    end do process_m
    !$omp end do nowait
    !$omp end parallel
    call TimerStop('Laser rev n (par M)')
  end subroutine pt_rev_laser_n_parallel_m
  !
  subroutine pt_rev_laser_n_parallel_l(psi,a,dt)
    type(sd_wfn), intent(inout) :: psi ! Wavefunction we act upon
    real(xk), intent(in)        :: a   ! Vector potential along the Z axis
    complex(xk), intent(in)     :: dt  ! Timestep, in atomic units
    !
    integer(ik) :: lval, mval, sval
    ! Intel Fortran serializes on the memory allocator. To get reasonable
    ! parallel performance, we are forced to pull scratch allocation upstream
    ! from the worker routines :(
    complex(rk) :: lang_c1c(sd_nradial), lang_c2c(sd_nradial), lang_c3c(sd_nradial)
    real(rk)    :: lang_c1r(sd_nradial), lang_c2r(sd_nradial), lang_c3r(sd_nradial)
    complex(rk) :: lang_xl(sd_nradial), lang_xp(sd_nradial)
    complex(rk) :: lmix_zp(sd_nradial,2), lmix_zm(sd_nradial,2)
    complex(rk) :: lmix_xp(sd_nradial,2), lmix_xm(sd_nradial,2)
    complex(rk) :: lmix_rhs(sd_nradial)
    complex(rk) :: lmix_leqp_c(3,sd_nradial), lmix_leqpf_c(3,sd_nradial)
    complex(rk) :: lmix_leqm_c(3,sd_nradial), lmix_leqmf_c(3,sd_nradial)
    real(rk)    :: lmix_leqp_r(3,sd_nradial), lmix_leqpf_r(3,sd_nradial)
    real(rk)    :: lmix_leqm_r(3,sd_nradial), lmix_leqmf_r(3,sd_nradial)
    !
    if (a==0._rk) return ! Nothing to do!
    call TimerStart('Laser rev n (par L)')
    !$omp parallel default(none) &
    !$omp& shared(sd_mmin,sd_mmax,sd_lmax,sd_nspin,a,dt,psi,pt_chunk) &
    !$omp& private(mval,lval,sval) &
    !$omp& private(lang_c1c,lang_c2c,lang_c3c,lang_c1r,lang_c2r,lang_c3r,lang_xl,lang_xp) &
    !$omp& private(lmix_zp,lmix_zm,lmix_xp,lmix_xm,lmix_rhs,lmix_leqp_c,lmix_leqpf_c) &
    !$omp& private(lmix_leqm_c,lmix_leqmf_c,lmix_leqp_r,lmix_leqpf_r,lmix_leqm_r,lmix_leqmf_r)
    process_m: do mval=sd_mmin,sd_mmax
      !$omp do schedule(dynamic,pt_chunk)
      process_diagonal_term: do lval=abs(mval),sd_lmax
        spin_part_3: do sval=1,sd_nspin
          call op_la2(a,dt,psi%wfn(:,sval,lval,mval))
        end do spin_part_3
      end do process_diagonal_term
      !$omp end do
      !
      !$omp do schedule(dynamic,pt_chunk)
      process_second_wave_l: do lval=abs(mval)+1,sd_lmax-1,2
        spin_part_2: do sval=1,sd_nspin
          call op_lmix_n(mval,lval,a,dt,psi%wfn(:,sval,lval,mval),psi%wfn(:,sval,lval+1,mval), &
                         lmix_zp,lmix_zm,lmix_xp,lmix_xm,lmix_rhs,lmix_leqp_c,lmix_leqpf_c, &
                         lmix_leqm_c,lmix_leqmf_c,lmix_leqp_r,lmix_leqpf_r,lmix_leqm_r,lmix_leqmf_r)
          call op_lang_n(mval,lval,a,dt,psi%wfn(:,sval,lval,mval),psi%wfn(:,sval,lval+1,mval), &
                         lang_c1c,lang_c2c,lang_c3c,lang_c1r,lang_c2r,lang_c3r,lang_xl,lang_xp)
        end do spin_part_2
      end do process_second_wave_l
      !$omp end do
      !
      !$omp do schedule(dynamic,pt_chunk)
      process_first_wave_l: do lval=abs(mval),sd_lmax-1,2
        spin_part_1: do sval=1,sd_nspin
          call op_lmix_n(mval,lval,a,dt,psi%wfn(:,sval,lval,mval),psi%wfn(:,sval,lval+1,mval), &
                         lmix_zp,lmix_zm,lmix_xp,lmix_xm,lmix_rhs,lmix_leqp_c,lmix_leqpf_c, &
                         lmix_leqm_c,lmix_leqmf_c,lmix_leqp_r,lmix_leqpf_r,lmix_leqm_r,lmix_leqmf_r)
          call op_lang_n(mval,lval,a,dt,psi%wfn(:,sval,lval,mval),psi%wfn(:,sval,lval+1,mval), &
                         lang_c1c,lang_c2c,lang_c3c,lang_c1r,lang_c2r,lang_c3r,lang_xl,lang_xp)
        end do spin_part_1
      end do process_first_wave_l
      !$omp end do nowait
    end do process_m
    !$omp end parallel
    call TimerStop('Laser rev n (par L)')
  end subroutine pt_rev_laser_n_parallel_l
  !
  !  Propagator for the left wavefunction; use operator transpose
  !
  subroutine pt_fwd_laser_t(psi,a,dt)
    type(sd_wfn), intent(inout) :: psi ! Wavefunction we act upon
    real(xk), intent(in)        :: a   ! Vector potential along the Z axis
    complex(xk), intent(in)     :: dt  ! Timestep, in atomic units
    !
    if (sd_mmax-sd_mmin+1>sd_lmax+1 .and. .not. pt_force_par_l) then
      call pt_fwd_laser_t_parallel_m(psi,a,dt)
    else
      call pt_fwd_laser_t_parallel_l(psi,a,dt)
    end if
  end subroutine pt_fwd_laser_t
  !
  subroutine pt_fwd_laser_t_parallel_m(psi,a,dt)
    type(sd_wfn), intent(inout) :: psi ! Wavefunction we act upon
    real(xk), intent(in)        :: a   ! Vector potential along the Z axis
    complex(xk), intent(in)     :: dt  ! Timestep, in atomic units
    !
    integer(ik) :: lval, mval, sval
    ! Intel Fortran serializes on the memory allocator. To get reasonable
    ! parallel performance, we are forced to pull scratch allocation upstream
    ! from the worker routines. This sucks :(
    complex(rk) :: lang_c1c(sd_nradial), lang_c2c(sd_nradial), lang_c3c(sd_nradial)
    real(rk)    :: lang_c1r(sd_nradial), lang_c2r(sd_nradial), lang_c3r(sd_nradial)
    complex(rk) :: lang_xl(sd_nradial), lang_xp(sd_nradial)
    complex(rk) :: lmix_zp(sd_nradial,2), lmix_zm(sd_nradial,2)
    complex(rk) :: lmix_xp(sd_nradial,2), lmix_xm(sd_nradial,2)
    complex(rk) :: lmix_rhs(sd_nradial), lmix_tmp(sd_nradial)
    complex(rk) :: lmix_leqp_c(3,sd_nradial), lmix_leqpf_c(3,sd_nradial)
    complex(rk) :: lmix_leqm_c(3,sd_nradial), lmix_leqmf_c(3,sd_nradial)
    real(rk)    :: lmix_leqp_r(3,sd_nradial), lmix_leqpf_r(3,sd_nradial)
    real(rk)    :: lmix_leqm_r(3,sd_nradial), lmix_leqmf_r(3,sd_nradial)
    !
    if (a==0._rk) return ! Nothing to do!
    call TimerStart('Laser fwd t (par M)')
    !$omp parallel default(none) &
    !$omp& shared(sd_mmin,sd_mmax,sd_lmax,sd_nspin,a,dt,psi,pt_chunk) &
    !$omp& private(mval,lval,sval) &
    !$omp& private(lang_c1c,lang_c2c,lang_c3c,lang_c1r,lang_c2r,lang_c3r,lang_xl,lang_xp) &
    !$omp& private(lmix_zp,lmix_zm,lmix_xp,lmix_xm,lmix_rhs,lmix_tmp,lmix_leqp_c,lmix_leqpf_c) &
    !$omp& private(lmix_leqm_c,lmix_leqmf_c,lmix_leqp_r,lmix_leqpf_r,lmix_leqm_r,lmix_leqmf_r)
    !$omp do schedule(dynamic,pt_chunk)
    process_m: do mval=sd_mmin,sd_mmax
      process_first_wave_l: do lval=abs(mval),sd_lmax-1,2
        spin_part_1: do sval=1,sd_nspin
          ! Changing sign of dt converts op_lang_n() to op_lang_t(); however, since 
          ! we need time running backwards, we change the sign the second time.
          call op_lang_n(mval,lval,a,dt,psi%wfn(:,sval,lval,mval),psi%wfn(:,sval,lval+1,mval), &
                         lang_c1c,lang_c2c,lang_c3c,lang_c1r,lang_c2r,lang_c3r,lang_xl,lang_xp)
          call op_lmix_t(mval,lval,a,dt,psi%wfn(:,sval,lval,mval),psi%wfn(:,sval,lval+1,mval), &
                         lmix_zp,lmix_zm,lmix_xp,lmix_xm,lmix_rhs,lmix_tmp,lmix_leqp_c,lmix_leqpf_c, &
                         lmix_leqm_c,lmix_leqmf_c,lmix_leqp_r,lmix_leqpf_r,lmix_leqm_r,lmix_leqmf_r)
        end do spin_part_1
      end do process_first_wave_l
      !
      process_second_wave_l: do lval=abs(mval)+1,sd_lmax-1,2
        spin_part_2: do sval=1,sd_nspin
          ! Changing sign of dt converts op_lang_n() to op_lang_t()
          call op_lang_n(mval,lval,a,dt,psi%wfn(:,sval,lval,mval),psi%wfn(:,sval,lval+1,mval), &
                         lang_c1c,lang_c2c,lang_c3c,lang_c1r,lang_c2r,lang_c3r,lang_xl,lang_xp)
          call op_lmix_t(mval,lval,a,dt,psi%wfn(:,sval,lval,mval),psi%wfn(:,sval,lval+1,mval), &
                         lmix_zp,lmix_zm,lmix_xp,lmix_xm,lmix_rhs,lmix_tmp,lmix_leqp_c,lmix_leqpf_c, &
                         lmix_leqm_c,lmix_leqmf_c,lmix_leqp_r,lmix_leqpf_r,lmix_leqm_r,lmix_leqmf_r)
        end do spin_part_2
      end do process_second_wave_l
      !
      process_diagonal_term: do lval=abs(mval),sd_lmax
        spin_part_3: do sval=1,sd_nspin
          ! op_la2() is diagonal, and coincides with its own transpose
          ! Then, we flip the time step to get time running backwards
          call op_la2(a,-dt,psi%wfn(:,sval,lval,mval))
        end do spin_part_3
      end do process_diagonal_term
    end do process_m
    !$omp end do nowait
    !$omp end parallel
    call TimerStop('Laser fwd t (par M)')
  end subroutine pt_fwd_laser_t_parallel_m
  !
  subroutine pt_fwd_laser_t_parallel_l(psi,a,dt)
    type(sd_wfn), intent(inout) :: psi ! Wavefunction we act upon
    real(xk), intent(in)        :: a   ! Vector potential along the Z axis
    complex(xk), intent(in)     :: dt  ! Timestep, in atomic units
    !
    integer(ik) :: lval, mval, sval
    ! Intel Fortran serializes on the memory allocator. To get reasonable
    ! parallel performance, we are forced to pull scratch allocation upstream
    ! from the worker routines :(
    complex(rk) :: lang_c1c(sd_nradial), lang_c2c(sd_nradial), lang_c3c(sd_nradial)
    real(rk)    :: lang_c1r(sd_nradial), lang_c2r(sd_nradial), lang_c3r(sd_nradial)
    complex(rk) :: lang_xl(sd_nradial), lang_xp(sd_nradial)
    complex(rk) :: lmix_zp(sd_nradial,2), lmix_zm(sd_nradial,2)
    complex(rk) :: lmix_xp(sd_nradial,2), lmix_xm(sd_nradial,2)
    complex(rk) :: lmix_rhs(sd_nradial), lmix_tmp(sd_nradial)
    complex(rk) :: lmix_leqp_c(3,sd_nradial), lmix_leqpf_c(3,sd_nradial)
    complex(rk) :: lmix_leqm_c(3,sd_nradial), lmix_leqmf_c(3,sd_nradial)
    real(rk)    :: lmix_leqp_r(3,sd_nradial), lmix_leqpf_r(3,sd_nradial)
    real(rk)    :: lmix_leqm_r(3,sd_nradial), lmix_leqmf_r(3,sd_nradial)
    !
    if (a==0._rk) return ! Nothing to do!
    call TimerStart('Laser fwd t (par L)')
    !$omp parallel default(none) &
    !$omp& shared(sd_mmin,sd_mmax,sd_lmax,sd_nspin,a,dt,psi,pt_chunk) &
    !$omp& private(mval,lval,sval) &
    !$omp& private(lang_c1c,lang_c2c,lang_c3c,lang_c1r,lang_c2r,lang_c3r,lang_xl,lang_xp) &
    !$omp& private(lmix_zp,lmix_zm,lmix_xp,lmix_xm,lmix_rhs,lmix_tmp,lmix_leqp_c,lmix_leqpf_c) &
    !$omp& private(lmix_leqm_c,lmix_leqmf_c,lmix_leqp_r,lmix_leqpf_r,lmix_leqm_r,lmix_leqmf_r)
    process_m: do mval=sd_mmin,sd_mmax
      !$omp do schedule(dynamic,pt_chunk)
      process_first_wave_l: do lval=abs(mval),sd_lmax-1,2
        spin_part_1: do sval=1,sd_nspin
          ! Changing sign of dt converts op_lang_n() to op_lang_t(); however, since 
          ! we need time running backwards, we change the sign the second time.
          call op_lang_n(mval,lval,a,dt,psi%wfn(:,sval,lval,mval),psi%wfn(:,sval,lval+1,mval), &
                         lang_c1c,lang_c2c,lang_c3c,lang_c1r,lang_c2r,lang_c3r,lang_xl,lang_xp)
          call op_lmix_t(mval,lval,a,dt,psi%wfn(:,sval,lval,mval),psi%wfn(:,sval,lval+1,mval), &
                         lmix_zp,lmix_zm,lmix_xp,lmix_xm,lmix_rhs,lmix_tmp,lmix_leqp_c,lmix_leqpf_c, &
                         lmix_leqm_c,lmix_leqmf_c,lmix_leqp_r,lmix_leqpf_r,lmix_leqm_r,lmix_leqmf_r)
        end do spin_part_1
      end do process_first_wave_l
      !$omp end do
      !
      !$omp do schedule(dynamic,pt_chunk)
      process_second_wave_l: do lval=abs(mval)+1,sd_lmax-1,2
        spin_part_2: do sval=1,sd_nspin
          ! Changing sign of dt converts op_lang_n() to op_lang_t()
          call op_lang_n(mval,lval,a,dt,psi%wfn(:,sval,lval,mval),psi%wfn(:,sval,lval+1,mval), &
                         lang_c1c,lang_c2c,lang_c3c,lang_c1r,lang_c2r,lang_c3r,lang_xl,lang_xp)
          call op_lmix_t(mval,lval,a,dt,psi%wfn(:,sval,lval,mval),psi%wfn(:,sval,lval+1,mval), &
                         lmix_zp,lmix_zm,lmix_xp,lmix_xm,lmix_rhs,lmix_tmp,lmix_leqp_c,lmix_leqpf_c, &
                         lmix_leqm_c,lmix_leqmf_c,lmix_leqp_r,lmix_leqpf_r,lmix_leqm_r,lmix_leqmf_r)
        end do spin_part_2
      end do process_second_wave_l
      !$omp end do
      !
      !$omp do schedule(dynamic,pt_chunk)
      process_diagonal_term: do lval=abs(mval),sd_lmax
        spin_part_3: do sval=1,sd_nspin
          ! op_la2() is diagonal, and coincides with its own transpose
          ! Then, we flip the time step to get time running backwards
          call op_la2(a,-dt,psi%wfn(:,sval,lval,mval))
        end do spin_part_3
      end do process_diagonal_term
      !$omp end do nowait
    end do process_m
    !$omp end parallel
    call TimerStop('Laser fwd t (par L)')
  end subroutine pt_fwd_laser_t_parallel_l
  !
  subroutine pt_rev_laser_t(psi,a,dt)
    type(sd_wfn), intent(inout) :: psi ! Wavefunction we act upon
    real(xk), intent(in)        :: a   ! Vector potential along the Z axis
    complex(xk), intent(in)     :: dt  ! Timestep, in atomic units
    !
    if (sd_mmax-sd_mmin+1>sd_lmax+1 .and. .not. pt_force_par_l) then
      call pt_rev_laser_t_parallel_m(psi,a,dt)
    else
      call pt_rev_laser_t_parallel_l(psi,a,dt)
    end if
  end subroutine pt_rev_laser_t
  !
  subroutine pt_rev_laser_t_parallel_m(psi,a,dt)
    type(sd_wfn), intent(inout) :: psi ! Wavefunction we act upon
    real(xk), intent(in)        :: a   ! Vector potential along the Z axis
    complex(xk), intent(in)     :: dt  ! Timestep, in atomic units
    !
    integer(ik) :: lval, mval, sval
    ! Intel Fortran serializes on the memory allocator. To get reasonable
    ! parallel performance, we are forced to pull scratch allocation upstream
    ! from the worker routines :(
    complex(rk) :: lang_c1c(sd_nradial), lang_c2c(sd_nradial), lang_c3c(sd_nradial)
    real(rk)    :: lang_c1r(sd_nradial), lang_c2r(sd_nradial), lang_c3r(sd_nradial)
    complex(rk) :: lang_xl(sd_nradial), lang_xp(sd_nradial)
    complex(rk) :: lmix_zp(sd_nradial,2), lmix_zm(sd_nradial,2)
    complex(rk) :: lmix_xp(sd_nradial,2), lmix_xm(sd_nradial,2)
    complex(rk) :: lmix_rhs(sd_nradial), lmix_tmp(sd_nradial)
    complex(rk) :: lmix_leqp_c(3,sd_nradial), lmix_leqpf_c(3,sd_nradial)
    complex(rk) :: lmix_leqm_c(3,sd_nradial), lmix_leqmf_c(3,sd_nradial)
    real(rk)    :: lmix_leqp_r(3,sd_nradial), lmix_leqpf_r(3,sd_nradial)
    real(rk)    :: lmix_leqm_r(3,sd_nradial), lmix_leqmf_r(3,sd_nradial)
    !
    if (a==0._rk) return ! Nothing to do!
    call TimerStart('Laser rev t (par M)')
    !$omp parallel default(none) &
    !$omp& shared(sd_mmin,sd_mmax,sd_lmax,sd_nspin,a,dt,psi,pt_chunk) &
    !$omp& private(mval,lval,sval) &
    !$omp& private(lang_c1c,lang_c2c,lang_c3c,lang_c1r,lang_c2r,lang_c3r,lang_xl,lang_xp) &
    !$omp& private(lmix_zp,lmix_zm,lmix_xp,lmix_xm,lmix_rhs,lmix_tmp,lmix_leqp_c,lmix_leqpf_c) &
    !$omp& private(lmix_leqm_c,lmix_leqmf_c,lmix_leqp_r,lmix_leqpf_r,lmix_leqm_r,lmix_leqmf_r)
    !$omp do schedule(dynamic,pt_chunk)
    process_m: do mval=sd_mmin,sd_mmax
      process_diagonal_term: do lval=abs(mval),sd_lmax
        spin_part_3: do sval=1,sd_nspin
          ! op_la2() is diagonal, and coincides with its own transpose
          ! Changing sign of -dt makes sure time is running backwards for the left solution
          call op_la2(a,-dt,psi%wfn(:,sval,lval,mval))
        end do spin_part_3
      end do process_diagonal_term
      !
      process_second_wave_l: do lval=abs(mval)+1,sd_lmax-1,2
        spin_part_2: do sval=1,sd_nspin
          call op_lmix_t(mval,lval,a,dt,psi%wfn(:,sval,lval,mval),psi%wfn(:,sval,lval+1,mval), &
                         lmix_zp,lmix_zm,lmix_xp,lmix_xm,lmix_rhs,lmix_tmp,lmix_leqp_c,lmix_leqpf_c, &
                         lmix_leqm_c,lmix_leqmf_c,lmix_leqp_r,lmix_leqpf_r,lmix_leqm_r,lmix_leqmf_r)
          ! Changing sign of dt converts op_lang_n() to op_lang_t(). Then we change the sign
          ! back to plus since time is running backwards for the left solution
          call op_lang_n(mval,lval,a,dt,psi%wfn(:,sval,lval,mval),psi%wfn(:,sval,lval+1,mval), &
                         lang_c1c,lang_c2c,lang_c3c,lang_c1r,lang_c2r,lang_c3r,lang_xl,lang_xp)
        end do spin_part_2
      end do process_second_wave_l
      !
      process_first_wave_l: do lval=abs(mval),sd_lmax-1,2
        spin_part_1: do sval=1,sd_nspin
          call op_lmix_t(mval,lval,a,dt,psi%wfn(:,sval,lval,mval),psi%wfn(:,sval,lval+1,mval), &
                         lmix_zp,lmix_zm,lmix_xp,lmix_xm,lmix_rhs,lmix_tmp,lmix_leqp_c,lmix_leqpf_c, &
                         lmix_leqm_c,lmix_leqmf_c,lmix_leqp_r,lmix_leqpf_r,lmix_leqm_r,lmix_leqmf_r)
          ! Changing sign of dt converts op_lang_n() to op_lang_t(), and again since it is
          ! the left solution
          call op_lang_n(mval,lval,a,dt,psi%wfn(:,sval,lval,mval),psi%wfn(:,sval,lval+1,mval), &
                         lang_c1c,lang_c2c,lang_c3c,lang_c1r,lang_c2r,lang_c3r,lang_xl,lang_xp)
        end do spin_part_1
      end do process_first_wave_l
    end do process_m
    !$omp end do nowait
    !$omp end parallel
    call TimerStop('Laser rev t (par M)')
  end subroutine pt_rev_laser_t_parallel_m
  !
  subroutine pt_rev_laser_t_parallel_l(psi,a,dt)
    type(sd_wfn), intent(inout) :: psi ! Wavefunction we act upon
    real(xk), intent(in)        :: a   ! Vector potential along the Z axis
    complex(xk), intent(in)     :: dt  ! Timestep, in atomic units
    !
    integer(ik) :: lval, mval, sval
    ! Intel Fortran serializes on the memory allocator. To get reasonable
    ! parallel performance, we are forced to pull scratch allocation upstream
    ! from the worker routines :(
    complex(rk) :: lang_c1c(sd_nradial), lang_c2c(sd_nradial), lang_c3c(sd_nradial)
    real(rk)    :: lang_c1r(sd_nradial), lang_c2r(sd_nradial), lang_c3r(sd_nradial)
    complex(rk) :: lang_xl(sd_nradial), lang_xp(sd_nradial)
    complex(rk) :: lmix_zp(sd_nradial,2), lmix_zm(sd_nradial,2)
    complex(rk) :: lmix_xp(sd_nradial,2), lmix_xm(sd_nradial,2)
    complex(rk) :: lmix_rhs(sd_nradial), lmix_tmp(sd_nradial)
    complex(rk) :: lmix_leqp_c(3,sd_nradial), lmix_leqpf_c(3,sd_nradial)
    complex(rk) :: lmix_leqm_c(3,sd_nradial), lmix_leqmf_c(3,sd_nradial)
    real(rk)    :: lmix_leqp_r(3,sd_nradial), lmix_leqpf_r(3,sd_nradial)
    real(rk)    :: lmix_leqm_r(3,sd_nradial), lmix_leqmf_r(3,sd_nradial)
    !
    if (a==0._rk) return ! Nothing to do!
    call TimerStart('Laser rev t (par L)')
    !$omp parallel default(none) &
    !$omp& shared(sd_mmin,sd_mmax,sd_lmax,sd_nspin,a,dt,psi,pt_chunk) &
    !$omp& private(mval,lval,sval) &
    !$omp& private(lang_c1c,lang_c2c,lang_c3c,lang_c1r,lang_c2r,lang_c3r,lang_xl,lang_xp) &
    !$omp& private(lmix_zp,lmix_zm,lmix_xp,lmix_xm,lmix_rhs,lmix_tmp,lmix_leqp_c,lmix_leqpf_c) &
    !$omp& private(lmix_leqm_c,lmix_leqmf_c,lmix_leqp_r,lmix_leqpf_r,lmix_leqm_r,lmix_leqmf_r)
    process_m: do mval=sd_mmin,sd_mmax
      !$omp do schedule(dynamic,pt_chunk)
      process_diagonal_term: do lval=abs(mval),sd_lmax
        spin_part_3: do sval=1,sd_nspin
          ! op_la2() is diagonal, and coincides with its own transpose
          ! Changing sign of -dt makes sure time is running backwards for the left solution
          call op_la2(a,-dt,psi%wfn(:,sval,lval,mval))
        end do spin_part_3
      end do process_diagonal_term
      !$omp end do
      !
      !$omp do schedule(dynamic,pt_chunk)
      process_second_wave_l: do lval=abs(mval)+1,sd_lmax-1,2
        spin_part_2: do sval=1,sd_nspin
          call op_lmix_t(mval,lval,a,dt,psi%wfn(:,sval,lval,mval),psi%wfn(:,sval,lval+1,mval), &
                         lmix_zp,lmix_zm,lmix_xp,lmix_xm,lmix_rhs,lmix_tmp,lmix_leqp_c,lmix_leqpf_c, &
                         lmix_leqm_c,lmix_leqmf_c,lmix_leqp_r,lmix_leqpf_r,lmix_leqm_r,lmix_leqmf_r)
          ! Changing sign of dt converts op_lang_n() to op_lang_t(). Then we change the sign
          ! back to plus since time is running backwards for the left solution
          call op_lang_n(mval,lval,a,dt,psi%wfn(:,sval,lval,mval),psi%wfn(:,sval,lval+1,mval), &
                         lang_c1c,lang_c2c,lang_c3c,lang_c1r,lang_c2r,lang_c3r,lang_xl,lang_xp)
        end do spin_part_2
      end do process_second_wave_l
      !$omp end do
      !
      !$omp do schedule(dynamic,pt_chunk)
      process_first_wave_l: do lval=abs(mval),sd_lmax-1,2
        spin_part_1: do sval=1,sd_nspin
          call op_lmix_t(mval,lval,a,dt,psi%wfn(:,sval,lval,mval),psi%wfn(:,sval,lval+1,mval), &
                         lmix_zp,lmix_zm,lmix_xp,lmix_xm,lmix_rhs,lmix_tmp,lmix_leqp_c,lmix_leqpf_c, &
                         lmix_leqm_c,lmix_leqmf_c,lmix_leqp_r,lmix_leqpf_r,lmix_leqm_r,lmix_leqmf_r)
          ! Changing sign of dt converts op_lang_n() to op_lang_t(), and again since it is
          ! the left solution
          call op_lang_n(mval,lval,a,dt,psi%wfn(:,sval,lval,mval),psi%wfn(:,sval,lval+1,mval), &
                         lang_c1c,lang_c2c,lang_c3c,lang_c1r,lang_c2r,lang_c3r,lang_xl,lang_xp)
        end do spin_part_1
      end do process_first_wave_l
      !$omp end do nowait
    end do process_m
    !$omp end parallel
    call TimerStop('Laser rev t (par L)')
  end subroutine pt_rev_laser_t_parallel_l
  !
  !          /         0            (L+1) C(L+1,M) / r \
  !  Hlang = |                                         | ((I e)/(m c))*A
  !          \ -(L+1) C(L+1,M) / r         0           /
  !
  !  where C(L,M) = ((L**2-M**2)/(4*L**2-1))**0.5
  !
  !  The actual operator we apply is:
  !
  !  X = [1 + I dt Hlang]^-1 [1 - I dt Hlang] Psi
  !
  !  Since the Hlang operator is block-diagonal, multiplying through in the propagation
  !  operator gives:
  !
  !  X(L)   = (r**2+b**2)**-1 [ Psi(L)  *(r**2-b**2) + (2*b*r) Psi(L+1) ]
  !  X(L+1) = (r**2+b**2)**-1 [ Psi(L+1)*(r**2-b**2) - (2*b*r) Psi(L)   ]
  !  b      = (dt)*(e A/(m c))*(L+1)*C(L+1,M)
  !
  !  Note that we do not require a separate routine for op_lang_t(): since Hlang is
  !  Hermitian, changing the sign of dt converts op_lang_n() to op_lang_t().
  !
  subroutine op_lang_n(mval,lval,a,dt,psi_l,psi_p, c1c,c2c,c3c,c1r,c2r,c3r,xl,xp)
    integer(ik), intent(in)    :: mval     ! M value of both sub-blocks
    integer(ik), intent(in)    :: lval     ! L value of the first sub-block; second sub-block is L+1
    real(xk), intent(in)       :: a        ! Vector-potential along the Z axis
    complex(xk), intent(in)    :: dt       ! Time step
    complex(rk), intent(inout) :: psi_l(:) ! Radial part of the L block of the wavefunction
    complex(rk), intent(inout) :: psi_p(:) ! Radial part of the L+1 block of the wavefunction
    ! Intel Fortran spends an excessively long time allocating our scratch arrays,
    ! so that we must pass them from the calling routine to get reasonable performance.
    complex(rk), intent(inout) :: c1c(:), c2c(:), c3c(:)
    real(rk), intent(inout)    :: c1r(:), c2r(:), c3r(:)
    complex(rk), intent(inout) :: xl(:), xp(:)
    !
    complex(rk) :: fac_bc
    real(rk)    :: fac_br
  ! complex(rk) :: c1c(sd_nradial), c2c(sd_nradial), c3c(sd_nradial)
  ! real(rk)    :: c1r(sd_nradial), c2r(sd_nradial), c3r(sd_nradial)
  ! complex(rk) :: xl(sd_nradial), xp(sd_nradial)
    !
    if (aimag(dt)==0._xk) then
      fac_br = real(dt*a*(electron_charge/electron_mass)*(lval+1)*sd_laser_clm(lval+1,mval),kind=rk)
      c3r(:) = 1._rk/(sd_rtab(:sd_nradial)**2 + fac_br**2)
      c1r(:) = c3r(:)*(sd_rtab(:sd_nradial)**2 - fac_br**2)
      c2r(:) = c3r(:)*2._rk*fac_br*sd_rtab(:sd_nradial)
      !
      xl(:) = c1r(:)*psi_l(:) + c2r(:)*psi_p(:)
      xp(:) = c1r(:)*psi_p(:) - c2r(:)*psi_l(:)
    else
      fac_bc = cmplx(dt*a*(electron_charge/electron_mass)*(lval+1)*sd_laser_clm(lval+1,mval),kind=rk)
      c3c(:) = 1._rk/(sd_rtab(:sd_nradial)**2 + fac_bc**2)
      c1c(:) = c3c(:)*(sd_rtab(:sd_nradial)**2 - fac_bc**2)
      c2c(:) = c3c(:)*2._rk*fac_bc*sd_rtab(:sd_nradial)
      !
      xl(:) = c1c(:)*psi_l(:) + c2c(:)*psi_p(:)
      xp(:) = c1c(:)*psi_p(:) - c2c(:)*psi_l(:)
    end if
    !
    psi_l(:) = xl(:)
    psi_p(:) = xp(:)
  end subroutine op_lang_n
  !
  !         I e A          /     0          M1(L+1)^{-1} D1(L+1) \  .
  !  Hmix = ----- C(L+1,M) |                                     |  .
  !          m c           \ M1(L)^-1 D1(L)          0           /  .
  !
  !  where C(L,M) = ((L**2-M**2)/(4*L**2-1))**0.5, and M1, D1 are tridiagonal
  !  matrices defining the radial gradient operator.
  !
  !  The actual operator we apply is:
  !
  !  X = [1 + I (dt) Hmix]^-1 [1 - I (dt) Hmix] Psi
  !
  !  or, in matrix form:
  !
  !   / 1                -b*M1(L+1)^-1*D1(L+1) \ /  X(L)  \     /  Y(L)  \  .
  !   |                                        | |        |  =  |        |  .
  !   \-b M1(L)^-1*D1(L)  1                    / \ X(L+1) /     \ Y(L+1) /  .
  !
  !   /  Y(L)  \     / 1                 b*M1(L+1)^-1*D1(L+1) \ /  Z(L)  \ .
  !   |        |  =  |                                        | |        | .
  !   \ Y(L+1) /     \ b*M1(L)^-1*D1(L)  1                    / \ Z(L+1) / .
  !
  !  where Z is the wavefunction on the right, and we are seeking X.
  !
  !  There are two distinct implementations of this operator, depending on the
  !  symmetry of M1/D1. When M1(L+1)/=M1(L) [ie L=0] the operator inverse on the 
  !  left has to be computed with dense matrix algebra, using Gaussian elimination.
  !  For all other cases [ie L>=1] we can instead solve an auxiliary tridiagonal
  !  linear system:
  !
  !    (M1-b*D1) X'(L)   = (M1+b*D1)(Z(L)+Z(L+1))
  !    (M1+b*A1) X'(L+1) = (M1-b*D1)(Z(L)-Z(L+1))
  !    X(L)   = 0.5*(X'(L)+X'(L+1))
  !    X(L+1) = 0.5*(X'(L)-X'(L+1))
  !    b = (dt)*(e/(m c))*A*C(L+1,M)
  !
  subroutine op_lmix_n(mval,lval,a,dt,psi_l,psi_p, &
                       zp,zm,xp,xm,rhs,leqp_c,leqpf_c,leqm_c,leqmf_c,leqp_r,leqpf_r,leqm_r,leqmf_r)
    integer(ik), intent(in)    :: mval     ! M value of both sub-blocks
    integer(ik), intent(in)    :: lval     ! L value of the first sub-block; second sub-block is L+1
    real(xk), intent(in)       :: a        ! Vector-potential along the Z axis
    complex(xk), intent(in)    :: dt       ! Time step
    complex(rk), intent(inout) :: psi_l(:) ! Radial part of the L block of the wavefunction
    complex(rk), intent(inout) :: psi_p(:) ! Radial part of the L+1 block of the wavefunction
    ! Intel Fortran serializes on memory allocation; to get decent parallel scaling,
    ! we are forced to pass all temp arrays from upstream. This sucks.
    complex(rk), intent(inout) :: zp(:,:), zm(:,:)
    complex(rk), intent(inout) :: xp(:,:), xm(:,:)
    complex(rk), intent(inout) :: rhs(:)
    complex(rk), intent(inout) :: leqp_c(:,:), leqpf_c(:,:)
    complex(rk), intent(inout) :: leqm_c(:,:), leqmf_c(:,:)
    real(rk), intent(inout)    :: leqp_r(:,:), leqpf_r(:,:)
    real(rk), intent(inout)    :: leqm_r(:,:), leqmf_r(:,:)
    !
    complex(rk) :: fac_b
    !
    fac_b = cmplx(dt*(electron_charge/electron_mass)*a*sd_laser_clm(lval+1,mval),kind=rk)
    !
    !  The brute force version does not rely on gradient operator blocks being the same.
    !  This version should work for all L, but is too expensive to use generally.
    !
    if (lval==0) then
      select case (pt_mix_solver)
        case default
          write (out,"('propagator_tools%op_lmix_n: pt_mix_solver=',a,' is not recognized.')") trim(pt_mix_solver)
          stop 'propagator_tools%op_lmix_n - bad pt_mix_solver'
        case ('lapack')
          call op_lmix_bruteforce_n(lval,fac_b,psi_l,psi_p)
        case ('bi-CG')
          call op_lmix_bicg_n(lval,fac_b,psi_l,psi_p)
        case ('SM','default')
          call op_lmix_sherman_morrison_n(lval,fac_b,psi_l,psi_p, &
                                          zp,zm,xp,xm,leqp_c,leqpf_c,leqm_c,leqmf_c)
      end select
    else
      call op_lmix_mixunmix_n  (lval,fac_b,psi_l,psi_p, &
                                zp(:,1),zm(:,1),xp(:,1),xm(:,1),rhs,leqp_c,leqpf_c,leqm_c,leqmf_c,leqp_r,leqpf_r,leqm_r,leqmf_r)
    end if
  end subroutine op_lmix_n
  !
  subroutine op_lmix_bruteforce_n(lval,fac_b,psi_l,psi_p)
    integer(ik), intent(in)    :: lval     ! L value of the first sub-block; second sub-block is L+1
    complex(rk), intent(in)    :: fac_b    ! Amplitude for the coupling operator in the propagator
    complex(rk), intent(inout) :: psi_l(:) ! Radial part of the L block of the wavefunction
    complex(rk), intent(inout) :: psi_p(:) ! Radial part of the L+1 block of the wavefunction
    !
    complex(rk) :: yl(sd_nradial), yp(sd_nradial) ! Result of the direct-operator application
    complex(rk) :: tmp(sd_nradial)
    complex(rk) :: leq(2*sd_nradial,2*sd_nradial), rhs(2*sd_nradial,1)
    real(rk)    :: op(sd_nradial,sd_nradial)
    integer(ik) :: ind
    !
    call TimerStart('Laser lmix n brutefore')
    !
    !  Apply the direct part of the operator [1-I(dt)Hmix]
    !
    if (lval==0) then
      call m3d_multiply(sd_d1n_l0,psi_l,tmp)
      call m3d_solve(sd_m1nf_l0,tmp,yp)
    else
      call m3d_multiply(sd_d1n_lx,psi_l,tmp)
      call m3d_solve(sd_m1nf_lx,tmp,yp)
    end if
    yp(:) = psi_p(:) + fac_b*yp(:)
    call m3d_multiply(sd_d1n_lx,psi_p,tmp)  ! L+1 can never be zero, so this operator is fixed!
    call m3d_solve(sd_m1nf_lx,tmp,yl)
    yl(:) = psi_l(:) + fac_b*yl(:)
    !
    !  Prepare linear system matrix for the implicit part
    !
    leq(1:sd_nradial,   1:sd_nradial   ) = 0
    leq(  sd_nradial+1:,  sd_nradial+1:) = 0
    fill_diagonal: do ind=1,2*sd_nradial
      leq(ind,ind) = 1
    end do fill_diagonal
    !
    call sd_expand_implicit_operator(sd_d1n_lx,sd_m1nf_lx,op)
    leq(1:sd_nradial,sd_nradial+1:) = -fac_b * op
    if (lval==0) then
      call sd_expand_implicit_operator(sd_d1n_l0,sd_m1nf_l0,op)
      leq(sd_nradial+1:,1:sd_nradial) = -fac_b * op
    else
      leq(sd_nradial+1:,1:sd_nradial) = leq(1:sd_nradial,sd_nradial+1:)
    end if
    !
    !  Solve for the operator inverse
    !
    rhs(1:sd_nradial   ,1) = yl
    rhs(  sd_nradial+1:,1) = yp
    call lapack_gesv(leq,rhs)
    !
    !  Copy out the solutions
    !
    psi_l = rhs(1:sd_nradial   ,1)
    psi_p = rhs(  sd_nradial+1:,1)
    call TimerStop('Laser lmix n brutefore')
  end subroutine op_lmix_bruteforce_n
  !
  subroutine op_lmix_bicg_n(lval,fac_b,psi_l,psi_p)
    integer(ik), intent(in)    :: lval     ! L value of the first sub-block; second sub-block is L+1
    complex(rk), intent(in)    :: fac_b    ! Amplitude for the coupling operator in the propagator
    complex(rk), intent(inout) :: psi_l(:) ! Radial part of the L block of the wavefunction
    complex(rk), intent(inout) :: psi_p(:) ! Radial part of the L+1 block of the wavefunction
    !
    complex(rk) :: yl(sd_nradial), yp(sd_nradial) ! Result of the direct-operator application
    complex(rk) :: tmp(sd_nradial)
    complex(rk) :: rhs(2*sd_nradial), x(2*sd_nradial)
    logical     :: fail
    !
    ! call TimerStart('Laser lmix n bicg')
    !
    !  Apply the direct part of the operator [1-I(dt)Hmix]
    !
    if (lval==0) then
      call m3d_multiply(sd_d1n_l0,psi_l,tmp)
      call m3d_solve(sd_m1nf_l0,tmp,yp)
    else
      call m3d_multiply(sd_d1n_lx,psi_l,tmp)
      call m3d_solve(sd_m1nf_lx,tmp,yp)
    end if
    if (pt_sense_real .and. aimag(fac_b)==0._rk) then
      yp(:) = psi_p(:) + real(fac_b,kind=rk)*yp(:)
    else
      yp(:) = psi_p(:) + fac_b*yp(:)
    end if
    call m3d_multiply(sd_d1n_lx,psi_p,tmp)  ! L+1 can never be zero, so this operator is fixed!
    call m3d_solve(sd_m1nf_lx,tmp,yl)
    if (pt_sense_real .and. aimag(fac_b)==0._rk) then
      yl(:) = psi_l(:) + real(fac_b,kind=rk)*yl(:)
    else
      yl(:) = psi_l(:) + fac_b*yl(:)
    end if
    !
    !  Solve the inverse part of the propagator
    !
    rhs(1:sd_nradial   ) = yl
    rhs(  sd_nradial+1:) = yp
    !
    call bicg_solve(rhs,x,(/lval,0/),(/fac_b/),mix_inverse,fail)
    if (fail) then
      call op_lmix_bruteforce_n(lval,fac_b,psi_l,psi_p)
      return
    end if
    !
    !  Copy out the solutions
    !
    psi_l = x(1:sd_nradial   )
    psi_p = x(  sd_nradial+1:)
    ! call TimerStop('Laser lmix n bicg')
  end subroutine op_lmix_bicg_n
  !
  !  The routine below is likely incomprehensible without reading the notes in:
  !  "hgm-revised-non-iterative-mixing-term.pdf" or the section 2.3.3 in the
  !  reference paper.
  !
  subroutine op_lmix_sherman_morrison_n(lval,fac_b,psi_l,psi_p, &
                                        zp,zm,xp,xm,mpd,mpd_f,mmd,mmd_f)
    integer(ik), intent(in)    :: lval     ! L value of the first sub-block; second sub-block is L+1
    complex(rk), intent(in)    :: fac_b    ! Amplitude for the coupling operator in the propagator
    complex(rk), intent(inout) :: psi_l(:) ! Radial part of the L block of the wavefunction
    complex(rk), intent(inout) :: psi_p(:) ! Radial part of the L+1 block of the wavefunction
    !
    complex(rk), intent(inout) :: zp(:,:), zm(:,:)      ! I'll be solving multiple linear systems
    complex(rk), intent(inout) :: xp(:,:), xm(:,:)
    complex(rk), intent(inout) :: mpd(:,:), mpd_f(:,:)  ! M+D and the corresponding factors
    complex(rk), intent(inout) :: mmd(:,:), mmd_f(:,:)  ! M-D and the corresponding factors
    ! Local temporaries; if compiler serializes on memory allocation, consider moving upstream
    ! complex(rk) :: zp(sd_nradial,2), zm(sd_nradial,2)      ! I'll be solving multiple linear systems
    ! complex(rk) :: xp(sd_nradial,2), xm(sd_nradial,2)
    ! complex(rk) :: mpd(3,sd_nradial), mpd_f(3,sd_nradial)  ! M+D and the corresponding factors
    ! complex(rk) :: mmd(3,sd_nradial), mmd_f(3,sd_nradial)  ! M-D and the corresponding factors
    !
    complex(rk) :: m1l0(2)        ! Non-zero elements of the first row of sd_m1n_l0
    complex(rk) :: m1lx(2)        !  ... ditto for sd_m1n_lx
    complex(rk) :: d1l0(2)        !  ... ditto for sd_d1n_l0
    complex(rk) :: d1lx(2)        !  ... ditto for sd_d1n_lx
    complex(rk) :: npe(2), nme(2) ! Non-zero elements of N+E and N-E
    complex(rk) :: smfact         ! Correction amplitude
    !
    !  First of all, make sure Sherman-Morrison formula applies to this term.
    !  At most we can handle a difference in the first row of the M_1 and D_1
    !  matrices; anything else should go to the iterative propagator.
    !  Identical M_1/D_1 matrices are OK, but a waste of effort.
    !
    if (lval==0) then
      if (any(sd_m1n_l0(:,2:)/=sd_m1n_lx(:,2:)) .or. sd_m1n_l0(2,1)/=sd_m1n_lx(2,1) .or. &
          any(sd_d1n_l0(:,2:)/=sd_d1n_lx(:,2:)) .or. sd_d1n_l0(2,1)/=sd_d1n_lx(2,1)) then
        stop 'propagator_tools%op_lmix_sherman_morrison_n - d1n/m1n is not in the expected form'
      end if
      ! Construct the average operator, and the two-point difference operator
      ! Don't try to be clever: this routine is usually not on the critical path
      mpd  = sd_m1n_l0 + fac_b*sd_d1n_l0
      mmd  = sd_m1n_l0 - fac_b*sd_d1n_l0
      ! Extract the first row of all matrices, and pre-multiply them by the factors we'll need
      m1l0 = sd_m1n_l0((/1,3/),1) * 0.5_rk
      m1lx = sd_m1n_lx((/1,3/),1) * 0.5_rk
      d1l0 = sd_d1n_l0((/1,3/),1) * 0.5_rk * fac_b
      d1lx = sd_d1n_lx((/1,3/),1) * 0.5_rk * fac_b
      ! Construct the difference operator
      npe  = (m1lx-m1l0) + (d1lx-d1l0)
      nme  = (m1lx-m1l0) - (d1lx-d1l0)
      ! Adjust the average operator
      mpd((/1,3/),1) = mpd((/1,3/),1) + npe
      mmd((/1,3/),1) = mmd((/1,3/),1) + nme
    else 
      ! lval>0, operators are the same, and the difference terms vanish
      mpd  = sd_m1n_lx + fac_b*sd_d1n_lx
      mmd  = sd_m1n_lx - fac_b*sd_d1n_lx
      nme  = 0
      npe  = 0
    end if
    !
    !  "Mix" the input wavefunction
    !
    xp(:,1) = psi_l(:) + psi_p(:)
    xm(:,1) = psi_l(:) - psi_p(:)
    !
    !  Construct the right-hand side:
    !
    !   zp = (M+D) . xp + (N-E) . xm
    !   zm = (M-D) . xm + (N+E) . xp
    !
    call m3d_multiply(mpd,xp(:,1),zp(:,1))
    call m3d_multiply(mmd,xm(:,1),zm(:,1))
    zp(1,1) = zp(1,1) + sum(nme*xm(1:2,1))
    zm(1,1) = zm(1,1) + sum(npe*xp(1:2,1))
    !
    !  Now construct the "unperturbed" linear system matrix for Sherman-Morrison.
    !  If you pay attention, it turns out to coinside with the laser-coupling 
    !  operator for L=lval.
    !
    mpd((/1,3/),1) = mpd((/1,3/),1) - npe
    mmd((/1,3/),1) = mmd((/1,3/),1) - nme
    call m3d_decompose(mpd,mpd_f)
    call m3d_decompose(mmd,mmd_f)
    !
    !  Construct the auxiliary right-hand side for the Sherman-Morrison correction
    !  This is vector "u"
    !
    zp(:,2) = 0 ; zp(1,2) = 1
    zm(:,2) = 0 ; zm(1,2) = 1
    !
    !  Solve the unperturbed linear system and the correction equation
    !
    call m3d_solve(mmd_f,zp,xp)
    call m3d_solve(mpd_f,zm,xm)
    !
    !  Correction factor: (v.y)/(1+(v.z))
    !  Vector v has only four non-zero elements
    !  y is the unperturbed solution; z is the correction-vector solution
    !
    smfact = (  sum(nme*xp(1:2,1))+sum(npe*xm(1:2,1))) &
           / (1+sum(nme*xp(1:2,2))+sum(npe*xm(1:2,2)))
    !
    !  Apply the correction
    !
    xp(:,1) = xp(:,1) - smfact*xp(:,2)
    xm(:,1) = xm(:,1) - smfact*xm(:,2)
    !
    !  Unmix, and we are done
    !
    psi_l(:) = 0.5_rk*(xp(:,1)+xm(:,1))
    psi_p(:) = 0.5_rk*(xp(:,1)-xm(:,1))
  end subroutine op_lmix_sherman_morrison_n
  !
  !  Linear matrix for which we need an inverse in op_lmix_bicg_n()
  !
  subroutine mix_inverse(op,iblob,cblob,vec,av)
    character(len=*), intent(in) :: op       ! 'N' or 'T'
    integer(ik), intent(in)      :: iblob(:) ! An opaque binary blob, passed along from the cg_solve() call
    complex(rk), intent(in)      :: cblob(:) ! An opaque binary blob, passed along from the cg_solve() call
    complex(rk), intent(in)      :: vec(:)   ! Vector
    complex(rk), intent(out)     :: av (:)   ! Matrix-vector product: av = a^op . vec
    !
    integer(ik)      :: lval
    character(len=1) :: mtype
    complex(rk)      :: fac_b
    complex(rk)      :: tmp(sd_nradial)
    !
    if (size(iblob)/=2 .or. size(cblob)/=1) stop 'propagator_tools%mix_inverse - bad blobs'
    if (size(vec)/=2*sd_nradial .or. size(av)/=2*sd_nradial) stop 'propagator_tools%mix_inverse - bad vectors'
    if (iblob(1)<0 .or. iblob(1)>sd_lmax-1) stop 'propagator_tools%mix_inverse - bad lval'
    if (iblob(2)<0 .or. iblob(2)>1) stop 'propagator_tools%mix_inverse - bad matrix type'
    !
    lval  = iblob(1) 
    mtype = 'N'
    if (iblob(2)==1) mtype = 'T'
    fac_b = cblob(1)
    !
    if (mtype==op) then
      ! Multiply by normal matrix
      if (lval==0) then
        call m3d_multiply(sd_d1n_l0,vec(1:sd_nradial),tmp)
        call m3d_solve(sd_m1nf_l0,tmp,av(sd_nradial+1:))
      else
        call m3d_multiply(sd_d1n_lx,vec(1:sd_nradial),tmp)
        call m3d_solve(sd_m1nf_lx,tmp,av(sd_nradial+1:))
      end if
      call m3d_multiply(sd_d1n_lx,vec(sd_nradial+1:),tmp)
      call m3d_solve(sd_m1nf_lx,tmp,av(1:sd_nradial))
    else 
      ! Multiply by matrix transpose
      if (lval==0) then
        call m3d_solve(sd_m1tf_l0,vec(sd_nradial+1:),tmp)
        call m3d_multiply(sd_d1t_l0,tmp,av(1:sd_nradial))
      else
        call m3d_solve(sd_m1tf_lx,vec(sd_nradial+1:),tmp)
        call m3d_multiply(sd_d1t_lx,tmp,av(1:sd_nradial))
      end if
      call m3d_solve(sd_m1tf_lx,vec(1:sd_nradial),tmp)
      call m3d_multiply(sd_d1t_lx,tmp,av(sd_nradial+1:))
    end if
    if (pt_sense_real .and. aimag(fac_b)==0._rk) then
      av = -real(fac_b,kind=rk) * av + vec
    else
      av = -fac_b * av + vec
    end if
  end subroutine mix_inverse
  !
  subroutine op_lmix_mixunmix_n(lval,fac_b,psi_l,psi_p, &
                                zp,zm,xp,xm,rhs,leqp_c,leqpf_c,leqm_c,leqmf_c,leqp_r,leqpf_r,leqm_r,leqmf_r)
    integer(ik), intent(in)    :: lval     ! L value of the first sub-block; second sub-block is L+1
    complex(rk), intent(in)    :: fac_b    ! Amplitude for the coupling operator in the propagator
    complex(rk), intent(inout) :: psi_l(:) ! Radial part of the L block of the wavefunction
    complex(rk), intent(inout) :: psi_p(:) ! Radial part of the L+1 block of the wavefunction
    ! Intel Fortran serializes on memory allocation; to get decent parallel scaling,
    ! we are forced to pass all temp arrays from upstream. This sucks.
    complex(rk), intent(inout) :: zp(:), zm(:)
    complex(rk), intent(inout) :: xp(:), xm(:)
    complex(rk), intent(inout) :: rhs(:)
    complex(rk), intent(inout) :: leqp_c(:,:), leqpf_c(:,:)
    complex(rk), intent(inout) :: leqm_c(:,:), leqmf_c(:,:)
    real(rk), intent(inout)    :: leqp_r(:,:), leqpf_r(:,:)
    real(rk), intent(inout)    :: leqm_r(:,:), leqmf_r(:,:)
    !
  ! complex(rk) :: zp(sd_nradial), zm(sd_nradial)
  ! complex(rk) :: xp(sd_nradial), xm(sd_nradial)
  ! complex(rk) :: rhs(sd_nradial)
  ! complex(rk) :: leqp_c(3,sd_nradial), leqpf_c(3,sd_nradial)
  ! complex(rk) :: leqm_c(3,sd_nradial), leqmf_c(3,sd_nradial)
  ! real(rk)    :: leqp_r(3,sd_nradial), leqpf_r(3,sd_nradial)
  ! real(rk)    :: leqm_r(3,sd_nradial), leqmf_r(3,sd_nradial)
    !
    if (lval==0) stop 'propagator_tools%op_lmix_mixunmix_n - lval=0 is not valid'
    !
    zp(:) = psi_l(:) + psi_p(:)
    zm(:) = psi_l(:) - psi_p(:)
    !
    !  Construct and factor linear system matrices.
    !
    if (pt_sense_real .and. aimag(fac_b)==0._rk) then
      leqp_r = sd_m1n_lx + real(fac_b,kind=rk)*sd_d1n_lx
      leqm_r = sd_m1n_lx - real(fac_b,kind=rk)*sd_d1n_lx
      call m3d_decompose(leqp_r,leqpf_r)
      call m3d_decompose(leqm_r,leqmf_r)
      !
      !  Solve for the mixed propagator
      !
      call m3d_multiply(leqp_r,zp,rhs)
      call m3d_solve(leqmf_r,rhs,xp)
      call m3d_multiply(leqm_r,zm,rhs)
      call m3d_solve(leqpf_r,rhs,xm)
    else
      leqp_c = sd_m1n_lx + fac_b*sd_d1n_lx
      leqm_c = sd_m1n_lx - fac_b*sd_d1n_lx
      call m3d_decompose(leqp_c,leqpf_c)
      call m3d_decompose(leqm_c,leqmf_c)
      !
      !  Solve for the mixed propagator
      !
      call m3d_multiply(leqp_c,zp,rhs)
      call m3d_solve(leqmf_c,rhs,xp)
      call m3d_multiply(leqm_c,zm,rhs)
      call m3d_solve(leqpf_c,rhs,xm)
    end if
    !
    !  Unmix, and we are done
    !
    psi_l(:) = 0.5_rk*(xp+xm)
    psi_p(:) = 0.5_rk*(xp-xm)
  end subroutine op_lmix_mixunmix_n
  !
  !  In matrix form, our operator is now:
  !
  !   / 1                         a*D1(L)^T*M1(L)^T^-1 \ /  X(L)  \     /  Y(L)  \  .
  !   |                                                | |        |  =  |        |  .
  !   \ a*D1(L+1)^T*M1(L+1)^T^-1  1                    / \ X(L+1) /     \ Y(L+1) /  .
  !
  !   /  Y(L)  \     /  1                         -a*D1(L)^T*M1(L)^T^-1 \ /  Z(L)  \ .
  !   |        |  =  |                                                  | |        | .
  !   \ Y(L+1) /     \ -a*D1(L+1)^T*M1(L+1)^T^-1  1                     / \ Z(L+1) / .
  !
  !  where Z is the wavefunction on the right, and we are seeking X.
  !
  !  For the symmetric case, we now have:
  !
  !    (M1^T+a*D1^T) Y'(L)   = (1-a*D1^T*M1^-T) (Z(L)+Z(L+1))
  !    (M1^T-a*D1^T) Y'(L+1) = (1+a*D1^T*M1^-T) (Z(L)-Z(L+1))
  !    X(L)   = 0.5*M1^T*(Y'(L)+Y'(L+1))
  !    X(L+1) = 0.5*M1^T*(Y'(L)-Y'(L+1))
  !
  subroutine op_lmix_t(mval,lval,a,dt,psi_l,psi_p, &
                       zp,zm,xp,xm,rhs,tmp,leqp_c,leqpf_c,leqm_c,leqmf_c,leqp_r,leqpf_r,leqm_r,leqmf_r)
    integer(ik), intent(in)    :: mval     ! M value of both sub-blocks
    integer(ik), intent(in)    :: lval     ! L value of the first sub-block; second sub-block is L+1
    real(xk), intent(in)       :: a        ! Vector-potential along the Z axis
    complex(xk), intent(in)    :: dt       ! Time step
    complex(rk), intent(inout) :: psi_l(:) ! Radial part of the L block of the wavefunction
    complex(rk), intent(inout) :: psi_p(:) ! Radial part of the L+1 block of the wavefunction
    ! Intel Fortran serializes on memory allocation; as the result, we must pass
    ! local temporaries in op_lmix_mixunmix_t() from upstream to get reasonable
    ! parallel performance. This sucks :O(
    complex(rk), intent(inout) :: zp(:,:), zm(:,:)
    complex(rk), intent(inout) :: xp(:,:), xm(:,:)
    complex(rk), intent(inout) :: rhs(:), tmp(:)
    complex(rk), intent(inout) :: leqp_c(:,:), leqpf_c(:,:)
    complex(rk), intent(inout) :: leqm_c(:,:), leqmf_c(:,:)
    real(rk), intent(inout)    :: leqp_r(:,:), leqpf_r(:,:)
    real(rk), intent(inout)    :: leqm_r(:,:), leqmf_r(:,:)
    !
    complex(rk) :: fac_b
    !
    !  Changing sign of dt makes time run backwards for the left eigenfunction
    !
    fac_b = cmplx(-dt*(electron_charge/electron_mass)*a*sd_laser_clm(lval+1,mval),kind=rk)
    !
    !  The brute force version does not rely on gradient operator blocks being the same.
    !  This version should work for all L, but is too expensive to use generally.
    !
    if (lval==0) then
      select case (pt_mix_solver)
        case default
          write (out,"('propagator_tools%op_lmix_t: pt_mix_solver=',a,' is not recognized.')") trim(pt_mix_solver)
          stop 'propagator_tools%op_lmix_t - bad pt_mix_solver'
        case ('lapack')
          call op_lmix_bruteforce_t(lval,fac_b,psi_l,psi_p)
        case ('bi-CG')
          call op_lmix_bicg_t(lval,fac_b,psi_l,psi_p)
        case ('SM','default')
          call op_lmix_sherman_morrison_t(lval,fac_b,psi_l,psi_p, &
                                          zp,zm,xp,xm,leqp_c,leqpf_c,leqm_c,leqmf_c)
      end select
    else
      call op_lmix_mixunmix_t(lval,fac_b,psi_l,psi_p, &
                              zp(:,1),zm(:,1),xp(:,1),xm(:,1),rhs,tmp,leqp_c,leqpf_c,leqm_c,leqmf_c,leqp_r,leqpf_r,leqm_r,leqmf_r)
    end if
  end subroutine op_lmix_t
  !
  subroutine op_lmix_bruteforce_t(lval,fac_b,psi_l,psi_p)
    integer(ik), intent(in)    :: lval     ! L value of the first sub-block; second sub-block is L+1
    complex(rk), intent(in)    :: fac_b    ! Amplitude for the coupling operator in the propagator
    complex(rk), intent(inout) :: psi_l(:) ! Radial part of the L block of the wavefunction
    complex(rk), intent(inout) :: psi_p(:) ! Radial part of the L+1 block of the wavefunction
    !
    complex(rk) :: yl(sd_nradial), yp(sd_nradial) ! Result of the direct-operator application
    complex(rk) :: tmp(sd_nradial)
    complex(rk) :: leq(2*sd_nradial,2*sd_nradial), rhs(2*sd_nradial,1)
    real(rk)    :: op(sd_nradial,sd_nradial)
    integer(ik) :: ind
    !
    call TimerStart('Laser lmix t bruteforce')
    !
    !  Apply the direct part of the operator [1-I(dt)Hmix]
    !
    if (lval==0) then
      call m3d_solve(sd_m1tf_l0,psi_p,tmp)
      call m3d_multiply(sd_d1t_l0,tmp,yl)
    else
      call m3d_solve(sd_m1tf_lx,psi_p,tmp)
      call m3d_multiply(sd_d1t_lx,tmp,yl)
    end if
    yl(:) = psi_l(:) + fac_b*yl(:)
    call m3d_solve(sd_m1tf_lx,psi_l,tmp)    ! L+1 can never be zero, so this operator is fixed!
    call m3d_multiply(sd_d1t_lx,tmp,yp)
    yp(:) = psi_p(:) + fac_b*yp(:)
    !
    !  Prepare linear system matrix for the implicit part
    !  It is less work to prepare direct matrix, then take the transpose
    !
    leq(1:sd_nradial,   1:sd_nradial   ) = 0
    leq(  sd_nradial+1:,  sd_nradial+1:) = 0
    fill_diagonal: do ind=1,2*sd_nradial
      leq(ind,ind) = 1
    end do fill_diagonal
    !
    call sd_expand_implicit_operator(sd_d1n_lx,sd_m1nf_lx,op)
    leq(1:sd_nradial,sd_nradial+1:) = -fac_b * op
    if (lval==0) then
      call sd_expand_implicit_operator(sd_d1n_l0,sd_m1nf_l0,op)
      leq(sd_nradial+1:,1:sd_nradial) = -fac_b * op
    else
      leq(sd_nradial+1:,1:sd_nradial) = leq(1:sd_nradial,sd_nradial+1:)
    end if
    leq = transpose(leq)
    !
    !  Solve for the operator inverse
    !
    rhs(1:sd_nradial   ,1) = yl
    rhs(  sd_nradial+1:,1) = yp
    call lapack_gesv(leq,rhs)
    !
    !  Copy out the solutions
    !
    psi_l = rhs(1:sd_nradial   ,1)
    psi_p = rhs(  sd_nradial+1:,1)
    call TimerStop('Laser lmix t bruteforce')
  end subroutine op_lmix_bruteforce_t
  !
  subroutine op_lmix_bicg_t(lval,fac_b,psi_l,psi_p)
    integer(ik), intent(in)    :: lval     ! L value of the first sub-block; second sub-block is L+1
    complex(rk), intent(in)    :: fac_b    ! Amplitude for the coupling operator in the propagator
    complex(rk), intent(inout) :: psi_l(:) ! Radial part of the L block of the wavefunction
    complex(rk), intent(inout) :: psi_p(:) ! Radial part of the L+1 block of the wavefunction
    !
    complex(rk) :: yl(sd_nradial), yp(sd_nradial) ! Result of the direct-operator application
    complex(rk) :: tmp(sd_nradial)
    complex(rk) :: rhs(2*sd_nradial), x(2*sd_nradial)
    logical     :: fail
    !
    ! call TimerStart('Laser lmix t bicg')
    !
    !  Apply the direct part of the operator [1-I(dt)Hmix]
    !
    if (lval==0) then
      call m3d_solve(sd_m1tf_l0,psi_p,tmp)
      call m3d_multiply(sd_d1t_l0,tmp,yl)
    else
      call m3d_solve(sd_m1tf_lx,psi_p,tmp)
      call m3d_multiply(sd_d1t_lx,tmp,yl)
    end if
    if (pt_sense_real .and. aimag(fac_b)==0._rk) then
      yl(:) = psi_l(:) + real(fac_b,kind=rk)*yl(:)
    else
      yl(:) = psi_l(:) + fac_b*yl(:)
    end if
    call m3d_solve(sd_m1tf_lx,psi_l,tmp)    ! L+1 can never be zero, so this operator is fixed!
    call m3d_multiply(sd_d1t_lx,tmp,yp)
    if (pt_sense_real .and. aimag(fac_b)==0._rk) then
      yp(:) = psi_p(:) + real(fac_b,kind=rk)*yp(:)
    else
      yp(:) = psi_p(:) + fac_b*yp(:)
    end if
    !
    rhs(1:sd_nradial   ) = yl
    rhs(  sd_nradial+1:) = yp
    !
    call bicg_solve(rhs,x,(/lval,1/),(/fac_b/),mix_inverse,fail)
    if (fail) then
      call op_lmix_bruteforce_t(lval,fac_b,psi_l,psi_p)
      return
    end if
    !
    !  Copy out the solutions
    !
    psi_l = x(1:sd_nradial   )
    psi_p = x(  sd_nradial+1:)
    ! call TimerStop('Laser lmix t bicg')
  end subroutine op_lmix_bicg_t
  !
  !  The routine below is likely incomprehensible without reading the notes in:
  !  "hgm-revised-non-iterative-mixing-term.pdf" or the section 2.3.3 in the
  !  reference paper.
  !
  subroutine op_lmix_sherman_morrison_t(lval,fac_b,psi_l,psi_p, &
                                        zp,zm,xp,xm,mpd,mpd_f,mmd,mmd_f)
    integer(ik), intent(in)    :: lval     ! L value of the first sub-block; second sub-block is L+1
    complex(rk), intent(in)    :: fac_b    ! Amplitude for the coupling operator in the propagator
    complex(rk), intent(inout) :: psi_l(:) ! Radial part of the L block of the wavefunction
    complex(rk), intent(inout) :: psi_p(:) ! Radial part of the L+1 block of the wavefunction
    !
    complex(rk), intent(inout) :: zp(:,:), zm(:,:)      ! I'll be solving multiple linear systems
    complex(rk), intent(inout) :: xp(:,:), xm(:,:)
    complex(rk), intent(inout) :: mpd(:,:), mpd_f(:,:)  ! M+D and the corresponding factors
    complex(rk), intent(inout) :: mmd(:,:), mmd_f(:,:)  ! M-D and the corresponding factors
    ! Local temporaries; if compiler serializes on memory allocation, consider moving upstream
    ! complex(rk) :: zp(sd_nradial,2), zm(sd_nradial,2)      ! I'll be solving multiple linear systems
    ! complex(rk) :: xp(sd_nradial,2), xm(sd_nradial,2)
    ! complex(rk) :: mpd(3,sd_nradial), mpd_f(3,sd_nradial)  ! M+D and the corresponding factors
    ! complex(rk) :: mmd(3,sd_nradial), mmd_f(3,sd_nradial)  ! M-D and the corresponding factors
    !
    complex(rk) :: m1l0(2)        ! Non-zero elements of the first column of sd_m1t_l0
    complex(rk) :: m1lx(2)        !  ... ditto for sd_m1t_lx
    complex(rk) :: d1l0(2)        !  ... ditto for sd_d1t_l0
    complex(rk) :: d1lx(2)        !  ... ditto for sd_d1t_lx
    complex(rk) :: npe(2), nme(2) ! Non-zero elements of N+E and N-E
    complex(rk) :: smfact         ! Correction amplitude
    !
    !  First of all, make sure Sherman-Morrison formula applies to this term.
    !  At most we can handle a difference in the first column of the M_1^T and D_1^T
    !  matrices; anything else should go to the iterative propagator.
    !  Identical M_1^T/D_1^T matrices are OK, but a waste of effort.
    !
    if (lval==0) then
      if (any(sd_m1t_l0(:,2:)/=sd_m1t_lx(:,2:)) .or. sd_m1t_l0(3,1)/=sd_m1t_lx(3,1) .or. &
          any(sd_d1t_l0(:,2:)/=sd_d1t_lx(:,2:)) .or. sd_d1t_l0(3,1)/=sd_d1t_lx(3,1)) then
        stop 'propagator_tools%op_lmix_sherman_morrison_t - d1t/m1t is not in the expected form'
      end if
      ! Construct the average operator, and the two-point difference operator
      ! Don't try to be clever: this routine is usually not on the critical path
      mpd  = sd_m1t_l0 + fac_b*sd_d1t_l0
      mmd  = sd_m1t_l0 - fac_b*sd_d1t_l0
      ! Extract the first column of all matrices, and pre-multiply them by the factors we'll need
      m1l0 = sd_m1t_l0((/1,2/),1) * 0.5_rk
      m1lx = sd_m1t_lx((/1,2/),1) * 0.5_rk
      d1l0 = sd_d1t_l0((/1,2/),1) * 0.5_rk * fac_b
      d1lx = sd_d1t_lx((/1,2/),1) * 0.5_rk * fac_b
      ! Construct the difference operator
      npe  = (m1lx-m1l0) + (d1lx-d1l0)
      nme  = (m1lx-m1l0) - (d1lx-d1l0)
      ! Adjust the average operator
      mpd((/1,2/),1) = mpd((/1,2/),1) + npe
      mmd((/1,2/),1) = mmd((/1,2/),1) + nme
    else 
      ! lval>0, operators are the same, and the difference terms vanish
      mpd  = sd_m1t_lx + fac_b*sd_d1t_lx
      mmd  = sd_m1t_lx - fac_b*sd_d1t_lx
      nme  = 0
      npe  = 0
    end if
    !
    !  Apply the inverse of the M^-T matrices to the input wavefunctions
    !
    if (lval==0) then
      call m3d_solve(sd_m1tf_lx,psi_l,zm(:,1)) ! ZM = M^-T_LX Psi_l
      call m3d_solve(sd_m1tf_l0,psi_p,zp(:,1)) ! ZP = M^-T_L0 Psi_p
    else
      call m3d_solve(sd_m1tf_lx,psi_l,zm(:,1)) ! ZM = M^-T_LX Psi_l
      call m3d_solve(sd_m1tf_lx,psi_p,zp(:,1)) ! ZP = M^-T_LX Psi_p
    end if
    !
    !  "Mix" the input wavefunction
    !
    xp(:,1) = zm(:,1) + zp(:,1)
    xm(:,1) = zm(:,1) - zp(:,1)
    !
    !  Construct the right-hand side:
    !
    !   zp = (M+D) . xp + (N+E) . xm
    !   zm = (M-D) . xm + (N-E) . xp
    !
    call m3d_multiply(mpd,xp(:,1),zp(:,1))
    call m3d_multiply(mmd,xm(:,1),zm(:,1))
    zp(1:2,1) = zp(1:2,1) + npe*xm(1,1)
    zm(1:2,1) = zm(1:2,1) + nme*xp(1,1)
    !
    !  Now construct the "unperturbed" linear system matrix for Sherman-Morrison.
    !  If you pay attention, it turns out to coinside with the laser-coupling 
    !  operator for L=lval.
    !
    mpd((/1,2/),1) = mpd((/1,2/),1) - npe
    mmd((/1,2/),1) = mmd((/1,2/),1) - nme
    call m3d_decompose(mpd,mpd_f)
    call m3d_decompose(mmd,mmd_f)
    !
    !  Construct the auxiliary right-hand side for the Sherman-Morrison correction
    !  This is vector "u"
    !
    zp(:,2) = 0 ; zp(1:2,2) = nme
    zm(:,2) = 0 ; zm(1:2,2) = npe
    !
    !  Solve the unperturbed linear system and the correction equation
    !
    call m3d_solve(mmd_f,zp,xp)
    call m3d_solve(mpd_f,zm,xm)
    !
    !  Correction factor: (v.y)/(1+(v.z))
    !  Vector v has only two non-zero elements (both unity)
    !  y is the unperturbed solution; z is the correction-vector solution
    !
    smfact = (xp(1,1)+xm(1,1))/(1+xp(1,2)+xm(1,2))
    !
    !  Apply the correction
    !
    xp(:,1) = xp(:,1) - smfact*xp(:,2)
    xm(:,1) = xm(:,1) - smfact*xm(:,2)
    !
    !  Unmix
    !
    zp(:,1) = 0.5_rk*(xp(:,1)+xm(:,1))
    zm(:,1) = 0.5_rk*(xp(:,1)-xm(:,1))
    !
    !  Scale by the appropriate M^T matrix, and we are done
    !
    if (lval==0) then
      call m3d_multiply(sd_m1t_lx,zp(:,1),psi_l)  !  psi_l = M^T_LX zp
      call m3d_multiply(sd_m1t_l0,zm(:,1),psi_p)  !  psi_p = M^T_L0 zm
    else
      call m3d_multiply(sd_m1t_lx,zp(:,1),psi_l)  !  psi_l = M^T_LX zp
      call m3d_multiply(sd_m1t_lx,zm(:,1),psi_p)  !  psi_p = M^T_LX zm
    end if
  end subroutine op_lmix_sherman_morrison_t
  !
  subroutine op_lmix_mixunmix_t(lval,fac_b,psi_l,psi_p, &
                                zp,zm,xp,xm,rhs,tmp,leqp_c,leqpf_c,leqm_c,leqmf_c,leqp_r,leqpf_r,leqm_r,leqmf_r)
    integer(ik), intent(in)    :: lval     ! L value of the first sub-block; second sub-block is L+1
    complex(rk), intent(in)    :: fac_b    ! Amplitude for the coupling operator in the propagator
    complex(rk), intent(inout) :: psi_l(:) ! Radial part of the L block of the wavefunction
    complex(rk), intent(inout) :: psi_p(:) ! Radial part of the L+1 block of the wavefunction
    ! Intel Fortran serializes on memory allocation; as the result, we must pass
    ! local temporaries in op_lmix_mixunmix_t() from upstream to get reasonable
    ! parallel performance. This sucks :O(
    complex(rk), intent(inout) :: zp(:), zm(:)
    complex(rk), intent(inout) :: xp(:), xm(:)
    complex(rk), intent(inout) :: rhs(:), tmp(:)
    complex(rk), intent(inout) :: leqp_c(:,:), leqpf_c(:,:)
    complex(rk), intent(inout) :: leqm_c(:,:), leqmf_c(:,:)
    real(rk), intent(inout)    :: leqp_r(:,:), leqpf_r(:,:)
    real(rk), intent(inout)    :: leqm_r(:,:), leqmf_r(:,:)
    !
  ! complex(rk) :: zp(sd_nradial), zm(sd_nradial)
  ! complex(rk) :: xp(sd_nradial), xm(sd_nradial)
  ! complex(rk) :: rhs(sd_nradial), tmp(sd_nradial)
  ! complex(rk) :: leqp_c(3,sd_nradial), leqpf_c(3,sd_nradial)
  ! complex(rk) :: leqm_c(3,sd_nradial), leqmf_c(3,sd_nradial)
  ! real(rk)    :: leqp_r(3,sd_nradial), leqpf_r(3,sd_nradial)
  ! real(rk)    :: leqm_r(3,sd_nradial), leqmf_r(3,sd_nradial)
    !
    if (lval==0) stop 'propagator_tools%op_lmix_mixunmix_t - lval=0 is not valid'
    !
    zp(:) = psi_l(:) + psi_p(:)
    zm(:) = psi_l(:) - psi_p(:)
    !
    !  Construct and factor linear system matrices
    !
    if (pt_sense_real .and. aimag(fac_b)==0._rk) then
      leqp_r = sd_m1t_lx + real(fac_b,kind=rk)*sd_d1t_lx
      leqm_r = sd_m1t_lx - real(fac_b,kind=rk)*sd_d1t_lx
      call m3d_decompose(leqp_r,leqpf_r)
      call m3d_decompose(leqm_r,leqmf_r)
      !
      !  Solve for the mixed propagator
      !
      call m3d_solve(sd_m1tf_lx,zp,tmp)
      call m3d_multiply(leqp_r,tmp,rhs)
      call m3d_solve(leqmf_r,rhs,tmp)
      call m3d_multiply(sd_m1t_lx,tmp,xp)
      !
      call m3d_solve(sd_m1tf_lx,zm,tmp)
      call m3d_multiply(leqm_r,tmp,rhs)
      call m3d_solve(leqpf_r,rhs,tmp)
      call m3d_multiply(sd_m1t_lx,tmp,xm)
    else
      leqp_c = sd_m1t_lx + fac_b*sd_d1t_lx
      leqm_c = sd_m1t_lx - fac_b*sd_d1t_lx
      call m3d_decompose(leqp_c,leqpf_c)
      call m3d_decompose(leqm_c,leqmf_c)
      !
      !  Solve for the mixed propagator
      !
      call m3d_solve(sd_m1tf_lx,zp,tmp)
      call m3d_multiply(leqp_c,tmp,rhs)
      call m3d_solve(leqmf_c,rhs,tmp)
      call m3d_multiply(sd_m1t_lx,tmp,xp)
      !
      call m3d_solve(sd_m1tf_lx,zm,tmp)
      call m3d_multiply(leqm_c,tmp,rhs)
      call m3d_solve(leqpf_c,rhs,tmp)
      call m3d_multiply(sd_m1t_lx,tmp,xm)
    end if
    !
    !  Unmix, and we are done
    !
    psi_l(:) = 0.5_rk*(xp+xm)
    psi_p(:) = 0.5_rk*(xp-xm)
  end subroutine op_lmix_mixunmix_t
  !
  !
  !  Hla2 = (e^2/(m c^2)) A^2
  !
  !  This is a diagonal operator; it is the same for normal and transpose forms.
  !
  !  The actual operator we apply is:
  !
  !  [1 + I (dt) Hla2]^-1 [1 - I (dt) Hla2]
  !
  !  In principle, since this is a scalar, we could also use exp(-I (dt) Hla2),
  !  but given that all the remaining terms in the propagator are of lower accuracy,
  !  what would be the point?
  !
  subroutine op_la2(a,dt,psi_l)
    real(xk), intent(in)       :: a        ! Vector-potential along the Z axis
    complex(xk), intent(in)    :: dt       ! Time step
    complex(rk), intent(inout) :: psi_l(:) ! Radial part of the L block of the wavefunction
    !
    complex(xk) :: hla2, scl
    !
    hla2  = 0.5_xk * (electron_charge**2 / electron_mass) * a**2
    scl   = (1-(0,1)*dt*hla2)/(1+(0,1)*dt*hla2)
    !scl   = exp(-(0,1)*dt*hla2)
    psi_l = cmplx(scl,kind=rk)*psi_l
  end subroutine op_la2
end module propagator_tools
