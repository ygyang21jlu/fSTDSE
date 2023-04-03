!
!  Solving atomic TDSE in the presence of a strong laser field.
!  The structure of this code is strongly inspired by the paper:
!
!   H.G. Muller, An Efficient Propagation Scheme for the Time-Dependent 
!   Schroedinger Equation in the Velocity Gauge, Laser Physics, 9, 138-148 (1999)
!
!  However, the capabilities of the code, and some of the numerical 
!  details are quite different from HGM's.
!
!  In parallel runs, this code may benefit significantly from large page support,
!  especially on systems with 2K default pages and small TLB (x86_64!). Under
!  Linux, enabling large-page support may requre adjusting the hugepage limits,
!  e.g.:
!
!   sysctl vm.nr_hugepages=128
!   sysctl vm.nr_overcommit_hugepages=4096
!
!  and linking with libhugepage. One would also need to set HUGETLB_MORECORE=yes
!  or HUGETLB_MORECORE=thp (if the kernel supports it) environment variable.
!  Executing:
!
!   hugeedit --data spherical_tdse.x
!
!  may also be a good idea.
!
!  Note that linking with libhugetlbfs prevents a possibility of static linking
!
module spherical_tdse
  use accuracy
  use bicg_tools
  use cap_tools
  use composition_analysis
  use constants
  use math
  use potential_tools
  use propagator_tools
  use rotation_tools
  use spherical_data
  use test_tools
  use timer
  use tridiagonal_tools
  use vectorpotential_tools
  use wavefunction_tools
  implicit none
  private
  public start
  !
  integer, parameter       :: iu_detail             = 29           ! Unit for detailed output; remains open during the entire run
  integer, parameter       :: iu_temp               = 22           ! An arbitrary unit number, which can be used here
  !                                                 
  integer(ik)              :: verbose               = 2_ik         ! How verbose do we need to be?
  integer(ik)              :: omp_num_threads       = 0_ik         ! Non-zero value will cause number of OpenMP threads
                                                                   ! to be set explicitly. It looks like some libraries
                                                                   ! mess with the environment variables, and stop
                                                                   ! OMP_NUM_THREADS from working.
  character(len=clen)      :: comment               = ' '          ! Descriptive string, to be copied to the output
  character(len=clen)      :: initial_wfn           = 'atomic'     ! Choice of the initial wavefunction. One of:
                                                                   ! 'random' = Use random initial guess
                                                                   ! 'unit'   = Use unit vector as the initial guess
                                                                   ! 'atomic' = Use atomic field-free solution (full matrix)
                                                                   !            Must supply initial_wfn_index.
                                                                   ! 'single' = Use atomic filed-free solution (single solution)
                                                                   !            Must supply initial_wfn_index and initial_wfn_energy
                                                                   ! 'read'   = Read initial wfn from a file
  integer(ik)              :: initial_wfn_index(3)  = (/0,0,1/)    ! Initial wavefunction quantum numbers; L,M,I
  complex(rk)              :: initial_wfn_energy    = -1._rk       ! Initial guess for wavefunction energy
  character(len=clen)      :: initial_wfn_file      = ' '          ! File containing initial wavefunction
  character(len=clen)      :: task                  = 'real time'  ! What to do after the initialization. One of:
                                                                   ! 'real time'      = Propagate in real time
                                                                   ! 'imaginary time' = Propagate in imaginary time
  logical                  :: skip_tests            = .false.      ! Skip all non-essential tests, and assume user knows
                                                                   ! where she is going
  real(xk)                 :: dt                    = 0.01_xk      ! Time step, in atomic units of time
  integer(ik)              :: timesteps             = 50000_ik     ! Number of time steps to perform
  character(len=clen)      :: rotation_mode         = 'auto'       ! Choice of the frame rotation implementation. One of:
                                                                   ! 'none'        - Use fixed frame. Vector-potential must be along Z;
                                                                   !                 other vector-potentials will be accepted, but 
                                                                   !                 (silently) produce incorrect results.
                                                                   !                 Choosing 'none' will also force field_unwrap=.false.
                                                                   ! 'brute force' - Explicit finite-angle rotation, using full-matrix
                                                                   !                 implementation through Wigner rotation matrices.
                                                                   !                 This version becomes numerically unstable for large
                                                                   !                 angular momenta, and is intended primarily for 
                                                                   !                 debugging.
                                                                   ! 'sparse'      - Small-angle rotation using sparse matrices. Will
                                                                   !                 implement large-angle rotation by breaking it up
                                                                   !                 in small steps (see rt_max_rot below, and the code)
                                                                   ! 'auto'        - Switches between 'none' and 'sparse' depending 
                                                                   !                 on the vector potential.
  logical                  :: field_unwrap          = .true.       ! Try to unwrap vector-potential in spherical coordinates
  real(rk)                 :: unwrap_threshold      = 1e-8_rk      ! Do not try to manipulate vector potential magnitude smaller 
                                                                   ! than unwrap_threshold * max(vp)
  character(len=clen)      :: field_preview         = 'field.table'! File containing laser field during the simulation.
  integer(ik)              :: output_each           = 100_ik       ! Reduce summary output by a factor output_each
  character(len=clen)      :: detail_output         = 'detail.table'! File containing full-accuracy results from the simulation
  real(rk)                 :: composition_threshold = 1e-10_rk     ! Threshold for reporting the field-free amplitudes
                                                                   ! Set to a negative number to disable composition analysis
  character(len=clen)      :: initial_wf_dump_prefix= ' '          ! Initial wavefunction dump in human-readable from
  character(len=clen)      :: final_wf_dump_prefix  = ' '          ! Final radial wavefunction in human-readable form is dumped in
                                                                   ! files with final_wf_dump_prefix used as a prefix. Empty
                                                                   ! string suppresses the dump
  character(len=clen)      :: vp_table              = 'vp.table'   ! Relevant if vp_shape='table'. The contents of the vpot_table() 
                                                                   ! array, in the Fortran free-form input format, are expected.
  type(sd_wfn)             :: wfn_r                                ! Our wavefunction (right)
  type(sd_wfn)             :: wfn_l                                ! Ditto (left)
  real(xk), allocatable    :: vpot_table(:,:)                      ! Vector-potential parameters table for the entire simulation
                                                                   ! First index: (0) = time, (1) = signed magnitude, (2) = theta, (3) = phi
                                                                   ! Second index: timestep, from -1 to 2*timesteps+1; even numbers are whole
                                                                   ! time steps; odd numbers are half-steps. 
                                                                   ! These data must be handled in at least double precision, or bad things
                                                                   ! will happen in time propagation
  real(xk), allocatable    :: efield_table(:,:)                    ! -(d A/d t) in the laboratory frame, otherwise known as the electric
                                                                   ! field. This quantity is needed for computing dipole acceleration.
                                                                   ! First index: 1=X, 2=Y, 3=Z components of the electric field
                                                                   ! Second index: timestep, from 0 to 2*timesteps; note that this is NOT
                                                                   ! the same as in vpot_table.
                                                                   ! efield_table is derived from the contents of vpot_table in
                                                                   ! fill_efield_table() below.
  !
  !  Simulation parameters; we are collecting variables from many modules.
  !
  namelist /sph_tdse/ &
                      ! Parameters defined locally
                      verbose, comment, &
                      omp_num_threads, &
                      initial_wfn, initial_wfn_index, initial_wfn_energy, initial_wfn_file, &
                      task, dt, timesteps, rotation_mode, &
                      field_unwrap, unwrap_threshold, field_preview, skip_tests, &
                      output_each, detail_output, &
                      composition_threshold, final_wf_dump_prefix, &
                      initial_wf_dump_prefix, &
                      vp_table, &
                      ! Parameters from spherical_data
                      sd_nradial, sd_nspin, sd_lmax, sd_mmin, sd_mmax, &
                      sd_rgrid, sd_rgrid_zeta, sd_rgrid_dr, sd_rgrid_r0, sd_rgrid_scale, sd_rgrid_npad, &
                      sd_rgrid_file, sd_rgrid_grad_delta, &
                      ! Parameters from potential_tools
                      pot_name, pot_param, &
                      ! Parameters from propagator_tools
                      pt_mix_solver, pt_chunk, pt_sense_real, pt_eps_dt, pt_force_par_l, &
                      ! Parameters from vectorpotential_tools
                      vp_scale, vp_scale_x, vp_shape, vp_param, vp_param_x, &
                      ! Parameters from timer
                      timer_disable, &
                      ! Parameters from bicg_tools
                      bicg_failtrace, bicg_epsilon, bicg_maxiter, &
                      ! Parameters from rotation_tools
                      rt_max_rot, rt_blocksize, rt_sense_real, &
                      ! Parameters from wavefunction_tools
                      wt_atomic_cache_prefix, wt_iterative_improvement, &
                      wt_disable_orthogonalization, wt_max_solution_iterations, &
                      ! Parameters from cap_tools
                      cap_name, cap_param
  !
  contains
  !
  subroutine random_initial_wfn(tag,wfn)
    character(len=*), intent(in) :: tag
    type(sd_wfn), intent(inout)  :: wfn
    integer(ik)                  :: lval, mval, sval, ipt
    real(rk)                     :: hr(2)
    !
    call random_seed()
    loop_m: do mval=sd_mmin,sd_mmax
      loop_l: do lval=abs(mval),sd_lmax
        loop_s: do sval=1,sd_nspin
          loop_pt: do ipt=1,sd_nradial
            call random_number(hr)
            wfn%wfn(:,sval,lval,mval) = cmplx(hr(1),hr(2),kind=rk)
          end do loop_pt
        end do loop_s
      end do loop_l
    end do loop_m
    write (out,"('Initial ',a,' wavefunction set to random values')") trim(tag)
  end subroutine random_initial_wfn
  !
  subroutine atomic_initial_wfn
    integer(ik)              :: ipt, lval, mval, ind, nvec, alloc
    complex(rk), allocatable :: evec(:,:,:)   ! Eigenvectors
    complex(rk), allocatable :: eval(:)       ! Eigenvalues
    !
    lval = initial_wfn_index(1)
    mval = initial_wfn_index(2)
    ind  = initial_wfn_index(3)
    nvec = sd_nradial
    if (initial_wfn=='single') then
      ind  = 1
      nvec = 1
    end if
    !
    write (out,"('Choosing atomic solution with L=',i0,' M=',i0,' I=',i0)") lval, mval, ind
    if (lval<0 .or. lval>sd_lmax .or. mval<sd_mmin .or. mval>sd_mmax .or. &
        abs(mval)>lval .or. ind<1 .or. ind>sd_nradial) then
      stop 'spherical_tdse%atomic_initial_wfn - bad initial_wfn_index'
    end if
    !
    allocate (evec(sd_nradial,nvec,2),eval(nvec),stat=alloc)
    if (alloc/=0) stop 'spherical_tdse%atomic_initial_wfn - allocation failed'
    select case (initial_wfn)
      case default
        stop 'spherical_tdse%atomic_initial_wfn - bad initial_wfn'
      case ('atomic')
        call wt_atomic_solutions(verbose,lval,eval,evec)
      case ('single')
        eval(ind) = initial_wfn_energy
        call wt_one_atomic_solution(verbose,lval,eval(ind),evec(:,ind,:))
    end select
    !
    write (out,"('Energy of the atomic solution is ',2(g28.16e3,1x),'Hartree')") eval(ind)
    wfn_l%wfn = 0
    wfn_l%wfn(:,1,lval,mval) = evec(:,ind,1)
    wfn_r%wfn = 0
    wfn_r%wfn(:,1,lval,mval) = evec(:,ind,2)
    !
    if (verbose>=2) then 
      write (out,"(/'Initial radial wavefunction was:'/)")
      write (out,"((1x,a6,5(1x,a24)))") &
             ' I ', '  R,Bohr ', '  Re(psi_l)  ', '  Im(psi_l)  ', '  Re(psi_r)  ', '  Im(psi_r)  ', &
             '---', '---------', '-------------', '-------------', '-------------', '-------------'
      print_radial: do ipt=1,sd_nradial
        write (out,"(1x,i6,5(1x,g24.13e3))") ipt, sd_rtab(ipt), evec(ipt,ind,:)
      end do print_radial
      write (out,"()")
    end if
    !
    deallocate (evec,eval)
  end subroutine atomic_initial_wfn
  !
  subroutine prepare_initial_wavefunction
    real(rk)    :: ram
    complex(rk) :: norm(2)
    integer(ik) :: alloc
    !
    call TimerStart('Initial wavefunction')
    ram = 2*2*rk_bytes()*real(sd_nradial,kind=rk)*sd_nspin*real(sd_lmax,kind=rk)*real(sd_mmax-sd_mmin+1,kind=rk)
    write (out,"('Wavefunction requires ',f12.3,' Mbytes of memory')") ram/(1024._rk**2)
    call flush_wrapper(out)
    allocate (wfn_l%wfn(sd_nradial,sd_nspin,0:sd_lmax,sd_mmin:sd_mmax), &
              wfn_r%wfn(sd_nradial,sd_nspin,0:sd_lmax,sd_mmin:sd_mmax),stat=alloc)
    if (alloc/=0) then
      write (out,"('Allocation of wavefunction array failed, code = ',i0)") alloc
      stop 'spherical_tdse%prepare_initial_wavefunction - out of memory'
    end if
    !
    select case (initial_wfn)
      case default
        write (out,"('Initial wavefunction choice ',a,' is not recognized')") trim(initial_wfn)
        stop 'spherical_tdse%prepare_initial_wavefunction - bad initial_wfn'
      case ('random')
        call random_initial_wfn('left',wfn_l)
        call random_initial_wfn('right',wfn_r)
      case ('unit')
        wfn_l%wfn = 1
        wfn_r%wfn = 1
      case ('atomic','single')
        call atomic_initial_wfn
      case ('read')
        call fetch_wavefunctions(initial_wfn_file)
    end select
    call wt_normalize(wfn_l,wfn_r,norm)
    write (out,"('        Initial wavefunction norm was ',g24.13e3,1x,g24.13e3)") norm(1)
    write (out,"('Right wavefunction Cartesian norm was ',g24.13e3)") real(norm(2),kind=rk)
    call TimerStop('Initial wavefunction')
  end subroutine prepare_initial_wavefunction
  !
  subroutine fill_vpot_table
    integer(ik) :: its, alloc
    real(xk)    :: time, vp, th, ph
    !
    call TimerStart('Prepare VP table')
    allocate (vpot_table(0:3,-1:2*timesteps+1),stat=alloc)
    if (alloc/=0) stop 'spherical_tdse%fill_vpot_table - allocation failed'
    if (vp_shape=='table') then
      open (iu_temp,form='formatted',recl=256,action='read',position='rewind',status='old',file=trim(vp_table))
      read (iu_temp,*) vpot_table
      close (iu_temp)
    else
      vp = vp_apot(0._xk) ! Initialize internal structures of vp_apot
      !
      !  Start by filling the table 
      !
      time_steps: do its=-1,2*timesteps+1
        ! Full step
        time = dt*0.5_xk*its
        vp   = vp_apot(time,theta=th,phi=ph)
        vpot_table(:,its) = (/ time, vp, th, ph /)
      end do time_steps
    end if
    call TimerStop('Prepare VP table')
  end subroutine fill_vpot_table
  !
  subroutine unwrap_vpot_table
    integer(ik) :: its
    real(xk)    :: vp, th, ph, vp_extrap, th_ref, ph_ref
    real(xk)    :: th_v1, ph_v1, th_v2, ph_v2, r2_v1, r2_v2
    real(xk)    :: vp_max
    !
    if (.not.field_unwrap .or. rotation_mode=='none') then
      write (out,"(/'Skipping vector-potential unwrapping'/)")
      return
    end if
    call TimerStart('Unwrap VP table')
    !
    !  Spherical representation is not necessarily continuous in the normal parameter domain
    !  vp [0:inf), th [0:pi], ph [0:2pi)
    !  Our propagator is much happier when all three parameters to be continuos though the 
    !  first order (althogh it should be able to handle rapid change in theta and phi), 
    !  so that we need to unwrap the v.p.
    !
    !  Step 1: make sure V.P. magnitude's derivative is continuous. Points -1 and 0 are
    !          assumed to be "good" already, and will be used to start the process.
    ! 
    vp_max = maxval(abs(vpot_table(1,:)))
    unwrap_magnitude: do its=1,2*timesteps+1
      vp_extrap = 2*vpot_table(1,its-1) - vpot_table(1,its-2)
      vp        =   vpot_table(1,its)
      ! Try to avoid messing with vector-potential when it is nearly zero
      if (abs(vp_extrap+vp)>abs(vp_extrap-vp) .or. (abs(vp)+abs(vp_extrap))<unwrap_threshold*vp_max) cycle unwrap_magnitude
      ! We get smoother vector-potential by flipping the magnitude; do it!
      vpot_table(1,its) = -vpot_table(1,its)
      vpot_table(2,its) = -pi + vpot_table(2,its)
    end do unwrap_magnitude
    !
    !  Step 2: keeping VP magnitude constant, try to choose the smoothest possible version of (theta,phi)
    !
    !   We have two main possibilities:
    !   1.  th =  th0 + 2pi*n ; ph = ph0 + 2pi*k
    !   2   th = -th0 + 2pi*n ; ph = ph0 + 2pi*k + pi
    !
    unwrap_theta_phi: do its=0,2*timesteps+1
      th     = vpot_table(2,its) 
      ph     = vpot_table(3,its) 
      th_ref = vpot_table(2,its-1)
      ph_ref = vpot_table(3,its-1)
      th_v1  = nearest_2pi( th,   th_ref)
      ph_v1  = nearest_2pi( ph,   ph_ref)
      th_v2  = nearest_2pi(-th,   th_ref)
      ph_v2  = nearest_2pi( ph+pi,ph_ref)
      r2_v1  = (th_v1-th_ref)**2 + (ph_v1-ph_ref)**2
      r2_v2  = (th_v2-th_ref)**2 + (ph_v2-ph_ref)**2
      if (r2_v1<=r2_v2) then
        vpot_table(2,its) = th_v1
        vpot_table(3,its) = ph_v1
      else
        vpot_table(2,its) = th_v2
        vpot_table(3,its) = ph_v2
      end if
    end do unwrap_theta_phi
    call TimerStop('Unwrap VP table')
    contains
    function nearest_2pi(x,ref) result(v)
      real(xk), intent(in) :: x    
      real(xk), intent(in) :: ref 
      real(xk)             :: v
      !
      v = x + 2*pi_xk*nint((ref-x)/(2*pi_xk),kind=ik)
    end function nearest_2pi
  end subroutine unwrap_vpot_table
  !
  !  We calculate the electric field using numerical differentiation of the 
  !  vector-potential. Since we require high numerical accuracy, we use
  !  implicit derivative expressions.
  !
  subroutine fill_efield_table
    integer(ik)            :: its, alloc
    real(xk), allocatable  :: vpot_xyz(:,:)              ! Vector-potential in Cartesian laboratory coordinates
    real(xk), allocatable  :: dt(:)                      ! dt(i) is the timestep leading to the time at point (i)
    real(xk), allocatable  :: d1(:,:), m1(:,:), m1f(:,:) ! tri-diagonal matrices defining the implict gradient
    real(xk), allocatable  :: tmp1(:,:), tmp2(:,:)       ! Temporary arrays
    !
    call TimerStart('Prepare efield table')
    !
    !  efield_table will remain until the end of the run. Remaining arrays must be
    !  deallocated before leaving fill_efied_table.
    !
    allocate (efield_table(3,0:2*timesteps), &
              vpot_xyz(0:2*timesteps,3), dt(0:2*timesteps+1), &
              d1(3,0:2*timesteps), m1(3,0:2*timesteps), m1f(3,0:2*timesteps), &
              tmp1(0:2*timesteps,3), tmp2(0:2*timesteps,3), &
              stat=alloc)
    if (alloc/=0) stop 'spherical_tdse%fill_efield_table - allocation failed'
    !
    !  Calculate Cartesian vector-potential using spherical vector-potential in vp_apot()
    !  Also build the table of the timesteps.
    !
    fill_vpot_xyz: do its=0,2*timesteps
      vpot_xyz(its,:) = vpot_sph2xyz(vpot_table(1:3,its))
      dt(its)         = vpot_table(0,its) - vpot_table(0,its-1)
    end do fill_vpot_xyz
    dt(2*timesteps+1) = vpot_table(0,2*timesteps+1)-vpot_table(0,2*timesteps)
    !
    !  We need to prepare the Delta_1 and M_1 matrices for the first derivative.
    !  The code snipped below is lifted from initialize_radial_gradient() in spherical_data.f90
    !  We can't just reuse initialize_radial_gradient, since the boundary conditions are different
    !
    grad_tables: do its=0,2*timesteps
      ! Diagonal
      d1(1,its) = 1/dt(its) - 1/dt(its+1) + (-dt(its) + dt(its+1))/ (dt(its)**2 + dt(its)*dt(its+1) + dt(its+1)**2)
      m1(1,its) = (dt(its) + dt(its+1))**2/ (2*(dt(its)**2 + dt(its)*dt(its+1) + dt(its+1)**2))
      if (its>=2*timesteps) cycle grad_tables
      ! Sub-diagonal
      d1(2,its) = -((dt(its+2)**2*(2*dt(its+1) + dt(its+2)))/ (dt(its+1)*(dt(its+1) + dt(its+2))* &
                         (dt(its+1)**2 + dt(its+1)*dt(its+2) + dt(its+2)**2)))
      m1(2,its) = dt(its+2)**2/(2*(dt(its+1)**2 + dt(its+1)*dt(its+2) + dt(its+2)**2))
      ! Super-diagonal
      d1(3,its) = (dt(its)**2*(dt(its) + 2*dt(its+1)))/(dt(its+1)* (dt(its) + dt(its+1))*(dt(its)**2 + &
                         dt(its)*dt(its+1) + dt(its+1)**2))
      m1(3,its) = dt(its)**2/(2*(dt(its)**2 + dt(its)*dt(its+1) + dt(its+1)**2))
    end do grad_tables
    call m3d_decompose_x(m1,m1f)
    !
    !  We are done with all the matrices; the derivative is now given by:
    !
    !    M1^-1 Delta_1 A
    !
    call m3d_multiply_x(d1,vpot_xyz,tmp1)
    call m3d_solve_x(m1f,tmp1,tmp2)
    efield_table = -transpose(tmp2)
    !
    deallocate (vpot_xyz,dt,d1,m1,m1f,tmp1,tmp2)
    call TimerStop('Prepare efield table')
    contains
    function vpot_sph2xyz(atp)
      real(xk), intent(in) :: atp(3)
      real(xk)             :: vpot_sph2xyz(3)
      !
      vpot_sph2xyz(1) = atp(1)*sin(atp(2))*cos(atp(3))
      vpot_sph2xyz(2) = atp(1)*sin(atp(2))*sin(atp(3))
      vpot_sph2xyz(3) = atp(1)*cos(atp(2))
    end function vpot_sph2xyz
  end subroutine fill_efield_table
  !
  subroutine preview_laser_field
    integer(ik) :: its
    real(xk)    :: tatp(0:3) ! (time, a, theta, phi)
    real(xk)    :: efield(3) ! (Ex,Ey,Ez)
    !
    if (field_preview==' ') return
    call TimerStart('Dump VP table')
    !
    open (iu_temp,form='formatted',action='write',position='rewind',status='replace',file=trim(field_preview))
    !
    write (iu_temp,"(('#',a11,1x,a14,3(1x,a20,4x),1x,3(1x,a20,4x)))") &
           ' step ', ' Time, au[t] ', ' V.P. magnitude ', ' V.P. theta ', ' V.P. phi ', ' Ex ', ' Ey ', ' Ez ', &
           '------', '-------------', '----------------', '------------', '----------', '----', '----', '----'
    time_steps: do its=0,2*timesteps
      tatp   = vpot_table(:,its)
      efield = efield_table(:,its)
      write (iu_temp,"(1x,f11.1,1x,f14.6,3(1x,g24.13e3),1x,3(1x,g24.13e3))") 0.5_rk*its, tatp, efield
    end do time_steps
    !
    close (iu_temp)
    call TimerStop('Dump VP table')
  end subroutine preview_laser_field
  !
  subroutine dump_wavefunctions(prefix)
    character(len=*), intent(in) :: prefix
    integer(ik)                  :: lval, mval, sval, ir
    character(len=clen)          :: filename
    !
    call TimerStart('Dump wavefunction')
    dump_l_channels: do lval=0,sd_lmax
      dump_m_channels: do mval=max(sd_mmin,-lval),min(lval,sd_mmax)
        if (mval>=0) then
          write (filename,"(a,'-L',i0.3,'-M+',i0.4)") trim(prefix), lval, mval
        else
          write (filename,"(a,'-L',i0.3,'-M-',i0.4)") trim(prefix), lval, -mval
        end if
        open(iu_temp,form='formatted',recl=256,action='write',position='rewind',status='replace',file=trim(filename))
        write (iu_temp,"('#',a2,1x,a24,2(2x,a34,1x,a34))") &
              'S', ' R, Bohr ', ' Re[left wfn] ', ' Im[left wfn] ', ' Re[right wfn] ', ' Im[right wfn] '
        dump_s_channels: do sval=1,sd_nspin
          dump_radial: do ir=1,sd_nradial
            write (iu_temp,"(1x,i2,1x,g24.13,2(2x,g34.21e3,1x,g34.21e3))") &
                   sval, sd_rtab(ir), wfn_l%wfn(ir,sval,lval,mval), wfn_r%wfn(ir,sval,lval,mval)
          end do dump_radial
        end do dump_s_channels
        close(iu_temp)
      end do dump_m_channels
    end do dump_l_channels
    call TimerStop('Dump wavefunction')
  end subroutine dump_wavefunctions
  !
  subroutine fetch_wavefunctions(prefix)
    character(len=*), intent(in) :: prefix
    integer(ik)                  :: lval, mval, sval, ir
    character(len=clen)          :: filename
    character(len=1)             :: buf
    integer(ik)                  :: line, ios
    integer(ik)                  :: stmp
    real(rk)                     :: rtmp(5)
    !
    call TimerStart('Fetch wavefunction')
    fetch_l_channels: do lval=0,sd_lmax
      fetch_m_channels: do mval=max(sd_mmin,-lval),min(lval,sd_mmax)
        if (mval>=0) then
          write (filename,"(a,'-L',i0.3,'-M+',i0.4)") trim(prefix), lval, mval
        else
          write (filename,"(a,'-L',i0.3,'-M-',i0.4)") trim(prefix), lval, -mval
        end if
        open(iu_temp,form='formatted',recl=256,action='read',position='rewind',status='old',file=trim(filename))
        line = 1
        read (iu_temp,"(a1)",iostat=ios) buf
        if (ios/=0) then
          write (out,"('Error ',i0,' skipping line ',i0,' of file ',a)") ios, line, trim(filename)
          stop 'spherical_tdse%fetch_wavefunctions - skip error'
        end if
        fetch_s_channels: do sval=1,sd_nspin
          fetch_radial: do ir=1,sd_nradial
            line = line + 1
            read (iu_temp,*,iostat=ios) stmp, rtmp
            if (ios/=0) then
              write (out,"('Error ',i0,' reading line ',i0,' of file ',a)") ios, line, trim(filename)
              stop 'spherical_tdse%fetch_wavefunctions - read error'
            end if
            if (stmp/=sval .or. abs(sd_rtab(ir)-rtmp(1))>1e-10_rk) then
              write (out,"('Grid mismatch reading line ',i0,' of file ',a)") line, trim(filename)
              stop 'spherical_tdse%fetch_wavefunctions - data error'
            end if
            wfn_l%wfn(ir,sval,lval,mval) = cmplx(rtmp(2),rtmp(3),kind=rk)
            wfn_r%wfn(ir,sval,lval,mval) = cmplx(rtmp(4),rtmp(5),kind=rk)
          end do fetch_radial
        end do fetch_s_channels
        close(iu_temp)
      end do fetch_m_channels
    end do fetch_l_channels
    call TimerStop('Fetch wavefunction')
  end subroutine fetch_wavefunctions
  !
  subroutine lab_dipole(wfn_l,wfn_r,th,ph,efield,norm,dipole,acceleration)
    type(sd_wfn), intent(in) :: wfn_l           ! Left wavefunction 
    type(sd_wfn), intent(in) :: wfn_r           ! Right wavefunction 
    real(rk), intent(in)     :: th              ! Theta angle
    real(rk), intent(in)     :: ph              ! Phi angle
    real(rk), intent(in)     :: efield(3)       ! Electric field, needed for the dipole acceleration
    complex(rk), intent(in)  :: norm            ! Wavefunction norm, ditto
    complex(rk), intent(out) :: dipole(3)       ! <L|q r|R> expectation value in the lab frame
    complex(rk), intent(out) :: acceleration(3) ! (d^2/d t^2) <L|q r|R> in the lab frame
    !
    complex(rk) :: loc_dipole(3)       ! Dipole in the local coordinate system
    complex(rk) :: loc_acceleration(3) ! Dipole acceleration in the local coordinate system
    real(rk)    :: rm(3,3)             ! Transformation matrix
    !
    call wt_dipole(wfn_l,wfn_r,loc_dipole,loc_acceleration)
    call MathRotationMatrix((/ph,th,0._rk/),rm)
    dipole       = matmul(transpose(rm),loc_dipole)
    acceleration = matmul(transpose(rm),loc_acceleration)
    !
    !  Add correction term to the dipole acceleration: (e^2/m) E <L|R>
    !
    acceleration = acceleration + (electron_charge**2/electron_mass) * efield * norm
  end subroutine lab_dipole
  !
  subroutine write_detail_header
    if (detail_output==' ') return
    !
    open(iu_detail,file=trim(detail_output),form='formatted',status='replace',position='rewind',recl=800,pad='no')
    write (iu_detail,"('# Field Columns Data')")
    write (iu_detail,"('#  1    2,13    Timestep')")
    write (iu_detail,"('#  2   15,46    Time, au[t]')")
    write (iu_detail,"('#  3   48,79    Vector-potential magnitude')")
    write (iu_detail,"('#  4   81,112   Vector-potential, lab theta')")
    write (iu_detail,"('#  5  114,145   Vector-potential, lab phi')")
    write (iu_detail,"('#  6  147,178   Re[<Psi_L|Psi_R>]')")
    write (iu_detail,"('#  7  180,211   Im[<Psi_L|Psi_R>]')")
    write (iu_detail,"('#  8  213,244   Re[<Psi_L|H_at+H_L+V_cap|Psi_R>], Hartree')")
    write (iu_detail,"('#  9  246,277   Im[<Psi_L|H_at+H_L+V_cap|Psi_R>], Hartree')")
    write (iu_detail,"('# 10  279,310   Re[<Psi_L|H_at+H_L|Psi_R>], Hartree')")
    write (iu_detail,"('# 11  312,343   Im[<Psi_L|H_at+H_L|Psi_R>], Hartree')")
    write (iu_detail,"('# 12  345,376   Re[<Psi_L|e . x|Psi_R>], e-Bohr')")
    write (iu_detail,"('# 13  378,409   Im[<Psi_L|e . x|Psi_R>], e-Bohr')")
    write (iu_detail,"('# 14  411,442   Re[<Psi_L|e . y|Psi_R>], e-Bohr')")
    write (iu_detail,"('# 15  444,475   Im[<Psi_L|e . y|Psi_R>], e-Bohr')")
    write (iu_detail,"('# 16  477,508   Re[<Psi_L|e . z|Psi_R>], e-Bohr')")
    write (iu_detail,"('# 17  510,541   Im[<Psi_L|e . z|Psi_R>], e-Bohr')")
    write (iu_detail,"('# 18  543,574   Re[(d^2/d t^2) <Psi_L|e . x|Psi_R>], e-Bohr/jiffy^2')")
    write (iu_detail,"('# 19  576,607   Im[(d^2/d t^2) <Psi_L|e . x|Psi_R>], e-Bohr/jiffy^2')")
    write (iu_detail,"('# 20  609,640   Re[(d^2/d t^2) <Psi_L|e . y|Psi_R>], e-Bohr/jiffy^2')")
    write (iu_detail,"('# 21  642,673   Im[(d^2/d t^2) <Psi_L|e . y|Psi_R>], e-Bohr/jiffy^2')")
    write (iu_detail,"('# 22  675,706   Re[(d^2/d t^2) <Psi_L|e . z|Psi_R>], e-Bohr/jiffy^2')")
    write (iu_detail,"('# 23  708,739   Im[(d^2/d t^2) <Psi_L|e . z|Psi_R>], e-Bohr/jiffy^2')")
    write (iu_detail,"('# WARNING: Expectation value of the Hamiltonian does not include the non-adiabatic term')")
    if (sd_pot_nonlocal) then
      write (iu_detail,"('# WARNING: Non-local potential detected. Dipole acceleration results are incorrect')")
      write (out,      "(/ 'WARNING: Non-local potential detected. Dipole acceleration results are incorrect'/)")
    end if
    write (iu_detail,"('#',1x,a12,22(1x,a32))") &
           ' i ', ' time ', ' vp ', ' theta ', ' phi ', ' re(norm) ', ' im(norm) ', &
           ' re(energy) ', ' im(energy) ', ' re(en-nocap) ', ' im(en-nocap) ', &
           ' re(dip_x) ', ' im(dip_x) ', ' re(dip_y) ', ' im(dip_y) ', ' re(dip_z) ', ' im(dip_z) ', &
           ' re(acc_x) ', ' im(acc_x) ', ' re(acc_y) ', ' im(acc_y) ', ' re(acc_z) ', ' im(acc_z) '
    write (iu_detail,"('#',1x,a12,22(1x,a32))") &
           ' 1 ', ' 2 ', ' 3 ', ' 4 ', ' 5 ', ' 6 ', ' 7 ', ' 8 ', ' 9 ', ' 10 ', ' 11 ', &
           ' 12 ', ' 13 ', ' 14 ', ' 15 ', ' 16 ', ' 17 ', &
           ' 18 ', ' 19 ', ' 20 ', ' 21 ', ' 22 ', ' 23 '
  end subroutine write_detail_header
  !
  subroutine propagation
    integer(ik) :: its
    real(xk)    :: time, vp, th, ph  ! Handle time and vector-potential in at least double precision
    real(rk)    :: efield(3)         ! The electric field
    complex(xk) :: cdt
    !
    !  Propagation routine can handle varying time steps, at least in principle
    !
    call write_detail_header
    !
    time_steps: do its=0,timesteps-1
      time   = vpot_table(0,2*its)
      vp     = vpot_table(1,2*its)
      th     = vpot_table(2,2*its)
      ph     = vpot_table(3,2*its)
      efield = real(efield_table(:,2*its),kind=kind(efield))
      !
      call write_output(force=.false.)
      !
      !  First half-step; vector-potential is at the current time step
      !  The time step is a half-step.
      !
      cdt = 0.5_rk*(vpot_table(0,2*its+1) - vpot_table(0,2*its))
      call pt_fwd_laser_n(wfn_r,vp,cdt)
      call pt_fwd_laser_t(wfn_l,vp,cdt)
      !
      !  We are at the mid-step; apply atomic propagator for the first half-step
      !
      call pt_fwd_atomic_n(wfn_r,cdt)
      call pt_rev_atomic_n(wfn_r,cdt)
      call pt_fwd_atomic_t(wfn_l,cdt)
      call pt_rev_atomic_t(wfn_l,cdt)
      !
      !  Rotation; nominally a full step. It will be broken into smaller sub-steps
      !  if this turns out to be necessary.
      !
      call rotate(from=vpot_table(2:3,2*its),to=vpot_table(2:3,2*its+2))
      !
      !  Atomic propagator for the second half-step
      !
      cdt = 0.5_rk*(vpot_table(0,2*its+2) - vpot_table(0,2*its+1))
      call pt_fwd_atomic_n(wfn_r,cdt)
      call pt_rev_atomic_n(wfn_r,cdt)
      call pt_fwd_atomic_t(wfn_l,cdt)
      call pt_rev_atomic_t(wfn_l,cdt)
      !
      !  Second half of the time step; vector potential is at the <i>next</i> time-step
      !
      vp = vpot_table(1,2*its+2)
      call pt_rev_laser_n(wfn_r,vp,cdt)
      call pt_rev_laser_t(wfn_l,vp,cdt)
    end do time_steps
    !
    time = vpot_table(0,2*its)
    vp   = vpot_table(1,2*its)
    th   = vpot_table(2,2*its)
    ph   = vpot_table(3,2*its)
    !
    call write_output(force=.true.)
    !
    if (detail_output/=' ') then
      close(iu_detail)
    end if
    !
    !  Rotate wavefunction back into the lab orientation for analysis
    !
    call rotate(from=vpot_table(2:3,2*its),to=(/0._xk,0._xk/))
    if (rotation_mode/='none') then
      write (out,"(/'Wavefunction restored back to laboratory frame'/)")
    end if
    !
    contains 
    !
    subroutine write_output(force)
      logical, intent(in) :: force
      !
      real(rk)            :: abs_dipole
      complex(rk)         :: energy(2), norm
      complex(rk)         :: dipole(3)
      complex(rk)         :: acceleration(3)
      !
      call wt_energy(wfn_l,wfn_r,vp,energy,norm)
      call lab_dipole(wfn_l,wfn_r,real(th,kind=rk),real(ph,kind=rk),efield,norm,dipole,acceleration)
      !
      if (detail_output/=' ') then
        write (iu_detail,"(1x,i12,22(1x,g32.22e3))") its, time, vp, th, ph, norm, energy, dipole, acceleration
      end if
      !
      abs_dipole = sqrt(sum(abs(dipole)**2))
      if (force .or. mod(its,output_each)==0) then
        write (out,"('@ ',i9,' t= ',f12.4,' a= ',f14.6,' th= ',f10.3,' ph= ',f10.3,' <l|r>= ',g23.12,1x,g17.6," // &
                   "' <l|h|r>= ',g23.12,1x,g17.6,' |d|= ',g20.9)") &
               its, time, vp, th, ph, norm, energy(1), abs_dipole
        call flush_wrapper(out)
      end if
    end subroutine write_output
    !
    subroutine rotate(from,to)
      real(xk), intent(in) :: from(:), to(:) ! Initial and final field orientation
      !
      select case (rotation_mode)
        case default
          write (out,"('spherical_tdse%propagation: Rotation mode ',a,' is not recognized')") trim(rotation_mode)
          stop 'spherical_tdse%propagation - bad rotation_mode'
        case ('none')
        case ('sparse')
          call rt_rotate(wfn_r,from=from,to=to,left=.false.)
          call rt_rotate(wfn_l,from=from,to=to,left=.true.)
        case ('brute force')
          call rt_rotate_bruteforce(wfn_r,from=from,to=to,left=.false.)
          call rt_rotate_bruteforce(wfn_l,from=from,to=to,left=.true.)
      end select
    end subroutine rotate
  end subroutine propagation
  !
  subroutine choose_rotation_code
    select case (rotation_mode)
      case default
        write (out,"('spherical_tdse%choose_rotation_mode: rotation mode ',a,' is not recognized')") trim(rotation_mode)
        stop 'spherical_tdse%choose_rotation_mode - bad rotation_mode'
      case ('none','brute force','sparse')
      case ('auto')
        ! If we use a range of M values, presumably we want rotation
        if (sd_mmin/=sd_mmax) rotation_mode = 'sparse'
        ! If we know that VP is along Z, we do not need rotation
        ! This includes 'z Gaussian', 'z Sin2', etc
        if (vp_shape(1:2)=='z ') rotation_mode = 'none'
        ! If we still do not know, we choose the general case (rotation)
        if (rotation_mode=='auto') rotation_mode = 'sparse'
        write (out,"(/'Chosen rotation_mode = ''',a,''''/)") trim(rotation_mode)
    end select 
    !
    !  Issue warnings based on the final choice of rotation_mode
    !
    select case (rotation_mode)
      case default
        stop 'spherical_tdse%choose_rotation_mode - logic error'
      case ('none')
        if (sd_mmin/=sd_mmax) then
          write (out,"(/'WARNING: Multiple, uncoupled M values are propagated.'/)") 
        end if
        if (vp_shape(1:2)/='z ') then
          write (out,"(/'WARNING: Vector-potential may contain components along X/Y, but the propagator')")
          write (out,"( 'WARNING: assumes it is along Z. The results are likely incorrect.'/)")
        end if
      case ('brute force','sparse')
        if (sd_mmin>-sd_lmax .or. sd_mmax<sd_lmax) then
          write (out,"(/'WARNING: Not all possible projections of angular momentum are included in propagation.')")
          write (out,"( 'WARNING: The results are likely incorrect.'/)")
        end if
    end select 
  end subroutine choose_rotation_code
  !
  subroutine math_init
    real(rk) :: dummy
    !
    if (log10(huge(1._rk))>=4930._rk) then
      ! Quad precision; don't forget "factorial_slack"
      dummy = MathFactorial(1750_ik-5_ik)
      dummy = MathLogFactorial(2000_ik)
    else if (log10(huge(1._rk))>=308._rk) then
      ! Double precision
      dummy = MathFactorial(170_ik-5_ik)
      dummy = MathLogFactorial(1000_ik)
    else
      ! Single precision
      dummy = MathFactorial(34_ik-5_ik)
      dummy = MathLogFactorial(500_ik)
    end if
  end subroutine math_init
  !
  subroutine start
    !$ use OMP_LIB
    logical :: have_openmp
    !
    write (out,"('Version: ',a/)") __BUILD_ID__
    write (out,"('  Integer kind = ',i0,' (',i0,' decimals)')") kind(1_ik), int(log10(huge(1_ik)+1._rk))
    write (out,"('     Real kind = ',i0,' (',i0,' decimals)')") kind(1._rk), precision(1._rk)
    write (out,"('Aux. real kind = ',i0,' (',i0,' decimals)')") kind(1._xk), precision(1._xk)
    write (out,"()")
    !
    call TimerStart('start')
    !
    call TimerStart('Initialization')
    !
    read (input,nml=sph_tdse)
    write (out,"(' ===== begin simulation parameters ===== ')")
    write (out,nml=sph_tdse)
    write (out,"(' ====== end simulation parameters ====== ')")
    !
    write (out,"(/a/)") trim(comment)
    !
    have_openmp = .false.
    !$ have_openmp = .true.
    !$ if (omp_num_threads/=0) then
    !$   write (out,"('Forcing number of OpenMP threads to ',i0)") omp_num_threads
    !$   call omp_set_num_threads(omp_num_threads)
    !$ end if
    !$ write (out,"('Maximum number of OpenMP threads for this run is ',i0)") omp_get_max_threads()
    if (omp_num_threads>0 .and. .not.have_openmp) then
      write (out,"(/'WARNING: omp_num_threads input parameter has been specified, but the code was built')")
      write (out,"( 'WARNING: without OpenMP support. Execution will continue on a single CPU core.'/)")
    end if
    !
    call math_init
    !
    !  Set up potential can CAP evaluation
    !
    call pt_initialize
    call cap_initialize
    !
    !  Set up grids and operators
    !
    call sd_initialize
    !
    !  Check accuracy of numetical derivatives
    !
    if (.not.skip_tests) call derivatives_test(verbose)
    !
    !  Construct field-free states using direct diagonalization
    !
    if (.not.skip_tests) call fieldfree_test(verbose)
    !
    !  We need a wavefunction to start from
    !
    call prepare_initial_wavefunction
    !
    if (initial_wf_dump_prefix/=' ') then
      write (out,"(/'Dumping initial wavefunction to disk, prefix = ',a/)") trim(initial_wf_dump_prefix)
      call flush_wrapper(out)
      call dump_wavefunctions(initial_wf_dump_prefix)
    end if
    !
    call TimerStop('Initialization')
    write (out,"(/'Done with initialization'/)")
    call TimerReport
    !
    select case (task)
      case default
        write (out,"('Task ',a,' is not recognized')") trim(task)
        stop 'spherical_tdse - bad task'
      case ('real time')
        write (out,"(/'Real-time propagation'/)")
        call choose_rotation_code
        call fill_vpot_table
        call unwrap_vpot_table
        call fill_efield_table
        call preview_laser_field
        call propagation
      case ('imaginary time')
        write (out,"(/'Imaginary-time propagation, A = ',g24.13/)") vp_scale
        call imaginary_propagation_test(verbose,wfn_l,wfn_r,apot=vp_scale,dt=dt,nstep=timesteps)
    end select
    !
    if (pt_mix_solver=='bi-CG') then
      write (out,"(/'Total number of failures in bi-CG solver = ',i0/)") bicg_failure_count
    end if
    call flush_wrapper(out)
    !
    write (out,"(/'Analyzing wavefunction composition'/)")
    call flush_wrapper(out)
    call ca_analyze(verbose,composition_threshold,wfn_l,wfn_r)
    !
    if (final_wf_dump_prefix/=' ') then
      write (out,"(/'Dumping final wavefunction to disk'/)")
      call flush_wrapper(out)
      call dump_wavefunctions(final_wf_dump_prefix)
    end if
    !
    call TimerStop('start')
    call TimerReport
  end subroutine start
end module spherical_tdse
!
program main
  use spherical_tdse
  call start
end program main
