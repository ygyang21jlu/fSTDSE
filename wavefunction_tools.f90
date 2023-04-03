!
!  Sundry tools for manipulating wavefunction
!
module wavefunction_tools
  use accuracy
  use constants
  use timer
  use spherical_data
  use tridiagonal_tools
  use sort_tools
  use lapack
  use math
  implicit none
  private
  public wt_atomic_cache_prefix, wt_iterative_improvement, wt_disable_orthogonalization
  public wt_max_solution_iterations
  public wt_atomic_solutions, wt_one_atomic_solution
  public wt_normalize, wt_energy, wt_dipole
  !
  integer, parameter  :: iu_temp                = 24  ! An arbitrary unit number, which can be used here
  character(len=clen) :: wt_atomic_cache_prefix = ' ' ! Atomic solution for each L will be cached
                                                      ! here; these quantities are expensive to
                                                      ! recompute, especially if we run a 2D simulation
                                                      ! Blank means no caching
  logical             :: wt_iterative_improvement = .true.
                                                      ! Set to .true. if the atomic eigenstates are to be 
                                                      ! iteratively refined
  logical             :: wt_disable_orthogonalization = .false.
                                                      ! Set to .true. to skip explicit Gram-Schmidt orthogonalization
                                                      ! for the atomic solutions
  integer(ik)         :: wt_max_solution_iterations = 20_ik
                                                      ! Max number of iterations in trying to find the eigenvectors
                                                      ! in wt_one_atomic_solution()
  !
  contains
  !
  !  Calculate a block of atomic solutions. We return the non-zero part of the
  !  eigenvectors only, instead of the much larger total wavefunction.
  !
  subroutine wt_atomic_solutions(verbose,lval,eval,evec)
    integer(ik), intent(in)  :: verbose     ! Reporting level
    integer(ik), intent(in)  :: lval        ! Desired angular momentum
    complex(rk), intent(out) :: eval(:)     ! Eigenvalues of the Hamiltonian matrix block for L=lval
    complex(rk), intent(out) :: evec(:,:,:) ! Left (:,:,1) and right (:,:,2) eigenvectors  
    !
    character(len=clen) :: cache
    integer(ik)         :: ios
    !
    call TimerStart('Atomic solutions')
    if (size(eval)/=sd_nradial .or. any(ubound(evec)/=(/sd_nradial,sd_nradial,2/)) .or. lval<0 .or. lval>sd_lmax) then
      stop 'wavefunction_tools%wt_atomic_solutions - bad arguments'
    end if
    !
    !  Try to load from cache first
    !
    if (wt_atomic_cache_prefix/=' ') then
      !
      !  This function may be called from a parallel region; make sure it
      !  does not try to access the same file from multiple threads!
      !
      !$omp critical
      write (cache,"(a,'-L=',i3.3)") trim(wt_atomic_cache_prefix), lval
      open (iu_temp,action='read',status='old',form='unformatted',iostat=ios,file=trim(cache))
      if (ios==0) then
        read (iu_temp) eval, evec
        close (iu_temp)
      end if
      !$omp end critical
      if (ios==0) then
        call TimerStop('Atomic solutions')
        return
      end if
      !
      ! Open was unsuccessful; proceed through calculation, and try to save the 
      ! results once we have them
      !
    end if
    !
    !  No cached result available; compute it from scratch. 
    !
    !  Begin by computing the guess eigendecomposition using LAPACK.
    !  Unfortunately, LAPACK sometimes produces very inaccurate results for general
    !  complex matrices, so that these eigenvalues/eigenvectors cannot be used as is
    !
    call wt_atomic_solution_guess(verbose,lval,eval,evec)
    !
    call wt_atomic_solution_biorthogonalize(evec)
    !
    if (wt_iterative_improvement) then
      call wt_atomic_solution_improve(verbose,lval,eval,evec)
    end if
    !
    call wt_atomic_solution_verify(lval,evec)
    !
    if (wt_atomic_cache_prefix/=' ') then
      !$omp critical
      open (iu_temp,action='write',status='new',form='unformatted',iostat=ios,file=trim(cache))
      if (ios/=0) then
        write (out,"('Error ',i0,' creating new cache file ',a)") ios, trim(cache)
        write (out,"('Skipping wavefunction save for L=',i0,' and continuing')") lval
      else
        write (iu_temp) eval, evec
        close (iu_temp)
      end if
      !$omp end critical
    end if
    call TimerStop('Atomic solutions')
  end subroutine wt_atomic_solutions
  !
  !  Compute a single atomic solution using inverse iterations.
  !  Requires approximate energy of the solution as the starting guess.
  !
  subroutine wt_one_atomic_solution(verbose,lval,eval,evec)
    integer(ik), intent(in)    :: verbose   ! Reporting level
    integer(ik), intent(in)    :: lval      ! Desired angular momentum
    complex(rk), intent(inout) :: eval      ! In: guess eigenvalue
                                            ! Out: improved eigenvalue
    complex(rk), intent(out)   :: evec(:,:) ! Left (:,1) and right (:,2) eigenvectors  
    !
    real(rk)    :: guess_buf(sd_nradial)
    complex(rk) :: eval_new, delta
    integer(ik) :: iter
    logical     :: converged
    !
    call TimerStart('Atomic solution (single)')
    if (any(ubound(evec)/=(/sd_nradial,2/)) .or. lval<0 .or. lval>sd_lmax) then
      stop 'wavefunction_tools%wt_one_atomic_solution - bad arguments'
    end if
    !
    !  Start with a random real guess. Don't forget to normalize!
    !
    call random_number(guess_buf)
    evec(:,1) = guess_buf
    call random_number(guess_buf)
    evec(:,2) = guess_buf
    call wt_normalize_one_atomic_solution(evec)
    !
    write (out,"('Solving for atomic eigenvector L=',i0,' E(guess)=',2(g26.12,1x),'using inverse iterations.')") lval, eval
    converged = .false.
    inverse_iteration_passes: do iter=1,wt_max_solution_iterations
      eval_new = eval
      call wt_improve_one_atomic_eigenstate(verbose,lval,eval_new,evec)
      delta = eval_new - eval
      eval  = eval_new
      if (verbose>=0) then
        write (out,"(1x,'Iteration ',i0,' E= ',2(g35.24e3,1x),' change = ',2(g15.6e3,1x))") iter, eval, delta
      end if
      if (abs(delta)<=100._rk*spacing(abs(eval_new))) then
        !
        !  We are converged, but do couple more iterations for a good measure
        !
        eval_new = eval
        call wt_improve_one_atomic_eigenstate(verbose,lval,eval_new,evec)
        call wt_improve_one_atomic_eigenstate(verbose,lval,eval_new,evec)
        delta = eval_new - eval
        eval  = eval_new
        if (verbose>=0) then
          write (out,"(1x,'Final iteration E= ',2(g35.24e3,1x),' change = ',2(g15.6e3,1x))") eval, delta
        end if
        converged = .true.
        exit inverse_iteration_passes
      end if
    end do inverse_iteration_passes
    !
    write (out,"('Final energy ',2(g35.24e3,1x))") eval
    if (.not.converged) then
      write (out,"('Inverse iterations failed to converge after ',i0,' passes.')") wt_max_solution_iterations
      write (out,"('Continuing with possibly unconverged eigenvectors.')")
    end if
    !
    call TimerStop('Atomic solution (single)')
  end subroutine wt_one_atomic_solution
  !
  subroutine wt_atomic_solution_guess(verbose,lval,eval,evec)
    integer(ik), intent(in)  :: verbose     ! Reporting level
    integer(ik), intent(in)  :: lval        ! Desired angular momentum
    complex(rk), intent(out) :: eval(:)     ! Eigenvalues of the Hamiltonian matrix block for L=lval
    complex(rk), intent(out) :: evec(:,:,:) ! Left (:,:,1) and right (:,:,2) eigenvectors  
    !
    real(rk), allocatable :: tmat(:,:)
    integer(ik)           :: order(sd_nradial)
    integer(ik)           :: alloc, ipt
    !
    call TimerStart('Atomic solutions: Guess')
    allocate (tmat(sd_nradial,sd_nradial),stat=alloc)
    if (alloc/=0) stop 'wavefunction_tools%wt_atomic_solution_guess - no memory for tmat'
    !
    if (lval==0) then
      call sd_expand_implicit_operator(sd_d2n_l0,sd_m2nf_l0,tmat)
    else
      call sd_expand_implicit_operator(sd_d2n_lx,sd_m2nf_lx,tmat)
    end if
    evec(:,:,1) = (-1._rk/(2._rk*electron_mass)) * tmat
    !
    deallocate (tmat)
    !
    add_potential: do ipt=1,sd_nradial
      evec(ipt,ipt,1) = evec(ipt,ipt,1) + sd_pottab(ipt,1,lval)
    end do add_potential
    if (sd_capped) then
      add_cap: do ipt=sd_cap_start,sd_nradial
        evec(ipt,ipt,1) = evec(ipt,ipt,1) + sd_captab(ipt)
      end do add_cap
    end if
    !
    if (verbose>=2) then
      write (out,"('       Largest matrix element of H for L=',i0,' is ',g24.13)") &
             lval, maxval(abs(evec(:,:,1)))
      write (out,"('Largest deviation from hermiticity for L=',i0,' is ',g24.13)") &
             lval, maxval(abs(evec(:,:,1)-conjg(transpose(evec(:,:,1)))))
    end if
    !
    call lapack_geev(evec,eval)
    call order_keys(real(eval,kind=rk),order)
    eval = eval(order)
    !
    !  LAPACK _GEEV conjugates the left eigenvectors. Let's begin by undoing the damage.
    !
    evec(:,:,1) = conjg(evec(:,order,1))
    evec(:,:,2) =       evec(:,order,2)
    call TimerStop('Atomic solutions: Guess')
  end subroutine wt_atomic_solution_guess
  !
  subroutine wt_atomic_solution_biorthogonalize(evec)
    complex(rk), intent(inout) :: evec(:,:,:) ! Left (:,:,1) and right (:,:,2) eigenvectors  
    !
    integer(ik)  :: imo, jmo
    complex(rk)  :: scl
    ! complex(rk)  :: fact_l(sd_nradial), fact_r(sd_nradial)
    !
    if (wt_disable_orthogonalization) return
    call TimerStart('Atomic solutions: Orthogonalize')
    !
    !  LAPACK _GEEV routines do not give biorthogonal eigenvectors
    !  in the present of multiple roots. Although this should not
    !  occur in our case, it is better to be safe than sorry. Since
    !  we also do not like LAPACK's normalization convention for the
    !  eigenvectors (Cartesian norm = 1 separately for the left and
    !  right vectors), we have to do a bit of work here
    !
    normalize_eigenvectors: do imo=1,sd_nradial
      !
      !  Enforce bi-orthogonality to all lower eigenvectors,
      !  using a variation of the Gram-Schmidt process.
      !  We assume that all lower eigenvectors are already normalized to
      !  our conventions (L.R = 1)
      !
      orthogonalize_lower: do jmo=1,imo-1
        ! Right eigenvector
        scl = sum(evec(:,jmo,1)*evec(:,imo,2))
        evec(:,imo,2) = evec(:,imo,2) - scl*evec(:,jmo,2)
        ! Left eigenvector
        scl = sum(evec(:,jmo,2)*evec(:,imo,1))
        evec(:,imo,1) = evec(:,imo,1) - scl*evec(:,jmo,1)
      end do orthogonalize_lower
      ! The matrix-vector code below should be faster; in fact, it is
      ! actually (much) slower. Oops.
      ! fact_l(:imo-1) = matmul(evec(:,imo,2),evec(:,:imo-1,1))
      ! fact_r(:imo-1) = matmul(evec(:,imo,1),evec(:,:imo-1,2))
      ! evec(:,imo,1) = evec(:,imo,1) - matmul(evec(:,:imo-1,1),fact_r(:imo-1))
      ! evec(:,imo,2) = evec(:,imo,2) - matmul(evec(:,:imo-1,2),fact_l(:imo-1))
      !
      call wt_normalize_one_atomic_solution(evec(:,imo,:))
    end do normalize_eigenvectors
    call TimerStop('Atomic solutions: Orthogonalize')
  end subroutine wt_atomic_solution_biorthogonalize
  !
  !  Normalize right eigenvector to Euclidian norm of 1,
  !  simultaneously making the largest coefficient positive real.
  !
  subroutine wt_normalize_one_atomic_solution(vecx)
    complex(rk), intent(inout) :: vecx(:,:)
    !
    real(rk)    :: nrm
    complex(rk) :: scl
    integer(ik) :: imax
    !
    if (any(ubound(vecx)/=(/sd_nradial,2/))) stop 'wavefunction_tools%wt_normalize_one_atomic_solution - bad arguments'
    !
    ! Normalize right eigenvector
    !
    nrm       = sum(abs(vecx(:,2))**2)
    imax      = maxloc(abs(vecx(:,2)),dim=1)
    scl       = abs(vecx(imax,2)) / (vecx(imax,2) * sqrt(nrm))
    vecx(:,2) = scl * vecx(:,2)
    !
    ! Normalize left eigenvector
    !
    scl       = sum(vecx(:,1)*vecx(:,2))
    scl       = 1._rk/scl
    vecx(:,1) = scl * vecx(:,1)
  end subroutine wt_normalize_one_atomic_solution
  !
  subroutine wt_atomic_solution_improve(verbose,lval,eval,evec)
    integer(ik), intent(in)    :: verbose     ! Reporting level
    integer(ik), intent(in)    :: lval        ! Desired angular momentum
    complex(rk), intent(inout) :: eval(:)     ! Eigenvalues of the Hamiltonian matrix block for L=lval
    complex(rk), intent(inout) :: evec(:,:,:) ! Left (:,:,1) and right (:,:,2) eigenvectors  
    !
    integer(ik) :: ipass
    integer(ik) :: iev
    complex(rk) :: eval_initial(size(eval))
    !
    call TimerStart('Atomic solutions: Improve')
    eval_initial = eval
    improvement_passes: do ipass=1,1
      scan_eigenvalues: do iev=1,sd_nradial
        call wt_improve_one_atomic_eigenstate(verbose,lval,eval(iev),evec(:,iev,:))
      end do scan_eigenvalues
      call wt_atomic_solution_biorthogonalize(evec)
    end do improvement_passes
    !
    if (verbose>=1) then
      eval_initial = eval_initial - eval
      iev = maxloc(abs(eval_initial),dim=1)
      write (out,"(/'Iterative update of L=',i0,' solutions complete')") lval
      write (out,"('        ground-state energy is ',g34.22,1x,g34.22,' change = ',g18.7,1x,g18.7)") eval(1), eval_initial(1)
      write (out,"('   most affected eigenvalue is ',g34.22,1x,g34.22,' change = ',g18.7,1x,g18.7)") eval(iev), eval_initial(iev)
    end if
    !
    call TimerStop('Atomic solutions: Improve')
  end subroutine wt_atomic_solution_improve
  !
  !  Improve a single eigenstate using inverse iterations
  !
  subroutine wt_improve_one_atomic_eigenstate(verbose,lval,eval,evec)
    integer(ik), intent(in)    :: verbose   ! Reporting level
    integer(ik), intent(in)    :: lval      ! Desired angular momentum
    complex(rk), intent(inout) :: eval      ! Eigenvalue to improve
    complex(rk), intent(inout) :: evec(:,:) ! Left (:,1) and right (:,2) eigenvectors to improve
    !
    integer(ik) :: iter_evec, iter_eval
    logical     :: fail(2)
    real(rk)    :: mn   (3,sd_nradial), mt   (3,sd_nradial)
    complex(rk) :: leqnf(3,sd_nradial), leqtf(3,sd_nradial)
    complex(rk) :: tmp (sd_nradial)
    complex(rk) :: vecx(sd_nradial,2)
    complex(rk) :: delta_m1  ! Correction to the eigenvalue^{-1}
    !
    improve_eval: do iter_eval=1,5
      !
      call build_right_linear_system(eval,leqnf,mn,fail(2))
      call build_left_linear_system (eval,leqtf,mt,fail(1))
      if (any(fail).and.(iter_eval<=1 .or. verbose>2)) then
        write (out,"('WARNING: wt_improve_one_atomic_eigenstate: Update iteration ',i0,' failed, leave solutions unchanged')") &
               iter_eval
        return
      end if
      !
      !  There is non-negligible cost to building linear system solutions;
      !  perform a few eigenvector updates before updating the eigenvalue
      !
      improve_evec: do iter_evec=1,3
        ! Inverse iteration: right vector
        call m3d_multiply(mn,evec(:,2),tmp)
        call m3d_solve(leqnf,tmp,vecx(:,2))
        ! Inverse iteration: left vector
        call m3d_solve(leqtf,evec(:,1),tmp)
        call m3d_multiply(mt,tmp,vecx(:,1))
        ! Compute correction to the eigenvalue
        delta_m1  = sum(evec(:,1)*vecx(:,2))
        if (verbose>2) then
          write (out,"('Update iteration ',i0,',',i0,' eval = ',2(g32.20,1x),' correction = ',2(g32.20,1x))") &
                 iter_eval, iter_evec, eval, 1/delta_m1
        end if
        call wt_normalize_one_atomic_solution(vecx)
        ! Store the updated eigenvectors
        evec(:,:) = vecx(:,:)
      end do improve_evec
      ! Update eigenvalue
      eval = eval + 1/delta_m1
    end do improve_eval
    !
    contains
    subroutine build_right_linear_system(eps,leqnf,mn,fail)
      complex(rk), intent(in)  :: eps        ! Current eigenvalue guess
      complex(rk), intent(out) :: leqnf(:,:) ! Factorization of the linear system
      real(rk), intent(out)    :: mn   (:,:) ! Scaling factor for the RHS
      logical, intent(out)     :: fail
      !
      complex(rk) :: vtmp(  sd_nradial)
      complex(rk) :: leqn(3,sd_nradial)
      !
      !  Prepare linear system matrices for solving right linear system: 
      !
      !     [-(1/(2m))M^{-1}Delta+V-eps] Y = Z.
      !
      !  This is done as:
      !
      !     [-(1/(2m))Delta+M(V-eps)] Y = M Z
      !
      vtmp(:) = sd_pottab(:,1,lval) - eps
      if (sd_capped) then
        vtmp(sd_cap_start:) = vtmp(sd_cap_start:) + sd_captab(sd_cap_start:)
      end if
      if (lval==0) then
        call m3d_right_scale(sd_m2n_l0,vtmp,leqn)
        leqn(:,:) = leqn + (-1._rk/(2._rk*electron_mass))*sd_d2n_l0
        mn  (:,:) = sd_m2n_l0
      else
        call m3d_right_scale(sd_m2n_lx,vtmp,leqn)
        leqn(:,:) = leqn + (-1._rk/(2._rk*electron_mass))*sd_d2n_lx
        mn  (:,:) = sd_m2n_lx
      end if
      call m3d_decompose(leqn,leqnf,fail)
    end subroutine build_right_linear_system
    !
    subroutine build_left_linear_system(eps,leqtf,mt,fail)
      complex(rk), intent(in)  :: eps        ! Current eigenvalue guess
      complex(rk), intent(out) :: leqtf(:,:) ! Factorization of the linear system
      real(rk), intent(out)    :: mt   (:,:) ! Scaling factor for the solution
      logical, intent(out)     :: fail
      !
      complex(rk) :: vtmp(  sd_nradial)
      complex(rk) :: leqt(3,sd_nradial)
      !
      !  Prepare linear system matrices for solving left linear system: 
      !
      !     [-(1/(2m))Delta^T M^{-T}+V-eps] Y = Z.
      !
      !  This is done as:
      !
      !     [-(1/(2m))Delta^T+(V-eps)M^T] Q = Z
      !     Y = M^T Q
      !
      vtmp(:) = sd_pottab(:,1,lval) - eps
      if (sd_capped) then
        vtmp(sd_cap_start:) = vtmp(sd_cap_start:) + sd_captab(sd_cap_start:)
      end if
      if (lval==0) then
        call m3d_left_scale(vtmp,sd_m2t_l0,leqt)
        leqt(:,:) = leqt(:,:) + (-1._rk/(2._rk*electron_mass))*sd_d2t_l0
        mt  (:,:) = sd_m2t_l0
      else
        call m3d_left_scale(vtmp,sd_m2t_lx,leqt)
        leqt(:,:) = leqt(:,:) + (-1._rk/(2._rk*electron_mass))*sd_d2t_lx
        mt  (:,:) = sd_m2t_lx
      end if
      call m3d_decompose(leqt,leqtf,fail)
    end subroutine build_left_linear_system
  end subroutine wt_improve_one_atomic_eigenstate
  !
  subroutine wt_atomic_solution_verify(lval,evec)
    integer(ik), intent(in) :: lval
    complex(rk), intent(in) :: evec(:,:,:) ! Left (:,:,1) and right (:,:,2) eigenvectors  
    !
    integer(ik)              :: imo, jmo, alloc
    real(rk)                 :: nrm
    complex(rk), allocatable :: norm(:,:)
    !
    call TimerStart('Atomic solutions: Verify')
    allocate (norm(sd_nradial,sd_nradial),stat=alloc)
    if (alloc/=0) stop 'wavefunction_tools%wt_atomic_solution_verify - allocation failed'
    !
    norm = matmul(transpose(evec(:,:,1)),evec(:,:,2))
    test_norm_j: do jmo=1,sd_nradial
      test_norm_i: do imo=1,sd_nradial
        nrm = 0
        if (imo==jmo) nrm = 1
        if (abs(norm(imo,jmo)-nrm)>1e5_rk*spacing(1._rk)) then
          write (out,"('WARNING: _GEEV eigenvectors ',i0,',',i0,' L=',i0,' are not biorthogonal. Error =',2(1x,g24.13))") &
                 imo, jmo, lval, norm(imo,jmo)-nrm
        end if
      end do test_norm_i
    end do test_norm_j
    deallocate (norm)
    call TimerStop('Atomic solutions: Verify')
  end subroutine wt_atomic_solution_verify
  !
  !  Reset wavefunction norm.
  !  We will set Cartesian norm of the right wavefunction to 1, while maintaining its phase.
  !  The left wavefunction will be adjusted so that <l|r> product is also 1.
  !  Keep in mind that the left wavefunction does not need to be conjugated!
  !
  subroutine wt_normalize(wfn_l,wfn_r,norm)
    type(sd_wfn), intent(inout) :: wfn_l   ! Left wavefunction
    type(sd_wfn), intent(inout) :: wfn_r   ! Right wavefunction
                                           ! Since we are dealing with (potentially) non-Hermitian
                                           ! operators, our wavefunctions always come as a pair of
                                           ! the left and right vectors.
    complex(rk), intent(out)    :: norm(2) ! norm(1) - Input wavefunction norm <l|r> 
                                           ! norm(2) - Cartesian norm of the input right wavefunction <r|r>
                                           ! norm(2) is guaranteed to be real.
    !
    integer(ik) :: lval, mval, sval
    real(rk)    :: sum_rr, scl_r
    complex(rk) :: sum_lr, scl_l
    !
    call TimerStart('WF normalize')
    sum_rr = 0
    sum_lr = 0
    !$omp parallel default(none) &
    !$omp& shared(sd_mmin,sd_mmax,sd_lmax,sd_nspin,wfn_l,wfn_r) &
    !$omp& private(mval,lval,sval) &
    !$omp& shared(sum_rr,sum_lr,scl_r,scl_l)
    sense_loop_m: do mval=sd_mmin,sd_mmax
      !$omp do reduction(+:sum_rr,sum_lr)
      sense_loop_l: do lval=abs(mval),sd_lmax
        sense_loop_s: do sval=1,sd_nspin
          sum_rr = sum_rr + sum(abs(wfn_r%wfn(:,sval,lval,mval))**2)
          sum_lr = sum_lr + sum(wfn_l%wfn(:,sval,lval,mval)*wfn_r%wfn(:,sval,lval,mval))
        end do sense_loop_s
      end do sense_loop_l
      !$omp end do nowait
    end do sense_loop_m
    !$omp barrier
    !$omp single
    scl_r = 1._rk / sqrt(sum_rr)
    scl_l = sqrt(sum_rr) / sum_lr
    !$omp end single
    scale_loop_m: do mval=sd_mmin,sd_mmax
      !$omp do
      scale_loop_l: do lval=abs(mval),sd_lmax
        scale_loop_s: do sval=1,sd_nspin
          wfn_r%wfn(:,sval,lval,mval) = scl_r*wfn_r%wfn(:,sval,lval,mval)
          wfn_l%wfn(:,sval,lval,mval) = scl_l*wfn_l%wfn(:,sval,lval,mval)
        end do scale_loop_s
      end do scale_loop_l
      !$omp end do nowait
    end do scale_loop_m
    !$omp end parallel
    norm(1) = sum_lr
    norm(2) = sum_rr
    call TimerStop('WF normalize')
  end subroutine wt_normalize
  !
  !  Calculate expectation value of the Hamiltonian
  !
  subroutine wt_energy(wfn_l,wfn_r,apot,energy,norm)
    type(sd_wfn), intent(in) :: wfn_l     ! Left wavefunction 
    type(sd_wfn), intent(in) :: wfn_r     ! Right wavefunction 
    real(xk), intent(in)     :: apot      ! Current value of the vector-potential; assumed to be along Z
    complex(rk), intent(out) :: energy(2) ! Expectation value of the total energy;
                                          ! [1] = including all terms in the Hamiltonian
                                          ! [2] = excluding the CAP term; only collected if sd_capped is .true.
    complex(rk), intent(out) :: norm      ! Norm of the input wavefunction; <psi|psi>
    !
    integer(ik) :: lval, mval
    complex(rk) :: hpsi(sd_nradial), tmp(sd_nradial)
    complex(rk) :: ap_energy
    real(rk)    :: vrminus(sd_nradial)
    real(rk)    :: a2_potential, a_factor
    !
    call TimerStart('WF energy')
    !
    if (sd_nspin/=1) stop 'wavefunction_tools%wt_energy - spinorbit not implemented'
    !
    !  Coupling coefficients and prefactors due to the vector potential
    !
    a_factor     = real((electron_charge/electron_mass) * apot,kind=rk)
    a2_potential = real(0.5_rk * (electron_charge**2 / electron_mass) * apot**2,kind=rk)
    vrminus(:)   = 1._rk / sd_rtab(1:sd_nradial)
    !
    energy = 0
    norm = 0
    !$omp parallel default(none) &
    !$omp& shared(sd_mmin,sd_mmax,sd_lmax,sd_d2n_l0,sd_m2nf_l0,sd_d2n_lx,sd_m2nf_lx) &
    !$omp& shared(sd_pottab,sd_capped,sd_captab,sd_cap_start) &
    !$omp& shared(a2_potential,a_factor,vrminus,sd_d1n_l0,sd_m1nf_l0,sd_d1n_lx,sd_m1nf_lx) &
    !$omp& shared(wfn_l,wfn_r) &
    !$omp& private(mval,lval,tmp,hpsi,ap_energy) &
    !$omp& reduction(+:energy,norm)
    energy = 0
    norm = 0
    sense_loop_m: do mval=sd_mmin,sd_mmax
      !$omp do
      sense_loop_l: do lval=abs(mval),sd_lmax
        !
        !  L-diagonal part of the Hamiltonian: field-free and A^2 terms
        !
        if (lval==0) then
          call m3d_multiply(sd_d2n_l0,wfn_r%wfn(:,1,lval,mval),tmp)
          call m3d_solve(sd_m2nf_l0,tmp,hpsi)
        else
          call m3d_multiply(sd_d2n_lx,wfn_r%wfn(:,1,lval,mval),tmp)
          call m3d_solve(sd_m2nf_lx,tmp,hpsi)
        end if
        hpsi = (-0.5_rk/electron_mass)*hpsi
        hpsi = hpsi + (a2_potential+sd_pottab(:,1,lval))*wfn_r%wfn(:,1,lval,mval)
        !
        if (sd_capped) then
          energy(2) = energy(2) + sum(wfn_l%wfn(:,1,lval,mval)*hpsi)
          hpsi(sd_cap_start:) = hpsi(sd_cap_start:) &
                              + sd_captab(sd_cap_start:)*wfn_r%wfn(sd_cap_start:,1,lval,mval)
        end if
        energy(1) = energy(1) + sum(wfn_l%wfn(:,1,lval,mval)*hpsi)
        norm      = norm   + sum(wfn_l%wfn(:,1,lval,mval)*wfn_r%wfn(:,1,lval,mval))
        !
        !  A.p terms; these only need to be calculated if vector potential is not zero
        !
        if (a_factor==0) cycle sense_loop_l
        !
        !  Act on the RHS with the radial-gradient part of Hmix
        !
        if (lval==0) then
          call m3d_multiply(sd_d1n_l0,wfn_r%wfn(:,1,lval,mval),tmp)
          call m3d_solve(sd_m1nf_l0,tmp,hpsi)
        else
          call m3d_multiply(sd_d1n_lx,wfn_r%wfn(:,1,lval,mval),tmp)
          call m3d_solve(sd_m1nf_lx,tmp,hpsi)
        end if
        !
        !  Act of the RHS with the (1/r) term in Hang
        !
        tmp = vrminus(:)*wfn_r%wfn(:,1,lval,mval)
        !
        !  Assemble the field-coupling term
        !
        ap_energy = 0
        if (lval>abs(mval)) then ! We can couple down
          ap_energy = ap_energy + real(sd_laser_clm(lval,  mval),kind=rk)*sum(wfn_l%wfn(:,1,lval-1,mval)*(hpsi+(lval  )*tmp))
        endif
        if (lval<sd_lmax) then ! We can couple up
          ap_energy = ap_energy + real(sd_laser_clm(lval+1,mval),kind=rk)*sum(wfn_l%wfn(:,1,lval+1,mval)*(hpsi-(lval+1)*tmp))
        end if
        energy = energy + (0,1)*a_factor*ap_energy
      end do sense_loop_l
      !$omp end do nowait
    end do sense_loop_m
    !$omp barrier
    !$omp end parallel
    call TimerStop('WF energy')
  end subroutine wt_energy
  !
  !  Calculate dipole expectation and dipole acceleration 
  !
  !  WARNING: Dipole acceleration is valid ONLY for multiplicative potentials
  !
  !  The calculated acceleration term does not include the time derivative of the
  !  vector potential; this term is best evaluated later, in the laboratory frame.
  !  Specifically, one needs to add:
  !
  !      -(e^2/m) (d A/d t) <L|R> = (e^2/m) E <L|R>
  !
  !  where E is the instantaneous value of the laser electric field.
  !
  subroutine wt_dipole(wfn_l,wfn_r,dipole,acceleration)
    type(sd_wfn), intent(in) :: wfn_l           ! Left wavefunction 
    type(sd_wfn), intent(in) :: wfn_r           ! Right wavefunction 
    complex(rk), intent(out) :: dipole(3)       ! <L|q r|R> expectation value
    complex(rk), intent(out) :: acceleration(3) ! [(d^2/d t^2) <L|q r|R>] + (e^2/m) (d A/d t) <L|R>
    !
    integer(ik) :: l_left, m_left, l_right, m_right, m_op, ispin
    complex(rk) :: dip_sph(-1:1), acc_sph(-1:1), wgt_dip, wgt_acc
    real(rk)    :: wgt_ang
    !
    call TimerStart('WF dipole')
    !
    !  Compute spherical-tensor form for the dipole
    !
    dip_sph(:) = 0
    acc_sph(:) = 0
    !$omp parallel default(none) &
    !$omp& shared(sd_mmin,sd_mmax,sd_lmax,sd_nradial,sd_nspin,wfn_l,wfn_r,sd_rtab,sd_dvdrtab) &
    !$omp& private(m_left,l_left,m_right,l_right,m_op,ispin,wgt_ang,wgt_dip,wgt_acc) &
    !$omp& reduction(+:dip_sph,acc_sph)
    dip_sph(:) = 0
    acc_sph(:) = 0
    loop_m_left: do m_left=sd_mmin,sd_mmax
      !$omp do
      loop_l_left: do l_left=abs(m_left),sd_lmax
        loop_m_right: do m_right=max(sd_mmin,m_left-1),min(sd_mmax,m_left+1)
          loop_l_right: do l_right=l_left-1,l_left+1
            if (l_right<abs(m_right)) cycle loop_l_right
            if (l_right>sd_lmax) cycle loop_l_right
            m_op    = m_left - m_right
            wgt_ang = angular_term(l_left,m_left,1_ik,m_op,l_right,m_right)
            wgt_dip = 0
            wgt_acc = 0
            loop_spin: do ispin=1,sd_nspin
              wgt_dip = wgt_dip + sum(wfn_l%wfn(:,ispin,l_left,m_left)*wfn_r%wfn(:,ispin,l_right,m_right)*sd_rtab(1:sd_nradial))
              wgt_acc = wgt_acc + sum(wfn_l%wfn(:,ispin,l_left,m_left)*wfn_r%wfn(:,ispin,l_right,m_right)*sd_dvdrtab)
            end do loop_spin
            dip_sph(m_op) = dip_sph(m_op) + electron_charge*wgt_ang*wgt_dip
            acc_sph(m_op) = acc_sph(m_op) - electron_charge*wgt_ang*wgt_acc/electron_mass
          end do loop_l_right
        end do loop_m_right
      end do loop_l_left
      !$omp end do nowait
    end do loop_m_left
    !$omp end parallel
    !
    !  Transform to Cartesian coordinates
    !
    call sph2cart(dip_sph,dipole)
    call sph2cart(acc_sph,acceleration)
    call TimerStop('WF dipole')
    !
    contains
    real(rk) function angular_term(l1,m1,l2,m2,l3,m3)
      integer(ik), intent(in) :: l1, m1, l2, m2, l3, m3
      ! Brute-force expression using 3J symbols; Mathematica phase convenstions
      angular_term = ((-1)**m1) * sqrt((2*l1+1)*(2*l2+1)*(2*l3+1)/(4*pi)) &
                   * Math3J(l1,l2,l3,0_ik,0_ik,0_ik) * Math3J(l1,l2,l3,-m1,m2,m3)
    end function angular_term
    subroutine sph2cart(sph,cart)
      complex(rk), intent(in)  :: sph (-1:1) ! Spherical tensor
      complex(rk), intent(out) :: cart(3)    ! Equivalent Cartesian tensor
      !
      cart(1) =       sqrt((2*pi)/3) * (sph(-1)-sph(+1))
      cart(2) = (0,1)*sqrt((2*pi)/3) * (sph(-1)+sph(+1))
      cart(3) =       sqrt((4*pi)/3) *  sph( 0)
    end subroutine sph2cart
  end subroutine wt_dipole
  !
end module wavefunction_tools
