!
!  Analysis aid: projects wavefunction on the field-free solutions,
!  and reports the results.
!
module composition_analysis
  use accuracy
  use constants
  use timer
  use spherical_data
  use wavefunction_tools
  implicit none
  private
  public ca_analyze
  !
  contains
  !
  !  Calculation of the field-free eigenspectrum
  !
  subroutine ca_analyze(verbose,threshold,wfn_l,wfn_r)
    integer(ik), intent(in)  :: verbose   ! Debugging level
    real(rk), intent(in)     :: threshold ! Reporting threshold for the amplitudes of the field-free solutions
    type(sd_wfn), intent(in) :: wfn_l     ! Left/right wavefunction pair to analyze
    type(sd_wfn), intent(in) :: wfn_r
    !
    integer(ik)              :: mval, lval, ispin, iev, mmin, mmax
    integer(ik)              :: alloc
    complex(rk), allocatable :: block_evec(:,:,:)         ! Eigenvectors (unused here)
    complex(rk)              :: block_eval(sd_nradial)    ! Eigenvalues
    complex(rk), allocatable :: amp(:,:,:,:,:), en(:,:)   ! Energies and amplitudes
    complex(rk)              :: wgt
    complex(rk)              :: pop_total, pop_bound
    complex(rk)              :: lm_norms(sd_mmin:sd_mmax,0:sd_lmax)
    !
    call TimerStart('Field-free analysis')
    !
    !  Begin by calculating norms of each L,M channel contribution; this is useful check
    !  on the quality of the analysis later.
    !
    !$omp parallel do default(none) &
    !$omp& private(lval,mmin,mmax,mval,ispin) &
    !$omp& shared(sd_lmax,sd_mmin,sd_mmax,sd_nspin,lm_norms,wfn_l,wfn_r)
    lm_norm_l: do lval=0,sd_lmax
      mmin = max(sd_mmin,-lval)
      mmax = min(lval,sd_mmax)
      lm_norm_m: do mval=mmin,mmax
        lm_norms(mval,lval) = 0
        lm_norm_s: do ispin=1,sd_nspin
          lm_norms(mval,lval) = lm_norms(mval,lval) + sum(wfn_l%wfn(:,ispin,lval,mval)*wfn_r%wfn(:,ispin,lval,mval))
        end do lm_norm_s
      end do lm_norm_m
    end do lm_norm_l
    !$omp end parallel do
    !
    write (out,"(/t5,'Final norm, by total angular momentum'/)")
    write (out,"(1x,a3,1x,a4,2x,a34,2x,a34)") ' L ', '  M ', '      Re(norm)    ', '     Im(norm)     ', &
                                              '---', ' ---', '------------------', '------------------'
    lm_print_l: do lval=0,sd_lmax
      mmin = max(sd_mmin,-lval)
      mmax = min(lval,sd_mmax)
      lm_print_m: do mval=mmin,mmax
        write (out,"(1x,i3,1x,i4,2x,g34.23e3,2x,g34.23e3)") lval, mval, lm_norms(mval,lval)
      end do lm_print_m
    end do lm_print_l
    write (out,"()")
    call flush_wrapper(out)
    !
    !  The rest is expensive; skip it if threshold is negative
    !
    if (threshold<0) then
      call TimerStop('Field-free analysis')
      return
    end if
    write (out,"(/t5,'Analyzing the final wavefunction in terms of field-free eigenstates'/)")
    !
    allocate (amp(sd_nradial,2,sd_nspin,sd_mmin:sd_mmax,0:sd_lmax),en(sd_nradial,0:sd_lmax),stat=alloc)
    if (alloc/=0) then
      write (out,"('ca_analyze: No memory for analysis. Error code ',i0)") alloc
      stop 'composition_analysis%ca_analyze - no memory for analysis'
    end if
    !$omp parallel default(none) &
    !$omp& private(lval,mval,mmin,mmax,ispin,block_evec,block_eval,alloc) &
    !$omp& shared(sd_nradial,sd_lmax,sd_mmin,sd_mmax,sd_nspin,verbose) &
    !$omp& shared(wfn_l,wfn_r,amp,en)
    allocate (block_evec(sd_nradial,sd_nradial,2),stat=alloc)
    if (alloc/=0) then
      write (out,"('ca_analyze: No per-thread memory for analysis. Error code ',i0)") alloc
      stop 'composition_analysis%ca_analyze - no per-thread memory for analysis'
    end if
    !$omp do schedule(guided)
    scan_l_channels: do lval=0,sd_lmax
      call wt_atomic_solutions(verbose,lval,block_eval,block_evec)
      en(:,lval) = block_eval(:)
      mmin = max(sd_mmin,-lval)
      mmax = min(lval,sd_mmax)
      scan_spin: do ispin=1,sd_nspin
        amp(:,1,ispin,mmin:mmax,lval) = matmul(transpose(block_evec(:,:,1)),wfn_r%wfn(:,ispin,lval,mmin:mmax))
        amp(:,2,ispin,mmin:mmax,lval) = matmul(transpose(block_evec(:,:,2)),wfn_l%wfn(:,ispin,lval,mmin:mmax))
      end do scan_spin
    end do scan_l_channels
    !$omp end do
    deallocate (block_evec)
    !$omp end parallel
    !
    !  Reporting part, in increasingly excruciating detail
    !
    write (out,"(/t5,'Final populations, by total angular momentum'/)")
    write (out,"((1x,a3,3(2x,a24,1x,a24)))") &
           ' L ', ' Total population ', ' ', ' Bound population ', ' ', ' Continuum population ', ' ', &
           '---', '------------------', ' ', '------------------', ' ', '----------------------', ' '
    l_resolved_l_channels: do lval=0,sd_lmax
      pop_total = 0
      pop_bound = 0
      !$omp parallel do default(none) &
      !$omp& private(mval,iev,wgt) reduction(+:pop_total,pop_bound) &
      !$omp& shared(sd_mmin,sd_mmax,sd_nradial,lval,amp,en)
      l_resolved_m_channels: do mval=max(sd_mmin,-lval),min(lval,sd_mmax)
        l_resolved_evs: do iev=1,sd_nradial
          ! wgt = abs(sum(amp(iev,1,:,mval,lval)*amp(iev,2,:,mval,lval)))
          wgt = sum(amp(iev,1,:,mval,lval)*amp(iev,2,:,mval,lval))
          pop_total = pop_total + wgt
          if (real(en(iev,lval),kind=rk)<0) pop_bound = pop_bound + wgt
        end do l_resolved_evs
      end do l_resolved_m_channels
      !$omp end parallel do
      write (out,"((1x,i3,3(2x,g24.13e3,1x,g24.13e3)))") lval, pop_total, pop_bound, pop_total-pop_bound
    end do l_resolved_l_channels
    write (out,"()")
    !
    write (out,"(/t5,'Final populations, by total angular momentum and angular momentum projection'/)")
    write (out,"((1x,a3,1x,a4,3(2x,a24,1x,a24)))") &
           ' L ', ' M ', ' Total population ', ' ', ' Bound population ', ' ', ' Continuum population ', ' ', &
           '---', '---', '------------------', ' ', '------------------', ' ', '----------------------', ' '
    lm_resolved_l_channels: do lval=0,sd_lmax
      lm_resolved_m_channels: do mval=max(sd_mmin,-lval),min(lval,sd_mmax)
        pop_total = 0
        pop_bound = 0
        !$omp parallel do default(none) &
        !$omp& private(iev,wgt) reduction(+:pop_total,pop_bound) &
        !$omp& shared(mval,sd_nradial,lval,amp,en)
        lm_resolved_evs: do iev=1,sd_nradial
          ! wgt = abs(sum(amp(iev,1,:,mval,lval)*amp(iev,2,:,mval,lval)))
          wgt = sum(amp(iev,1,:,mval,lval)*amp(iev,2,:,mval,lval))
          pop_total = pop_total + wgt
          if (real(en(iev,lval),kind=rk)<0) pop_bound = pop_bound + wgt
        end do lm_resolved_evs
        !$omp end parallel do
        write (out,"((1x,i3,1x,i4,3(2x,g24.13e3,1x,g24.13e3)))") lval, mval, pop_total, pop_bound, pop_total-pop_bound
      end do lm_resolved_m_channels
    end do lm_resolved_l_channels
    write (out,"()")
    !
    write (out,"(/t5,'Large amplitudes of individual field-free states'/)")
    write (out,"(1x,a3,1x,a4,1x,a5,2x,a24,1x,a24,2x,a24,1x,a24,2x,a24,1x,a24,2x,a24,1x,a24)") &
      ' L ', ' M ', ' I ', ' Re[E(i)], H ', ' Im[E(i)] ', '  Re[Wgt]  ', '  Im[Wgt]  ', &
                            ' Re[<I|W>] ', ' Im[<I|W>] ', ' Re[<W|I>] ', ' Im[<W|I>] ', &
      '---', '---', '---', '-------------', '----------', '-----------', '-----------', &
                              '-----------', '-----------', '-----------', '---------'
    state_resolved_l_channels: do lval=0,sd_lmax
      state_resolved_m_channels: do mval=max(sd_mmin,-lval),min(lval,sd_mmax)
        state_resolved_evs: do iev=1,sd_nradial
          if (all(abs(amp(iev,:,:,mval,lval))<threshold)) cycle state_resolved_evs
          ! wgt = abs(sum(amp(iev,1,:,mval,lval)*amp(iev,2,:,mval,lval)))
          wgt = sum(amp(iev,1,:,mval,lval)*amp(iev,2,:,mval,lval))
          write (out,"(1x,i3,1x,i4,1x,i5,2(2x,g24.13e3,1x,g24.13e3),2x,(t120,2(g24.13e3,1x,g24.13e3)))") &
                 lval, mval, iev, en(iev,lval), wgt, amp(iev,:,:,mval,lval)
        end do state_resolved_evs
      end do state_resolved_m_channels
    end do state_resolved_l_channels
    write (out,"()")
    !
    deallocate (amp,en)
    call flush_wrapper(out)
    !
    call TimerStop('Field-free analysis')
  end subroutine ca_analyze
end module composition_analysis
