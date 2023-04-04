!
!  Library of simple atomic potentials
!  These routines should not be used within a critical path; if their output is
!  needed repeatedly, then cache it.
!
module potential_tools
  use accuracy
  use constants
  use timer
  implicit none
  private
  public pot_name, pot_param
  public pt_initialize
  public pt_evaluate_potential
  !
  !  We recognize
  !
  character(len=clen), save :: pot_name     = 'hydrogenic'   ! Name of the potential; can be one of:
                                                             ! 'hydrogenic' - Z/r
                                                             ! 'harmonic'   - (1/2)*k*r**2
                                                             ! 'argon'      - Effective potential from HGM PRA 60, 1341 (1999)
                                                             ! 'argon 2P'   - Valence-only effective potential. This is NOT the
                                                             !                same potential given in HGM; I can't reproduce
                                                             !                any results claimed for that potential.
  real(rk), save            :: pot_param(5) = 1._rk          ! Parameters of the potential
  !
  contains
  !
  subroutine pt_initialize
    !
    select case (pot_name)
      case default
        write (out,"('potential_tools%pt_setup_potential: potential ',a,' is not recognized')") trim(pot_name)
        stop 'potential_tools%pt_setup_potential - unknown potential'
      case ('hydrogenic')
        write (out,"('Coulombic nuclear potential. Nuclear charge = ',g24.13)") pot_param(1)
      case ('harmonic')
        write (out,"('Harmonic nuclear potential. Force constant = ',g24.13)") pot_param(1)
      case ('argon')
        write (out,"('""Argon"" potential from H.G. Muller, PRA 60, 1341 (1999)')")
      case ('argon 2P')
        write (out,"('""Argon"" valence-only potential, inspired by H.G. Muller, PRA 60, 1341 (1999).')")
        write (out,"('WARNING: This is NOT the same potential as HGM Eq. 2!')")
      case ('scp')
        write (out,"('lamda = ',g24.13)") pot_param(1)
    end select 
  end subroutine pt_initialize
  !
  function pt_evaluate_potential(l,j,r,centrifugal) result(v)
    integer(ik), intent(in) :: l           ! Orbital angular momentum channel
    integer(ik), intent(in) :: j           ! Total angular momentum sub-channel (currently unused)
                                           ! -1: l-1/2 sub-channel
                                           !  0: weighed-average of the sub-channels
                                           ! +1: l+1/2 sub-channel
    real(rk), intent(in)    :: r           ! Radial position where potential is needed
    logical, intent(in)     :: centrifugal ! Add centrifugal term
    real(rk)                :: v           ! Potential at a grid point, including the centrifugal term
    !
    real(rk) :: zeff, arg
    !
    !  Sanity check: L must be non-negative; J must be in the +1,-1 range, and R must be positive
    !
    if (l<0 .or. abs(j)>1 .or. r<=0._rk) then
      stop 'potential_tools%pt_evaluate_potential - domain error'
    end if
    !
    !  First, the system-specific part of the potential
    !
    select case (pot_name)
      case default
        write (out,"('potential_tools%pt_evaluate_potential: potential ',a,' is not recognized')") trim(pot_name)
        stop 'potential_tools%pt_evaluate_potential - unknown potential'
      case ('hydrogenic')
        v = electron_charge*pot_param(1)/r
      case ('harmonic')
        v = 0.5_rk*pot_param(1)*r**2
      case ('argon')
        zeff = 1._rk + 5.4_rk*exp(-r) + (17._rk-5.4_rk)*exp(-3.682_rk*r)
        v    = electron_charge*zeff/r
      case ('argon 2P')
        zeff = 1._rk + 7.625195_rk*exp(-1.02557_rk*r)
        arg  = 10.00_rk*(r-0.37110_rk)
        if (arg<0.9_rk*log(huge(1._rk))) then
          zeff = zeff - 124.55_rk/(1._rk + exp(arg))
        end if
        v = electron_charge*zeff/r
      case ('scp') 
        v = electron_charge*exp(-pot_param(1)*r)/r
    end select 
    ! 
    !  Add the universal centrifugal term. Nominally, that's a kinetic energy contribution, but ...
    !
    if (centrifugal) then
      v = v + (1._rk/(2._rk*electron_mass)) * l*(l+1) / r**2
    end if
  end function pt_evaluate_potential
  !
end module potential_tools
