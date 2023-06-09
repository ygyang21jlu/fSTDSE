!
!  Vector-potential shape.
!  This module is never on the critical path. 
!  It is safe to call from a parallel region, but probably should not be.
!  This module operates in at least double precision, even if propagation
!  is in single precision.
!
!  We define vector-potential as a direction and a signed magnitude in the spherical 
!  polar coordinate system, consisting of:
!
!    - Signed magnitude of the vector-potential
!    - polar angle theta (angle between the direction of the VP and lab Z axis)
!    - polar angle phi (angle between the X axis and XY projection of the VP)
!
!  This is the choice of parameters most convenient for HGM-style propagation.
!  If your VP is defined in terms of some other set of parameters, please convert them
!  to apot/theta/phi here. Do not try to hack downstream routines - that way lies madness!
!
!  (partial) List of vector-potentials supported here:
!
!  'z Gaussian'  - finite-time Gaussian envelope (see code), polarized along lab Z
!  'zy Gaussian' - finite-time Gaussian envelope, independent components along lab Z and Y
!  'xy Gaussian' - finite-time Gaussian envelope, independent components along lab Z and Y
!  'z Sin2'      - Sin^2 envelope, polarized along lab Z
!  'z Flat-Sin2' - Flat-top pulse with sin^2 raise/fall
!  'zy CW'       - CW field, circularly polarized in the ZY plane
!
!  Note that these are NOT the final angle values which will be used in the propagation:
!  we'll still need to perform the unwrap on the angles, to ensure that all of r, theta,
!  and phi are smooth in time.
!
module vectorpotential_tools
  use accuracy
  use constants
  implicit none
  private
  public vp_scale, vp_shape, vp_param, vp_param_x, vp_scale_x
  public vp_apot
  !
  !  List of vector-potential names. These are used in a couple different places; 
  !  I would prefer any typos to cause a compile-time error!
  ! 
  character(len=*), parameter :: vpn_zGaussian  = 'z Gaussian'
  character(len=*), parameter :: vpn_zSin2      = 'z Sin2'
  character(len=*), parameter :: vpn_zFlatSin2  = 'z Flat-Sin2'
  character(len=*), parameter :: vpn_zyGaussian = 'zy Gaussian'
  character(len=*), parameter :: vpn_xyGaussian = 'xy Gaussian'
  character(len=*), parameter :: vpn_zyCW       = 'zy CW'
  !
  !add by ygyang
  character(len=*), parameter :: vpn_SFA        = 'z SFA'
  !
  !
  real(xk), save            :: vp_scale   = 0.00_xk     ! Overall magnitude of the vector potential
  character(len=clen), save :: vp_shape   = ' '         ! Shape of the vector-potential. See vp_apot()
                                                        ! below for the list of possibilites. There is
                                                        ! a special magic value: vp_shape='table'. This
                                                        ! value is handled elsewhere.
  real(xk), save            :: vp_param(20)             ! Parameters with meaning depending on vp_shape;
                                                        ! First 10 values are intended for somewhat generic
                                                        ! parameters, which apply in more than one case.
                                                        ! Values from 11 on correspond to specific vp_shape choices.
                                                        ! See vp_apot() and below. 
  real(xk), save            :: vp_scale_x = 0.00_xk     ! Additional scaling parameter for the extra polatization component
  real(xk), save            :: vp_param_x(20)           ! Same as vp_param, for pulse shapes defined by two orthogonal
                                                        ! polarizations. ("x" stands for "extra", not the X axis!)
  !
  !  Parameters below this line are unpacked from vp_param() upon the first call to vp_apot()
  !
                                                        ! polarizations. ("x" stands for "extra", not the X axis!)
  !
  !  Parameters below this line are unpacked from vp_param() upon the first call to vp_apot()
  !
  logical, save             :: first      = .true.      ! Initialization is necessary
  !
  !  Values relevant for more than one vp_shape
  !
  real(xk), save            :: omega      = 0._xk       ! vp_param(1) : Carrier frequency [atomic units]
  real(xk), save            :: phase      = 0._xk       ! vp_param(2) : Carrier phase at pulse origin [radians]
  real(xk), save            :: origin     = 0._xk       ! vp_param(3) : Pulse origin [atomic units]
  real(xk), save            :: width      = 0._xk       ! vp_param(4) : Pulse width; meaning and units differ
  real(xk), save            :: x_omega    = 0._xk       ! vp_param_x(1) : Carrier frequency [atomic units]
  real(xk), save            :: x_phase    = 0._xk       ! vp_param_x(2) : Carrier phase at pulse origin [radians]
  real(xk), save            :: x_origin   = 0._xk       ! vp_param_x(3) : Pulse origin [atomic units]
  real(xk), save            :: x_width    = 0._xk       ! vp_param_x(4) : Pulse width; meaning and units differ
  !
  !  Values relevant for vp_shape=='Gaussian':
  !
  !  (width) corresponds to full width at half maximum of power, in atomic units of time
  !
  real(xk), save            :: gau_toff1    = 0._xk     ! vp_param(11) : Beginning of the hard turn-off, relative to origin
  real(xk), save            :: gau_toff2    = 0._xk     ! vp_param(12) : End of the hard turn-off, relative to origin
                                                        !                gau_toff2 must be greater than gau_toff1
  real(xk), save            :: gau_alpha    = 0._xk     ! Pulse width parameter; derived from width
  real(xk), save            :: x_gau_toff1  = 0._xk     ! vp_param_x(11) : Beginning of the hard turn-off, relative to origin
  real(xk), save            :: x_gau_toff2  = 0._xk     ! vp_param_x(12) : End of the hard turn-off, relative to origin
                                                        !                gau_toff2 must be greater than gau_toff1
  real(xk), save            :: x_gau_alpha  = 0._xk     ! Pulse width parameter; derived from width
  !
  !  Values relevant for vp_shape='Flat-Sin2'
  !
  real(xk), save            :: flat_raise   = 0._xk     ! vp_param(11) : Duration of the raising and falling edges, in 
                                                        !                atomic units of time
  !
  contains
  !
  function vp_apot(t,theta,phi) result(apot)
    real(xk), intent(in)            :: t     ! Time [atomic units] at which vector-potential is needed
    real(xk)                        :: apot  ! Vector-potential
    real(xk), intent(out), optional :: theta ! Polar angles
    real(xk), intent(out), optional :: phi   ! 
    !
    real(xk) :: th, ph  ! Local copies of theta and phi to simplify code logic
    real(xk) :: az, ay, ax
    !
    call initialize
    select case (vp_shape)
      case default
        write (out,"('Vector-potential shape ',a,' is not recognized')") trim(vp_shape)
        stop 'vectorpotential_tools%vp_apot - bad vp_shape'
      case (vpn_zyCW)
        apot = vp_scale
        th   = omega*t
        ph   = 0._xk
      case (vpn_zGaussian)
        apot = vp_scale * GaussianVP(t,omega,phase,origin,width,gau_toff1,gau_toff2,gau_alpha)
        th   = 0._xk
        ph   = 0._xk
      case (vpn_zSin2)
        apot = vp_scale * Sin2VP(t,omega,phase,origin,width)
        th   = 0._xk
        ph   = 0._xk
      case (vpn_zFlatSin2)
        apot = vp_scale * FlatSin2VP(t,omega,phase,origin,width,flat_raise)
        th   = 0._xk
        ph   = 0._xk
      case (vpn_zyGaussian)
        az = vp_scale   * GaussianVP(t,  omega,  phase,  origin,  width,  gau_toff1,  gau_toff2,  gau_alpha)
        ay = vp_scale_x * GaussianVP(t,x_omega,x_phase,x_origin,x_width,x_gau_toff1,x_gau_toff2,x_gau_alpha)
        apot = sqrt(az**2+ay**2)
        if (apot>0) then
          th = acos(az/apot)
        else
          th = 0._xk
        end if
        if (ay>=0) then
          ph   = 0._xk
        else
          ph   = pi_xk
        end if
      case (vpn_xyGaussian)
        ax = vp_scale   * GaussianVP(t,  omega,  phase,  origin,  width,  gau_toff1,  gau_toff2,  gau_alpha)
        ay = vp_scale_x * GaussianVP(t,x_omega,x_phase,x_origin,x_width,x_gau_toff1,x_gau_toff2,x_gau_alpha)
        apot = sqrt(ax**2+ay**2)
        if (apot>0) then
          th = 0.5_xk*pi_xk
        else
          th = 0._xk
        end if
        ph = atan2(ay,ax)
      !
      !add by ygyang
      case (vpn_SFA)
        apot = vp_scale * SFAVP(t,omega,phase,origin,width)
        th   = 0._xk
        ph   = 0._xk
      !
    end select
    if (present(theta)) theta = th
    if (present(phi)  ) phi   = ph
  end function vp_apot
  !
  subroutine initialize
    if (.not.first) return
    !
    !$omp critical
    if (first) then ! Flag may have been set by the time we acquired the lock
      select case (vp_shape)
        case default
          write (out,"('Vector-potential shape ',a,' is not recognized')") trim(vp_shape)
          stop 'vectorpotential_tools%initialize - bad vp_shape'
        case ('table')
          stop 'vectorpotential_tools%initialize - logic error: vp_shape=table must be handled elsewhere'
        case (vpn_zGaussian)
          call init_GaussianVP('along lab Z',vp_param,omega,phase,origin,width,gau_toff1,gau_toff2,gau_alpha)
        case (vpn_zSin2)
          call init_Sin2VP('along lab Z',vp_param,omega,phase,origin,width)
        case (vpn_zFlatSin2)
          call init_FlatSin2VP('along lab Z',vp_param,omega,phase,origin,width,flat_raise)
        case (vpn_zyGaussian)
          call init_GaussianVP('along lab Z',vp_param,    omega,  phase,  origin,  width,  gau_toff1,  gau_toff2,  gau_alpha)
          call init_GaussianVP('along lab Y',vp_param_x,x_omega,x_phase,x_origin,x_width,x_gau_toff1,x_gau_toff2,x_gau_alpha)
        case (vpn_xyGaussian)
          call init_GaussianVP('along lab X',vp_param,    omega,  phase,  origin,  width,  gau_toff1,  gau_toff2,  gau_alpha)
          call init_GaussianVP('along lab Y',vp_param_x,x_omega,x_phase,x_origin,x_width,x_gau_toff1,x_gau_toff2,x_gau_alpha)
        case (vpn_zyCW)
          omega = vp_param(1)
        case (vpn_SFA)
          call init_Sin2VP('along lab Z',vp_param,omega,phase,origin,width)
      end select
      first = .false.
      !$omp flush(first)
    end if
    !$omp end critical
  end subroutine initialize
  !
  !  Implementation of truncated-Gaussian vector-potential
  !
  function GaussianVP(time,omega,phase,origin,width,gau_toff1,gau_toff2,gau_alpha)
    real(xk), intent(in) :: time        ! Time
    real(xk)             :: GaussianVP  ! Vector-potential
    real(xk), intent(in) :: omega, phase, origin, width, gau_toff1, gau_toff2, gau_alpha
    !
    real(xk)             :: t_carrier   ! Time used for the carrier part
    real(xk)             :: t_envelope  ! Time used for the envelope part
    real(xk)             :: t_red       ! Time reduced into the [0:1] range for stretching
    real(xk)             :: exp_arg     ! Argument for the exponential envelope
    real(xk)             :: ph          ! Carrier phase at this time
    !
    GaussianVP = 0._xk
    t_carrier  = time - origin
    if (gau_toff2<=0._xk .or. abs(t_carrier)<=gau_toff2) then
      !
      !  The field is not zero; decide the effective envelope time
      !
      t_envelope = t_carrier
      if (gau_toff1<gau_toff2 .and. abs(t_envelope)>gau_toff1) then
        !
        !  We are in the envelope stretch zone; scale the envelope time.
        !  We do not need the sign for the envelope part, so do not bother with that
        !
        t_red      = (abs(t_envelope) - gau_toff1)/(gau_toff2-gau_toff1)
        t_red      = min(1._xk-spacing(10._xk),t_red)
        t_envelope = gau_toff1 + (2._xk/pi_xk) * (gau_toff2-gau_toff1) * tan((pi_xk/2._xk) * t_red)
      end if
      exp_arg    = gau_alpha*t_envelope**2
      ph         = omega * t_carrier + phase
      if (exp_arg<-log(spacing(1._xk))) then
        GaussianVP = exp(-exp_arg) * cos(ph)
      end if
    end if
  end function GaussianVP
  !
  subroutine init_GaussianVP(tag,vp_param,omega,phase,origin,width,gau_toff1,gau_toff2,gau_alpha)
    character(len=*), intent(in) :: tag
    real(xk), intent(in)         :: vp_param(:)
    real(xk), intent(out)        :: omega, phase, origin, width, gau_toff1, gau_toff2, gau_alpha
    !
    !  Unpack parameters from vp_param()
    !
    omega     = vp_param(1)
    phase     = vp_param(2)
    origin    = vp_param(3)
    width     = vp_param(4)
    gau_toff1 = vp_param(11)
    gau_toff2 = vp_param(12)
    !
    write (out,"(/'Using truncated Gaussian vector-potential ',a/)") trim(tag)
    write (out,"(t5,'         Carrier frequency = ',g24.13,' [(a.u. time)^-1]')") omega
    write (out,"(t5,'   Carrier phase at origin = ',g24.13,' [Radian]')") phase
    write (out,"(t5,'              Pulse origin = ',g24.13,' [a.u. time]')") origin
    write (out,"(t5,'  Full width at half power = ',g24.13,' [a.u. time]')") width
    write (out,"(t5,'Beginning of hard turn-off = ',g24.13,' [a.u. time]')") gau_toff1
    write (out,"(t5,'      End of hard turn-off = ',g24.13,' [a.u. time]')") gau_toff2
    !
    call prepare_alpha
    !
    contains
    !
    !  Calculate Gaussian envelope parameter from the given width.
    !  The desired value is the fixed point of a rapidly converging recursion.
    !
    subroutine prepare_alpha
      real(xk) :: alpha0
      real(xk) :: alpha
      !
      alpha0   = alpha_iteration(0._xk)
      alpha    = alpha_iteration(alpha0)
      gau_alpha = alpha0
      fixed_point: do while(abs(gau_alpha-alpha)>=100._xk*spacing(alpha))
        gau_alpha = alpha
        alpha    = alpha_iteration(alpha)
      end do fixed_point
      !
      write (out,"(/'Gaussian vector-potential envelope exponent = ',g19.11)") gau_alpha
      write (out,"( '    The long-pulse approximate GVP exponent = ',g19.11/)") alpha0
    end subroutine prepare_alpha
    !
    real(xk) function alpha_iteration(alpha)
      real(xk), intent(in) :: alpha
      !
      alpha_iteration = (2._xk/width**2) * (log(2._xk) + log(1._xk + (width*alpha/omega)**2))
    end function alpha_iteration
  end subroutine init_GaussianVP
  !
  !  Sin2 vector-potential envelope implementation
  !
  function Sin2VP(time,omega,phase,origin,width)
    real(xk), intent(in) :: time        ! Time
    real(xk)             :: Sin2VP      ! Vector-potential
    real(xk), intent(in) :: omega, phase, origin, width
    !
    real(xk)             :: t_carrier   ! Time relative to the pulse origin
    !
    t_carrier  = time - origin
    if (abs(t_carrier)>=0.5_xk*width) then
      Sin2VP = 0._xk
    else
      ! The envelope becomes cos^2 relative to the pulse centre
      Sin2VP = (cos(pi_xk*t_carrier/width)**2) * sin(omega*t_carrier+phase)
    end if
  end function Sin2VP
  !
  !
  subroutine init_Sin2VP(tag,vp_param,omega,phase,origin,width)
    character(len=*), intent(in) :: tag
    real(xk), intent(in)         :: vp_param(:)
    real(xk), intent(out)        :: omega, phase, origin, width
    !
    !  Unpack parameters from vp_param()
    !
    omega     = vp_param(1)
    phase     = vp_param(2)
    origin    = vp_param(3)
    width     = vp_param(4)
    !
    write (out,"(/'Using Sin^2 vector-potential envelope ',a/)") trim(tag)
    write (out,"(t5,'         Carrier frequency = ',g24.13,' [(a.u. time)^-1]')") omega
    write (out,"(t5,'   Carrier phase at origin = ',g24.13,' [Radian]')") phase
    write (out,"(t5,'              Pulse origin = ',g24.13,' [a.u. time]')") origin
    write (out,"(t5,'        Full width at zero = ',g24.13,' [a.u. time]')") width
    !
  end subroutine init_Sin2VP
  !
  !  Cos2 E-potential
  !
  function SFAVP(time,omega,phase,origin,width)
    real(xk), intent(in) :: time           
    real(xk)             :: SFAVP    
    real(xk), intent(in) :: omega, phase, origin, width
    !
    real(xk)             :: t_carrier
    !
    t_carrier = time - origin*2*pi_xk/omega
    if (abs(t_carrier)>=(0.5_xk*width*2._xk*pi_xk/omega)) then
      SFAVP = 0._xk
    else
      SFAVP = 0.5_xk*sin((omega/width+omega)*t_carrier+phase)/(omega/width+omega)&
            + 0.5_xk*sin((omega/width-omega)*t_carrier-phase)/(omega/width-omega)&
            + sin(omega*t_carrier+phase)/omega
    end if
  end function SFAVP 
    !
    !
  !
  !  Flat-top vector-potential with Sin^2 raising/falling edge implementation
  !
  function FlatSin2VP(time,omega,phase,origin,width,flat_raise)
    real(xk), intent(in) :: time        ! Time
    real(xk)             :: FlatSin2VP  ! Vector-potential
    real(xk), intent(in) :: omega, phase, origin, width, flat_raise
    !
    real(xk)             :: t_carrier   ! Time relative to the pulse origin
    real(xk)             :: t_edge      ! Edge relative to the pulse origin
    !
    t_carrier  = time - origin
    if (abs(t_carrier)>=0.5_xk*width) then
      FlatSin2VP = 0._xk
    else 
      t_edge = 0.5_xk*width-flat_raise
      !
      !  Pulse is non-zero
      !
      if (abs(t_carrier)>=t_edge) then
        ! The envelope becomes cos^2 relative to the edge
        FlatSin2VP = cos(pi_xk*(abs(t_carrier)-t_edge)/(2*flat_raise))**2
      else 
        ! Flat-top portion
        FlatSin2VP = 1._xk
      end if
      FlatSin2VP = FlatSin2VP * sin(omega*t_carrier+phase)
    end if
  end function FlatSin2VP
  !
  subroutine init_FlatSin2VP(tag,vp_param,omega,phase,origin,width,flat_raise)
    character(len=*), intent(in) :: tag
    real(xk), intent(in)         :: vp_param(:)
    real(xk), intent(out)        :: omega, phase, origin, width, flat_raise
    !
    !  Unpack parameters from vp_param()
    !
    omega      = vp_param(1)
    phase      = vp_param(2)
    origin     = vp_param(3)
    width      = vp_param(4)
    flat_raise = vp_param(11)
    !
    write (out,"(/'Using flat-top vector-potential witn Sin^2 raising/falling edges',a/)") trim(tag)
    write (out,"(t5,'            Carrier frequency = ',g24.13,' [(a.u. time)^-1]')") omega
    write (out,"(t5,'      Carrier phase at origin = ',g24.13,' [Radian]')") phase
    write (out,"(t5,'                 Pulse origin = ',g24.13,' [a.u. time]')") origin
    write (out,"(t5,'           Full width at zero = ',g24.13,' [a.u. time]')") width
    write (out,"(t5,'Raising/falling edge duration = ',g24.13,' [a.u. time]')") flat_raise
    !
    if (width<=2*flat_raise) stop 'vectorpotential_tools%init_FlatSin2VP - bad pulse parameters'
  end subroutine init_FlatSin2VP
end module vectorpotential_tools
