!
!  Type and precision definitions
!  Theoretically, changing the value of rk and/or ik here should be sufficient to 
!  switch the entire program to the desired accuracy.
!
module accuracy
  implicit none
  private
  public ik, rk, xk
  public input, out
  public clen
  public srk, drk, qrk
  public rk_bytes, ik_bytes
  public flush_wrapper
  !
  integer, parameter :: ik     = selected_int_kind(8)       ! "Normal" integers
! integer, parameter :: ik     = selected_int_kind(15)      ! "Long" integers
! integer, parameter :: rk     = selected_real_kind(6,17)   ! Low-precision
  integer, parameter :: rk     = selected_real_kind(14,17)  ! "Normal" reals and complex (complexi? :-)
! integer, parameter :: rk     = selected_real_kind(28,50)  ! quad-precision
  !
  !  ifort 15.0.1 mistakenly issues a warning for the following line, complaining about 
  !  the mismatch in the type/kind of precision(1._rk) and 14. Disregard the warning.
  integer, parameter :: xk     = selected_real_kind(max(precision(1._rk),14),max(range(1._rk),17))
                                                            ! Real kind used for manipulating time intervals
                                                            ! and field. Will have at least the range and
                                                            ! accuracy of the (rk) real kind, but not less
                                                            ! than the range and accuracy of double-precision type
                                                            ! This real kind should never appear on time-critical
                                                            ! path.
  integer, parameter :: input  = 5                          ! Standard input I/O channel
  integer, parameter :: out    = 6                          ! Output I/O channel
  integer, parameter :: clen   = 255                        ! Standard character length; enough for most
                                                            ! keywords and file names
  !
  !  System kinds; should only be used where we interface to external libraries
  !
  integer, parameter :: srk    = kind(1.)                   ! System default real kind
  integer, parameter :: drk    = kind(1.d0)                 ! System double-precision real kind
!*nq  integer, parameter :: qrk    = kind(1.d0)             ! No quad-precision real kind; use double dummy
!*qd  integer, parameter :: qrk    = kind(1.q0)             ! System quad-precision real kind
  !
  contains

  integer(ik) function rk_bytes ()
    rk_bytes = real_bytes(radix(1._rk) ,digits(1._rk ),maxexponent(1._rk )-minexponent(1._rk ))
  end function rk_bytes
 
  integer(ik) function ik_bytes ()
    ik_bytes = int_bytes(radix(1_ik ),digits(1_ik ))
  end function ik_bytes

  integer(ik) function int_bytes(radix,digits)
    integer(ik), intent(in) :: radix, digits
    !
    real(rk) :: bits
    !
    bits = digits*(log(real(radix,kind=rk))/log(2._rk))
    int_bytes = ceiling((1+bits)/8._rk)
  end function int_bytes

  integer(ik) function real_bytes(radix,digits,exp_range)
    integer(ik), intent(in) :: radix, digits, exp_range
    !
    real(rk) :: exp_bits, mant_bits
    !
    exp_bits   = log(real(exp_range,kind=rk))/log(2._rk)
    mant_bits  = digits*(log(real(radix,kind=rk))/log(2._rk))
    real_bytes = ceiling((exp_bits+mant_bits)/8._rk)
  end function real_bytes

  subroutine flush_wrapper(unit)
    integer, intent(in) :: unit
    !
    flush (unit)  ! Fortran-2003 form
    ! intrinsic flush
    ! call flush(unit) ! Common extension; uncomment if Fortran-2003 statement is not available
  end subroutine flush_wrapper

end module accuracy
