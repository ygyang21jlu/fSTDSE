!
!  Sundry mathematical definitions.
!  This is a subset of the routines exported by multigrid's math.f90
!
module math
  use accuracy
  use constants
  !$ use OMP_LIB
  implicit none
  private
  public MathFactorial, MathLogFactorial, MathLegendrePn, MathLegendrePnm
  public MathJacobiPn
  public MathDoubleFactorial, MathLogDoubleFactorial
  public MathBinomial, MathPochhammer
  public Math3J, MathLogSum
  public MathRotationMatrix, MathYJMRotationMatrix
  public MathYLM, MathAllYLM
!
!
  integer(ik), parameter       :: factorial_slack = 5         ! Extra factorials to produce while filling the cache
  integer(ik), save            :: factorial_max = -1          ! Largest value of the factorials cached in the table
  real(rk), allocatable, save  :: factorial_table(:)
  integer(ik), save            :: log_factorial_max = -1      ! Largest value of the factorials cached in the table
  real(rk), allocatable, save  :: log_factorial_table(:)
  integer(ik), save            :: dfactorial_max = -2         ! Largest value of the double factorials cached in the table
  real(rk), allocatable, save  :: dfactorial_table(:)
  integer(ik), save            :: log_dfactorial_max = -2     ! Largest value of the double factorials cached in the table
  real(rk), allocatable, save  :: log_dfactorial_table(:)
!
  contains
  !
  !  External interfaces
  !
  function MathFactorial(n) result(v)
    integer(ik), intent(in)        :: n
    real(rk)                       :: v
    !
    if (n<0) stop 'math%MathFactorialReal - domain error'
    if (n>factorial_max) call fill_factorial_table(n+factorial_slack)
    v = factorial_table(n)
  end function MathFactorial
  !
  function MathLogFactorial(n) result(v)
    integer(ik), intent(in)        :: n
    real(rk)                       :: v
    !
    if (n<0) stop 'math%MathLogFactorial - domain error'
    if (n>log_factorial_max) call fill_log_factorial_table(n+factorial_slack)
    v = log_factorial_table(n)
  end function MathLogFactorial
  !
  function MathDoubleFactorial(n) result(v)
    integer(ik), intent(in)        :: n
    real(rk)                       :: v
    !
    if (n<-1) stop 'math%MathDoubleFactorial - domain error'
    if (n>dfactorial_max) call fill_dfactorial_table(n+factorial_slack)
    v = dfactorial_table(n)
  end function MathDoubleFactorial
  !
  function MathLogDoubleFactorial(n) result(v)
    integer(ik), intent(in) :: n
    real(rk)                :: v
    !
    if (n<-1) stop 'math%MathLogDoubleFactorial - domain error'
    if (n>log_dfactorial_max) call fill_log_dfactorial_table(n+factorial_slack)
    v = log_dfactorial_table(n)
  end function MathLogDoubleFactorial
  !
  !  Pochhammer function: product of (n) integers starting at (a)
  !
  function MathPochhammer(a,n) result(v)
    integer(ik), intent(in) :: a, n
    real(rk)                :: v
    !
    integer(ik) :: aa     ! Starting integer of the equivalent positive sequence
    logical     :: minus
    !
    if (n<0) stop 'math%MathPochhammer - domain error'
    !
    if (n==0) then
      v = 1._rk
      return
    end if
    !
    if (a<=0) then
      !
      !  Catch sequences containing zero factors: these are always zero
      !
      if (a+n-1>=0) then
        v = 0._rk
        return
      end if
      aa    = -(a+n-1)
      minus = mod(n,2)==1
    else
      aa    = a
      minus = .false.
    end if
    !
    v = MathLogFactorial(aa+n-1)-MathLogFactorial(aa-1)
    v = exp(v)
    if (minus) v = -v
  end function MathPochhammer
  !
  !  Binomial coefficients
  !
  function MathBinomial(n,m) result(cnm)
    integer(ik), intent(in)        :: n, m
    real(rk)                       :: cnm
    !
    if (n<0 .or. m<0 .or. m>n) stop 'MathBinomialIntegerReal - domain error'
    cnm = exp(MathLogFactorial(n)-MathLogFactorial(n-m)-MathLogFactorial(m))
  end function MathBinomial
  !
  !  Evaluate log(exp(a)+exp(b)), where a and b may be too large to fit 
  !  in the exponent range.
  !
  function MathLogSum(a,b) result(c)
    real(rk), intent(in) :: a, b
    real(rk)             :: c
    !
    if (a>=b) then
      c = a + log(1._rk + exp(b-a))
    else
      c = b + log(1._rk + exp(a-b))
    end if
  end function MathLogSum
  !
  !  If all values of Pn or Pnm up to a given order are required, it is much
  !  more efficient to compute them all at once!
  !
  !  Accuracy of the recursions used to calculate Ledendre functions
  !  deteriorates in the vicinity of the polynomial roots, especially
  !  for very high orders (n,m)
  !
  function MathLegendrePn(n,x) result(pn)
    integer(ik), intent(in) :: n    ! Order of the Legendre polynomial
    real(rk), intent(in)    :: x    ! Coordinate
    real(rk)                :: pn   ! Value of the Legenre polynomial
    !
    real(rk) :: tab(0:n)
    !
    call legendrePn_table(n,x,tab)
    pn = tab(n)
  end function MathLegendrePn
  !
  function MathLegendrePnm(n,m,x) result(pnm)
    integer(ik), intent(in) :: n, m ! Order of the associated Legendre polynomial
    real(rk), intent(in)    :: x    ! Coordinate, abs(x)<=1
    real(rk)                :: pnm  ! Value of the associated Legenre polynomial
    !
    real(rk) :: tab(0:n,0:m)
    !
    call legendrePnm_table(n,m,x,tab)
    pnm = tab(n,m)
  end function MathLegendrePnm
  !
  !  Use recurrence with respect to degree, see Abramowitz & Stegun 22.7.1
  !
  function MathJacobiPn(n,alp,bet,x) result(pn)
    integer(ik), intent(in) :: n        ! Order of the polynomial
    real(rk), intent(in)    :: alp, bet ! Powers of the weight function, must be > -1
    real(rk), intent(in)    :: x        ! Coordinate, abs(x)<=1
    real(rk)                :: pn
    !
    real(rk) :: tab(0:n)
    !
    call jacobiPn_table(n,alp,bet,x,tab)
    pn = tab(n)
  end function MathJacobiPn
  !
  !  Computes Wigner 3J symbols. The code below is a direct implementation
  !  of L&L 3J formulae. The accuracy of this routine is reduced relative to
  !  that is theoretically possible, due to the use of logarithms. The routine
  !  I had in MNDO99 is more accurate and can handle broader range of J values.
  !
  function Math3J(j1,j2,j3,m1,m2,m3) result(v)
    integer(ik), intent(in) :: j1, j2, j3  ! / J1 J2 J3 \ 3-j
    integer(ik), intent(in) :: m1, m2, m3  ! \ M1 M2 M3 / 
    real(rk)                :: v
    !
    integer(ik) :: ij0, ij1, ij2, ij3, im1a, im1b, im2a, im2b, im3a, im3b
    integer(ik) :: t1, t2, t3
    integer(ik) :: z, minz, maxz
    real(rk)    :: logscale, logterm
    !
    !  Before we do anything, check whether this 3J symbol satisfies the
    !  vector addition constraints
    !
    ij0  =   j1 + j2 + j3 + 1
    ij1  =   j1 + j2 - j3
    ij2  =   j1 - j2 + j3
    ij3  = - j1 + j2 + j3
    im1a =   j1 - m1 ; im1b = j1 + m1
    im2a =   j2 - m2 ; im2b = j2 + m2
    im3a =   j3 - m3 ; im3b = j3 + m3
    if (ij1<0 .or. ij2<0 .or. ij3<0 .or. im1a<0 .or. im1b<0 .or. im2a<0 .or. im2b<0 .or. im3a<0 .or. im3b<0 .or. m1+m2+m3/=0) then
      v = 0
      return
    end if
    !
    logscale = MathLogFactorial(ij1)  + MathLogFactorial(ij2)  + MathLogFactorial(ij3)  &
             + MathLogFactorial(im1a) + MathLogFactorial(im1b) + MathLogFactorial(im2a) &
             + MathLogFactorial(im2b) + MathLogFactorial(im3a) + MathLogFactorial(im3b) &
             - MathLogFactorial(ij0)
    logscale = 0.5_rk * logscale
    !
    t1   = j2 - j3 - m1
    t2   = j1 + m2 - j3
    t3   = j1 - j2 - m3
    minz = max(0_ik,t1,t2)
    maxz = min(ij1,im1a,im2b)
    v = 0
    sum_terms: do z=minz,maxz,1
      logterm = logscale - MathLogFactorial(z)      - MathLogFactorial(ij1-z)  - MathLogFactorial(im1a-z) &
                         - MathLogFactorial(im2b-z) - MathLogFactorial(z-t1)   - MathLogFactorial(z-t2)
      if (abs(logterm)>=0.9_rk*log(huge(1._rk))) then
        write (out,"('Math3J: Intermediate logarithm ',g12.5,' exceeds the real(rk) dynamic range.')") logterm
        write (out,"('Math3J: The 3J arguments were: ',6i10)") j1, j2, j3, m1, m2, m3
        stop 'math%Math3J - exceeded dynamic range'
      end if
      if (mod(z+t3,2)==0) then
        v = v + exp(logterm)
      else
        v = v - exp(logterm)
      end if
    end do sum_terms
    !
  end function Math3J
  !
  !  Given Euler angles, construct rotation matrix for coordinate axes 
  !  in 3D space.  The Euler angles are defined as follows:
  !   1. Rotate coordinate axes by alpha around the Z axis
  !   2. Rotate axes by beta around the new Y axis
  !   3. Rotate axes by gamma around the new Z axis
  !  This definition of the Euler angles matches the definition from
  !  section 58 of L&L III.
  !
  !  Note that prior to June 10th, 2010 the definition of all three Euler 
  !  angles in this routine used a wrong sign, corresponding to anti-screw
  !  rotation sense. Thanks, Mike.
  !
  !  If you would like to rotate an object instead of the axes, take the
  !  transpose or (equivalently) replace (alpha,beta,gamma) with
  !  (-gamma,-beta,-alpha).
  !
  !  Some useful special cases are:
  !
  !    a     b     c  
  !   ---   ---   --- 
  !    x     0     0   Rotate coordinate system by x around the Z axis
  !    0     0     x   Rotate coordinate system by x around the Z axis
  !    0     x     0   Rotate coordinate system by x around the Y axis
  !  -pi/2   x   pi/2  Rotate coordinate system by x around the X axis
  !      
  subroutine MathRotationMatrix(euler_angles,mat)
    real(rk), intent(in)  :: euler_angles(3) ! Euler rotation angles: alpha, beta, and gamma
    real(rk), intent(out) :: mat(3,3)        ! Rotation matrix
    !
    real(rk) :: a, b, g, rma(3,3), rmb(3,3), rmg(3,3)
    !
    a = euler_angles(1)
    b = euler_angles(2)
    g = euler_angles(3)
    !
    rma(1,:) = (/  cos(a),  sin(a),   0._rk /)
    rma(2,:) = (/ -sin(a),  cos(a),   0._rk /)
    rma(3,:) = (/  0._rk,    0._rk,   1._rk /)
    !
    rmb(1,:) = (/  cos(b),   0._rk, -sin(b) /)
    rmb(2,:) = (/  0._rk,    1._rk,   0._rk /)
    rmb(3,:) = (/  sin(b),   0._rk,  cos(b) /)
    !
    rmg(1,:) = (/  cos(g),  sin(g),  0._rk  /)
    rmg(2,:) = (/ -sin(g),  cos(g),  0._rk  /)
    rmg(3,:) = (/  0._rk,    0._rk,  1._rk  /)
    !
    mat = matmul(rmg,matmul(rmb,rma))
  end subroutine MathRotationMatrix
  !
  !  Rotation matrix for angular momentum eigenfunctions, following L&L III Eq. 58.10
  !  Both integer and half-integer J values are OK.
  !
  !  The resulting rotation matrix is accurate to 2ulp for multiplicities up to 6,
  !  with error increasing to 4ulp for multiplicity 20. It loses about 11 decimal places
  !  of accuracy for multiplicity 81, and overflows IEEE double at higher multiplicities.
  !
  !  Note that the rotation matrix uses somewhat weird conventions: it rotates transposed
  !  harmonics from the primed coordinate system defined by the Euler angles back into
  !  the lab system:
  !  
  !    Y(L,M) = Sum Y(L,M') D(M',M)
  !
  !  Furthermore, it looks like the expression for the Wigner matrix in the 5th Russian
  !  edition of L&L is actually incorrect. To get the correct expression, it is necessary
  !  to change the sign of the Euler beta angle. The code below is a literal implementation
  !  of L&L 58.10, so don't forget to flip the sign of beta when calling it!
  !
  subroutine MathYJMRotationMatrix(euler_angles,mult,mat)
    real(rk), intent(in)     :: euler_angles(3) ! Euler rotation angles: alpha, beta, and gamma
                                                 ! See comments in MathRotationMatrix
    integer(ik), intent(in)   :: mult            ! Multipliplicity of the angular-momentum state,
                                                 ! mult = 2*j+1
    complex(rk), intent(out) :: mat(:,:)        ! Rotation matrix
    !
    real(rk)    :: a, b, g
    real(rk)    :: cosb2, sinb2
    complex(rk) :: expa2, expg2
    integer(ik)  :: j2, m2, mp2
    integer(ik)  :: im, imp
    !
    if (mult<1) then
      stop 'math%MathYJMRotationMatrix - multiplicity: domain error'
    end if
    if (size(mat,dim=1)/=mult .or. size(mat,dim=2)/=mult) then
      stop 'math%MathYJMRotationMatrix - rotation matrix dimensions do not match multiplicity'
    end if
    !
    a = euler_angles(1)
    b = euler_angles(2)
    g = euler_angles(3)
    !
    !  We need to take special care when angle beta approaches n*pi. For these angles,
    !  
    !
    sinb2 = sin(0.5_rk*b)
    cosb2 = cos(0.5_rk*b)
    !
    expa2 = exp(cmplx(0._rk,0.5_rk*a,kind=rk))
    expg2 = exp(cmplx(0._rk,0.5_rk*g,kind=rk))
    !
    j2  = mult - 1
    mp2 = -j2
    row_mp: do imp=1,mult
      m2 = -j2
      column_m: do im=1,mult
        mat(imp,im) = sqrt(hf(j2+mp2)*hf(j2-mp2)/(hf(j2+m2)*hf(j2-m2))) &
                    * modJacobiPn((j2-mp2)/2,(mp2-m2)/2,(mp2+m2)/2,sinb2,cosb2) &
                    * expa2**m2 * expg2**mp2
        m2 = m2 + 2
      end do column_m
      mp2 = mp2 + 2
    end do row_mp
    !
    contains
    !
    !  (n/2)! where n must be even
    !
    function hf(n) result(f)
      integer(ik), intent(in) :: n ! Must be even
      real(rk)               :: f
      !
      if (mod(n,2)/=0) stop 'math%MathYJMRotationMatrix%hf - domain error'
      f = MathFactorial(n/2)
    end function hf
    !
    !  Specialized derivative of Jacobi P polynomial:
    !
    !    y^a x^b JacobiP(n,a,b,x^2-y^2)
    !
    !  where x^2+y^2 is equal to 1. Care should be taken in evaluating this function
    !  when either a or b are negative: standard recursion with respect to degree 
    !  becomes undefined in this case.
    !
    !  As the result, we have to use series expansion around +1/-1 argument to
    !  evaluate this function.
    !
    function modJacobiPn(n,a,b,y,x) result(jp)
      integer(ik), intent(in) :: n, a, b ! Parameters of the Jacobi P polynomial
      real(rk), intent(in)   :: y, x    ! Coordinate
      real(rk)               :: jp
      !
      integer(ik) :: k
      !
      if (abs(y)<abs(x)) then
        !
        !  Small y, expand JacobiP around z=+1
        !
        jp = 0
        expand_plus1: do k=max(0,-a),n
          jp = jp + MathPochhammer(a+k+1,n-k) * MathPochhammer(-n,k) * MathPochhammer(a+b+n+1,k) &
                  * y**(2*k+a) / MathFactorial(k)
        end do expand_plus1
        jp = x**b * jp / MathFactorial(n)
      else 
        !
        !  Small x, expand JacobiP around z=-1
        !
        jp = 0
        expand_minus1: do k=max(0,-b),n
          jp = jp + MathPochhammer(b+k+1,n-k) * MathPochhammer(-n,k) * MathPochhammer(a+b+n+1,k) &
                  * x**(2*k+b) / MathFactorial(k)
        end do expand_minus1
        jp = y**a * jp / MathFactorial(n)
        if (mod(n,2)==1) jp = -jp
      endif
!       !
!       !  The general case; do not use
!       !
!       jp = y**a * x**b * MathJacobiPn(n,real(a,kind=rk),real(b,kind=rk),x**2-y**2)
    end function modJacobiPn
  end subroutine MathYJMRotationMatrix
  !
  !  Evaluate a given spherical harmonic. If you need to evaluate the same spherical
  !  harmonic at multiple points, it is preferable to use FLharmonics() in fields.f90
  !
  function MathYLM(l,m,dir) result(ylm)
    integer(ik), intent(in) :: l, m
    real(rk), intent(in)    :: dir(3)   ! (Unnormalized) direction vector
    complex(rk)             :: ylm
    !
    complex(rk) :: harm_fact
    real(rk)    :: r, ct, xymod
    complex(rk) :: xy
    !
    harm_fact = (-1)**((m+abs(m))/2)
    harm_fact = harm_fact * (0._rk,1._rk)**l
    harm_fact = harm_fact * sqrt((2*l+1)/(4*pi))
    harm_fact = harm_fact * sqrt(MathFactorial(l-abs(m))/MathFactorial(l+abs(m)))
    !
    r  = sqrt(sum(dir**2))
    !
    ct = dir(3)/r
    if (m>0) then
      xy = cmplx(dir(1), dir(2),kind=rk)
    else
      xy = cmplx(dir(1),-dir(2),kind=rk)
    end if
    xymod = abs(xy)
    if (xymod>0._rk) then
      xy = xy / xymod
    else
      xy = 1._rk
    end if
    ylm = harm_fact * MathLegendrePnm(l,abs(m),ct) * xy**abs(m)
  end function MathYLM
  !
  !  Calculate all spherical harmonics up to angular momentum l_max
  !
  subroutine MathAllYLM(l_max,dir,ylm)
    integer(ik), intent(in)  :: l_max   ! Desired maximum angular momentum
    real(rk), intent(in)     :: dir(3)  ! (Unnormalized) direction vector
    complex(rk), intent(out) :: ylm(-l_max:l_max,0:l_max)
    !
    integer(ik) :: l, m
    complex(rk) :: harm_fact
    real(rk)    :: r, ct, xymod
    complex(rk) :: xy_plus, xy_minus
    real(rk)    :: legendre_tab(0:l_max,0:l_max)
    complex(rk) :: xy_powm(-l_max:l_max)
    !
    r  = sqrt(sum(dir**2))
    ct = dir(3)/r
    call legendrePnm_table(l_max,l_max,ct,legendre_tab)
    !
    xy_plus  = cmplx(dir(1), dir(2),kind=rk)
    xy_minus = cmplx(dir(1),-dir(2),kind=rk)
    xymod    = abs(xy_plus)
    if (xymod>0._rk) then
      xy_plus  = xy_plus  / xymod
      xy_minus = xy_minus / xymod
    else
      xy_plus  = 1._rk
      xy_minus = 1._rk
    end if
    !
    xy_pow1: do m=-l_max,-1
      xy_powm(m) = xy_minus ** abs(m)
    end do xy_pow1
    xy_powm(0) = 1._rk
    xy_pow2: do m=1,l_max
      xy_powm(m) = xy_plus ** abs(m)
    end do xy_pow2
    !
    tab_l: do l=0,l_max
      tab_m: do m=-l,l
        harm_fact = (-1)**((m+abs(m))/2) * (0._rk,1._rk)**l * sqrt((2*l+1)/(4*pi)) &
                  * sqrt(MathFactorial(l-abs(m))/MathFactorial(l+abs(m)))
        ylm(m,l) = harm_fact * legendre_tab(l,abs(m)) * xy_powm(m)
      end do tab_m
    end do tab_l
  end subroutine MathAllYLM
  !
  !  Auxiliary functions
  !
  subroutine legendrePn_table(nmax,x,pn)
    integer(ik), intent(in) :: nmax  ! Maximum order of the Legendre polynomials desired
    real(rk), intent(in)    :: x     ! Coordinate at which P values are needed
    real(rk), intent(out)   :: pn(:) ! Values of LegendreP from n=0 to n=nmax
    !
    integer(ik) :: n
    real(rk)    :: invn
    !
    if (nmax<0) stop 'math%legendreP_table - negative-order polynomial requested'
    pn(1) = 1._rk
    if (nmax<1) return
    pn(2) = x
    !
    n_recursion: do n=2,nmax
      invn = 1._rk / n
      pn(1+n) = (2._rk-invn)*x*pn(1+(n-1)) - (1._rk-invn)*pn(1+(n-2))
    end do n_recursion
  end subroutine legendrePn_table
  !
  subroutine legendrePnm_table(nmax,mmax,x,pnm)
    integer(ik), intent(in) :: nmax     ! Maximum order of the Legendre polynomials desired
    integer(ik), intent(in) :: mmax     ! Maximum order of the Legendre polynomials desired
    real(rk), intent(in)    :: x        ! Coordinate at which P values are needed, abs(x) must be <=1
    real(rk), intent(out)   :: pnm(:,:) ! Values of LegendreP from n,m=0 to n,m=nmax,mmax
                                        ! n is the first subscript; m is the second subscript
    !
    integer(ik) :: n, m
    real(rk)    :: sqfac ! sqrt(1-x**2)
    !
    sqfac = 1._rk - x**2
    if (sqfac<0) stop 'math%legendrePnm_table - domain error'
    sqfac = sqrt(sqfac)
    !
    call legendrePn_table(nmax,x,pnm(:,1))
    if (mmax<1) return
    !
    !  Special case for m=1: recursion is truncated for n=1, m=1
    !
    pnm(1+0,1+1) = 0._rk
    if (nmax>=1) pnm(1+1,1+1) = -sqfac
    m = 1
    n1_recursion: do n=2,nmax
        pnm(1+n,1+m) = pnm(1+(n-2),1+m) - (2*n-1)*sqfac*pnm(1+(n-1),1+(m-1))
    end do n1_recursion
    !
    m_recursion: do m=2,mmax
      pnm(1+0:1+min(nmax,(m-1)),1+m) = 0._rk
      nm_recursion: do n=m,nmax
        pnm(1+n,1+m) = pnm(1+(n-2),1+m) - (2*n-1)*sqfac*pnm(1+(n-1),1+(m-1))
      end do nm_recursion
    end do m_recursion
  end subroutine legendrePnm_table
  !
  !  Recurrence from A&S 22.7.1. These recurrences are not completely numerically
  !  stable, and lose up to 3 decimal digits for n>10.
  !  Unfortunately, these recurrences do not work for negative integer a and b:
  !  sooner or later, we always hit division by zero. Oops.
  !
  subroutine jacobiPn_table(nmax,a,b,x,pn)
    integer(ik), intent(in) :: nmax  ! Maximum order of the Legendre polynomials desired
    real(rk), intent(in)    :: a, b  ! alpha and beta parameters of the weight function
    real(rk), intent(in)    :: x     ! Coordinate at which P values are needed
    real(rk), intent(out)   :: pn(:) ! Values of LegendreP from n=0 to n=nmax
    !
    integer(ik) :: n
    real(rk)    :: a1n, a2n, a3n, a4n
    !
    if (nmax<0) stop 'math%jacobiPn_table - negative-order polynomial requested'
    pn(1) = 1._rk
    if (nmax<1) return
    pn(2) = 0.5_rk*((a-b) + (2._rk+a+b)*x)
    !
    !  A&S 22.7 recursions are written for a polynomial of degree (n+1).
    !  To avoid confusion, we'll do the same. Note that the a3n coefficient in A&S 22.7.1
    !  is written using non-standard series notation, which is likely to cause confusion ...
    !  Keep in mind that n-th degree polynomial is at pn(n+1)
    !
    n_recursion: do n=1,nmax-1
      a1n = 2 * (n+1)*(n+a+b+1) * (2*n+a+b)
      a2n = (2*n+a+b+1) * (a**2-b**2)
      a3n = (2*n+a+b) * (2*n+a+b+1) * (2*n+a+b+2)
      a4n = 2 * (n+a) * (n+b) * (2*n+a+b+2)
      !
      pn(n+2) = ( (a2n+a3n*x)*pn(n+1) - a4n*pn(n) ) / a1n
    end do n_recursion
  end subroutine jacobiPn_table
  !
  subroutine fill_factorial_table(nmax)
    integer(ik), intent(in) :: nmax
    integer(ik)             :: n, alloc
    real(rk)                :: fac
    !
    !$ if (omp_in_parallel()) then
    !$   stop 'math%fill_factorial_table - unsafe call to MathFactorial'
    !$ end if
    !
    if (factorial_max>=0) then
      deallocate (factorial_table,stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i10,' deallocating factorial table')") alloc
        stop 'math%fill_factorial_table - deallocate'
      end if
    end if
    !
    n   = 0
    fac = 1._rk
    !
    allocate (factorial_table(0:nmax),stat=alloc)
    if (alloc/=0) then
      write (out,"('Error ',i10,' allocating ',i10,'-element factorial table')") & 
             alloc, nmax
      stop 'math%fill_factorial_table - allocate'
    end if
    !
    fill_factorials: do while(n<=nmax-1)
      factorial_table(n) = fac
      n = n + 1
      !
      if (huge(fac)/n<=fac) then
        write (out,"(1x,i10,'! would exceed dynamic range of the chosen real kind')") n
        stop 'math%fill_factorial_table - range exceeded'
      end if
      fac = fac * n
    end do fill_factorials
    factorial_table(n) = fac
    !
    factorial_max = nmax
  end subroutine fill_factorial_table
  !
  subroutine fill_log_factorial_table(nmax)
    integer(ik), intent(in) :: nmax
    integer(ik)             :: n, alloc
    real(rk)                :: fac
    !
    !$ if (omp_in_parallel()) then
    !$   stop 'math%fill_factorial_table - unsafe call to MathLogFactorial'
    !$ end if
    !
    if (log_factorial_max>=0) then
      deallocate (log_factorial_table,stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i10,' deallocating log-factorial table')") alloc
        stop 'math%fill_factorial_table - deallocate'
      end if
    end if
    !
    n   = 0
    fac = 0._rk
    !
    allocate (log_factorial_table(0:nmax),stat=alloc)
    if (alloc/=0) then
      write (out,"('Error ',i10,' allocating ',i10,'-element log-factorial table')") & 
             alloc, nmax
      stop 'math%fill_log_factorial_table - allocate'
    end if
    !
    fill_factorials: do while(n<=nmax-1)
      log_factorial_table(n) = fac
      n = n + 1
      !
      fac = fac + log(real(n,kind=rk))
    end do fill_factorials
    log_factorial_table(n) = fac
    !
    log_factorial_max = nmax
  end subroutine fill_log_factorial_table
  !
  subroutine fill_dfactorial_table(nmax)
    integer(ik), intent(in) :: nmax
    integer(ik)             :: n, alloc
    !
    !$ if (omp_in_parallel()) then
    !$   stop 'math%fill_dfactorial_table - unsafe call to MathDoubleFactorial'
    !$ end if
    !
    if (dfactorial_max>=0) then
      deallocate (dfactorial_table,stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i10,' deallocating double factorial table')") alloc
        stop 'math%fill_dfactorial_table - deallocate'
      end if
    end if
    !
    allocate (dfactorial_table(-1:nmax),stat=alloc)
    if (alloc/=0) then
      write (out,"('Error ',i10,' allocating ',i10,'-element double factorial table')") & 
             alloc, nmax
      stop 'math%fill_dfactorial_table - allocate'
    end if
    !
    dfactorial_table(-1:1) = 1._rk
    n = 2
    fill_factorials: do while(n<=nmax)
      if (huge(1._rk)/n<=dfactorial_table(n-2)) then
        write (out,"(1x,i10,'!! would exceed dynamic range of the chosen real kind')") n
        stop 'math%fill_dfactorial_table - range exceeded'
      end if
      dfactorial_table(n) = dfactorial_table(n-2) * n
      n = n + 1
    end do fill_factorials
    !
    dfactorial_max = nmax
  end subroutine fill_dfactorial_table
  !
  subroutine fill_log_dfactorial_table(nmax)
    integer(ik), intent(in) :: nmax
    integer(ik)             :: n, alloc
    !
    !$ if (omp_in_parallel()) then
    !$   stop 'math%fill_dfactorial_table - unsafe call to MathLogDoubleFactorial'
    !$ end if
    !
    if (log_dfactorial_max>=0) then
      deallocate (log_dfactorial_table,stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i10,' deallocating log-double factorial table')") alloc
        stop 'math%fill_dfactorial_table - deallocate'
      end if
    end if
    !
    allocate (log_dfactorial_table(-1:nmax),stat=alloc)
    if (alloc/=0) then
      write (out,"('Error ',i10,' allocating ',i10,'-element log-double factorial table')") & 
             alloc, nmax
      stop 'math%fill_log_dfactorial_table - allocate'
    end if
    !
    log_dfactorial_table(-1:1) = 0._rk
    n = 2
    fill_factorials: do while(n<=nmax)
      log_dfactorial_table(n) = log_dfactorial_table(n-2) + log(real(n,kind=rk))
      n = n + 1
    end do fill_factorials
    !
    log_dfactorial_max = nmax
  end subroutine fill_log_dfactorial_table
  !
end module math
