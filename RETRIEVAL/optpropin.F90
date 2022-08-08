
!> Module to calculate VLIDORT input parameters in each layer
!>
!> author:  Martina M. Friedrich
!> date:    Dec 2013
!> module to calculate the VLIDORT input parameters in each layer from the basic 
!> atomospheric properties
!>
!> contains:
!>
!> (1) calc_cross:          calculates the cross section at a given lambda and T
!> (2) calc_gasabs_coef:    calculates the gas absorption coefficient from C, rho and cross section
!> (4) calc_gasscat_coef:   calculates the gas scattering coefficient from rho and Qray
!> (3) calc_rayscatt_cross: calculates the rayleigh scattering crossec from lambda and Delta
!> (0) calc_ref_index_air:  calculates the refractive index of air from lambda,T,P, Pwater, x(CO2)
!> (-) function dDelta:     converts King factor F to depolarization factor
!> (12) calc_Delta:          calculates the depolarization factor at lambda 
!> 7-9) calc_vlidort_input: calculates finally the input desired for vlidort 

module optpropin

use constants

implicit none
real(kind=16) :: one=1.0, two=2.0, three=3.0, six=6.0, seven=7.0
real(kind=16) :: zero=0.0

 contains

subroutine calc_cross(lambda, cross_sectab,cross_sec)
  !> implements (1)
  !> calculates the cross section at wavelength lamba and temperature T
  !> for the moment, the temperature dependence is not implemented
  !> so the only thing that is done is linear interpolating from the table

  real(kind=16), intent (in), allocatable,dimension(:,:)  :: cross_sectab
  real(kind=16), intent (in)  :: lambda     !> wavelength in nm
  real(kind=16), intent (out) :: cross_sec  !> cross section in cm²
  real(kind=16)  :: odpos, resid
  integer           :: ipos,iposp1

  ipos = minloc(abs(cross_sectab(:,1)-lambda),1)
  odpos  = cross_sectab(ipos,1)
  resid  = lambda - odpos

  if (resid .gt. 0.0) then
    iposp1 = ipos +1
  else
    iposp1 = ipos
    ipos   = ipos - 1
    resid  = (cross_sectab(iposp1,1)-cross_sectab(ipos,1))+resid
  endif
  resid = resid / (cross_sectab(iposp1,1)-cross_sectab(ipos,1))
  cross_sec = cross_sectab(ipos,2)+(cross_sectab(iposp1,2)-cross_sectab(ipos,2))*resid
  if (lambda.gt. maxval(cross_sectab(:,1))) then
     write(*,*) 'warning: wavelength is outside range of the cross section file:'
     write(*,*) 'supplied wavelength:   ', lambda
     write(*,*) 'max in cross sec file: ', maxval(cross_sectab(:,1))
  elseif (lambda.lt. minval(cross_sectab(:,1))) then
     write(*,*) 'warning: wavelength is outside range of the cross section file:'
     write(*,*) 'supplied wavelength:   ', lambda
     write(*,*) 'min in cross sec file: ', minval(cross_sectab(:,1))
  endif
end subroutine calc_cross


subroutine calc_gasabs_coef(rho,Cq,crosssec,abs_coef)
  !> implements (2)
  real(kind=16), intent (in)  :: rho        !> gas density in layer in /m³
  real(kind=16), intent (in)  :: Cq         !> volume mixing ratio of gas in ppb
  real(kind=16), intent (in)  :: crosssec   !> cross section in cm²
  real(kind=16), intent (out) :: abs_coef   !> gas absorption coefficient /cm
  real(kind=16)               :: rhoincmp3  !> rho in /cm³
#ifdef TG
  rhoincmp3 = rho/dble(1000000.0)
  abs_coef = crosssec * rhoincmp3 * Cq/dble(1000000000.0)
#endif
#ifdef AEROSOL
  rhoincmp3 = rho/dble(1000000.0)*rho/dble(1000000.0) ! careful, this is now rho²
  abs_coef = crosssec * rhoincmp3 * (Cq/dble(1000000000.0)*Cq/dble(1000000000.0))
#endif
end subroutine calc_gasabs_coef



subroutine calc_gasscat_coef(rho,QRay,scatt_coef)
  !> implements (4)
  real(kind=16), intent (in)  :: rho         !> gas density in layer in /m³
  real(kind=16), intent (in)  :: QRay        !> Rayleigh cross section in cm²
  real(kind=16), intent (out) :: scatt_coef  !> gas scattering coefficient
  real(kind=16)               :: rhoincmp3   !> rho in /cm³
  rhoincmp3 = rho/dble(1000000.0)
  scatt_coef = Qray * rhoincmp3
end subroutine calc_gasscat_coef



subroutine calc_ref_index_air(lambda,Tt,Pp,hh,xCO2,nn,rho)
  !> calculates the air refractive index at a given temperature, pressure, humidity
  !> CO2 content. It is needed to calculate the Rayleigh scattering cross section
  !> I implement the calcuatlions of Ciddor, P.E (1996)
  !> Tested against the values given in the paper: ok
  !> I also use it to calculate the density from the pressure. I return the
  !> number density of molecules of dry air
  !> the densities are in kg/m³
  !> implements (0)
  real(kind=16), intent (in)  :: lambda  !> wavelength in nm
  real(kind=16), intent (in)  :: Tt      !> temperature in K
  real(kind=16), intent (in)  :: Pp      !> pressure in Pa
  !real(kind=16), intent (in)  :: Pw      !> partial pressure water vapor in Pa
  real(kind=16), intent (in)  :: hh      !> fractional humidity (0 --1) (replaces Pw)
  real(kind=16), intent (in)  :: xCO2    !> volume mixing ratio CO2 in ppm
  real(kind=16), intent (out) :: nn      !> refractive index of air
  real(kind=16), intent (out) :: rho     !> number density at given Pp and Tt in 1/ m³
  real(kind=16)  :: Pw                        !> partial pressure water vapor in Pa
  real(kind=16)  :: TtC                       !> temperature in C
  real(kind=16)  :: k0,k1,k2,k3               !> constants needed for Eq. 1,2
  real(kind=16)  :: a0,a1,a2,b0,b1,c0,c1,dd,ee!> constants needed for Eq. 12
  real(kind=16)  :: w0,w1,w2,w3               !> constants needed for Eq. 3
  !real(kind=16)  :: naxs, nas     !> n at 15 deg, 101.325 Pa, 0% hum, CO2= xc/ 450 ppm. 
  !real(kind=16)  :: nws           !> n for water vaport at standart conditions 
  real(kind=16)  :: nasm1         !> (nas-1) x 10⁸
  real(kind=16)  :: naxsm1        !> (naxs -1) x 10⁸
  real(kind=16)  :: nwsm1         !> (nws-1) x 10⁸
  real(kind=16)  :: sigma, sigma2 !> wave number in 1/micro, wave number squared
  real(kind=16)  :: Ma, Za, Zw    !> molar mass of air xCO2, compressibility dry air/water vapor
  real(kind=16)  :: rhoaxs, rhows !> densities of standart air and water vapor
  real(kind=16)  :: rhodry, rhowv !> density of air component and water vapor component
  real(kind=16)  :: nnm1          !> (nn -1)*1e8
  real(kind=16)  :: cf = 1.022d0    !> correction factor for nws
  real(kind=16)  :: xw, Zz        ! xw=ff*Pw/Pp Zz=comperssibility at actual conditions
  real(kind=16)  :: Aa,Bb,Cc,D2,svp ! constants for svp (saturation vapor pressure)
  real(kind=16)  :: alpha, beta, gamma, ff  ! constants for ff (enhancement factor of water vapor in air)
  
  TtC    = Tt -dble(273.15d0)
  k0     = dble(238.0185d0)     ! /micro²
  k1     = dble(5792105.0d0)    ! /micro²
  k2     = dble(57.362d0)       ! /micro²
  k3     = dble(167917.0d0)     ! /micro²
  a0     = dble(1.58123d-6)   ! K/ Pa
  a1     = dble(-2.9331d-8)   ! / Pa
  a2     = dble(1.1043d-10)   ! /(K*Pa)
  b0     = dble(5.707d-6)     ! K/Pa
  b1     = dble(-2.051d-8)    ! /Pa
  c0     = dble(1.9898d-4)    ! K/Pa
  c1     = dble(-2.376d-6)    ! /Pa
  dd     = dble(1.83d-11)     ! K²/Pa²
  ee     = dble(-0.765d-8)    ! K²/Pa²
  w0     = dble(295.235d0)      ! /micro²  ! this actually should not have a unit
  w1     = dble(2.6422d0)       ! /micro²  ! and this should be micro² not 1/micro²
  w2     = dble(-0.032380d0)    ! /micro⁴  ! accordingly, this should be micro⁴
  w3     = dble(0.004028d0)     ! /micro⁶  ! and this micro⁶
  Aa     = dble(1.2378847d-5) ! /K²
  Bb     = dble(-1.9121316d-2)! /K
  Cc     = dble(33.93711047d0)  !
  D2     = dble(-6.3431645d3) ! K
  alpha  = dble(1.00062d0)      !
  beta   = dble(3.14d-8)      ! /Pa
  gamma  = dble(5.6d-7)       ! /C²
  
  svp    = exp(Aa*Tt*Tt+Bb*Tt+Cc+D2/Tt) !(checked)
  Pw     = hh *svp  ! new, changed from direct input.
  ff     = alpha + beta*Pp + gamma*TtC*TtC !(checked)
  xw     = ff* Pw/Pp  ! see page 3 !(checked)
  sigma  = dble(1.0d0)/(lambda/dble(1000d0))
  sigma2 = sigma*sigma
  nasm1  = k1/(k0-sigma2) + k3/(k2-sigma2) !(checked)
  naxsm1 = nasm1*(dble(1.0d0)+dble(0.534d-6)*(xCO2-dble(450.0d0)))!*1.0e-8 ! ??? factor 1e-8?
  nwsm1  = cf *(w0+w1*sigma2+w2*sigma2*sigma2+w3*sigma2*sigma2*sigma2) !checked
  Ma     = 1.0d-3*(dble(28.9635d0)+dble(12.011d-6)*(xCO2-dble(400.0d0)))!checked
  Za     = dble(1.0d0)-dble(101325.0d0)/dble(288.15d0)*(a0+a1*dble(15.0d0)+a2*dble(225.0d0))+ &
               (dble(101325.0d0) /dble(288.15d0)* dble(101325.0d0) /dble(288.15d0))*dd !(checked)
  Zw     = dble(1.0d0)-dble(1333.0d0) /dble(293.15d0) *(a0+a1*dble(20.0d0)+a2*dble(400.0d0)+ &
                b0+b1*dble(20.0d0)+c0+c1*dble(20.0d0))+ &
                (dble(1333.0d0) /dble(293.15d0)* dble(1333.0d0) /dble(293.15d0))*(dd+ee) !(checked)
  rhoaxs =  dble(101325.0d0)*Ma  /(Za*Rgas*dble(288.15d0))   ! in kg/m³ !(checked)
  rhows  =  dble(1333.0d0 )*Mmwv/(Zw*Rgas*dble(293.15d0))    ! in kg/m³ !(checked)
  Zz     = dble(1.0d0)- (Pp/Tt)*   &
          (a0 + a1*TtC + a2*TtC*TtC+ &
          (b0 + b1*TtC) *xw +        &
          (c0 + c1*TtC) *xw*xw) +     &
          (Pp/Tt)*(Pp/Tt)*(dd+ee*xw*xw)  !(checked)
  rhodry = Pp*Ma/(Zz*Rgas*Tt)*(dble(1.0d0)-xw) !(checked)
  rhowv  = Pp*Mmwv*         xw /(Zz*Rgas*Tt)  !(checked)
  rho    = rhodry/Ma*Av   ! number density in /m³
  nnm1   =  (rhodry/rhoaxs) * naxsm1 + (rhowv/rhows) * nwsm1
  !nn     = nnm1*dble(1.0e-8)+dble(1.0)  ! MMF 18.11.2014 I think the factor 1.0e-8 should not be there, therefore I change this to: 
  nn     = nnm1+dble(1.0)
end subroutine calc_ref_index_air


real(kind=16) function dDelta(F)
  real(kind=16) :: F
    dDelta = dble(6.0)* (dble(1.0)-F)/(dble(-7.0)*F-dble(3.0))  ! 9.2 instead of 6?
end function dDelta

subroutine calc_Delta(lambda,Delta)
  !> needed to calculate the rayleigh cross section and the beta coefficients
  !> This subroutine implements the wavelength dependence of the
  !> depolarization factor of N2, O2 given by Bates 1984
  !> maybe later, I can implement that it is also a function of xCO2,  
  !> for the moment, I ignore the change 
  !> implements (12)
  real(kind=16), intent (in) :: lambda  !> wavelength in nm
  real(kind=16), intent (out):: Delta   !> depolarization ratio1
  real(kind=16)              :: lambdainmicro  !> lambda in micro meter
  real(kind=16)              :: lambdainmicro2 !> lambda * lambda
  real(kind=16)              :: xN2  = 780840.0 !> VMR of N2  in ppm
  real(kind=16)              :: xO2  = 209460.0 !> VMR of O2  in ppm
  real(kind=16)              :: xAr  =   9340.0 !> VMR of Ar  in ppm
  real(kind=16)              :: xCO2 =    400.0 !> VMR of CO2 in ppm
  real(kind=16)              :: mil  =1000000.0 !> one million in dp
  real(kind=16)              :: FN2             !> King Factors for N2, O2,Ar, CO2
  real(kind=16)              :: FO2             !> from Bates 1984
  real(kind=16)              :: DeltaN2         !> can be transformed to Delta via
  real(kind=16)              :: DeltaO2         !> F := (6+3*Delta)/(6-7*Delta)
  real(kind=16)              :: DeltaAr = 0.0   !> Ar and CO2 are given not lambda 
  real(kind=16)              :: DeltaCO2= 0.0814! dependend

  lambdainmicro  = lambda/ dble(1000.0) 
  lambdainmicro2 = lambdainmicro*lambdainmicro
  
  FN2     = dble(1.034)+dble(3.17e-4) /lambdainmicro2 
  FO2     = dble(1.096)+dble(1.385e-3)/lambdainmicro2 + &
            dble(1.448e-4)/(lambdainmicro2*lambdainmicro2) 
  DeltaO2 = dDelta(FO2)
  DeltaN2 = dDelta(FN2)
  Delta   = DeltaN2 * xN2 / mil + &
            DeltaO2 * xO2 / mil + &
            DeltaAr * xAr / mil + &
            DeltaCO2* xCO2/ mil 
  Delta = 0.035 ! according to Goody, Yung 1989 just after eq. 7.38
  ! something goes wrong in the calculation above. May 2022
end subroutine calc_Delta

subroutine calc_rayscatt_cross(lambda, m, N2, Delta, rayscatt, beta)
  !> check the influence of xCO2
  !> implements (3)
  !> tested the rayleigh cross section roughly against the approximation
  !> formula in Nicolet,M. 1984
  !> just to get an idea if I do the right thing.

  real(kind=16), intent (in)  :: lambda  !> in nm
  real(kind=16), intent (out) :: rayscatt !> in cm²
  real(kind=16), intent (out) :: beta(0:2)
  real(kind=16), intent (in)  :: m !> refractive index of air, Ciddor, P.E. 1996
  real(kind=16), intent (in)  :: N2 !> number density of air in /m³
  real(kind=16), intent (in)  :: Delta !> = 0.0350
  ! effective depolarization index of air, Goody&Yung 1989, cpt7,p11 instead,
  ! use the wavelength dependence from Bates 1984
  real(kind=16) :: Nincmp3       !> N2 in /cm³
  real(kind=16) :: lambincm      !> lambda in cm
  real(kind=16) :: c1, c2, c3, c4, c5, c6, zzz, yyy

  lambincm = lambda/dble(10000000.0)
  Nincmp3 = N2/dble(1000000.0)

  ! according to Goody, Yung 1989 Eq. 7.37 or Platt et al radiation optics in
  ! atmosphere Eq. 19.5
  ! with correct units from Frohlich and Glenn E. Shaw
  rayscatt = dble(24.0)*pi*pi*pi*(m*m-dble(1.0))*(m*m-dble(1.0)) / &
    (lambincm*lambincm*lambincm*lambincm*Nincmp3*Nincmp3*(m*m+dble(2.0))* &
    (m*m+dble(2.0))) &
            *(dble(6.0)+dble(3.0)*Delta)/(dble(6.0)-dble(7.0)*Delta)
  ! according to Spurr et al 2001
  beta = (/one, zero,(one-Delta)/(two+Delta)/)
  ! There seem to be an issue with rayscatt. Use instead the simplified formula
  ! by Nicolet,M. 1984 (see also radiations_optics_atmos.pdf page 1166
  if (lambda < 550) then
    rayscatt = 4.02d-28/((lambda/1000)**(4+(0.000389*lambda + 94.26/lambda - 0.3328)))
  else
    rayscatt = 4.02d-28/((lambda/1000)**4.04)
  endif
  
  ! implementation as in disort
  c1 = 2.318239D-8
  c2 = 1.798522D-4
  c3 = 1.495366d-5
  c4 = -9.886386d-7
  c5 = 1.078786d-2
  c6 = 3.97839d-38
  zzz = 1.D+6/(lambda*lambda)
  yyy = c1*zzz
  yyy = c2 + zzz*(c3+zzz*(c4+yyy))
  yyy = (1.d0 + zzz*(c5 + zzz*yyy))   
  rayscatt = c6*zzz*zzz*yyy*1.D+10
  !    write(*,*) "Raleigh xs, ref idx ", rayscatt, lambda
end subroutine calc_rayscatt_cross
  
subroutine calc_vlidort_input(ngreek_moments_input,gas_scat, &
                              gas_abs,aer_abs,aer_scat,layerz, gas_beta, &
                              aer_beta,omega, tau,beta,dtaudC, domegadC,dbetadC)
  !> implements (7), (8) and (9)
  !>
  !> If I include aerosols, I should include more than three moments.
  !> also, the moments should start from 0
  !
  integer, intent(in) :: ngreek_moments_input
  real(kind=16), intent (in)  :: gas_scat    !> gas scattering coefficient in /cm
  real(kind=16), intent (in)  :: gas_abs     !> gas absorption coefficient in /cm
  real(kind=16), intent (in)  :: aer_abs     !> aerosol absorption coefficient
  real(kind=16), intent (in)  :: aer_scat    !> aerosol scattering coefficient
  real(kind=16), intent (in)  :: layerz      !> layer thickness in km
  real(kind=16), intent (in)  :: gas_beta(0:2) !> phase function expansion coeff gas
  real(kind=16), intent (in), dimension(:)  :: aer_beta(0:ngreek_moments_input)
  !> phase function expansion coeff aer
  real(kind=16), intent (out) :: omega       !> total single scattering albedo
  real(kind=16), intent (out) :: tau         !> total optical depth 
  real(kind=16), intent (out), dimension(:) :: beta(0:ngreek_moments_input)
  !> total phase function moments
  real(kind=16), intent (out) :: dtaudC      !> derivative C/tau   * d tau   / d C
  real(kind=16), intent (out) :: domegadC    !> derivative C/omega * d omega / d C
  real(kind=16), intent (out),  dimension(:) :: dbetadC(0:ngreek_moments_input)
  !> derivative C/beta  * d beta  / d C
  real(kind=16)               :: extinc      !> total extinction coefficient
  real(kind=16)               :: layerzincm  !> z in cm
  integer                        :: i 


  !write(*,*) 'gas, aerosol scatt', gas_scat, aer_scat

  layerzincm = layerz * dble(100000.0)
  extinc = gas_scat + gas_abs +aer_abs + aer_scat
  tau = layerzincm * extinc
  omega = min((gas_scat +aer_scat) / extinc, 1-1.0d-15)
  ! MMF July 2017: to prevent vlidort input error
  ! Force zeroth moment to be one?
  beta(0) = 1.0

  do i = 1,2
    beta(i) = (gas_beta(i)*gas_scat+aer_beta(i)*aer_scat)/(gas_scat+aer_scat)
  enddo

  do i = 3,ngreek_moments_input
    beta(i) = (aer_beta(i)*aer_scat)/(gas_scat+aer_scat)
  enddo

  ! the following quantities are for the weighting function output with respect
  ! to the concentration of trace gas in the layer, see VLIDORT user guide page 40
#ifdef TG
  dtaudC   =  gas_abs / extinc  ! was tau, but I'm pretty sure should be extinc
  domegadC = -gas_abs / extinc  ! ""
  dbetadC  = 0.0
#endif
#ifdef AEROSOL
  ! In fact, this should be changes induced by changing the absorption coefficient only.
  ! to the changes in AOD in each layer , see VLIDORT user guide page 40
  ! These equations were re-checked (July 2017)
  dtaudC   =  (aer_scat+aer_abs) ! was the same
  !domegadC = ((aer_scat/dtaudC-omega)*dtaudC)/extinc/omega
  ! this is equal to:
  domegadC = (aer_scat - omega * (aer_scat + aer_abs)) / (extinc * omega)

  do i = 0, ngreek_moments_input
   dbetadC(i)  =  aer_scat*(aer_beta(i)-beta(i))/(gas_scat+aer_scat)/beta(i)
   ! checked Aug. 2017
  enddo
  dtaudC = dtaudC/extinc
#endif
end subroutine calc_vlidort_input
end module


