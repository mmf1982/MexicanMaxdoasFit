
 !> module to initialize the grid and inputs for vlidort that don't change
 !>
 !> author:         Martina M. Friedrich \n
 !> date:           April, 2014 \n
 !> latest change:  Oct, 2019, added FORW option \n

module initialize


 use VLIDORT_pars,        only: zero, vlidort_success
 use VLIDORT_io_defs,     only: VLIDORT_Modified_Inputs, &
                                VLIDORT_Input_Exception_Handling, &
                                VLIDORT_Fixed_Inputs, VLIDORT_Sup_InOut
 use VLIDORT_lin_io_defs, only: VLIDORT_LinSup_InOut, VLIDORT_Fixed_LinInputs, &
                                VLIDORT_Modified_LinInputs
 use VLIDORT_aux

 implicit none

 TYPE(VLIDORT_Fixed_Inputs)             :: VLIDORT_FixIn
 TYPE(VLIDORT_Sup_InOut)                :: VLIDORT_Sup
 TYPE(VLIDORT_LinSup_InOut)             :: VLIDORT_LinSup
 TYPE(VLIDORT_Fixed_LinInputs)          :: VLIDORT_LinFixIn
 TYPE(VLIDORT_Input_Exception_Handling) :: VLIDORT_InputStatus
 TYPE(VLIDORT_Modified_LinInputs)       :: VLIDORT_LinModIn
 TYPE(VLIDORT_Modified_Inputs)          :: VLIDORT_ModIn

 double precision,dimension(:)  , allocatable   :: scatt_coef !> gas scattering coefficient
 double precision,dimension(:,:),allocatable    :: gas_beta   !> gas phase function moments
 double precision,dimension(:)  , allocatable   :: rho        !> total density in 1/m³
 double precision,dimension(:),allocatable :: Ppgrid !> pressure in Pa
 double precision,dimension(:),allocatable :: Ttgrid !> temperature in K
 double precision,dimension(:)  , allocatable   :: abs_aer    !> aerosol absorption coefficient
 double precision,dimension(:)  , allocatable   :: scatt_aer  !> aerosol scattering coefficient
 double precision,dimension(:,:),allocatable    :: aer_beta   !> aerosol phase function moments
 double precision :: cross_sec       !> molecular cross section in cm²

 contains

subroutine initialize_gas(tpfile, height)
  use constants, only: C0inK, cross_secsetup, m2cm, km2m
  use handle_setupfile, only: crosssecfile, lambda, &
                        groundoffset, hh, xCO2, nlayers
  use optpropin, only: calc_cross, calc_ref_index_air, &
                calc_gasabs_coef, calc_rayscatt_cross, calc_gasscat_coef, &
                calc_Delta
  use atmosconditions, only: temppress, temppress2

  implicit none
  character(len=500), intent(in) :: tpfile
  double precision, dimension(:), allocatable, intent(in):: height ! this is above surface
  !double precision,dimension(:),allocatable :: Ppgrid !> pressure in Pa
  !double precision,dimension(:),allocatable :: Ttgrid !> temperature in K
  double precision,dimension(:),allocatable :: nn     !> refractive index
  double precision :: rayscatt                        !> rayleigh scattering cross section
  double precision :: Delta                           !> depolarization ratio of air
  integer :: ii
  double precision, dimension(:,:),allocatable  :: cross_sectable!> to read in

  allocate (Ttgrid(nlayers))
  allocate (Ppgrid(nlayers))
  allocate (rho(nlayers))
  allocate (nn(nlayers))
  allocate (gas_beta(0:2,nlayers))
  allocate (scatt_coef(nlayers))

  call cross_secsetup(crosssecfile, cross_sectable) ! load the molecule cross section table
  call calc_Delta(lambda,Delta)                     ! calculate depolarizaion factor of air
  Ttgrid = 0.0
  ! the height passed to temppress2 should be the height above m.s.l. and hence +groundoffset
  call temppress2(tpfile, groundoffset, height, Ttgrid, Ppgrid)
  Ttgrid = Ttgrid + C0inK         ! Tt in Kelvin

#ifdef FORW
    cross_sec = 0  ! I just need to set something here, in FORW, it is supplied as input
#else
  call calc_cross(lambda, cross_sectable, cross_sec)
! In theory, it should be largely independed of the exact cross section
! below instead of above could be used
!#ifdef AEROSOL
!    cross_sec = 5.0d-46
!#else
!    cross_sec = 5.0d-20
!#endif

  if (cross_sec.le.0.0) then
    write(*,*) "ERROR"
    write(*,*) "The cross section at ", lambda, "nm, is invalid:", cross_sec
    write(*,*) "Check the file: ", crosssecfile
    !cross_sec = 2.0d-46
    ! write(*,*) "Since cross section should not matter, use default of ", cross_sec
    stop
  endif
#endif
  ! write(*,*) "cross sec", cross_sec
  do ii = 1, nlayers
    ! calc refraction index and density
    call calc_ref_index_air(lambda,Ttgrid(ii),Ppgrid(ii),hh(ii),xCO2(ii),nn(ii),rho(ii))
    ! calc Raleigh crosssec & beta
    call calc_rayscatt_cross(lambda, nn(ii), rho(ii), Delta, rayscatt, gas_beta(:,ii))
    ! calc gas scattering coefficient
    call calc_gasscat_coef(rho(ii),rayscatt,scatt_coef(ii))
  enddo

end subroutine initialize_gas

subroutine initialize_aerosol(aerprofshape,taer,waer,gaer,ngreek_moments_input, deltaz)
    use constants, only: m2cm, km2m
    use handle_setupfile, only: nlayers
    implicit none
    double precision, intent(in) :: taer, gaer, waer
    integer, intent(in) :: ngreek_moments_input
    double precision, dimension(:), allocatable, intent(in):: aerprofshape
    double precision, dimension(:), allocatable, intent(in):: deltaz
    integer :: ii
    double precision :: parcel
    allocate (aer_beta(0:ngreek_moments_input,nlayers))
    allocate (abs_aer(nlayers))
    allocate (scatt_aer(nlayers))
    aer_beta(0,:) = 1.0d0
    do ii = 1, ngreek_moments_input
        aer_beta(ii,:) = dble(2*ii+1) * gaer ** dble(ii)
    enddo
    abs_aer = 0.0
    scatt_aer = 0.0
    if (taer /= 0.0) then
        parcel = taer/(sum(deltaz(1:nlayers)*km2m*m2cm*aerprofshape(1:nlayers)))
    else
        parcel = 0.0
    endif
    ! 15.01.2015: The reason why the layer thickness is not included here is
    ! because it is only included when calculating the total tau value in
    ! calc_vlidort_input in optpropin(2).
    do ii = 1, nlayers
        scatt_aer(ii) = waer*parcel*aerprofshape(ii)
        abs_aer(ii) = parcel*(dble(1.0)-waer)*aerprofshape(ii)
    enddo
        ! to check if the total extinction adds up correctly
#ifdef DETAIL
    write(*,*) "aerosol check"
    write(*,*) sum((abs_aer(1:nlayers) + &
                   scatt_aer(1:nlayers))*deltaz(1:nlayers)*km2m*m2cm)
    write(*,*) taer
#endif

end subroutine initialize_aerosol


subroutine initialize_vlidort(vlidortlogfile)
  use handle_setupfile, only: nlayers, height_grid
  character(len=500), intent(in) :: vlidortlogfile
  integer :: n, nf, na

    VLIDORT_FixIn%Cont%TS_NLAYERS_NOMS              = 0
    VLIDORT_FixIn%Cont%TS_NLAYERS_CUTOFF            = 0
    VLIDORT_FixIn%Chapman%TS_PRESSURE_GRID          = ZERO
    VLIDORT_FixIn%Chapman%TS_TEMPERATURE_GRID       = ZERO
    VLIDORT_FixIn%Chapman%TS_FINEGRID               = 0
    !VLIDORT_FixIn%Optical%TS_LTE_DELTAU_VERT_INPUT  = ZERO
    !VLIDORT_FixIn%Optical%TS_LTE_THERMAL_BB_INPUT   = ZERO
    VLIDORT_Sup%BRDF%TS_EXACTDB_BRDFUNC             = ZERO
    VLIDORT_Sup%BRDF%TS_BRDF_F_0                    = ZERO
    VLIDORT_Sup%BRDF%TS_BRDF_F                      = ZERO
    VLIDORT_Sup%BRDF%TS_USER_BRDF_F_0               = ZERO
    VLIDORT_Sup%BRDF%TS_USER_BRDF_F                 = ZERO
    VLIDORT_Sup%BRDF%TS_EMISSIVITY                  = ZERO
    VLIDORT_Sup%BRDF%TS_USER_EMISSIVITY             = ZERO
    VLIDORT_Sup%SLEAVE%TS_SLTERM_ISOTROPIC          = ZERO
    VLIDORT_Sup%SLEAVE%TS_SLTERM_USERANGLES         = ZERO
    VLIDORT_Sup%SLEAVE%TS_SLTERM_F_0                = ZERO
    VLIDORT_Sup%SLEAVE%TS_USER_SLTERM_F_0           = ZERO
    VLIDORT_LinSup%BRDF%TS_LS_EXACTDB_BRDFUNC       = ZERO
    VLIDORT_LinSup%BRDF%TS_LS_BRDF_F_0              = ZERO
    VLIDORT_LinSup%BRDF%TS_LS_BRDF_F                = ZERO
    VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F_0         = ZERO
    VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F           = ZERO
    VLIDORT_LinSup%BRDF%TS_LS_EMISSIVITY            = ZERO
    VLIDORT_LinSup%BRDF%TS_LS_USER_EMISSIVITY       = ZERO
    VLIDORT_ModIn%MBool%TS_DO_FO_CALC = .false. ! 22.06.2015 added. Not sure what this is
    VLIDORT_FixIn%Optical%TS_ATMOS_WAVELENGTH = 0.333d0 ! 22.06.2015 This should not be used
    VLIDORT_LinModIn%MCont%TS_DO_SLEAVE_WFS         = .FALSE.
    VLIDORT_LinFixIn%Cont%TS_N_SLEAVE_WFS           = 0
    VLIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_ISOTROPIC  = ZERO
    VLIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_USERANGLES = ZERO
    VLIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_F_0        = ZERO
    VLIDORT_LinSup%SLEAVE%TS_LSSL_USER_SLTERM_F_0   = ZERO
    VLIDORT_FixIn%Chapman%TS_height_grid(0:NLAYERS)  = height_grid

    VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG(1:NLAYERS)   = .true.
    VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER(1:NLAYERS) = 1

    IF ( VLIDORT_InputStatus%TS_status_inputread .ne. vlidort_success ) then
        open(1,file = vlidortlogfile, status = 'unknown')
        WRITE(1,*)' FATAL:   Wrong input from VLIDORT input file-read'
        WRITE(1,*)'  ------ Here are the messages and actions '
        write(1,'(A,I3)')'    ** Number of messages = ',&
          VLIDORT_InputStatus%TS_NINPUTMESSAGES
        DO N = 1, VLIDORT_InputStatus%TS_NINPUTMESSAGES
          NF = LEN_STRING(VLIDORT_InputStatus%TS_INPUTMESSAGES(N))
          NA = LEN_STRING(VLIDORT_InputStatus%TS_INPUTACTIONS(N))
          write(1,'(A,I3,A,A)')'Message # ',N,' : ',&
            VLIDORT_InputStatus%TS_INPUTMESSAGES(N)(1:NF)
          write(1,'(A,I3,A,A)')'Action  # ',N,' : ',&
            VLIDORT_InputStatus%TS_INPUTACTIONS(N)(1:NA)
        ENDDO
        close(1)
        STOP 'Read-input fail: Look at file 2p6_VLIDORT_ReadInput.log'
    ENDIF

end subroutine initialize_vlidort

end module
