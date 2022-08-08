
!> This module contains methods to read the setup file
!>
!> author: Martina M. Friedrich
!> date:  March, 2014
!
!> module to read in the setup file for each file that should be analyzed.
!> the setup file should contain information about the number of layers, 
!> the grid boundaries, directory and file names, apriori profile, etc
module handle_setupfile

implicit none

  character (len=500) :: input_filename         ! input filename
  character (len=500) :: apriori_file
  character (len=500) :: height_file
  character (len=500) :: tpfile
  character (len=500) :: Sa_file
  character (len=500) :: vlidortlogfile
  character (len=500) :: output_write_dir
  character (len=500) :: crosssecfile
  character(len=500) :: vlidort_config
  logical             :: do_read_Sa
  real(kind=16)     :: lambda                 ! wavelength in nm
  integer             :: nlayers                ! number of layers for simulation
  integer             :: nlayers_retrieve       ! number of layers for retrieval
  character (len=10)  :: S1
  real(kind=16)    :: groundoffset           ! altitude in km 
  real(kind=16)    :: l_albedo               ! surface albedo
  real(kind=16)    :: Saparam1               ! parameters for Sa matrix, scaling!
  real(kind=16)    :: scaleheight            ! scaleheight for off-diag elements of Sa in km 
  real(kind=16), dimension(:), allocatable :: height_grid ! height grid in km 
  real(kind=16), dimension(:), allocatable :: hh          ! fractional humidity
  real(kind=16), dimension(:), allocatable :: xCO2        ! VMR of CO2 in ppm
  real(kind=16), dimension(:), allocatable :: aerprofshape ! aerosol profile shape
  real(kind=16) :: taer, gaer, waer, real0 ! aerosol ext, unisotropy, omega, zenith
  character (len=200) ::  SaM1file

contains

subroutine readsetupfile(filename)

use VLIDORT_AUX

implicit none
  character (len=*), intent(in) :: filename
  logical             :: isthere
  character (len=160) :: PAR_STR
  integer             :: ii, io
  character (len=1)   :: prefix ='*'
  logical             :: error
  integer             :: FILUNIT = 12

  open(FILUNIT,file =filename,action='read', status='old',iostat=io)
  if(io.ne.0) then
    write(0,'(A)')'  Error opening setupfile, aborting...', filename
    stop
  end if

  PAR_STR = 'number of layers'
  IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) then
  ! the lowest most NLAYERS_RETRIEVE will be used for the retrieval
        READ (FILUNIT,*) NLAYERS
        READ (FILUNIT, *) NLAYERS_RETRIEVE
        allocate(HEIGHT_GRID(0:NLAYERS))
        allocate(hh(NLAYERS))
        allocate(xCO2(NLAYERS))
        allocate(aerprofshape(nlayers))
  ENDIF


  do_read_Sa = .false.
  PAR_STR = "Sa read in?"
  IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) then
        READ (FILUNIT,*) S1
        S1 = to_upper(S1)
        if (S1.eq.'.TRUE.'.or.S1.eq.'.T.'.or.S1.eq.'TRUE'.or.S1.eq."T") then
            do_read_Sa=.true.
        else
            do_read_Sa=.false.
  ENDIF

  PAR_STR = 'upper layer boundaries above ground in km'
  IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) then
        do ii = 0, NLAYERS-1
            READ (FILUNIT,*) HEIGHT_GRID(ii)
        enddo
        HEIGHT_GRID(NLAYERS) = 0.0
  ENDIF

  PAR_STR = 'input filename'
  IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
        READ (FILUNIT,*) INPUT_FILENAME

  PAR_STR = 'output directory'
  IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
        READ (FILUNIT,*) output_write_dir

  PAR_STR = 'apriori profile filename'
  IF (GFINDPAR (FILUNIT, PREFIX, ERROR, PAR_STR)) then
        READ (FILUNIT,*)apriori_file
  ENDIF

  PAR_STR = 'Sa filename'
  IF (GFINDPAR (FILUNIT, PREFIX, ERROR, PAR_STR)) then
        READ (FILUNIT,*) Sa_file
  ENDIF

  PAR_STR = 'heights for apriori profile filename'
  IF (GFINDPAR (FILUNIT, PREFIX, ERROR, PAR_STR)) then
        READ (FILUNIT,*) height_file
  ENDIF

  l_albedo = 0.06
  PAR_STR = 'surface albedo'
  IF (GFINDPAR (FILUNIT, PREFIX, ERROR, PAR_STR)) then
        READ (FILUNIT,*) l_albedo
  ENDIF

  PAR_STR = 'vlidort config filename'
  IF (GFINDPAR (FILUNIT, PREFIX, ERROR, PAR_STR)) then
        READ (FILUNIT,*) vlidort_config
  ENDIF

  PAR_STR = 'temperature pressure filename'
  IF (GFINDPAR (FILUNIT, PREFIX, ERROR, PAR_STR)) then
        READ (FILUNIT,*) tpfile
  ENDIF

  hh = 0.0
  PAR_STR = 'fractional humidity'
  IF (GFINDPAR (FILUNIT, PREFIX, ERROR, PAR_STR)) then
         do ii = 1, NLAYERS
           READ (FILUNIT,*) hh(ii)
         enddo
  ENDIF

  xCO2 = 400.0
  PAR_STR = 'VMR of CO2 in ppm'
  IF (GFINDPAR (FILUNIT, PREFIX, ERROR, PAR_STR)) then
         do ii = 1, NLAYERS
           READ (FILUNIT,*) xCO2(ii)
         enddo
  ENDIF

  aerprofshape = 0.00001
  PAR_STR = 'aerosol profile shape'
  !write(*,*) NLAYERS_RETRIEVE+1, NLAYERS
  !write(*,*) aerprofshape
  IF (GFINDPAR (FILUNIT, PREFIX, ERROR, PAR_STR)) then
         do ii = NLAYERS-NLAYERS_RETRIEVE+1, NLAYERS
           READ (FILUNIT,*) aerprofshape(ii)
         enddo
  ENDIF

  PAR_STR = 'wavelength in nm'
  IF (GFINDPAR (FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*) lambda

  scaleheight = 0.2
  PAR_STR = 'Sa correlation length in km'
  IF (GFINDPAR (FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*) scaleheight

  groundoffset = 0.0
  PAR_STR = 'altitude in km'
  IF (GFINDPAR (FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*) groundoffset

  PAR_STR = 'cross section folder'
  IF (GFINDPAR (FILUNIT, PREFIX, ERROR, PAR_STR)) then
        READ (FILUNIT,*) crosssecfile
        inquire(FILUNIT , EXIST=isthere)
        if (.not.isthere) then
            write(0,'(A)')'  Error, cross section dir not existent'
            write(0,'(A)') 'tried opening ', crosssecfile
            stop
        endif
  endif

  Saparam1 = 1.0
  PAR_STR = 'scaling parameter for Sa matrix'
  IF (GFINDPAR (FILUNIT, PREFIX, ERROR, PAR_STR)) then
           READ (FILUNIT,*) Saparam1
  endif

  PAR_STR = 'aerosol ext t, unisotropy g, omega w'
  IF (GFINDPAR (FILUNIT, PREFIX, ERROR, PAR_STR)) then
           READ (FILUNIT,*) taer
           READ (FILUNIT,*) gaer
           READ (FILUNIT,*) waer
  ENDIF

  real0 = 0.01
  PAR_STR = 'real elevation angle for zenith direction'
  IF (GFINDPAR (FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*) real0
  endif

  close(FILUNIT)
  !write(*,*) "aerosol:", aerprofshape
  !write(*,*) "HEIGHT_GRID", HEIGHT_GRID

end subroutine readsetupfile

Pure Function to_upper (str) Result (string)

!   ==============================
!   Changes a string to upper case
!   ==============================

    Implicit None
    Character(*), Intent(In) :: str
    Character(LEN(str))      :: string

    Integer :: ic, i

    Character(26), Parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    Character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

!   Capitalize each letter if it is lowecase
    string = str
    do i = 1, LEN_TRIM(str)
        ic = INDEX(low, str(i:i))
        if (ic > 0) string(i:i) = cap(ic:ic)
    end do

End Function to_upper

end module


