
!> this module contains some physical constants
!>
module constants
 

implicit none 

real(kind=16), parameter :: pi     = 3.1415926359  ! pi
real(kind=16), parameter :: C0inK  = 273.15d0         ! 0 deg C in Kelvin

! the following two parameters are from Ciddor,P.E. (1996)
real(kind=16), parameter :: Rgas = 8.314462618d0     ! gas constant in J/mol/K  
real(kind=16), parameter :: Mmwv = 0.018015d0      ! molar mass of water vapor in kg/mol 
real(kind=16), parameter :: Av   = 6.02214129e23 ! avogadro number #/mol
real(kind=16), parameter :: g0   = 9.80665       ! standart grav acc in m/s²

real(kind=16), parameter :: m2cm = 100.          ! m in cm
real(kind=16), parameter :: km2m = 1000.        ! km in m
real(kind=16), parameter :: ppb  = 1000000000.   ! 1e9




 ! Format specifiers
 CHARACTER(len=70) :: crossecformat =  '(e8.4,5x,es12.6E3))'

contains

subroutine cross_secsetup(infile, crossec_table)

  character(len=*), intent (in)  :: infile
  real(kind=16), dimension(:,:), intent (out), allocatable  :: crossec_table
  integer :: ii, i2, jj=0
  real(kind=16) :: t1, t2
  character(len=200) :: trash
  ! this is solely to determine the length of the table to read in 
  open(1,file=infile,status='old' )
  ii= 0
  do while (1.eq.1)
    read(1,'(A)',END=200) trash
    ii = ii +1
    if (SCAN(trash,";").gt.0) then
        jj = ii
    endif
  enddo
  200 continue
  !close(1)
  ! read in the table
  allocate(crossec_table(1:ii-jj,2))
  rewind(1)
  !write(*,*) "ii, jj", ii, jj
  !open(1,file=infile,status='old' )
  do i2 = 1,jj
    read(1,'(A)') trash
    !write(*,*) trash
  enddo
  
  do i2 = 1, ii-jj-1
    !write(*,*) i2
    read(1,*) t1,t2
    crossec_table(i2,1) = t1
    crossec_table(i2,2) =t2
  enddo
  !write(*,*) "first entry", crossec_table(1,1),crossec_table(i2,1)
  close(1)
  
end subroutine cross_secsetup

end module 
