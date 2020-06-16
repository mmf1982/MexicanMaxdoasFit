!> This module contains routines that take in a height grid and return temperature and pressure grids.
!>
!> One uses the ISA, the other measured data. \n
!> Author: Martina M. Fredrich

module atmosconditions




use constants 
use strings, only: parse
use sorting, only:  Sort

implicit none

  !double precision, allocatable, dimension(:) :: pressure1, temperature1, height1

  double precision, dimension(8) :: heightlow= dble((/  0.0,11.0,20.0,32.0,&
                                                        47.0,51.0,71.0,84.852/))
  double precision, dimension(8) :: Templow  = dble((/ 15.0,-56.5,-56.5,&
                                                      -44.5,-2.5,-2.5,-58.5,-86.28/))
  double precision, dimension(8) :: Tslope   = dble((/ -6.5,0.0,1.0,2.8,0.0,&
                                                       -2.8,-2.0,0.0/))
  double precision, dimension(8) :: Presslow = dble((/101325.0,22632.0,&
                                                      5474.89,868.02,110.91,66.94,3.96,0.3734/))

  
 contains
  subroutine temppress(height, temperature, pressure)
!> author:   Martina M. Friedrich
!> date:     Dec, 2013
! 
!> This module takes in a height grid and calculates the temperatures and pressures in the other layers 
!> assuming ISA (international standard atmosphere) 
!> http://en.wikipedia.org/wiki/Barometric_formula
!> in:  height (scalar) in km
!> out: pressure (scalar) in Pa, temperature (scalar) in C
!

    double precision, dimension(:), allocatable, intent (out) :: temperature !> in C
    double precision, dimension(:), allocatable, intent (out) :: pressure !> in Pa 
    double precision, dimension(:), allocatable, intent (in)  :: height  !> in km
    double precision, dimension(:), allocatable :: hinm
    double precision               :: odpos, resid
    integer                        :: ipos, ii
    double precision               :: hlowinm, TlowinK,slopinm
    double precision               :: M = 0.0289644 ! molar mass of earth's air 
    !in kg/ mol as function of xCO2, see optpropin, calc_ref_index_air

    allocate(hinm(size(height)))
    allocate(temperature(size(height)))  ! WHY DO I NEED TO ALLOCATE T AND P HERE
    allocate(pressure(size(height)))     ! AS WELL AS IN INITIALIZE.F90, FROM WHERE I CALL?
    hinm = height * dble(1000.0)
    do ii = 1, size(height)
       ipos  = minloc(abs(heightlow - height(ii)),1)
       odpos = heightlow(ipos)
       resid = height(ii) - odpos
       if (resid .lt. 0.0) then 
           ipos   = ipos -1
       endif
       temperature(ii) = Templow(ipos) + Tslope(ipos) * (height(ii) - heightlow(ipos))
       hlowinm = heightlow(ipos) *dble(1000.0)
       TlowinK = Templow(ipos) + dble(273.15)
       slopinm = Tslope(ipos)/ dble(1000.0)
       if (abs(Tslope(ipos)).le.0.1) then
           pressure(ii) = Presslow(ipos)*exp(-g0*M*(hinm(ii)-hlowinm)/(Rgas*TlowinK))
       else 
           pressure(ii) = Presslow(ipos)* &
       (TlowinK/(TlowinK+slopinm*(hinm(ii)-hlowinm)))**(g0*M/(Rgas*slopinm))
       endif
    enddo
  end subroutine temppress

  subroutine readtemppress(stationheight,inputfile,height1,temperature1,pressure1)
  ! this subroutine replaces the one below. It assumes a simple ascii table with
  ! the following structure: the altitude is a.s.l
  ! altitude (km) pressure (hPa) temp (K)
  character (len=*), intent(in) :: inputfile
  character (len=200) :: firstline
  integer :: ii=0, jj=1, io=0, kk=0, columns !jj=2 accounts for 1 line of header ii was 2
  double precision :: testheight
  double precision, intent(in) :: stationheight
  double precision, allocatable, dimension(:), intent(out) :: pressure1  ! in hPa
  double precision, allocatable, dimension(:), intent(out) :: temperature1  ! in C
  double precision, allocatable, dimension(:), intent(out) :: height1  ! in km now
  character (len=70) :: margs(12)
  open(1,file=inputfile,status='old',iostat=io)
  if(io.ne.0) then
    write(0,'(A)')'  Error opening temperature/ pressure file, aborting...'
    write(0,*) inputfile
    stop
  end if
  read(1,'(A)') firstline
  do while (io.eq.0)
        read(1,'(A)',iostat=io) firstline
        !call parse(firstline,' ', margs,columns)
        !read(margs(1),*) testheight
        !if (testheight.lt.stationheight) jj = ii
        ii = ii + 1
  enddo
  !ii = ii - jj
  ii = ii - 1
  allocate(pressure1(ii))
  allocate(temperature1(ii))
  allocate(height1(ii))
  rewind(1)
  !do kk = 1,jj
  read(1,'(A)') firstline
  !enddo
  do kk=1,ii
     read(1,'(A)',iostat=io) firstline
     call parse(firstline,' ', margs,columns)
     read(margs(1),*)  height1(kk)
     read(margs(2),*)  pressure1(kk)
     read(margs(3),*)  temperature1(kk)
     ! write(*,*) height1(kk), temperature1(kk)
  end do
  temperature1 = temperature1 - C0inK  ! to convert to degree Celsius
  end subroutine readtemppress

  subroutine temppress2(tpfile, stationheight, height2, temperature2, pressure2) ! km, C, Pa
  character(len=500), intent(in) :: tpfile
  double precision, allocatable, dimension(:), intent(out) :: pressure2, temperature2
  double precision, allocatable, dimension(:), intent(in) :: height2  ! height above ground
  double precision, intent(in) :: stationheight
  double precision, allocatable, dimension(:) :: pressure1,temperature1, height1  ! this height is a.s.l
  integer :: mydim,ipos,iposp1, ii,mysize,ip
  double precision :: odpos,resid
  double precision, allocatable, dimension(:) :: h2_1,h2_2,p2_1,p2_2,t2_1,t2_2
  call readtemppress(stationheight, tpfile, height1,temperature1,pressure1) ! km, C, hPa
  pressure1 = pressure1*dble(100.0) ! Need to convert P in hPA to PA
  height1 = height1 - stationheight ! Need to remove groudoffset
  mydim = size(height2)
  allocate(temperature2(mydim))
  allocate(pressure2(mydim))
  ! MMF July 2018: I think I need to sort T and P according to H first: No
  ! call Sort(height1, pressure1, temperature1, size(height1)) not needed
  
  
  
  ! disable the possibility to have a too small TP height grid.
  !if (maxval(height2).gt.maxval(height1)) then
  !   ip=minloc(abs(height2(:)-maxval(height1)),1)
  !   if (height2(ip)-maxval(height1).lt.0.0) ip=ip-1
  !   mysize=mydim-ip
  !   allocate(h2_1(mysize))
  !   allocate(t2_1(mysize))
  !   allocate(p2_1(mysize))
  !   allocate(h2_2(ip))
  !   allocate(t2_2(ip))
  !   allocate(p2_2(ip))
  !   h2_1=height2(ip+1:)
  !   h2_2=height2(1:ip)
  !   !do ii=1,mysize
  !   call temppress(h2_2,t2_2,p2_2)
  !   !enddo
  !   temperature2(1:ip) = t2_2
  !   pressure2(1:ip) = p2_2
  !else
     allocate(h2_1(mydim))
     allocate(t2_1(mydim))
     allocate(p2_1(mydim))
     h2_1 = height2
     mysize=mydim
     ip = 0
  !endif
  ! Here comes the interpolation for the temperature2 and pressure2 on the height2 grid
  do ii = 1,mysize
      ipos = minloc(abs(height1(:)-h2_1(ii)),1)
      odpos  = height1(ipos)
      resid  = h2_1(ii) - odpos
      if (resid .gt. 0.0) then
        iposp1 = ipos +1
        resid = resid / (height1(iposp1)-height1(ipos))
      elseif (resid .lt. 0.0) then
        iposp1 = ipos
        ipos   = ipos - 1
        resid  = (height1(iposp1)-height1(ipos))+resid
        resid = resid / (height1(iposp1)-height1(ipos))
      else
        iposp1 = ipos
        resid = 0
      endif
      p2_1(ii) = pressure1(ipos)+(pressure1(iposp1)-pressure1(ipos))*resid
      t2_1(ii) = temperature1(ipos)+(temperature1(iposp1)-temperature1(ipos))*resid
  enddo
  temperature2(ip+1:) = t2_1
  pressure2(ip+1:) = p2_1
  end subroutine temppress2
end module


