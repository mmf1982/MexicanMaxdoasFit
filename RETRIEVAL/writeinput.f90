
!> Module to read and write DSCDs that resamble QDOAS output
!>
!> author:   Martina M. Friedrich
!> date:     March 2014
!> Module to read and write DSCDs and angle data to mimic those produced by QDOAS
!> it needs a QDOAS blue-print output file.
!> Program assumes a "," as delimiter. 
!> uses module strings from Dr. George Benthien, http://gbenthien.net/strings/
module writeinput

 use strings

 implicit none

 integer :: dataset(0:300)
 double precision, dimension(:,:), allocatable :: dscd, dscderr
 double precision, dimension(:)  , allocatable :: szangle, elevangle, date, &
                                              time, solazim, telazim,rel_azimang
 
 contains
 
 subroutine readwritefile(filename,trace_gas,SCD,usedEA)
  double precision, allocatable, dimension(:)           :: usedEA
 double precision,allocatable,dimension(:),intent(in)   :: SCD
 CHARACTER(len=*), intent(in)                           :: filename 
 CHARACTER(len=650)                                     :: firstline, otherlines
 integer                                                :: io, ii, columns, jj,kk,ll,mm(1:10), rows
 character(len=70),dimension(1:70)                      :: margs
 character(len=100), intent(in)                         :: trace_gas
 character(len=30), dimension(1:10)                     :: gases
 character(len=7), dimension(1:10,1:10)                 :: spectrum
 character(len=50), dimension(:,:), allocatable         :: datamatrix_c
 integer*8, dimension(:),allocatable                    :: specnum,datetime,scannum
 integer                                                :: szangleidx=0,elevangleidx=0,dateidx=0,timeidx=0, &
        specnumidx=0,solazimidx=0,telazimidx=0,datetimeidx=0,dscdidx(10)=0, &
        dscderridx(10)=0, scannumidx=0
 double precision                                       :: test
 double precision                                       :: errorfrac=0.1 !> what fraction of dscd should be written as measurement error
 character(len=200)                                     :: newfilename
 newfilename ="_new"
 newfilename = trim(filename)//newfilename

 open(2,file = trim(newfilename),action='write', status='unknown',iostat=io)
 open(1,file = trim(filename),action='read', status='old',iostat=io)
 if(io.ne.0) then
    write(0,'(A)')'  Error opening DSCD file, aborting...'
    stop
 end if
 read(1,'(A)',iostat=io) firstline
 call parse(firstline,',',margs,columns)

 mm = 0
 kk = 0
 do ii = 1,columns
    select case (margs(ii))
    case ("Date (DD/MM/YYYY)")
        dateidx = ii
    case ("Time (hh:mm:ss)")
        timeidx = ii
    case ("Day number")
        !write(*,*) 'Day Number is', ii
    case ("Scan")
        scannumidx = ii
    case ("Elev. viewing angle")
        elevangleidx = ii
    case ("Azim. viewing angle")
        telazimidx = ii
    case ("Spec No")
        specnumidx = ii
    case ("# Date & time (YYYYMMDDhhmmss)")
        datetimeidx = ii
    case ("Solar Azimuth angle","Solar Azimuth Angle")
        solazimidx = ii
    case("SZA")
        szangleidx = ii
    end select

    if (INDEX(margs(ii), 'Chi', .true.).gt.0) then
       kk = kk+1
       jj = INDEX(margs(ii), 'Chi', .true.)
       gases(kk) = margs(ii)(1:jj-2)
       !write(*,*) gases(kk),'Chi is', ii
    endif
    if (INDEX(margs(ii), 'RMS', .true.).gt.0) then
       !write(*,*) gases(kk), 'RMS is', ii
    endif
    if (INDEX(margs(ii), 'SlCol', .true.).gt.0) then
       ll = INDEX(margs(ii), 'SlCol', .false.)
       mm(kk) =mm(kk)+1
       spectrum(kk,mm(kk)) = margs(ii)(ll+6:len(trim(margs(ii)))-1)
       !write(*,*) gases(kk), 'SlCol ', spectrum(kk,mm(kk)), ' is', ii
       if (INDEX(gases(kk),trim(spectrum(kk,mm(kk))),.true.).gt.0) then
          dscdidx(kk) = ii
       endif
    endif
    if (INDEX(margs(ii),'SlErr', .true.).gt.0) then
       dscderridx(kk) = ii
    endif
 enddo
    write(2,'(A)',advance="no") trim(margs(datetimeidx))
    write(2,'(A)',advance="no") ","
    write(2,'(A)',advance="no") trim(margs(specnumidx))
    write(2,'(A)',advance="no") ","
    write(2,'(A)',advance="no") trim(margs(elevangleidx))
    write(2,'(A)',advance="no") ","
    write(2,'(A)',advance="no") trim(margs(telazimidx))
    write(2,'(A)',advance="no") ","
    write(2,'(A)',advance="no") trim(margs(solazimidx))
    write(2,'(A)',advance="no") ","
    write(2,'(A)',advance="no") trim(margs(szangleidx))
    write(2,'(A)',advance="no") ","
    write(2,'(A)',advance="no") trim(trace_gas) // ".Chi"
    write(2,'(A)',advance="no") ","
    write(2,'(A)',advance="no") trim(trace_gas) // ".RMS"
    write(2,'(A)',advance="no") ","
    write(2,'(A)',advance="no") trim(trace_gas) // ".SlCol"
    write(2,'(A)',advance="no") ","
    write(2,'(A)',advance="yes") trim(trace_gas) // ".SlErr"
!just to find out how many lines of data there are
 rows = 0
 do while (io.eq.0) 
    rows= rows + 1
    read(1,*,iostat=io) otherlines
    if(io.lt.0) exit
    if(io.gt.0) then
       write(0,'(A,I4,A)')'  Error reading ', filename ,', line',ii,' aborting...'
       stop
    end if
 end do
 close(1)
rows=rows-1
! allocate the data vectors
 allocate(datamatrix_c(rows,columns))
 allocate(szangle(rows))
 allocate(elevangle(rows))
 allocate(date(rows))
 allocate(time(rows))
 allocate(specnum(rows))
 allocate(solazim(rows))
 allocate(telazim(rows))
 allocate(datetime(rows))
 allocate(dscd(kk,rows))
 allocate(dscderr(kk,rows))
 allocate(scannum(rows))

! Now, open again to read in the data
 open(1,file = trim(filename),&
      action='read', status='old',iostat=io)
 read(1,'(A)',iostat=io) firstline
 do ii = 1,rows
    read(1,'(A)',iostat=io) otherlines
    call parse(otherlines,',',datamatrix_c(ii,:),columns)
    if (szangleidx.gt.0)   Read(datamatrix_c(ii,szangleidx)  , * ) szangle(ii)
    if (elevangleidx.gt.0) Read(datamatrix_c(ii,elevangleidx), * ) elevangle(ii)
    if (elevangle(ii).gt. 89.9999) elevangle(ii) = dble( 89.9999)
    if (elevangle(ii).le.-89.9999) elevangle(ii) = dble(-89.9999)
    !if (dateidx.gt.0)      Read(datamatrix_c(ii,dateidx)     , * ) date(ii)  ! I don't need that
    !if (timeidx.gt.0)      Read(datamatrix_c(ii,timeidx)     , * ) time(ii)  ! I don't need that
    if (solazimidx.gt.0)   then
       Read(datamatrix_c(ii,solazimidx)  , * ) solazim(ii)
       !solazim(ii) = solazim(ii) +dble(180.0) ! the output from QDOAS seems to be with respect to south
    else
       write(*,*) "no solar azimuth angle found"
       stop
    endif
    if (telazimidx.gt.0) then
       Read(datamatrix_c(ii,telazimidx)  , * ) telazim(ii)
    else 
       write(*,*) "no telescope azimuth angle found"
       stop
    endif
    if (specnumidx.gt.0) then
       Read(datamatrix_c(ii,specnumidx)  , * ) specnum(ii)
    else 
       write(*,*) "no spectrum number found"
       stop
    endif
    if (datetimeidx.gt.0) then
       Read(datamatrix_c(ii,datetimeidx) , * ) datetime(ii)
    else 
       write(*,*) "no no date/ time index"
       stop
    endif
    if (scannumidx.gt.0)   Read(datamatrix_c(ii,scannumidx)  , * ) scannum(ii)
    
    do jj = 1,kk
       Read(datamatrix_c(ii,dscdidx(jj)),* ) dscd(jj,ii)
       Read(datamatrix_c(ii,dscderridx(jj)),*) dscderr(jj,ii)
    enddo
 enddo
 close(1)
 !write(*,*) rows
 do ii = 1, rows !-1
   do jj=1,rows !-1
      if (abs(elevangle(ii)-usedEA(jj)).lt.1.0e-5) then
        !write(*,*) jj
        exit
      endif
   enddo
   if (jj.ge.rows) jj=rows
   write(2,'(i14)',advance="no") datetime(ii)
   write(2,'(A)',advance="no") ","
   write(2,'(i8)',advance="no") specnum(ii)
   write(2,'(A)',advance="no") ","
   write(2,'(f12.4)',advance="no") elevangle(ii)
   write(2,'(A)',advance="no") ","
   write(2,'(f12.4)',advance="no") telazim(ii)
   write(2,'(A)',advance="no") ","
   write(2,'(f12.4)',advance="no") solazim(ii)
   write(2,'(A)',advance="no") ","
   write(2,'(f12.4)',advance="no") szangle(ii)
   write(2,'(A)',advance="no") ","
   write(2,'(A)',advance="no") ' '
   write(2,'(A)',advance="no") ","
   write(2,'(A)',advance="no") ' '
   write(2,'(A)',advance="no") ","
   write(2,'(d12.4)',advance="no") SCD(jj)   ! why is this jj?
   write(2,'(A)',advance="no") ","
   write(2,'(d12.4)',advance="yes") errorfrac*abs(SCD(jj))
 enddo
 
 
! Now, use the scan number to bundle the data into different packages.
! szangle 	holds SZA of each measurement
! solazim 	holds the solar azimuth angle
! elevangle holds the telescope elevation angle
! telazim 	holds the solar azimuth angle
! specnum 	holds the number of the spectrum (the same for each position of the telescope)
! datetime 	holds a unique date/time integer YYYYMMDDhhmmss
! gases		holds the different gases that are analysed
! dscd 		is a matrix, each column vec holds the dscd of a gas, corresponding to gases
! dscderr	is a matrix, each column vec holds the error on dscd of the gas corresponding to gases
! rel_azimang	holds the relative azimuth angle
 
 !I should not do this here
 !allocate(rel_azimang(rows))
 !telazim = (180.0+telazim) ! because telazim is measured against positive telescope elev)
 !rel_azimang = abs(telazim-solazim)
 ! This is to organzie the data, each entry in dataset is the startindex of a new
 ! scan
 !ii = 0
 !jj = 0
 !dataset(0) = 1
 !do while (ii.lt.size(specnum))
 !   jj = jj +1
 !   ii = ii +1
 !   if (ii.lt.size(specnum)) then
 !      do while (specnum(ii+1).gt.specnum(ii))
 !         ii = ii +1
 !         dataset(jj) = ii+1
 !      enddo
 !   endif
 !enddo
  close(2)
end subroutine readwritefile
end module


