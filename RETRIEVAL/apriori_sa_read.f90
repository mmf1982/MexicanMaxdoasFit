
!> This module contains routines to handle data to read a-priori and Sa matrix
!>
!> author:  Martina M. Friedrich
!> date:    February 2018
!> 
!> 
!> 

module apriori_sa_read
  use sorting, only:  Sort, Sortmat
  use constants, only: km2m
  
  implicit none

  contains
  ! need to implement 2D interpolation for Sa if needed
  subroutine read_apriori(apriori_file, height_file, profile, hgrid)
    character(len=*), intent(in) :: apriori_file
    character(len=*), intent(in) :: height_file
    real(kind=16), dimension(:), allocatable, intent(inout) :: profile
    real(kind=16), dimension(:), allocatable, intent(in) :: hgrid  ! in km
    real(kind=16), dimension(:), allocatable :: hgrid_external  ! in m
    real(kind=16), dimension(:), allocatable :: profile_external
    integer :: rows, io, rows2, sizehgrid, ipos, iposp1, ii
    real(kind=16) :: hipos, resid, pipos, piposp1, hiposp1

    sizehgrid=size(hgrid)
    open(1,file= trim(apriori_file), action='read',status='old',iostat=io)
    if(io.ne.0) then
        write(0,'(A)')'  Error opening apriori file: ', apriori_file
        stop
    endif
    read(1,*, iostat=io) rows
    if(io.ne.0) then
        write(0,'(A)')'Error reading number of layers in: ', apriori_file
        stop
    endif
    open(2,file= height_file, action='read',status='old',iostat=io)
    if(io.ne.0) then
        write(0,'(A)')'  Error opening height file: ', height_file
        stop
    endif
    read(2,*, iostat=io) rows2
    if(io.ne.0) then
        write(0,'(A)')'Error reading number of layers in: ', height_file
        stop
    end if
    if (rows.ne.rows2) then
        write(0,'(A)') 'Error apriori and height files unequal rows: '
        write(0,'(A)') apriori_file, height_file, rows, rows2
    endif

    allocate(hgrid_external(rows))
    allocate(profile_external(rows))
    read(1,*, iostat=io) profile_external
    if(io.ne.0) then
        write(0,'(A)')'Error reading apriori file:', apriori_file
        stop
    endif
    read(2,*, iostat=io) hgrid_external
    if(io.ne.0) then
        write(0,'(A)')'Error reading height file:', height_file
        stop
    endif
    close(1)
    close(2)
    
    !write(*,*) "external profile"
    !do ii = 1, rows
    !    write(*,*) hgrid_external(ii), profile_external(ii)
    !enddo
    ! write(*,*) "hgrid, prof"
    hgrid_external = hgrid_external/km2m
    ! asure that the externally provided grid is ordered from ground to top:
    call Sort(hgrid_external, profile_external, profile_external, rows)
    do ii = 1, sizehgrid
        if (hgrid_external(rows) <= hgrid(ii)) then
            profile(ii) = profile_external(rows)
        else
            ipos  = minloc(abs(hgrid_external-hgrid(ii)),1)
            hipos = hgrid_external(ipos)
            resid  = (hgrid(ii) - hipos)
            if (resid.lt.0.0) then
                iposp1 = ipos -1
            else
                iposp1 = ipos +1
            endif
            hiposp1 = hgrid_external(iposp1)
            pipos = profile_external(ipos)
            piposp1 = profile_external(iposp1)
            profile(ii) = pipos + (piposp1-pipos) * resid/(hiposp1-hipos)
       endif
       !write(*,*) hgrid(ii), profile(ii) !, ii, ipos, iposp1, minval(abs(hgrid_external-hgrid(ii)))
    enddo
  end subroutine read_apriori

  subroutine read_sa(Sa_file, height_file, Sa, hgrid)
    character(len=500), intent(in) :: Sa_file
    character(len=500), intent(in) :: height_file
    real(kind=16), dimension(:,:), allocatable, intent(inout) :: Sa
    real(kind=16), dimension(:), intent(in) :: hgrid  ! in km
    real(kind=16), dimension(:), allocatable :: hgrid_external  ! in m
    real(kind=16), dimension(:,:), allocatable :: Sa_external
    integer :: rows, io, rows2, sizehgrid, ipos, iposp1, ipos2, iposp1_2, ii, jj
    real(kind=16) :: hipos, resid, pipos, piposp1, hipos2, resid2, temp1, temp2
    real(kind=16) :: pipos11, pipos12, pipos21, pipos22, hiposp1_2, hiposp1
    sizehgrid=size(hgrid)

    open(1,file= Sa_file, action='read',status='old',iostat=io)
    if(io.ne.0) then
        write(0,'(A)')'  Error opening Sa file: ', Sa_file
        stop
    endif
    read(1,*, iostat=io) rows
    if(io.ne.0) then
        write(0,'(A)')'Error reading number of layers in: ', Sa_file
        stop
    endif
    open(2,file= height_file, action='read',status='old',iostat=io)
    if(io.ne.0) then
        write(0,'(A)')'  Error opening height file: ', height_file
        stop
    endif
    read(2,*, iostat=io) rows2
    if(io.ne.0) then
        write(0,'(A)')'Error reading number of layers in: ', height_file
        stop
    end if
    if (rows.ne.rows2) then
        write(0,'(A)') 'Error Sa and height files unequal rows: '
        write(0,'(A)') Sa_file, height_file, rows, rows2
    endif

    allocate(hgrid_external(rows))
    allocate(Sa_external(rows,rows))
    read(1,*, iostat=io) Sa_external
    if(io.ne.0) then
        write(0,'(A)')'Error reading Sa file:', Sa_file
        stop
    endif
    read(2,*, iostat=io) hgrid_external
    if(io.ne.0) then
        write(0,'(A)')'Error reading height file:', height_file
        stop
    endif
    close(1)
    close(2)
    hgrid_external = hgrid_external/km2m
    ! asure that the externally provided grid is ordered from ground to top:
    call Sortmat(hgrid_external, Sa_external, rows)
    do ii = 1, sizehgrid
        if (hgrid_external(rows) <= hgrid(ii)) then
            write(0,'(A)')'Error: height of provided Sa not sufficient'
            stop
        else
            do jj = 1, sizehgrid
                ipos  = minloc(abs(hgrid_external-hgrid(ii)),1)
                hipos = hgrid_external(ipos)
                resid  = (hgrid(ii) - hipos)
                ipos2 = minloc(abs(hgrid_external-hgrid(jj)),1)
                hipos2 = hgrid_external(ipos2)
                resid2 = (hgrid(jj)-hipos2)
                if (resid.lt.0.0) then
                    iposp1 = ipos -1
                else
                    iposp1 = ipos +1
                endif
                if (resid2.lt.0.0) then
                    iposp1_2 = ipos2 -1
                else
                    iposp1_2 = ipos2 +1
                endif
                hiposp1 = hgrid_external(iposp1)
                hiposp1_2 = hgrid_external(iposp1_2)
                pipos11 = Sa_external(ipos,  ipos2)
                pipos12 = Sa_external(ipos,  iposp1_2)
                pipos21 = Sa_external(iposp1,ipos2)
                pipos22 = Sa_external(iposp1,iposp1_2)
                temp1   = pipos11 + (pipos21-pipos11) * resid/(hiposp1-hipos)
                temp2   = pipos12 + (pipos22-pipos12) * resid/(hiposp1-hipos)
                Sa(ii,jj) = temp1 + (temp2 - temp1) * resid2/(hiposp1_2-hipos2)
            enddo
       endif
    enddo
  end subroutine read_sa
end module
