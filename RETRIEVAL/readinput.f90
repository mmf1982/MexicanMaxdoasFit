
!> Module to read DSCDs and angle data table
!>
!> author:   Martina M. Friedrich \n
!> date:     February 2018

module readinput
  contains
  subroutine readfile(filename, dscd, dscderr, elevangle, raa, szangle, rows, sn, skippedrows)

    implicit none


    integer :: io
    character(len=*), intent(in) :: filename
    double precision, dimension(:,:), allocatable:: datamatrix
    double precision, dimension(:,:), allocatable:: datamatrix_internal
    integer, dimension(:),allocatable,  intent(out):: skippedrows
    double precision, dimension(:), allocatable, intent(out) :: dscd
    double precision, dimension(:), allocatable, intent(out) :: dscderr
    double precision, dimension(:), allocatable, intent(out) :: elevangle
    double precision, dimension(:), allocatable, intent(out) :: raa
    double precision, dimension(:), allocatable, intent(out) :: szangle
    integer, intent(out) :: rows
    integer :: rows_internal, ii
    character(len=*), intent(in) :: sn

    open(1,file = trim(filename),action='read', status='old',iostat=io)
    if(io.ne.0) then
        write(0,'(A)')'Error opening DSCD file: ', filename
        write(0,'(A)')'   for scan: ', sn
        stop
    end if
    read(1,*, iostat=io) rows_internal
    if(io.ne.0) then
        write(0,'(A)')'Error reading number of angles in DSCD file: ', filename
        write(0,'(A)')'   for scan: ', sn
        stop
    end if

    allocate(datamatrix(5, rows_internal))
    allocate(datamatrix_internal(5, rows_internal))
    allocate(skippedrows(rows_internal))
    skippedrows = -1  ! if skippedrows is -1, all ok. If 1, need to insert nan

    read(1,*, iostat=io) datamatrix_internal
    if(io.ne.0) then
        write(0,'(A)')'Error reading  DSCD file:', trim(filename)
        write(0,'(A)')'   for scan: ', sn
        stop
    end if
    close(1)
    rows = 0
    do ii = 1, rows_internal
        if ((isnan(datamatrix_internal(4,ii))).or.(isnan(datamatrix_internal(5,ii)))) then
            skippedrows(ii) = 1
        else
            rows = rows + 1
            datamatrix(:, rows)=datamatrix_internal(:, ii)
        endif
    enddo
    if (rows.eq.0) then
        write(0,*) 'scan ', sn, 'is empty. '
        stop
    endif
    allocate(dscd(rows))
    allocate(dscderr(rows))
    allocate(elevangle(rows))
    allocate(raa(rows))
    allocate(szangle(rows))
    szangle = datamatrix(1, 1:rows)
    raa = datamatrix(2, 1:rows)
    elevangle = datamatrix(3, 1:rows)
    if (minval(elevangle).eq.0.0) then
        ! write(*,*) "elevation angle of 0.0 set to 0.01"
        elevangle = max(0.01, elevangle)
    endif
    ! the elevation angle is defined as measured against the zenith
    ! however, in vlidort, it is defined as from the ground, therefore:
    elevangle = 90. - elevangle
    dscd = datamatrix(4, 1:rows)
    dscderr =datamatrix(5, 1:rows)

  end subroutine readfile
end module
