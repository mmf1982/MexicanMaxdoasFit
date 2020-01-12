
!> This module contains routines to handle the aerosol data from ceilometer measurements
!>
!> author:  Martina M. Friedrich
!> date:    June, 2014
!> 
!> contains:
!> 
!>   readaerosol: reads in the aerosol data. This data is assumed to be written in a single line
!>                where the first entry is a dummy, the next 500 entries are the height in metres
!>                then there are 287 "lines" each with 501 enties. The first entry of
!> each line is the    time with minutes/ seconds as fractions of hour, the next 500 are
!> a measure for reflectancy.
!>
!>   interpaerosol: interpolates the aerosol profile at a specific time to a new height grid
!
!

module aerosolread

implicit none 
 double precision, dimension(:,:) :: aedata(500,287)
 double precision, dimension(:)   :: height(500), time(287)
 
contains

subroutine readaerosol(infilename,actime,hgrid,aerosolprof)
 double precision, dimension(:),allocatable,intent(in) :: hgrid
 double precision, dimension(:),allocatable :: hgrid2
 character(len=500), intent(in) :: infilename
 double precision, intent(in)   :: actime
 integer                     :: io,ii
 double precision, dimension(:) :: aervec(500)
 double precision, dimension(:),allocatable, intent(out) :: aerosolprof
 double precision :: dummy
 
 allocate(hgrid2(size(hgrid)))
 
 hgrid2(:) = hgrid(:) * 1000.0
 height(:)=0.0
 open(1,file=infilename,action='read',status='old',iostat=io)
 if(io.ne.0) then
    write(0,'(A)')'  Error opening aerosol file, aborting...'
    write(*,*)'  tried to open ',infilename 
    stop
 end if
 read(1,'(1(f12.6))',iostat=io,advance="no") dummy
 read(1,'(500(f12.6))',iostat=io,advance="no") height(1:500)
 do ii=1,287
    read(1,'(1(f12.6))',iostat=io,advance="no") time(ii)
    read(1,'(500(f12.6))',iostat=io,advance="no") aedata(1:500,ii)
 enddo
 close(1)
 
 call findtime(actime,aervec)
 call interpaerosol(aervec,height,hgrid2,aerosolprof)
end subroutine readaerosol

subroutine findtime(actime,aervec)
 !> to find the time position and to return a correct aerosol vector
 !> which need to be interpolated

 double precision, intent(in) :: actime
 double precision, intent(out), dimension(:) ::aervec(500)
 integer :: ipos , iposp1
 double precision :: odpos, resid
 
 ipos  = minloc(abs(time - actime),1)
 iposp1= ipos+1
 odpos = time(ipos)
 resid = actime - odpos

 if (resid .lt. 0.0) then 
   iposp1   = ipos -1
 endif
 aervec = aedata(:,ipos) + (aedata(:,iposp1) - aedata(:,ipos)) * resid

end subroutine findtime

subroutine interpaerosol(mdata,h1,h2,mdata2)
 double precision, intent(in), dimension(:) :: mdata(500)  !> aerosol profile at specified time at the heights as in h1
 double precision, intent(in), dimension(:) :: h1(500)     !> heights for original aerosol grid
 double precision, intent(in),  dimension(:),allocatable :: h2 !> new height grid at which data is needed.
 double precision, intent(out), dimension(:),allocatable :: mdata2
 integer :: ii, ipos,iposp1
 double precision :: resid, odpos
 
 allocate(mdata2(size(h2)))
 do ii=1,size(h2)
    if (h2(ii).gt.maxval(h1)) then
       mdata2(ii) =0.0001
    else
       ipos  = minloc(abs(h1-h2(ii)),1)
       iposp1 =ipos+1
       odpos = h1(ipos)
       resid  = h2(ii) -odpos
       if (resid .lt. 0.0) iposp1 = ipos -1
       mdata2(ii) = mdata(ipos) + (mdata(iposp1)-mdata(ipos))*resid
    endif
 enddo
end subroutine interpaerosol


end module
