program test_readininput
! compile with:
! gfortran readinput.f90 test_readinput.f90
! A file test.txt should exist with the following content:
!
! test.txt:
!
! 3
! 1,10,100,1000,10000
! 2.0,20.0,200.0,2000.0,20000.0
! 3e0,3e1,3.e2,3.e3,3.e4
!
! run as:
! ./a.out

use readinput, only: readfile, dscd, dscderr, elevangle, rel_azimang, szangle

implicit none
character(len=20) :: filename="test.txt"



call readfile(filename)

write(*,*) "sum(dscd - expected)       ", sum(dscd - (/1,2,3/))
write(*,*) "sum(dscderr - expected)    ", sum(dscderr - (/10,20,30/))
write(*,*) "sum(elevagle - expected)   ", sum(elevangle - (/100,200,300/))
write(*,*) "sum(rel_azimang - expected)", sum(rel_azimang - (/1000,2000,3000/))
write(*,*) "sum(szangle - expected)    ", sum(szangle - (/10000,20000,30000/))


end program test_readininput
