program test_sort

use sorting

implicit none
double precision :: x(10),y(10), z(10)

x = (/ 10.,9.,8.,7.,6.,5.,4.,3.,2.,1./)

y = (/ 1.,2.,3.,4.,5.,6.,7.,8.,9.,10./)

z = y

call Sort(x,y,z, 10)

write(*,*) x

write(*,*) y

write(*,*) z

end program