! --------------------------------------------------------------------
!> MODULE  Sorting:
!>    This module can sort a set of numbers.  The method used is
!> usually referred to as "selection" method.
! --------------------------------------------------------------------

MODULE  Sorting
  IMPLICIT  NONE
  PRIVATE   :: FindMinimum, Swap
CONTAINS

! --------------------------------------------------------------------
!> INTEGER FUNCTION  FindMinimum():
!>    This function returns the location of the minimum in the section
!> between Start and End.
! --------------------------------------------------------------------

  INTEGER FUNCTION  FindMinimum(x, Start, myend)
    IMPLICIT  NONE
    double precision, DIMENSION(1:), INTENT(IN) :: x
    INTEGER, INTENT(IN)                :: start, myend
    double precision                   :: Minimum
    INTEGER                            :: Location
    INTEGER                            :: i

    Minimum  = x(start)          ! assume the first is the min
    Location = start             ! record its position
    DO i = Start+1, myend        ! start with next elements
      IF (x(i) < Minimum) THEN   !   if x(i) less than the min?
        Minimum  = x(i)          !      Yes, a new minimum found
        Location = i             !      record its position
      END IF
    END DO
    FindMinimum = Location            ! return the position
  END FUNCTION  FindMinimum

! --------------------------------------------------------------------
! SUBROUTINE  Swap():
!    This subroutine swaps the values of its two formal arguments.
! --------------------------------------------------------------------

  SUBROUTINE  Swap(a, b)
    IMPLICIT  NONE
    double precision, INTENT(INOUT) :: a, b
    double precision                :: temp
    temp = a
    a    = b
    b    = temp
  END SUBROUTINE  Swap

! --------------------------------------------------------------------
! SUBROUTINE  Sort():
!    This subroutine receives two arrays x() and y() and sorts y according to
! ascending order of x.
! --------------------------------------------------------------------

   SUBROUTINE  Sort(x, y, z, mysize)
      IMPLICIT  NONE
      double precision, dimension(:), INTENT(inout) :: x    ! sort after
      double precision, dimension(:), intent(inout) :: y ! gets sorted
      double precision, dimension(:), intent(inout) :: z
      INTEGER, INTENT(IN)                   :: mysize
      INTEGER                               :: i
      INTEGER                               :: Location

      DO i = 1, mysize             ! for a strange reason, was mysize-1
         Location = FindMinimum(x, i, mysize)  ! find min from this to last
         call  Swap(x(i), x(location))  ! swap this and the minimum
         call  swap(y(i), y(location))
         call  swap(z(i), z(location))
      END DO
   END SUBROUTINE  Sort

   SUBROUTINE Sortmat(x, ymat, mysize)
    IMPLICIT None
    double precision, dimension(:), intent(inout)   :: x
    double precision, dimension(:,:), intent(inout) :: ymat
    integer, intent(in) :: mysize
    integer :: ii, jj, location
    Do ii = 1, mysize
        Location = FindMinimum(x, ii, mysize)
        call Swap(x(ii), x(location))
        do jj = 1, mysize
            call Swap(ymat(ii,jj),ymat(location, jj))
        enddo
        do jj = 1, mysize
            call Swap(ymat(jj,ii),ymat(jj,location))
        enddo
    enddo
  end subroutine Sortmat
END MODULE  Sorting

