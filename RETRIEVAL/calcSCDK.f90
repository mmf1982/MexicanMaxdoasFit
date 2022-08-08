
!> Module to calculate the SCD and the weighting function from VLIDORT output
!>
!> author:   Martina M. Friedrich \n
!> date:     Dec, 2013 \n
!>
!>
!> module to calculate the SCD and the weighting function of the SCD from the 
!> VLIDORT output which is intensity and intensity weighting function
!>
!> in:  I1, I2, I3, I4, crossec \n
!> out: SCD \n

module calcSCDK

implicit none 
 contains

subroutine calcscd(I1,I2,I3,I4,crossec,SCD)

  real(kind=16), intent (in),dimension(:),allocatable  :: I1  !> I_0(alpha)
  real(kind=16), intent (in)  :: I2  !> I  (90)
  real(kind=16), intent (in),dimension(:),allocatable  :: I3  !> I  (alpha)
  real(kind=16), intent (in)  :: I4  !> I_0(90)
  real(kind=16), intent (in)  :: crossec !> what is the unit here? /cm² I think. 
  real(kind=16), intent (inout),dimension(:),allocatable :: SCD !> only I/O for allocate reason

  SCD = log(I1*I2/(I3*I4))/crossec

end subroutine calcscd


subroutine calckernel(I1,I2,I3,I4,K1,K2,K3,K4,crossec,KS, nstart)

  real(kind=16), intent (in),   dimension(:),allocatable  :: I1
  real(kind=16), intent (in)                              :: I2
  real(kind=16), intent (in),   dimension(:),allocatable  :: I3
  real(kind=16), intent (in)                              :: I4
  real(kind=16), intent (in),   dimension(:,:),allocatable:: K1
  real(kind=16), intent (in),   dimension(:),allocatable  :: K2
  real(kind=16), intent (in),   dimension(:,:),allocatable:: K3
  real(kind=16), intent (in),   dimension(:),allocatable  :: K4
  real(kind=16), intent (inout),dimension(:,:),allocatable:: KS
  real(kind=16), intent (in)  :: crossec
  integer, intent(in) :: nstart
  integer :: layernumber,ii,anglenumber, jj, idx

  
  
  layernumber = size(K1,1)
  anglenumber = size(K1,2)
  KS = 0.0
  do ii=1,layernumber-nstart+1
    !KS(ii,:)  = ( K1(ii,:)*I2*I3*I4+I1*K2(ii)*I3*I4-I1*I2*K3(ii,:)*I4-I1*I2*I3*K4(ii) ) / ( I1*I2*I3*I4 ) /crossec
    idx  = ii+nstart-1
    do jj=1,anglenumber
        KS(ii,jj)  = (K1(idx,jj)/I1(jj) + K2(idx)/I2 - K3(idx,jj)/I3(jj) - K4(idx)/I4) /crossec
    enddo
  enddo

  ! TEST
  !KS = 0.0
  !do ii=1,layernumber
  !  KS(ii,:)  = ( I2*I3+K2(ii)*I3-I2*K3(ii,:)-I2*I3 ) / ( I2*I3 ) /crossec
  !enddo
  
  
  ! since K1 and K4 should be 0 for gas, the following calculation replaces the former one for gas
  !do ii=1,layernumber
  !  KS(ii,:)  = (I1*K2(ii)*I3*I4-I1*I2*K3(ii,:)*I4) / ( I1*I2*I3*I4 )/crossec
  !enddo

  ! which is the same as: 
  !KS=0.0
  !do ii = 1,layernumber
  !  do jj =1,anglenumber
  !    KS(ii,jj) = (K2(ii)/I2-K3(ii,jj)/I3(jj))/crossec
  !  enddo 
  !enddo
end subroutine calckernel

end module 
