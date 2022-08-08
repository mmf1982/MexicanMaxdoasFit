!> This module contains routines to perform the inversion and write inversion products
!>
!> author:  Martina M. Friedrich
!> date:    June, 2014 -- Sept. 2018
!> 
module gaussnewton
  use handle_setupfile, only: output_write_dir, SaM1file
  use constants, only: km2m, m2cm
  implicit none
  contains

  subroutine inversion(&
        Kmat, SCDmeasured, SCD, X_apriori_orig, Sm_inv, Sa_inv, nlayers, &
        nangles,X_orig, scale_fact_plus1)
  !> K is actually reversed, therfore my K is in reality K_trans and vice versa
  real(kind=16), dimension(:,:),allocatable, intent(in) :: Kmat     !> layers_retr x angles_retr
  real(kind=16), dimension(:),  allocatable, intent(in) :: SCDmeasured   !> angles x 1
  real(kind=16), dimension(:),  allocatable, intent(in) :: SCD           !> angles x 1
  real(kind=16), dimension(:),  allocatable, intent(in) :: X_apriori_orig!> layers x 1
  real(kind=16), dimension(:),  allocatable, intent(inout) :: X_orig     !> layers x 1
  real(kind=16), dimension(:,:),allocatable, intent(in) :: Sm_inv         !> angles x angles
  real(kind=16), dimension(:,:),allocatable, intent(in) :: Sa_inv    !> layers_retr x layers_retr
  real(kind=16), intent(in) :: scale_fact_plus1
  real(kind=16), dimension(:), allocatable :: X_apriori
  real(kind=16), dimension(:), allocatable :: X
  integer, intent(in) :: nlayers
  integer, intent(in) :: nangles
  integer :: ii,kk,io, nstart, nlayers_orig
  logical :: SaM1_switch=.false.  ! this used to be read in handle_setupfile
  integer :: upperhalf=0  ! this used to be read in handle_setupfile
  real(kind=16) :: upper_p=1.0, lower_p=1.0   ! this used to be read in handle_setupfile
  real(kind=16), dimension(:,:),  allocatable :: K_trans
  real(kind=16), dimension(:,:),  allocatable :: Gmat
  real(kind=16), dimension(:,:),  allocatable :: G_inv
  real(kind=16), dimension(:),    allocatable :: temp
  real(kind=16), dimension(:),    allocatable :: temp2
  real(kind=16), dimension(:,:),    allocatable :: deltazM
  real(kind=16), dimension(:,:), allocatable :: debug1

  nlayers_orig = size(X_orig)
  nstart = nlayers_orig - nlayers + 1

  allocate(K_trans(nangles,nlayers))
  allocate(Gmat(nlayers,nlayers))
  allocate(G_inv(nlayers,nlayers))
  allocate(temp(nlayers))
  allocate(temp2(nangles))
  allocate(deltazM(nlayers,nlayers))
  allocate(X_apriori(nlayers))
  allocate(X(nlayers))

  X         = X_orig(nstart:)
  X_apriori = X_apriori_orig(nstart:)
  temp(:)   = X(:)-X_apriori(:)
  temp2(:)  = (SCDmeasured(:)-SCD(:))
  K_trans   = transpose(Kmat)

   !-----------------------------------------------------------------------------
   ! 23.11.2014: the below seems to be wrong: (1) it uses L2 instead of L1 
   ! (2) why upper half/lower half? -> should only be different in alpha. 
!   if (SaM1_switch) then
!       deltazM(:,:)=0.0
!       do ii=1,nlayers
!           deltazM(ii,ii)=1.0/X_apriori(ii)!/deltaz(ii)
!       enddo
!       if (upperhalf.gt.nlayers) upperhalf=nlayers
!       Sa_inv(:,:) = 0.0
!       do ii=2,upperhalf-1
!          Sa_inv(ii,ii) = 2.0
!          Sa_inv(ii,ii-1) = -1.0
!          Sa_inv(ii,ii+1) = -1.0
!       enddo
!       Sa_inv(1,1) = 1.0
!       Sa_inv(1,2) = -1.0
!       Sa_inv(upperhalf,upperhalf) = 1.0
!       Sa_inv(upperhalf,upperhalf-1) = -1.0
!       Sa_inv=Sa_inv *upper_p
!       if (upperhalf.lt.nlayers) then
!          if (upperhalf.lt.nlayers-1) then
!             do ii=upperhalf+2,nlayers-1
!                Sa_inv(ii,ii) =2.0
!                Sa_inv(ii,ii-1) = -1.0
!                Sa_inv(ii,ii+1) = -1.0
!             enddo
!             Sa_inv(upperhalf+1,upperhalf+1) = 1.0
!             Sa_inv(upperhalf+1,upperhalf+2) = -1.0
!          endif
!          Sa_inv(nlayers,nlayers)   =  1.0
!          Sa_inv(nlayers,nlayers-1) = -1.0
!          Sa_inv = Sa_inv*lower_p
!       endif
!       Sa_inv=matmul(matmul(deltazM,Sa_inv),deltazM)
!       ! 03.02.2015: Wolfgang wants a posibility to read in an R matrix directly. 
!       if (SaM1file.ne.'zero') then 
!          open(11,file=SaM1file,action='read',status='old',iostat=io)
!          if(io.ne.0) then
!             write(0,'(A)')'  Error opening regulation matrix file, aborting...'
!             write(*,*)'  tried to open ',SaM1file
!             stop
!          end if
!          do kk=1,nlayers
!             read(11,*), Sa_inv(:,kk)
!          enddo
!          close(11)
!       endif 
!   !    Sa_inv=matmul(transpose(Sa_inv),Sa_inv)
!   endif  
  !-----------------------------------------------------------------------------
  ! optimal estimation
  !============================================================
  !Gmat          = matmul(matmul(Kmat,Sm_inv),K_trans)+Sa_inv     ! layers x layers
  Gmat  = matmul(matmul(Kmat,Sm_inv),K_trans)+Sa_inv*(scale_fact_plus1) ! gmat as needed by Rodgers 5.36
  call inverse(Gmat,G_inv,nlayers)
  !X         = X_apriori + matmul(matmul(matmul(G_inv,Kmat),Sm_inv),(temp2+matmul(K_trans,temp)))
  ! this is eq. 5.9 from Rodgers
  X         = X + matmul(G_inv,(matmul(Kmat,matmul(Sm_inv,temp2))-matmul(Sa_inv,temp)))
  ! this is eq. 5.8 from Rodgers.
  ! Rodgers 5.36 is equivalent to 5.8 with the change in the definition of Gmat as implemented above
  X_orig(nstart:) = X

 end subroutine inversion

subroutine write_inversion(Kmat, Sm, Sa_o, nlayers, nangles, xx_o, x_apriori_o, nstart, heights)
  real(kind=16), dimension(:,:),allocatable, intent(inout) :: Kmat     !> layers_retr x angles_retr
  real(kind=16), dimension(:,:),allocatable, intent(in) :: Sm       !> angles x angles
  real(kind=16), dimension(:,:),allocatable, intent(in) :: Sa_o     !> layers_retr x layers_retr
  real(kind=16), dimension(:),allocatable, intent(in) :: xx_o
  real(kind=16), dimension(:),allocatable, intent(in) :: x_apriori_o
  real(kind=16), dimension(:),allocatable, intent(in) :: heights !> layer thickness in km
  integer, intent(in) :: nlayers
  integer, intent(in) :: nangles
  integer, intent(in) :: nstart
  real(kind=16), dimension(:),allocatable :: xx
  real(kind=16), dimension(:),allocatable :: x_apriori
  real(kind=16), dimension(:,:),allocatable :: Sa
  real(kind=16), dimension(:,:),allocatable :: errormat  !> layers x layers
  real(kind=16), dimension(:,:),allocatable :: AvK       !> layers x layers Averaging Kernel
  real(kind=16), dimension(:,:),allocatable :: GAIN      !> layers x angles ???
  real(kind=16), dimension(:,:),allocatable :: Gmat
  real(kind=16), dimension(:,:),allocatable :: G_inv
  real(kind=16), dimension(:,:),allocatable :: K_trans
  real(kind=16) :: dof                  !> degree of freedom 
  character(len=200) :: convspec2,fmtnangles
  real(kind=16), dimension(:,:),  allocatable :: Sm_inv
  real(kind=16), dimension(:,:),  allocatable :: Sa_inv
  integer :: ii

  write(convspec2,*) NLAYERS !"'(", N_USER_VZANGLES,"(f14.5))'"
  write(convspec2,*) trim(convspec2),"es14.5e3,2x)"
  convspec2="("//adjustl(trim(convspec2))
  write(fmtnangles,*) nangles !"'(", N_USER_VZANGLES,"(f14.5))'" !!!
  write(fmtnangles,*) trim(fmtnangles),"e14.5,2x)"
  fmtnangles="("//adjustl(trim(fmtnangles))
  !============================================================
  
  allocate(x_apriori(nlayers))
  allocate(xx(nlayers))
  allocate(Avk(nlayers,nlayers))
  allocate(GAIN(nlayers,nangles))
  allocate(errormat(nlayers,nlayers))
  allocate(Sa(nlayers,nlayers))
  allocate(K_trans(nangles,nlayers))
  allocate(Sm_inv(nangles,nangles))
  allocate(Sa_inv(nlayers,nlayers))
  allocate(Gmat(nlayers,nlayers))
  allocate(G_inv(nlayers,nlayers))
  x_apriori = x_apriori_o(nstart:)
  xx = xx_o(nstart:)
  Sa = Sa_o
  ! Jan 2018: If in log, I need to convert the averaging kernel to lin.
  ! Also, I probably do not want the gamma factor in there, so re-calculate

#ifdef LOGSPACE
  do ii = 1, nlayers   ! ----< is this re-calculatio
    Sa(ii,:) = Sa(ii,:)*exp(x_apriori) !/exp(xx)
    Kmat(ii,:) = Kmat(ii,:)/exp(xx(ii))
  enddo
  do ii = 1, nlayers
    Sa(:,ii) = Sa(:,ii)*exp(x_apriori)!/exp(xx)
  enddo
#endif


#ifdef DETAIL
  open(1,file = trim(output_write_dir) //'Kmat.out',status = 'unknown') !,position="append") ! layers x angles
  write(1,*) "matrix,", nangles
  write(1,'(A)') "Kernel, --, K"
  do ii = 1, nlayers
    write(1,fmtnangles) Kmat(ii,:)
  enddo
  close(1)
  open(1,file = trim(output_write_dir) //'Smmat.out',status = 'unknown') ! layers x angles
  write(1,*) "matrix,", nangles
  write(1,'(A)') "Sm, --, measurement error matrix"
  do ii = 1, nangles
    write(1,fmtnangles) Sm(ii,:)
  enddo
  close(1)
  open(1,file = trim(output_write_dir) //'Samat.out',status = 'unknown') ! layers x angles
  write(1,*) "matrix,", nlayers
  write(1,'(A)') "Sa, --, a-priori covariance matrix"
  do ii = 1, nlayers
#ifdef LOGSPACE
    write(1,convspec2) Sa(ii,:) !*exp(x_apriori)*exp(x_apriori(ii))
#else
    write(1,convspec2) Sa(ii,:)
#endif
  enddo
  close(1)
#endif




  call inverse(Sa, Sa_inv, nlayers)
  call inverse(Sm, Sm_inv, nangles)
  !! For diagnostics
  !------------------------
  !! This first step is the gain matrix, that is "G" in Rodgers
  K_trans = transpose(Kmat)
  Gmat  = matmul(matmul(Kmat,Sm_inv),K_trans)+Sa_inv
  call inverse(Gmat,G_inv,nlayers)
  GAIN = matmul(matmul(G_inv,Kmat),Sm_inv)
  AvK = matmul(GAIN,K_trans)


  dof = 0.0
  open(1,file = trim(output_write_dir) //'avk.dat',status = 'unknown')
  write(1,*) "matrix,", nlayers
  write(1,'(A)') "averaging kernel, --, calculated as Gain * K"
  do ii=1,nlayers
    write(1,*) Avk(:,ii)
    dof = dof + Avk(ii,ii)
  enddo
  close(1)

  open(1,file = trim(output_write_dir) //'dof.dat',status = 'unknown')
  write(1,*)  dof
  close(1)

  open(1,file = trim(output_write_dir) //'GAIN.dat',status = 'unknown')
  write(1,*) "matrix,", nangles
  write(1,'(A)') "Gain matrix, --, calculated as G_inv * K * Sm_inv"
  do ii=1,nlayers
    write(1,fmtnangles) GAIN(ii,:)
  enddo
  close(1)

  errormat=matmul(matmul(GAIN,Sm), transpose(GAIN))
  open(1,file = trim(output_write_dir) //'noiseerror.dat',status = 'unknown')
  write(1,*) "retrieval error"
  write(1,'(A)') "retrieval error, ??, calculated as sqrt(Gain * Sm * Gain)"
  do ii=1,nlayers
#ifdef AEROSOL
    ! in case of AE, I want aod/1km.
    write(1,*) sqrt(errormat(ii,ii))/heights(ii+nstart-1)
#else
    ! in case of TG, I want concentration in cm**3
    write(1,*) sqrt(errormat(ii,ii))/(heights(ii+nstart-1)*km2m* m2cm)
#endif
  enddo
  close(1)

!#ifdef DETAIL
  open(1,file = trim(output_write_dir) //'n_err_cov.dat',status = 'unknown')
  write(1,*) "matrix,", nlayers
  write(1,'(A)') "retrieval error, ??, calculated as Gain * Sm * Gain"
  do ii=1,nlayers
    write(1,*) errormat(ii,:)
  enddo
  close(1)
!#endif

  do ii = 1, nlayers
    Avk(ii,ii) = Avk(ii,ii)-1.
  enddo

  errormat = matmul(matmul(Avk,Sa),transpose(Avk))
  open(1,file = trim(output_write_dir) //'smootherror.dat',status = 'unknown')
  write(1,*) "retrieval error"
  write(1,'(A)') "retrieval error, ??, calculated as sqrt((avk-1) * Sa * (avk-1)T)"
  do ii=1,nlayers
#ifdef AEROSOL
    ! in case of AE, I want aod/1km.
     write(1,*) sqrt(errormat(ii,ii))/heights(ii+nstart-1)
#else
    ! in case of TG, I want concentration in cm**3
    write(1,*) sqrt(errormat(ii,ii))/(heights(ii+nstart-1)*km2m* m2cm)
#endif
  enddo
  close(1)

  !#ifdef DETAIL # this is in partial col
  open(1,file = trim(output_write_dir) //'sm_err_cov.dat',status = 'unknown')
  write(1,*) "matrix,", nlayers
  write(1,'(A)') "retrieval cov smooth error, ??, calculated as (avk-1) * Sa * (avk-1)"
  do ii=1,nlayers
    write(1,*) errormat(ii,:)
  enddo
  close(1)
!#endif
  

  end subroutine write_inversion

  
  subroutine calccost(Sa_inv,Sm_inv,temp,temp2,cost)
    real(kind=16), dimension(:,:), allocatable, intent(in) :: Sa_inv
    real(kind=16), dimension(:,:), allocatable, intent(in) :: Sm_inv
    real(kind=16), dimension(:), allocatable, intent(in) :: temp
    real(kind=16), dimension(:), allocatable, intent(in) :: temp2
    real(kind=16), intent(out):: cost
    cost = dot_product(temp2,matmul(temp2,Sm_inv)) + dot_product(matmul(Sa_inv,temp),temp)
  end subroutine calccost

  subroutine inverse(a2,c,n)
!============================================================
!> Inverse matrix
!> Method: Based on Doolittle LU factorization for Ax=b
!> Alex G. December 2009
!-----------------------------------------------------------
!> input ...
!> a(n,n) - array of coefficients for matrix A
!> n      - dimension
!> output ...
!> c(n,n) - inverse matrix of A
!> comments ...
!> the original matrix a(n,n) will be destroyed 
!> during the calculation
!===========================================================
implicit none 
integer n
real(kind=16) a(n,n), c(n,n), a2(n,n)
real(kind=16) L(n,n), U(n,n), b(n), d(n), x(n)
real(kind=16) coeff
integer i, j, k

a = a2   ! MMF July 2017 fix to keep original matrix as it is
! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse

function eye(Nn)
 integer Nn, ii
 real(kind=16), dimension(Nn,Nn):: eye
 eye = 0.0
   do ii = 1, Nn
      eye (ii, ii) = 1.0
   enddo 
   return
end function



end module
