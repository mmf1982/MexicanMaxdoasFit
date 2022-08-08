
!> Module to write profile and dscd SCD_table output

module write_output
  use constants, only: m2cm, km2m, ppb
  ! USE, INTRINSIC :: IEEE_ARITHMETIC
  implicit none
  contains


  subroutine writelog(logfile, SCD, SCDmeasured, SCDerrmeasured, nangle, nlayers, convcrit, counter )
    character(len=*), intent(in):: logfile
    real(kind=16), dimension(:), allocatable, intent(in):: SCD, SCDmeasured, SCDerrmeasured
    integer, intent(in) :: nangle, nlayers, counter
    real(kind=16), intent(in) :: convcrit
    open(17,file=logfile,status="unknown",position="append")
    write(17,*) 'error', sum(abs(SCD-SCDmeasured))
    write(17,*) 'relative error', sum(abs(SCD-SCDmeasured)/SCDerrmeasured)/nangle
    write(17,*) 'converged after', counter, ' steps.'
    write(17,*) '--------------------------------------------------------'
  end subroutine writelog

  
   ! name was SCD_table.dat
  subroutine writescd(scdfile, usedEA, SCD, SCDmeasured, SCDerrmeasured, skippedrows)
    character(len=*), intent(in):: scdfile
    real(kind=16), dimension(:), allocatable, intent(in):: SCD, &
            SCDmeasured, usedEA, SCDerrmeasured
    integer, dimension(:), allocatable, intent(in) :: skippedrows
    integer :: ii, jj
    real(kind=16) :: temp
    open(101,FILE= scdfile, STATUS='UNKNOWN')
    write(101,'(A)', advance="no") "theta/ degree, "
    write(101,'(A)', advance="no") "SCD/ (cm**-2), "
    write(101,'(A)', advance="no") "measured SCD/ (cm**-2), "
    write(101,'(A)') "error SCD/ (cm**-2)"
    write(101,'(A)') "  SCD, (cm**-2), summary of O4 SCDs"
    jj = size(usedEA) + 1
    do ii=1,size(skippedrows)
        if (skippedrows(ii).eq.-1) then
            jj = jj - 1
            write(101,'(f14.4,4(e14.5,1x))') 90.-usedEA(jj), &
                                                SCD(jj), &
                                        SCDmeasured(jj), &
                                        SCDerrmeasured(jj)
        else
            ! temp = IEEE_VALUE(1.0, IEEE_QUIET_NAN)
            write(101,'(A)')  "nan nan nan nan"
        endif
    enddo
    close(101)
  end subroutine writescd

  subroutine writetau(taufile, xx)
    character(len=*), intent(in) :: taufile
    real(kind=16), dimension(:), allocatable, intent(in) :: xx
    open(1,file=taufile)
#ifdef LOGSPACE
    write(1,*) sum(exp(xx))
#else
    write(1,*) sum(xx)
#endif
    close(1)
  end subroutine writetau

  subroutine writevcd(vcdfile,xx)
  character(len=*), intent(in) :: vcdfile
    real(kind=16), dimension(:), allocatable, intent(in) :: xx
    open(1,file=vcdfile)
#ifdef LOGSPACE
    write(1,*) sum(exp(xx))
#else
    write(1,*) sum(xx)
#endif
    close(1)
  end subroutine writevcd

  subroutine writevcdo4(vcdfile, xx)
    ! o4 vcd is the same fore lin or log processing
    character(len=*), intent(in) :: vcdfile
    real(kind=16), dimension(:), allocatable, intent(in) :: xx
    open(1,file= vcdfile)
    write(1,*) sum(xx)
    close(1)
  end subroutine writevcdo4

  ! name was height_profiles.dat
  subroutine writeprofile(profilefile, xx, deltaz, xx_apriori, height_grid, rho, &
        Cq_apriori, middleheights, nlayers, nstart, groundoffset, temp, press)
    character(len=*), intent(in) :: profilefile
    real(kind=16), dimension(:),allocatable, intent(in) :: xx, deltaz,&
        xx_apriori, height_grid, rho, Cq_apriori, middleheights, temp, press
    integer, intent(in):: nlayers, nstart
    real(kind=16), intent(in) :: groundoffset
    integer :: ii
  ! probably I do not need the height_grid if I use middleheights
    OPEN(102, FILE= profilefile, STATUS='UNKNOWN') !, position="append")
#ifdef TG
    write(102,'(A)', advance="no") 'retrieved VMR/ppb, '
    write(102,'(A)', advance="no") 'apriori conentration/ (cm**-3), ' !VMR/ppb
    write(102,'(A)', advance="no") 'height/km, '
    write(102,'(A)', advance="no") 'partial VCD/ (cm**-2), '
    write(102,'(A)', advance="no") 'concentration/ (cm**-3), '
    write(102,'(A)', advance="no") 'delta h/km, '
    write(102,'(A)', advance="no") 'rho/(m**-3), '
    write(102,'(A)') 'air column / (cm**-2)'
    write(102,'(A)') 'VMR , ppb, retrieved VMR '
#endif
#ifdef AEROSOL
    write(102,'(A)', advance="no") 'partial tau/ km, '
    write(102,'(A)', advance="no") 'apriori partial tau/ km, '
    write(102,'(A)', advance="no") 'height/km, '
    write(102,'(A)', advance="no") 'partial tau, '
    write(102,'(A)', advance="no") 'layer height/km, '
    write(102,'(A)', advance="no") 'air column / (cm**-2), '
    write(102,'(A)', advance="no") 'pressure/ hPa, '
    write(102,'(A)', advance="no") 'temperature/ K, '
    write(102,'(A)') 'middle height above sea level/km, '
    write(102,'(A)') "extinction coefficient, --, aerosol optical depth in each layer"
#endif
    do ii = NSTART,NLAYERS
#ifdef AEROSOL
#ifdef LOGSPACE
        write(102,'(9(es14.5e3))') exp(xx(ii))/deltaz(ii), &
                                   exp(xx_apriori(ii))/deltaz(ii),&
                                   (height_grid(ii-1)-height_grid(ii))/2.0+height_grid(ii) ,&
                                   exp(xx(ii)), &
                                   (height_grid(ii-1)-height_grid(ii)), &
                                   rho(ii) * (height_grid(ii-1)-height_grid(ii))*0.1, &
                                   !MMF July 2017: 1/km * 1/m**3 to 1/ cm**2 gives factor 0.1
                                   press(ii), &
                                   temp(ii), &
                                   (height_grid(ii-1)-height_grid(ii))/2.0+height_grid(ii)+groundoffset
#else
        write(102,'(9(es14.5e3))') xx(ii)/deltaz(ii), &
                                   xx_apriori(ii)/deltaz(ii),&
                                   (height_grid(ii-1)-height_grid(ii))/2.0+height_grid(ii) ,&
                                   xx(ii), &
                                   (height_grid(ii-1)-height_grid(ii)), &
                                   rho(ii) * (height_grid(ii-1)-height_grid(ii))*0.1, &
                                   !MMF July 2017: 1/km * 1/m**3 to 1/ cm**2 gives factor 0.1
                                   press(ii), &
                                   temp(ii), &
                                   (height_grid(ii-1)-height_grid(ii))/2.0+height_grid(ii)+groundoffset
#endif
#endif

#ifdef TG
#ifdef LOGSPACE
        write(102,'(8(es14.5e3))') &
        exp(xx(ii))/rho(ii)*(m2cm)**3*ppb/ deltaz(ii)/(km2m*m2cm), &
        Cq_apriori(ii)*rho(ii)/(m2cm**3)/ppb,  &
        middleheights(ii) , &
        exp(xx(ii)), &
        exp(xx(ii))/(deltaz(ii)*km2m*m2cm), &
        (height_grid(ii-1)-height_grid(ii)), &
        rho(ii), &
        rho(ii) * (height_grid(ii-1)-height_grid(ii))*0.1
        !MMF July 2017: 1/km * 1/m**3 to 1/ cm**2 gives factor 0.1
#else
        write(102,'(8(es14.5e3))') &
        xx(ii)/rho(ii)*(m2cm)**3*ppb/ deltaz(ii)/(km2m*m2cm), &
        Cq_apriori(ii)*rho(ii)/(m2cm**3)/ppb,  &
        middleheights(ii) , &
        xx(ii), &
        xx(ii)/(deltaz(ii)*km2m*m2cm), &
        (height_grid(ii-1)-height_grid(ii)), &
        rho(ii), &
        rho(ii) * (height_grid(ii-1)-height_grid(ii))*0.1
#endif
#endif
    enddo
    close(102)
    !write(*,*) xx
    !    pause
  end subroutine writeprofile


 subroutine writedscdini(mfolder, usedela, scd)
 character(len=*), intent(in):: mfolder
 real(kind=16), allocatable, dimension(:), intent(in):: usedela, scd
 integer :: ii
 open(101,FILE= trim(mfolder) //'SCD_ini.dat', STATUS='UNKNOWN')
    write(101,'(a)') "theta/ degree,  scd/ (cm**-2)"
    write(101,'(a)') 'scd, cm**-2, scd that results from a-priori profile'
    do ii=1,size(usedela)
        write(101,'(f14.4,1(e14.5,1x))') usedela(ii), scd(ii)
    enddo
close(101)
end subroutine writedscdini


 end module
