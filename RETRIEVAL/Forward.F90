!> program to run vlidort as desired for the trace gas retrieval
!> trace gas retrieval program
!> author: Martina M. Friedrich
!> date:   Dec 2013
!> last update:  Oct 2017

program profiling

    use vlidort_lps_masters, only: vlidort_lps_master
    use vlidort_pars,        only: MAX_ATMOSWFS, MAXLAYERS, MAXSTOKES_SQ, &
                                   MAXMOMENTS_INPUT, omega_smallnum, one
    use vlidort_io_defs,     only: VLIDORT_Modified_Inputs, VLIDORT_Outputs
    use vlidort_lin_io_defs, only: VLIDORT_LinOutputs
    use vlidort_aux,         only: vlidort_write_status
    use vlidort_l_inputs,    only: vlidort_l_input_master
    use handle_setupfile,    only: readsetupfile, input_filename, NLAYERS,  &
                                   aerprofshape, height_grid, &
                                   taer,gaer,waer,real0, &
                                   nlayers_retrieve, & 
                                   apriori_file, height_file, tpfile, &
                                   do_read_Sa, Saparam1, output_write_dir, &
                                   vlidort_config, groundoffset, scaleheight, &
                                   Sa_file, l_albedo
    use readinput,           only: readfile
    use optpropin,           only: calc_vlidort_input
    use initialize,          only: VLIDORT_FixIn, VLIDORT_Sup, &
                                   VLIDORT_LinSup, VLIDORT_LinFixIn, &
                                   VLIDORT_LinModIn, VLIDORT_InputStatus, &
                                   rho, gas_beta, aer_beta, abs_aer, &
                                   scatt_aer, scatt_coef, &
                                   initialize_aerosol, initialize_gas, &
                                   initialize_vlidort, Ttgrid, Ppgrid
    use constants,           only: m2cm, km2m, ppb
    use calcSCDK,            only: calcscd
    use apriori_sa_read,     only: read_apriori
    use sorting,             only: Sort
    implicit none
    ! variables for setting up vlidort
    real(kind=16),dimension(:,:,:),allocatable :: greekmat_total_input
    real(kind=16) :: l_omega_total_input (max_atmoswfs, maxlayers) = 0.0
    real(kind=16) :: l_deltau_vert_input (max_atmoswfs, maxlayers) = 0.0
    real(kind=16) :: l_greekmat_total_input(&
        max_atmoswfs, 0:maxmoments_input, maxlayers, maxstokes_sq) = 0.0
    integer          :: layer_vary_number(maxlayers)
    real(kind=16) :: omega_total_input(maxlayers) = 0.0
    logical          :: layer_vary_flag( maxlayers )
    logical :: openfileflag
    integer :: isgood = 1
    type(vlidort_modified_inputs)          :: vlidort_modin
    type(vlidort_outputs)                  :: vlidort_out
    type(vlidort_linoutputs)               :: vlidort_linout
    real(kind=16),dimension(:), allocatable    :: deltau_vert_input
    ! measurement variables
    real(kind=16),allocatable,dimension(:)   :: SCDerrmeasured
    real(kind=16),allocatable,dimension(:)   :: SCDmeasured
    real(kind=16),allocatable,dimension(:)   :: elevangle
    real(kind=16),allocatable,dimension(:)   :: szangle
    real(kind=16),allocatable,dimension(:)   :: rel_azimang
    !variables for the inversion step
    real(kind=16) :: I90 = 0.0, I090 = 0.0
    real(kind=16),allocatable,dimension(:)     :: I0alpha
    real(kind=16),allocatable,dimension(:)     :: Ialpha
    real(kind=16), allocatable, dimension(:)   :: middleheights, deltaz
    real(kind=16),dimension(:), allocatable    :: temp
    real(kind=16),dimension(:), allocatable    :: dbetadC, beta
    real(kind=16),dimension(:), allocatable    :: temp2
    real(kind=16),dimension(:), allocatable    :: SCD_gas, SCD_O4
    real(kind=16),allocatable,dimension(:)     :: abs_coef, abs_coef2, abs_gas, abs_O4
    real(kind=16),allocatable,dimension(:)     :: Cq_gas, Cq_O4
    integer, allocatable, dimension(:) :: skippedrows
    real(kind=16) :: omega
    real(kind=16) :: h1, h2, fact
    real(kind=16) :: tau, temps, temps2
    real(kind=16) :: dtaudC
    real(kind=16) :: cseco4
    character(len=100) :: cseco4_as_str
    real(kind=16) :: csectg
    character(len=100) :: csectg_as_str
    real(kind=16) :: domegadC
    real(kind=16) :: start, finish, convcrit_scd
    integer :: num_args, ix, ngr_moms=80, ii, nangle, jj
    integer :: n_user_vzangles, nstart
    character(len=500), dimension(:), allocatable :: args
    character(len=500) :: vlidortlogfile = '2p6_VLIDORT_Execution.log'
    character(len=500) :: setupname
    character(len=500) :: mfolder
    logical :: isthere
    character(len=3) :: scannumber
    num_args = command_argument_count() ! extract the MMF folder name
    allocate(args(num_args))
    do ix = 1, num_args
         call get_command_argument(ix, args(ix))
    end do
    if (num_args > 0) then
        mfolder = trim(args(1)) //"/"
    else
        write(0,"(A)") "first argument passed needs to be folder MMF folder"
        stop
    endif
    if (num_args > 1) then
        setupname = trim(args(2))
    else
        write(0,"(A)") "second argument passed needs to be setup file path"
        stop
    endif

    if (num_args > 2) then
        cseco4_as_str = trim(args(3))
        read(cseco4_as_str , *) cseco4
    else
        write(0,"(A)") "third argument passed needs to be O4 xs"
        stop
    endif

    if (num_args > 3) then
        csectg_as_str = trim(args(4))
        read(csectg_as_str , *) csectg
    else
        write(0,"(A)") "forth argument passed needs to be TG xs"
        stop
    endif
    !write(*,*) csectg, cseco4
    
    
    call readsetupfile(trim(setupname))
    if ((gaer<-1) .or. (gaer>1)) then
        write(*,*) "ERROR in scan ", scannumber, " asym out of range [-1,1]: ", gaer
        stop
    elseif ((gaer>-0.01) .and. (gaer<0.01)) then
        gaer = 0.01
        write(*,*) "WARNING in scan ", scannumber, " asym almost 0, set to 0.01"
    elseif (gaer>0.99) then
        write(*,*) "WARNING in scan ", scannumber, " asym almost 1, set to 0.99"
        gaer = 0.99
    elseif (gaer<-0.99) then
        write(*,*) "WARNING in scan ", scannumber, " asym almost -1, set to -0.99"
        gaer = -0.99
    endif
    if ((waer < 0) .or. (waer >1.0)) then
        write(*,*) "ERROR in scan ", scannumber, "waer out of range [0,1]: ", waer
        stop
    elseif (waer<0.01) then
        write(*,*) "WARNING in scan ", scannumber, " waer almost 0, set to 0.01"
        waer = 0.01
    elseif (waer>0.99) then
        write(*,*) "WARNING in scan ", scannumber, " waer almost 1, set to 0.99"
        waer = 0.99
    endif
    inquire(FILE= vlidort_config, EXIST=isthere)
    if (.not.isthere) then
        write(0,'(A)', advance="no")'Error, vlidort config file not found: ', vlidort_config
        stop
    endif

    nstart = nlayers - nlayers_retrieve + 1
    allocate(middleheights(nlayers))
    allocate(deltaz(nlayers))
    allocate(greekmat_total_input(0:ngr_moms,nlayers,MAXSTOKES_SQ))
    allocate(DELTAU_VERT_INPUT(nlayers))
    allocate(temp(nlayers_retrieve))
    allocate(abs_coef(nlayers))
    allocate(abs_coef2(nlayers))
    allocate(abs_gas(nlayers))
    allocate(abs_O4(nlayers))
    allocate(dbetadC(0:ngr_moms))
    allocate(beta(0:ngr_moms))
    allocate(CQ_gas(nlayers))
    allocate(CQ_O4(nlayers))
    ! initialize variables to 0
    real0                   = abs(real0)
    abs_coef                = 0.0
    abs_coef2               = 0.0
    middleheights           = 0.0
    deltaz                  = 0.0
    greekmat_total_input    = 0.0
    DELTAU_VERT_INPUT       = 0.0
    temp                    = 0.0
    dbetadC                 = 0.0
    beta                    = 0.0
    CQ_gas                  = 0.0
    CQ_O4                   = 0.0
    call VLIDORT_L_INPUT_MASTER ( &
      vlidort_config,     & ! Input
      VLIDORT_FixIn,      & ! Outputs
      VLIDORT_ModIn,      & ! Outputs
      VLIDORT_LinFixIn,   & ! Outputs
      VLIDORT_LinModIn,   & ! Outputsq
      VLIDORT_InputStatus ) ! Outputs
    VLIDORT_FixIn%Cont%TS_NLAYERS = NLAYERS       ! to overwrite read input
    VLIDORT_ModIn%MUserVal%TS_USER_LEVELS=NLAYERS ! to overwrite read input
    VLIDORT_ModIn%MCont%TS_ngreek_moments_input = ngr_moms
    call readfile(trim(mfolder)//input_filename, SCDmeasured, SCDerrmeasured, elevangle, &
                  rel_azimang, szangle, nangle, scannumber, skippedrows)
    do ii=1,NLAYERS
        deltaz(ii) = height_grid(ii-1)-height_grid(ii)
        middleheights(ii) = height_grid(ii)+deltaz(ii)/2.0
    enddo
    ! a-priori file needs to be the filename of the a-priori file
    ! height_file needs to be the filename of the height file.
    ! both of them should be read in in the setup file
    call read_apriori(trim(trim(height_file) //"/aerosol_apriori.txt"), &
        trim(trim(height_file) //'/apriori_heights.txt'), aerprofshape, middleheights)
    Cq_O4 = 0.20946
    call read_apriori(trim(trim(height_file) //"/gas_apriori.txt"), &
        trim(trim(height_file) //'apriori_heights.txt'), Cq_gas, middleheights)
    ! this sets aer_beta, abs_aer, scatt_aer
    call initialize_aerosol(aerprofshape, taer, waer, gaer, ngr_moms, deltaz)
    ! tpfile needs to be the TP filename. should be read in via setup file
    call initialize_gas(tpfile, middleheights)
#ifdef VMRIN
    abs_gas = csectg * rho/dble(1000000.0) * Cq_gas/dble(1000000000.0)
#else
    ! in this case, the read in Cq_apriori is concentration in cm**3.
    abs_gas = csectg * Cq_gas
    !write(*,*) "xs ",csectg
    !write(*,*) "Cq", Cq_gas
#endif
    allocate (Ialpha (nangle))
    allocate (I0alpha(nangle))
    allocate (SCD_gas(nangle))
    allocate (SCD_O4(nangle))
    allocate (temp2(nangle))
    Ialpha      = 0.0
    I0alpha     = 0.0
    temp2       = 0.0
    abs_O4      = cseco4 * rho/dble(1000000.0)*rho/dble(1000000.0) * (Cq_O4*Cq_O4)
    !write(*,*) abs_O4
    !pause
    ! May 2018 add surface albedo
    VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO = l_albedo
    ! setup the angles, they can be outside the convergence loop now
    VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT = 0.0
    n_user_vzangles                               = nangle+1
    VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES     = n_user_vzangles
    VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(1) = real0
    call sort(elevangle, SCDmeasured, SCDerrmeasured,nangle)
    VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(2:n_user_vzangles) =  &
                        elevangle
    VLIDORT_ModIn%MSunRays%TS_SZANGLES = sum(szangle)/nangle
        ! below, I should make sure that I only have rel_azimang gt 0!?
        ! also, I should make sure that this is actually the right one now
        ! w.r.t all the adding and subtracting of 180 degree
    VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS = sum(rel_azimang)/nangle

    ! NO GAS ABSORPTION
    ! *********************************************************************
    layersloop: do ii = 1, NLAYERS
        call calc_vlidort_input(ngr_moms, &
                scatt_coef(ii),abs_coef(ii),abs_aer(ii), &
                scatt_aer(ii),deltaz(ii),gas_beta(:,ii), &
                aer_beta(:,ii),omega,tau,beta, &
                dtaudC, domegadC,dbetadC)
        deltau_vert_input(ii)                     = tau
        if (omega .gt. one-omega_smallnum) omega  = one-omega_smallnum
        if (omega .lt. omega_smallnum) omega      = omega_smallnum
        omega_total_input(ii)                     = omega
        greekmat_total_input(0:ngr_moms,ii,1)     = beta
        L_deltau_vert_input(1,ii)                 = dtaudC
        L_greekmat_total_input(1,0:ngr_moms,ii,1) = dbetadC
        L_omega_total_input(1,ii)                 = domegadC
	!write(*,*) tau, omega
    enddo layersloop
    VLIDORT_FixIn%Optical%TS_deltau_vert_input(1:NLAYERS) = deltau_vert_input
    VLIDORT_FixIn%Optical%TS_greekmat_total_input(0:ngr_moms,1:NLAYERS,:) = &
            greekmat_total_input
    VLIDORT_ModIn%MOptical%TS_omega_total_input(1:NLAYERS) = &
        omega_total_input(1:NLAYERS)
    VLIDORT_LinFixIn%Optical%TS_L_deltau_vert_input  = L_deltau_vert_input
    VLIDORT_LinFixIn%Optical%TS_L_greekmat_total_input(:,0:ngr_moms,1:NLAYERS,:) = &
        L_greekmat_total_input(:,0:ngr_moms,1:NLAYERS,:)
    VLIDORT_LinFixIn%Optical%TS_L_omega_total_input   = L_omega_total_input
    ! vlidortlogfile should be a name read in by the setup file
    call initialize_vlidort(vlidortlogfile)
    call VLIDORT_LPS_MASTER (VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup, &
                VLIDORT_Out, VLIDORT_LinFixIn, VLIDORT_LinModIn, VLIDORT_LinSup, &
                VLIDORT_LinOut )
    call VLIDORT_WRITE_STATUS (vlidortlogfile,35, OPENFILEFLAG, VLIDORT_Out%Status)
    I0alpha = VLIDORT_Out%Main%TS_Stokes(1,(/(jj,jj=2,N_USER_VZANGLES)/),1,2)
    I090    = VLIDORT_Out%Main%TS_Stokes(          1,1          ,1,2)
    ! GAS ABSORPTION
    ! *********************************************************************
    do ii = 1,2
        if (ii == 1) then
            abs_coef2 = abs_gas
        else
            abs_coef2 = abs_O4
        endif
        !write(*,*) "IN GAS ABS"
        do jj = 1, nlayers
            call calc_vlidort_input(ngr_moms, scatt_coef(jj), abs_coef2(jj), abs_aer(jj), &
                    scatt_aer(jj), deltaz(jj), gas_beta(:,jj), aer_beta(:,jj),omega, tau, beta, &
                    dtaudC, domegadC, dbetadC)
            !if (ii == 2) then
            !    write(*,*) jj, tau
            !endif
            deltau_vert_input(jj) = tau
            if (omega .gt. one-omega_smallnum) omega  = one-omega_smallnum*dble(10.0)
            if (omega .lt. omega_smallnum) omega  = omega_smallnum*dble(10.0)
            omega_total_input(jj) = omega
            greekmat_total_input(0:ngr_moms,jj,1) = beta
            L_deltau_vert_input(1,jj) = dtaudC
            L_omega_total_input(1,jj) = domegadC
            L_greekmat_total_input(1,0:ngr_moms,jj,1) = dbetadC
            layer_vary_number(jj) = 1
            layer_vary_flag(jj)   = .true.
        enddo
        !pause
        VLIDORT_FixIn%Optical%TS_deltau_vert_input(1:NLAYERS) = deltau_vert_input(1:NLAYERS)
        VLIDORT_FixIn%Optical%TS_greekmat_total_input(0:ngr_moms,1:NLAYERS,:)=greekmat_total_input(0:ngr_moms,1:NLAYERS,:)
        VLIDORT_ModIn%MOptical%TS_omega_total_input(1:NLAYERS)=omega_total_input(1:NLAYERS)
        VLIDORT_LinFixIn%Optical%TS_L_deltau_vert_input(:,1:NLAYERS)=L_deltau_vert_input(:,1:NLAYERS)
        VLIDORT_LinFixIn%Optical%TS_L_greekmat_total_input(:,0:ngr_moms,1:NLAYERS,:) = &
            L_greekmat_total_input(:,0:ngr_moms,1:NLAYERS,:)
        VLIDORT_LinFixIn%Optical%TS_L_omega_total_input(:,1:NLAYERS) = &
            L_omega_total_input(:,1:NLAYERS)
        call initialize_vlidort(vlidortlogfile)
        CALL VLIDORT_LPS_MASTER (VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup, &
            VLIDORT_Out, VLIDORT_LinFixIn, VLIDORT_LinModIn, VLIDORT_LinSup, &
            VLIDORT_LinOut )
        call VLIDORT_WRITE_STATUS (vlidortlogfile,35, OPENFILEFLAG, VLIDORT_Out%Status )
        Ialpha =  VLIDORT_Out%Main%TS_Stokes(1,(/(jj,jj=2,N_USER_VZANGLES)/),1,2)
        I90  = VLIDORT_Out%Main%TS_Stokes(1,1 ,1,2)
        if (ii == 1) then
            call calcscd(I0alpha,I90,Ialpha,I090,csectg,SCD_gas)
        else
            call calcscd(I0alpha,I90,Ialpha,I090,cseco4,SCD_O4)
        endif
    enddo
    open(101,FILE= trim(mfolder)//"/dscd.txt", STATUS='UNKNOWN')
    !write(101,'(A)', advance="no") "theta/ degree, "
    !write(101,'(A)', advance="no") "tg SCD/ (molec/ cm**-2) "
    !write(101,'(A)', advance="yes") "O4 SCD/ (molec**2/ cm**-5) "
    ! this is now a specific format so it works for my forward model table
    ! sza, raa and then elevangle measurements
    jj = size(elevangle) + 1
    ! write(*,*) "raa", rel_azimang
    write(101, *) szangle(1)
    write(101, *) rel_azimang(1)
    do ii=1,size(elevangle)
        jj = jj - 1
        write(101,'(3(e14.5,1x))') 90.-elevangle(jj), SCD_gas(jj), SCD_O4(jj)
    enddo
    close(101)
    open(101,FILE= trim(mfolder)//"/profile.txt", STATUS='UNKNOWN')
    write(101,'(A)', advance="no") "height, "
    write(101,'(A)', advance="no") "partial aod "
    write(101,'(A)', advance="no") "aod per km "
    write(101,'(A)', advance="yes") "concentration"
    do ii=1, nlayers
        jj = jj - 1
        write(101,'(4(e14.5,1x))') middleheights(ii), &
            (abs_aer(ii)+scatt_aer(ii))*deltaz(ii)*100000.0, &
            (abs_aer(ii)+scatt_aer(ii))*100000.0, Cq_gas(ii)
    enddo
    close(101)
    !write(*,*) "tau", sum((abs_aer+scatt_aer)*deltaz*100000.0)
    !do ii = NLAYERS,1, -1
    !    write(*,'(2(e12.4))') (abs_aer(ii)+scatt_aer(ii))*100000.0, (height_grid(ii-1)-height_grid(ii))/2.0+height_grid(ii)
    !enddo
    !write(*,*) "I90",I90,I090
    !write(*,*) "Ialpha", Ialpha
    !write(*,*) "I0alpha",I0alpha
    
    deallocate (Ialpha)
    deallocate (I0alpha)
    deallocate (SCD_gas)
    deallocate (SCD_O4)
    deallocate (SCDmeasured)
    deallocate (temp2)
    deallocate (SCDerrmeasured)
    deallocate (rho)
    deallocate (deltaz)
    deallocate (gas_beta)
    deallocate (scatt_coef)
    deallocate (abs_aer)
    deallocate (scatt_aer)
    deallocate (aer_beta)
    deallocate (greekmat_total_input)
    deallocate (DELTAU_VERT_INPUT)
    deallocate (abs_coef)
    deallocate (middleheights)
    deallocate (beta)
    deallocate (dbetadC)
    deallocate (abs_coef2)
    deallocate (temp)
end program profiling
