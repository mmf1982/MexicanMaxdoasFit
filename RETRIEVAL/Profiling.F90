!> Program to run vlidort as desired for the trace gas retrieval
!> ############################################################################
!> Main starting point for MMF fortran part
!> ############################################################################
!> trace gas/ aerosol  retrieval program
!> @author: Martina M. Friedrich
!> @date:   Dec 2013
!> last update:  Oct 2018
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
    use optpropin,           only: calc_gasabs_coef, calc_vlidort_input
    use initialize,          only: VLIDORT_FixIn, VLIDORT_Sup, &
                                   VLIDORT_LinSup, VLIDORT_LinFixIn, &
                                   VLIDORT_LinModIn, VLIDORT_InputStatus, &
                                   rho, gas_beta, aer_beta, abs_aer, &
                                   scatt_aer, cross_sec, scatt_coef, &
                                   initialize_aerosol, initialize_gas, &
                                   initialize_vlidort, Ttgrid, Ppgrid
    use gaussnewton,         only: inversion, calccost, inverse, write_inversion
    use constants,           only: m2cm, km2m, ppb
    use calcSCDK,            only: calcscd, calckernel
    use apriori_sa_read,     only: read_apriori, read_sa
    use sorting,             only: Sort
    use write_output,        only: writescd, writeprofile, writetau, writevcd
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
    integer :: isgood = 2
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
    integer, allocatable, dimension(:) :: skippedrows
    !variables for the inversion step
    real(kind=16) :: I90 = 0.0, I090 = 0.0
    real(kind=16),allocatable,dimension(:)     :: K90
    real(kind=16),allocatable,dimension(:)     :: K090
    real(kind=16),allocatable,dimension(:)     :: I0alpha
    real(kind=16),allocatable,dimension(:)     :: Ialpha
    real(kind=16),allocatable,dimension(:,:)   :: Kalpha
    real(kind=16),allocatable,dimension(:,:)   :: K0alpha
    real(kind=16), allocatable, dimension(:,:) :: Sa
    real(kind=16),allocatable,dimension(:,:)   :: KS, Sm, KS_previous
    real(kind=16), allocatable, dimension(:,:) :: Sa_inv
    real(kind=16), allocatable, dimension(:,:) :: Sm_inv
    real(kind=16), allocatable, dimension(:)   :: middleheights, deltaz
    real(kind=16),dimension(:), allocatable    :: xx, xxold, xx_interm, temp
    real(kind=16),dimension(:), allocatable    :: xx_apriori, dbetadC, beta
    real(kind=16),dimension(:), allocatable    :: temp2
    real(kind=16),dimension(:), allocatable    :: SCD, SCD_previous, SCD_old
    real(kind=16),allocatable,dimension(:)     :: abs_coef, abs_coef2
    real(kind=16),allocatable,dimension(:)     :: Cq, Cq_apriori
    real(kind=16), allocatable, dimension(:)   :: interm
    real(kind=16) :: omega
    real(kind=16) :: costfuncfrac
    real(kind=16) :: h1, h2, fact
    real(kind=16) :: tau, temps, temps2
    real(kind=16) :: dtaudC
    real(kind=16) :: domegadC
    real(kind=16) :: scale_fact_plus1
    real(kind=16) :: convcrit, convfrac, cost = 0.0, cost_old = 1.e12
    real(kind=16) :: start, finish, convcrit_scd
    integer :: num_args, ix, ngr_moms=80, ii, nangle, jj, maxcounter=25
    integer :: n_user_vzangles, nstart
    character(len=500), dimension(:), allocatable :: args
    character(len=500) :: vlidortlogfile = '2p6_VLIDORT_Execution.log'
    character(len=500) :: setupname
    character(len=500) :: mfolder
    logical :: isthere
    logical :: isnotconverged = .true.
    logical :: switch = .true.
    integer :: counter = 0
    character(len=3) :: scannumber
    real(kind=16):: version_number = 1.0
#ifdef TG
    convfrac = 0.015
    costfuncfrac = 1.0001
    scale_fact_plus1 = 256.
#else
    costfuncfrac = 1.0001
    convfrac = 0.015
    scale_fact_plus1 = 4096.
#endif
    
    ! START OF EXECUTABLE SECTION
    ! initialize some quantities:
    call cpu_time(start)
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
        scannumber = trim(args(3))
    else
        scannumber = "-1"
    endif

    open(111,FILE= trim(mfolder) // "version.dat", status="unknown")
    write(111,*) version_number
    close(111)
    !vlidort_config = trim(mfolder) // trim(vlidort_config)
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
        write(0,'(A)')'   for scan: ', scannumber
        stop
    endif

    nstart = nlayers - nlayers_retrieve + 1
    allocate(middleheights(nlayers))
    allocate(deltaz(nlayers))
    allocate(Sa(nlayers_retrieve,nlayers_retrieve))
    allocate(Sa_inv(nlayers_retrieve,nlayers_retrieve))
    allocate(greekmat_total_input(0:ngr_moms,nlayers,MAXSTOKES_SQ))
    allocate(DELTAU_VERT_INPUT(nlayers))
    allocate(xx(nlayers))
    allocate(xxold(nlayers))
    allocate(xx_interm(nlayers))
    allocate(temp(nlayers_retrieve))
    allocate(K90 (nlayers))
    allocate(K090(nlayers))
    allocate(abs_coef(nlayers))
    allocate(abs_coef2(nlayers))
    allocate(xx_apriori(nlayers))
    allocate(dbetadC(0:ngr_moms))
    allocate(beta(0:ngr_moms))
    allocate(CQ(nlayers))
    allocate(CQ_apriori(nlayers))
    allocate(interm(nlayers))
    ! initialize variables to 0
    real0                   = abs(real0)
    abs_coef                = 0.0
    abs_coef2               = 0.0
    middleheights           = 0.0
    deltaz                  = 0.0
    Sa                      = 0.0
    Sa_inv                  = 0.0
    greekmat_total_input    = 0.0
    DELTAU_VERT_INPUT       = 0.0
    xx                      = 0.0
    xxold                   = 0.0
    temp                    = 0.0
    K90                     = 0.0
    K090                    = 0.0
    xx_apriori              = 0.0
    dbetadC                 = 0.0
    beta                    = 0.0
    CQ                      = 0.0
    CQ_apriori              = 0.0

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
    call readfile(input_filename, SCDmeasured, SCDerrmeasured, elevangle, &
                  rel_azimang, szangle, nangle, scannumber, skippedrows)
    do ii=1,NLAYERS
        deltaz(ii) = height_grid(ii-1)-height_grid(ii)
        middleheights(ii) = height_grid(ii)+deltaz(ii)/2.0
    enddo
    ! a-priori file needs to be the filename of the a-priori file
    ! height_file needs to be the filename of the height file.
    ! both of them should be read in in the setup file

#ifdef AEROSOL
    call read_apriori(apriori_file, height_file, aerprofshape, middleheights)
    Cq_apriori = 0.20946*ppb
    Cq = Cq_apriori
#else
    ! if this is gas retrieval aerprofshape comes from the setup_file
    ! The aerosol profile is written directly inside the setup
    call read_apriori(apriori_file, height_file, Cq_apriori, middleheights)
#endif
    ! this sets aer_beta, abs_aer, scatt_aer
    call initialize_aerosol(aerprofshape, taer, waer, gaer, ngr_moms, deltaz)
    ! tpfile needs to be the TP filename. should be read in via setup file
    call initialize_gas(tpfile, middleheights)
#ifdef TG
#ifdef VMRIN
    Cq = Cq_apriori
#else
    xx = Cq_apriori * deltaz * km2m * m2cm
    Cq_apriori = xx/rho*m2cm**3*ppb/ deltaz/(km2m*m2cm)
    Cq = Cq_apriori
#endif
#endif
    allocate (Ialpha (nangle))
    allocate (I0alpha(nangle))
    allocate (SCD(nangle))
    allocate (temp2(nangle))
    allocate (SCD_previous(nangle))
    allocate (SCD_old(nangle))
    allocate (Kalpha (nlayers,nangle))
    allocate (K0alpha(nlayers,nangle))
    allocate (KS(nlayers_retrieve,nangle))
    allocate (KS_previous(nlayers_retrieve,nangle))
    allocate (Sm(nangle,nangle))
    allocate (Sm_inv(nangle,nangle))
    Ialpha      = 0.0
    I0alpha     = 0.0
    temp2       = 0.0
    SCD_previous= 0.0
    SCD_old     = 1.0
    Kalpha      = 0.0
    K0alpha     = 0.0
    SCD         = 0.0
    KS          = 0.0
    Sm          = 0.0
    Sm_inv      = 0.0
    KS_previous = 0.0
    do ii = 1, nlayers
        call calc_gasabs_coef(rho(ii),Cq(ii),cross_sec,abs_coef2(ii))
    enddo

#ifdef TG
        ! As var to vary, use the VCD/ cm² in each layer
    xx = abs_coef2/cross_sec*deltaz*km2m*m2cm
#endif
#ifdef AEROSOL
    xx = (abs_aer+scatt_aer)*deltaz*km2m*m2cm
#endif
    xx_apriori = xx

    if (do_read_Sa) then
      call read_sa(Sa_file, height_file, Sa, middleheights(nstart:nlayers))
      Sa = Sa * Saparam1
#ifdef TG
#ifdef VMRIN
      do jj = 1, nlayers_retrieve
        temps  = rho(jj+nstart-1) / m2cm /m2cm  * deltaz(jj + nstart -1) * km2m/ppb
        do ii = 1, nlayers_retrieve
            temps2 = rho(ii+nstart-1) / m2cm /m2cm  * deltaz(ii + nstart -1) * km2m/ppb
            Sa(ii, jj) = Sa(ii, jj)* temps *temps2
#ifdef LOGSPACE
            Sa(ii, jj) = Sa(ii, jj) / xx_apriori(ii) / xx_apriori(jj)
            if (Sa(ii,jj) < 0) then
                Sa(ii,jj) = 0.0
            endif
#endif
        enddo
      enddo
#endif
#endif
        ! This needs implementing
        write(*,*) "read Sa"
    else
    ! This assumes a correlation length equal to the layer height
#ifdef LOGSPACE
        !write(*,*) middleheights
        do jj =1,nlayers_retrieve
            h1 = middleheights(nstart+jj-1)
            do ii = 1, nlayers_retrieve
                h2 = middleheights(nstart+ii-1)
                fact = ((h2-h1)/scaleheight)**2
                if (ii.eq.jj) then
                    Sa(ii,jj) = Saparam1*Saparam1
                else
                    Sa(ii,jj) = (exp(-log(2.)*fact))*Saparam1*Saparam1  ! Saparam1 replaces 0.5 now
                    !Sa(ii,jj) = (exp(-fact))*Saparam1*Saparam1
                endif
            enddo
        enddo
#else
    ! THIS IS NOT TESTED YET
        do jj = 1, nlayers_retrieve
            interm(jj) = xx_apriori(jj+nstart-1)*xx_apriori(jj+nstart-1)*Saparam1*Saparam1
        enddo
        do jj = 1, nlayers_retrieve
            h1 = middleheights(nstart+jj-1)
            do ii = 1, nlayers_retrieve
                h2 = middleheights(nstart+ii-1)
                fact = ((h2-h1)/scaleheight)**2
                Sa(ii,jj) = sqrt(interm(ii)*interm(jj))*exp(-log(2.)*fact)
            enddo
        enddo
#endif
        !write(*,*) "construct Sa"
    endif

#ifdef LOGSPACE
    do ii = 1,nlayers
        if (xx_apriori(ii).eq.0.0) then
            xx_apriori(ii) = 1.0e-30
        endif
    enddo
    xx_apriori = log(xx_apriori)
    xx = xx_apriori
#endif
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
    do jj = 1, nangle
        Sm(jj,jj) = SCDerrmeasured(jj)*SCDerrmeasured(jj)
    enddo
    call inverse(Sm, Sm_inv, nangle)
    call inverse(Sa,Sa_inv,nlayers_retrieve)

    convergenceloop: do while (isnotconverged.and.(counter.le.maxcounter))
        counter = counter + 1
#ifdef TG
      if (counter.eq.1) then
#endif
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
        enddo layersloop
        !pause
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
        call VLIDORT_LPS_MASTER ( &
                    VLIDORT_FixIn, &
                    VLIDORT_ModIn, &
                    VLIDORT_Sup, &
                    VLIDORT_Out, &
                    VLIDORT_LinFixIn, &
                    VLIDORT_LinModIn, &
                    VLIDORT_LinSup, &
                    VLIDORT_LinOut )
        call VLIDORT_WRITE_STATUS (vlidortlogfile,35, OPENFILEFLAG, VLIDORT_Out%Status )
        I0alpha = VLIDORT_Out%Main%TS_Stokes(1,(/(jj,jj=2,N_USER_VZANGLES)/),1,2)
        I090    = VLIDORT_Out%Main%TS_Stokes(          1,1          ,1,2)
#ifdef LOGSPACE
        K090    = VLIDORT_LinOut%prof%ts_profilewf(1,1:nlayers,1,1,1,2)
        K0alpha = VLIDORT_LinOut%prof%ts_profilewf(1,1:nlayers,1,2:N_USER_VZANGLES,1,2)
#else
        do ii = 1, nlayers
            K090(ii)    = VLIDORT_LinOut%prof%ts_profilewf(1,ii,1,1,1,2)/xx(ii)
            K0alpha(ii,:) = VLIDORT_LinOut%prof%ts_profilewf(1,ii,1,2:N_USER_VZANGLES,1,2)/xx(ii) !(nstart:)
        enddo
#endif

#ifdef TG
      endif
#endif
        ! GAS ABSORPTION
        ! *********************************************************************
        do jj = 1, nlayers
            call calc_vlidort_input(ngr_moms, &
                    scatt_coef(jj), abs_coef2(jj), abs_aer(jj), &
                    scatt_aer(jj), deltaz(jj), gas_beta(:,jj), &
                    aer_beta(:,jj),omega, tau, beta, &
                    dtaudC, domegadC, dbetadC)
            deltau_vert_input(jj) = tau
            if (omega .gt. one-omega_smallnum) omega  = one-omega_smallnum*dble(10.0)
            if (omega .lt. omega_smallnum) omega  = omega_smallnum*dble(10.0)
            omega_total_input(jj) = omega
            greekmat_total_input(0:ngr_moms,jj,1) = beta
            L_deltau_vert_input(1,jj) = dtaudC
            L_omega_total_input(1,jj) = domegadC
            L_greekmat_total_input(1,0:ngr_moms,jj,1) = dbetadC !this was ii[26.01.2015]
            layer_vary_number(jj) = 1
            layer_vary_flag(jj)   = .true.
            !write(*,*) tau, omega
        enddo
        VLIDORT_FixIn%Optical%TS_deltau_vert_input(1:NLAYERS) = &
            deltau_vert_input(1:NLAYERS)
        VLIDORT_FixIn%Optical%TS_greekmat_total_input(0:ngr_moms,1:NLAYERS,:) = &
            greekmat_total_input(0:ngr_moms,1:NLAYERS,:)
        VLIDORT_ModIn%MOptical%TS_omega_total_input(1:NLAYERS)  = &
            omega_total_input(1:NLAYERS)
        VLIDORT_LinFixIn%Optical%TS_L_deltau_vert_input(:,1:NLAYERS)  = &
            L_deltau_vert_input(:,1:NLAYERS)
        VLIDORT_LinFixIn%Optical%TS_L_greekmat_total_input(:,0:ngr_moms,1:NLAYERS,:) = &
            L_greekmat_total_input(:,0:ngr_moms,1:NLAYERS,:)
        VLIDORT_LinFixIn%Optical%TS_L_omega_total_input(:,1:NLAYERS) = &
            L_omega_total_input(:,1:NLAYERS)
        call initialize_vlidort(vlidortlogfile)
        CALL VLIDORT_LPS_MASTER ( &
                VLIDORT_FixIn, &
                VLIDORT_ModIn, &
                VLIDORT_Sup, &
                VLIDORT_Out, &
                VLIDORT_LinFixIn, &
                VLIDORT_LinModIn, &
                VLIDORT_LinSup, &
                VLIDORT_LinOut )
        call VLIDORT_WRITE_STATUS ( &
            vlidortlogfile,35, OPENFILEFLAG, VLIDORT_Out%Status )
        Ialpha =  VLIDORT_Out%Main%TS_Stokes(1,(/(jj,jj=2,N_USER_VZANGLES)/),1,2)
        I90  = VLIDORT_Out%Main%TS_Stokes(1,1 ,1,2)
#ifdef LOGSPACE
        K90(:) = VLIDORT_LinOut%prof%ts_profilewf(1,1:NLAYERS,1,1,1,2)
        Kalpha = VLIDORT_LinOut%prof%ts_profilewf( &
                1,1:nlayers,1,2:N_USER_VZANGLES,1,2)
#else
        do ii = 1, nlayers
            K90(ii)      = VLIDORT_LinOut%prof%ts_profilewf(1,ii,1,1,1,2)/xx(ii)
            Kalpha(ii,:) = VLIDORT_LinOut%prof%ts_profilewf(1,ii,1,2:N_USER_VZANGLES,1,2)/xx(ii) !(nstart:)
        enddo
#endif
    !write(*,*) "I90: ", I90
    !write(*,*) "Ialpha: ", Ialpha
    !write(*,*) "I090: ", I090
    !write(*,*) "I0alpha: ", I0alpha

#ifdef LOGSPACE
        convcrit = sum(abs(exp(xxold)-exp(xx)))/sum(0.5*(exp(xx)+exp(xxold)))
#else
        convcrit = sum(abs(xxold-xx))/sum(0.5*(xx+xxold))
#endif
        convcrit_scd = sum(abs(SCD_old-SCD)/(0.5*(abs(SCD_old)+abs(SCD))))
        if ((((convcrit.le.convfrac).and.convcrit_scd.le.convfrac).and.(scale_fact_plus1<2.5)) &
              .or.counter.ge.maxcounter) then
            isnotconverged = .false.
            if (((convcrit.le.convfrac).and.convcrit_scd.le.convfrac).and.(scale_fact_plus1<2.5)) then
                isgood = 0
            elseif (((convcrit.le.convfrac*5.).and.convcrit_scd.le.convfrac*5.) .and.(scale_fact_plus1<16.5)) then
                isgood = 1
            else
                isgood = 2
            endif
        else

            SCD_old = SCD
            call calcscd(I0alpha,I90,Ialpha,I090,cross_sec,SCD)
            call calckernel(I0alpha, I90, Ialpha, I090, K0alpha, K90, &
                            Kalpha, K090, cross_sec, KS, nstart)
            temp2 = (SCDmeasured-SCD)
            temp = xx(NLAYERS-NLAYERS_RETRIEVE+1:) - &
                   xx_apriori(NLAYERS-NLAYERS_RETRIEVE+1:)
            !call calccost(Sa_inv,Sm_inv,temp,temp2,cost)
            !write(*,*) "conv", convcrit_scd
            !write(*,*) "original:", cost, scale_fact_plus1
            !write(*,*) "   temp2", temp2
            temp2 = temp2
            call calccost(Sa_inv,Sm_inv,temp,temp2,cost) ! JAN2018, Sa was KS_new
            !write(*,*) "convcrit", convcrit_scd
            !write(*,*) "old factor", scale_fact_plus1
            !write(*,*) "new:", cost, "old:", cost_old, "fraction", cost/cost_old
            !write(*,*) "x", xx(NLAYERS-NLAYERS_RETRIEVE+1:)
            !write(*,*) "temp2", temp2
            ! this is the LM check
            if ((cost/cost_old.le.costfuncfrac).or.counter.eq.1.or.scale_fact_plus1>65535) then
                xxold = xx
                cost_old = cost
                KS_previous = KS
                SCD_previous = SCD
                if (scale_fact_plus1>1.5) then
                        scale_fact_plus1 = scale_fact_plus1/4.
                else
                    scale_fact_plus1 = 1.0
                endif
            else
                scale_fact_plus1 = scale_fact_plus1*16.
                xx = xxold
                SCD = SCD_previous
                KS = KS_previous
            endif
            ! write(*,*) "new factor", scale_fact_plus1
            call inversion(KS,SCDmeasured,SCD,xx_apriori,Sm_inv,Sa_inv, &
                    NLAYERS_RETRIEVE,nangle,xx, scale_fact_plus1)

#ifdef AEROSOL
#ifdef LOGSPACE
            !do ii= 1,nlayers  ! vlidort crashes on too large values
            !    xx(ii) = min(xx(ii),3.5)
            !enddo
#else
            do ii=1,nlayers ! vlidort crashes on too large negative values
                if (xx(ii).le.0.0) then
                    xx(ii) = xx(ii)*1.0e-5
                    write(0,*, advance="no") 'warning:',ii,'of xx is negative, set to almost 0'
                    write(0,'(A)')'   for scan: ', scannumber
                endif
            enddo
#endif
#endif
            !write(*,*) "SCD_old", SCD_old
            !write(*,*) "xx: ", xx
            !write(*,*) "SCD: ", SCD
            !write(*,*) "cross_sec", cross_sec
            !pause
            !if (counter.eq.1) then
            !    write(*,*) I90
            !    write(*,*) I090
                ! call writing of scd_ini
            !endif
#ifdef TG
#ifdef LOGSPACE
            abs_coef2 = exp(xx)*cross_sec/deltaz/(km2m*m2cm)
            !do ii=1,nlayers
            !    if (abs_coef2(ii)*deltaz(ii)*km2m*m2cm.gt.2.) then
            !        write(0,*) "warning, big abs_coef2"
            !        write(0,'(A)')'   for scan: ', scannumber
            !        abs_coef2(ii) = 2./km2m/m2cm/deltaz(ii)
            !        xx(ii) = log(abs_coef2(ii)*deltaz(ii)*km2m*m2cm/cross_sec)
            !    endif
            !enddo
#else
            abs_coef2= xx*cross_sec/deltaz/(km2m*m2cm)
#endif
#endif
#ifdef AEROSOL
            scatt_aer=xx/deltaz/(km2m*m2cm)*waer
            abs_aer= xx/deltaz/(km2m*m2cm)*(1.-waer)
#ifdef LOGSPACE
            scatt_aer=exp(xx)/deltaz/(km2m*m2cm)*waer
            abs_aer= exp(xx)/deltaz/(km2m*m2cm)*(1.-waer)
#endif
#endif
        endif
    enddo convergenceloop
    !pause
    !write(*,*) "cross_sec", cross_sec
    call write_inversion(KS, Sm, Sa, nlayers_retrieve, nangle, xx, xx_apriori, &
                         nstart, deltaz)
    call writeprofile(trim(output_write_dir) //"/"// "height_profiles.dat", &
                      xx, deltaz, xx_apriori, height_grid, &
                      rho, Cq_apriori, middleheights, nlayers, nstart, groundoffset, &
                      Ttgrid, Ppgrid)

    call writescd(trim(output_write_dir) //"/"// "SCD_table.dat", elevangle, &
                  SCD, SCDmeasured, SCDerrmeasured, skippedrows)
#ifdef AEROSOL
    call writetau(trim(output_write_dir) //"/"//"aod.dat",xx)
    ! write(*,*) "vcd", sum(abs_coef2/cross_sec*deltaz*km2m*m2cm)
    open(1,file= trim(output_write_dir) //"/"//"vcd.dat")
    write(1,*) sum(abs_coef2/cross_sec*deltaz*km2m*m2cm)
    close(1)
#else
    call writevcd(trim(output_write_dir) //"/"//"vcd.dat",xx)
#endif
    deallocate (Ialpha)
    deallocate (I0alpha)
    deallocate (SCD)
    deallocate (SCDmeasured)
    deallocate (temp2)
    deallocate (Kalpha )
    deallocate (K0alpha)
    deallocate (KS)
    deallocate (Sm)
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
    deallocate (xx)
    deallocate (K90 )
    deallocate (K090)
    deallocate (Sa)
    deallocate (abs_coef)
    deallocate (xx_apriori)
    deallocate (middleheights)
    deallocate (xxold)
    deallocate (beta)
    deallocate (dbetadC)
    deallocate (abs_coef2)
    deallocate (temp)
    call cpu_time(finish)
#ifdef TG
    open(111,FILE= trim(mfolder) //'timing_TG.dat', STATUS='REPLACE') !'UNKNOWN',  position="append")
#endif
#ifdef AEROSOL
    open(111,FILE= trim(mfolder) //'timing_AE.dat', STATUS='REPLACE') !UNKNOWN',  position="append")
#endif
    !write(*,*) "done ", scannumber, counter, scale_fact_plus1, convcrit, convcrit_scd
    !write(*,*) "-----"
    write(111,*) finish-start, counter, scale_fact_plus1, convcrit, convcrit_scd
    close(111)
    open(111,FILE=trim(mfolder)//"flag.dat", STATUS="REPLACE") !UNKNOWN",position="append")
    write(111, *) isgood
    close(111)
    open(111, FILE=trim(mfolder)//"counter.dat", STATUS="REPLACE") !UNKNOWN", position="")
    write(111,*) counter
    close(111)
end program profiling
