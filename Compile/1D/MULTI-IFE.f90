!======================================================================!
!                                                                      !
! M    M  U    U  LL      TTTTTT  IIIIII        IIIIII  FFFFFF  EEEEEE !
! MM  MM  U    U  LL      TTTTTT  IIIIII        IIIIII  FFFFFF  EEEEEE !
! MMMMMM  U    U  LL        TT      II            II    F       E      !
! M MM M  U    U  LL        TT      II            II    F       E      !
! M    M  U    U  LL        TT      II   ======   II    FFFFF   EEEEE  !
! M    M  U    U  LL        TT      II            II    F       E      !
! M    M  U    U  LL        TT      II            II    F       E      !
! M    M  UUUUUU  LLLLLL    TT    IIIIII        IIIIII  F       EEEEEE !
! M    M   UUUU   LLLLLL    TT    IIIIII        IIIIII  F       EEEEEE !
!                                                                      !
!======================================================================!
!                                                                      !
!     MULTI-IFE - A One-Dimensional Computer Code for Inertial Fusion  ! 
!     Energy (IFE) Target Simulations                                  !
!                                                                      !
!     R.Ramis and J. Meyer-ter-vehn                                    !
!                                                                      !
!     to be submitted to Computer Physics Communications               !
!                                                                      !
!     Version March 4, 2015                                            !
!                                                                      !
!======================================================================!
!                                                                      !
!     unidimensional radiation hydrodynamic simulation code including  !
!                                                                      !
!     - planar, cylidrical, or spherical geometry                      !
!     - multi-layer multi-material shells                              !
!     - laser deposition by either:                                    !
!         - WKB aproximation with normal incidence and a fraction of   !
!           energy coupled to the critical density                     !
!         - Maxwell equation solver (only for planar geometry)         !
!         - WKB with 3D ray tracing (only for spherical geometry)      !
!     - implicit hydrodynamics with separate ion and electron temp.    !
!     - multigroup non-lte radiation transport                         !
!     - flux limited electron thermal diffusion                        !
!     - ion thermal diffusion                                          !
!     - tabulated matter properties (EOS, opacities, emissivities,     !
!       and ionization)                                                !
!     - DT fusion package with non-local deposition (alpha particle    !
!       diffusion)                                                     !
!     - electron collisionality modelled by either:                    !
!         - classical plasma model                                     !
!         - electron-phonon interaction                                !
!         - Drude-Sommerfeld model                                     !
!                                                                      !
!======================================================================!
      module multidata
!======================================================================!
!                                                                      !
!     definitions of data structures used in all parts of the code     !
!                                                                      !
!======================================================================!
!-----mathematical constants--------------------------------------------
      real*8, parameter :: cpi   = 3.1415926535897931d0 ! pi 
!-----physical constants--(from S.Atzeni and J.Meyer-ter-Vehn book)-----
      real*8, parameter :: pm    = 1.6726e-24  ! proton mass (g)
      real*8, parameter :: bolz  = 1.6022e-12  ! Boltzmann C. (cgs/eV)
      real*8, parameter :: em    = 9.1096e-28  ! electron mass (g)
      real*8, parameter :: eq    = 4.8033e-10  ! electron charge (e.u)
      real*8, parameter :: c     = 2.9979e+10  ! light velocity (cm/s)
      real*8, parameter :: hb    = 1.0546e-27  ! Planck constant (cgs)
!-----derived constants-------------------------------------------------
      real*8, parameter :: hh    = 2*cpi*hb    ! Planck constant (cgs)
      real*8, parameter :: sigma = (2*cpi**5*bolz**4)/(15*c**2*hh**3) 
                                               ! Stefan-Boltzmann C.
!-----memory control----------------------------------------------------
      integer           :: memory = 0          ! Allocated memory
!=======================================================================
!-----equation of state table (inverted: P(e,rho) and T(e,rho))
      type eostable
      integer                       :: nr=0
      integer                       :: ne=0
      real*8, dimension(:), pointer :: rho => null()
      real*8, dimension(:), pointer :: e0 => null()
      real*8, dimension(:), pointer :: de => null()
      real*8, dimension(:), pointer :: p => null()
      real*8, dimension(:), pointer :: t => null()
      end type eostable
!=======================================================================
!-----group opacity table or Z table
      type grouptable
      integer                       :: nr=0
      integer                       :: nt=0
      real*8                        :: fmin=0,fmax=0
      real*8, dimension(:), pointer :: rho => null()
      real*8, dimension(:), pointer :: t => null()
      real*8, dimension(:), pointer :: o => null()
      end type grouptable
!=======================================================================
!-----multi group opacity table
      type multitable
      integer                       :: ng=0
      type(grouptable), dimension(:), pointer :: tables => null()
      end type multitable
!=======================================================================
!-----material struture: atomic number, mass numbers, and tables
      type typematerial
      character(len=80) :: name
      real*8            :: zi,ai
      character(len=80) :: eeos_file
      integer           :: eeos_id
      type(eostable)    :: eeos
      character(len=80) :: ieos_file
      integer           :: ieos_id
      type(eostable)    :: ieos
      character(len=80) :: z_file
      integer           :: z_id
      type(grouptable)  :: z
      character(len=80) :: planck_file
      integer           :: planck_id
      type(multitable)  :: planck
      character(len=80) :: ross_file
      integer           :: ross_id
      type(multitable)  :: ross
      character(len=80) :: eps_file
      integer           :: eps_id
      type(multitable)  :: eps
      end type typematerial
!=======================================================================
!-----Definition of a WKB laser pulse (parallel/convergent rays)
      type typesimlaserwkb
      integer :: enable=0  ! main flag
      integer :: inter     ! laser direction (+1/-1 from RHS/LHS)
      real*8  :: pimax     ! pulse power
      real*8  :: pitime    ! pulse time
      real*8  :: wl        ! laser wavelength
      integer :: itype     ! pulse type 1 - sin**2
                           !            2 - rectangular
                           !            4 - tabulated
      real*8  :: delta     ! fraccional absorption at crit. density
      end type typesimlaserwkb
!=======================================================================
!-----Definition of a WKB laser pulse (spherical 3D with refraction)
      type typesimlaser3d
      integer :: enable=0  ! main flag
      real*8  :: pimax     ! pulse power
      real*8  :: pitime    ! pulse time
      real*8  :: wl        ! laser wavelength
      integer :: itype     ! pulse type 1 - sin**2
                           !            2 - rectangular
                           !            4 - tabulated
      integer :: nr        ! number of rays
      real*8  :: rmax      ! beam radius
      real*8  :: fwhm      ! beam fwhm
      real*8  :: bexp      ! supergauss exponent
      end type typesimlaser3d
!=======================================================================
!-----Definition of a laser wave (Maxwell solver)
      type typesimlasermaxwell
      integer :: enable=0  ! main flag
      integer :: inter     ! laser direction (+1/-1 from RHS/LHS)
      real*8  :: pimax     ! pulse intensity (in propagation dir.)
      real*8  :: pitime    ! pulse time
      real*8  :: wl        ! laser wavelength
      integer :: itype     ! pulse type 1 - sin**2
                           !            2 - rectangular
                           !            4 - tabulated
      integer :: idep=1    ! number of cell subdivisions
      real*8  :: angle     ! incidence angle
      integer :: pol       ! polarization: 1/2 for p/s
      end type typesimlasermaxwell
!=======================================================================
!-----definition of a external pulse (used when itype==4)
      type typesimpulse
      integer :: ntab=0        ! number of points
      real*8  :: ttab(100)     ! times
      real*8  :: ptab(100)     ! values
      integer :: mode          ! 1 - linear interp. in value**expo
                               ! 2 - linear interp. in log(value)
                               ! 3 - linear interp. in exp(value)
      real*8  :: expo          ! used for mode=1
      end type typesimpulse
!=======================================================================
!-----constant simulation parameters.
      type typesimparam
  
!-----wfyuan,2020/10/27, index for genetic algorithm input/output file	  
      character(32)   :: igafile
      real*8          :: rhorDT
      real*8          :: alphaDT
      real*8          :: timestag
      real*8          :: timean
      real*8          :: vmean
      real*8          :: ifarmax
      real*8          :: rhomax
      real*8          :: pmaxdt
      real*8          :: Ekimp
      real*8          :: Edelas
      real*8, pointer :: gax(:)
!-----wfyuan,2020/10/27, index of file number for genetic algorithm	
  
!-----geometry type: 1/2/3 for planar/cylindrical/spherical
      integer :: igeo 
!-----boundary condition on the first/last interface 0/1 for free/rigid
      integer :: iright, ileft
!-----number of cells with DT fuel (0 inhibes nuclear reactions)
      integer :: nfuel            
!-----maximum number of time steps
      integer :: nexit        
!-----partition of input energy between processes 1/0 for optimum/fixed
      integer :: iwctrl 
!-----enable/disable radiation transport, hydrodynamics, and ion heat (1/0)
      integer :: iradia,ihydro,iheation   
!-----enable/disable negative pressures (1/0)
      integer :: inegpre
!-----model for electron collisions and heat condutivity
!     model=0  - classical plasma (Braginskii)
!     model=1  - electron-phonon interaction for low temperatures and
!                classical plasma for high temperatures with a smooth
!                interpolation
!     model=2  - Drude-Sommerfeld model with smooth interpolation to
!                classical plasma
!     fheat, fei, flaser - multipliers used in model=2, for heat
!                conduction, e-i transfer, and laser deposition, resp.
!     zmin     - minimum value of ionization
!     flf      - flux limit factor
      integer :: model               
      real*8  :: fheat,fei,flaser,zmin,flf
!-----boundary conditions for radiation
      real*8  :: alphal,alphar   ! reflection factors on L/R boundaries
      real*8  :: betal,betar     ! incoming power factor on L/R Bound.
      integer :: itype           ! pulse type 1 - sin**2
                                 !            2 - rectangular
                                 !            4 - tabulated
      real*8  :: tau             ! radiation pulse time
      real*8  :: trad            ! radiation pulse temperature
!-----initial position of leftmost interface
      real*8  :: xmin  
!-----number of cells
      integer :: n
!-----mass of cells and atomic mass number
      real*8, dimension(:), pointer :: dm,ai
!-----for each cell number, material index in array 'mat'
      integer, dimension(:), pointer :: mid     
!-----frequency boundary for radiation groups
      real*8, dimension(:), pointer :: frei => null()
!-----number of radiation groups
      integer :: ng  
!-----number of substeps in a full time step 
      integer :: nsplit,msplit        ! required value and used value
!-----interval to be simulated
      real*8  :: texit
!-----output control
      real*8  :: dt_bout, dt_aout     ! binary/ascii output time interval
      integer :: ns_bout, ns_aout     ! binary/ascii output step interval
      integer :: nreduce              ! output one of each nreduce cell 
                                      ! or iterface values
!-----write radiation spectrum at these interfaces
      integer :: irad_left            ! interface number
      integer :: irad_right           ! interface number
!-----time step control
      real*8  :: dtmin                ! minimum dt
      real*8  :: dtmax                ! maximum dt
      real*8  :: dtinit               ! initial dt
      real*8  :: dtrvar               ! relative change of density
      real*8  :: dttevar              ! relative change of elec. temp.
      real*8  :: dttivar              ! relative change of ion temp.
      real*8  :: dtbreak              ! break step change
      real*8  :: dtfactor             ! dt increment factor
!-----minimum temperature
      real*8  :: tempmin = 0.0000861  ! 1 Kelvin
!-----material properties. array of fixed size
      integer :: nmat=0            
      type(typematerial), dimension(100) :: mat
!-----laser and pulse
      type (typesimlaserwkb)      :: laser_wkb 
      type (typesimlaser3d)       :: laser_3d 
      type (typesimlasermaxwell)  :: laser_maxwell 
      type (typesimpulse)         :: pulse 
      end type typesimparam
!=======================================================================
!-----state vector (independent values)
      type typesimstate
      real*8 :: time                     ! time
      real*8, pointer :: x(:)  => null() ! positions of interfaces
      real*8, pointer :: v(:)  => null() ! velocities of interfaces
      real*8, pointer :: r(:)  => null() ! mass density
      real*8, pointer :: ee(:) => null() ! spec. electron internal energy
      real*8, pointer :: ei(:) => null() ! spec. ion internal energy
      real*8, pointer :: ea(:) => null() ! energy density of alpha part.
      real*8, pointer :: f(:)  => null() ! molar fraction of T (or D)
      real*8, pointer :: wa(:) => null() ! group frac. depos. (1+ng=heat) 
      end type typesimstate
!=======================================================================
!-----thermodynamic variables (can be computed from state vector)
      type typesimthermo
      real*8, pointer :: te(:)     => null() ! electron temperature
      real*8, pointer :: pe(:)     => null() ! electron pressure
      real*8, pointer :: dtedee(:) => null() ! derivative Dte/Dee
      real*8, pointer :: dtedr(:)  => null() ! derivative Dte/Dr
      real*8, pointer :: dpedee(:) => null() ! derivative Dpe/Dee
      real*8, pointer :: dpedr(:)  => null() ! derivative Dpe/Dr
      real*8, pointer :: ti(:)     => null() ! ion temperature
      real*8, pointer :: pi(:)     => null() ! ion pressure
      real*8, pointer :: dtidei(:) => null() ! derivative Dti/Dei
      real*8, pointer :: dtidr(:)  => null() ! derivative Dti/Dr
      real*8, pointer :: dpidei(:) => null() ! derivative Dpi/Dei
      real*8, pointer :: dpidr(:)  => null() ! derivative Dpi/Dr
      real*8, pointer :: zi(:)     => null() ! ion number (actual value)
      end type typesimthermo
!=======================================================================
!-----auxiliar quantities (temporal storage and output)
      type typesimenergy
!-----time step couters        
      integer :: istep         ! number of time steps
      integer :: itry          ! number of tried time steps
!-----laser deposition (power/mass): d(n)
      real*8, dimension(:), pointer :: d => null()
!-----radiation field: strahl(n+1,ng,2)
!-----                 strahl(i,j,1)=S+ at interface i for group j
!-----                 strahl(i,j,2)=S- at interface i for group j
      real*8, dimension(:,:,:), pointer :: strahl => null() 
!-----variables used for energy balance
      real*8 :: enie0          ! electron energy (initial)
      real*8 :: enii0          ! ion energy (initial)
      real*8 :: enia0          ! alpha particle energy (initial)
      real*8 :: enki0          ! kinetic energy (initial)
      real*8 :: enie           ! electron energy
      real*8 :: enii           ! ion energy
      real*8 :: enia           ! alpha particle energy
      real*8 :: enki           ! kinetic energy
      real*8 :: delas          ! absorbed laser energy
      real*8 :: derad          ! absorbed radiation
      real*8 :: defus          ! deposited fusion energy
      real*8 :: radoutl        ! Rad. output energy (left boundary)
      real*8 :: radoutr        ! Rad. output energy (right boundary)
      real*8 :: radinl         ! Rad. input energy (left boundary)
      real*8 :: radinr         ! Rad. input energy (right boundary)
      real*8 :: powoutl        ! Rad. output power (left boundary)
      real*8 :: powoutr        ! Rad. output power (right boundary)
      real*8 :: powinl         ! Rad. input power (left boundary)
      real*8 :: powinr         ! Rad. input power (right boundary)
      real*8 :: power          ! incident laser power
      real*8 :: alphas1        ! energy lost by alphas (last step)
      real*8 :: fusion1        ! fusion energy (last step)
      real*8 :: defus1         ! coupled fusion energy (last step)
      real*8 :: alphas         ! energy lost by alphas
      real*8 :: fusion         ! fusion energy 
      real*8 :: laser          ! incident laser energy
      real*8 :: neutrons       ! number of fusion neutrons
      real*8 :: ref            ! laser reflection fraction
      real*8 :: tra            ! laser transmision fraction
      end type typesimenergy
!=======================================================================
!-----definition of layered configuration
      type typesimlayers
      integer                 :: nl            ! number of layers
      character(len=80)       :: material(100) ! material names
      integer, dimension(100) :: nc            ! number of cells
      real*8,  dimension(100) :: thick         ! thickness
      real*8,  dimension(100) :: zonpar        ! griding parameter
      real*8,  dimension(100) :: r0,te0,ti0    ! initial values
      end type typesimlayers
!=======================================================================
      contains
!======================================================================!
!                                                                      !
!     control routines                                                 !
!                                                                      !
!======================================================================!
!-----main subroutine
      subroutine main
      implicit none
      type(typesimparam)  :: p
      type(typesimstate)  :: s
      type(typesimenergy) :: e 
      call init(p,s,e)
      call integrate(p,s,e)
      call finish(p,s,e)
      end subroutine main
!=======================================================================
!-----numerical integration
      subroutine integrate(p,s,e)
      implicit none
      type(typesimparam),  intent(in)    :: p
      type(typesimstate),  intent(inout) :: s
      type(typesimenergy), intent(inout) :: e
      real*8                             :: dt,taout,tbout
      character(len=32)                  :: filename
   
      !wfyuan, 2020/10/27, swtich off the output of fort.10,fort.11 and print_int_head
!-----print heading of run section
      !call print_int_head
!-----open files
      !open(10,file='fort.10',form='formatted')
      filename = "fort_"//trim(p%igafile)//".10"
      open(10,file=filename,form='formatted')
      !open(11,file='fort.11',form='formatted')
!-----initialize time step
      dt               = p%dtinit
!-----time of the next ascii and binary dumps
      taout            = 0
      tbout            = 0
!-----loop 
      do 
         e%istep       = e%istep+1
!-----integration step
         call step(p,s,e,dt) 
!-----unformatted output on unit 10, wfyuan mycall
         if(s%time>=tbout.or.mod(e%istep,p%ns_bout)==0)then
            tbout      = s%time+p%dt_bout
            call write_output_record(p,s,e) !wfyuan 2020/10/27 switch off for gaout

            !wfyuan, 2020/10/27, calculate fitness--rhoDT,Vimp
            call ga_calculate_fitness(p,s,e)
         end if
!-----formatted output on unit 11
         if(s%time>=taout.or.mod(e%istep,p%ns_aout)==0)then
            taout      = s%time+p%dt_aout
            !call ascii_output(p,s,e)
         end if
!-----exit loop
         if(s%time>=p%texit)exit
         if(e%istep>=p%nexit)exit
      end do
!-----close files
      close(10)
      !close(11)
!-----print tail of run section
      !call print_int_tail
      end subroutine integrate
!=======================================================================
!-----numerical integration
      subroutine step(p,s1,e,dt)
      implicit none
      type(typesimparam),  intent(in)    :: p
      type(typesimstate),  intent(inout) :: s1
      type(typesimenergy), intent(inout) :: e
      real*8,              intent(inout) :: dt
      real*8                             :: var
      integer                            :: icomp
      type(typesimstate)                 :: s2
!-----allocate memory
      call allocsimstate(p%n,p%ng,s2)
!-----loop for time step attempts 
      do  
         e%itry        = e%itry+1
!-----integration routine
         call schritt(p,s1,s2,dt,icomp,e,var)
!-----negative density, try again
         if(icomp<0.and.dt>p%dtmin)then
            call print_int_line(e%istep,s1%time,dt,var,-1)
            dt         = max(p%dtmin,dt/2)
!-----allowed maximum variation exceeded, try again
         else if(var>=p%dtbreak.and.dt>p%dtmin) then
            !call print_int_line(e%istep,s1%time,dt,var,-2)
            dt         = max(p%dtmin,dt/2)
!-----succesfull step
		 else
			! wfyuan,2020/10/27, switch off print for ga
            !call print_int_line(e%istep,s1%time,dt,var,1)
            exit
         end if
      end do
!-----replace state
      call copysimstate(s2,s1)
!-----update energy balance
      call update_energies(p,s1,dt,e)
!-----compute the next integration step
      if(var<=0)then
         dt            = dt*p%dtfactor
      else
         dt            = dt*min(p%dtfactor,1/var)
      end if
      dt               = max(p%dtmin,min(p%dtmax,dt))
!-----free memory
      call freesimstate(s2)
      end subroutine step
!=======================================================================
!-----perform one time step
      subroutine schritt(p,s1,s2,dt,icomp,e,var)
      implicit none 
      type(typesimparam),  intent(in)    :: p
      type(typesimstate),  intent(in)    :: s1
      type(typesimstate),  intent(inout) :: s2
      real*8,              intent(in)    :: dt
      integer,             intent(out)   :: icomp
      type(typesimenergy), intent(inout) :: e
      real*8,              intent(out)   :: var
!---- local variables
      real*8  :: wc,wb(p%ng+1)     ! for energy partition between groups
      real*8  :: d(p%n)            ! specific power laser deposition
      real*8  :: dae(p%n)          ! alpha->electron (power per volume unit)
      real*8  :: dai(p%n)          ! alpha->ion (power per unit of volume)
      real*8  :: qe(p%n)           ! specific power to electrons
      real*8  :: qi(p%n)           ! specific power to ions
      integer :: ictrl             ! cell with maximum deposition
      integer :: ictrl1(1)         ! auxiliar array needed by "maxloc"
      type(typesimstate)  :: s3    ! auxiliar states
      type(typesimthermo) :: t2,t3 ! auxiliar thermodynamic struct.
      integer :: n,ng              ! number of cells and groups
      integer :: isplit,kgroup     ! do loop counters
      real*8  :: dtf               ! sucycle time step
      real*8  :: aux               ! auxiliar storage
      real*8  :: varr,varte,varti  ! max variations in a subcycle
      real*8  :: timem             ! subcycle time 
!-----number  of cells and groups
      n        = p%n 
      ng       = p%ng
!-----sub-cycle time step
      dtf      = dt/p%msplit
!---- allocate local structures
      call allocsimstate(n,ng,s3)
      call allocsimthermo(n,t2)
      call allocsimthermo(n,t3)
!-----relative variation of variables / requested variation
      var      = 0
!-----powers: radiation fluxes, alpha particle losses, and fusion energy
!-----and laser deposition
      e%powoutr= 0
      e%powoutl= 0
      e%powinr = 0
      e%powinl = 0
      e%alphas1= 0
      e%fusion1= 0
      e%defus1 = 0
      e%d      = 0
!-----copy state (1-->2)
      call copysimstate(s1,s2)
!-----loop for subcycling 
      do isplit=1,p%msplit
!-----compute pressure, temp., ionization, and derivatives from tables
         call eosall(p,s2,t2)
!-----copy state (2-->3)
         call copysimstate(s2,s3)
         call copysimthermo(t2,t3)
!-----compute laser deposition
         timem = s1%time+(real(isplit)-0.5)/real(p%msplit)*dt
         call laser(timem,d,t3,s3,e,p)
!-----upgrade laser deposition on cells
         e%d           = e%d+d*s1%wa(p%ng+1)/p%msplit
!-----fusion 
         call fusion(p,s3,t3,e,dae,dai,dtf)
!-----electron heat transport (and e-i interchange)
         qi            = dai/s3%r
         qe            = dae/s3%r+d*s1%wa(p%ng+1)
         call heat(p,dtf,s3,t3,qe,qi)
!-----ion heat transport
         if(p%iheation>0) call heation(p,dtf,s3,t3)
!-----hydrodynamics
         if(p%ihydro>0) call hydro(p,dtf,s3,t3)
         if(minval(s3%r)<0)then
            icomp=-1
            go to 1
         end if
!-----variation of temperature at the cell with maximum deposition
         ictrl1        = maxloc(d)
         ictrl         = ictrl1(1)
         wc            = d(ictrl)*dt*t2%dtedee(ictrl)
         wb(p%ng+1)    = t3%te(ictrl)-t2%te(ictrl)
!-----radiation transport: loop for groups assigned to sub-cycle
         if(p%iradia==1)then
            do kgroup = isplit,p%ng,p%msplit
               if(kgroup<=p%ng)then
                  qe      = d*s1%wa(kgroup)
                  e%d     = e%d+d*s1%wa(kgroup)
                  aux     = t3%te(ictrl)
                  call group(p,e,dt,s3,t3,kgroup,qe)
                  wb(kgroup)=t3%te(ictrl)-aux
               end if
            end do  
         end if
!-----compute variation of variables
         call variation(n,s2%r,s3%r,0,varr)
         call variation(n,t2%te,t3%te,1,varte)
         call variation(n,t2%ti,t3%ti,1,varti)
         var=max(var,varr/p%dtrvar,varte/p%dttevar,varti/p%dttivar)
!-----copy state (3-->2)
         call copysimstate(s3,s2)
      end do  
!-----advance time
      s2%time  = s1%time+dt
!-----recompute the energy partition
      call wctrl(p%ng+1,s2%wa,wb,wc,p%iwctrl,p%iradia)
!-----set return indicator as OK
      icomp    = 1
!-----free memory
1     continue
      call freesimstate(s3)
      call freesimthermo(t2)
      call freesimthermo(t3)
      end subroutine schritt
!=======================================================================
!-----maximum relative variation of a vector
!-----if mode==1, use maximum of 3 neighbour elements as initial value
      subroutine variation(n,xi,xf,mode,vmax)
      implicit none
      integer,intent(in)  :: n,mode
      real*8, intent(in)  :: xi(n),xf(n)
      real*8, intent(out) :: vmax
      real*8              :: xr(n)
      integer             :: i
      if(mode==0)then
         vmax     = maxval(abs(xf-xi)/xi)
      else
         do i=2,n-1
            xr(i) = max(xi(i-1),xi(i),xi(i+1))
         end do
         xr(1)    = xr(2)
         xr(n)    = xr(n-1)
         vmax     = maxval(abs(xf-xi)/xr)
      end if
      end subroutine variation
!=======================================================================
!-----compute energy partition coefficients
      subroutine wctrl(ng,wa,wb,wc,iw,irad)
      implicit none
      integer, intent(in)    :: ng,iw,irad
      real*8,  intent(inout) :: wa(ng)
      real*8,  intent(in)    :: wb(ng),wc
      real*8                 :: waux
!-----if radiation has been switched out, all deposition to hydro
      if(irad==0)then
         wa    = 0
         wa(ng)= 1
!-----if deposition is zero, or if selected by the user, set equal parts
      else if(wc<=0.or.iw==0)then
         wa    = 1/real(ng)
!-----in standard situation, try to minimize peak temperature variations
      else
         waux  = sum(wb)/real(ng)
         wa    = max(0.d0,wa+(waux-wb)/wc)
         waux  = sum(wa)
         wa    = wa/waux
      end if
      end subroutine wctrl
!======================================================================!
!                                                                      !
!     D+T ignition and burning package with non-local deposition and   !
!     a simple model for alpha particle diffusion                      !
!                                                                      !
!======================================================================!
!-----advance fusion reactions (s%ea and s%f)
!-----depe and depi are the power to elect. ions per unit of volume
      subroutine fusion(p,s,t,e,depe,depi,dt)
      implicit none
      type(typesimparam), intent(in)      :: p
      type(typesimstate), intent(inout)   :: s
      type(typesimthermo),intent(in)      :: t
      type(typesimenergy),intent(inout)   :: e
      real*8,             intent(in)      :: dt
      real*8, dimension(p%n), intent(out) :: depe,depi
      real*8, dimension(p%n)   :: ea,f      ! initial values of ea and f
      real*8, dimension(p%n)   :: vol,den   ! cell volume and number density
      real*8, dimension(p%n+1) :: a         ! interface area
      real*8, dimension(p%n)   :: sigm      ! cross section
      real*8, dimension(p%n)   :: tauae,tau ! relaxation times
      real*8, dimension(p%n)   :: coulog    ! Coulomb logarithm
      real*8, dimension(p%n+1) :: q         ! diffusion coefficient
      real*8, dimension(3,p%n) :: as        ! tridiagonal matrix
      real*8, dimension(p%n)   :: aux       ! auxiliar storage 
      integer                  :: n,i
      real*8                   :: ffick,flimit,factor,egen,eloss
      real*8                   :: delta_ea,delta_ee,delta_ei
      real*8, dimension(p%n)   :: tk,xi,zeta
!-----physical constants
      real*8, parameter        :: eps0 = 5.632E-6
      real*8, parameter        :: v0   = 1.297E+9
!-----reactivity coefficients from S. Atzeni and J. Meyer-ter-Vehn
!-----'The Physics of Inertial Fusion' Table 1.3
      real*8, parameter        :: c0 =   6.6610d+00
      real*8, parameter        :: c1 =   643.41d-16
      real*8, parameter        :: c2 =   15.136d-03
      real*8, parameter        :: c3 =   75.189d-03
      real*8, parameter        :: c4 =   4.6064d-03
      real*8, parameter        :: c5 =   13.500d-03
      real*8, parameter        :: c6 = -0.10675d-03
      real*8, parameter        :: c7 =  0.01366d-03
!-----inmediate return
      if(p%nfuel==0)then
         depe  = 0
         depi  = 0
         return
      end if
!-----number of cells and initial values
      n        = p%n
      ea       = s%ea
      f        = s%f
!-----area of the interfaces and cell volumen
      select case(p%igeo)
      case(1)
         a     = 1
         vol   = s%x(2:n+1)-s%x(1:n)
      case(2)
         a     = 2.*cpi*s%x
         vol   = cpi*(s%x(2:n+1)**2-s%x(1:n)**2)
      case(3)
         a     = 4.*cpi*s%x**2
         vol   = 4.d0/3.d0*cpi*(s%x(2:n+1)**3-s%x(1:n)**3)
      end select
!-----compute electron number density, cross sections and relaxation times
      den      = s%r*t%zi/pm/p%ai
      tk       = t%ti/1000.
      xi       = c0/tk**(1.d0/3.d0)
      zeta     = 1-(c2*tk+c4*tk**2+c6*tk**3)/(1+c3*tk+c5*tk**2+c7*tk**3)
      sigm     = c1*xi**2/zeta**(5.d0/6.d0)*exp(-3*xi*zeta**(1.d0/3.d0))
      coulog   = max(2.d0,min(23.46-log(t%zi*sqrt(den)/ &
                 t%te**1.5),25.26-log(sqrt(den)/t%te)))
      tauae    = 3./16.*sqrt(2./cpi)*(4.*pm)*(t%te*bolz) &
                 **(3./2.)/(coulog*(2.*eq**2)**2*den*sqrt(em))
      tau      = tauae*33000./(33000.+t%te)
!-----compute diffusion coefficient for alpha particles
      aux      = v0**2*tauae/(9+1.251d-6*t%zi/p%ai*t%te**1.5)
      do i=2,n
         q(i)  = (aux(i-1)+aux(i))/(s%x(i+1)-s%x(i-1))
         ffick = abs(q(i)*(ea(i)-ea(i-1)))
         flimit= v0*max(ea(i),ea(i-1))/3.
         if(ffick<=0)then
            factor     = 1.0
         else
            factor     = min(1.d0,flimit/ffick)
         end if
         q(i)  = factor*q(i)*a(i)
      end do
      q(1)     = 0
      q(n+1)   = 0
!-----solve rate equation
      s%f      = f/(1+f*den*sigm*dt)
!-----solve diffusion equation for alfa particles
      as(1,2:n)= -q(2:n)/vol(1:n-1)
      as(2,1:n)= 1/dt+1/tau+(q(1:n-1)+q(2:n))/vol
      as(3,1:n-1)= -q(2:n)/vol(2:n)
      as(2,n)  = as(2,n)+v0*a(n+1)/2./vol(n)
      as(3,n-1)= as(3,n-1)-v0*a(n+1)/6./vol(n)
      s%ea     = (ea-eps0*den*(s%f-f))/dt
      call tridia(n,as,s%ea)
!-----energy deposition to electron and ions
      aux      = 33000./(33000.+t%te)
      !depe     = s%ea/tau*aux
      !depi     = s%ea/tau*(1.-aux)
      !wfyuan gafusion
      depe     = 0*s%ea/tau*aux
      depi     = 0*s%ea/tau*(1.-aux)
!-----generated power
      egen     = sum((f-s%f)*den*vol*eps0/dt)
      e%fusion1= e%fusion1+5*egen*dt
!-----power loss by flying alpha particles
      eloss    = (1.5*s%ea(n)-0.5*s%ea(n-1))*v0/3.*a(n)
      e%alphas1= e%alphas1+eloss*dt
!-----fusion energy coupled to the fluid
      delta_ea = sum((s%ea-ea)*p%dm/s%r)        !--on alphas
      delta_ee = sum(depe/s%r*p%dm)*dt          !--to electrons
      delta_ei = sum(depi/s%r*p%dm)*dt          !--to ions
      e%defus1 = e%defus1+delta_ea+delta_ee+delta_ei
      end subroutine fusion
!======================================================================!
!                                                                      !
!     electron heat transport and electron-ion interchange             !
!                                                                      !
!======================================================================!
!-----heat conduction sub-step
      subroutine heat(p,dt,s,t,de,di)
      implicit none
      real*8,                 intent(in)    :: dt
      type(typesimparam),     intent(in)    :: p
      type(typesimthermo),    intent(inout) :: t
      type(typesimstate),     intent(inout) :: s
      real*8, dimension(p%n), intent(in)    :: de,di
!-----local variables
      real*8, dimension(p%n)   :: ai,zi,te,ti,rho,dens,ee,ei
      real*8, dimension(p%n)   :: ecf_heat,ecf_qei,coefic,coefis
      real*8, dimension(p%n)   :: te32,kei,transf
      real*8, dimension(p%n+1) :: a,grt52,qsat,qclas,coeff,phi,cc
      real*8, dimension(p%n)   :: dte,dte_te_m,dte_te,dte_te_p,dte_ti
      real*8, dimension(p%n)   :: dti,dti_te,dti_ti
      real*8, dimension(p%n)   :: alpha,beta,bs
      real*8, dimension(3,p%n) :: as
      integer                  :: n
!-----copy variables
      n        = p%n
      ai       = p%ai
      rho      = s%r
      zi       = t%zi
      ee       = s%ee 
      ei       = s%ei 
      te       = t%te 
      ti       = t%ti
!-----area of interfaces
      select case (p%igeo)
      case (1)
         a     = 1
      case (2)
         a     = 2*cpi*s%x
      case (3)
         a     = 4*cpi*s%x**2
      end select
      dens     = zi*rho/(pm*ai)
      te32     = sqrt(abs(te))*te
!-----collision frequencies for ion-electron interchange
      call collfreq(p,n,dens,te,ti,zi,0.0d0,ecf_qei,3)
!-----heat conductivity
      if(p%flf>0)then
!-----collision frequencies for electron thermal flux interchange
         call collfreq(p,n,dens,te,ti,zi,0.0d0,ecf_heat,2)
!-----specific power transfer e-i: kei(i)*(te(i)-ti(i))
         kei      = 3.*(zi*bolz/ai/pm)/(ai*pm/em/ecf_qei)
!-----classical transport coefficient: q=-coefic*te**2.5*grad(te)
         coefic   = 3.22554*(zi+0.24)/(1.+0.24*zi)*dens*bolz**2/ &
                    (ecf_heat*em)/te32
!-----saturated transport coefficient: |q|=coefis*te**1.5*rho
         coefis   = (p%flf*bolz**1.5/pm/sqrt(em))*zi/ai
!-----temperature gradient times te**1.5
         grt52(2:n)= (te32(2:n)+te32(1:n-1))*(te(2:n)-te(1:n-1))/ &
                      (s%x(3:n+1)-s%x(1:n-1))
!-----saturated flux
         qsat(2:n)= (coefis(2:n)+coefis(1:n-1))*(te32(2:n)+te32(1:n-1))* &
                    (rho(2:n)+rho(1:n-1))/8
!-----classical flux
         qclas(2:n)= -0.5*(coefic(2:n)*te(2:n)+coefic(1:n-1)*te(1:n-1))* &
                     grt52(2:n)
!-----flux coefficient
         coeff(2:n)= -0.5*qsat(2:n)/(qsat(2:n)+abs(qclas(2:n)))* &
                   (coefic(2:n)*te(2:n)+coefic(1:n-1)*te(1:n-1))
!-----heat transfer: cc(i)*(te(i-1)-te(i)) for 1 < i < n+1
         cc(1)    = 0
         cc(2:n)  = -(te32(2:n)+te32(1:n-1))/(s%x(3:n+1)-s%x(1:n-1))* &
                     a(2:n)*coeff(2:n)
         cc(n+1)  = 0
      else
!-----no heat transfer
         cc       = 0
      end if
!-----power transfer e-i
      transf   = kei*(ti-te)
!-----heat flux
      phi(1)   = 0
      phi(2:n) = cc(2:n)*(te(1:n-1)-te(2:n))
      phi(n+1) = 0
!-----time derivatives
      dte      = (de+transf+(phi(1:n)-phi(2:n+1))/p%dm)*t%dtedee
      dti      = (di-transf)*t%dtidei
!-----derivatives of the derivatives
      dte_te_p = cc(2:n+1)*t%dtedee/p%dm;
      dte_te_m = cc(1:n)*t%dtedee/p%dm;
      dte_ti   = kei*t%dtedee
      dte_te   = -(dte_ti+dte_te_p+dte_te_m)
      dti_te   = kei*t%dtidei
      dti_ti   = -dti_te
!-----assemble system
      beta     = dt/(1-dt*dti_ti)
      alpha    = beta*dti_te
      beta     = beta*dti
      as(1,2:n)= dte_te_p(1:n-1)
      as(2,1:n)= dte_te-1/dt+dte_ti*alpha
      as(3,1:n-1) = dte_te_m(2:n)
      bs       = -dte-dte_ti*beta
!-----solve system
      call tridia(n,as,bs)
!-----advance variables
      t%te     = t%te+bs
      t%ti     = t%ti+alpha*bs+beta
      s%ee     = s%ee+(t%te-te)/t%dtedee
      s%ei     = s%ei+(t%ti-ti)/t%dtidei
      t%pe     = t%pe+t%dpedee*(s%ee-ee)
      t%pi     = t%pi+t%dpidei*(s%ei-ei)
!-----limit temperatures
      t%te     = max(t%te,p%tempmin)
      t%ti     = max(t%ti,p%tempmin)
      end subroutine heat
!======================================================================!
!                                                                      !
!     ion heat transport                                               !
!                                                                      !
!======================================================================!
!-----heat conduction sub-step
      subroutine heation(p,dt,s,t)
      implicit none
      real*8,                 intent(in)    :: dt
      type(typesimparam),     intent(in)    :: p
      type(typesimthermo),    intent(inout) :: t
      type(typesimstate),     intent(inout) :: s
!-----local variables
      real*8, dimension(p%n)   :: ti,ei,coef,bs
      real*8, dimension(p%n)   :: dti,dti_ti_m,dti_ti,dti_ti_p
      real*8, dimension(p%n+1) :: a,phi,cc
      real*8, dimension(3,p%n) :: as
      integer                  :: n,i
      real*8                   :: te,omp,vth,bxbn,clog,dens,delta,ai,zi
!-----copy variables
      n        = p%n
      ei       = s%ei 
      ti       = t%ti
!-----area of interfaces
      select case (p%igeo)
      case (1)
         a     = 1
      case (2)
         a     = 2*cpi*s%x
      case (3)
         a     = 4*cpi*s%x**2
      end select
!-----classical transport coefficient: q=-coefic*grad(te)
      do i=1,n
         ai    = p%ai(i)
         zi    = t%zi(i)
         te    = t%te(i)
         dens  = zi*s%r(i)/(pm*ai)
         omp   = sqrt(4*cpi*dens*eq*eq/em)
         vth   = min(sqrt(bolz*te*em)/hb,bolz*te/(zi*eq*eq))
         bxbn  = max(1.d0,2*vth*sqrt(te*bolz/em)/omp)
         clog  = log(1.+bxbn)
         delta = 0.095*(zi+0.24)/(1.+0.24*zi)
         coef(i) = delta*20.0*(2.0/cpi)**1.5*bolz* &
                   abs(bolz*ti(i))**2.5/ (clog*eq**4*sqrt(pm*ai)*zi**4)
      end do
!-----heat transfer: cc(i)*(ti(i-1)-ti(i)) for 1 < i < n+1
      cc(1)    = 0
      cc(2:n)  = 2.0*max(coef(2:n),coef(1:n-1))/ &
                 (s%x(3:n+1)-s%x(1:n-1))*a(2:n)
      cc(n+1)  = 0
!-----heat flux
      phi(1)   = 0
      phi(2:n) = cc(2:n)*(ti(1:n-1)-ti(2:n))
      phi(n+1) = 0
!-----time derivatives
      dti      = (phi(1:n)-phi(2:n+1))/p%dm*t%dtidei
!-----derivatives of the derivatives
      dti_ti_p = cc(2:n+1)*t%dtidei/p%dm;
      dti_ti_m = cc(1:n)*t%dtidei/p%dm;
      dti_ti   = -(dti_ti_p+dti_ti_m)
!-----assemble system
      as(1,2:n)= dti_ti_p(1:n-1)
      as(2,1:n)= dti_ti-1/dt
      as(3,1:n-1) = dti_ti_m(2:n)
      bs       = -dti
!-----solve system
      call tridia(n,as,bs)
!-----advance variables
      t%ti     = t%ti+bs
      s%ei     = s%ei+(t%ti-ti)/t%dtidei
      t%pi     = t%pi+t%dpidei*(s%ei-ei)
!-----limit temperatures
      t%ti     = max(t%ti,p%tempmin)
      end subroutine heation
!======================================================================!
!                                                                      !
!     implicit hydrodynamics                                           !
!                                                                      !
!======================================================================!
!-----hydrodynamic sub-step
      subroutine hydro(p,dt,s1,t1)
      implicit none
      type(typesimparam),  intent(in)    :: p
      real*8,              intent(in)    :: dt
      type(typesimthermo), intent(inout) :: t1
      type(typesimstate),  intent(inout) :: s1
!-----local variables
      real*8, dimension(p%n)     :: ee,ei,ea,r
      real*8, dimension(p%n)     :: vol,rho,pe,pi,pa,pt,c2e,c2i,c2a,c2
      real*8, dimension(p%n)     :: zeta,qa,qb,dv2,aux,f,g,dv
      real*8, dimension(p%n)     :: rcoef,tcoef,qcoef,scoef,aeff,beff
      real*8, dimension(p%n+1)   :: x,area,volume,b
      real*8, dimension(3,p%n+1) :: a
      integer                    :: n
!-----numerical parameters for ideal method
!     real*8, parameter        :: av1    = 2.0
!     real*8, parameter        :: av2    = 1.0
!     real*8, parameter        :: nu     = 1.0
!     real*8, parameter        :: auto   = 1.0
!     real*8, parameter        :: alpha  = 7.0/12.0
!     real*8, parameter        :: beta   = 1.0/12.0
!     real*8, parameter        :: lambda = 1.0/12.0
!-----numerical parameters for classical explicit method
!     real*8, parameter        :: av1    = 2.0
!     real*8, parameter        :: av2    = 0.0
!     real*8, parameter        :: nu     = 0.0
!     real*8, parameter        :: auto   = 0.0
!     real*8, parameter        :: alpha  = 0.0
!     real*8, parameter        :: beta   = 0.0
!     real*8, parameter        :: lambda = 0.0
!-----numerical parameters for robust implicit method
      real*8, parameter        :: av1    = 2.0
      real*8, parameter        :: av2    = 2.0
      real*8, parameter        :: nu     = 1.0
      real*8, parameter        :: auto   = 0.0
      real*8, parameter        :: alpha  = 1.0 
      real*8, parameter        :: beta   = 0.0
      real*8, parameter        :: lambda = 0.0
!-----copy variables
      n        = p%n
      r        = s1%r
      ee       = s1%ee
      ei       = s1%ei
      ea       = s1%ea
!-----intial pressures of electrons, ions, alphas (as ideal gas), and total
!-----(at initial time n-1/2)
      pe       = t1%pe
      pi       = t1%pi
      pa       = 2.0d0/3.0d0*s1%ea
      pt       = pe+pi+pa
!-----compute isentropic sound speed squared
!-----c**2=dP/dr, de=-Pd(1/r), dP=dPde*de+dPdr*dr
!-----c**2=dPdr+P/r**2*dPde 
!-----for monoatomic ideal gas c**2=5/3*p/r
!-----(assumed constant during the hydro step)
      c2e      = t1%dpedr+pe/s1%r**2*t1%dpedee
      c2i      = t1%dpidr+pi/s1%r**2*t1%dpidei
      c2a      = 5.0d0/3.0d0*pa/s1%r
      c2       = c2e+c2i+c2a
!-----advance node position (at intermediate time n)
      x        = s1%x+s1%v*(dt/2)
!-----areas and volumes of interfaces (at intermediate time n)
      select case(p%igeo)
      case(1)
         area  = 1
         volume= x
      case(2)
         area  = 2*cpi*x
         volume= area*x/2
      case(3)
         area  = 4*cpi*x**2
         volume= area*x/3
      end select
!-----cell volumes (at intermediate time n)
      vol      = volume(2:n+1)-volume(1:n)   
!-----mass density at t+dt/2 (at intermediate time n)
      rho      = p%dm/vol 
!-----time derivative of cell volumes (at time n-eps)
      dv       = area(2:n+1)*s1%v(2:n+1)-area(1:n)*s1%v(1:n)
!-----effective values of alpha and beta in each cell
      aeff     = alpha
      beff     = beta
      if(dt>0.and.auto>0)then
         where(c2>0)
            aux        = (x(2:n+1)-x(1:n))/(auto*dt) 
            aux        = aux**2/c2
            aux        = 1-(1-4*lambda)*aux
            aeff       = max((alpha+beta)*aux-beta,beta)
         end where 
      end if
!-----the effective thermodynamic pressure can be expressed as:
!-----pt-rcoef*c2-tcoef*c2*increment(dv)
!-----(at time n)
      rcoef    = rho*dt/vol*dv*(0.5+aeff-beff)
      tcoef    = rho*dt/vol*aeff
!-----effective viscous pressure can be expressed as:
!-----qcoef+scoef*increment(dv)
!-----time derivative of cell volumes per unit of area
      zeta     = 2*dv/(area(1:n)+area(2:n+1))
!-----viscous pressure qa : traditional expression
      qa       = av1*rho*zeta**2
!-----viscous pressure qb : on oscillating cells
      qb       = -av2*rho*sqrt(abs(c2))*zeta
      where(dv(1:n-2)<0.or.dv(3:n)<0.or.c2(2:n-1)<0)qb(2:n-1) = 0
      qb(1)    = 0
      qb(n)    = 0
!-----apply viscous pressure on collapsing cells
      where(dv<0)
         qcoef = max(qa,qb)
         scoef = nu*qcoef/dv
      elsewhere
         qcoef = 0
         scoef = 0
      end where
!-----effective total pressure times dt can be expressed as:
!-----f + g * increment(dv)
      f        = (pt-c2*rcoef+qcoef)*dt
      g        = (scoef-c2*tcoef)*dt
!-----assemble system of equations
      b(1)     = -f(1)*area(1)
      b(2:n)   = (f(1:n-1)-f(2:n))*area(2:n)
      b(n+1)   = f(n)*area(n+1)
      a(2,1)   = (0.5-lambda)*p%dm(1)-area(1)**2*g(1)
      a(2,2:n) = (0.5-lambda)*(p%dm(1:n-1)+p%dm(2:n)) &
                 -area(2:n)**2*(g(1:n-1)+g(2:n))
      a(2,n+1) = (0.5-lambda)*p%dm(n)-area(n+1)**2*g(n)
      a(1,2:n+1)= p%dm*lambda+area(1:n)*area(2:n+1)*g
      a(3,1:n) = a(1,2:n+1)
!-----apply boundary conditions
      if(p%ileft==1)then
          b(1) = 0
          a(2,1) = 1
          a(1,2) = 0
      end if
      if(p%iright==1)then
          b(n+1) = 0
          a(2,n+1) = 1
          a(3,n) = 0
      end if
!-----solve tridiagonal system: b contains increments of velocity
      call tridia(n+1,a,b)
!-----advance velocities and positions (at time n+1)
      s1%v = s1%v+b
      s1%x = x+s1%v*dt/2
!-----advance densities (at time n+1)
      select case (p%igeo)
      case (1)
         s1%r=p%dm/((s1%x(2:n+1)-s1%x(1:n)))
      case (2)
         s1%r=p%dm/(cpi*(s1%x(2:n+1)**2-s1%x(1:n)**2))
      case (3)
         s1%r=p%dm/(4.d0/3.d0*cpi*(s1%x(2:n+1)**3-s1%x(1:n)**3))
      end select
!-----increment(dv)
      dv2      = area(2:n+1)*b(2:n+1)-area(1:n)*b(1:n)
!-----effective pressures at pressures (at time n)
      aux      = rcoef+tcoef*dv2
      pe       = pe-c2e*aux
      pi       = pi-c2i*aux+qcoef+scoef*dv2
      pa       = pa-c2a*aux
!-----upgrade energies
      aux      = (dv+0.5*dv2)*dt/p%dm
      s1%ee    = ee - pe*aux
      s1%ei    = ei - pi*aux
      s1%ea    = (ea/r-pa*aux)*s1%r
!-----compute the new temperature and pressure (overwrite r, ee, and ei)
      r        = s1%r-r
      ee       = s1%ee-ee
      ei       = s1%ei-ei
      t1%te    = t1%te+t1%dtedr*r+t1%dtedee*ee
      t1%pe    = t1%pe+t1%dpedr*r+t1%dpedee*ee
      t1%ti    = t1%ti+t1%dtidr*r+t1%dtidei*ei
      t1%pi    = t1%pi+t1%dpidr*r+t1%dpidei*ei
!-----limit temperatures
      t1%te    = max(t1%te,p%tempmin)
      t1%ti    = max(t1%ti,p%tempmin)
      end subroutine hydro
!======================================================================!
!                                                                      !
!     radiation transport                                              !
!                                                                      !
!======================================================================!
!-----solve one group equation
      subroutine group(p,e,dt,ss,t,kgroup,qe)
      implicit none
      type(typesimparam),       intent(in)    :: p
      type(typesimenergy),      intent(inout) :: e
      type(typesimstate),       intent(inout) :: ss
      type(typesimthermo),      intent(inout) :: t
      integer,                  intent(in)    :: kgroup
      real*8,                   intent(in)    :: dt
      real*8, dimension(p%n),   intent(in)    :: qe
!-----Planck and Rosseland opacities, emissivity, Up and dUp/dT
      real*8, dimension(p%n)   :: op,or,up,upd,eps
!-----U (at cells), and increment of T
      real*8, dimension(p%n)   :: v,w
!-----U and S (at interfaces)
      real*8, dimension(p%n+1) :: u,s
!-----system of linear equations
      real*8                   :: z1,z2,z3,z4
      real*8, dimension(p%n+1) :: eta
      real*8, dimension(p%n)   :: b1,b2,a12,a21,a22
      real*8, dimension(3,p%n) :: a11
!-----other local variables
      integer                  :: n,igeo,i
      real*8                   :: alphal,alphar,betal,betar,fa,fb,f
      real*8                   :: areal,arear
      real*8                   :: tbath(1),ue(1),due(1),epse(1),flux
      n        = p%n
      igeo     = p%igeo
      alphal   = p%alphal
      alphar   = p%alphar
      betal    = p%betal
      betar    = p%betar
      fa       = p%frei(kgroup)
      fb       = p%frei(kgroup+1)
!-----obtain opacities and normalized emissivity at cells from tables
      call opacit(p,ss%r,t%te,(fa+fb)/2,op,or,eps)
!-----optionally compute an incident flux
      if(p%trad>0)then
         call pulse(p%itype,p%tau,p%trad,ss%time,tbath(1),p%pulse)
      else
         tbath = 0
      end if
      if(tbath(1)>0)then
         epse  = 1
         call planck(1,tbath,fa,fb,ue,due,epse)
         flux  = 0.25*c*ue(1)
      else
         flux  = 0
      end if
!-----compute planckian energy density 'up' and its derivative 'upd'
      call planck(n,t%te,fa,fb,up,upd,eps)
!-----solve the group equations
      call radia(n,igeo,ss%x,ss%r,op,or,up,upd,alphar,alphal, &
                 a11,a12,b1,eta,z1,z2,z3,z4,flux*betar,flux*betal)
      call mee(n,dt,t%dtedee,op,up,upd,qe,a21,a22,b2)
      call solve3(n,a11,a12,a21,a22,b1,b2,v,w)
!-----upgrade energy
      ss%ee    = ss%ee+w/t%dtedee
!-----limit electron temperature
      t%te     = max(t%te+w,p%tempmin)
!-----compute energy density
      u(1)     = z1*v(1)-z2*v(2)
      do i=2,n
         f     = (ss%x(i)-ss%x(i-1))/(ss%x(i+1)-ss%x(i-1))
         u(i)  = f*v(i)+(1-f)*v(i-1)
      end do
      u(n+1)   = -z3*v(n-1)+z4*v(n)
!-----compute the radiation flux
      s(1)     = -((1-alphal)*u(1)*c/2-2*flux*betal)/(1+alphal)
      s(2:n)   = eta(2:n)*(v(1:n-1)-v(2:n))
      s(n+1)   = ((1-alphar)*u(n+1)*c/2-2*flux*betar)/(1+alphar)
!-----compute S+ and S-
      e%strahl(:,kgroup,1) = 0.25*c*u+0.5*s
      e%strahl(:,kgroup,2) = 0.25*c*u-0.5*s
!-----energy balance
      select case(igeo)
      case(1)
         areal = 1
         arear = 1
      case(2)
         areal = 2*cpi*ss%x(1)
         arear = 2*cpi*ss%x(n+1)
      case(3)
         areal = 4*cpi*ss%x(1)**2
         arear = 4*cpi*ss%x(n+1)**2
      end select
      e%powinl         = e%powinl +areal*e%strahl(1,kgroup,1)
      e%powoutr        = e%powoutr+arear*e%strahl(n+1,kgroup,1)
      e%powoutl        = e%powoutl+areal*e%strahl(1,kgroup,2)
      e%powinr         = e%powinr +arear*e%strahl(n+1,kgroup,2)
      end subroutine group
!=======================================================================
!-----Planck function and its temperature derivative
      subroutine planck(n,t,fra,frb,up,upd,eps)
      implicit none
      integer,              intent(in)  :: n
      real*8,               intent(in)  :: fra,frb
      real*8, dimension(n), intent(in)  :: t,eps
      real*8, dimension(n), intent(out) :: up,upd
      real*8, dimension(n)              :: xa,xb,fa,fb,dfadx,dfbdx
      real*8, parameter                 :: utot=4*sigma/c
      xa       = fra/t
      xb       = frb/t
      call pdstrb(n,xa,fa,dfadx)
      call pdstrb(n,xb,fb,dfbdx)
      up       = eps*utot*t**4*(fb-fa)
      upd      = eps*utot*4*t**3*((fb-fa)-0.25*(xb*dfbdx-xa*dfadx))
      end subroutine planck
!=======================================================================
!-----compute the Planck function and its integral
!-----estimated precision 1/1000
      subroutine pdstrb(n,x,f,dfdx)
      implicit none
      integer,                intent(in)  :: n
      real*8,  dimension(n),  intent(in)  :: x
      real*8,  dimension(n),  intent(out) :: f,dfdx
      real*8,  dimension(n)               :: x1,x2,f1,f2
      integer, parameter                  :: nt=200
      real*8,  dimension(nt), save        :: xt,ft
      integer, save                       :: ivez=1
      real*8,  parameter                  :: dintgr=cpi**4/15.
      integer                             :: i,ind
!-----the first time create a table
      if(ivez==1)then
         xt(1)         = 1.e-12
         do i=2,nt
            xt(i)      = 0.03*real(i-1)+2.3e-8*real(i-1)**4.5
         end do
         ft(1)         = 0
         do i=2,nt
            ft(i)      = ft(i-1)+0.5/dintgr*(xt(i)**3/(exp(xt(i))-1)+ &
                         xt(i-1)**3/(exp(xt(i-1))-1.))*(xt(i)-xt(i-1))
         end do
         ft            = ft/ft(nt)
         ivez          = 0
      end if
!-----linear interpolation in table
      do i=1,n
         do ind=2,nt-1
            if(xt(ind)>x(i)) exit
         end do
         x1(i)         = xt(ind-1)
         x2(i)         = xt(ind)
         f1(i)         = ft(ind-1)
         f2(i)         = ft(ind)
      end do
      f                = f1+(f2-f1)/(x2-x1)*(x-x1)
      dfdx             = x**3/(dintgr*(exp(min(x,200.d0))-1))
      end subroutine pdstrb
!=======================================================================
!-----equation of radiation transfer
      subroutine radia(n,igeo,x,r,op,or,us,usd,alphar,alphal,au,aw,b, &
                       eta,z1,z2,z3,z4,fluxr,fluxl)
      implicit none
      integer,                intent(in)  :: n,igeo
      real*8,                 intent(in)  :: alphar,alphal,fluxr,fluxl
      real*8, dimension(n+1), intent(in)  :: x
      real*8, dimension(n),   intent(in)  :: r,op,or,us,usd
      real*8,                 intent(out) :: z1,z2,z3,z4
      real*8, dimension(n+1), intent(out) :: eta
      real*8, dimension(n),   intent(out) :: aw,b
      real*8, dimension(3,n), intent(out) :: au
      real*8                              :: theta(n+1),beta(n)
      real*8                              :: g,opamat,opageo,hl,hr
      integer                             :: i
!-----compute material and 'geometric' opacities at interfaces
      g        = 1./2.
      do i=2,n
         opamat=2/(1/(or(i-1)*r(i-1))+1/(or(i)*r(i)))
         opageo=((2*x(i)/(x(i)+x(i-1)))**(igeo-1)- &
                (2*x(i)/(x(i)+x(i+1)))**(igeo-1))/(x(i+1)-x(i-1))
         eta(i)=2*c*g**2/((x(i+1)-x(i-1))*(opamat+opageo))
      end do
!-----extrapolation of energy density at the boundaries
!-----U(left)  =  z1*U(1)  -z2*U(2)
!-----U(right) = -z3*U(n-1)+z4*U(n)
      z2       = (x(2)-x(1))/(x(3)-x(1))
      z3       = (x(n+1)-x(n))/(x(n+1)-x(n-1))
      z1       = 1+z2
      z4       = 1+z3
      hl       = c/2*(1-alphal)/(1+alphal)
      hr       = c/2*(1-alphar)/(1+alphar)
!-----compute intermediate quantities
      select case(igeo)
      case(1)
         theta = 1
         beta  = c*op*r*(x(2:n+1)-x(1:n))
      case(2)
         theta = 2*cpi*x
         beta  = 2*cpi*c*op*r*(x(2:n+1)**2-x(1:n)**2)/real(2)
      case(3)
         theta = 4*cpi*x**2
         beta  = 4*cpi*c*op*r*(x(2:n+1)**3-x(1:n)**3)/real(3)
      end select
!-----build banded system of equations
      au(2,1)          = theta(2)*eta(2)+theta(1)*z1*hl+beta(1)
      au(1,2)          = -theta(2)*eta(2)-theta(1)*z2*hl
      do i=2,n-1
         au(3,i-1)     = -theta(i)*eta(i)
         au(2,i)       = theta(i)*eta(i)+theta(i+1)*eta(i+1)+beta(i)
         au(1,i+1)     = -theta(i+1)*eta(i+1)
      end do
      au(3,n-1)        = -z3*hr*theta(n+1)-theta(n)*eta(n)
      au(2,n)          = z4*hr*theta(n+1)+theta(n)*eta(n)+beta(n)
      aw               = -beta*usd
      b                = beta*us
      b(1)             = b(1)+2/(1+alphal)*fluxl*theta(1)
      b(n)             = b(n)+2/(1+alphar)*fluxr*theta(n+1)
      end subroutine radia
!=======================================================================
!-----matter energy equation
      subroutine mee(n,dt,dtde,op,us,usd,q,au,aw,b)
      implicit none
      integer,              intent(in)  :: n
      real*8,               intent(in)  :: dt
      real*8, dimension(n), intent(in)  :: dtde,op,us,usd,q
      real*8, dimension(n), intent(out) :: au,aw,b
      real*8, dimension(n)              :: phi
      phi      = c*dt*op*dtde
      au       = -phi
      aw       = 1+phi*usd
      b        = -phi*us+dt*dtde*q
      end subroutine mee
!=======================================================================
!-----solve radiative transport equations
!-----    a11 * v + a12 * w = b1
!-----    a21 * v + a22 * w = b2
!-----v,w vectors, a12,a21,a22 diagonal matrices, a11 tridiagonal matrix
      subroutine solve3(n,a11,a12,a21,a22,b1,b2,v,w)
      implicit none
      integer,                intent(in)  :: n
      real*8, dimension(3,n), intent(in)  :: a11
      real*8, dimension(n),   intent(in)  :: a12,a21,a22,b1,b2
      real*8, dimension(n),   intent(out) :: v,w
      real*8                              :: a(3,n),alpha(n)
      a(1,:)   = a11(1,:)
      a(3,:)   = a11(3,:)
      alpha    = a12/a22
      a(2,:)   = a11(2,:)-alpha*a21
      v        = b1-alpha*b2
      call tridia(n,a,v)
      w       = (b2-a21*v)/a22
      end subroutine solve3
!======================================================================!
!                                                                      !
!     Upgrade energy balance                                           !
!                                                                      !
!======================================================================!
!-----if dt<0,  the structure 'e' is created and initialized
!-----if dt>=0  update energy balance structure 'e'
      subroutine update_energies(p,s,dt,e)
      implicit none
      type(typesimparam),  intent(in)    :: p
      type(typesimstate),  intent(in)    :: s
      real*8,              intent(in)    :: dt
      type(typesimenergy), intent(inout) :: e
      real*8, dimension(p%n+1)           :: dmi
      integer                            :: n,ng,nfuel
      n        = p%n
      ng       = p%ng
      nfuel    = p%nfuel
      if(dt<0)then
         allocate(e%d(n),e%strahl(n+1,ng,2))
         memory = memory+2
      end if
      dmi(1)   = p%dm(1)/2
      dmi(2:n) = (p%dm(2:n)+p%dm(1:n-1))/2
      dmi(n+1) = p%dm(n)/2
      e%enii   = sum(s%ei*p%dm)
      e%enie   = sum(s%ee*p%dm)
      e%enki   = sum(0.5*dmi*s%v**2)
      e%enia   = sum(p%dm/s%r*s%ea)
      if(dt<0)then
         e%d           = 0
         e%strahl      = 0
         e%tra         = 0
         e%ref         = 0
         e%radoutl     = 0
         e%radoutr     = 0
         e%radinl      = 0
         e%radinr      = 0
         e%alphas      = 0
         e%fusion      = 0
         e%neutrons    = 0
         e%laser       = 0
         e%defus       = 0
         e%derad       = 0
         e%delas       = 0
         e%enii0       = e%enii
         e%enie0       = e%enie
         e%enia0       = e%enia
         e%enki0       = e%enki
         e%istep       = 0
         e%itry        = 0
      else
         e%radoutl     = e%radoutl+e%powoutl*dt
         e%radoutr     = e%radoutr+e%powoutr*dt
         e%radinl      = e%radinl +e%powinl*dt
         e%radinr      = e%radinr +e%powinr*dt
         e%alphas      = e%alphas+e%alphas1
         e%fusion      = e%fusion+e%fusion1
         e%neutrons    = sum((0.5-s%f(1:nfuel))*p%dm(1:nfuel)/ &
                         (p%ai(1:nfuel)*pm))
         e%laser       = e%laser+e%power*dt
         e%defus       = e%defus+e%defus1
         e%derad       = e%radinl+e%radinr-e%radoutl-e%radoutr
         e%delas       = e%delas+sum(e%d*p%dm*dt)
      end if
      end subroutine update_energies
!======================================================================!
!                                                                      !
!     laser deposition (and pulse definition)                          !
!                                                                      !
!======================================================================!
      subroutine laser(time,d,t,s,e,p)
      implicit none
      real*8,              intent(in)     :: time
      type(typesimparam),  intent(in)     :: p
      type(typesimthermo), intent(in)     :: t
      type(typesimstate),  intent(in)     :: s
      type(typesimenergy), intent(inout)  :: e
      real*8,              intent(out)    :: d(p%n)
      integer                             :: ic
      ic       = 0
      if(p%laser_maxwell%enable==1)then
         if(p%igeo /= 1)then
             print *,'ERROR in subroutine laser'
             print *,'      Maxwell solver only available for igeo=1'
             stop
         end if
         call laser_maxwell(time,d,t,s,e,p,p%laser_maxwell,p%pulse)
         ic    = ic+1
      end if
      if(p%laser_wkb%enable==1)then
         call laser_wkb(time,d,t,s,e,p,p%laser_wkb,p%pulse)
         ic    = ic+1
      end if
      if(p%laser_3d%enable==1)then
         if(p%igeo /= 3)then
             print *,'ERROR in subroutine laser'
             print *,'      3D solver only available for igeo=3'
             stop
         end if
         call laser_3d(time,d,t,s,e,p,p%laser_3d,p%pulse)
         ic    = ic+1
      end if
      if(ic==0)then
         e%power= 0
         d     = 0
      end if
      if(ic>1) then
         print *,'ERROR in subroutine laser'
         print *,'      more than one active laser'
         stop
      end if
      end subroutine laser
!=======================================================================
! laser deposition routine using the WKB approximation
      subroutine laser_wkb(time,extern,t,st,e,p,laserwkb,puls)
      implicit none
      real*8,                intent(in)    :: time
      type(typesimparam),    intent(in)    :: p
      type(typesimthermo),   intent(in)    :: t
      type(typesimstate),    intent(in)    :: st
      type(typesimlaserwkb), intent(in)    :: laserwkb
      type(typesimpulse),    intent(in)    :: puls
      type(typesimenergy),   intent(inout) :: e
      real*8,                intent(out)   :: extern(p%n)
      real*8, dimension(p%n+1) :: s
      real*8, dimension(p%n)   :: a,dene,ecf
      integer                  :: n,itype,inter
      real*8                   :: pimax,pitime,wl,delta,freq,omega
      real*8                   :: dncrt
!-----copy variables
      n        = p%n
      pimax    = laserwkb%pimax
      pitime   = laserwkb%pitime
      wl       = laserwkb%wl
      itype    = laserwkb%itype
      inter    = laserwkb%inter
      delta    = laserwkb%delta
!-----compute critical density
      freq     = c/wl
      omega    = 2*cpi*freq
      dncrt    = cpi*freq**2*em/eq**2
!-----compute incident power
      call pulse(itype,pitime,pimax,time,e%power,puls)
!-----if incident power is zero, return inmediatelly
      if(e%power<=0)then
         extern= 0
         return
      end if
!-----compute electron number density and relative density
      dene     = (st%r*t%zi)/(p%ai*pm)
      s(1)     = 0
      s(2:n)   = (dene(2:n)+dene(1:n-1))/(2*dncrt)
      s(n+1)   = 0
!-----compute electron collision frequency
      call collfreq(p,n,dene,t%te,t%ti,t%zi,omega,ecf,1)
!-----compute absortion coefficient
      a=(dncrt/c)*(ecf/dene)*(st%x(2:n+1)-st%x(1:n))
      if(inter>0)then !-----laser comes from the right
         call atenua(n,s,a,delta,extern,e%tra,e%ref)
      else            !-----laser comes from the left
         s     = s(n+1:1:-1)
         a     = a(n:1:-1)
         call atenua(n,s,a,delta,extern,e%tra,e%ref)
         extern= extern(n:1:-1)
      end if
!-----deposition per unit of mass
      extern   = extern*e%power/p%dm
      end subroutine laser_wkb
!=======================================================================
!-----compute deposition profile for normaliced incidence
      subroutine atenua(n,s,a,delta,dep,tra,ref)
      implicit none
      integer, intent(in)  :: n
      real*8,  intent(in)  :: s(n+1),a(n),delta
      real*8,  intent(out) :: dep(n),tra,ref
      real*8               :: gamma(n+1),phi(n+1),depth(n)
      real*8               :: aux,factor,phicrt
      integer              :: i,icrt
!-----find critical cell
      do i=n,1,-1
         if(s(i)>=1)exit
      end do
      icrt     = i
!-----compute optical depth of the underdense cells
      do i=icrt+1,n+1
         gamma(i)=-sqrt(1-s(i))*(6*s(i)**2+8*s(i)+16)/15
      end do
      if(icrt<n)then
         do i=icrt+1,n
            if(s(i)==s(i+1))then
               depth(i) = a(i)*s(i)**2/sqrt(1-s(i))
            else
               depth(i) = a(i)*(gamma(i+1)-gamma(i))/(s(i+1)-s(i))
            end if
         end do
      end if
!-----compute optical depth of critical cell, if one exists
      if(icrt/=0)depth(icrt) = a(icrt)*gamma(icrt+1)/(s(icrt+1)-1)
!-----compute the incident flux
      phi(n+1) = 1
      do i=n,1,-1
         if(i>icrt)then
            phi(i) = phi(i+1)*exp(-min(100.d0,depth(i)))
         else
            phi(i) = 0
         end if
      end do
!-----incident flux on the critical surface and total flux
      if(icrt/=0)then
         phicrt = phi(icrt+1)*exp(-min(100.d0,depth(icrt)))
         if(phicrt>0)then
            do i=icrt+1,n+1
               phi(i) = phi(i)-phicrt**2*(1-delta)/phi(i)
            end do
         end if
      end if
!-----energy deposition
      do i=1,n
         dep(i) = phi(i+1)-phi(i)
      end do
!-----smoothing at the critical surface
      if(icrt>1.and.icrt<n)then
         factor      = (s(icrt)-1)/(s(icrt)-s(icrt+1))
         aux         = dep(icrt)
         dep(icrt-1) = (1-factor)*0.5*aux
         dep(icrt)   = 0.5*aux
         dep(icrt+1) = dep(icrt+1)+factor*0.5*aux
      end if
!-----transmission and reflection
      tra            = phi(1)
      if(icrt/=0)then
         ref         = phicrt**2*(1-delta)
      else
         ref         = 0
      end if
      end subroutine atenua
!=======================================================================
!-----laser deposition routine solving the 1D Maxwell's equations
      subroutine laser_maxwell(time,extern,t,st,e,p,lasermaxwell,puls)
      implicit none
      real*8,                    intent(in)    :: time
      type(typesimparam),        intent(in)    :: p
      type(typesimstate),        intent(in)    :: st
      type(typesimthermo),       intent(in)    :: t
      type(typesimlasermaxwell), intent(in)    :: lasermaxwell
      type(typesimpulse),        intent(in)    :: puls
      type(typesimenergy),       intent(inout) :: e
      real*8,                    intent(out)   :: extern(p%n)
      real*8, dimension(p%n+1)                 :: g,h,z,dene
      real*8, dimension(p%n*lasermaxwell%idep) :: as,dx,ro,dens,temp
      real*8, dimension(p%n*lasermaxwell%idep) :: zet,tion,dep,ecf
      character(len=1) :: pol
      integer          :: n,inter,itype,idep,nm,l,lk,i
      real*8           :: pimax,pitime,wl,angle,freq,omega,wvec,dncrt
      real*8           :: fk,hkl
      n                = p%n
      inter            = lasermaxwell%inter
      pimax            = lasermaxwell%pimax
      pitime           = lasermaxwell%pitime
      wl               = lasermaxwell%wl
      itype            = lasermaxwell%itype
      idep             = lasermaxwell%idep
      angle            = lasermaxwell%angle
      if(lasermaxwell%pol==1)pol='P'
      if(lasermaxwell%pol==2)pol='S'
      pimax            = pimax*cos(angle/180*cpi)
!-----compute critical density
      freq             = c/wl
      omega            = 2*cpi*freq
      wvec             = 2*cpi/wl
      dncrt            = cpi*freq**2*em/eq**2
!-----compute incident power
      call pulse(itype,pitime,pimax,time,e%power,puls)
!-----if incident power is zero, return inmediatelly
      if(e%power<=0)then
         extern        = 0
         return
      end if
!-----compute interface centered temperatures 'h' and 'g', ionization 'z',
!-----and electron number density 'dene'
      h(1)             = t%te(1)
      g(1)             = t%ti(1)
      z(1)             = t%zi(1)
      dene(1)          = 0
      do i=2,n
         h(i)          = (t%te(i)+t%te(i-1))/2
         g(i)          = (t%ti(i)+t%ti(i-1))/2
         z(i)          = (t%zi(i)+t%zi(i-1))/2
         dene(i)       = 2/((pm*p%ai(i))/(st%r(i)*t%zi(i))+ &
                         (pm*p%ai(i-1))/(st%r(i-1)*t%zi(i-1)))
      end do  
      h(n+1)           = t%te(n)
      g(n+1)           = t%ti(n)
      z(n+1)           = t%zi(n)
      dene(n+1)        = 0
!-----subdivide each cell in 'idep' subcells. for each subcell compute:
!     dx:   cell size times wavevector
!     temp: electron temperature
!     tion: ion temperature
!     dens: electron number density
!     ro:   electron number density over critical density
!     zet:  ionization 
      nm              = idep*n
      fk              = 1/real(idep)
      do l=1,nm
         lk           = (l-1)/idep + 1
         hkl          = mod(l-1,idep)+0.5
         dx(l)        = fk*wvec*(st%x(lk+1)-st%X(lk))
         dens(l)      = dene(lk)+fk*hkl*(dene(lk+1)-dene(lk))
         ro(l)        = dens(l)/dncrt
         temp(l)      = h(lk)+fk*hkl*(h(lk+1)-h(lk))
         tion(l)      = g(lk)+fk*hkl*(g(lk+1)-g(lk))
         zet(l)       = z(lk)+fk*hkl*(z(lk+1)-z(lk))
      end do
!-----compute electron collision frequency
      call collfreq(p,nm,dens,temp,tion,zet,omega,ecf,1)
!-----compute normalized collision frequency
      as              = ecf/omega
!-----laser comes from the right hand side
      if(inter>0)then
         dx           = dx(nm:1:-1)
         ro           = ro(nm:1:-1)
         as           = as(nm:1:-1)
         call depos(nm,angle,ro,dx,as,pol,dep,e%tra,e%ref)
         dep          = dep(nm:1:-1)
!-----laser comes from the left hand side
      else
         call depos(nm,angle,ro,dx,as,pol,dep,e%tra,e%ref)
      end if
!-----compute cell deposition from subdivided Lagrangean cells
      do i=1,n
         extern(i)    = sum(dep((i-1)*idep+1:i*idep))
      end do
!-----deposition per unit of mass
      extern          = extern*e%power/p%dm
      end subroutine laser_maxwell
!=======================================================================
!-----laser depostion using planar Maxwell's equations
!-----compute deposition profile for normalized incidence laser 
!-----deposition over a profile composed of successive density steps 
!-----(Born/Wolf). Routine originally written by K.Eidmann and S. Hueller 
!-----Aug. 1991/Feb. 1992. Modified by R. Ramis Jul. 2009/ Jul. 2010/ 
!-----Jan. 2014.
      subroutine depos(nm,angle,ro,dx,ab,pol,dep,tra,ref)
      implicit none
      integer,               intent(in)  :: nm
      real*8,                intent(in)  :: angle
      real*8, dimension(nm), intent(in)  :: ro,dx,ab
      character(len=1),      intent(in)  :: pol
      real*8, dimension(nm), intent(out) :: dep
      real*8,                intent(out) :: tra,ref
      real*8  :: alf,pp0,trans,flux(nm+1)
      integer :: l,nmb
      complex :: uright,b11,b12,b21,b22,ep,br,ar,si,co,pp
      complex :: a(nm,2,2),u(nm+1),v(nm+1)
!-----sin of angle of incidence
      alf              = sin(cpi*angle/180)**2
!-----complex refractive index imaginary part of ep (conductivity)
      pp0              = sqrt(1-alf)
!-----characteristic matrix of the l-th layer
      do l=1,nm
         ep            = 1-ro(l)/cmplx(1,ab(l))
         br            = sqrt(ep-alf)
         if(pol=='P'.or.pol=='p')then
            pp         = br/ep
         else
            pp         = br
         end if
         ar            = br*dx(l)
         co            = cos(ar)
         si            = (0,1)*sin(ar)
         b11           = co
         b12           =-si/pp
         b21           =-si*pp
         b22           = co
         if(l==1)then
            a(1,1,1)   = b11
            a(1,1,2)   = b12
            a(1,2,1)   = b21
            a(1,2,2)   = b22
         else
            a(l,1,1)   = a(l-1,1,1)*b11+a(l-1,1,2)*b21
            a(l,1,2)   = a(l-1,1,1)*b12+a(l-1,1,2)*b22
            a(l,2,1)   = a(l-1,2,1)*b11+a(l-1,2,2)*b21
            a(l,2,2)   = a(l-1,2,1)*b12+a(l-1,2,2)*b22
         end if
!-----monitor transmision and cut off computations if it becomes too small
         trans         = abs(2*pp0/(a(l,2,1)+(a(l,1,1)+a(l,2,2))*pp0+ &
                         a(l,1,2)*pp0*pp0))**2
         if(trans<1e-6)then
            nmb        = l
            go to 1
         end if
      end do  
      nmb=nm
1     continue
!-----electric and magnetic field of the laser 
!-----uright is u(n+1)/sqrt(8*pi*phi(inc)/c*cos(theta))
      uright           = 2*pp0/(a(nmb,2,1)+(a(nmb,1,1)+a(nmb,2,2)+ &
                         a(nmb,1,2)*pp0)*pp0)
!-----for p-polarization u=b_x, for s-polarization u=e_y (divided by
!-----sqrt(8*pi*phi(inc)/c*cos(theta))
      u(1)             = uright*(a(nmb,1,1)+a(nmb,1,2)*pp0)
      v(1)             = uright*(a(nmb,2,1)+a(nmb,2,2)*pp0)
      do l=1,nmb
         u(l+1)        = a(l,2,2)*u(1)-a(l,1,2)*v(1)
         v(l+1)        =-a(l,2,1)*u(1)+a(l,1,1)*v(1)
      end do
      if(nmb<nm)then
         do l=nmb+1,nm
            u(l+1)     = u(nmb+1)
            v(l+1)     = v(nmb+1)
         end do
      end if
!-----energy flux (divided by phi(inc))
      flux             = (real(u)*real(v)+aimag(u)*aimag(v))/pp0
      do l=1,nm
         dep(l)        = flux(l)-flux(l+1)
      end do
!-----reflection and transmission
      ref              = 1-flux(1)
      tra              = flux(nm+1)
      end subroutine depos
!=======================================================================
!-----laser deposition routine solving 3D ray tracing         
      subroutine laser_3d(time,extern,t,st,e,p,laser3d,puls)
      implicit none
      real*8,               intent(in)    :: time
      type(typesimparam),   intent(in)    :: p
      type(typesimthermo),  intent(in)    :: t
      type(typesimstate),   intent(in)    :: st
      type(typesimlaser3d), intent(in)    :: laser3d
      type(typesimpulse),   intent(in)    :: puls
      type(typesimenergy),  intent(inout) :: e
      real*8,               intent(out)   :: extern(p%n)
      real*8, dimension(p%n)   :: dene,ecf,dray
      real*8, dimension(p%n+1) :: s
      integer                  :: n,itype,nr,i
      real*8                   :: pimax,pitime,wl,rmax,fwhm,bexp
      real*8                   :: freq,omega,dncrt,wt,pi,w
!-----copy variables
      n        = p%n
      pimax    = laser3d%pimax
      pitime   = laser3d%pitime
      wl       = laser3d%wl
      itype    = laser3d%itype
      nr       = laser3d%nr
      rmax     = laser3d%rmax
      fwhm     = laser3d%fwhm
      bexp     = laser3d%bexp 
!-----compute critical density
      freq     = c/wl
      omega    = 2*cpi*freq
      dncrt    = cpi*freq**2*em/eq**2
!-----compute incident power
      call pulse(itype,pitime,pimax,time,e%power,puls)
!-----if incident power is zero, return inmediatelly
      if(e%power<=0)then
         extern=0
         return
      end if
!-----compute electron number density
      dene     = (st%r*t%zi)/(p%ai*pm)
!-----number density over critical density
      s(1)     = 0
      s(2:n)   = (dene(2:n)+dene(1:n-1))/(2*dncrt)
      s(n+1)   = 0
!-----compute electron collision frequency 'ecf'
      call collfreq(p,n,dene,t%te,t%ti,t%zi,omega,ecf,1)
!-----add deposition for all impact parameters
      wt=0
      extern=0
      do i=1,nr
!-----impact parameter is 'pi' and has weight 'w'
         pi    = real(i)/real(nr)*rmax
         w     = pi*exp(-(2.*pi/fwhm)**bexp*log(2.0))
         wt    = wt+w
!-----compute normalized deposition
         call laser_3d_ray(n,st%x,s,pi,ecf,dray)
         extern= extern+dray*w
      end do
!-----normalize
      extern   = extern/wt
!-----compute reflection
      e%ref = 1
      do i=1,n
         e%ref = e%ref-extern(i)
      end do
!-----deposition per unit of mass
      extern   = extern*e%power/p%dm
      end subroutine laser_3d
!=======================================================================
!-----compute 3d laser deposition for an individual ray
      subroutine laser_3d_ray(n,r,s,p,ecf,dray)
      implicit none
      integer, intent(in)  :: n
      real*8,  intent(in)  :: r(n+1),s(n+1),ecf(n),p
      real*8,  intent(out) :: dray(n)
      real*8               :: p2,x1,x2,g1,g2,a,b,f,fluxcr,flux(n+1)
      integer              :: i,icrt
!-----check intersection
      if(p>=r(n+1))then
         dray=0
         return
      end if
!-----compute incident flux
      flux(n+1) = 1
      p2        = p*p
      do i=n,1,-1
!-----the optical depth of a cell is given by
!-----ecf/c*integral((1-n**2)/sqrt(n**2-p**2/x**2),dx)
         x1     = r(i)**2
         x2     = r(i+1)**2
!-----index of refraction squared = n**2 = 1-s = a*x+b
         a      = (s(i)-s(i+1))/(x2-x1)
         b      = 1-s(i+1)-a*x2
!-----the opticall depth of a cell is given by (c=-p**2)
!-----ecf/(2*c)*integral((1-b-a*x)/sqrt(a*x**2+b*x+c),dx)
         f=ecf(i)/(2*c)
         if(x1*(1-s(i))>p2)then
             call intlaser(a,b,-p2,x1,g1)
             call intlaser(a,b,-p2,x2,g2)
             flux(i)=flux(i+1)*exp(f*(g1-g2))
         else
             call intlaser(a,b,-p2,x2,g2)
             fluxcr=flux(i+1)*exp(-f*g2)
             exit
         end if
      end do
      icrt     = i
!-----take into acount outcoming flux (avoid division by zero)
      if(icrt>0)then
         do i=1,icrt
            flux(i)=0
         end do
         if(fluxcr>1.0e-36)then
            do i=icrt+1,n+1
               flux(i)=flux(i)-fluxcr*fluxcr/flux(i)
            end do
         end if
      end if
!-----compute normalized deposition
      dray=flux(2:n+1)-flux(1:n)
      end subroutine laser_3d_ray
!=======================================================================
!-----comp. integral of (1-b-a*x)/sqrt(a*x**2+b*x+c) between xmin and x
!-----xmin is the smaller one of the positive roots of a*x**2+b*x+c
      subroutine intlaser(a,b,c,x,g)
      implicit none
      real*8, intent(in)   :: a,b,c,x
      real*8, intent(out)  :: g
      real*8               :: d,u
!-----compute first the integral of 1/sqrt(a*x**2+b*x+c)
      if(a==0)then
         g=(2/b)*sqrt(abs(b*x+c))
      else
         d=sqrt(abs(b*b-4.*a*c))
         u=(2*a*x+b)/d
         if(a>0)then
            g=log(sqrt(abs(u*u-1))+u)/sqrt(a)
         else
            g=acos(u)/sqrt(-a)
         end if
      end if
!-----now compute full integral
      g=(1-b/2)*g-sqrt(abs(a*x*x+b*x+c))
      end subroutine intlaser
!=======================================================================
!-----computes power or radiation temperature as a function of time
      subroutine pulse(itype,pitime,pimax,time,power,puls)
      implicit none
      integer,            intent(in)  :: itype
      real*8,             intent(in)  :: pitime,pimax,time
      real*8,             intent(out) :: power
      type(typesimpulse), intent(in)  :: puls
      real*8, dimension(puls%ntab)    :: ttab,ptab
      real*8                          :: factor,p1,p2
      integer                         :: i
      select case(itype)
      case(1)
         if(time<2*pitime)then
            power=pimax*(sin(cpi/2*(time/pitime)))**2
         else
            power=0
         end if
      case(2)
         if(time<0)then
            power=0
         else if(time<=pitime)then
            power=pimax
         else
            power=0
         end if
      case(4)
         if(puls%ntab==0) then
            print *,'ERROR in subroutine pulse'
            print *,'      table not available'
            stop
         end if
         ttab=puls%ttab(1:puls%ntab)*pitime
         ptab=puls%ptab(1:puls%ntab)*pimax
         do i=1,puls%ntab
            if (ttab(i)>time) exit
         end do
         if(i==1)then
            power=ptab(1)
         else if(i<=puls%ntab)then
            factor     = (time-ttab(i-1))/(ttab(i)-ttab(i-1))
            p1         = ptab(i-1)
            p2         = ptab(i)
            select case (puls%mode)
            case(1)
               p1      = p1**puls%expo
               p2      = p2**puls%expo
               power   = p1+(p2-p1)*factor
               power   = power**(1/puls%expo)
            case(2)
               p1      = log(p1)
               p2      = log(p2)
               power   = p1+(p2-p1)*factor
               power   = exp(power)
            case(3)
               p1      = exp(p1)
               p2      = exp(p2)
               power   = p1+(p2-p1)*factor
               power   = log(power)
            end select
         else
            power=ptab(puls%ntab)
         end if
      case default
         print *,'ERROR in subroutine pulse'
         print *,'      unknown pulse type'
         stop
      end select
      end subroutine pulse
!======================================================================!
!                                                                      !
!     models for electron collision frequency                          !
!                                                                      !
!======================================================================!
!-----compute collision frequency
      subroutine collfreq(p,n,dens,tele,tion,zet,omega,ecf,mode)
      implicit none
      type(typesimparam),   intent(in)  :: p
      integer,              intent(in)  :: n,mode
      real*8,               intent(in)  :: omega
      real*8, dimension(n), intent(in)  :: dens,tele,tion,zet
      real*8, dimension(n), intent(out) :: ecf
      real*8                            :: vorf
      select case(p%model)
      case(0) !-----classical plasma
         call ecf_class(n,dens,tele,zet,omega,ecf)
      case(1) !-----electron phonon interaction
         if(mode==1)vorf=p%fheat
         if(mode==2)vorf=p%fei
         if(mode==3)vorf=p%flaser
         call ecf_kse(n,dens,tele,tion,zet,omega,ecf,vorf)
      case(2) !-----Drude-Sommerfeld model 
         call stossfrequenz(n,dens,tele,zet,omega,ecf)
      case default
         print *,'ERROR in subroutine collfreq'
         print *,'      bad value of model'
         stop
      end select
      end subroutine collfreq
!=======================================================================
!-----compute collision frequency for a classical plasma
      subroutine ecf_class(n,dens,tele,zet,omega,ecf)
      implicit none
      integer,              intent(in)  :: n
      real*8,               intent(in)  :: omega
      real*8, dimension(n), intent(in)  :: dens,tele,zet
      real*8, dimension(n), intent(out) :: ecf
      real*8                            :: coef,omp,vth,bxbn,clog
      integer                           :: i
!-----coefficient for classical collision frequency
      coef= sqrt(2*cpi)*4./3. *eq**4 *em/(em*bolz)**1.5
      do i=1,n
!-----plasma frequency
         omp =sqrt(4*cpi*dens(i)*eq*eq/em)
!-----thermal wavenumber (1/pmin)
         vth = min(sqrt(bolz*tele(i)*em)/hb,bolz*tele(i)/(zet(i)*eq*eq))
!-----classical Coulomb logarithm
         bxbn=max(1.d0,2*vth*sqrt(tele(i)*bolz/em)/max(omp,omega))
         clog= log(1.+bxbn)
!-----classical collision frequency
         ecf(i) = clog*coef*dens(i)*zet(i)/(abs(tele(i))**1.5)
      end do
      end subroutine ecf_class
!=======================================================================
!-----compute collision frequency by Eidmann-Hueller model described in
!-----K. Eidmann et al, Phys. Rew. E (2000) 1202
      subroutine ecf_kse(n,dens,tele,tion,zet,omega,ecf,vorf)
      implicit none
      integer,              intent(in)  :: n
      real*8,               intent(in)  :: omega,vorf
      real*8, dimension(n), intent(in)  :: dens,tele,tion,zet
      real*8, dimension(n), intent(out) :: ecf
      real*8                            :: coef,vth,bxbn,clog,ecfcla
      real*8                            :: vf,ecfeph,rwi,ecflim,omp
      integer                           :: i
!-----coefficient for classical collision frequency
      coef= sqrt(2*cpi)*4./3.*eq**4*em/(em*bolz)**1.5
      do i=1,n
!-----plasma frequency
         omp =sqrt(4*cpi*dens(i)*eq*eq/em)
!-----thermal wavenumber (1/pmin)
         vth = min(sqrt(bolz*tele(i)*em)/hb,bolz*tele(i)/(zet(i)*eq*eq))
!-----classical Coulomb logarithm
         bxbn=max(1.d0,2*vth*sqrt(tele(i)*bolz/em)/max(omp,omega))
         clog= log(1+bxbn)
!-----classical collision frequency
         ecfcla = clog*coef*dens(i)*zet(i)/(abs(tele(i))**1.5)
!-----collision frequency governed by scattering of electrons by phonons
!-----or lattice vibrations (Eidmann et AL. PRE 62 (2000) 1202, pag 1203,
!-----from Yakolev and Urpin, Astron. Zh. 57 (1980) 526) (vf=Fermi
!-----velocity)
         vf  = (hb/em)*(3*cpi*cpi*dens(i))**(1./3.)
         ecfeph = (vorf*2*eq**2*bolz/hb**2)*tion(i)/vf
!-----Wigner-Seitz radius
         rwi  = (0.75*zet(i)/dens(i)/cpi)**(1./3.)
!-----frequency upper limit
         ecflim=sqrt(vf**2+(bolz/em)*tele(i))/rwi
!-----average and limit collision frequency
         ecf(i) = min(1/(1/ecfcla+1/ecfeph),ecflim)
      end do
      end subroutine ecf_kse
!=======================================================================
!-----compute collision frequency by Meyer-ter-Vehn/Tronier model
      subroutine stossfrequenz(n,dens,tele,zet,omega,stoss)
      implicit none
      integer,              intent(in)  :: n
      real*8,               intent(in)  :: omega
      real*8, dimension(n), intent(in)  :: dens,tele,zet
      real*8, dimension(n), intent(out) :: stoss
      real*8, parameter :: ab=hb**2/em/eq**2
      real*8, parameter :: e0=eq**2/ab/bolz
      integer           :: i
      real*8            :: densab,telat,tf,tft,omegaev,omp,omegatf
      real*8            :: y,z,pauli,t0
      t0=ab**1.5*em**0.5/eq
      do i=1,n
         densab=dens(i)*ab**3
         telat=tele(i)/e0
!-----Fermi temperature
         tf=0.5*(3*cpi**2*densab)**(2./3.)
!-----'Theta' in Comp. Phys. Comm. 183 (2012) 637
         tft=telat/tf
!-----plasma frequency
         omegaev=omega*t0
         omp=sqrt(4*cpi*densab)
         omegatf=omegaev/tf
!-----chemical potential
         y=-log(3*sqrt(cpi)*tft**(3./2.)/4)+(0.25054*tft**(-1.858) &
         +0.0720*tft**(-0.929))/(1+0.25054*tft**(-0.858))
!-----Pauli blocking ('F' in Comp. Phys. Comm. 183 (2012) 637)
         z=omegatf/tft
         if(z<=0.000001)then
            pauli=1/(1+exp(-y))
         else
            if(y>100)then
               if(y-z>100)then
                  pauli=z/(1-exp(-z))
               else
                  pauli=(y-log(1+exp(y-z)))/(1-exp(-z))
               end if
            else
               pauli=log((1+exp(y))/(1+exp(y-z)))/(1-exp(-z))
            end if
         end if
         pauli=3*sqrt(cpi)/4*tft**1.5*pauli
!-----collision frequency
         stoss(i)=(2*sqrt(2*cpi)/t0)*zet(i)*densab/telat**1.5 &
         *log(1+1*(1.32/sqrt(2*cpi))*telat/(zet(i)**2 &
         *(omegaev**2+omp**2))**(1./3.))*pauli
      end do
      end subroutine stossfrequenz
!======================================================================!
!                                                                      !
!     eos and opacity interpolation routines                           !
!                                                                      !
!======================================================================!
!-----interpolate temperature, pressure, and ionization
      subroutine eosall(p,s,t)
      implicit none
      type(typesimparam),  intent(in)    :: p
      type(typesimstate),  intent(in)    :: s
      type(typesimthermo), intent(inout) :: t
      integer                            :: n,i,i1,i2,nn
      n        = p%n
!-----loop over material layers
      i1       = 1
      do i=2,n
         if(p%mid(i1)/=p%mid(i).or.i==n)then
!-----cells from i1 to i2 have the same material
            if(i==n)then
               i2      = i
            else
               i2      = i-1
            end if
            nn         = i2-i1+1
!-----electron temperature and pressure
            call eos(p%mat(p%mid(i1))%eeos,nn,s%r(i1:i2),s%ee(i1:i2), &
            t%pe(i1:i2),t%dpedr(i1:i2),t%dpedee(i1:i2), &
            t%te(i1:i2),t%dtedr(i1:i2),t%dtedee(i1:i2))
            t%te(i1:i2)     = max(t%te(i1:i2),p%tempmin)
!-----ion temperature and pressure
            call eos(p%mat(p%mid(i1))%ieos,nn,s%r(i1:i2),s%ei(i1:i2), &
            t%pi(i1:i2),t%dpidr(i1:i2),t%dpidei(i1:i2), &
            t%ti(i1:i2),t%dtidr(i1:i2),t%dtidei(i1:i2))
            t%ti(i1:i2)     = max(t%ti(i1:i2),p%tempmin)
!-----ionization
            if(p%mat(p%mid(i1))%z%nr>0)then
                call zbr(p%mat(p%mid(i1))%z,nn,s%r(i1:i2),t%te(i1:i2), &
                t%zi(i1:i2))
                t%zi(i1:i2) = min(t%zi(i1:i2),p%mat(p%mid(i1))%zi)
            else
                t%zi(i1:i2) = p%mat(p%mid(i1))%zi
            end if
!-----next layer
            i1=i
         end if
      end do 
      t%zi     = max(p%zmin,t%zi)
      end subroutine eosall
!======================================================================
!-----interpolate initial internal energy 
      subroutine eosenergy(eostab,r,t,e)
      implicit none
      type(eostable),      intent(in)    :: eostab
      real*8,              intent(in)    :: r,t
      real*8,              intent(out)   :: e
      real*8, dimension(eostab%ne)       :: ttab,ttab1,ttab2
      integer                            :: nr,ne,k
      real*8                             :: factor,e0
      nr       = eostab%nr
      ne       = eostab%ne
!-----locate indices for density
      k        = isearch(nr,eostab%rho,r)-1
      factor   = (r-eostab%rho(k))/(eostab%rho(k+1)-eostab%rho(k))
!-----cold energy
      e0       = eostab%e0(k)+factor*(eostab%e0(k+1)-eostab%e0(k))
!-----temperatures at given density
      ttab1    = eostab%t(k:nr*(ne-1)+k:nr)
      ttab2    = eostab%t(k+1:nr*(ne-1)+k+1:nr)
      ttab     = ttab1+(ttab2-ttab1)*factor
!-----locate indices for temperature
      k        = isearch(ne,ttab,t)-1
      factor   = (t-ttab(k))/(ttab(k+1)-ttab(k))
!-----obtain energy
      e        = eostab%de(k)+factor*(eostab%de(k+1)-eostab%de(k))+e0
      end subroutine eosenergy
!======================================================================
!-----interpolate pressure and temperature 
      subroutine eos(eostab,n,r,e,p,dpdr,dpde,t,dtdr,dtde)
      implicit none
      type(eostable),        intent(in)   :: eostab
      integer,               intent(in)   :: n
      real*8, dimension(n),  intent(in)   :: r,e
      real*8, dimension(n),  intent(out)  :: p,dpdr,dpde,t,dtdr,dtde
      real*8, dimension(n)                :: e0,de0dr,de
!-----obtain cold energy for the given densities
      call eoslin(eostab%nr,eostab%rho,eostab%e0,n,r,e0,de0dr)
!-----energy above the cold energy (must be positive)
      de       = max(0.d0,e-e0)
!-----obtain pressure and its derivatives
      call eosbin(eostab%nr,eostab%ne,eostab%rho,eostab%de,eostab%p, &
                  n,r,de,p,dpdr,dpde)
      dpdr     = dpdr-dpde*de0dr
!-----obtain temperature and its derivatives
      call eosbin(eostab%nr,eostab%ne,eostab%rho,eostab%de,eostab%t, &
                  n,r,de,t,dtdr,dtde)
      dtdr     = dtdr-dtde*de0dr
      end subroutine eos
!=======================================================================
!-----interpolate ionization
      subroutine zbr(grouptab,n,r,t,z)
      implicit none
      type(grouptable),      intent(in)   :: grouptab
      integer,               intent(in)   :: n
      real*8, dimension(n),  intent(in)   :: r,t
      real*8, dimension(n),  intent(out)  :: z
      real*8, dimension(n)                :: rlog,tlog,zlog
      real*8, dimension(n)                :: dummy1,dummy2
      rlog             = log10(r)
      tlog             = log10(t)-3.
      call eosbin(grouptab%nr,grouptab%nt,grouptab%rho,grouptab%t, &
                  grouptab%o,n,rlog,tlog,zlog,dummy1,dummy2)
      z                = 10.**zlog
      end subroutine zbr
!=======================================================================
!-----interpolate opacity
      subroutine opacit(p,r,t,f,op,or,eps)
      implicit none
      type(typesimparam),     intent(in)  :: p
      real*8, dimension(p%n), intent(in)  :: r,t
      real*8,                 intent(in)  :: f
      real*8, dimension(p%n), intent(out) :: op,or,eps
      integer                             :: n,i,i1,i2,nn
      n        = p%n
!-----loop over material layers
      i1       = 1
      do i=2,n
         if(p%mid(i1)/=p%mid(i).or.i==n)then
!-----cells from i1 to i2 have the same material
            if(i==n)then
               i2      = i
            else
               i2      = i-1
            end if
            nn         = i2-i1+1
!-----get Planck opacity
            call opa(p%mat(p%mid(i1))%planck,nn,r(i1),t(i1),f,op(i1))
!-----get Rosseland opacity
            call opa(p%mat(p%mid(i1))%ross,nn,r(i1),t(i1),f,or(i1))
!-----get emissivity 
            if(p%mat(p%mid(i1))%eps%ng>0)then
              call opa(p%mat(p%mid(i1))%eps,nn,r(i1),t(i1),f,eps(i1))
            else
               eps(i1:i2) = 1.0
            end if
!-----next layer
            i1         = i
         end if
      end do
      end subroutine opacit
!=======================================================================
!-----interpolate in multigroup table to obtain opacity
      subroutine opa(multitab,n,r,t,f,o)
      implicit none
      type(multitable),      intent(in)   :: multitab
      integer,               intent(in)   :: n
      real*8, dimension(n),  intent(in)   :: r,t
      real*8,                intent(in)   :: f
      real*8, dimension(n),  intent(out)  :: o
      type(grouptable)                    :: grouptab
      real*8, dimension(n)                :: rlog,tlog,olog
      real*8, dimension(n)                :: dummy1,dummy2
      real*8                              :: fa,fb
      integer                             :: i
      do i=1,multitab%ng
         fa            = multitab%tables(i)%fmin
         fb            = multitab%tables(i)%fmax
         if(fa==fb.or.(fa<=f.and.f<=fb))then
            rlog       = log10(r)
            tlog       = log10(t)
            grouptab   = multitab%tables(i)
            call eosbin(grouptab%nr,grouptab%nt,grouptab%rho, &
            grouptab%t,grouptab%o,n,rlog,tlog,olog,dummy1,dummy2)
            o          = 10.**olog
            return
         end if
      end do
      print *,'ERROR in subroutine opa'
      print *,'      frequency not found in tables'
      stop
      end subroutine opa
!=======================================================================
!-----linear interpolation: y(x), dydx(x) from table yt(xt)
      subroutine eoslin(nt,xt,yt,n,x,y,dydx)
      implicit none
      integer,               intent(in)  :: nt,n
      real*8, dimension(nt), intent(in)  :: xt,yt
      real*8, dimension(n),  intent(in)  :: x
      real*8, dimension(n),  intent(out) :: y,dydx
      integer                            :: ind(n)
      real*8, dimension(n)               :: x1,x2,y1,y2
      integer                            :: i
      do i=1,n
         ind(i)        = isearch(nt,xt,x(i))
      end do
      x1               = xt(ind-1)
      x2               = xt(ind)
      y1               = yt(ind-1)
      y2               = yt(ind)
      dydx             = (y2-y1)/(x2-x1)
      y                = y1+dydx*(x-x1)
      end subroutine eoslin
!=======================================================================
!-----bilinear interpolation: z(x,y), dzdx(x,y), and dzdy(x,y) from 
!-----table zt(xt,xt)
      subroutine eosbin(nx,ny,xt,yt,zt,n,x,y,z,dzdx,dzdy)
      implicit none
      integer,                   intent(in)  :: n,nx,ny
      real*8,  dimension(nx),    intent(in)  :: xt
      real*8,  dimension(ny),    intent(in)  :: yt
      real*8,  dimension(nx*ny), intent(in)  :: zt
      real*8,  dimension(n),     intent(in)  :: x,y
      real*8,  dimension(n),     intent(out) :: z,dzdx,dzdy
      integer, dimension(n)                  :: ix,iy
      real*8,  dimension(n)                  :: x1,x2,y1,y2,z1,z2,z3,z4
      real*8,  dimension(n)                  :: alpha,beta,aa,bb,cc
      integer                                :: i
      do i=1,n
         ix(i)         = isearch(nx,xt,x(i))
         iy(i)         = isearch(ny,yt,y(i))
      end do   
      x1               = xt(ix-1)
      x2               = xt(ix)
      y1               = yt(iy-1)
      y2               = yt(iy)
      z1               = zt(ix-1+(iy-2)*nx)
      z2               = zt(ix  +(iy-2)*nx)
      z3               = zt(ix-1+(iy-1)*nx)
      z4               = zt(ix  +(iy-1)*nx)
      alpha            = (x-x1)/(x2-x1)
      beta             = (y-y1)/(y2-y1)
      aa               = z2-z1
      bb               = z3-z1
      cc               = z4-z3-z2+z1
      z                = z1+aa*alpha+bb*beta+cc*alpha*beta
      dzdx             = (aa+cc*beta)/(x2-x1)
      dzdy             = (bb+cc*alpha)/(y2-y1)
      end subroutine eosbin
!=======================================================================
!-----find x(isearch-1) <= xt <= x(isearch), forcing 2 <= isearch <= n
      integer function isearch(n,x,xt)
      implicit none
      integer, intent(in) :: n
      real*8, intent(in)  :: x(n),xt
      integer             :: i
      do i=2,n
         if(x(i)>xt) exit
      end do  
      isearch=min(i,n)
      end function isearch
!======================================================================!
!                                                                      !
!     management of 'typesimstate' and 'typesimthermo' structures      !
!                                                                      !
!======================================================================!
!-----allocate memory for a typesimstate structure
      subroutine allocsimstate(n,ng,s)
      implicit none
      integer,            intent(in)    :: n,ng
      type(typesimstate), intent(inout) :: s
      allocate(s%x(n+1),s%v(n+1))
      allocate(s%r(n),s%f(n),s%ee(n),s%ei(n),s%ea(n))
      allocate(s%wa(ng+1))
      memory     = memory+8
      end subroutine allocsimstate
!=======================================================================
!-----copy members of a typesimstate structure
      subroutine copysimstate(s1,s2)
      implicit none
      type(typesimstate), intent(in)    :: s1
      type(typesimstate), intent(inout) :: s2
      s2%time   = s1%time
      s2%x      = s1%x
      s2%v      = s1%v
      s2%r      = s1%r
      s2%ee     = s1%ee
      s2%ei     = s1%ei
      s2%ea     = s1%ea
      s2%f      = s1%f
      s2%wa     = s1%wa
      end subroutine copysimstate
!=======================================================================
!-----free memory for a typesimstate structure
      subroutine freesimstate(s)
      implicit none
      type(typesimstate), intent(inout) :: s
      deallocate(s%x,s%v)
      deallocate(s%r,s%f,s%ee,s%ei,s%ea)
      deallocate(s%wa)
      memory     = memory-8
      end subroutine freesimstate
!=======================================================================
!-----allocate memory for a typesimthermo structure
      subroutine allocsimthermo(n,t)
      implicit none
      integer,             intent(in)    :: n
      type(typesimthermo), intent(inout) :: t
      allocate(t%te(n),t%pe(n))
      allocate(t%dtedee(n),t%dtedr(n),t%dpedee(n),t%dpedr(n))
      allocate(t%ti(n),t%pi(n))
      allocate(t%dtidei(n),t%dtidr(n),t%dpidei(n),t%dpidr(n))
      allocate(t%zi(n))
      memory     = memory+13
      end subroutine allocsimthermo
!=======================================================================
!-----copy members of a typesimthermo structure
      subroutine copysimthermo(t1,t2)
      implicit none
      type(typesimthermo), intent(in)    :: t1
      type(typesimthermo), intent(inout) :: t2
      t2%te    = t1%te
      t2%pe    = t1%pe
      t2%dtedee= t1%dtedee
      t2%dtedr = t1%dtedr
      t2%dpedee= t1%dpedee
      t2%dpedr = t1%dpedr
      t2%ti    = t1%ti
      t2%pi    = t1%pi
      t2%dtidei= t1%dtidei
      t2%dtidr = t1%dtidr
      t2%dpidei= t1%dpidei
      t2%dpidr = t1%dpidr
      t2%zi    = t1%zi
      end subroutine copysimthermo
!=======================================================================
!-----free memory for a typesimthermo structure
      subroutine freesimthermo(t)
      implicit none
      type(typesimthermo), intent(inout) :: t
      deallocate(t%te,t%pe)
      deallocate(t%dtedee,t%dtedr,t%dpedee,t%dpedr)
      deallocate(t%ti,t%pi)
      deallocate(t%dtidei,t%dtidr,t%dpidei,t%dpidr)
      deallocate(t%zi)
      memory     = memory-13
      end subroutine freesimthermo
!======================================================================!
!                                                                      !
!     mathematical routine                                             !
!                                                                      !
!======================================================================!
!-----Solve a tridiagonal system of equations [a]{x}={b}
!-----a(1,i) - upper diagonal (from a(1,2) to (1,n))
!-----a(2,i) - main diagonal
!-----a(3,i) - lower diagonal (from a(3,1) to (3,n-1))
!-----x is returned in the place of b
      subroutine tridia(n,a,b)
      implicit none
      integer, intent(in)    :: n
      real*8,  intent(inout) :: a(3,n),b(n)
      integer                :: i
      real*8                 :: aux
      do i=1,n-1
         if(a(2,i)==0) call tridia_error(1)
         aux           = a(3,i)/a(2,i)
         a(2,i+1)      = a(2,i+1)-aux*a(1,i+1)
         b(i+1)        = b(i+1)-aux*b(i)
      end do
      if(a(2,n)==0)    call tridia_error(2)
      b(n)=b(n)/a(2,n)
      do i=n-1,1,-1
         if(a(2,i)==0) call tridia_error(3)
         b(i)          = (b(i)-b(i+1)*a(1,i+1))/a(2,i)
      end do
      end subroutine tridia
!=======================================================================
!-----print error message
      subroutine tridia_error(i)
      implicit none
      integer, intent(in) :: i
      print *,'ERROR in subroutine tridia'
      print *,'      singular matrix found in ',i
      stop
      end subroutine tridia_error
!======================================================================!
!                                                                      !
!     initialization routines                                          !
!                                                                      !
!======================================================================!
!-----initialization routine
      subroutine init(p,s,e)
      implicit none
      type(typesimparam),  intent(out) :: p
      type(typesimstate),  intent(out) :: s
      type(typesimenergy), intent(out) :: e
      type(typesimlayers)              :: layers
      integer                          :: i
  
      call get_command_argument(1,p%igafile)
  
      call read_parameters(p)
      call ga_read_laser_3d(p%laser_3d,p%igafile)
      call ga_read_laser_wkb(p%laser_wkb,p%igafile)
      call ga_read_laser_maxwell(p%laser_maxwell,p%igafile)
      call ga_read_pulse(p%pulse,p%igafile)
      !call read_laser_3d(p%laser_3d)
      !call read_laser_wkb(p%laser_wkb)
      !call read_laser_maxwell(p%laser_maxwell)
      !call read_pulse(p%pulse)
      call ga_read_material_definitions(p)
      !call read_material_definitions(p)
      call read_material_data(p)
      call ga_read_layer_definition(layers,p%igafile)  
     ! call read_layer_definition(layers)
  
      call init_state(p,layers,s)    
      call update_energies(p,s,-1.0d0,e)
      !wfyuan, 2020/10/27,switch off print for ga
      !call print_heading
      !call print_parameters(p)
      !call print_layer_definition(layers)
      !do i=1,p%nmat
      !   call print_material(p%mat(i),1)
      !end do
      !call print_laser_3d(p%laser_3d)
      !call print_laser_wkb(p%laser_wkb)
      !call print_laser_maxwell(p%laser_maxwell)
      !call print_pulse(p%pulse)
      end subroutine init
!=======================================================================  
!-----finishing routine
      subroutine finish(p,s,e)
      implicit none
      type(typesimparam),  intent(inout) :: p
      type(typesimenergy), intent(inout) :: e
      type(typesimstate),  intent(inout) :: s
      !-----wfyuan, 2020/10/27, switch off print for ga
!-----energy balance
      !call print_energies(e)
  
!-----wfyuan, 2020/10/27, print fitness
      call ga_print_fitness(p)
  
!-----free memory
      call deallocate_parameters(p)
      call freesimstate(s)
      call deallocate_energies(e)
!-----final message
      !call print_tail
      end subroutine finish
!=======================================================================
!-----read pulse shape
      subroutine read_pulse(p)
      implicit none
      type(typesimpulse), intent(out) :: p
      integer                         :: ntab,mtab,icode
      real*8, dimension(size(p%ttab)) :: ttab,ptab
      real*8                          :: expo
      namelist /pulse_shape/ ntab,mtab,expo,ttab,ptab
      call namelistblock('pulse_shape',1,icode)
      if(icode==0)return
      open(3,file='block',form='formatted')
      read(3,pulse_shape)
      close(3)
      p%ntab   = ntab
      p%mode   = mtab
      p%expo   = expo
      p%ttab   = ttab
      p%ptab   = ptab
      end subroutine read_pulse
!=======================================================================
!-----read laser parameters (ortogonal rays)
      subroutine read_laser_wkb(p)
      implicit none
      type(typesimlaserwkb), intent(out) :: p
      integer :: inter,itype,icode
      real*8  :: pimax,pitime,wl,delta
      namelist /pulse_wkb/ inter,pimax,pitime,wl,delta,itype
      call namelistblock('pulse_wkb',1,icode)
      if(icode==0)then
         p%enable = 0
         return
      end if
      open(3,file='block',form='formatted')
      read(3,pulse_wkb)
      close(3)
      p%enable = 1
      p%inter  = inter
      p%pimax  = pimax
      p%pitime = pitime
      p%wl     = wl
      p%itype  = itype
      p%delta  = delta
      end subroutine read_laser_wkb
!=======================================================================
!-----read laser parameters (3D ray tracing)
      subroutine read_laser_3d(p)
      implicit none
      type(typesimlaser3d), intent(out) :: p
      integer :: itype,nr,icode
      real*8  :: pimax,pitime,wl,rmax,fwhm,bexp
      namelist /pulse_3d/ pimax,pitime,wl,itype,nr,rmax,fwhm,bexp
      call namelistblock('pulse_3d',1,icode)
      if(icode==0)then
         p%enable = 0
         return
      end if
      open(3,file='block',form='formatted')
      read(3,pulse_3d)
      close(3)
      p%enable = 1
      p%pimax  = pimax
      p%pitime = pitime
      p%wl     = wl
      p%itype  = itype
      p%nr     = nr
      p%rmax   = rmax
      p%fwhm   = fwhm
      p%bexp   = bexp 
      end subroutine read_laser_3d
!=======================================================================
!-----read laser parameters (Maxwell's equations solver)
      subroutine read_laser_maxwell(p)
      implicit none
      type(typesimlasermaxwell), intent(out) :: p
      integer          :: inter,itype,idep,icode
      real*8           :: pimax,pitime,wl,angle
      character(len=1) :: pol
      namelist /pulse_maxwell/ inter,pimax,pitime,wl,itype, &
                               idep,angle,pol
      call namelistblock('pulse_maxwell',1,icode)
      if(icode==0)then
         p%enable = 0
         return
      end if
      open(3,file='block',form='formatted')
      read(3,pulse_maxwell)
      close(3)
      p%enable = 1
      p%pimax  = pimax
      p%inter  = inter
      p%pitime = pitime
      p%wl     = wl
      p%itype  = itype
      p%idep   = idep
      p%angle  = angle
      p%pol    = 1
      if(pol=='s'.or.pol=='S') p%pol=2
      end subroutine read_laser_maxwell
!=======================================================================
!-----read simulation parameters
      subroutine read_parameters(p)
      implicit none
      type(typesimparam), intent(inout) :: p
      integer :: igeo,nfuel,iright,ileft,itype,iradia,ihydro,inegpre
      integer :: model,nsplit,iwctrl,nexit,ns_bout,ns_aout,irad_left
      integer :: irad_right,iheation,icode,nreduce
      real*8  :: texit,xmin,alphal,alphar,betal,betar,tau,trad
      real*8  :: fheat,fei,flaser,flf,zmin,dtmin,dtmax,dtinit,dtrvar
      real*8  :: dttevar,dttivar,dtbreak,dtfactor,dt_bout,dt_aout
      namelist /parameters/ igeo,xmin,texit,nfuel, &
       iright,ileft,alphal,alphar,betal,betar,itype,tau,trad, &
       iradia,ihydro,model,fheat,fei,flaser,zmin,flf,inegpre, &
       nsplit,dtmin,dtmax,dtinit,dtrvar,dttevar,dttivar,dtbreak, &
       dtfactor,nexit,iwctrl,iheation,nreduce, &
       dt_bout,ns_bout,dt_aout,ns_aout,irad_left,irad_right
      !call namelistblock('parameters',1,icode)
      call ga_namelistblock('parameters',1,icode,p%igafile)
      if(icode==0)then
         print *,'ERROR in subroutine read_parameters'
         print *,'      namelist ''parameters'' not found'
         stop
      end if
      open(3,file='block_'//trim(p%igafile),form='formatted')
      !open(3,file='block',form='formatted')
      read(3,parameters)
      close(3)
      p%igeo           = igeo
      p%xmin           = xmin
      p%texit          = texit 
      p%nfuel          = nfuel
      p%iright         = iright
      p%ileft          = ileft
      p%alphal         = alphal
      p%alphar         = alphar
      p%betal          = betal
      p%betar          = betar
      p%itype          = itype
      p%tau            = tau
      p%trad           = trad
      p%iradia         = iradia
      p%ihydro         = ihydro
      p%iheation       = iheation
      p%model          = model
      p%fheat          = fheat
      p%fei            = fei
      p%flaser         = flaser
      p%zmin           = zmin
      p%flf            = flf
      p%inegpre        = inegpre
      p%nsplit         = nsplit
      p%dtmin          = dtmin
      p%dtmax          = dtmax
      p%dtinit         = dtinit
      p%dtrvar         = dtrvar
      p%dttevar        = dttevar
      p%dttivar        = dttivar
      p%dtbreak        = dtbreak
      p%dtfactor       = dtfactor
      p%nexit          = nexit 
      p%iwctrl         = iwctrl
      p%dt_aout        = dt_aout
      p%ns_aout        = ns_aout
      p%dt_bout        = dt_bout
      p%ns_bout        = ns_bout
      p%nreduce        = nreduce
      p%irad_left      = irad_left
      p%irad_right     = irad_right
      end subroutine read_parameters
!=======================================================================
!-----read material definitions
      subroutine read_material_definitions(p)
      implicit none
      type(typesimparam), intent(inout) :: p
      integer                           :: icode
      p%nmat=0
      do
         if(p%nmat>size(p%mat))then
            print *,'ERROR in subroutine read_material_definitions'
            print *,'      more than ',size(p%mat),' materials'
            stop
         end if
         call namelistblock('material',p%nmat+1,icode)
         if(icode==0)return
         p%nmat        = p%nmat+1
         call read_one_material_definition(p%mat(p%nmat))
      end do
      end subroutine read_material_definitions
!=======================================================================
!-----read one material definition
      subroutine read_one_material_definition(m)
      implicit none
      type(typematerial),    intent(out) :: m
      character(len=80)   :: name
      real*8              :: ai,zi
      character(len=80)   :: eeos_file,ieos_file,z_file
      character(len=80)   :: planck_file,ross_file,eps_file
      integer             :: eeos_id,ieos_id,z_id
      integer             :: planck_id,ross_id,eps_id
      namelist /material/ name,ai,zi, &
       eeos_file,eeos_id,ieos_file,ieos_id,z_file,z_id, &
       planck_file,planck_id,ross_file,ross_id,eps_file,eps_id
      eeos_file   = "fort.13"
      ieos_file   = "fort.13"
      z_file      = "fort.13"
      planck_file = "fort.13"
      ross_file   = "fort.13"
      eps_file    = "fort.13"
      eeos_id     = 0
      ieos_id     = 0
      z_id        = 0
      planck_id   = 0
      ross_id     = 0
      eps_id      = 0
      open(3,file='block',form='formatted')
      read(3,material)
      close(3)
      m%name           = name
      m%ai             = ai
      m%zi             = zi
      m%eeos_file      = eeos_file
      m%eeos_id        = eeos_id
      m%ieos_file      = ieos_file
      m%ieos_id        = ieos_id
      m%z_file         = z_file
      m%z_id           = z_id
      m%planck_file    = planck_file
      m%planck_id      = planck_id
      m%ross_file      = ross_file
      m%ross_id        = ross_id
      m%eps_file       = eps_file
      m%eps_id         = eps_id
      end subroutine read_one_material_definition
!=======================================================================
!-----read layered initial structure
      subroutine read_layer_definition(l)
      implicit none
      type(typesimlayers), intent(out) :: l
      character(len=80)                :: material
      integer                          :: nc,icode
      real*8                           :: thick,r0,te0,ti0,zonpar
      namelist /layer/  nc,thick,r0,te0,zonpar,ti0,material
      l%nl=0
      do
         call namelistblock('layer',l%nl+1,icode)
         if(icode==0)return
         open(3,file='block',form='formatted')
         read(3,layer)
         close(3)
         l%nl           = l%nl+1
         if(l%nl>size(l%nc))then
            print *,'ERROR in subroutine read_layer_definition'
            print *,'      more than ',size(l%nc),' layers'
            stop
         end if
         l%nc     (l%nl)  = nc
         l%thick  (l%nl)  = thick
         l%zonpar (l%nl)  = zonpar
         l%r0     (l%nl)  = r0
         l%te0    (l%nl)  = te0
         l%ti0    (l%nl)  = ti0
         l%material(l%nl) = material
      end do
      end subroutine read_layer_definition
!=======================================================================
!-----tranfer namelist block 'name' from file 'fort.12' to file 'block_index'
      subroutine ga_namelistblock(name,number,icode,fileindex)
      implicit none
      character(len=*), intent(in)  :: name
      integer,          intent(in)  :: number
      integer,          intent(out) :: icode
      character(len=32)             :: fileindex
      character(len=80)             :: inputline
      integer                       :: i,j
      open(1,file='fort.12',form='formatted')
      open(2,file='block_'//trim(fileindex),form='formatted')
      i    = 0
      do
         read(1,'(a)',end=2) inputline
         j = locblock(name,inputline)
         if(j>0.and.i==number)exit
         if(j==2) i=i+1
         if(i==number)then
            j=scan(inputline,"!")
            if(j>0) inputline=inputline(:j-1)
            write(2,'(a)')trim(inputline)
         end if
      end do  
2     continue
      close(1)
      close(2)
      icode=0
      if(i==number)icode=1
      end subroutine ga_namelistblock
!=======================================================================
!-----tranfer namelist block 'name' from file 'fort.12' to file 'block'
      subroutine namelistblock(name,number,icode)
      implicit none
      character(len=*), intent(in)  :: name
      integer,          intent(in)  :: number
      integer,          intent(out) :: icode
      character(len=80)             :: inputline
      integer                       :: i,j
      open(1,file='fort.12',form='formatted')
      open(2,file='block',form='formatted')
      i    = 0
      do
         read(1,'(a)',end=2) inputline
         j = locblock(name,inputline)
         if(j>0.and.i==number)exit
         if(j==2) i=i+1
         if(i==number)then
            j=scan(inputline,"!")
            if(j>0) inputline=inputline(:j-1)
            write(2,'(a)')trim(inputline)
         end if
      end do  
2     continue
      close(1)
      close(2)
      icode=0
      if(i==number)icode=1
      end subroutine namelistblock
!=======================================================================
!-----check if 'line' contains namelist block 'name'
!-----0 - no
!-----1 - namelist begin line
!-----2 - namelist begin line with name 'name'
      integer function locblock(name,line)
      implicit none
      character(len=*), intent(in) :: name,line
      integer                      :: i,i1,nc,k1,k2
      integer, parameter           :: idelta=ichar('A')-ichar('a')
      locblock = 0
      i1       = scan(line,'&')
      if(i1==0)return
      locblock = 1
      nc       = scan(line(i1:),' ')-2
      if(nc/=len(name))return
      do i=1,nc
         k1    = ichar(line(i1+i:i1+i))
         k2    = ichar(name(i:i))
         if(k1/=k2.and.k1/=k2+idelta)return
      end do
      locblock = 2
      end function locblock
!=======================================================================
!-----read material data from table files
      subroutine read_material_data(p)
      implicit none
      type(typesimparam), intent(inout) :: p
      type(typematerial)                :: m
      integer                           :: i
      do i=1,p%nmat
         m=p%mat(i)
         call load_table(1,m%eeos_file,  m%eeos_id,  m)
         call load_table(2,m%ieos_file,  m%ieos_id,  m)
         call load_table(3,m%z_file,     m%z_id,     m)
         call load_table(4,m%planck_file,m%planck_id,m)
         call load_table(5,m%ross_file,  m%ross_id,  m)
         call load_table(6,m%eps_file,   m%eps_id,   m)
         if(p%inegpre/=1)then
            m%eeos%p = max(m%eeos%p,0.d0)
            m%ieos%p = max(m%ieos%p,0.d0)
         end if
         p%mat(i)=m
      end do
      end subroutine read_material_data
!=======================================================================
!-----load one table
      subroutine load_table(mode,file,id,m)
      implicit none
      integer, intent(in)               :: mode,id
      character (len=*)                 :: file
      type(typematerial), intent(inout) :: m
      integer                           :: ind,dummy
      if(mode==3.and.id<0)return
      if(mode==6.and.id<0)return
      ind      = 1
      open(3,file=trim(file),status='old',form='formatted',err=1)
      select case(mode)
      case(1)
         call load_table_eos(3,id,ind,m%eeos)
      case(2)
         call load_table_eos(3,id,ind,m%ieos)
      case(3)
         call load_table_group(3,id,dummy,m%z)
      case(4)
         call load_table_multi(3,id,dummy,m%planck)
         if(m%planck%ng==0)ind = 0
      case(5)
         call load_table_multi(3,id,dummy,m%ross)
         if(m%ross%ng==0)ind = 0
      case(6)
         call load_table_multi(3,id,dummy,m%eps)
      end select
      close(3)
      if(ind==0)then
          print *,'ERROR in subroutine read_material_data'
          print *,'      table not found'
          print *,'      file = ',trim(file)
          print *,'      id   = ',id
          stop
      end if
      return
1     continue
      print *,'ERROR in subroutine read_material_data'
      print *,'      file not found'
      print *,'      file = ',trim(file)
      stop
      end subroutine load_table
!=======================================================================
!-----load equation-of-state table from table file
!-----unit   - fortran I/O unit
!-----idin   - table code (0 means the first EOS found)
!-----ind    - 1 - table read
!-----         0 - EOF found
!-----eostab - equation-of-state table
      subroutine load_table_eos(unit,idin,ind,eostab)
      implicit none
      integer,          intent(in)    :: unit,idin
      integer,          intent(out)   :: ind
      type(eostable),   intent(out)   :: eostab
      real*8,           allocatable   :: buffer(:)
      integer                         :: nr,ne,id,itype,numrec,n
      do
!-----read table heading
         call load_table_heading(unit,id,itype,nr,ne,numrec)
!-----EOF reached
         if(itype==0)then
            ind         = 0
            return
         end if
!-----found table with requested id and correct type
         if((id==idin.or.idin==0).and.itype==1)then
            allocate(buffer(numrec*4))
            memory      = memory+1
            call load_table_read(unit,numrec,buffer)
            eostab%nr   = nr
            eostab%ne   = ne
            allocate(eostab%rho(nr))
            allocate(eostab%de(ne))
            allocate(eostab%e0(nr))
            allocate(eostab%p(nr*ne))
            allocate(eostab%t(nr*ne))
            memory      = memory+5
            n           = 0
            eostab%rho  = buffer(n+1:n+nr)
            n           = n+nr
            eostab%de   = buffer(n+1:n+ne)*1.e12
            n           = n+ne
            eostab%e0   = buffer(n+1:n+nr)*1.e12
            n           = n+nr
            eostab%p    = buffer(n+1:n+ne*nr)*1.e12
            n           = n+ne*nr
            eostab%t    = buffer(n+1:n+ne*nr)/11605.
            deallocate(buffer)
            memory      = memory-1
            ind         = 1
            return
         else 
!-----skip this table
            call load_table_skip(unit,numrec)
         end if
      end do
      end subroutine load_table_eos
!=======================================================================
!-----load group table (opac., emiss., or ionization) from table-file
!-----unit   - fortran I/O unit
!-----idin   - table code (0 means the first group-table found)
!-----ind    - 1 - table read
!-----         0 - EOF found
!-----grouptab - one group table
      subroutine load_table_group(unit,idin,ind,grouptab)
      implicit none
      integer,          intent(in)    :: unit,idin
      integer,          intent(out)   :: ind
      type(grouptable), intent(out)   :: grouptab
      real*8,           allocatable   :: buffer(:)
      integer                         :: id,itype,nr,nt,numrec,n
      do
!-----read table heading
         call load_table_heading(unit,id,itype,nr,nt,numrec)
!-----EOF reached
         if(itype==0)then
            ind        = 0
            return
         end if
!-----found table with requested id and correct type
         if((id==idin.or.idin==0).and.itype>1)then
            allocate(buffer(numrec*4))
            memory      = memory+1
            call load_table_read(unit,numrec,buffer)
            grouptab%nr= nr
            grouptab%nt= nt
            allocate(grouptab%rho(nr))
            allocate(grouptab%t(nt))
            allocate(grouptab%o(nr*nt))
            memory      = memory+3
            n          = 0
            if(itype>=6)then
               grouptab%fmin = buffer(1)
               grouptab%fmax = buffer(2)
               n       = n+4
            else
               grouptab%fmin = 0.0
               grouptab%fmax = 0.0
            end if
            grouptab%rho= buffer(n+1:n+nr)
            n           = n+nr
            grouptab%t  = buffer(n+1:n+nt)
            n           = n+nt
            grouptab%o  = buffer(n+1:n+nr*nt)
            deallocate(buffer)
            memory      = memory-1
            ind         = 1
            return
         else 
!-----skip this table
            call load_table_skip(unit,numrec)
         end if
      end do  
      end subroutine load_table_group
!=======================================================================
!-----load multigroup table (opacity, emissivity) from table file
!-----unit   - fortran I/O unit
!-----idin   - table code (0 means load all tables found)
!-----ind    - 1 - single group table read
!-----         0 - EOF reached (zero or more single group tables loaded)
!-----multitab - multi group table
      subroutine load_table_multi(unit,idin,ind,multitab)
      implicit none
      integer,          intent(in)    :: unit,idin
      integer,          intent(out)   :: ind
      type(multitable), intent(out)   :: multitab
      type(grouptable)                :: grouptab
      type(grouptable), dimension(:), pointer :: aux => null()
      multitab%ng = 0
      do
!-----read one group table
         call load_table_group(unit,idin,ind,grouptab)
!-----EOF reached
         if(ind==0)return
!-----put group table in multi table
         if(multitab%ng==0) then
            allocate(multitab%tables(1))
            memory     = memory+1
         else
            allocate(aux(multitab%ng))
            aux = multitab%tables
            deallocate(multitab%tables)
            allocate(multitab%tables(multitab%ng+1))
            multitab%tables(1:multitab%ng) = aux
            deallocate(aux)
         end if
         multitab%ng = multitab%ng+1
         multitab%tables(multitab%ng) = grouptab
!-----for a single group table, return inmediatelly
         if(grouptab%fmax==grouptab%fmin)then
            ind=1
            return
         end if
      end do
      end subroutine load_table_multi
!=======================================================================
!-----read table heading from table file
!-----unit   - fortran I/O unit
!-----id     - table code
!-----itype  - table type
!-----         0 - EOF reached
!-----         1 - inverted EOS
!-----         2 - ionization
!-----         3 - Planck opacity (frequency averaged)
!-----         4 - Rosseland opacity (frequency averaged)
!-----         5 - emissivity (frequency averaged)
!-----         6 - Planck opacity (frequency averaged)
!-----         7 - Rosseland opacity (frequency averaged)
!-----         8 - emissivity (frequency averaged)
!-----nr,nt  - dimensions (number of dens., number of temp./energies)
!-----numrec - number of lines (excluding heading)
      subroutine load_table_heading(unit,id,itype,nr,nt,numrec)
      implicit none
      integer,         intent(in)    :: unit
      integer,         intent(out)   :: id,itype,nr,nt,numrec
      real*8                         :: a1,a3,a4
      character(len=15)              :: a2
      read(unit,'(e15.0,a,2e15.8)',end=1,err=1)a1,a2,a3,a4
      id       = int(a1)
      nr       = int(a3)
      nt       = int(a4)
      itype    = 1
      if(a2==' 0.60000000E+01')itype=2
      if(a2==' 0.60000000e+01')itype=2
      if(a2=='PLANCK 1       ')itype=3
      if(a2==' 0.70000000E+01')itype=3
      if(a2==' 0.70000000e+01')itype=3
      if(a2=='ROSSELAND 1    ')itype=4
      if(a2==' 0.40000000E+01')itype=4
      if(a2==' 0.40000000e+01')itype=4
      if(a2=='EPS 1          ')itype=5
      if(a2=='PLANCK M       ')itype=6
      if(a2=='ROSSELAND M    ')itype=7
      if(a2=='EPS M          ')itype=8
      if(itype==1)then
         numrec     = (2*nr+nt+2*nr*nt+3)/4
      else if(itype<6)then
         numrec     = (nr+nt+nt*nr+3)/4
      else
         numrec     = (nr+nt+nt*nr+3)/4+1
      end if
      return
1     continue
      itype    = 0         
      end subroutine load_table_heading
!=======================================================================
!-----read numrec records from table file
      subroutine load_table_read(unit,numrec,buffer)
      implicit none
      integer,         intent(in)    :: unit,numrec
      real*8,          intent(out)   :: buffer(numrec*4)
      integer                        :: i
      do i=1,numrec
         read(unit,'(4e15.8)',end=1,err=2)buffer(i*4-3:i*4)
      end do
      return
1     continue
      print *,'ERROR in subroutine load_table_read'
      print *,'      premature eof'
      stop
2     continue
      print *,'ERROR in subroutine load_table_read'
      print *,'      input error'
      stop
      end subroutine load_table_read
!=======================================================================
!-----skip numrec records in table file
      subroutine load_table_skip(unit,numrec)
      implicit none
      integer,         intent(in)    :: unit,numrec
      integer                        :: i
      character(len=15)              :: a2
      do i=1,numrec
          read(unit,'(a)',end=1)a2
      end do
      return
1     continue
      print *,'ERROR in subroutine load_table_skip'
      print *,'      premature eof'
      stop
      end subroutine load_table_skip
!=======================================================================
!-----deallocate material data
      subroutine deallocate_parameters(p)
      implicit none
      type(typesimparam), intent(inout) :: p
      type(typematerial)                :: m
      integer                           :: i
      if(associated(p%dm)) then
         deallocate(p%dm)
         memory= memory-1
      end if
      if(associated(p%ai)) then
         deallocate(p%ai)
         memory= memory-1
      end if
      if(associated(p%mid)) then
         deallocate(p%mid)
         memory= memory-1
      end if
      if(associated(p%frei)) then
         deallocate(p%frei)
         memory= memory-1
      end if
      do i=1,p%nmat
         m=p%mat(i)
         call deallocate_table_eos(m%eeos)
         call deallocate_table_eos(m%ieos)
         call deallocate_table_group(m%z)
         call deallocate_table_multi(m%planck)
         call deallocate_table_multi(m%ross)
         call deallocate_table_multi(m%eps)
      end do
      end subroutine deallocate_parameters
!=======================================================================
      subroutine deallocate_table_eos(eostab)
      implicit none
      type(eostable),   intent(inout)   :: eostab
      deallocate(eostab%rho,eostab%e0,eostab%de,eostab%p,eostab%t)
      memory   = memory-5
      end subroutine deallocate_table_eos
!=======================================================================
      subroutine deallocate_table_multi(multitab)
      implicit none
      type(multitable), intent(inout)   :: multitab
      integer                           :: i
      do i=1,multitab%ng
         call deallocate_table_group(multitab%tables(i))
      end do
      if(associated(multitab%tables))then
         deallocate(multitab%tables)
         memory   = memory-1
      end if
      end subroutine deallocate_table_multi
!=======================================================================
      subroutine deallocate_table_group(grouptab)
      implicit none
      type(grouptable), intent(inout)   :: grouptab
      if(grouptab%nr/=0)then
         deallocate(grouptab%rho,grouptab%t,grouptab%o)
         memory= memory-3
      end if
      end subroutine deallocate_table_group
!=======================================================================
      subroutine deallocate_energies(e)
      implicit none
      type(typesimenergy), intent(inout) :: e
      if(associated(e%d)) then
         deallocate(e%d)
         memory= memory-1
      end if
      if(associated(e%strahl)) then
         deallocate(e%strahl)
         memory= memory-1
      end if
      end subroutine deallocate_energies
!=======================================================================
!-----initialize state vector 's' and upgrade structure 'p'
      subroutine init_state(p,l,s)
      implicit none
      type(typesimparam),  intent(inout) :: p
      type(typesimstate),  intent(out)   :: s
      type(typesimlayers), intent(inout) :: l    ! wfyuan 2020/10/27, for ga
      !type(typesimlayers), intent(in)    :: l
      real*8                             :: ee,ei,wb(1)
      integer                            :: n,i,j,mid
  
      integer*4                          :: narg
      character(len=32)                  :: varg
      character(len=32)                  :: filename
      real*8, dimension(10)              :: mylayerthick
      real*8, pointer                    :: gaxacc(:)
      integer*4                          :: ntab
      real*8                             :: elaser = 0
  
!-----get size and reallocate memory
      n        = sum(l%nc(1:l%nl))
      p%n      = n
      if(associated(p%dm)) then
         deallocate(p%dm)
         memory= memory-1
      end if
      if(associated(p%ai)) then
         deallocate(p%ai)
         memory= memory-1
      end if
      if(associated(p%mid))then
         deallocate(p%mid)
         memory= memory-1
      end if
      allocate(p%dm(n),p%ai(n),p%mid(n)) 
      memory   = memory+3
        
!-----get group frequencies
      call frequencies(p)
!-----adjust number of subcycles
      p%msplit = max(1,min(p%ng,p%nsplit))
!-----allocate memory
      call allocsimstate(p%n,p%ng,s)
        
!-----wfyuan 2020/10/27, read parameters sety by genetic algorithm
      !call get_command_argument(1,p%igafile)
      filename = "inp_"//trim(p%igafile)//".dat"
      write(*,*) filename
      call ga_read_inputfile(filename,p)
      
      ntab         = (size(p%gax)-3)/2
      p%pulse%ntab = ntab
      allocate(gaxacc(p%pulse%ntab))
      gaxacc(1) = p%gax(1)
      do i=2,p%pulse%ntab
         gaxacc(i) = gaxacc(i-1)+p%gax(i)
         elaser    = elaser+0.5*p%gax(i)*(p%gax(i-1+ntab)+p%gax(i+ntab))
      end do

      !The following lines are for the pulse shape
      do i=1,p%pulse%ntab
         p%pulse%ttab(i)= gaxacc(i)
         p%pulse%ptab(i)= p%gax(i+p%pulse%ntab)
      end do

      !JingyunYang 20230323 for laser plateau
      p%pulse%ptab(p%pulse%ntab-1)=p%gax(2*p%pulse%ntab-3)
      p%pulse%ptab(p%pulse%ntab-2)=p%gax(2*p%pulse%ntab-3)

      !The following lines are for the target thickness
      l%thick(1)= p%gax(2*ntab+1)
      l%thick(2)= p%gax(2*ntab+2)
  
      ! check input laser energy, gapulse
      !if(elaser>27*5.58 .or. elaser<1*5.6) then
        !p%rhorDT    = 0
        !p%alphaDT   = 0
        !p%timean    = 0
        !p%vmean     = 0
        !p%timestag  = 0
        !p%ifarmax   = 0
        !p%rhomax    = 0
        !p%pmaxdt    = 0
        !p%Edelas    = elaser*1e10
        !p%Ekimp     = 0
        !call ga_print_fitness(p)
        !stop
      !end if

!-----fill in spatial variables
      p%n      = 0
      s%x(1)   = p%xmin
      do i=1,l%nl
         mid        = 0
         do j=1,p%nmat
            if(l%material(i)==p%mat(j)%name) mid=j
         end do
         if(mid==0)then
            print *,'ERROR in subroutine init_state'
            print *,'      unknown material ',l%material(i)
            stop 
         end if
         call eosenergy(p%mat(mid)%eeos,l%r0(i),l%te0(i),ee)
         call eosenergy(p%mat(mid)%ieos,l%r0(i),l%ti0(i),ei)
         do j=p%n+1,p%n+l%nc(i)
            s%r(j)     = l%r0(i)         ! density
            s%ee(j)    = ee              ! elec. spec. energy
            s%ei(j)    = ei              ! elec. spec. energy
            p%ai(j)    = p%mat(mid)%ai   ! atomic mass
            s%f(j)     = 0               ! fuel fraction
            s%ea(j)    = 0               ! alpha energy density
            p%mid(j)   = mid             ! material code
         end do
         call zoning(l%nc(i),l%zonpar(i),l%thick(i), &
                     s%x(p%n+1:p%n+l%nc(i)+1))
         p%n   = p%n+l%nc(i)
      end do
      s%v              = 0               ! velocity
      s%f(1:p%nfuel)   = 0.5             ! T fraction in fuel
!-----compute mass of cells
      select case (p%igeo)
      case (1)
         p%dm=(s%x(2:n+1)-s%x(1:n))*s%r
      case (2)
         p%dm=cpi*(s%x(2:n+1)**2-s%x(1:n)**2)*s%r
      case (3)
         p%dm=4./3.*cpi*(s%x(2:n+1)**3-s%x(1:n)**3)*s%r
      end select
!-----set the initial input energy partition (equal parts)
      call wctrl(p%ng+1,s%wa,wb,0.d0,0,p%iradia)
!-----set the initial time
      s%time           = 0
      end subroutine init_state
!=======================================================================
!-----wfyuan, 2020/10/27, read variables from a input file
      subroutine ga_read_inputfile(filename,p)
      implicit none
      type(typesimparam),  intent(inout) :: p
      character(len=32)                  :: filename
      integer                            :: nline, ios, i

      ! open the input file 
      open(4,file=filename,status='old')

      ! get the number of lines in this file
      nline=0
      do 
          read(4,*,iostat=ios)
          if (ios/=0) exit 
          nline=nline+1          
      end do

      ! allocate space
      allocate(p%gax(nline))
      rewind(4)

      ! read the variables and put them into p%gax(:)
      do i=1,nline
          read(4,*) p%gax(i)
      end do
      close(4)

      end subroutine ga_read_inputfile
  
!=======================================================================
!-----wfyuan, 2020/10/27, print fitness for ga	  
      subroutine ga_print_fitness(p)
      implicit none
      type(typesimparam),  intent(inout) :: p
      character(len=32)                  :: filename

      filename   = "fit_"//trim(p%igafile)//".dat"
      open(4,file=filename,status='new')

      !if(p%timestag*1e9<0.93*p%pulse%ttab(p%pulse%ntab)) p%rhorDT = 0

      write(4,100) p%rhorDT, p%timean*1e-3,p%vmean*1e-5, p%ifarmax, &
                   p%alphaDT,p%rhomax,     p%pmaxdt*1e-15,          &
                   p%Ekimp/5.6*1e-10,      p%Edelas/5.6*1e-10
100   Format(9(F9.3))
      close(4)
      end subroutine ga_print_fitness
  
!=======================================================================
!-----routine to initializate coordinates
!-----every cell is a factor 'zonpar' greater than their left neighbour,
!-----for negative 'zonpar' abs(zonpar) and 1/abs(zonpar) are applied
!-----to half intervals
      subroutine zoning(n,zonpar,delta,x)
      implicit none
      integer, intent(in)    :: n
      real*8,  intent(in)    :: zonpar,delta
      real*8,  intent(inout) :: x(n+1)
      real*8                 :: fact,dx(n)
      integer                :: i
      if(zonpar>0)then
         dx(1)=1 
         do i=2,n
            dx(i)=abs(zonpar)*dx(i-1)
         end do
      else
         dx(1)=1 
         do i=2,n/2
            dx(i)=abs(zonpar)*dx(i-1)
         end do
         dx(n)=1 
         do i=n-1,n/2+1,-1
            dx(i)=abs(zonpar)*dx(i+1)
         end do
      end if
      fact=delta/sum(dx)
      do i=1,n
         x(i+1)=x(i)+fact*dx(i)
      end do
      end subroutine zoning
!======================================================================
!-----obtain group definition from loaded tables
      subroutine frequencies(p)
      implicit none
      type(typesimparam),  intent(inout) :: p
      type(multitable)                   :: multitab
      integer, parameter                 :: nfrei=10000
      real*8                             :: frei(nfrei),fa,fb
      integer                            :: ng,i,j,k
!-----frequency boundaries of all stored tables
      ng       = 0
      do i=1,p%nmat
         do j=1,3
            select case (j)
            case (1)
               multitab= p%mat(i)%planck
            case (2)
               multitab= p%mat(i)%ross
            case (3)
               multitab= p%mat(i)%eps
            end select
            do k=1,multitab%ng
               fa      = multitab%tables(k)%fmin
               fb      = multitab%tables(k)%fmax
               if(fa/=fb)then
                  call glist(nfrei,ng,frei,fa)
                  call glist(nfrei,ng,frei,fb)
               end if
            end do
         end do
      end do
!-----no multigroup tables found, use a wide spectrum
      if(ng==0)then
         ng            = 2
         frei(1)       = 1
         frei(2)       = 1000000
      end if
!-----store data
      if(associated(p%frei))then
         deallocate(p%frei)
         memory        = memory-1
      end if
      allocate(p%frei(ng))
      memory   = memory+1
      p%frei   = frei(1:ng)
      p%ng     = ng-1
      end subroutine frequencies
!=======================================================================
!-----insert xele in sorted list x(n)
      subroutine glist(nmax,n,x,xele)
      implicit none
      integer, intent(in)    :: nmax
      integer, intent(inout) :: n
      real*8,  intent(inout) :: x(nmax)
      real*8,  intent(in)    :: xele
      integer                :: i
      do i=1,n
         if(x(i)==xele)return
      end do
      if(n==nmax)then
         print *,'ERROR in subroutine glist'
         print *,'      n==nmax (',nmax,')'
         stop
      end if
      do i=1,n
         if(xele<x(i))then
            x(i+1:n+1) = x(i:n)
            x(i)       = xele
            n          = n+1
            return
         end if
      end do
      x(n+1)           = xele
      n                = n+1
      end subroutine glist
!======================================================================!
!                                                                      !
!     output routines                                                  !
!                                                                      !
!======================================================================!
!-----write ascii output
      subroutine ascii_output(p,s,e)
      implicit none
      type(typesimparam),  intent(in) :: p
      type(typesimstate),  intent(in) :: s
      type(typesimenergy), intent(in) :: e
      integer                         :: i
      type(typesimthermo)             :: t
      call allocsimthermo(p%n,t)
      call eosall(p,s,t)
      write(11,*)repeat('=',78)
      write(11,*)' '
      write(11,*)'subroutine  ascii_ouput'
      write(11,*)' '
      write(11,'(a,i5,a,e14.6)')' step=',e%istep,' time=',s%time
      write(11,*)' '
      write(11,'(a4,1x,5a14)')'i','x','v','rho','te','depo'
      write(11,*)' '
      do i=1,p%n
         write(11,'(i4,1x,5e14.6)')i,s%x(i),s%v(i),s%r(i), &
         t%te(i),e%d(i)*s%r(i)
      end do
      write(11,'(i4,1x,5e14.6)')p%n+1,s%x(p%n+1),s%v(p%n+1)
      write(11,*)' '
      call freesimthermo(t)
      end subroutine ascii_output
!=======================================================================
!-----wfyuan, 2020/10/27, calculate fitness for genetic algorithm
      subroutine ga_calculate_fitness(p,s,e)
      implicit none
      type(typesimparam)              :: p
      type(typesimstate),  intent(in) :: s
      type(typesimenergy), intent(in) :: e
      type(typesimthermo)             :: t
      real*8                          :: vimplo,rhorDT,rhorHS
      real*8                          :: Timean
      real*8                          :: ifar
      real*8                          :: alphaDT
      real*8                          :: rhomax
      real*8                          :: pmaxdt,ptemp
      real*8                          :: Ekimp
      real*8, dimension(p%n)          :: vcell
      integer                         :: i,nfuel
  
      !-----get thermodynamic quantities
      call allocsimthermo(p%n,t)
      call eosall(p,s,t)
  
      call compute_vimplo(p,s,vimplo)
      call compute_rhor(p,s,t,0.0d0,rhorDT)
      call compute_ifar(p,s,ifar)
      call compute_alpha(p,s,t,alphaDT)

      nfuel  = p%nfuel
      rhomax = 0;
      do i=61,nfuel
         if(s%r(i)>rhomax)then
            rhomax = s%r(i)
         end if

          ptemp=t%pe(i)+t%pi(i)
         if(ptemp>pmaxdt)then
            pmaxdt = ptemp
         end if

      end do


      if(p%vmean < vimplo)then
         p%vmean = vimplo
         !p%alphaDT  = alphaDT
      end if

      ifar = s%x(161)/(s%x(161)-s%x(61))
      if(p%ifarmax<ifar)then 
          p%ifarmax = ifar
          !p%alphaDT  = alphaDT
      end if
  
      if(p%rhorDT<rhorDT)then
         p%rhorDT   = rhorDT
         p%pmaxdt   = pmaxdt
         p%timestag = s%time
         p%alphaDT  = alphaDT
      end if

      if(p%rhomax<rhomax)then
         p%rhomax   = rhomax
      end if

      !Timean = sum(t%ti(p%nfuel/2+1:p%nfuel))/p%nfuel*2
      Timean = sum(t%ti(61:nfuel))/(nfuel-60)
      if(p%timean < Timean) p%timean = Timean
  
      ! calculate implosion kinentic energy
      vcell(1:p%n) = (s%v(1:p%n)+s%v(2:p%n+1))/2
      where(vcell>0) vcell = 0
      Ekimp        = sum(0.5*p%dm(1:p%n)*vcell(1:p%n)**2)
      if(p%Ekimp < Ekimp) p%Ekimp = Ekimp

      p%Edelas = e%delas

      call freesimthermo(t)
  
      end subroutine ga_calculate_fitness
!=======================================================================
!-----write binary output record
      subroutine write_output_record(p,s,e)
      implicit none
      type(typesimparam),  intent(in) :: p
      type(typesimstate),  intent(in) :: s
      type(typesimenergy), intent(in) :: e
      integer,                   save :: mode=0
      integer                         :: i,i1,icount,n,ng,nr
      real*8                          :: astep,vimplo,rhorDT,rhorHS,ifar
      real*8                          :: alpha
      real*8, dimension(p%n)          :: cmc,nc,xc,dene,ad,pt
      real*8, dimension(p%n+1)        :: ni,cmi,tr,u,sp,sm,a,sr
      real*8, dimension(p%ng)         :: frec,radl,radr
      type(typesimthermo)             :: t
      n        = p%n
      ng       = p%ng
      nr       = p%nreduce
!-----get thermodynamic quantities
      call allocsimthermo(n,t)
      call eosall(p,s,t)
!-----compute auxiliar quantities
      if(mode==0)then
         forall(i=1:n)   nc(i)=i
         forall(i=1:n+1) ni(i)=i
         cmi(1)=0
         do i=2,n+1
            cmi(i)=cmi(i-1)+p%dm(i-1)
         end do
         cmc   = 0.5*(cmi(1:n)+cmi(2:n+1))
         frec  = 0.5*(p%frei(2:ng+1)+p%frei(1:ng))
      end if
      astep    = e%istep
      xc       = 0.5*(s%x(1:n)+s%x(2:n+1))
      dene     = s%r*t%zi/p%ai/pm
      ad       = e%d*s%r
      select case(p%igeo)
      case(1)
         a     = 1
      case(2)
         a     = 2*cpi*s%x
      case(3)
         a     = 4*cpi*s%x**2
      end select
      call compute_radia(n,ng,e%strahl,tr,u,sp,sm)
      sp       = sp*a    
      sm       = sm*a    
      sr       = sp-sm
      pt       = t%pe+t%pi+2.0d0/3.0d0*s%ea
      call compute_spectra(n,ng,e%strahl,a,p%frei,p%irad_left,radl,2)
      call compute_spectra(n,ng,e%strahl,a,p%frei,p%irad_right,radr,1)
      call compute_vimplo(p,s,vimplo)
      call compute_rhor(p,s,t,0.0d0,rhorDT)
      call compute_rhor(p,s,t,5.0d3,rhorHS)
      call compute_ifar(p,s,ifar)
      call compute_alpha(p,s,t,alpha)
!-----first section
      if(mode==0)then
         do i=0,3
             icount=0
             call write_vect(i,icount,nr,'CMC     ',n,cmc)
             call write_vect(i,icount,nr,'NC      ',n,nc)
             call write_vect(i,icount,-nr,'CMI     ',n+1,cmi)
             call write_vect(i,icount,-nr,'NI      ',n+1,ni)
             call write_vect(i,icount,0,'FREC    ',ng,frec)
             call write_vect(i,icount,0,'FREI    ',ng+1,p%frei)
             if(i==0)write(10,'(i6)')icount
         end do
      end if
!-----second section
      i1=0
      if(mode==1)i1=3
      do i=i1,3
         icount=0
         call write_item(i,icount,'TIME    ',s%time)
         call write_item(i,icount,'ISTEP   ',astep)
         call write_item(i,icount,'LASER   ',e%laser)
         call write_item(i,icount,'LPOWER  ',e%power)
         call write_item(i,icount,'ALPHAS  ',e%alphas)
         call write_item(i,icount,'FUSION  ',e%fusion)
         call write_item(i,icount,'NEUTRONS',e%neutrons)
         call write_item(i,icount,'ENIA    ',e%enia)
         call write_item(i,icount,'ENIE    ',e%enie)
         call write_item(i,icount,'ENII    ',e%enii)
         call write_item(i,icount,'ENKI    ',e%enki)
         call write_item(i,icount,'DELAS   ',e%delas)
         call write_item(i,icount,'DERAD   ',e%derad)
         call write_item(i,icount,'DEFUS   ',e%defus)
         call write_item(i,icount,'ENRL1   ',e%radoutl)
         call write_item(i,icount,'ENRL4   ',e%radoutr)
         call write_item(i,icount,'ENRG1   ',e%radinl)
         call write_item(i,icount,'ENRG4   ',e%radinr)
         call write_item(i,icount,'REF     ',e%ref)
         call write_item(i,icount,'TRA     ',e%tra)
         call write_item(i,icount,'VIMPLO  ',vimplo)
         call write_item(i,icount,'RHORDT  ',rhorDT)
         call write_item(i,icount,'RHORHS  ',rhorHS)
         call write_item(i,icount,'IFAR    ',ifar)
         call write_item(i,icount,'ALPHA   ',alpha)
         call write_vect(i,icount,nr,'R       ',n,s%r)
         call write_vect(i,icount,nr,'PT      ',n,pt)
         call write_vect(i,icount,nr,'TE      ',n,t%te)
         call write_vect(i,icount,nr,'PE      ',n,t%pe)
         call write_vect(i,icount,nr,'ZI      ',n,t%zi)
         call write_vect(i,icount,nr,'XC      ',n,xc)
         call write_vect(i,icount,nr,'D       ',n,e%d)
         call write_vect(i,icount,nr,'DENE    ',n,dene)
         call write_vect(i,icount,nr,'EA      ',n,s%ea)
         call write_vect(i,icount,nr,'F       ',n,s%f)
         call write_vect(i,icount,nr,'TI      ',n,t%ti)
         call write_vect(i,icount,nr,'PI      ',n,t%pi)
         call write_vect(i,icount,nr,'AD      ',n,ad)
         call write_vect(i,icount,-nr,'X       ',n+1,s%x)
         call write_vect(i,icount,-nr,'V       ',n+1,s%v)
         call write_vect(i,icount,-nr,'TR      ',n+1,tr)
         call write_vect(i,icount,-nr,'S+      ',n+1,sp)
         call write_vect(i,icount,-nr,'S-      ',n+1,sm)
         call write_vect(i,icount,-nr,'S       ',n+1,sr)
         call write_vect(i,icount,0,'I-      ',ng,radl)
         call write_vect(i,icount,0,'I+      ',ng,radr)
         if(i==0)write(10,'(i6)')icount
      end do
      mode=1
      call freesimthermo(t)
      end subroutine write_output_record
!=======================================================================
!-----write binary record fragment for a scalar
      subroutine write_item(mode,icount,label,value)
      implicit none
      integer,              intent(in)    :: mode
      character(len=8),     intent(in)    :: label
      integer,              intent(inout) :: icount
      real*8,               intent(in)    :: value
      real*8                              :: vector(1)
      vector(1)= value
      call write_vect(mode,icount,0,label,1,vector)
      end subroutine write_item
!=======================================================================
!-----write binary record fragment for a vector
      subroutine write_vect(mode,icount,nr,label,n,value)
      implicit none
      integer,              intent(in)    :: mode,n,nr
      character(len=8),     intent(in)    :: label
      integer,              intent(inout) :: icount
      real*8, dimension(n), intent(in)    :: value
      integer                             :: i,i1,i2,m
!-----number of values computed as function of n and nr
      if(nr==0)then                           ! all values
         m=n
      else if(nr>0) then                      ! cell values
         m=(n+nr-1)/nr
      else                                    ! interface values
         m=(n-nr-2)/(-nr)+1
      end if
      select case(mode)
      case(0)                                 ! count number of elements
         icount=icount+m
      case(1)                                 ! name of variables
         write(10,'(a)')(label,i=1,m)
      case(2)                                 ! indices 
         if(n==1) then
            write(10,'(i6)')0
         else
            write(10,'(i6)')(i,i=1,m)   
         end if
      case(3)
         if(nr==0)then                        ! all values
            write(10,'(e12.5)')(value(i),i=1,n)
         else if(nr>0)then                    ! cell values
            do i=1,m
               i1=(i-1)*nr+1
               i2=min(n,i1+nr-1)
               write(10,'(e12.5)')sum(value(i1:i2))/(i2-i1+1)
            end do
         else                                 ! interface values
            do i=1,m
               i1=min((i-1)*(-nr)+1,n)
               write(10,'(e12.5)')value(i1)
            end do
         end if         
      end select
      end subroutine write_vect
!=======================================================================
      subroutine compute_spectra(n,ng,strahl,area,frei,i,s,idir)
      implicit none
      integer,                intent(in)  :: n,ng,i,idir
      real*8,                 intent(in)  :: strahl(n+1,ng,2)
      real*8,                 intent(in)  :: area(n+1),frei(ng+1)
      real*8,                 intent(out) :: s(ng)
      integer                             :: ii
      real*8                              :: dfre(ng)
      ii       = max(1,min(i,n+1))
      dfre     = frei(2:ng+1)-frei(1:ng)
      if(idir==1)then
         s        = strahl(ii,:,1)/(cpi*dfre)*area(ii)
      else
         s        = strahl(ii,:,2)/(cpi*dfre)*area(ii)
      end if
      end subroutine compute_spectra
!=======================================================================
!-----compute relevant radiation quantities
!-----n        - number of cells
!-----ng       - number of groups
!-----strahl   - radiation field: strahl(n+1,ng,2)
!-----           strahl(i,j,1)=S+ at interface i for group j
!-----           strahl(i,j,2)=S- at interface i for group j
!-----tr       - radiation temperature
!-----u        - radiation density of energy
!-----sp       - S+
!-----sm       - S-
      subroutine compute_radia(n,ng,strahl,tr,u,sp,sm)
      implicit none
      integer,                intent(in)  :: n,ng
      real*8,                 intent(in)  :: strahl(n+1,ng,2)
      real*8, dimension(n+1), intent(out) :: tr,u,sp,sm
      integer                             :: i
      sp       = 0
      sm       = 0
      do i=1,ng
         sp    = sp+strahl(:,i,1)
         sm    = sm+strahl(:,i,2)
      end do
      u        = 2*(sp+sm)/c
      tr       = (0.25*c*abs(u)/sigma)**0.25
      end subroutine compute_radia
!=======================================================================
!-----compute implosion velocity of DT fuel
      subroutine compute_vimplo(p,s,vimplo)
      implicit none
      type(typesimparam), intent(in)  :: p
      type(typesimstate), intent(in)  :: s
      real*8,             intent(out) :: vimplo
      integer                         :: nf
      real*8                          :: ms,m(p%nfuel+1)
      nf               = p%nfuel
      m                = 0
      m(1:nf)          = p%dm(1:nf) 
      m(2:nf+1)        = m(2:nf+1)+p%dm(1:nf)
      where(s%v(1:nf+1)>0)  m = 0
      ms               = sum(m)
      if(ms>0)then
         vimplo        = -sum(m*s%v(1:nf+1))/ms
      else
         vimplo        = 0
      end if
      end subroutine compute_vimplo
!=======================================================================
!-----compute rhoR
      subroutine compute_rhor(p,s,t,tmin,rhor)
      implicit none
      type(typesimparam),  intent(in)  :: p
      type(typesimstate),  intent(in)  :: s
      type(typesimthermo), intent(in)  :: t
      real*8,              intent(in)  :: tmin
      real*8,              intent(out) :: rhor
      integer                          :: nf
      real*8                           :: rhodr(p%nfuel)
      nf               = p%nfuel
      rhodr            = s%r(1:nf)*(s%x(2:nf+1)-s%x(1:nf))
      where(t%ti(1:nf)<tmin)  rhodr = 0
      rhor             = sum(rhodr)
      end subroutine compute_rhor
!=======================================================================
!-----compute ifar (in flight aspect ratio)
      subroutine compute_ifar(p,s,ifar)
      implicit none
      type(typesimparam), intent(in)  :: p
      type(typesimstate), intent(in)  :: s
      real*8,             intent(out) :: ifar
      real*8,             parameter   :: factor=0.3678794   ! 1/e
      real*8                          :: rho(p%n+1),f,rhomax,xmin,xmax
      integer                         :: i,n,imax
      n                = p%n
      rho(1)           = 0
      rho(2:n)         = s%r(1:n-1)+s%r(2:n)
      rho(n+1)         = 0
      rhomax           = 0
      do i=2,n
         if(rho(i)>rhomax)then
            imax       = i 
            rhomax     = rho(i)
         end if
      end do
      rhomax           = rhomax*factor
      do i=imax,n
         if((rho(i)-rhomax)*(rho(i+1)-rhomax)<0)then
            f          = (rhomax-rho(i+1))/(rho(i)-rho(i+1))
            xmax       = s%x(i+1)*(1-f)+s%x(i)*f
            exit
         end if
      end do
      do i=imax,2,-1
         if((rho(i)-rhomax)*(rho(i-1)-rhomax)<0)then
            f          = (rhomax-rho(i-1))/(rho(i)-rho(i-1))
            xmin       = s%x(i-1)*(1-f)+s%x(i)*f
            exit
         end if
      end do
      ifar             = (xmax+xmin)/(xmax-xmin)/2
      end subroutine compute_ifar
!=======================================================================
!-----compute average of fuel pressure/fermi pressure
      subroutine compute_alpha(p,s,t,alpha)
      implicit none
      type(typesimparam),  intent(in)  :: p
      type(typesimstate),  intent(in)  :: s
      type(typesimthermo), intent(in)  :: t
      real*8,              intent(out) :: alpha
      real*8,               parameter  :: factor=0.1
      real*8                           :: dene,fermi,rcut,mass
      integer                          :: i,nfuel
      alpha    = 0
      nfuel    = p%nfuel
      if(nfuel==0)return
      rcut     = factor*maxval(s%r(1:nfuel))
      mass     = 0
      do i=1,nfuel
          if(s%r(i)>rcut)then
             dene  = s%r(i)*t%zi(i)/p%ai(i)/pm
             fermi = hb**2/(5*em)*(3*cpi**2*dene)**(2d0/3d0)*dene
             mass  = mass+p%dm(i)
             alpha = alpha+p%dm(i)*(t%pe(i)+t%pi(i))/fermi
          end if
      end do
      alpha    = alpha/mass
      end subroutine compute_alpha
!=======================================================================
!-----print heading
      subroutine print_heading
      implicit none
      integer                           :: i
      character (len=5*8), dimension(9) :: a1
      character (len=4*8), dimension(9) :: a2
      a1(1)    = 'M    M  U    U  LL      TTTTTT  IIIIII'
      a1(2)    = 'MM  MM  U    U  LL      TTTTTT  IIIIII'
      a1(3)    = 'MMMMMM  U    U  LL        TT      II'
      a1(4)    = 'M MM M  U    U  LL        TT      II'
      a1(5)    = 'M    M  U    U  LL        TT      II'
      a1(6)    = 'M    M  U    U  LL        TT      II'
      a1(7)    = 'M    M  U    U  LL        TT      II'
      a1(8)    = 'M    M  UUUUUU  LLLLLL    TT    IIIIII'
      a1(9)    = 'M    M   UUUU   LLLLLL    TT    IIIIII'
      a2(1)    = '        IIIIII  FFFFFF  EEEEEE'
      a2(2)    = '        IIIIII  FFFFFF  EEEEEE'
      a2(3)    = '          II    F       E'
      a2(4)    = '          II    F       E'
      a2(5)    = '======    II    FFFFF   EEEEE'
      a2(6)    = '          II    F       E'
      a2(7)    = '          II    F       E'
      a2(8)    = '        IIIIII  F       EEEEEE'
      a2(9)    = '        IIIIII  F       EEEEEE'
      print *,repeat('=',78)
      print *,' '
      print *,' '
      do i=1,size(a1)
         print '(4x,a)',a1(i)//a2(i)
      end do
      print *,' '
      print *,' '
      print *,'MULTI-IFE - A One-Dimensional Computer Code for'
      print *,'Inertial Fusion Energy (IFE) Target Simulations'
      print *,' '
      print *,'R.Ramis and J. Meyer-ter-Vehn'
      print *,' '
      print *,'Version March 4, 2015'
      print *,' '
      end subroutine print_heading
!=======================================================================
!-----print tail
      subroutine print_tail
      print *,repeat('=',78)
      print *,' '
      print *,'normal program termination'
      print *,' '
      print *,'used memory check'
      print *,' ' 
      call pri('memory',  memory,    'allocations - deallocations')
      print *,' ' 
      print *,repeat('=',78)
      end subroutine print_tail
!=======================================================================
!-----print head of integration section
      subroutine print_int_head
      print *,repeat('=',78)
      print *,' ' 
      print *,'subroutine integrate'
      print *,' ' 
      print '(a)','   step      time          dt        max-var'
      print *,' '
      end subroutine print_int_head
!=======================================================================
!-----print line of integration section
      subroutine print_int_line(step,time,dt,var,mode)
      integer, intent(in) :: step,mode
      real*8,  intent(in) :: time,dt,var
      select case(mode)
      case(-1)
          print '(i6,2e13.5,a)',step,time,dt,' neg. density'
      case(-2)
          print '(i6,3e13.5,a)',step,time,dt,var,' var fail'
      case(1)
          print '(i6,3e13.5,a)',step,time,dt,var
      end select
      end subroutine print_int_line
!=======================================================================
!-----print tail of integration section
      subroutine print_int_tail
      print *,' '
      end subroutine print_int_tail
!=======================================================================
!-----print significant values in structure 'multitab'
      subroutine print_multitable(multitab,ind)
      implicit none
      type(multitable),   intent(in) :: multitab
      integer,            intent(in) :: ind
      integer                        :: nr,nt,ng,i
      type(grouptable)               :: grouptab
      real*8, dimension(multitab%ng) :: xt,yt
      real*8                         :: rmin,rmax,tmin,tmax,fmin,fmax
      ng       = multitab%ng
      call prc('ng',ng,0,'table not loaded')
      if(ng==0)return
      grouptab = multitab%tables(1)
      nr       = grouptab%nr  
      nt       = grouptab%nt  
      do i=1,ng
         xt(i)    = multitab%tables(i)%fmin
         yt(i)    = multitab%tables(i)%fmax
      end do
      rmin     = 10**minval(grouptab%rho)
      rmax     = 10**maxval(grouptab%rho)
      tmin     = 10**minval(grouptab%t)
      tmax     = 10**maxval(grouptab%t)
      fmin     = minval(xt)
      fmax     = maxval(yt)
      call pri('nr'  ,nr,  'number of tabulated densities')
      call prf('rmin',rmin,'min. density')
      call prf('rmax',rmax,'max. density')
      call pri('nt'  ,nt,  'number of tabulated temperatures')
      call prf('tmin',tmin,'min. temperature')
      call prf('tmax',tmax,'max. temperature')
      call pri('ng'  ,ng,  'number of groups')
      call prf('fmin',fmin,'min. frequency')
      call prf('fmax',fmax,'max. frequency')
      if(ind==1)then
         print *,' '
         call prt(ng,'fmin','fmax',xt,yt)
      end if
      end subroutine print_multitable
!======================================================================
!-----print significant values in structure 'eostab'
      subroutine print_eostable(eostab)
      implicit none
      type(eostable),   intent(in) :: eostab
      integer                      :: nr,ne
      real*8                       :: rmin,rmax,tmin,tmax,pmin,pmax
      nr       = eostab%nr   
      ne       = eostab%ne   
      call prc('nr',nr,0,'table not loaded')
      if(nr==0)return
      rmin     = minval(eostab%rho)
      rmax     = maxval(eostab%rho)
      tmin     = maxval(eostab%t(1:nr))
      tmax     = minval(eostab%t(ne*nr-nr+1:nr*ne))
      pmin     = minval(eostab%p)
      pmax     = maxval(eostab%p)
      call pri('nr'  ,nr,  'number of tabulated densities')
      call pri('ne'  ,ne,  'number of tabulated energies')
      call prf('rmin',rmin,'min. density')
      call prf('rmax',rmax,'max. density')
      call prf('tmin',tmin,'min. temperature')
      call prf('tmax',tmax,'max. temperature')
      call prf('pmin',pmin,'min. pressure')
      call prf('pmax',pmax,'max. pressure')
      end subroutine print_eostable
!======================================================================
!-----print significant values in structure 'grouptab'
      subroutine print_grouptable(grouptab)
      implicit none
      type(grouptable), intent(in) :: grouptab
      integer                      :: nr,nt
      real*8                       :: rmin,rmax,tmin,tmax,omin,omax
      nr       = grouptab%nr   
      nt       = grouptab%nt   
      call prc('nr',nr,0,'table not loaded')
      if(nr==0)return
      rmin     = 10**minval(grouptab%rho)
      rmax     = 10**maxval(grouptab%rho)
      tmin     = 10**minval(grouptab%t)
      tmax     = 10**maxval(grouptab%t)
      omin     = 10**minval(grouptab%o)
      omax     = 10**maxval(grouptab%o)
      call pri('nr'  ,nr,  'number of tabulated densities')
      call prf('rmin',rmin,'min. density')
      call prf('rmax',rmax,'max. density')
      call pri('nt'  ,nt,  'number of tabulated temperatures')
      call prf('tmin',tmin,'min. temperature')
      call prf('tmax',tmax,'max. temperature')
      call prf('omin',omin,'min. value')
      call prf('omax',omax,'max. value')
      end subroutine print_grouptable
!=======================================================================
!-----print significant values in structure 'typematerial'
      subroutine print_material(m,ind)
      implicit none
      type(typematerial), intent(in) :: m
      integer,            intent(in) :: ind
      print *,repeat('=',78)
      print *,' '
      print *,'subroutine print_material'
      print *,' '
      call pra('name',       m%name,       'material name')
      call prf('zi',         m%zi,         'atomic number')
      call prf('ai',         m%ai,         'atomic weight')
      print *,' '
      call prn('eeos_file',  m%eeos_file,  'electron EOS file')
      call pri('eeos_id',    m%eeos_id,    'electron EOS id')
      call print_eostable(m%eeos)
      print *,' '
      call prn('ieos_file',  m%ieos_file,  'ion EOS file')
      call pri('ieos_id',    m%ieos_id,    'ion EOS id')
      call print_eostable(m%ieos)
      print *,' '
      call prn('z_file',     m%z_file,     'Z file')
      call pri('z_id',       m%z_id,       'Z id')
      call print_grouptable(m%z)
      print *,' '
      call prn('planck_file',m%planck_file,'Planck file')
      call pri('planck_id',  m%planck_id,  'Planck id')
      call print_multitable(m%planck,ind)
      print *,' '
      call prn('ross_file',  m%ross_file,  'Rosseland file')
      call pri('ross_id',    m%ross_id,    'Rosseland id')
      call print_multitable(m%ross,ind)
      print *,' '
      call prn('eps_file',   m%eps_file,   'Emissivity file')
      call pri('eps_id',     m%eps_id,     'Emissivity id')
      call print_multitable(m%eps,ind)
      print *,' '
      end subroutine print_material
!=======================================================================
!-----print significant values in structure 'typesimparam'
      subroutine print_parameters(p)
      implicit none
      type(typesimparam), intent(in) :: p
      print *,repeat('=',78)
      print *,' '
      print *,'subroutine print_parameters'
      print *,' '
      print *,'problem definition (general):'
      print *,' '
      call prc('igeo',      p%igeo,1,'planar geometry')
      call prc('igeo',      p%igeo,2,'cylindrical geometry')
      call prc('igeo',      p%igeo,3,'spherical geometry')
      call prf('xmin',      p%xmin,'initial value x at left bound.')
      call prf('texit',     p%texit,'end of simulation')
      call pri('nfuel',     p%nfuel,'cells with DT mixture')
      print *,' '
      print *,'problem definition (boundary conditions):'
      print *,' '
      call prc('iright',    p%iright,0,'free expansion at RHS')
      call prc('iright',    p%iright,1,'rigid wall at RHS')
      call prc('ileft',     p%ileft,0,'free expansion at LHS')
      call prc('ileft',     p%ileft,1,'rigid wall at LHS')
      call prf('alphal',    p%alphal,'reflected rad. frac. at L')
      call prf('alphar',    p%alphar,'reflected rad. frac. at R')
      call prf('betal',     p%betal,'rad. pulse coupling at L')
      call prf('betar',     p%betar,'rad. pulse coupling at R')
      call prc('itype',     p%itype,1,'rad. pulse shape: sin**2')
      call prc('itype',     p%itype,2,'rad. pulse shape: constant')
      call prc('itype',     p%itype,4,'rad. pulse shape: tabulated')
      call prf('tau',       p%tau,'rad. pulse duration')
      call prf('trad',      p%trad,'rad. pulse peak temperature')
      print *,' '
      print *,'problem definition (modeling options):'
      print *,' '
      call prc('iradia',    p%iradia,0,'disable radiation')
      call prc('iradia',    p%iradia,1,'enable radiation')
      call prc('ihydro',    p%ihydro,0,'disable hydrodynamics')
      call prc('ihydro',    p%ihydro,1,'enable hydrodynamics')
      call prc('iheation',  p%iheation,0,'disable ion heat conduction')
      call prc('iheation',  p%iheation,1,'enable ion heat conduction')
      call prc('model',     p%model,0,'classical electron collisions')
      call prc('model',     p%model,1,'electron-phonon collisions')
      call prc('model',     p%model,2,'Drude-Sommerfeld collisions')
      if(p%model==1)then
         call prf('fheat',  p%fheat,'heat cond. factor')
         call prf('fei',    p%fei,'e-i interchange factor')
         call prf('flaser', p%flaser,'bremsstrahlung factor')
      end if
      call prf('zmin',      p%zmin,'minimum ionization')
      call prf('flf',       p%flf,'flux limit factor')
      call prc('inegpre',   p%inegpre,0,'force positive Pe')
      call prc('inegpre',   p%inegpre,1,'allow negative Pe')
      print *,' '
      print *,'numerical options:'
      print *,' '
      call pri('n',         p%n,'number of cells')
      call pri('ng',        p%ng,'number of groups')
      call pri('nsplit',    p%nsplit,'requested num. of subcycles')
      call pri('msplit',    p%msplit,'actual num. of subcycles')
      call prf('dtmin',     p%dtmin,'minimum dt')
      call prf('dtmax',     p%dtmax,'maximum dt')
      call prf('dtinit',    p%dtinit,'initial dt')
      call prf('dtrvar',    p%dtrvar,'nominal density var.')
      call prf('dttevar',   p%dttevar,'nominal elec. temp. var.')
      call prf('dttivar',   p%dttivar,'nominal ion temp. var.')
      call prf('dtbreak',   p%dtbreak,'break step factor')
      call prf('dtfactor',  p%dtfactor,'dt increment factor')
      call pri('nexit',     p%nexit,'maximum number of steps')
      call prc('iwctrl',    p%iwctrl,0,'fixed deposition factors')
      call prc('iwctrl',    p%iwctrl,1,'optimize deposition factors')
      print *,' '
      print *,'group distribution:'
      print *,' '
      call prt(p%ng,'fmin','fmax',p%frei(1:p%ng),p%frei(2:p%ng+1))
      print *,' '
      print *,'output control:'
      print *,' '
      call prf('dt_aout',   p%dt_aout,'ascii output time interval')
      call pri('ns_aout',   p%ns_aout,'ascii output step interval')
      call prf('dt_bout',   p%dt_bout,'bin. output time interval')
      call pri('ns_bout',   p%ns_bout,'bin. output step interval')
      call pri('irad_left', p%irad_left,'left spectra interf.')
      call pri('irad_right',p%irad_right,'right spectra interf.')
      call pri('nreduce',   p%nreduce,'reduce spatial resolution')
      print *,' '
      end subroutine print_parameters
!=======================================================================
!-----print significant values in structure 'typesimlayers'
      subroutine print_layer_definition(l)
      implicit none
      type(typesimlayers), intent(in) :: l
      integer                         :: i
      print *,repeat('=',78)
      print *,' '
      print *,'subroutine print_layer_definition'
      print *,' '
      do i=1,l%nl
         call pri('i',       i,          'layer number')
         call pra('material',l%material(i),'name')
         call prf('thick',   l%thick(i), 'thickness')
         call prf('r0',      l%r0(i),    'initial density')
         call prf('te0',     l%te0(i),   'initial electron temperature')
         call prf('ti0',     l%ti0(i),   'initial ion temperature')
         call pri('nc',      l%nc(i),    'number of cells')
         call prf('zonpar',  l%zonpar(i),'griding parameter')
         print *,' '
      end do
      end subroutine print_layer_definition
!=======================================================================
!-----print significant values in structure 'typesimlaserwkb'
      subroutine print_laser_wkb(p)
      implicit none
      type(typesimlaserwkb), intent(in) :: p
      real*8                            :: freq,dncrt
      if(p%enable==0)return
      freq     = c/p%wl
      dncrt    = cpi*freq**2*em/eq**2
      print *,repeat('=',78)
      print *,' '
      print *,'subroutine print_laser_wkb'
      print *,' '
      call prf('pimax', p%pimax,   'nominal laser power')
      call prf('pitime',p%pitime,  'nominal laser time')
      call prf('wl',    p%wl,      'laser wavelenght')
      call prc('inter', p%inter,1, 'laser comes from RHS')
      call prc('inter', p%inter,-1,'laser comes from LHS')
      call prf('delta', p%delta,   'frac. abs. at crit. dens.')
      call prc('itype', p%itype,1, 'sin**2 shape')
      call prc('itype', p%itype,2, 'rectangular shape')
      call prc('itype', p%itype,4, 'tabulated shape')
      call prf('dncrt', dncrt,     'critical density')
      print *,' '
      end subroutine print_laser_wkb
!=======================================================================
!-----print significant values in structure 'typesimlasermaxwell'
      subroutine print_laser_maxwell(p)
      implicit none
      type(typesimlasermaxwell), intent(in) :: p
      real*8                                :: freq,dncrt
      if(p%enable==0)return
      freq  = c/p%wl
      dncrt = cpi*freq**2*em/eq**2
      print *,repeat('=',78)
      print *,' '
      print *,'subroutine print_laser_maxwell'
      print *,' '
      call prf('pimax', p%pimax,   'nominal laser intensity')
      call prf('pitime',p%pitime,  'nominal laser time')
      call prf('wl',    p%wl,      'laser wavelenght')
      call prc('inter', p%inter,1, 'laser comes from RHS')
      call prc('inter', p%inter,-1,'laser comes from LHS')
      call prc('itype', p%itype,1, 'sin**2 shape')
      call prc('itype', p%itype,2, 'rectangular shape')
      call prc('itype', p%itype,4, 'tabulated shape')
      call pri('idep',  p%idep,    'cell subdivisions')
      call prf('angle', p%angle,   'incidence angle')
      call prc('pol',   p%pol,1,   'p-polarization')
      call prc('pol',   p%pol,2,   's-polarization')
      call prf('dncrt', dncrt,     'critical density')
      print *,' '
      end subroutine print_laser_maxwell
!=======================================================================
!-----print significant values in structure 'typesimlaser3d'
      subroutine print_laser_3d(p)
      implicit none
      type(typesimlaser3d), intent(in) ::  p
      real*8                                :: freq,dncrt
      if(p%enable==0)return
      freq  = c/p%wl
      dncrt = cpi*freq**2*em/eq**2
      print *,repeat('=',78)
      print *,' '
      print *,'subroutine print_laser_3d'
      print *,' '
      call prf('pimax', p%pimax,  'nominal laser power')
      call prf('pitime',p%pitime, 'nominal laser time')
      call prf('wl',    p%wl,     'laser wavelenght')
      call prc('itype', p%itype,1,'sin**2 shape')
      call prc('itype', p%itype,2,'rectangular shape')
      call prc('itype', p%itype,4,'tabulated shape')
      call pri('inr',   p%nr,     'number of rays')
      call prf('rmax',  p%rmax,   'maximum radius')
      call prf('fwhm',  p%fwhm,   'full width at half maximum')
      call prf('bexp',  p%bexp,   'super-Gaussian exponent')
      call prf('dncrt', dncrt,    'critical density')
      print *,' '
      end subroutine print_laser_3d
!=======================================================================
!-----print significant values in structure 'typesimpulse'
      subroutine print_pulse(p)
      implicit none
      type(typesimpulse), intent(in) ::  p
      if(p%ntab==0)return
      print *,repeat('=',78)
      print *,' '
      print *,'subroutine print_pulse'
      print *,' '
      call prc('mode',p%mode,1,'linear interp. in value**expo')
      call prc('mode',p%mode,2,'linear interp. in log(value)')
      call prc('mode',p%mode,3,'linear interp. in exp(value)')
      if(p%mode==1)call prf('expo',p%expo,'exponent')
      print *,' '
      call prt(p%ntab,'time','value',p%ttab,p%ptab)
      print *,' '
      end subroutine print_pulse
!=======================================================================
!-----print energies
      subroutine print_energies(e)
      implicit none
      type(typesimenergy), intent(in) :: e
      real*8                          :: energy(11),epos=0,eneg=0,error
      integer                         :: i
!-----energies incorporated to the system 
      energy(1)        = e%enie0
      energy(2)        = e%enii0
      energy(3)        = e%enia0
      energy(4)        = e%enki0
      energy(5)        = e%delas
      energy(6)        = e%defus
      energy(7)        = e%derad
!-----energies extracted from the system
      energy(8)        = -e%enie
      energy(9)        = -e%enii
      energy(10)       = -e%enki
      energy(11)       = -e%enia
!-----positive and negative energies
      do i=1,11
         if(energy(i)>0)then
             epos      = epos+energy(i)
         else
             eneg      = eneg-energy(i)
         end if
      end do
!-----relative error
      error            = abs(epos-eneg)/((epos+eneg)/2)
      print *,repeat('=',78)
      print *,' ' 
      print *,'subroutine print_energies'
      print *,' ' 
      print *,'integration steps'
      print *,' ' 
      call pri('itry',    e%itry,    'total')
      call pri('istep',   e%istep,   'valid')
      print *,' ' 
      print *,'initial energies'
      print *,' ' 
      call prf('enie0',   e%enie0,   'electron energy')
      call prf('enii0',   e%enii0,   'ion energy')
      call prf('enia0',   e%enia0,   'alpha energy')
      call prf('enki0',   e%enki0,   'kinetic energy')
      print *,' ' 
      print *,'actual energies'
      print *,' ' 
      call prf('enie',    e%enie,    'electron energy')
      call prf('enii',    e%enii,    'ion energy')
      call prf('enia',    e%enia,    'alpha energy')
      call prf('enki',    e%enki,    'kinetic energy')
      print *,' ' 
      print *,'energy deposition'
      print *,' ' 
      call prf('delas',   e%delas,   'laser')
      call prf('derad',   e%derad,   'thermal radiation')
      call prf('defus',   e%defus,   'fusion')
      print *,' ' 
      print *,'radiation input/output'
      print *,' ' 
      call prf('radoutl', e%radoutl, 'left boundary output')
      call prf('radinl',  e%radinl,  'left boundary input')
      call prf('radoutr', e%radoutr, 'right boundary output')
      call prf('radinr',  e%radinr,  'right boundary input')
      print *,' ' 
      print *,'miscellaneous values'
      print *,' ' 
      call prf('laser',   e%laser,   'laser energy')
      call prf('alphas',  e%alphas,  'alpha particles losses')
      call prf('fusion',  e%fusion,  'fusion energy')
      call prf('neutrons',e%neutrons,'fusion neutrons')
      call prf('error',   error,     'relative error in energy')
      print *,' ' 
      end subroutine print_energies
!======================================================================
!-----print character variable
      subroutine pra(name,value,comment)
      implicit none
      character(len=*), intent(in) :: value
      character(len=*), intent(in) :: name,comment
      character(len=12)     :: name1
      character(len=16)     :: value1
      name1    = name
      value1   = adjustl(value)
      name1    = adjustl(name1)
      value1   = adjustr(value1)
      print '(a)',' '//name1//'='//value1//' ('//comment//')'
      end subroutine pra
!======================================================================
!-----print filename
      subroutine prn(name,value,comment)
      character(len=*), intent(in) :: value
      character(len=*), intent(in) :: name,comment
      character(len=12)     :: name1
      character(len=80)     :: value1
      name1    = name
      value1   = adjustl(value)
      name1    = adjustl(name1)
      value1   = adjustr(value1)
      if(value1(65:65)/=' ')then
         value1(65:66)='..'
      end if
      print '(a)',' '//name1//'='//value1(65:80)//' ('//comment//')'
      end subroutine prn
!======================================================================
!-----print real variable
      subroutine prf(name,value,comment)
      implicit none
      real*8,           intent(in) :: value
      character(len=*), intent(in) :: name,comment
      character(len=16)            :: buffer
      write(buffer,'(g16.6)')value
      call pra(name,buffer,comment)
      end subroutine prf
!======================================================================
!-----print integer variable
      subroutine pri(name,value,comment)
      implicit none
      integer,          intent(in) :: value
      character(len=*), intent(in) :: name,comment
      character(len=16)            :: buffer
      write(buffer,'(i16)')value
      call pra(name,buffer,comment)
      end subroutine pri
!======================================================================
!-----print integer variable (if value==vc)
      subroutine prc(label,value,vc,comment)
      implicit none
      integer,          intent(in) :: value,vc
      character(len=*), intent(in) :: label,comment
      if(vc==value) call pri(label,value,comment)
      end subroutine prc
!======================================================================
!-----print table
      subroutine prt(n,labelx,labely,x,y)
      implicit none
      integer,              intent(in) :: n
      character(len=*),     intent(in) :: labelx,labely
      real*8, dimension(n), intent(in) :: x,y
      character(len=16)     :: lx,ly
      integer               :: i
      lx       = labelx
      ly       = labely
      print '(4x,a,a)',adjustr(lx),adjustr(ly)
      print *,' '
      do i=1,n
          print '(i8,2g16.6)',i,x(i),y(i)
      end do
      end subroutine prt
  
!=======================================================================
      end module multidata
!======================================================================!
!                                                                      !
!     main program                                                     !
!                                                                      !
!======================================================================!
      program multi_ife
      !use ga_routines
      use multidata
      call main
      end program multi_ife
!======================================================================
!module ga_routines
!	implicit none
       !contains
      !-----read layered initial structure
      subroutine ga_read_layer_definition(l,fileindex)
      use multidata
      implicit none
      type(typesimlayers), intent(out) :: l
      character(len=80)                :: material
      integer                          :: nc,icode
      real*8                           :: thick,r0,te0,ti0,zonpar
      namelist /layer/  nc,thick,r0,te0,zonpar,ti0,material
      character(len=32)                :: fileindex
  
      !write(*,*) fileindex
      l%nl=0
      do
         !call namelistblock('layer',l%nl+1,icode)
         call ga_namelistblock('layer',l%nl+1,icode,fileindex)
 
         if(icode==0)return
         !open(3,file='block',form='formatted')
         open(3,file='block_'//trim(fileindex),form='formatted')
         read(3,layer)
         close(3)
         l%nl           = l%nl+1
         if(l%nl>size(l%nc))then
            print *,'ERROR in subroutine read_layer_definition'
            print *,'      more than ',size(l%nc),' layers'
            stop
         end if
         l%nc     (l%nl)  = nc
         l%thick  (l%nl)  = thick
         l%zonpar (l%nl)  = zonpar
         l%r0     (l%nl)  = r0
         l%te0    (l%nl)  = te0
         l%ti0    (l%nl)  = ti0
         l%material(l%nl) = material
      end do
      end subroutine ga_read_layer_definition
!=======================================================================
!-----read material definitions
      subroutine ga_read_material_definitions(p)
      use multidata
      implicit none
      type(typesimparam), intent(inout) :: p
      integer                           :: icode
      p%nmat=0
      do
         if(p%nmat>size(p%mat))then
            print *,'ERROR in subroutine read_material_definitions'
            print *,'      more than ',size(p%mat),' materials'
            stop
         end if
         call ga_namelistblock('material',p%nmat+1,icode,p%igafile)
         if(icode==0)return
         p%nmat        = p%nmat+1
         call ga_read_one_material_definition(p%mat(p%nmat),p%igafile)
      end do
      end subroutine ga_read_material_definitions
!=======================================================================
!-----read one material definition
      subroutine ga_read_one_material_definition(m,fileindex)
      use multidata
      implicit none
      type(typematerial),    intent(out) :: m
      character(len=80)   :: name
      real*8              :: ai,zi
      character(len=80)   :: eeos_file,ieos_file,z_file
      character(len=80)   :: planck_file,ross_file,eps_file
      integer             :: eeos_id,ieos_id,z_id
      integer             :: planck_id,ross_id,eps_id
      character(len=32)   :: fileindex
      namelist /material/ name,ai,zi, &
       eeos_file,eeos_id,ieos_file,ieos_id,z_file,z_id, &
       planck_file,planck_id,ross_file,ross_id,eps_file,eps_id
      eeos_file   = "fort.13"
      ieos_file   = "fort.13"
      z_file      = "fort.13"
      planck_file = "fort.13"
      ross_file   = "fort.13"
      eps_file    = "fort.13"
      eeos_id     = 0
      ieos_id     = 0
      z_id        = 0
      planck_id   = 0
      ross_id     = 0
      eps_id      = 0
      !open(3,file='block',form='formatted')
      open(3,file='block_'//trim(fileindex),form='formatted')
      read(3,material)
      close(3)
      m%name           = name
      m%ai             = ai
      m%zi             = zi
      m%eeos_file      = eeos_file
      m%eeos_id        = eeos_id
      m%ieos_file      = ieos_file
      m%ieos_id        = ieos_id
      m%z_file         = z_file
      m%z_id           = z_id
      m%planck_file    = planck_file
      m%planck_id      = planck_id
      m%ross_file      = ross_file
      m%ross_id        = ross_id
      m%eps_file       = eps_file
      m%eps_id         = eps_id
      end subroutine ga_read_one_material_definition
!=======================================================================
!-----read pulse shape
      subroutine ga_read_pulse(p,fileindex)
      use multidata
      implicit none
      type(typesimpulse), intent(out) :: p
      integer                         :: ntab,mtab,icode
      real*8, dimension(size(p%ttab)) :: ttab,ptab
      !real*8, dimension(:)            :: gax
      real*8                          :: expo
      character(len=32)               :: fileindex
      namelist /pulse_shape/ ntab,mtab,expo,ttab,ptab
      call ga_namelistblock('pulse_shape',1,icode,fileindex)
      !call namelistblock('pulse_shape',1,icode)
      if(icode==0)return
      open(3,file='block_'//trim(fileindex),form='formatted')
      !open(3,file='block',form='formatted')
      read(3,pulse_shape)
      close(3)
      p%ntab   = ntab
      p%mode   = mtab
      p%expo   = expo
      p%ttab   = ttab
      p%ptab   = ptab
      !write(*,*) size(p%ttab)
      !write(*,*) size(gax)
      end subroutine ga_read_pulse
!=======================================================================
!-----read laser parameters (ortogonal rays)
      subroutine ga_read_laser_wkb(p,fileindex)
      use multidata
      implicit none
      type(typesimlaserwkb), intent(out) :: p
      integer :: inter,itype,icode
      real*8  :: pimax,pitime,wl,delta
      character(len=32)  :: fileindex
      namelist /pulse_wkb/ inter,pimax,pitime,wl,delta,itype
   
      call ga_namelistblock('pulse_wkb',1,icode,fileindex)
      !call namelistblock('pulse_wkb',1,icode)
      if(icode==0)then
         p%enable = 0
         return
      end if
      open(3,file='block_'//trim(fileindex),form='formatted')
      !open(3,file='block',form='formatted')
      read(3,pulse_wkb)
      close(3)
      p%enable = 1
      p%inter  = inter
      p%pimax  = pimax
      p%pitime = pitime
      p%wl     = wl
      p%itype  = itype
      p%delta  = delta
      end subroutine ga_read_laser_wkb
!=======================================================================
!-----read laser parameters (3D ray tracing)
      subroutine ga_read_laser_3d(p,fileindex)
      use multidata
      implicit none
      type(typesimlaser3d), intent(out) :: p
      integer :: itype,nr,icode
      real*8  :: pimax,pitime,wl,rmax,fwhm,bexp
      character(len=32)  :: fileindex
      namelist /pulse_3d/ pimax,pitime,wl,itype,nr,rmax,fwhm,bexp
  
      call ga_namelistblock('pulse_3d',1,icode,fileindex)
      !call namelistblock('pulse_3d',1,icode)
      if(icode==0)then
         p%enable = 0
         return
      end if
      open(3,file='block_'//trim(fileindex),form='formatted')
      !open(3,file='block',form='formatted')
      read(3,pulse_3d)
      close(3)
      p%enable = 1
      p%pimax  = pimax
      p%pitime = pitime
      p%wl     = wl
      p%itype  = itype
      p%nr     = nr
      p%rmax   = rmax
      p%fwhm   = fwhm
      p%bexp   = bexp 
      end subroutine ga_read_laser_3d
!=======================================================================
!-----read laser parameters (Maxwell's equations solver)
      subroutine ga_read_laser_maxwell(p,fileindex)
      use multidata
      implicit none
      type(typesimlasermaxwell), intent(out) :: p
      integer          :: inter,itype,idep,icode
      real*8           :: pimax,pitime,wl,angle
      character(len=1) :: pol
      character(len=32)  :: fileindex
      namelist /pulse_maxwell/ inter,pimax,pitime,wl,itype, &
                               idep,angle,pol
  
      call ga_namelistblock('pulse_maxwell',1,icode,fileindex)
      !call namelistblock('pulse_maxwell',1,icode)
      if(icode==0)then
         p%enable = 0
         return
      end if
      open(3,file='block_'//trim(fileindex),form='formatted')
      !open(3,file='block',form='formatted')
      read(3,pulse_maxwell)
      close(3)
      p%enable = 1
      p%pimax  = pimax
      p%inter  = inter
      p%pitime = pitime
      p%wl     = wl
      p%itype  = itype
      p%idep   = idep
      p%angle  = angle
      p%pol    = 1
      if(pol=='s'.or.pol=='S') p%pol=2
     end subroutine ga_read_laser_maxwell

!=======================================================================
!-----computes power or radiation temperature as a function of time
      subroutine ga_pulse(itype,pitime,pimax,time,power,puls,dtgap)
      use multidata
      implicit none
      integer,            intent(in)  :: itype
      real*8,             intent(in)  :: pitime,pimax,time
      real*8,             intent(out) :: power
      type(typesimpulse), intent(in)  :: puls
      real*8, dimension(puls%ntab)    :: ttab,ptab
      real*8                          :: factor,p1,p2
      real*8                          :: dtgap   ! for laser pulse shape
      integer                         :: i
      select case(itype)
      case(1)
         if(time<2*pitime)then
            power=pimax*(sin(cpi/2*(time/pitime)))**2
         else
            power=0
         end if
      case(2)
         if(time<0)then
            power=0
         else if(time<=pitime)then
            power=pimax
         else
            power=0
         end if
      case(4)
         if(puls%ntab==0) then
            print *,'ERROR in subroutine pulse'
            print *,'      table not available'
            stop
         end if
         ttab=puls%ttab(1:puls%ntab)*pitime
         ptab=puls%ptab(1:puls%ntab)*pimax
         do i=1,puls%ntab
            if (ttab(i)>time) exit
         end do
         if(i==1)then
            power=ptab(1)
         else if(i<=puls%ntab)then
            factor     = (time-ttab(i-1))/(ttab(i)-ttab(i-1))
            p1         = ptab(i-1)
            p2         = ptab(i)
            select case (puls%mode)
            case(1)
               p1      = p1**puls%expo
               p2      = p2**puls%expo
               power   = p1+(p2-p1)*factor
               power   = power**(1/puls%expo)
            case(2)
               p1      = log(p1)
               p2      = log(p2)
               power   = p1+(p2-p1)*factor
               power   = exp(power)
            case(3)
               p1      = exp(p1)
               p2      = exp(p2)
               power   = p1+(p2-p1)*factor
               power   = log(power)
            end select
         else
            power=ptab(puls%ntab)
        end if
!--------wfyuan, 2020/10/27, for prepulse
         if(time-23e-9+dtgap*pitime>=0 .and. time-23e-9+dtgap*pitime<=1e-9)then 
           !power = power+1*3.5e24*(sin(cpi/2*(time-1e-12)/25e-15))**2
            power = power+1*pimax
         end if
      case default
         print *,'ERROR in subroutine pulse'
         print *,'      unknown pulse type'
         stop
      end select
      end subroutine ga_pulse
!end module ga_routines	  
