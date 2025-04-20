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
