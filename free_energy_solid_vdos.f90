! gibbs-solid.f90
! Brian's modified version of Amanuel's code
! modified heavily by Ronald E. Cohen
implicit none
INTEGER :: status
integer :: ierr,arguments,iarg,l
!logical :: convertCray=.false.,convertfromLAMMPS=.false.,convertfromxyz=.false.
character(:), allocatable:: command_line,argument


integer i, j, npts, natom
real*8 k_B,h,T,P,dos,w_A,w_S,w_E,w_ZP,x,v,v1,v2,hv,kBT,n,Lx,Ly,Lz,E
real*8, allocatable :: dosXw_A(:),dosXw_S(:),dosXw_E(:),dosXw_ZP(:)
real*8 :: WCor_A,WCor_S,WCor_E
real*8, allocatable :: dosWCor_A(:),dosWCor_S(:),dosWCor_E(:)
real*8 int_dosWCor_A,int_dosWCor_S,int_dosWCor_E
real*8 int_dosXw_A,int_dosXw_S,int_dosXw_E,int_dosXw_ZP
real*8 CorE,CorAH,CorS,eV2Hartrees,kBeV
real*8 q,enth,G,sum_1,sum_2,sum_3,D,sum_4,sum_5,sum_6,sum_7
character*12 fdos, fpos
! This is joule per kelvin so SI 1.6021764
! q is conversion from J to eV

parameter (k_B=1.380649e-23,h=6.626070E-34,q=1.60217646E-19)  
!kB J/K h Js 
parameter (eV2Hartrees=0.0367493090027428)
parameter (kBeV=8.617333262145e-5)
!inputfile is VDOS.dat and is normalized to 3N write(*,*) "name of VDOS.dat file in THZ normalized to 3N:"

! read command line arguments
arguments=COMMAND_ARGUMENT_COUNT()
IF (arguments /= 0) THEN
   do iarg=1,arguments,2 !assumes come in pairs
      CALL GET_COMMAND_ARGUMENT(iarg, LENGTH=l, STATUS=status)
      if(allocated(command_line))deallocate(command_line)
      ALLOCATE(CHARACTER(l) :: command_line)
      CALL GET_COMMAND_ARGUMENT(iarg, VALUE=command_line, STATUS=status)
      IF (status /= 0) THEN
         STOP 'Sorry, but GET_COMMAND_ARGUMENT failed unexpectedly'
      END IF
      select case (command_line)
      case ("-t") !temperature
            CALL GET_COMMAND_ARGUMENT(iarg+1, LENGTH=l, STATUS=status)
            if(allocated(argument))deallocate(argument)
            ALLOCATE(CHARACTER(l) :: argument)
            CALL GET_COMMAND_ARGUMENT(iarg+1, VALUE=argument, STATUS=status)
            IF (status /= 0) THEN
               STOP 'Sorry, but GET_COMMAND_ARGUMENT failed unexpectedly'
            END IF
            READ (argument, *) T 
            write(*,*)' Temperature K:',T
         case("-n") !natoms
            CALL GET_COMMAND_ARGUMENT(iarg+1, LENGTH=l, STATUS=status)
            if(allocated(argument))deallocate(argument)
            ALLOCATE(CHARACTER(l) :: argument)
            CALL GET_COMMAND_ARGUMENT(iarg+1, VALUE=argument, STATUS=status)
            IF (status /= 0) THEN
               STOP 'Sorry, but GET_COMMAND_ARGUMENT failed unexpectedly'
            END IF
            READ (argument, *) natom
            write(*,*)' Atoms per cell:',natom
         case default
            write(*,*)' free_energy_solid_vdos -n natomsPerCell -t temperature_K'
            stop 'usage'
         end select
      enddo
   else
      write(*,*)' Command line input required. Type -h for help'
      stop 'input'
   endif



! Determine size of VDOS file, 3 columns. Column 2 Freq THz. Col3 VDOS normalized to 3N
open(10,file='VDOS.dat',status='old')
npts=0
100 read(10,*,END=200) 
npts = npts + 1
goto 100
200 continue
REWIND (10)

! Allocate arrays accordingly
allocate(dosXw_A(npts),dosXw_S(npts),dosXw_E(npts),dosXw_ZP(npts))
allocate(dosWCor_A(npts),dosWCor_S(npts),dosWCor_E(npts))
! Do the calculation
kBT = k_B*T
do i = 1, npts
  read(10,*)  j,v, dos
  if (i==1) v1=v
  if (i==2) v2=v  
  hv = h*v*1.0E+12
  x= hv/KBT
  w_A = -dlog(dexp(-x/2.0d0)/(1.0d0-dexp(-x)))
  WCor_A = w_A - log(x)
! Was  w_S = x/(dexp(x)-1.0d0)-log(1.0d0-dexp(-x)) agrees with Berens et alk. Maybe equivalent? Yes should be same
  w_S = x*dexp(-x)/(1.0d0-dexp(-x))-log(1.0d0-dexp(-x))
  WCor_S = x/(dexp(x)-1d0)-log(1.0-dexp(-x))+log(x)-1d0
  w_E = x/2.0d0 + x/(dexp(x)-1.0d0)
  WCor_E = w_E - 1d0
  w_ZP = x/2.0d0
  if (v==0.0d0) then
    dosXw_A(i) = 0.0d0
    dosXw_S(i) = 0.0d0
    dosXw_E(i) = 0.0d0
    dosXw_ZP(i) = 0.0d0
    dosWcor_A(i) = 0.0d0
    dosWcor_S(i) = 0.0d0
    dosWcor_E(i) = 0.0d0
 else
   dosXw_A(i) = dos*w_A
   dosXw_S(i) = dos*w_S
   dosXw_E(i) = dos*w_E
   dosXw_ZP(i) = dos*w_ZP
   dosWcor_A(i) = dos*Wcor_A
   dosWcor_S(i) = dos*Wcor_S
   dosWcor_E(i) = dos*Wcor_E
  endif
enddo 
close(10)

sum_1=dosXw_A(1)+dosXw_A(npts)
sum_2=dosXw_S(1)+dosXw_S(npts)
sum_3=dosXw_E(1)+dosXw_E(npts)
sum_4=dosXw_ZP(1)+dosXw_ZP(npts)
sum_5=dosWcor_A(1)+dosWcor_A(npts)
sum_6=dosWcor_S(1)+dosWcor_S(npts)
sum_7=dosWcor_E(1)+dosWcor_E(npts)
D=4.0d0
do i = 2,npts-1
   sum_1 = sum_1 + D*dosXw_A(i)
   sum_2 = sum_2 + D*dosXw_S(i)
   sum_3 = sum_3 + D*dosXw_E(i)
   sum_4 = sum_4 + D*dosXw_ZP(i)
   sum_5 = sum_5 + D*dosWcor_A(i)
   sum_6 = sum_6 + D*dosWcor_S(i)
   sum_7 = sum_7 + D*dosWcor_E(i)
   D=6.0d0-D
enddo

int_dosXw_A=sum_1*(v2-v1)/3.0d0
int_dosXw_S=sum_2*(v2-v1)/3.0d0
int_dosXw_E=sum_3*(v2-v1)/3.0d0
int_dosXw_ZP=sum_4*(v2-v1)/3.0d0
int_dosWcor_A=sum_5*(v2-v1)/3.0d0
int_dosWcor_S=sum_6*(v2-v1)/3.0d0
int_dosWcor_E=sum_7*(v2-v1)/3.0d0
kBT = kBT/q

! Print the results to screen
write(*,*) "Entropy per atom in k_B:", int_dosXw_S/dfloat(natom)
write(*,*) "Entropy per atom in eV/K:", int_dosXw_S/dfloat(natom)&
     &*kBeV
write(*,*) "Entropy per atom in H/K:", int_dosXw_S/dfloat(natom)&
     &*kBeV*ev2Hartrees 
write(*,*) "Entropy per cell in k_B:", int_dosXw_S
write(*,*) "Entropy per cell in eV/K:", int_dosXw_S*kBeV
write(*,*) "Entropy per cell in H/K:", int_dosXw_S*kBeV&
     &*ev2Hartrees 
write(*,*) "TS per atom in eV:", T*int_dosXw_S/dfloat(natom)*kBeV
write(*,*) "TS per atom in H:", T*int_dosXw_S/dfloat(natom)*kBeV*ev2Hartrees 
write(*,*) "TS per cell in eV:", T*int_dosXw_S*kBeV
write(*,*) "TS per cell in H:", T*int_dosXw_S*kBeV*ev2Hartrees
write(*,*) "Zero Point Energy per atom in eV:", kbT*int_dosXw_ZP&
     &/dfloat(natom) 
write(*,*) "Zero Point Energy per atom in H:", kbT*int_dosXw_ZP&
     &/dfloat(natom)*ev2Hartrees 
write(*,*) "Zero Point Energy per cell in eV:", kbT*int_dosXw_ZP
write(*,*) "Zero Point Energy per cell in H:", kbT*int_dosXw_ZP*ev2Hartrees
write(*,*) "Internal Energy per atom in eV:", kbT*int_dosXw_E/dfloat(natom)
write(*,*) "Internal Energy per atom in H:", kbT*int_dosXw_E&
     &/dfloat(natom)*ev2Hartrees 
write(*,*) "Internal Energy per cell in eV:", kbT*int_dosXw_E
write(*,*) "Internal Energy per cell in H:", kbT*int_dosXw_E*ev2Hartrees
write(*,*) "Free Energy per atom in eV:", kbT*int_dosXw_E/dfloat(natom)
write(*,*) "Free Energy per atom in H:", kbT*int_dosXw_E&
     &/dfloat(natom)*ev2Hartrees 
write(*,*) "Free Energy per cell in eV:", kbT*int_dosXw_E
write(*,*) "Free Energy per cell in H:", kbT*int_dosXw_E*ev2Hartrees
! quantum corrections to classical MD
write(*,*)' Quantum corrections to classical MD. H/cell: E(H), A, S '
write(*,*)int_dosWcor_E*kbT*ev2Hartrees,int_dosWcor_A*kbT*ev2Hartrees&
     &,int_dosWcor_S*kBeV*ev2Hartrees 



END program
