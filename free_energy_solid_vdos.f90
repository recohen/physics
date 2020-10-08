! gibbs-solid.f90
! Brian's modified version of Amanuel's code

implicit none
implicit none
INTEGER :: status
integer :: ierr,arguments,iarg,l
!logical :: convertCray=.false.,convertfromLAMMPS=.false.,convertfromxyz=.false.
character(:), allocatable:: command_line,argument


integer i, j, npts, natom
real*8 k_B,h,T,P,dos,w_A,w_S,w_E,w_ZP,x,v,v1,v2,hv,kBT,n,Lx,Ly,Lz,E
real*8 dosXw_A(:),dosXw_S(:),dosXw_E(:),dosXw_ZP(:)
real*8 int_dosXw_A,int_dosXw_S,int_dosXw_E,int_dosXw_ZP
real*8 vol,q,enth,G,sum_1,sum_2,sum_3,D,sum_4
real*8 alat,ax,ay,az,bx,by,bz,cx,cy,cz
character*12 fdos, fpos
! This is joule per kelvin so SI1.6021764
! q is conversion from J to eV

parameter (k_B=1.3806503E-23,h=6.626068E-34,q=1.60217646E-19)  
parameter (eV2Hartrees=0.0367493090027428)
allocatable :: dosXw_A,dosXw_S,dosXw_E,dosXw_ZP,

! Retrieve user input
!inputfile is VDOS.dat and is normalized to 3N write(*,*) "name of VDOS.dat file in THZ normalized to 3N:"
!read(*,*) fdos
!write(*,*) "name of POSCAR file for volume calculation:"
!read(*,*) fpos
!write(*,*) "temperature in K:"
!read(*,*) T

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
         case("-v") !volume
            CALL GET_COMMAND_ARGUMENT(iarg+1, LENGTH=l, STATUS=status)
            if(allocated(argument))deallocate(argument)
            ALLOCATE(CHARACTER(l) :: argument)
            CALL GET_COMMAND_ARGUMENT(iarg+1, VALUE=argument, STATUS=status)
            IF (status /= 0) THEN
               STOP 'Sorry, but GET_COMMAND_ARGUMENT failed unexpectedly'
            END IF
            READ (argument, *) vol
            write(*,*)' Cell Volume au:',vol
         case default
            write(*,*)' free_energy_solid_vdos -n natomsPerCell -v cellvolume_au -t temperature_K'
            stop 'usage'
         end select
      enddo
   else
      write(*,*)' Command line input required. Type -h for help'
   endif



! Determine size of VDOS file
open(10,file='VDOS.dat',status='old')
npts=0
100 read(10,*,END=200) 
npts = npts + 1
goto 100
200 continue
REWIND (10)

! Allocate arrays accordingly
allocate(dosXw_A(npts),dosXw_S(npts),dosXw_E(npts))mdosXw_ZP(npts))

! Do the calculation
kBT = k_B*T
do i = 1, npts
  read(10,*)  j,v, dos
  if (i==1) v1=v
  if (i==2) v2=v  
  hv = h*v*1.0E+12
  x= hv/KBT
  w_A = -dlog(dexp(-x/2.0d0)/(1.0d0-dexp(-x)))
  W_S = x/(dexp(x)-1.0d0)-log(1.0d0-dexp(-x))
  w_E = x/2.0d0 + x/(dexp(x)-1.0d0)
  w_ZP = x/2.0d0
  if (v==0.0d0) then
    dosXw_A(i) = 0.0d0
    dosXw_S(i) = 0.0d0
    dosXw_E(i) = 0.0d0
    dosXw_ZP(i) = 0.0d0
  else
   dosXw_A(i) = dos*w_A
   dosXw_S(i) = dos*w_S
   dosXw_E(i) = dos*w_E
   dosXw_ZP(i) = dos*w_ZP
   CHECK UNITS!
  endif
enddo 
close(10)

sum_1=dosXw_A(1)+dosXw_A(npts)
sum_2=dosXw_S(1)+dosXw_S(npts)
sum_3=dosXw_E(1)+dosXw_E(npts)
sum_4=dosXw_ZP(1)+dosXw_ZP(npts)
D=4.0d0
do i = 2,npts-1
   sum_1 = sum_1 + D*dosXw_A(i)
   sum_2 = sum_2 + D*dosXw_S(i)
   sum_3 = sum_3 + D*dosXw_E(i)
   sum_4 = sum_4 + D*dosXw_ZP(i)
   D=6.0d0-D
enddo

int_dosXw_A=sum_1*(v2-v1)/3.0d0
int_dosXw_S=sum_2*(v2-v1)/3.0d0
int_dosXw_E=sum_3*(v2-v1)/3.0d0
int_dosXw_ZP=sum_4*(v2-v1)/3.0d0

kBT = kBT/q

! Print the results to screen
write(*,*) "Entropy per atom in k_B:", int_dosXw_S/dfloat(natom)
write(*,*) "Entropy per atom in eV/K:", int_dosXw_S/dfloat(natom)*8.617343d-05
write(*,*) "Entropy per atom in H/K:", int_dosXw_S/dfloat(natom)*8.617343d-05*ev2Hartrees
write(*,*) "Entropy per cell in k_B:", int_dosXw_S
write(*,*) "Entropy per cell in eV/K:", int_dosXw_S*8.617343d-05
write(*,*) "Entropy per cell in H/K:", int_dosXw_S*8.617343d-05*ev2Hartrees
write(*,*) "TS per atom in eV:", T*int_dosXw_S/dfloat(natom)*8.617343d-05
write(*,*) "TS per atom in H:", T*int_dosXw_S/dfloat(natom)*8.617343d-05*ev2Hart
write(*,*) "TS per cell in eV:", T*int_dosXw_S*8.617343d-05
write(*,*) "TS per cell in H:", T*int_dosXw_S*8.617343d-05*ev2Hartrees
write(*,*) "Zero Point Energy per atom in eV:", int_dosXw_ZP/dfloat(natom)/kbT
write(*,*) "Zero Point Energy per atom in H:", int_dosXw_ZP/dfloat(natom)/kbT*ev2Hartree
write(*,*) "Zero Point Energy per cell in eV:", int_dosXw_ZP/kbT
write(*,*) "Zero Point Energy per cell in H:", int_dosXw_ZP/kbT*ev2Hartree
write(*,*) "Internal Energy per atom in eV:", int_dosXw_E/dfloat(natom)/kbT
write(*,*) "Internal Energy per atom in H:", int_dosXw_E/dfloat(natom)/kbT*ev2Hartree
write(*,*) "Internal Energy per cell in eV:", int_dosXw_E/kbT
write(*,*) "Internal Energy per cell in H:", int_dosXw_E/kbT*ev2Hartree
write(*,*) "Free Energy per atom in eV:", int_dosXw_E/dfloat(natom)/kbT
write(*,*) "Free Energy per atom in H:", int_dosXw_E/dfloat(natom)/kbT*ev2Hartree
write(*,*) "Free Energy per cell in eV:", int_dosXw_E/kbT
write(*,*) "Free Energy per cell in H:", int_dosXw_E/kbT*ev2Hartree
! Just numbers to file (what for?)
write(*,*)&
&int_dosXw_S/dfloat(natom),&
&int_dosXw_S/dfloat(natom)*8.617343d-05,&
&int_dosXw_S/dfloat(natom)*8.617343d-05*ev2Hartrees,&
&int_dosXw_S,&
&int_dosXw_S*8.617343d-05,&
&int_dosXw_S*8.617343d-05*ev2Hartrees,&
&T*int_dosXw_S/dfloat(natom)*8.617343d-05,&
&T*int_dosXw_S/dfloat(natom)*8.617343d-05*ev2Hart,&
&T*int_dosXw_S*8.617343d-05,&
&T*int_dosXw_S*8.617343d-05*ev2Hartrees,&
&int_dosXw_ZP/dfloat(natom)/kbT,&
&int_dosXw_ZP/dfloat(natom)/kbT*ev2Hartree,&
&int_dosXw_ZP/kbT,&
&int_dosXw_ZP/kbT*ev2Hartree,&
&int_dosXw_E/dfloat(natom)/kbT,&
&int_dosXw_E/dfloat(natom)/kbT*ev2Hartree,&
&int_dosXw_E/kbT,&
&int_dosXw_E/kbT*ev2Hartree,&
&int_dosXw_E/dfloat(natom)/kbT,&
&int_dosXw_E/dfloat(natom)/kbT*ev2Hartree,&
&int_dosXw_E/kbT,&
&int_dosXw_E/kbT*ev2Hartree



END program
