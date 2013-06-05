! program coordination.f90
!****************************************************
! Calculate the coordination of all atoms in a
! simulation over time and the system average
!****************************************************

implicit none

integer maxSteps, maxAtoms, i, j, k, l
parameter (maxSteps=100000,maxAtoms=1024)
integer natom, tstep, pbc_round(3), coord(maxAtoms), max_coord
real*8 x(maxAtoms), y(maxAtoms), z(maxAtoms), input_value(3)
real*8 alat, blat, clat, dx, dy, dz, r, cutoff, avg(maxAtoms)
character*2 typat(maxAtoms)
character*128 fxyz

! Retrieve user input
write(*,*) "Name of xyz file:"
read(*,*) fxyz
write(*,*) "alat, blat, clat (angstroms):"
read(*,*) alat, blat, clat
write(*,*) "cutoff (angstroms):"
read(*,*) cutoff

! Open necessary files
open(1,file=fxyz,status='old',ERR=90)
open(2,file='coordination.xyz',status='unknown')
open(3,file='coordination.dat',status='unknown')

5 format(A2,x,i3)

! Perform the calculation
do i = 1, maxSteps

  ! Read from xyz file
  read(1,*,END=100) natom
  read(1,*) tstep

  ! Write to coordination.xyz file
  write(2,*) natom
  write(2,*) tstep

  ! Read snapshot coordinates
  do j = 1, natom
    read(1,*) typat(j), x(j), y(j), z(j)
  enddo

  ! Loop over all atom pairs
  do j = 1, natom

    ! Reset coordination list to zero
    coord(j) = 0

    do k = 1, natom
      if (j.ne.k) then

        ! Calculate the distance between atom j & k
        dx = X(k) - X(j)
        dy = Y(k) - Y(j)
        dz = Z(k) - Z(j)
        input_value(1) = dx/alat
        input_value(2) = dy/blat
        input_value(3) = dz/clat
        do l = 1, 3
          pbc_round(l) = int(input_value(l))
          if (abs(input_value(l)-pbc_round(l)).ge.0.5) then
            if (input_value(l).gt.0) pbc_round(l) = pbc_round(l) + 1
            if (input_value(l).lt.0) pbc_round(l) = pbc_round(l) - 1
          endif
        enddo
        dx = dx - alat*pbc_round(1)
        dy = dy - blat*pbc_round(2)
        dz = dz - clat*pbc_round(3)
        r = ( dx**2 + dy**2 + dz**2 )**0.5

        ! If k is within cutoff distance of j
        if (r.lt.cutoff) then
          coord(j) = coord(j) + 1
        endif
      endif
    enddo

    ! Write each atom's coordination to xyz style file
    write(2,5) typat(j), coord(j)

    if (coord(j).gt.max_coord) then
      max_coord = coord(j)
    endif

    avg(coord(j)+1) = avg(coord(j)+1) + 1

  enddo

enddo

100 continue

tstep = i - 1
6 format(i3,x,f12.8)

do i = 1, max_coord + 2
  write(3,6) i-1, avg(i) / (tstep*natom)
enddo

close(1)
close(2)

90 continue

END
