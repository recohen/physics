!       program xyz_from_vasp_3types.f90
!***********************************************************
!       Create an xyz file from vasp output files
!       Made from XDATCAR file of a mixture of 3 types
!***********************************************************
        implicit none
        
        ! n = max # of timesteps, m = max # of atoms
        integer maxSteps, maxAtoms
        parameter (maxSteps=10000000,maxAtoms=100000)
        integer i, j
        integer natom1, natom2, natom3
        real x, y, z, alat, a, b, c, d0, d1, d2
        character*2 typat1, typat2, typat3

        write(*,*) 'Atom type 1'
        read(*,*) typat1
        write(*,*) 'Atom type 2'
        read(*,*) typat2
        write(*,*) 'Atom type 3'
        read(*,*) typat3

        ! Open the input and output files
        open(1,file='XDATCAR',status='old',ERR=90)
        open(2,file='POSCAR',status='old',ERR=80)
        open(3,file='TRAJEC.xyz')

        ! Get cell and natom from POSCAR file
        read(2,*) 
        read(2,*) alat
        read(2,*) a, d1, d2
        read(2,*) d0, b, d2
        read(2,*) d0, d1, c
! VASP 5 ONLY        read(2,*) typat1, typat2, typat3
        read(2,*) natom1, natom2, natom3
        close(2)

        a = a*alat
        b = b*alat
        c = c*alat

        ! Read off XDATCAR header
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*) 

        ! Loop over timesteps
        do i=1,maxSteps

          ! Read off the whitespace line
          read(1,*,END=100)

          ! Write natom and timestep to xyz file
          write(3,*) natom1+natom2+natom3
          write(3,*) i

          ! Loop over natom1
          do j=1,natom1

            ! Read in
            read(1,*,END=100) x, y, z
            
            ! Convert coords from reduced to angstroms
            x = x*a
            y = y*b
            z = z*c

            ! Write converted coords to TRAJEC.xyz
            write(3,*) typat1, x, y, z

          enddo

          ! Loop over natom2
          do j=1,natom2

            ! Read in
            read(1,*,END=100) x, y, z
            
            ! Convert coords from reduced to angstroms
            x = x*a
            y = y*b
            z = z*c

            ! Write converted coords to TRAJEC.xyz
            write(3,*) typat2, x, y, z

          enddo

          ! Loop over natom3
          do j=1,natom3

            ! Read in
            read(1,*,END=100) x, y, z
            
            ! Convert coords from reduced to angstroms
            x = x*a
            y = y*b
            z = z*c

            ! Write converted coords to TRAJEC.xyz
            write(3,*) typat3, x, y, z

          enddo

        enddo

 100    continue

        close(3)

 80     continue

        close(1)

 90     continue

        END
