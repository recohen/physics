!       program xyz_from_mol.f90
!***********************************************************
!       Convert mol file back to xyz file.
!***********************************************************
        implicit none
        
        integer maxSteps, maxAtoms, i, j, k
        parameter (maxSteps=100000,maxAtoms=1024)
        integer natom
        real x, y, z
        character*8 typat,d1,d2,d3
        character*128 fmol

        ! Retrieve user input
        write(6,*) 'Name of .mol file'
        read(5,*) fmol

        ! Open input and output files
        open(1,file=fmol,status='old',ERR=90)
        open(2,file='TRAJEC.xyz')

        ! Read in .mol header
        read(1,*) d1, d2, d3, natom
        read(1,*) 
        read(1,*) 
        read(1,*) 
        read(1,*) 
        read(1,*) 
        read(1,*) 
        read(1,*) 
        read(1,*) 
        read(1,*) 

        do i=1,maxSteps

          read(1,*,END=100) typat, x, y, z
          write(2,'(i3)') natom
          write(2,'(i6)') i
          write(2,'(A1,x,f10.6,x,f10.6,x,f10.6)') typat, x, y, z

          do j=1,natom-1

            read(1,*,END=100) typat, x, y, z
            write(2,'(A1,x,f10.6,x,f10.6,x,f10.6)') typat, x, y, z

          enddo
          
        enddo

 100    continue

        close(1)
        close(2)

 90     continue

        END
