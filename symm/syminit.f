        subroutine syminit()
        implicit none
c
c a subroutine to set up square periodic boundary condition
c a b & c are external parameters, providing translation in the x,y,z axes
c The output of this routine are the vectors
c iblock(13) and psym(ipermax)
c iblock(i) is a pointer to psym and it is the last atom of the i-th
c symmetry operation, psym is a pointer from the list of atoms related
c by symmetry to the list of neighbors lista,
c Let (x,y,z) denote the transformation along the three principal axes
c and "+" a positive transformation, "0" no transformation and "-" negative
c transformation
c the 13 transformations are:
c (+,0,+) ; (+,0,0) ; (+,0,-) ; (0,0,+) ; (+,+,+) ; (+,+,0)
c (+,+,-) ; (0,+,+) ; (0,+,0) ; (0,+,-) ; (+,-,+) ; (+,-,0) ; (+,-,-)
c  where the original system is denoted by (0,0,0)
c Note that for metal plates only the first ten transformations are used
c Noe also that it is crucial that the cutoff distance will be no more than
c half the box size.
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/SYMM.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'

        integer i
        logical first

        data first/.true./

        save first

        ipwat1 = 0
        ipwat2 = 0

        if (prll_on_off) then
              nlist1 = point1(dpoipt(monp(my_pe+1)))
              nlist2 = point2(dpoipt(monp(my_pe+1)))
              nlist3 = point3(dpoipt(monp(my_pe+1)))
              nwlist1 = poinwat1(my_nwat+1)
              nwlist2 = poinwat2(my_nwat+1)
              psymwt1(0) = nwlist1
              psymwt2(0) = nwlist2
              iblckwt1(0) = 0
              iblckwt2(0) = 0
        else 
              nwlist1 = poinwat1(nwaters)
              nwlist2 = poinwat2(nwaters)
              nlist1  = point1(npt-1)
              nlist2  = point2(npt-1)
              nlist3  = point3(npt-1)
        end if  


        indxsym1 = 0
        indxsym2 = 0
        indxsym3 = 0

        if (first) then
                first  = .false.
                if (metalyes) then 
                 symnum = 10
                else
                 symnum = 13
                end if
                do 1 i=1,symnum
                        rotop(i) = .false.
1               continue
c
c Here a new and more manual insertion of transformations is used
c in order to get no-Y transformation first. This will make metal
c walls calculations more efficient
c

c (+,0,+) (1)
        
                symop(1,1)=1
                symop(2,1)=0
                symop(3,1)=1

c (+,0,0) (2)

                symop(1,2)=1
                symop(2,2)=0
                symop(3,2)=0

c (+,0,-) (3)

                symop(1,3)=1
                symop(2,3)=0
                symop(3,3)=-1
        
c (0,0,+) (4)

                symop(1,4)=0
                symop(2,4)=0
                symop(3,4)=1

C End of transformation that remains the same even in the presence
c of metal walls

c (+,+,+) (5)

                symop(1,5)=1
                symop(2,5)=1
                symop(3,5)=1

c (+,+,0) (6)

                symop(1,6)=1
                symop(2,6)=1
                symop(3,6)=0

c (+,+,-) (7)

                symop(1,7)=1
                symop(2,7)=1
                symop(3,7)=-1

c (0,+,+) (8)

                symop(1,8)=0
                symop(2,8)=1
                symop(3,8)=1

c (0,+,0) (9)

                symop(1,9)=0
                symop(2,9)=1
                symop(3,9)=0

c (0,+,-) (10)

                symop(1,10)=0
                symop(2,10)=1
                symop(3,10)=-1

         if (.not.metalyes) then
c metal (images charges) doe snot require
c the last 3 transformations
c (+,-,+) (11)

                symop(1,11)=1
                symop(2,11)=-1
                symop(3,11)=1

c
c (+,-,0) (12)

                symop(1,12)=1
                symop(2,12)=-1
                symop(3,12)=0

c (+,-,-) (13)

                symop(1,13)=1
                symop(2,13)=-1
                symop(3,13)=-1
         end if
c end of symop initialization
        end if

c
c now generate symmetry neighbour lists for
c the different transformations
c
                do i = 1, symnum
                  call nbmsym(symop(1,i),symop(2,i),symop(3,i),i)
                end do
                

        return
        end
