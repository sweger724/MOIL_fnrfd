        subroutine rchain(urcrd,styl)
c
c a subrotuine to read initial chain coordinates. Three styles
c are supported: 
c (iii) INIT (initialize calculation reading initial and final crd files
        implicit none
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/CCRD.BLOCK'
        include 'COMMON/PATH2.BLOCK'
        include 'COMMON/SEARCH.BLOCK'
        integer urcrd
        logical find
        integer of
        character*4 styl
c local
        integer i,j,k,l,namel
        integer u1,u2
        double precision rms
        character*6 name

        name  = 'rchain'
        namel = 6

       if (styl.eq.'INIT') then
        
         if (first) then
          if (urcrd.ne.stdi) rewind urcrd
c getting unit number of reactants (CHARM format)
          call rline(name,namel,urcrd)
          if (find('file')) then
            u1 = of()
          else
            level = 1
            call alert(name,namel,'Missing file name',17,level)
          end if
c getting unit number of products (CHARM format)
          call rline(name,namel,urcrd)
          if (find('file')) then
            u2 = of()
          else
            level = 1
            call alert(name,namel,'Missing file name',17,level)
          end if

          call getcrd(u1,'CHARM')
          do  k=1,npt
            do l = 1,3
              r_initial(l,k)  = coor(l,k)
            end do
          end do
C       center r_initial to 0,0,0        
          !call rmsd_weight(npt,coor,r_initial(1,1),rms,.false.,ptms)

          call getcrd(u2,'CHARM')
          do  k=1,npt
            do l = 1,3
              r_final(l,k)  = coor(l,k)
            end do
          end do
C       center r_final to 0,0,0        
          !call rmsd_weight(npt,coor,r_final(1,1),rms,.false.,ptms)

c overlapping products with respect to reactants
      !call rmsd_weight(npt,r_initial(1,1),r_final(1,1),rms,.false.,ptms)

        end if ! first



         else if (styl.eq.'PATH') then
             rewind urcrd

C            Read very first structure
             if (first) then
               write(stdo,*)' Reading path format i = ',1
               call rpath_seq(urcrd,1)
               do  k=1,npt
                 do l = 1,3
                    r_initial(l,k) = coor(l,k)
                    centers(l,k) = coor(l,k)
                 end do
               end do
           !call rmsd_weight(npt,centers(1,1),r_initial,rms,.false.,ptms)

               do  i=2,igrid
                   write(stdo,*)' Reading path format i = ',i

                   call rpath_seq(urcrd,i)
                   write(stdo,*)' Storing path format i = ',i
                   !call rmsd_weight(npt,r_initial,coor,rms,.false.,ptms)
                   j = (i-1)*npt
                   do  k=1,npt
                      l = j + k
                      centers(1,l)  = coor(1,k)
                      centers(2,l)  = coor(2,k)
                      centers(3,l)  = coor(3,k)
                   end do
               end do

C   store final structure in the first processor
               do  k=1,npt
                 r_final(1,k)  = coor(1,k)
                 r_final(2,k)  = coor(2,k)
                 r_final(3,k)  = coor(3,k)
               end do

            end if ! first
         else 
           write(6,*)"Style not supported."
           stop
         end if
         
         return
         end
