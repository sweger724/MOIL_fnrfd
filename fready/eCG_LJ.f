        Subroutine eCG_LJ()

C       Lennard Jones CG energy terms
        implicit none

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/SPECL.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/FREADY.BLOCK'

        double precision rx,ry,rz,r2,s,s2,s3,s4,s6,s12,E,df,r
        double precision xi,yi,zi,dxi,dyi,dzi, x(0:9)
        integer k,i,j, T, loop

        e_vdw = 0.d0

c
c first is the list of particle pairs up to cutvdw
c included charged pairs so we calculate BOTH vdw and ele
c
        do i=1,npt-1

                xi  = coor(1,i)
                yi  = coor(2,i)
                zi  = coor(3,i)
                dxi = 0.d0
                dyi = 0.d0
                dzi = 0.d0

                do k=LJ_list1(i-1)+1,LJ_list1(i)
                  
                  j=LJ_list2(k)
                  rx = xi - coor(1,j)
                  ry = yi - coor(2,j)
                  rz = zi - coor(3,j)
                  r2=rx*rx+ry*ry+rz*rz
                  r = dsqrt(r2)
                  T = LJ_Type(k)
                  call eCG_NB(r,T,E,df)
C                  if ( moname(poimon(i)).eq. "CGTR" .or. 
C     &                 moname(poimon(j)).eq. "CGTR" ) then
C               write(stdo,*)'~~~~~~~~~~~~~~~~~~~~~~~~'
C               write(stdo,*)'i and j are ',ptnm(i),ptnm(j),T
C          write(stdo,*)'Mi, Mj:',moname(poimon(i))," ",moname(poimon(j))
C               write(stdo,*)'r, eLJ is ',r,E
C               write(*,*)"df is ",df
C               write(stdo,*)'~~~~~~~~~~~~~~~~~~~~~~~~~'
C               end if

                if (E*LJ_14(k) .gt. E_CG_max) then
      write(stdo,'(a,i4,1x,a5,a1,a3,l,i4,1x,a5,a1,a3,f6.2,1x,f8.4)')
     &   ,"Enb:",poimon(i),moname(poimon(i)),"-",ptnm(i),smooth_hardcore
     &   ,poimon(j),moname(poimon(j)),"-",ptnm(j),r,E*LJ_14(k)
                endif
                
                E = E * LJ_14(k)
                df = df * LJ_14(k)

                e_vdw = e_vdw + E 

                rx = df*rx
                ry = df*ry
                rz = df*rz
                
                dxi = dxi + rx
                dyi = dyi + ry
                dzi = dzi + rz
                
                dpot(1,j) = dpot(1,j) - rx
                dpot(2,j) = dpot(2,j) - ry
                dpot(3,j) = dpot(3,j) - rz

c               write (*,*) "index j is ",j
c               write (*,*) "dpot1 is",dpot(1,j)
c               write (*,*) "dpot2 is",dpot(2,j)
c               write (*,*) "dpot3 is",dpot(3,j)        

             enddo                
                
                dpot(1,i) = dpot(1,i) + dxi
                dpot(2,i) = dpot(2,i) + dyi
                dpot(3,i) = dpot(3,i) + dzi 
        end do

        return 
        end
