      subroutine CG_nb_list()
      
      implicit none
            
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/FREADY.BLOCK'

c      local variables

      integer i,j,l,counter, id1, id2, istart,mi,mj,irand
      integer ipt,jpt,lpt,i1pt,j1pt,l1pt,k1,k
      double precision rx,ry,rz,r2,cutoff2, small_cutoff
      logical CYSCYSbond

       integer TMP_list1(0:maxmono),TMP_list2(maxCGcontacts)
     &          ,TMP_type(maxCGcontacts),minK1,minK2, m, T
       double precision ENE(maxCGcontacts),fA,minE,E
       double precision rrj(3),pi(3),pj(3),pl(3),di(3),dj(3),dl(3)
       double precision cosA,cosB,cosC,dEdr,dEdA,dEdB,dEdC,si,sj,sij

     

      cutoff2= CG_cutoff*CG_cutoff
      small_cutoff = cutoff2
      
      counter =0
      LJ_list1(0)=0 
c
c.....generate the LJ lists
c
      do i = 1, npt-1
        mi = poimon(i)
        id1 = CGid(i)
        istart = npt+1
        if (ptnm(i).eq."CA  ") then
          if (mi+1 .le. totmon) istart = poipt(mi+1) + 2
        else
           istart = i + 2
        endif 
C        write(6,*)"istart: ",istart
        do j =istart, npt

            rx = coor(1,i) - coor(1,j)
            ry = coor(2,i) - coor(2,j)
            rz = coor(3,i) - coor(3,j)
            
            r2 = rx*rx + ry*ry + rz*rz
            mj = poimon(j)

       if (.FALSE.) then
        
         if (r2.lt.(15.0**2) .and .mj.gt.mi+9 ) then
           id2 = CGid(j)
           if (id1.eq.10 .and. id2.eq.10) then
             if (CYSCYSbond(mi,mj)) goto 123
           endif
           if (ptnm(i).eq."CM  " .and. ptnm(j).eq."CM  " ) then
             write(6,*)"Distance: ",mj-mi,
     &                  moname(mi)," ",moname(mj),sqrt(r2)
           endif
           if (ptnm(i).eq."CA  " .and. ptnm(j).eq."CA  " ) then
             write(6,*)"DistCACA: ",mj-mi,
     &                  moname(mi)," ",moname(mj),sqrt(r2)
           endif
           if (ptnm(i).eq."CA  " .and. ptnm(j).eq."CM  ") then
             write(6,*)"DistCACM: ",mj - mi,
     &                  moname(mi)," ",moname(mj),sqrt(r2)
           endif
           if(ptnm(i).eq."CM  " .and. ptnm(j).eq."CA  ") then
             write(6,*)"DistCACM: ",mj-mi,
     &                  moname(mj)," ",moname(mi),sqrt(r2)
           endif
         endif
       endif

        
            if (r2 .lt. small_cutoff .or. (r2.lt.cutoff2 .and. 
     &       ptnm(i).eq."CM  " .and. ptnm(j).eq."CM  " )) then
                id2 = CGid(j)
C               check whether it isn't a CYS - CYS bond         
                if (ptnm(i).eq."CM  " .and. ptnm(j).eq."CM  ") then
                  if (id1.eq.10 .and. id2.eq.10) then
                     if (CYSCYSbond(mi,mj)) goto 123
                  endif
                endif 
                counter = counter + 1
                LJ_list2(counter) = j
C               write(6,*)"Particles: ", ptnm(i), ptnm(j)
                LJ_Type(counter) = Namino*(id1-1) + id2
                if(mj .eq. mi + 1) then 
                    LJ_14(counter) = 0.3d0
                else
                    LJ_14(counter) = 1.0d0
                endif
                if (ptnm(i).eq."CA  " .and. ptnm(j).eq."CA  ") then
                  LJ_Type(counter) = LJ_Type(counter) + 2*Namino**2
                  if(mj .eq. mi + 3) LJ_14(counter) = 0.3d0
                endif
                
                if (ptnm(i).eq."CM  " .and. ptnm(j).eq."CA  ") then
                  LJ_Type(counter) = LJ_Type(counter) + Namino**2
                  if(mj .eq. mi + 2) LJ_14(counter) = 0.3d0
                endif
                if (ptnm(i).eq."CA  " .and. ptnm(j).eq."CM  ") then
                  LJ_Type(counter) = Namino*(id2-1) + id1 + Namino**2
                  if(mj .eq. mi + 2) LJ_14(counter) = 0.3d0 
                endif
C               LJ_14(counter) = 0.0d0
C                write(6,*)'Pair i,j',i,j
            endif
123         continue            
         end do 
         
         LJ_list1(i)=counter
C         write(6,*)'LJ_list1(i) ',i,LJ_list1(i)-LJ_list1(i-1) 
      end do
C      write(6,*)'Contacts in the structure: ',counter, CG_cutoff
c
c.....generate the HB lists
c
      counter = 0
      TMP_list1(0)=0

       do 43 i = 1, totmon
        ipt = poipt(i-1)+1
        i1pt = poipt(i)+1

        do m = 1,3
          pi(m) = 0.5d0 * (coor(m,ipt) + coor(m,i1pt))
          di(m) =coor(m,i1pt) - coor(m,ipt)
        end do
        
        if(moname(i).eq."CGTR" .or. moname(i+1).eq."CGTR") then
          TMP_list1(i)=counter 
          goto 43
        endif

         do 45 j =1,totmon 
          jpt = poipt(j-1)+1
          j1pt = poipt(j)+1

          do m =1,3
            pj(m) = 0.5d0 * (coor(m,jpt) + coor(m,j1pt))
            dj(m) = coor(m,j1pt) - coor(m,jpt)
            rrj(m) = pj(m) - pi(m)
          end do

          if(j.lt.i+3 .and. j.gt.i-3) goto 45
          if(moname(j).eq."CGTR" .or. moname(j+1).eq."CGTR") goto 45
            rx = 0.5d0*(coor(1,jpt) + coor(1,j1pt) 
     &                - coor(1,ipt) - coor(1,i1pt))
            ry = 0.5d0*(coor(2,jpt) + coor(2,j1pt) 
     &                - coor(2,ipt) - coor(2,i1pt))
            rz = 0.5d0*(coor(3,jpt) + coor(3,j1pt) 
     &                - coor(3,ipt) - coor(3,i1pt))

            r2 = rx*rx + ry*ry + rz*rz

            if (r2 .lt. HB_cutoff**2 ) then
                counter = counter + 1
                TMP_list2(counter) = j
                TMP_type(counter) = 1
                if (moname(i).eq."PROZ") 
     &             TMP_type(counter) = TMP_type(counter) + 1
                if (moname(j).eq."PROZ") 
     &             TMP_type(counter) = TMP_type(counter) + 1

C  calculate cos Alpha, Beta, Gamma

                 si = dsqrt(di(1)**2 + di(2)**2 + di(3)**2)
                 sj = dsqrt(dj(1)**2 + dj(2)**2 + dj(3)**2)
                 sij = dsqrt(rrj(1)**2 + rrj(2)**2 + rrj(3)**2)

                cosA = (di(1)*dj(1) + di(2)*dj(2) + di(3)*dj(3))/(si*sj)
            cosB = (di(1)*rrj(1) + di(2)*rrj(2) + di(3)*rrj(3))/(si*sij)
            cosC = (dj(1)*rrj(1) + dj(2)*rrj(2) + dj(3)*rrj(3))/(sj*sij)
                T = TMP_Type(counter)

                call eCG_backboneNB(T,sij,cosA,cosB,cosC,E,dEdr,
     &                              dEdA,dEdB,dEdC,.FALSE.)
                ENE(counter) = E
            endif
45     continue


        TMP_list1(i)=counter
C        write(6,*)'TMP_list1(i) ',i,TMP_list1(i)-TMP_list1(i-1)
43     continue

C        write(6,*) "HB contacts:", counter, HB_cutoff

        counter = 0

C         Find the actuall bonds that maximize HB energy
          do i=1,totmon
               minE=1.d10
               minK1 = 0
               minK2 = 0

               ipt = poipt(i-1)+1
               i1pt = poipt(i)+1

               do m=1,3
                 pi(m)  = 0.5d0 * (coor(m,ipt) + coor(m,i1pt))
               end do

               if (moname(i).ne."PROZ") then
               do k=TMP_list1(i-1)+1,TMP_list1(i)-1
                 j = TMP_list2(k)

                 jpt = poipt(j-1)+1
                 j1pt = poipt(j)+1

                 do m=1,3
                   pj(m)  = 0.5d0 * (coor(m,jpt) + coor(m,j1pt))
                 end do

                 do k1 = k+1,TMP_list1(i)
                  l = TMP_list2(k1)

C    this IF is a heuristic to speed it up the calculation
                  if ((l-j).gt.5) then

                    lpt = poipt(l-1)+1
                    l1pt = poipt(l)+1

                    do m=1,3
                      pl(m)  = 0.5d0 * (coor(m,lpt) + coor(m,l1pt))
                    end do

                    do m=1,3
                      dj(m) = pj(m) - pi(m)
                      di(m) = pi(m) - pl(m)
                    end do

                    si = dsqrt(di(1)**2 + di(2)**2 + di(3)**2)
                    sj = dsqrt(dj(1)**2 + dj(2)**2 + dj(3)**2)

                    cosA = (di(1)*dj(1) + di(2)*dj(2) + di(3)*dj(3))
     &                   /(si*sj)

                    E = ( ENE(k)+ENE(k1) ) * fA(cosA)

                    if (E.lt.minE) then
                      minK1 = k
                      minK2 = k1
                      minE = E
                    endif

                  endif ! if ((l-j).gt.5)

                 end do
               end do

               endif !  if PROZ

C   now check only a single HB from atom i
               do k=TMP_list1(i-1)+1,TMP_list1(i)
                  if (ENE(k).lt.minE) then
                    minK1 = k
                    minK2 = 0
                    minE = ENE(k)
                  endif
               end do

C  generate the final lists:

               if (minK1 .ne. 0) then
                 HBbond1(i)=TMP_list2(minK1)
                 HBtype1(i) = TMP_Type(minK1)
               else
                 HBbond1(i) = 0
                 HBtype1(i) = 0
               endif

               if (minK2 .ne. 0) then
                 HBbond2(i)=TMP_list2(minK2)
                 HBtype2(i) = TMP_Type(minK2)
               else
                 HBbond2(i) = 0
                 HBtype2(i) = 0
               endif 

         end do

      end


C*******************************************************
      function CYSCYSbond(i,j)
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/FREADY.BLOCK'
        
        logical CYSCYSbond
        integer i,j,k

        do k=1,Ncyscys
           if (SSbond(k,1).eq.i .and. SSbond(k,2).eq.j) then
              CYSCYSbond = .true.
C             write(6,*)"CYS - CYS bond detected."
              return 
           end if
           if (SSbond(k,1).eq.j .and. SSbond(k,2).eq.i) then
              CYSCYSbond = .true.
C             write(6,*)"CYS - CYS bond detected."
              return
           end if
        enddo
        
        CYSCYSbond = .false.
        return
        
      end
