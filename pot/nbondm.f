        subroutine nbondm()
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/LINE2.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/SPECL.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        integer nlist

C-------local variables--------------------------
        integer i,j,k,temp
        double precision rxyz,r2
        integer pointm(maxmono*maxdiv)
        integer nblistm(ihugemo)
        integer nwlist1,nwlist2,iwat

c watflag is true if the i-th monomer is water
c watflag2 is true if the i-th and the j-th monomers are waters
c
        logical sflag,watflag,watflag2,nflag


c -------------------------------------------------------------
c yael - new variables
        integer lstart,lend

        data sflag/.true./
        save sflag

        if (.not. (evdyes.or.eelyes.or.eCGyes)) return

C--------------------------------------------------------------------
C    if CG model:
        if (eCGyes) then
         call CG_nb_list()
         return
        endif


c--------------------------------------------------------------------
c initial partition 
        if (prll_on_off) then
           if (sflag) then 
              temp = totdmon / num_pes 
              do i=0,num_pes-1
                 monp(i) = temp*i
              end do   
              sflag = .false.
           else
              do i=0,num_pes-1
                 monp(i) = new_monp(i)
              end do   
           endif
           monp(num_pes) = totdmon
        endif

c If cutmono2 was not already set...
c set cutmono2 (cutoff distance^2 for monomer list) larger by 25A^2
c (trying now only 9 - 3^2A^2
c compared to cutebig2

        if (cutmono2.lt.0) cutmono2 = cutebig2 + 9.d0
c loop over first monomer
        if (debug) then
         write(stdo,*)' in nbondm npt =',npt
         write(stdo,1000)(dpoipt(i),i=1,totdmon)
1000     format(10(1x,i6))
        end if

c-----------------------------------------
C### Parallelization of cm of monomers left for future ref. when comm.
c### will be faster
c
	no_of_away =0
        do 2 i=1,totdmon
        if(vp_flag) call vp_locate()
        if (moname(realmono(i)).eq.'TIP3' .or.
     1        moname(realmono(i))(1:3).eq.'SPC') then
                j = poipt(realmono(i))-2
                mono_cent(1,i) = coor(1,j)
                mono_cent(2,i) = coor(2,j)
                mono_cent(3,i) = coor(3,j)
c RE 5/1/2011: A hack to find water molecules sufficiently far from
c the center of the box (and the solute. Here the minimal disance
c is hard cded to 25/30. Improvement possible and likely.
c
                if (scale_away) then
                 r2 = 0.d0
                 do k=1,3
                  r2 = r2 + mono_cent(k,i)*mono_cent(k,i)
                 end do
                 if(r2.gt.900) then
                  do k=0,2
                   no_of_away = no_of_away + 1
                   far_water(no_of_away) = j + k
                  end do
                 end if
                end if
                go to 2
         end if
         mono_cent(1,i) = 0.d0
         mono_cent(2,i) = 0.d0
         mono_cent(3,i) = 0.d0
        if (mdivlist(0,realmono(i)).eq.-1) then
           lstart = mdivlist(1,realmono(i)) + 1
           lend = mdivlist(2,realmono(i))
        else
           lstart = dpoipt(i-1)+1
           lend   = dpoipt(i)
        endif

         do 1 j=lstart,lend
          do 1 k=1,3
          mono_cent(k,i) = mono_cent(k,i) + coor(k,j)
1        continue
         r2 = 1.d0/dfloat(lend-lstart+1)
         mono_cent(1,i)  = mono_cent(1,i)*r2
         mono_cent(2,i)  = mono_cent(2,i)*r2
         mono_cent(3,i)  = mono_cent(3,i)*r2
2       continue
c nlist = index for monomer non-bonding list
        nlist = 0
        nwlist1 = 0
        nwlist2 = 0

        if (prll_on_off) then
           lstart = monp(my_pe)+1
           lend = monp(my_pe+1)
        else
           lstart = 1
           lend = totdmon
        endif   

        nflag = .false.
        iwat = 0
        poinwat1(0) = 0
        poinwat2(0) = 0

        do 100  i=lstart,lend
c set pointers, include self in list

c check if water.
                watflag = (moname(realmono(i)).eq.'TIP3')
     1          .or. (moname(realmono(i))(1:3).eq.'SPC')

                nlist          = nlist + 1
                pointm(i)      = nlist
                nblistm(nlist) = i
                if (watflag) then

                 !if all neighbors of previous monomer (i-1) were waters and
                 ! monomer i-1 was water itself it can be removed from the list
                 if (nflag) then
                  nlist         = nlist - 1
                  pointm(i)     = nlist
                  nblistm(nlist) = i
                 end if

                 iwat          = iwat + 1
                 indxwat(iwat) = i
                 poinwat1(iwat)= nwlist1
                 poinwat2(iwat)= nwlist2
                 watptr(i)     = iwat   
                 nflag = .true.
                else
                 nflag = .false.
                end if
        
                if (i.eq.totdmon) go to 100
c loop over second monomer
                do 200 j=i+1,totdmon
                        rxyz=mono_cent(1,i) - mono_cent(1,j)
                        r2 = rxyz*rxyz
                        if (r2.gt.cutmono2) go to 200
                        rxyz=mono_cent(2,i) - mono_cent(2,j)
                        r2 = r2 + rxyz*rxyz
                        if (r2.gt.cutmono2) go to 200
                        rxyz=mono_cent(3,i) - mono_cent(3,j)
                        r2 = r2 + rxyz*rxyz

c check distances^2. If smaller than cutmono2, include pair
c i,j in list, update pointers, go for next j.
c no need for extra cutoffs for waters. For waters use cutebig2 only
c
            watflag2 = watflag .and. ( moname(realmono(j)).eq.'TIP3'
     1          .or. (moname(realmono(j))(1:3).eq.'SPC') )

                        if (watflag2) then
                         if (r2.le.cutvbig2) then
                          nwlist1 = nwlist1 + 1
                          listwt1(nwlist1) = j
                         else if (r2.le.cutebig2) then
                          nwlist2 = nwlist2 + 1
                          listwt2(nwlist2) = j
                         end if

                        else if (r2.le.cutmono2) then

                         nlist=nlist+1
                         nblistm(nlist)=j
                         nflag = .false.

                        end if

200             continue
100     continue
        


         if (prll_on_off) then
          if (nflag) then
             pointm(monp(my_pe+1)+1) = nlist
          else
            pointm(monp(my_pe+1)+1) = nlist+1
          endif
          my_nwat = iwat
          poinwat1(iwat+1) = nwlist1
          poinwat2(iwat+1) = nwlist2
         else
          if (nflag) then
              pointm(totdmon+1) = nlist
           else
              pointm(totdmon+1) = nlist+1
           endif
           my_nwat = nwaters-1
        end if

        if (nwlist1.gt.wtrshrt) then
                write(*,*)' nwlist1 wtrshrt ',nwlist1,wtrshrt
                call alert(name,namel,'Water list1 too short',21,1)
        else if (nwlist2.gt.wtrlng) then
                write(*,*)' nwlist2 wtrshrt ',nwlist2,wtrshrt
                call alert(name,namel,'Water list2 too short',21,1)
        end if

 
        call nbond(pointm,nblistm)

c -------------------------------------------

        if (debug) write(*,*) ' leaving nbonm npt = ',npt
        return
        end

