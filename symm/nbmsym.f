      subroutine nbmsym(ia,jb,kc,lsym)
      implicit none
c     
c     NOTE: the generation of the non-bonded list is optimized for
c     the use of TIP3/SPC/SPCE (waters). Any complains should be directed
c     to /dev/null . You were told (just now)!

      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/NBLIST.BLOCK'
      include 'COMMON/SYMM.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/PARALLEL.BLOCK'
      integer nlist
      integer ia,jb,kc,lsym

C-------local variables--------------------------

      integer i,j,k,total_need_only
      double precision rx,r2
      double precision da,db,dc,ha,hb,hc
      double precision xch,ych,zch,rch2
      double precision barx,bary,barz,closex,closey,closez
      double precision large_cutmono
      double precision mono_csym(3,maxmono*maxdiv)
c     
c     not all molecules in the box need to be considered. Here we prepare a subset of them.
c     need_only is a pointer to these molecules
c     
      integer need_only(maxmono*maxdiv)
      integer pointm(0:maxmono*maxdiv+1)
      integer nblistm(ihugemo),my_part
      logical watflag,oneig1,oneig2
      logical lbarx,lbary,lbarz
      integer namel
      character*6 name

      logical metal_now

c     yael - variables for parallelization
      
      integer lstart,lend,temp,count,nforproc,watpt1,watpt2
      integer watptr1(maxmono*maxdiv),watptr2(maxmono*maxdiv)
      integer nbmlists(0:maxpe),tempvec(0:maxpe)
      integer sflag,lstart2,lend2
      save sflag
      save name,namel

      data sflag/1/
      data name,namel/'nbmsym',6/

      if (.not. (evdyes.or.eelyes)) return

      do i=0,maxpe
        tempvec(i)=0
      end do

      metal_now = .false.
c      write(6,*)' cutmono2 = ',cutmono2
      large_cutmono   = cutmono2
c     
c     get box side locations for use in loop
      ha=0.5*a
      hb=0.5*b
      hc=0.5*c

c     
c     set symmetry translations
c     
      da = ia*a
      db = jb*b
      dc = kc*c

c     the REAL monomers that need to be considered must be at least cutmono
c     distance from the edges of the box.
c     
      barx=0.5d0*da
      bary=0.5d0*db
      barz=0.5d0*dc
      lbarx = dabs(barx).gt.0
      lbary = dabs(bary).gt.0
      lbarz = dabs(barz).gt.0
c     write(6,*)' barx bary barz ', barx,bary,barz
      do 330 i=1,totdmon
         if (dabs(mono_cent(1,i)).gt.ha .or.
     $        dabs(mono_cent(2,i)).gt.hb .or.
     $        dabs(mono_cent(3,i)).gt.hc) then
            large_cutmono = cutmono2*4.d0
            go to 335
         end if
 330  continue
 335  continue

      total_need_only = 0
      do 340 i=1,totdmon
         if (lbarx) then
            closex = mono_cent(1,i)-barx
            if (closex*closex.gt.cutmono2) go to 340
c         else if (lbary) then
         end if
         if (lbary) then
            closey = mono_cent(2,i)-bary
            if (closey*closey.gt.cutmono2) go to 340
c         else if (lbarz) then
         end if
         if (lbarz) then
            closez = mono_cent(3,i)-barz
            if (closez*closez.gt.cutmono2) go to 340
         end if
         total_need_only = total_need_only + 1
         need_only(total_need_only) = i
 340  continue

C     write(6,*)'total, # of real monomers considered = ',totdmon,total_need_only
c     create the initial partition of monomers between the processors
c     or copy the partition calculated at the last step

      if (prll_on_off) then
         if (sflag.le.symnum) then
            temp = totdmon /num_pes
            do i=0,num_pes-1
               monsym(lsym,i) = temp*i
            end do 
            sflag = sflag + 1
         else
            do i=0,num_pes-1
               monsym(lsym,i) = new_monsym(lsym,i)
 360        end do
         end if
         monsym(lsym,num_pes) = totdmon
      end if
c     nlist       = index for monomer non-bonding list
      nlist = 0
      watflag = .false.
      
c     loop over real system monomers to create the relevant symmetry 
c     particles
c     
      if (metalyes .and. lsym.gt.4) metal_now = .true.
      if (metal_now) then
         do 1 i=1,totdmon
            mono_csym(1,i) = mono_cent(1,i) + da
            mono_csym(3,i) = mono_cent(3,i) + dc
            if (mono_cent(2,i).gt.0.d0) then
               mono_csym(2,i) = b - mono_cent(2,i)
            else
               mono_csym(2,i) = - mono_cent(2,i) - b
            end if
 1       continue
      else
         do 2 i=1,totdmon
            mono_csym(1,i) = mono_cent(1,i) + da
            mono_csym(2,i) = mono_cent(2,i) + db
            mono_csym(3,i) = mono_cent(3,i) + dc
 2       continue
      end if

      if (prll_on_off) then
         lstart = monsym(lsym,my_pe) + 1
         lend = monsym(lsym,my_pe+1)
      else
         lstart = 1
         lend = totdmon
      endif
      pointm(lstart-1) = 0


      do 10 i=lstart,lend
c     
c     check if the translated monomer is sufficiently closed to the box
c     
         watptr1(i) = 0
         watptr2(i) = 0

         xch = abs(mono_csym(1,i)) - ha
         ych = abs(mono_csym(2,i)) - hb
         zch = abs(mono_csym(3,i)) - hc

         if (xch.lt.0.d0) xch=0.d0
         if (ych.lt.0.d0) ych=0.d0
         if (zch.lt.0.d0) zch=0.d0
         rch2=xch*xch+ych*ych+zch*zch
         if (rch2.gt.large_cutmono) then
            pointm(i) = nlist
            go to 10
         endif
         watflag = .false.
         if (moname(realmono(i)).eq.'TIP3' .or.
     1       moname(realmono(i))(1:3).eq.'SPC') then 
            watflag = .true.
            oneig1  = .false.
            oneig2  = .false.
         end if

c     loop over second (real) monomers and compute distances
         do 6 k=1,total_need_only
            j = need_only(k)          
            rx = mono_csym(1,i)-mono_cent(1,j)
            r2 = rx*rx
            if (r2.gt.cutmono2) go to 6
            rx = mono_csym(2,i)-mono_cent(2,j)
            r2 = r2 + rx*rx
            if (r2.gt.cutmono2) go to 6
            rx = mono_csym(3,i)-mono_cent(3,j)
            r2 = r2 + rx*rx
c     
c     special treatment for water-water interactions
c     water - water are stored directly in monomer list
c     Metal charge images should be stored in the list for charged
c     interactions only
c     
            if (watflag.and. (moname(realmono(j)).eq.'TIP3' .or.
     1          moname(realmono(j))(1:3).eq.'SPC') ) then
               if (r2.lt.cutvbig2 .and. (.not.metal_now)) then
                  oneig1  = .true.
                  nwlist1 = nwlist1 + 1
                  listwt1(nwlist1) = j
               else if (r2.lt.cutebig2) then
                  oneig2  = .true.
                  nwlist2 = nwlist2 + 1
                  listwt2(nwlist2) = j
               end if

c     check distances. If smaller than cutmono2, include pair
c     i,j in list, update pointers, go for next j.
c     
            else if (r2.lt.cutmono2) then
               nlist=nlist+1
               nblistm(nlist)=j
            end if

            
 6       continue
         pointm(i) = nlist
         if (watflag) then
            if (oneig1) then
               ipwat1 = ipwat1 + 1      
               symmwat1(ipwat1) = i
               watptr1(i) = ipwat1
            end if
            if (oneig2) then
               ipwat2 = ipwat2 + 1
               symmwat2(ipwat2) = i
               watptr2(i) = ipwat2
            end if
            if (ipwat1.ne.0) psymwt1(ipwat1) = nwlist1
            if (ipwat2.ne.0) psymwt2(ipwat2) = nwlist2
         end if
 10   continue

      
      iblckwt1(lsym) = ipwat1
      iblckwt2(lsym) = ipwat2

      if (nwlist1.gt.wtrshrt) then
         write(*,*)' nwlist1 wtrshrt ',nwlist1,wtrshrt
         call alert(name,namel,'Water list1 too short',21,1)
      else if (nwlist2.gt.wtrlng) then
         write(*,*)' nwlist2 wtrlng ',nwlist2,wtrlng
         call alert(name,namel,'Water list2 too short',21,1)
      end if

c----------------------------------------------------------------
c     divide to monomer list between the processors.
c     first, send my monomer list size to all other processors

      if (prll_on_off) then

         lstart = iblckwt1(lsym-1)
         lend = iblckwt1(lsym)
         lstart2 = iblckwt2(lsym-1)
         lend2 = iblckwt2(lsym)

         nbmlists(0) = nlist + psymwt1(lend) - psymwt1(lstart) +
     $        psymwt2(lend2) - psymwt2(lstart2)
         call gather_nbmlists(nbmlists)
         nbmlists(num_pes) = 0

C       write(6,*)'load balancing for symmetry operation num ',lsym,' :'
C       do i=0,num_pes-1
C         write(6,*)'pes ',i,' has ',nbmlists(i)
C       end do 

         do i=1,num_pes
            nbmlists(i) = nbmlists(i) + nbmlists(i-1)
            tempvec(i) = 0
         end do

         tempvec(0) = 0

         nforproc = nbmlists(num_pes-1) / num_pes

         if (nforproc.ne.0) then
            if (my_pe.eq.0) then
               count = 0
               temp = 0
            else
               count = int(nbmlists(my_pe-1) / nforproc)
               temp =  nbmlists(my_pe-1) - (nforproc*count)
            endif

            do 380 i=monsym(lsym,my_pe)+1,monsym(lsym,my_pe+1)
               if (moname(realmono(i)).eq.'TIP3' .or.
     1             moname(realmono(i))(1:3).eq.'SPC' ) then
                  watpt1 = watptr1(i)
                  watpt2 = watptr2(i)
                  if (watpt1.ne.0) then
                     temp = temp + psymwt1(watpt1) - psymwt1(watpt1-1)
                  endif
                  if (watpt2.ne.0) then
                     temp = temp + psymwt2(watpt2) - psymwt2(watpt2-1)
                  endif
               endif
               temp = temp + (pointm(i) - pointm(i-1))
               if (temp.ge.nforproc) then
                  tempvec(count+1) = i
                  temp = temp - nforproc
                  count = count + 1
               endif
 380        continue

            my_part = int(nbmlists(my_pe) / nforproc)
            if (count.lt.my_part) then
               do 410 i=count,my_part-1
                  tempvec(i+1) = tempvec(i)+1
c     write(*,*)'!!',i+1,tempvec(i)+1
 410           continue
            endif

c     sum the work vector for all process
            call reduce_int(tempvec,num_pes)

         endif
      endif

      do i=0 , num_pes
         new_monsym(lsym,i) = tempvec(i)
      end do
C     @ ---------------------------------------------------

      if (metalyes .and. lsym.gt.4) then
         call nbsym_metal(pointm,nblistm,lsym,ia,jb,kc)
      else
         call nbsym(pointm,nblistm,lsym,ia,jb,kc)
      end if

      return
      end













