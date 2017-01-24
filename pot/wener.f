        subroutine wener(uwene)
        implicit none 
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/SPECL.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/SSBP.BLOCK'
        include 'COMMON/SYMM.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/EWALD.BLOCK'
        include 'COMMON/METAL.BLOCK'
        include 'COMMON/SGB.BLOCK'
        integer uwene
        double precision gradff
        integer i,j,n
        integer k1,k2,k3
        integer symidx,istart
        integer jbeg1,jbeg2,jbeg3,jend1,jend2,jend3,nbwatf

        nbwatf = 5

        
       if (num_pes .gt. 1) call reduce_energies()
        call force_norm(gradff)

c --------------------------------------------------

C        if (my_pe.eq.0) then
cc@yef
c       do 999 i=1,npt
c               write(*,*)i,(dpot(j,i),j=1,3)
c999    continue
         write(uwene,100)e_total,e_bond,e_theta,e_tors,e_imp,e_vdw,
     1e_el,e_el14,e_vdw14,e_cnst,e_vsym,e_lsym,e_cent,e_hyd
100      format(1x,/,' ENERGIES: E total = ',f12.3,/,
     11x,'E bond = ',f10.3,3x,'E angl = ',f10.3,3x,
     2'E tors = ',f10.3,/,1x,'E impr = ',f10.3,3x,
     3'E vdw  = ',f10.3,3x,'E elec = ',f10.1,/,1x,
     4'E 14el = ',f10.3,3x,'E 14vd = ',f10.3,/,1x,
     5'E cnst = ',f10.3,3x,'E evsym= ',f10.3,3x,
     6'E elsym= ',f10.3,/,1x,'E centr= ',f10.3,3x,
     7'E hydro= ',f10.3)

c YS
         if (gbsabool) then
             write(uwene,111) e_gbsa,e_gbsa - e_gbnonpol,e_gbnonpol
         end if
111      format(1x,'E gbsa =', f10.3, /,
     &         1x,'E pol =', f10.3, 3x,'E nonpol =', f10.3)

         if (qssbp) then
         write(uwene,101) enssbp,enkir,encav,enhsr,enangu,enharm
101      format(1x,'E ssbp = ',f10.3,3x,'E kirk = ',f10.3,3x,
     1'E cavi = ',f10.3,/,1x,
     2'E hsr  = ',f10.3,3x,'E angu = ',f10.3,3x,
     3'E harm = ',f10.3)
         end if

         if (ewaldyes) then
         write(uwene,*) 
         write(uwene,*) ' EWALD-type splitting of electrostatics'
         write(uwene,*) ' note - E elec and E elsym do not',
     1 ' contain recip. and corr. contr.'
         write(uwene,117) e_dir+e_ew_receip+e_self+e_corr,
     1                    e_dir,e_ew_receip,e_self,e_corr
         write(uwene,118) fpme
117      format(1x,/,' E elec TOTAL (incl. elsym) = ',f15.3,/,1x,
     1'E dir  = ',f12.3,3x,
     2'E recip= ',f12.3,/,1x,'E self = ',f12.3,3x,
     3'E corr = ',f12.3)
118      format(' Receip (without corrections) norm. force = ',f10.4)
         write(uwene,*) 
         end if

         if (metalyes) then
          write(uwene,119)e_wall,e_lctd,e_lmet
119       format(/,1x,' Metal E-s : ',/,1x,'E Wall = ',f10.3,3x,
     1'E Elect= ',f10.3,3x,'E image= ',f10.3,//)
         end if

        if (eballyes) then
         write(uwene,120)e_ball
120      format(/,1x,' E Solvation BALL = ',f10.3,/)
        end if

        if (esteeryes) then
         write(uwene,121)e_steer
121      format(/,1x,' E steer ',f10.3,/)
        end if
         write(uwene,102)gradff
102      format(1x,' Norm Force = ',f10.3)
        
c ---------------------------------------------------------
        if (eteth_yes) write(uwene,103)e_tether
103     format(1x,'E tether = ',f10.3)

c -----------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do 888 n=1,nmb
         if(repyes(n)) then
          write(uwene,661)
661       format(1x, ' repulsion was used ')
          write(uwene,777) e_repuls(n)
777       format(1x,' E repulsion = ',f10.3,/,1x)
          else if (emyes(n) .and. nmb .gt. 0) then
          write(uwene,669)
669       format(1x, ' Morse potential was used ')
           write(uwene,677) e_morseb(n)
677        format(1x,' E morse  = ',f10.3,/,1x)
          endif
888      continue
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C WRITE DOWN DATA OF SYMMETRY
c@
        if (esymyes) then
          if (prll_on_off) then
            k1 = nblists(1,my_pe)
            k2 = nblists(2,my_pe)
            k3 = nblists(3,my_pe)
          else
            k1 = point1(npt-1)
            k2 = point2(npt-1)
            k3 = point3(npt-1)
          end if
        
         write(uwene,104)indxsym1,indxsym2,indxsym3
104      format(//,1x,' Total number of symmetry related particles
     1  (symm)',/,1x,' short range:',i7,/,1x,' Uncharged vdw:',i7
     2         ,/,1x,' electrostatic: ',i7)
         if (iblock1(symnum).eq.0) then
          write(uwene,105) 0
         else
          write(uwene,105)psym1(iblock1(symnum))-k1
         end if
105      format(1x,' Number of symmetry neighbours (short rang) ',i7)
         if (iblock2(symnum).eq.0) then
          write(uwene,106) 0
         else
          write(uwene,106)psym2(iblock2(symnum))-k2
         end if
106      format(1x,' Number of symmetry neighbours: van der Waals ',i7)
         if (iblock3(symnum).eq.0) then
          write(uwene,107)0
         else
          write(uwene,107)psym3(iblock3(symnum))-k3
         end if
107      format(1x,' Number of symmetry neighbours: electrostatic ',i7)

         if (debug) then
         do 3 i=1,symnum
          if (iblock1(i).eq.0) then
           j = 0
          else
           j = psym1(iblock1(i)) - k1
          end if
          write(uwene,108)symop(1,i),symop(2,i),symop(3,i),j
          if (j.ne.0) k1 = psym1(iblock1(i))
108       format(1x,' symmetry operation (short range) ',3(i2,1x),
     1    'number of particle pairs ',i7)
          if (iblock2(i).eq.0) then
           j = 0
          else
            j = psym2(iblock2(i)) - k2
          end if
          write(uwene,109)symop(1,i),symop(2,i),symop(3,i),j
          if (j.ne.0) k2 = psym2(iblock2(i))
109       format(1x,' symmetry operation (van der Waals) ',3(i2,1x),
     1    'number of particle pairs ',i7)
          if (iblock3(i).eq.0) then
           j = 0
          else
           j = psym3(iblock3(i)) - k3
          end if
          write(uwene,110)symop(1,i),symop(2,i),symop(3,i),j
          if (j.ne.0) k3 = psym3(iblock3(i))
110       format(1x,' symmetry operation (electrostatic) ',3(i2,1x),
     1    'number of particle pairs ',i7)
3        continue
        end if
        end if
c
c write down information on water symmetry
c
        if (esymyes.and.nwaters.gt.0) then
          if (prll_on_off) then
            k1 = poinwat1(my_nwat+1)
            k2 = poinwat2(my_nwat+1)
          else
            k1 = poinwat1(nwaters)
            k2 = poinwat2(nwaters)
          end if
        
         write(uwene,202)ipwat1,ipwat2
202      format(//,1x,' Total number of symmetry waters
     1  (symm)',/,1x,' van der Waals range: ',i7
     2         ,/,1x,' electrostatic only : ',i7)
         if (iblckwt1(symnum).eq.0) then
          write(uwene,203) 0
         else
          write(uwene,203)psymwt1(iblckwt1(symnum))-k1
         end if
203      format(1x,' Num. of symm water neighbours (shrt rng) ',i7)
         if (iblckwt2(symnum).eq.0) then
          write(uwene,204) 0
         else
          write(uwene,204)psymwt2(iblckwt2(symnum))-k2
         end if
204      format(1x,' Num. of symm water neighbours: (electro) ',i7)

        if (debug) then
         do 4 i=1,symnum
          if (iblckwt1(i).eq.0) then
           j = 0
          else
           j = psymwt1(iblckwt1(i)) - k1
          end if
          write(uwene,206)symop(1,i),symop(2,i),symop(3,i),j
206       format(1x,' symmetry operation (short range) ',3(i2,1x),
     1    'number of particle pairs ',i7)
          if (j.ne.0) k1 = psymwt1(iblckwt1(i))
          if (iblckwt2(i).eq.0) then
           j = 0
          else
            j = psymwt2(iblckwt2(i)) - k2
          end if
          write(uwene,208)symop(1,i),symop(2,i),symop(3,i),j
208       format(1x,' symmetry operation (van der Waals) ',3(i2,1x),
     1    'number of particle pairs ',i7)
          if (j.ne.0) k2 = psymwt2(iblock2(i))
4        continue
        end if
        end if

         if (nwaters.eq.0) then
                k1 = 0
                k2 = 0
         else 
                k1 = poinwat1(nwaters)*nbwatf
                k2 = poinwat2(nwaters)*nbwatf
         end if
         if (prll_on_off) then 
           write(uwene,209) nblists(1,num_pes),nblists(2,num_pes),
     1          nblists(3,num_pes),nblists(4,num_pes)
209     format(//,1x,' Number of neighbours for short range int.'
     1 ,i10,/,1x, ' Number of uncharged vdW interactions     '
     2 ,i10,/,1x, ' Number of elec. only interactions        '
     3 ,i10,/,1x, ' Number of wat-wat neighbors  '
     4 ,i10)
         else
           write(uwene,210)point1(npt-1),point2(npt-1),
     1          point3(npt-1),k1,k2
         endif   
210     format(//,1x,' Number of neighbours for short range int.'
     1 ,i10,/,1x, ' Number of uncharged vdW interactions     '
     2 ,i10,/,1x, ' Number of elec. only interactions        '
     3 ,i10,/,1x, ' Number of wat-wat shrt. range neighbors  '
     4 ,i10,/,1x, ' Number of wat-wat long  range neighbors  '
     5 ,i10)
C        end if

        if (debug) then
        write(uwene,*)'dx,dy,dz'
        do 9 i=1,npt
         write(uwene,10) i,dpot(1,i),dpot(2,i),dpot(3,i)
9       continue
        jbeg1 = 1
        jbeg2 = 1
        jbeg3 = 1
        do 11 i=1,npt-1
        jend1 = point1(i)
        jend2 = point2(i)
        jend3 = point3(i)
        write(uwene,*)' i list1 = ',i,(list1(k1),k1=jbeg1,jend1)
        write(uwene,*)' i list2 = ',i,(list2(k1),k1=jbeg2,jend2)
        write(uwene,*)' i list3 = ',i,(list3(k1),k1=jbeg3,jend3)
        jbeg1 = jend1 + 1
        jbeg2 = jend2 + 1
        jbeg3 = jend3 + 1
10      format(1x,/,i4,3f12.5)
11      continue
        jbeg1 = point1(npt-1) + 1
        istart = 1
        do 13 symidx=1,symnum
         do 12 i=istart,iblock1(symidx)
          jend1 = psym1(i)
           write(uwene,*)' ireal ',symreal1(i)
           write(uwene,*)' list ',(list1(k1),k1=jbeg1,jend1)
12       continue
         jbeg1 = jend1 + 1
         istart = iblock1(symidx) + 1
13      continue
        jbeg2 = point2(npt-1) + 1
        istart = 1
        do 15 symidx=1,symnum
         do 14 i=istart,iblock2(symidx)
          jend2 = psym2(i)
           write(uwene,*)' ireal ',symreal2(i)
           write(uwene,*)' list ',(list2(k1),k1=jbeg2,jend2)
14       continue
         jbeg2 = jend2 + 1
         istart = iblock2(symidx) + 1
15      continue
        jbeg3 = point3(npt-1) + 1
        istart = 1
        do 17 symidx=1,symnum
         do 16 i=istart,iblock3(symidx)
          jend3 = psym3(i)
           write(uwene,*)' ireal ',symreal3(i)
           write(uwene,*)' list ',(list3(k1),k1=jbeg3,jend3)
16       continue
         jbeg3 = jend3 + 1
         istart = iblock3(symidx) + 1
17      continue

        end if
        return
        end
