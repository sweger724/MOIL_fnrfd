        subroutine shakinit(shakl,shakb,shakm,epsilon)
c
c initiating the bond (to be shaked) pointer arrays for shake
c ishak1/2(k) is the number of the first/second particle in
c the k-th bond constraint. dissq holds the ideal distance values^2

c input : shakl - if true, bonds with light (m<1.1) particles are shaked
c         shakb - all bonds are shaked
c also shaked are light particles that participate in angles even if
c not bonded i.e.
c                              O
c                            /   \
c                           H.....H
c The H....H distance is shaked

c output: ishak1(nshak) ishak2(nshak) nshak dissq(nshak)
c         (in SHAKE.BLOCK)
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/SHAKE.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/FREEZ.BLOCK'
        logical shakl,shakb,shakm
        double precision epsilon
c local
        character*8 name
        integer namel,i,i1,j1,level
        integer ibond1
        logical shak_accpt_bnd

         if (nshak.gt.0) then
          level = 0
          call alert('shakinit',8,' Shakinit called twice',22,level)
          return
         end if
        if (shakb) then
                write(stdo,100)
100             format(/,1x,' All bonds will be shaked ',/)
            else if (shakl) then
                 write(stdo,101)
101   format(/,1x,' All bonds with light particles will 
     1   be shaked',/)
                else
                    level = 0
                    call alert(name,namel,
     *                     'routine called but not needed',29,level)
                    return
        end if
      if (shakm) write (6,*) 'shakm is on, not shaking TIP3/SPC/SPCE'
         do 1 i=1,nb_all
          if(shak_accpt_bnd(shakl,shakb,shakm,i)) then
             nshak = nshak + 1
             ishak1(nshak) = ib1(i)
             ishak2(nshak) = ib2(i)
             dissq(nshak)  = req(i)*req(i)
             shak_diff(nshak) = dissq(nshak)*epsilon
          endif
1        continue

         write (6,*) 'this gives a total of ',nshak,
     &   ' actual bond constraints'


c
c Now shaking the angle of two hydrogens connected to a shared atom
c by fixing the distance between the two hydrogens
c since the bonds have identical properties (two hydroegns connected to the
c same atom), only the properties of one of the bonds is extracted.
c

         do 15 i=1,nangl_all
          if ((ptms(iangl1(i)).lt.1.1d0 .and.
     1                  ptms(iangl3(i)).lt.1.1d0) .and.
     2  ((.not.shakm).or.(moname(poimon(iangl1(i))).ne.'TIP3')) .and.
     2 ((.not.shakm).or.(moname(poimon(iangl1(i)))(1:3).ne.'SPC')) .and.
     3  (.not.(zerofrz(iangl1(i)).eq.0.and.zerofrz(iangl3(i)).eq.0)))
     4           then
           nshak     = nshak + 1 
           i1        = iangl1(i)
           j1        = iangl3(i) 
           ishak1(nshak) = i1
           ishak2(nshak) = j1
           call find_bond(iangl2(i),i1,ibond1)
           dissq(nshak)  = 2.d0*req(ibond1)*sin(angleq(i)*0.5d0)
           dissq(nshak)  = dissq(nshak)*dissq(nshak)
           shak_diff(nshak) = dissq(nshak)*epsilon
          end if
15       continue
         write(stdo,*)' Number of shake constraints = ',nshak
         if (nshak.gt.maxshak) then
          level = 1
          call alert(name,namel,'Max. # of constraints exceed',28,level)
         end if
c chen
         pt_start = 1
         pt_end = npt

!         if (prll_on_off) then 
c gather the full list of shaked particles 
!            call gather_shake(nshak,ishak1,ishak2,dissq,shak_diff)          
!            if (matshak) then 
!               call shake_balance()
!            else
!               sk_start = 1
!               sk_end = nshak
!            endif
!         else
            sk_start = 1
            sk_end = nshak
!         endif
C@
C        write(*,*)'!!!!!! ',pt_start,pt_end,matshak
         return
        end
        
c--------------------------------------------------------------------------
        logical function shak_accpt_bnd(shakl,shakb,shakm,i)
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/SHAKE.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/FREEZ.BLOCK'
        logical shakl,shakb,shakm
        integer i

        if(zerofrz(ib1(i)).eq.0 .or. zerofrz(ib2(i)).eq.0) then
                shak_accpt_bnd = .false.
                return
        end if
        if ((shakm) .and. (moname(poimon(ib1(i))).eq.'TIP3' .or.
     1     moname(poimon(ib1(i)))(1:3).eq.'SPC')) then
                shak_accpt_bnd = .false.
                return
        end if
        if (lesid(ib1(i)).ne.lesid(ib2(i))) then
                shak_accpt_bnd = .false.
                return
        end if
        if ((shakl) .and. 
     *      (ptms(ib1(i)).ge.1.1d0 .and. ptms(ib2(i)).ge.1.1d0)) then
                shak_accpt_bnd = .false.
                return
        else if ((shakl) .and. 
     *   (ptnm(ib1(i))(1:1).eq.'H'.or.ptnm(ib2(i))(1:1).eq.'H')) then
                shak_accpt_bnd = .true.
                return
        end if
        shak_accpt_bnd = .true.
        return
        end

c-------------------------------------------------------------------------
        subroutine shake_balance()

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/SHAKE.BLOCK'
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        
        integer lastmon,startmon        
        integer n_count,i,tnshak,istart,iend

        tnshak = nshak
        call load_balance(tnshak, my_pe,num_pes,istart,iend)

        startmon = poimon(ishak1(istart))
        lastmon = poimon(ishak1(iend))
        if (my_pe.eq.num_pes-1) lastmon = poimon(ishak2(iend))
c set the start end end of the bonds list in start end end of monomers

        if ((istart.ne.1).and.(poimon(ishak1(istart-1)).eq.
     1   startmon)) then
           do 40 i=istart,iend
              if (poimon(ishak1(i)).ne.startmon) then
                 istart = i
                 startmon = poimon(ishak1(istart))
                 goto 5
              endif
 40        continue
        endif  
        
5       continue

        shake_bal(0) = istart
        call gather_shake_balance(shake_bal)
        shake_bal(num_pes) = nshak + 1

        if (my_pe.ne.0) pt_start = ishak1(shake_bal(my_pe))
        if (my_pe.ne.num_pes-1) pt_end = ishak1(shake_bal(my_pe+1))-1

        sk_start = shake_bal(my_pe)
        sk_end = shake_bal(my_pe+1)-1
        
        shared_pt_num = 0

        do 70 i=sk_start,sk_end
           if ((ishak2(i).lt.pt_start).or.(ishak2(i).gt.pt_end)) then
              shared_pt_num = shared_pt_num + 1
              shared_pt(shared_pt_num) = ishak2(i)
              write(*,*)'!!',shared_pt_num,shared_pt(shared_pt_num)
           endif
 70     continue

        n_count = (pt_end - pt_start + 1)*3
        write(*,*)'I have',pt_start,pt_end,pt_end-pt_start
        call gather_int(n_count,nptsh)

        disp(0) = 0
        do 60 i=1,num_pes
           disp(i) = disp(i-1) + nptsh(i-1)
 60     continue

        return 
        end
