        subroutine nbond(pointm,nblistm)
        implicit none

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/FREEZ.BLOCK'
        include 'COMMON/SPECL.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/DYNA.BLOCK'

        integer pointm(*),nblistm(*)
        double precision rx,ry,rz,r2
C Check for EVERY monomer.

        logical flag2
        integer i,j,k,l,m,ki,kj,jbeg,jend
        integer nlist1,nlist2,nlist3
        integer jx
        integer jbegx,jendx

C# yael - new variables for parallelization

        integer totng(0:maxpe)
        integer nforproc,lstart,lend
        integer p1num,p2num,p3num,temp,count,iwat
        integer nbwatf
        integer namel
        character*5 name
        data namel,name/5,'nbond'/

        save namel,name
        
        nlist1 = 0
        nlist2 = 0
        nlist3 = 0

        nbwatf = 5
        if (prll_on_off) then
           lstart = monp(my_pe)+1
           lend = monp(my_pe+1)
        else
           lstart = 1
           lend = totdmon
        endif
       
        point1(dpoipt(lstart-1)) = 0
        point2(dpoipt(lstart-1)) = 0
        point3(dpoipt(lstart-1)) = 0
        !write(6,*)"SSS",lstart,lend

        do 100 i=lstart,lend
          ki = dpoipt(i-1)+1
          jbeg = pointm(i)
          jend = pointm(i+1) - 1

          do 200 l=ki,dpoipt(i)
             do 300 j=jbeg,jend
                k=nblistm(j)
                !write(6,*)"VVV",i,k,j
                if (i.ne.k) then
                 kj = dpoipt(k-1)+1
                else
                 kj = l+1
                end if

                do 400 m=kj,dpoipt(k)


c check if there are frozen particles at all
        if (inofrz.ne.npt) then
c skip test if both particles are frozen
         if (zerofrz(l).eq.0 .and. zerofrz(m).eq.0) go to 400
        end if

        flag2=((lesid(l).eq.lesid(m)).and.
     *               (cplbl(l).ne.cplbl(m))).and.
     *               (lesid(l).ne.0.and.lesid(m).ne.0)
                    if (flag2) goto 400
c ------------------------------------------------------------
                    if (repyes0) then
                     if((l.eq.imb1(1)) .or.(l.eq.imb1(2)).or.
     2                 (l.eq.imb1(3)) .or.(l.eq.imb1(4))) then
                      if((m.eq.imb2(1)) .or. (m.eq.imb2(2)).or.
     1                (m.eq.imb2(3)).or.(m.eq.imb2(4))) go to 400
                     else if((m.eq.imb1(1)) .or.(m.eq.imb1(2)).or.
     2                 (m.eq.imb1(3)) .or.(m.eq.imb1(4))) then
                      if((l.eq.imb2(1)) .or. (l.eq.imb2(2)).or.
     1                (l.eq.imb2(3)).or.(l.eq.imb2(4))) go to 400
                     end if
                    end if
c ------------------------------------------------------------------
                      jbegx=exc1(l-1)+1
                      jendx=exc1(l)
                      if (jbegx.le.jendx) then
                       do 1500 jx=jbegx,jendx
                        if(exc2(jx).eq.m) then
                                go to 400
                        end if
1500                   continue
                      end if
                        rx =coor(1,l)-coor(1,m)
                        ry =coor(2,l)-coor(2,m)
                        rz =coor(3,l)-coor(3,m)
                        r2 = rx*rx + ry*ry + rz*rz
c
c check for hydrogens. Hydrogens have zero van-der Waals rad if arith=.false.
c (default) and should be placed in the third (electrostatic only)
c list
c
c now the decision is made depending on use of shake

        if (hvdw0) then
         flag2 = (.not.arith) 
     2    .and. (r2.lt.cutebig2)
     1    .and. (epsgm12(l).lt.1.d-3 .or. epsgm12(m).lt.1.d-3) 
         if (flag2 .and. flagchr(l) .and. flagchr(m)) then
                  nlist3 = nlist3 + 1
                  list3(nlist3) = m
                  go to 400
         else if (flag2) then
                go to 400
         end if
        end if

c make now a decision to which of the three
c lists this pair belongs: list1 includes both
c van der Waals and electrostatic and the max.
c distance is cutvbig2. list2 uses the same
c cutoff (cutvbig2) but intends to uncharged particle only.
c list3 is for for distances larger than cutvdw2 and smaller
c than cutebig2 and includes only electrostatic.

                if(r2.lt.cutvbig2)then
                  if (flagchr(l).and.flagchr(m)) then
                   nlist1 = nlist1 + 1
                   list1(nlist1) = m
                  else
                   nlist2 = nlist2 + 1
                   list2(nlist2) = m
                  end if
                else if (r2.lt.cutebig2) then
                  if (flagchr(l).and.flagchr(m)) then
                   nlist3 = nlist3 + 1
                   list3(nlist3) = m
                  end if
                end if
                  
400             continue
300          continue
             point1(l)=nlist1
             point2(l)=nlist2
             point3(l)=nlist3
200     continue
100     continue


        nblists(1,0) = nlist1
        nblists(2,0) = nlist2
        nblists(3,0) = nlist3

        if (nlist1.gt.ichgvdw) then
                write(*,*)' nlist1 ichgvdw ',nlist1,ichgvdw
                call alert(name,namel,'List1 too short',15,1)
         else if (nlist2.gt.ivdw) then
                write(*,*)' nlist2 ivdw ',nlist2,ivdw
                call alert(name,namel,'List2 too short',15,1)
         else if (nlist3.gt.ichg) then
                write(*,*)' nlist3 ichg ',nlist3,ichg
                call alert(name,namel,'List3 too short',15,1)
         end if
         

       if (my_nwat.ne.0) then
          nblists(4,0) = (poinwat1(my_nwat+1) +
     1           poinwat2(my_nwat+1))*nbwatf
       else
          nblists(4,0) = 0
       endif

        if (.not.prll_on_off) return
c-------------------- yael---------------------------
c divide the work for the next iteration

c send the total number of neighbors (for this process) 
c
          
           call gather_nblists(nblists)    

           nblists(1,num_pes) = 0
           nblists(2,num_pes) = 0
           nblists(3,num_pes) = 0
           nblists(4,num_pes) = 0
           
           totng(0) = nblists(1,0) + nblists(2,0) +
     1           nblists(3,0) + nblists(4,0)
           do i=1,num_pes
              totng(i) = nblists(1,i)+nblists(2,i)+nblists(3,i)
     1           + nblists(4,i) + totng(i-1)
           end do   

           do i=0,num_pes-1
              nblists(1,num_pes) = nblists(1,num_pes) + nblists(1,i)
              nblists(2,num_pes) = nblists(2,num_pes) + nblists(2,i)
              nblists(3,num_pes) = nblists(3,num_pes) + nblists(3,i)
              nblists(4,num_pes) = nblists(4,num_pes) + nblists(4,i)
              new_monp(i) = 0
           end do   

           nforproc = totng(num_pes) / num_pes
          
           if (nforproc .gt. 0) then

           if (my_pe.eq.0) then
              count = 0
              temp = 0
           else
              count = int(totng(my_pe-1) / nforproc)
              temp =  totng(my_pe-1) - (nforproc*count)
           endif   
c
          do 610 j=monp(my_pe)+1,monp(my_pe+1)
              if (moname(realmono(j)).eq.'TIP3' .or.
     1             moname(realmono(j))(1:3).eq.'SPC') then
                 iwat = watptr(j)
                 p1num = (poinwat1(iwat+1) - poinwat1(iwat))*nbwatf
                 p3num = (poinwat2(iwat+1) - poinwat2(iwat))*nbwatf
                 temp = temp + p1num + p3num
                 if (temp.ge.nforproc) then
                    new_monp(count+1) = j
                    temp = temp - nforproc
                    count = count + 1
                 endif
             endif     
             do 605 i=dpoipt(j-1)+1,dpoipt(j)
             p1num = point1(i) - point1(i-1)
             p2num = point2(i) - point2(i-1)
             p3num = point3(i) - point3(i-1)
             temp = temp + p1num + p2num + p3num
             if (temp.ge.nforproc) then
                   new_monp(count+1) = j
                   temp = temp - nforproc
                   count = count + 1
             endif
605        continue
610     continue
c            
c sum the work vector for all processes
             call reduce_int(new_monp,num_pes)  
             do i=0,num_pes+1
               !write(6,*)"XXX",i,new_monp(i)
             end do
         else
           write(6,*) "ERROR!!!"
           stop
         end if ! numproc .eq. 0


        return
        end
