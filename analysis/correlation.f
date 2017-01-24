      program correlation
      implicit none

C sizes of things      
      include 'COMMON/LENGTH.BLOCK'
C need it for coor() array
      include 'COMMON/COORD.BLOCK'
C line intrpreter
      include 'COMMON/LINE.BLOCK'
C npt, etc
      include 'COMMON/CONNECT.BLOCK'

      include 'COMMON/CONVERT.BLOCK'

      integer i1,i,j,k,l
      integer stdi,urdyn
      integer namel,urcrd,strb,stre,size,uwdat,ucon,strd
      integer thismol,thatmol

      character*4 name
      character*4 cstyl

      double precision dotme
      double precision step
      double precision r(3,200000)
      double precision rdyn(3,200000)
      double precision corr(20000),normalize,getd
      double precision timea(20000)

C Routines/externals
      integer of,geti,igrid,mode,skno,count,urtim
      double precision avg,dx,dy,dz,avgvel,avgvel2,avgtime
      logical find
      

      jnkf=25
      stdi=5
      namel=4
      normalize = 1.0d0
      name = 'corr'


C Read our input file
      open(unit=jnkf,status='scratch')

 1    continue

      call rline(name,namel,stdi)

      strb = geti('strb',strb)
      stre = geti('stre',stre)
      mode = geti('mode',mode)

      skno = geti('skno',skno)
      step = getd('step',step)

      strd = geti('strd',strd)


      if (find('file')) then
         if (find('conn')) ucon = of()
         if (find('rcrd')) urcrd = of()
         if (find('rdyn')) urdyn = of()
         if (find('wdat')) uwdat = of()
C         if (find('rtim')) urtim = of()
      end if
      if (find('cpth')) cstyl = 'PATH'
      if (find('cdyn')) cstyl = 'DYNA'



      if (find('acti')) goto 2

      goto 1  
 2    continue

      close(jnkf)



C ---- Done reading input file -- start code -------------
C     Get connectivity

      call rconn(ucon)


C Now read coordinates
      
      if (cstyl.eq.'DYNA') then
         call rchain(urcrd,r,stre,'DYNA',0)
      else if (cstyl.eq.'PATH') then
         call rchain(urcrd,r,stre,'PATH',0)
      end if



      if (mode.eq.0) then 
C Now read in our times
c      avgtime = 0
c      do i = 1,stre
c         read(urtim,*) count,timea(i)
c         write (6,*) 'times array ',count,timea(i)*tconv
c         avgtime = avgtime + timea(i)*tconv
c      end do
c      write (6,*) 'total time ',avgtime
c      avgtime = avgtime/stre
c      write (6,*) 'our avg time is ',avgtime
c      stop
C Now lets do our correlation----------------------------------

C     init our array
         
         size = (stre-strb)/2

         do i = 1,size
            corr(i) = 0.0d0
         end do
         


c     (i = delta : <v(t+i)|v(t)> )
c     (j = first structure to last structure we can avg over (v(t) )

         do i = 1,size-1
            write (6,*) 'i ',i
            corr(i) = 0.0d0

C Now go over (size) structures (to make sure each entry has same
C number of entries for normalization
            count = 0
            do j = strb,stre-i

               count = count+1
               thismol = (j-1)*npt
               thatmol = (j+(i-1)-1)*npt
c     dot this molecule

               dotme = 0.0d0

               do k = 1,npt

                  dotme = dotme + 
     $                 r(1,thismol+k)*r(1,thatmol+k) + 
     $                 r(2,thismol+k)*r(2,thatmol+k) + 
     $                 r(3,thismol+k)*r(3,thatmol+k) 

               end do
               
               corr(i) = corr(i) + dotme

c     go over all structures

            end do

c     normalize by our degree's of freedom and structures we've avged over
            write (6,*) 'j ',i,count
            corr(i) = corr(i)/(3*(npt*1.0d0)*count)
            write (uwdat,*) (i-1)*1.0d0*step,corr(i)

c     calculate for all i 

         end do
C Mode != 0 (i.e. correlation )

       else 

C Lets read in the dynamics bits too
      call rchain(urdyn,rdyn,strd,'PATH',0)

      do i = 1,strd

         do j = 1,npt
            avgvel = avgvel + abs(rdyn(1,(i-1)*npt+j))
            avgvel = avgvel + abs(rdyn(1,(i-1)*npt+j))
            avgvel = avgvel + abs(rdyn(1,(i-1)*npt+j))
         end do
      end do

      count = 1
      avgvel2 = 0
      do i=1,strd,skno+1
         avg = 0
         
         do j=1,npt
            dx = r(1,(count-1)*npt+j) - rdyn(1,(i-1)*npt+j)
            dy = r(2,(count-1)*npt+j) - rdyn(2,(i-1)*npt+j)
            dz = r(3,(count-1)*npt+j) - rdyn(3,(i-1)*npt+j)
            avg = avg + dx*dx + dy*dy + dz*dz
            avgvel2 = avgvel2 + abs(r(1,(count-1)*npt+j))
            avgvel2 = avgvel2 + abs(r(2,(count-1)*npt+j))
            avgvel2 = avgvel2 + abs(r(3,(count-1)*npt+j))
         end do
         write (6,*) 'comparing ',i,' to ',count
         write (uwdat,*) (count-1)*step*1.0,avg/(3*npt)
         count = count+1

      end do
      write (6,*) 'avg velocity dyn',avgvel/(3*npt*strd)
      write (6,*) count,stre
      write (6,*) 'avgvel me ',avgvel2/(3*npt*stre)
      write (6,*) avg,npt,count

      end if
      
      stop

      end

