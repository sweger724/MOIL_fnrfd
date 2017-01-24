           program xdifferential
c
c Xtract coordinates from the three first structures in a  trajectory
C to compute a path that is solution of the verlet algorithm in length
C If original path and new path are the same convergence in a refinement
C procedure has been achieved
c
           implicit none
           include 'COMMON/LENGTH.BLOCK'
           include 'COMMON/CONNECT.BLOCK'
           include 'COMMON/UNITS.BLOCK'
           include 'COMMON/LINE.BLOCK'
      include 'COMMON/NBLIST.BLOCK'
      include 'COMMON/SYMM.BLOCK'
           include 'COMMON/COORD.BLOCK'
      include 'COMMON/CONVERT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/VELOC.BLOCK'
        include 'COMMON/SPECL.BLOCK'
        include 'COMMON/CONSTRAN.BLOCK'
        include 'COMMON/EWALD.BLOCK'
        include 'COMMON/METAL.BLOCK'
           include 'COMMON/DEBUG.BLOCK'
           include 'COMMON/FREEZ.BLOCK'
           include 'COMMON/CCRD.BLOCK'
      include 'COMMON/PDQ.BLOCK'
        include 'COMMON/MASSFAC.BLOCK'
            include 'COMMON/SCNDRV.BLOCK'
      include 'COMMON/GETVEC.BLOCK'
C SGB MODIFICATIONS :
        include 'COMMON/SGB.BLOCK'
C END SGB MODIFICATIONS
c ipick - pick subset of particles, a vector of length maxpt
c         value of zero if particle not selected one if it is
c npick  - number of picked particles
c of     - integer function for Opening a File, returned value is the
c          assigned unit number
c geti   - integer function to get integer value from a command line
c nstru  - number of structures in dynamics file
c namel  - length of program name
c inofrz - number of moving particles. used in reading the dynamics file
c level  - level of error found, level=0 only warning will be issued,
c          level=1 program stop
c urcrd,ucon - units of dynamics coord and connectivity files
c
           integer ipick(maxpt)
           integer npick 
           integer of,geti,nstru
           integer namel,i,j,l,level,i1,i2,i3,i4
           integer urcrd,uwtor,ucon,ucrd
           double precision dx2,dy2,dz2,dx3,dy3,dz3
           double precision ux,uy,uz,vx,vy,vz,uv,uu,vv,phi
           double precision d0(lgrid),r(3,3*maxpt),tmpx,tmpy,tmpz
           double precision e0,dv(3,maxpt),drdl(3,maxpt),dave2
           double precision dr22,sprod,d2rdl2an(3,maxpt)
           double precision getd,tst,d02(lgrid),e(3)
           integer sgbboolint,nlist,k

           

           character*6 name
           character*4 ctype
           logical find,fopen
           logical pickpt,empty
           data ucon,urcrd,uwtor/3*99/

           ctype='DYNA'
           i1    = 0
           i2    = 0
           i3    = 0
           i4    = 0
           lpstr = 1
           norew = .false.

           stdi=5
           stdo=6
           totmon=0
           npt=0
        nb     = 0
        nmb    = 0
        nangl  = 0
        ntors  = 0
        nimp   = 0
        lestyp = 0
        lpstr  = 1
           name='xtors'
           namel=5
           debug  = .false.
c  open junk file for rline
c
            jnkf=25
            open(unit=jnkf,status='scratch')
c default parameters
            nstru= 1
            pickpt=.false.
            sgbboolint= 0
            sgba = 0
        gbsu = 0
        gbsabool = .false.

C SGB Modifications :
        chainno = 1
        alphac(chainno) = 1
C END SGB MODIFICATIONS
        lcent = .false.
        ctrue = .true.
        shift = .false.
        hydro_scale = 1.d0
             call init_ef()
             nocut  = .true.


            call init_var()

1           continue
            call rline(name,namel,stdi)
            if (find('coor')) then
               call get4c(ctype,empty)
               write (6,*) 'file type ',ctype
            end if
            if (find('norw')) norew=.true.
            if (find('file')) then
               if (find ('conn')) then  
                ucon=of()
c get connectivity 
                call rconn(ucon)
c get coordinate file
               else if (find ('rcrd')) then
                if (npt.eq.0) then
                 level = 1
                 call alert(name,namel,'Must read con file first',
     1                  24,level)
                end if
                urcrd=of()
               else if (find('wtor')) then
                uwtor=of()
               end if
              else 
               nstru=geti('#str',nstru)
               ENERGYC=getd('ENER',ENERGYC)
        nlist   = geti('list',nlist)
        if (find('hvdw')) hvdw0 = .false.
        rmax     = getd('rmax',rmax)
        eps     = (getd('epsi',(eps)))
        if (find('cdie')) ctrue = .true.
        if (find ('sgbb')) then
           sgbbool=.true.
        end if
        if (find ('gbsa')) then
           gbsabool=.true.
        end if
        gbsu=geti('gbsu',gbsu)
               if (find('pick')) then
                call pick(ipick,npick)
                do 11 i=1,npt
                 if (i1.eq.0 .and. ipick(i).ne.0) then
                  i1=i
                 else if (i2.eq.0 .and. ipick(i).ne.0) then
                  i2=i
                 else if (i3.eq.0 .and. ipick(i).ne.0)then
                  i3=i
                 else if (i4.eq.0 .and. ipick(i).ne.0)then
                  i4=i
                 end if
11              continue
               end if
               if (find ('action')) goto  3
             end if
             go to 1
3            continue
        if (rmax.gt.0) then
                cutvdw2 = rmax
                cutele2 = rmax
        end if

        if (cutvbig2.lt.0.d0) cutvbig2 = cutvdw2 + 2.d0
        if (cutebig2.lt.0.d0) cutebig2 = cutele2 + 2.d0

        cutvdw2  = cutvdw2*cutvdw2
        cutvbig2 = cutvbig2*cutvbig2
        cutele2  = cutele2*cutele2
        cutebig2 = cutebig2*cutebig2

        if (cutmono2.lt.0) then
                cutmono2 = cutebig2*1.44
        else
                cutmono2 = cutmono2*cutmono2
        end if

c initialze no freez vecotr
        inofrz = npt
        do 21 i=1,npt
                zerofrz(i) = 1
21      continue



          if (.not. fopen(ucon)) then
             level=1
             call alert(name,namel,'ucon not opened',15,level)
           else if (.not. fopen(urcrd)) then
             level=1
             call alert(name,namel,'urcrd not opened',16,level)
           end if
             
C Compute average distance:
           do l=1,nstru-1
           d0(l)=0.0d0
           enddo
           rewind urcrd
           l=1
            call rpath_seq(urcrd,1)
            do i=1,npt
            r(1,i) = coor(1,i)*dsqrt(ptms(i))
                r(2,i) = coor(2,i)*dsqrt(ptms(i))
                r(3,i) = coor(3,i)*dsqrt(ptms(i))
            enddo
      if (gbsabool) then
         call make_rborn
       endif
c         call nbondm()
c           call eforce()
c           call wener(6)
            call rpath_seq(urcrd,2)
            do i=1,npt
            r(1,i+npt) = coor(1,i)*dsqrt(ptms(i))
            r(2,i+npt) = coor(2,i)*dsqrt(ptms(i))
            r(3,i+npt) = coor(3,i)*dsqrt(ptms(i))
            enddo
            do i=1,npt
            tmpx = r(1,i+npt) - r(1,i)
            tmpy = r(2,i+npt) - r(2,i)
            tmpz = r(3,i+npt) - r(3,i)
            d0(1)=d0(1)+tmpx*tmpx + tmpy*tmpy + tmpz*tmpz
            enddo
c           call eforce()
c           call wener(6)
            do l=3,nstru
            do i=1,npt
            r(1,i) = r(1,i+npt)
            r(2,i) = r(2,i+npt)
            r(3,i) = r(3,i+npt)
            enddo
            call rpath_seq(urcrd,l)
c           call eforce()
c           call wener(6)
            do i=1,npt
            r(1,i+npt) = coor(1,i)*dsqrt(ptms(i))
            r(2,i+npt) = coor(2,i)*dsqrt(ptms(i))
            r(3,i+npt) = coor(3,i)*dsqrt(ptms(i))
            enddo
            do i=1,npt
            tmpx = r(1,i+npt) - r(1,i)
            tmpy = r(2,i+npt) - r(2,i)
            tmpz = r(3,i+npt) - r(3,i)
            d0(l-1)=d0(l-1)+tmpx*tmpx + tmpy*tmpy + tmpz*tmpz
            enddo
            write(6,*) l-1,d0(l-1)
            enddo
            dave =0.0d0
            dave2 =0.0d0
            do l=2,2
c            do l=1,nstru-1
            d02(l) =d0(l)/10000.0
            d0(l)= dsqrt(d0(l))/100.0
c            write(6,*) l,d0(l)
            dave  = dave + d0(l)
            dave2  = dave2 + d02(l)
            enddo
c            dave = dave / (nstru-1)
c            dave2 = dave2 / (nstru-1)
            write(6,*) 'dave',dave,energyc,'dave2 ',dave2
   
            call init_wre(stdo)
c read dynamics structures
           rewind urcrd
           j = 0
c phi ala

           do 9 l=1,2

                 call rpath_seq(urcrd,l)

           if (l.eq.1) then
           do i=1,npt
                r(1,i) = coor(1,i)*dsqrt(ptms(i))
                r(2,i) = coor(2,i)*dsqrt(ptms(i))
                r(3,i) = coor(3,i)*dsqrt(ptms(i))
           enddo
           else if (l.eq.2) then
           do i=1,npt
                r(1,i+npt) = coor(1,i)*dsqrt(ptms(i))
                r(2,i+npt) = coor(2,i)*dsqrt(ptms(i))
                r(3,i+npt) = coor(3,i)*dsqrt(ptms(i))
            tmpx = r(1,i+npt) - r(1,i)
            tmpy = r(2,i+npt) - r(2,i)
            tmpz = r(3,i+npt) - r(3,i)           
          r(1,i+npt)=r(1,i)+tmpx/100.0
          r(2,i+npt)=r(2,i)+tmpy/100.0
          r(3,i+npt)=r(3,i)+tmpz/100.0
          coor(1,i)=r(1,i+npt)/dsqrt(ptms(i))
          coor(2,i)=r(2,i+npt)/dsqrt(ptms(i))
          coor(3,i)=r(3,i+npt)/dsqrt(ptms(i))
            enddo
         if (esymyes) call squeeze()
         if (esymyes) call syminit()
         call nbondm()
           call eforce()
c           write(6,*) evdyes, eelyes,rmax
           call wener(6)
           e0 = e_total
            write(6,*) 'etotal',e0
            do  i = 1,npt
               dv(1,i) = dpot(1,i)/dsqrt(ptms(i))
               dv(2,i) = dpot(2,i)/dsqrt(ptms(i))
               dv(3,i) = dpot(3,i)/dsqrt(ptms(i))
c       write(6,*) dv(1,i),dv(2,i),dv(3,i)
            end do
             endif
9          continue
        do l=1,10000
        dr22=0.0d0
C Compute intermediate structure:
        do i=1,npt
      if (l.ne.1) then
      r(1,i) = r(1,i+npt)
      r(2,i) = r(2,i+npt)
      r(3,i) = r(3,i+npt)
      r(1,i+npt) = r(1,i+2*npt)
      r(2,i+npt) = r(2,i+2*npt)
      r(3,i+npt) = r(3,i+2*npt)
      endif
        r(1,i+2*npt) = -r(1,i) + 2.0*r(1,i+npt)
        r(2,i+2*npt) = -r(2,i) + 2.0*r(2,i+npt)
        r(3,i+2*npt) = -r(3,i) + 2.0*r(3,i+npt)
C Calculate tangential vector:
        drdl(1,i)=r(1,i+npt)-r(1,i)
        drdl(2,i)=r(2,i+npt)-r(2,i)
        drdl(3,i)=r(3,i+npt)-r(3,i)
        dr22=dr22+(drdl(1,i)**2+drdl(2,i)**2
     & +drdl(3,i)**2)
      enddo
      write(6,*) 'dr22', dr22
       dr22=dsqrt(dr22)
        tst=0.0d0
        do i=1,npt
      drdl(1,i)=drdl(1,i)/dr22
      drdl(2,i)=drdl(2,i)/dr22
      drdl(3,i)=drdl(3,i)/dr22
        tst=tst+(drdl(1,i)**2+drdl(2,i)**2
     & +drdl(3,i)**2)
        enddo
        write(6,*) 'tst',dsqrt(tst)
C Compute analytical derivative:
        do i=1,npt
      sprod=dv(1,i)*drdl(1,i)+dv(2,i)*drdl(2,i)+
     & dv(3,i)*drdl(3,i)
      d2rdl2an(1,i)=dv(1,i)-sprod*drdl(1,i)
      d2rdl2an(2,i)=dv(2,i)-sprod*drdl(2,i)
      d2rdl2an(3,i)=dv(3,i)-sprod*drdl(3,i)
      d2rdl2an(1,i)=d2rdl2an(1,i)/2.0/(ENERGYC-e0)
      d2rdl2an(2,i)=d2rdl2an(2,i)/2.0/(ENERGYC-e0)
      d2rdl2an(3,i)=d2rdl2an(3,i)/2.0/(ENERGYC-e0)
      r(1,i+2*npt)=r(1,i+2*npt)-d2rdl2an(1,i)*dave2
      r(2,i+2*npt)=r(2,i+2*npt)-d2rdl2an(2,i)*dave2
      r(3,i+2*npt)=r(3,i+2*npt)-d2rdl2an(3,i)*dave2
      coor(1,i) =  r(1,i+2*npt)/dsqrt(ptms(i))
      coor(2,i) =  r(2,i+2*npt)/dsqrt(ptms(i))
      coor(3,i) =  r(3,i+2*npt)/dsqrt(ptms(i))
C double checking:
      e(1)=(r(1,i+2*npt)+r(1,i)-2.0*r(1,i+npt))/dave2
     & +d2rdl2an(1,i)
      e(2)=(r(2,i+2*npt)+r(2,i)-2.0*r(2,i+npt))/dave2
     & +d2rdl2an(2,i)
      e(3)=(r(3,i+2*npt)+r(3,i)-2.0*r(3,i+npt))/dave2
     & +d2rdl2an(3,i)
c       write(6,*) e(1),e(2),e(3)
      enddo
       call eforce()
           call wener(6)
      e0=e_total
            do  i = 1,npt
               dv(1,i) = dpot(1,i)/dsqrt(ptms(i))
               dv(2,i) = dpot(2,i)/dsqrt(ptms(i))
               dv(3,i) = dpot(3,i)/dsqrt(ptms(i))
            end do

        k=npt
          if (mod(l,100).eq.0) then
              write(uwtor) e0,
     $             (r(1,j)/dsqrt(ptms(j-k)),j=k+1,k+npt),
     $             (r(2,j)/dsqrt(ptms(j-k)),j=k+1,k+npt),
     $             (r(3,j)/dsqrt(ptms(j-k)),j=k+1,k+npt)
          endif
         enddo

      

         K=2*npt
              write(uwtor) e0,
     $             (r(1,j)/dsqrt(ptms(j-k)),j=k+1,k+npt),
     $             (r(2,j)/dsqrt(ptms(j-k)),j=k+1,k+npt),           
     $             (r(3,j)/dsqrt(ptms(j-k)),j=k+1,k+npt)        

         K=3*npt
              write(uwtor) e0,
     $             (r(1,j)/dsqrt(ptms(j-k)),j=k+1,k+npt),
     $             (r(2,j)/dsqrt(ptms(j-k)),j=k+1,k+npt), 
     $             (r(3,j)/dsqrt(ptms(j-k)),j=k+1,k+npt)
           stop
           end
      
