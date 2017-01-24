        subroutine init_press()
c

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/NBLIST.BLOCK'
        include 'COMMON/DYNA.BLOCK'


        character*8 name
        integer namel

        logical failure,flag

        integer j,k,l,iat1,iat2,istep,ip
        integer kl,km,nptid,i,m,kk

        double precision rcutv,ss,ss2,ss6
        double precision apc,bpc
        double precision epsg12id(maxunq),epsg6id(maxunq)
        double precision denstypid(maxunq)
        double precision mypi
         integer iii


c zero array and calc long range correct to short range (vdw) vir
c (Allen Tild eqn. 2.137)
c
        mypi = 4.d0*datan(1.d0)

c conversion factor from internal units to bars is 69480.d0
c (from 1kcal/molA = 69.48 pN ..)

        pconv = 69480.d0

c
       eninth=8.0d0/9.0d0
       fthird=4.0d0/3.0d0


c prepare num density of atm types
c for long range corr to press

       do i=1,npt
        numtyp(i) = 0.d0
       enddo

       do i=1,npt
        kk = ptid(i)
        numtyp(kk) = numtyp(kk)+1
       enddo
       do i=1,npt
        kk = ptid(i)
        denstyp(kk) = dble(numtyp(kk))/volbx
       enddo

c count how many atm types
c and pass info from atom id
c to atom type id

      nptid=1
        epsg12id(nptid) = epsgm12(1)
        epsg6id(nptid) = epsgm6(1)
        denstypid(nptid) = denstyp(ptid(1))

       do i=2,npt
        do l = 1,i-1
         if(ptid(i).eq.ptid(l)) goto 111
        enddo
        nptid=nptid+1
        epsg12id(nptid) = epsgm12(i)
        epsg6id(nptid) = epsgm6(i)
        denstypid(nptid) = denstyp(ptid(i))
111     continue
       enddo
c
        plrc = 0.d0
        elrc = 0.d0
        rcutv = dsqrt(cutvdw2)
        ss = 1.0d0/rcutv
        ss2 = ss*ss
        ss6 = ss2*ss2*ss2


c long range corrections (pressure and energy)

        do  l=1,nptid
          do   m=1,l

             apc= epsg12id(l)*epsg12id(m)
             bpc= epsg6id(l)*epsg6id(m)
        
             apc = apc*ss6*ss2*ss
             bpc = bpc*ss2*ss
             plrc = plrc + mypi*denstypid(l)*denstypid(m)*
     &             (eninth*apc - fthird*bpc)

             elrc = elrc + mypi*denstypid(l)*denstypid(m)*volbx*
     &             (eninth*apc/4.d0 - fthird*bpc/4.d0)
          enddo
        enddo
          
        write(6,'(a31,f10.4)') "L. r. correction to Pressure = ", 
     &  plrc*pconv," (bars)"
        write(6,'(a29,f10.4)') "L. r. correction to Energy = ", elrc

c   convert plrc to a viral term

      virlrc = -plrc*(3.d0*volbx)
c
       return
       end
