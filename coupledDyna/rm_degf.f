      subroutine rm_degf(tpo,tempi,tempf,tgroup,ntemp,
     >     nofreez,inofrz,numremdgf)
c     
c     Pick temperatures used in different velocity scaling in the
c     system. The defualt is that all particles belong to temperature 1.
c     Maximum number of temperature is itempg, a variable that is set up
c     in dynamics. Currently is 10. It is doubtfull that it will be modified
c     tpo - pointer to current temperature tpo(i) is the pointer to tempi & tempf
c     to get the assigned temperature of the i-th atom
c     tempi - initial assigned temperature
c     tempf - final assigned temperature
c     ntemp - actual number of temperatures found
c     
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/SHAKE.BLOCK'
      include 'COMMON/MSHAKE.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/SSBP.BLOCK'
c     
      integer tpo(*),nofreez(*),tgroup(*),ntemp,inofrz,numremdgf
      double precision tempi(ntemp),tempf(ntemp)
c     local
      integer i,j,namel,level
      logical next
      double precision number
      character*8 name


      name = 'picktemp'
      namel = 8
c     
      do 1 i=1,ntemp
         tgroup(i) = 0
 1    continue

      if (nshak.gt.0) then
	 do 2 i=1,nshak
            j = tpo(ishak1(i))
            tgroup(j) = tgroup(j) - 1
            j = tpo(ishak2(i))
            tgroup(j) = tgroup(j) - 1
 2       continue
      end if
c     divide by two since subtracted a total of two for each constraint above
c     
      if (nshak.gt.0) then
	 do 3 i=1,ntemp
            tgroup(i) = tgroup(i)/2
 3       continue
      end if
c     subtract matrix shake degrees of freedom
      if (nshakm.gt.0) then
         do 12 i=1,nwaters
            j = tpo(dpoipt(idxtip3(i)))
            tgroup(j) = tgroup(j) - 3
 12      continue
      end if

      do 4 i=1,inofrz
	 j = nofreez(i)
	 tgroup(tpo(j)) = tgroup(tpo(j)) + 3
 4    continue

c     
c     we assume below that only diatomic has two rotations
c     we are going to introduce error in strictly linear molecules
c     which are not in our agenda anyway.
c     
      tgroup(1) = tgroup(1) - numremdgf

      do 5 i=1,ntemp
	 if (tgroup(i).le.0) then
            level = 1
            call alert(name,namel,' subset with no particles ',26,level)
	 end if
 5    continue

      return
      end

