           program addion
           implicit none
c
c Starting from a system in a water box already add counter ions
c to make the system neutral. Changes that are made:
c 1 new coordinate file in which the last oxygens of the water molecules are
c        replaaced by the ions (the hydrogens are removed)
c 2 a new poly file (with the ion is produced
c
           include 'COMMON/LENGTH.BLOCK'
           include 'COMMON/COORD.BLOCK'
           include 'COMMON/UNITS.BLOCK'
           include 'COMMON/CONNECT.BLOCK'
           include 'COMMON/DEBUG.BLOCK'
           include 'COMMON/LINE.BLOCK'
           include 'COMMON/FREEZ.BLOCK'
           include 'COMMON/CONVERT.BLOCK'
           include 'COMMON/CCRD.BLOCK'
           include 'COMMON/MUTA.BLOCK'
c ipick - pick subset of particles, a vector of length maxpt
c         value of zero if particl not selected one if it is
c ipikpo - pointer to picked particles ipikpo(i) is the particle
c          number of the i-th selected particle
c npick  - number of picked particles
c of     - integer function for Opening a File, returned value is the
c          assigned unit number
c geti   - integer function to get integer value from a command line
c getd   - double precision function to get DP value from a line
c nstru  - number of structures in dynamics file
c namel  - length of program name
c inofrz - number of moving particles. used in reading the dynamics file
c level  - level of error found, level=0 only warning will be issued,
c          level=1 program stop
c urcrd,ucon - units of dynamics coord and connectivity files
c rcut   - cutoff distance to define a collision (rcut>distance)
c          input rcut, internally modified and used as rcut^2
c dd     - temporary variable storing the distance square between 2 pt
c xdiff ydiff zdiff - difference in x,y,z between pt (used in distance calc.)
c name - name of program (character) = contact
c find - find a charcter in line (logical function)
c fopen - check if file is open (logical function
c pickpt - true if pick instruction was found in input
c
           integer ipick(maxpt), ipikpo(maxpt)
           integer npick 
           integer of,geti,nstru
           integer namel,i,j,k,l,m,level
           integer urcrd,ucon,uwcrd,upoly
           double precision getd

           integer ion_name_length,ion_number,ion_monomer_length
           integer ncrd
           double precision step,cmx,cmy,cmz,cmm,tmpx,tmpy,tmpz,rg
           double precision time
           character*6 name
           character*5 getchar
           character*4 ion_name,ion_monomer
           character*80 tmp_char
           logical find,fopen
           logical pickpt
           logical lwxy
           data ucon,urcrd,uwcrd,upoly/4*99/

           integer placed,p,water(maxmono),Swaters, irand
           real rr

           irand = -1
           lpstr = 1

           stdi=5
           stdo=6
           totmon=0
           npt=0
           muta = .false.
           name='addion'
           namel=6
c  open junk file for rline
c
            jnkf=25
            open(unit=jnkf,status='scratch')
c default parameters
                ion_name_length = 4
                ion_monomer_length = 4
                ion_number = 0
                ion_name = 'NONE'
                ion_monomer = 'NONE'

1           continue
            call rline(name,namel,stdi)
            if (find('file')) then
               if (find ('conn')) then  
                if (find('muta')) muta = .true.
                ucon=of()
* get connectivity 
                call rconn(ucon)
* get coordinate file
               else if (find ('rcrd')) then
                 if (npt.eq.0) then
                  level = 1
                  call alert(name,namel,'Must read con file first',
     1                  24,level)
                  end if
                urcrd=of()
                call getcrd(urcrd,'CHARM')
                else if (find('wcrd')) then
                        uwcrd = of()
                else if (find('poly')) then
                        upoly = of()
                end if
               else if (find('acti')) then
                        go to 33
               else
               ion_name=getchar('iona',ion_name,ion_name_length)
              ion_monomer=getchar('ionm',ion_monomer,ion_monomer_length)
               ion_number=geti('#ion',ion_number)
               irand = geti('rand',irand)
                end if
             go to 1
33            continue

                call RLUXGO(223,irand,0,0) 

                if (ion_name.eq.'NONE' .or. ion_number.eq.0) then
                        write(*,*)' Nothing to do ! '
                        write(8,*)' ion name ',ion_name
                        write(*,*)' ion_number ',ion_number
                        call alert(name,namel,'no ions?',8,1)
                        stop
                end if
                
                Swaters = 0
                do i = 1,totmon
          if (moname(i).eq.'TIP3' .or. moname(i)(1:3).eq.'SPC') then
                    Swaters = Swaters + 1
                    water(Swaters) = i
                  end if
                end do
                if (Swaters .lt. ion_number) then
                  write(*,*)' not enough solvent to put ion'
                  write(*,*)' STOP '
                  stop
                end if
                
                placed = 0
                do while (placed .lt. ion_number)
                  write(6,*)"Drawing p..."
                  call RANLUX(rr,1)
                  p = int(rr*(Swaters-1))+1
                  write(6,*)"p=",p
                  if (moname(water(p)).eq. 'TIP3' .or. 
     &                moname(water(p))(1:3) .eq. 'SPC') then
                    moname(water(p)) = ion_monomer
                    k = poipt(water(p)-1)+1
                    ptnm(k) = ion_name
                  
                    do i = k+1,npt-2
                      poimon(i) = poimon(i+2)
                      ptnm(i) = ptnm(i+2)
                      ptid(i) = ptid(i+2)
                      do l =1,3
                        coor(l,i) = coor(l,i+2)
                      end do
                    end do
  
                    do i = water(p),totmon
                      poipt(i) = poipt(i)-2
                    end do
                  
                    npt = npt - 2
                    placed = placed + 1
                    write(6,*)"water number", water(p)," was replaced."
                  end if
                end do
                
                 call putcrd(uwcrd,'CHARM')
                 close(uwcrd)
             write(upoly,'(a/a6,a4,a9,i5.5/a)')
     $     "~","MOLC=(",BULK(1),")   #mon=",totmon,"~"
            write(upoly,'(10a5)') (moname(i), i=1,totmon)
            write(upoly,'(a4)') "*EOD"
            close(upoly)

           stop
           end
