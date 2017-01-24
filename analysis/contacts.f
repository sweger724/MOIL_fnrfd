           program contact
c
c calculate distance and collision numbers for a picked group
c For the selected group (e.g. a diatomic ligand), the program
c calculate the distance between all atoms and the selected group
c and send to output list of atoms that were in a distance smaller
c than rcut, this is done for a sequence of dynamics structures.
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
c ipick - pick subset of particles, a vector of length maxpt
c         value of zero if particl not selected one if it is
c ipikpo - pointer to picked particles ipikpo(i) is the particle
c          number of the i-th selected particle
c nopick - similar to ipikpo in reveres. nopick(i) is the particle
c          number of the i-th NOT-selected  particle
c coll   - collision counter, coll(i) the number of collisions of
c          the nopick particles with any of the picked particles
c          this is an average over the complete trajectory
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
           integer ipick(maxpt), ipikpo(maxpt), nopick(maxpt)
           integer col(maxpt)
           integer tot_col
           integer npick 
           integer of,geti,nstru
           integer namel,i,j,k,l,level
           integer urcrd,ucon,usum,uave
           integer rbin
           double precision rcut,dd,xdiff,ydiff,zdiff
           double precision getd
           character*7 name
           logical find,fopen
           logical pickpt
           data ucon,urcrd,usum,uave/4*99/


           lpstr = 1
           norew = .false.

           stdi=5
           stdo=6
           rbin = 1 

           totmon=0
           npt=0
           name='contact'
           namel=7
*  open junk file for rline
*
            jnkf=25
            open(unit=jnkf,status='scratch')
* default parameters
            nstru=1
            pickpt=.false.
            rcut = 5.d0

1           continue
            call rline(name,namel,stdi)
            if (find('norw')) norew=.true.
            if (find('file')) then
               if (find ('conn')) then  
                ucon=of()
* get connectivity 
                call rconn(ucon)
               end if
               if (find('wsum')) then
c write summary file
                usum = of()
               end if
                if (find('wave')) then
c write average collision number
                 uave = of()
                end if
* get coordinate file
               if (find ('rcrd')) then
                 if (npt.eq.0) then
                  level = 1
                  call alert(name,namel,'Must read con file first',
     1                  24,level)
                 end if
                 urcrd=of()
               end if
              else 
               rcut = getd('rcut',rcut)
               nstru=geti('#str',nstru)
               if (find('pick')) then
                call pick(ipick,npick)
                pickpt=.true.
               end if
               if (find ('action')) goto  3
             end if
             go to 1
3            continue

             rcut = rcut*rcut
c initialized the nofreez vector and the col (collision counters) vector
c
             inofrz=npt
             do 4 i=1,npt
                nofreez(i)=i
                col(i) = 0
4            continue

          if (.not. fopen(ucon)) then
             level=1
             call alert(name,namel,'ucon not opened',15,level)
            else if (.not. fopen(urcrd)) then
             level=1
             call alert(name,namel,'urcrd not opened',16,level)
            else if (.not. fopen(usum)) then
             level=1
             call alert(name,namel,'usum  not opened',16,level)
            else if (.not. fopen(uave)) then
             level =1
             call alert(name,namel,'uave  not opened',16,level)
            end if
             
           if (.not.pickpt) then
            level = 1
            call alert(name,namel,'Missing pick cmnd',17,level)
           end if

           npick = 0
           k = 0
           do 6 i=1,npt
            if (ipick(i) .eq. 1) then
             npick = npick + 1
             ipikpo(npick) = i
            else
             k = k + 1
             nopick(k) = i
            end if
6          continue

   
* read dynamics structures
           rewind urcrd
           do 8 l=1,nstru
                if (.not.norew) rewind urcrd
                call rdyncrd(urcrd,l,inofrz,nofreez,rbin)
                tot_col = 0
c
c calculate distances between the npick  and all other partciles
c
                write(stdo,100)l
100             format(1x,'Dynamics structure # ',i5,1x,' Collisions:')
                do 7 j=1,npt-npick
                 k = nopick(j)
                 do 7 i=1,npick
                  xdiff = coor(1,ipikpo(i)) - coor(1,k)
                  ydiff = coor(2,ipikpo(i)) - coor(2,k)
                  zdiff = coor(3,ipikpo(i)) - coor(3,k)
                  dd = xdiff*xdiff + ydiff*ydiff + zdiff*zdiff
                  if (dd.lt.rcut) then
                   col(k) = col(k) + 1
                   tot_col = tot_col + 1
                   write(stdo,101)ptnm(k),moname(poimon(k)),poimon(k)
101                format(1x,a4,1x,a4,1x,i7)
                  end if
7               continue
           write(uave,*)l,tot_col
8          continue
           write(usum,102)
102        format(1x,//,'Summary of collisions')
           do 9 j=1,npt-npick
            k = nopick(j)
            if (col(k).gt.0) then
            write(usum,103)ptnm(k),moname(poimon(k)),poimon(k),
     1       float(col(k))/nstru
c @TFB: the above gives error with gfortran; replace with ->
c       write(usum,103)ptnm(k),moname(poimon(k)),poimon(k),(col(k))/nstru
103          format(1x,' pt: ',a4,1x,a4,1x,i7,'  Collisions: ',f8.3)
            end if
9          continue
           stop
           end
