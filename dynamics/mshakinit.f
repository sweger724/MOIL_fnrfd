        subroutine mshakinit(debug) 
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/MSHAKE.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/UNITS.BLOCK'

        integer i,j
        integer namel,level
        character*9 name
        double precision dtnrm,detm

        logical debug

c
c
c consa() if first particle in constraint, consb() is second
c pointm() is last constraint coupled in current set,conseq() is equil. pos.
c
c creates list of particle numbers vlist() and another vector pointer vpart()
c giving the last position in the total list for a particular constraint set
c also creates the vector pointer (vcon) showing the last constraint 
c for this set
c
c first make the lists of constraints. this version currently
c finds all tip3 waters and shakes only them
c
        name = 'mshakinit'
        namel = 9
        write (stdo,*) 'initializing matrix shake on TIP3 water'
        write (stdo,*)
        write (stdo,*) ' General message: '
        write (stdo,*) '*** Freezing TIP3-s and MSHAKING them '
        write (stdo,*) '*** is not a good idea. '
        write (stdo,*)

        if (prll_on_off) then
         tip3_mono = nwaters
         if (tip3_mono.eq.0) then
                level = 0
                call alert(name,namel,' NO TIP3! ',10,level)
                return
         end if
         !write(6,*)"A:",my_pe,num_pes
         call load_balance(tip3_mono,my_pe,num_pes,tip3_st,tip3_en)
         write(stdo,*) ' No. of tip3 in p-e ',my_pe,' is ',tip3_mono
         write(stdo,*) ' Starting from ',tip3_st,' ending at ',tip3_en
        else
         tip3_st   = 1
         tip3_en   = nwaters
        end if
        nshakm = 3 * nwaters
c
c initiate the C matrix for TIP3/SPC/SCPE water
c It is assumed that if one of the water is TIP3
c then all water is TIP3 (and similarly SPC/SPCE)
c
        mreal = realmono(idxtip3(nwaters))
        if (moname(mreal).eq.'TIP3') then
          reqoh2 = 0.9572d0*0.9572d0
          reqhh2 = 1.5139d0*1.5139d0
        else if (moname(mreal)(1:3).eq.'SPC') then
          reqoh2 = 1.00
          reqhh2 = 1.63298086184d0**2
        end if
        
        invmh = 1.d0/ptms(dpoipt(mreal))
        invmo = 1.d0/ptms(dpoipt(mreal)-2)

        cmat(1,1)=reqoh2*(invmo+invmh)
        cmat(2,2)=cmat(1,1)
        cmat(3,3)=reqhh2*(invmh+invmh)
c
c calculate the cos directly from the distances so the numerical
c error will be minimized
c We use the cosine formula:
c rhh^2 = 2*roh^-2*roh^2*cos(hoh)
c cos(hoh)=(2roh^2-rhh^2)/(2roh^2)
c roh^2 = roh^2+rhh^2-2roh*rhh*cos(ohh)
c cos(ohh) = rhh/(2*roh)
c
c AND
c roh^2*cos(hoh)=2*roh^2-rhh^2
c roh*rhh*cos(ohh)=rhh^2/2
c
        cmat(1,2)   = invmo*(reqoh2-0.5d0*reqhh2)
        cmat(2,1)   = cmat(1,2)
        cmat(3,1)   = invmh*0.5d0*reqhh2
        cmat(3,2)   = cmat(3,1)
        cmat(1,3)   = cmat(3,1)
        cmat(2,3)   = cmat(3,1)
c
c invert the matrix as required by mshake
c
        call invmat(cmat,3,dtnrm,detm)

        return
        end

