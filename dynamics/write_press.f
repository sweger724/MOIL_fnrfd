        subroutine write_press()
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/DYNA.BLOCK'
        include 'COMMON/EWALD.BLOCK'


        double precision virele,vir
        double precision press,press2


c vir = short range vir (incl. long range corr) + constr vir
c       + coul vir (= - coul ene)
c 1-4 vdw interactions are virvdw5
c (see ener14)

c sum up virials

         vir=0.d0

         virele = - (e_dir+e_ew_receip+e_self+e_corr)

c W = - r*f (and virlrc = - plrc * 3V )
c
        if(pdebug) then
        write(6,*) "wat virvdw_blk=",-virvdw1,
     &  "wat virvdw_sym=",-virvdw2
        write(6,*) "wat virvdw=",-virvdw1-virvdw2
        write(6,*) "prot virvdw_blk=",-virvdw3,
     &  "prot virvdw_sym=",-virvdw4
        write(6,*) "prot virvdw=",-virvdw3-virvdw4
        write(6,*) "1-4 vdw virvdw=",-virvdw5
        write(6,*) "tot virvdw=",
     &  -virvdw1-virvdw2-virvdw3-virvdw4-virvdw5
        write(6,*) "virele=",virele
        write(6,*) "virele (+virele14)=",virele-e_el14
        write(6,*) "vircon=",-vircon
        write(6,*) "vircent=",-vircent   
        write(6,*) "enkin=",enkin
        endif

         vir = -virvdw1-virvdw2-virvdw3-virvdw4-virvdw5-vircon
     &  - vircent + (virele - e_el14) 

         vir = vir + virlrc

        write(6,'(a11,f16.4)') "Press = ",
     &  ((2.d0*enkin - vir)/(3.d0*volbx))*pconv


        return

        end
