       subroutine Read_CG_input()

        implicit none

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/VELOC.BLOCK'
        include 'COMMON/LINE.BLOCK'
        include 'COMMON/CONSTRAN.BLOCK'
        include 'COMMON/FREADY.BLOCK'

        integer uCGpar, i, of, urtet,lc
        integer p1,p2,p3,p4
        double precision r,rx,ry,rz, CGtorsion, angle, Cconst,rms
        logical find
        character*6 getchar
        character*5 crdtyp

        lc=5


        if (find('fix2')) then
c get coordinate file for tether and overlap calculations
c they are read into the velocities and immediately transfered
c to coor2  ... this was the general functionality
c
c here only filling the CGstr(i) datastructure  from more(i) is done
                urtet = of()
                crdtyp = getchar('ctyp','CHARM',lc)
C                write(6,*)"DDD:",crdtyp
                call getvel(urtet,crdtyp(1:4))
                close (urtet)
                do i=1,npt
                        coor2(1,i) = velo(1,i)
                        coor2(2,i) = velo(2,i)
                        coor2(3,i) = velo(3,i)
                        CGstr(i) = more(i)
                end do
             call rmsd_weight(npt,coor(1,1),coor2(1,1),rms,.false.,ptms)
            Fix_2nd_structure = .TRUE.

            ecnyes = .TRUE.
        Cconst = 20.d0

        do i = 1,ntors
          if ((CGstr(itor2(i)).ne.0.0) .and.
     &        (CGstr(itor3(i)).ne.0.0)) then
            ncnst = ncnst + 1
            icnst1(ncnst) = itor1(i)
            icnst2(ncnst) = itor2(i)
            icnst3(ncnst) = itor3(i)
            icnst4(ncnst) = itor4(i)
            kcns(ncnst)   = Cconst
            cnseq(ncnst)  = CGtorsion(itor1(i),itor2(i),
     &                                itor3(i),itor4(i))
          endif
        end do

        do i = 1,nangl
          ncnst = ncnst + 1
          icnst1(ncnst) = iangl1(i)
          icnst2(ncnst) = iangl2(i)
          icnst3(ncnst) = iangl3(i)
          icnst4(ncnst) = 0
          kcns(ncnst)   = Cconst
          cnseq(ncnst)  = angle(iangl1(i),iangl2(i),iangl3(i))
        enddo

        do i = 1,nb

          ncnst = ncnst + 1
          icnst1(ncnst) = ib1(i)
          icnst2(ncnst) = ib2(i)
          icnst3(ncnst) = 0
          icnst4(ncnst) = 0
          kcns(ncnst)   = Cconst
          rx = coor2(1,ib1(i)) - coor2(1,ib2(i))
          ry = coor2(2,ib1(i)) - coor2(2,ib2(i))
          rz = coor2(3,ib1(i)) - coor2(3,ib2(i))
          r = dsqrt(rx**2 + ry**2 + rz**2)
          cnseq(ncnst)  = r
        enddo
        
        end if

c       CG parameters file
         if (find('CGpr')) then
            uCGpar = of()
            call read_CG_parameters(uCGpar)

            eCGyes = .true.
            nocut  = .true.
            ebyes  = .false.
            ethyes = .false.
            etoyes = .false.
            eimyes = .false.
            evdyes = .false.
            eelyes = .false.
            e14el_yes  = .false.
            e14v_yes   = .false.
            eenmyes = .false.

         endif

      return

      end
