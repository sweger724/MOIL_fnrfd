        subroutine comc(r0,dmass,npt,pointr,
     1  nselec,grdcmx,grdcmy,grdcmz,grdlx,grdly,grdlz,sigma,debug)
        implicit none

        include 'COMMON/LENGTH.BLOCK'
c
c calculate gradient of the Center Of Mass Constraints.
c and then orthonormalize the constraints vectors.
c
c *** input
c r0       - the coordinates of the reference system for which
c            cenetr of mass translation and rotation are compared
c dmass    - double precision vector of masses
c debug    - logical. true=print a lot of debugging information
c work vectors - x1,y1,z1,dmass1,dmass2
c
c *** output
c constraints value at initial configuration, should remain
c               constants.
c grdcmx - gradient of the translation of center of mass (X)
c grdcmy                   and Y
c grdcmz                   and Z
c
c grdlx  - gradient of rotation (X)
c grdly                         (Y)
c grdlz                         (Z)
c
c NOTE THAT THE VECTORS ABOVE UNDERWENT ORTHONORMALIZTION WITH
c RESPCT TO EACH OTHER.
c
c A change as of 8/22 is the addition of selection array. pointr
c is a pointer to the selecetd atomsn and nselec is the number of
c selected atoms
c
        include 'COMMON/UNITS.BLOCK'
        double precision r0(3,*)
        double precision grdcmx(*),grdcmy(*),grdcmz(*)
        double precision grdlx(*),grdly(*),grdlz(*)
        double precision sigma(*),dmass(*)
        integer npt,nselec,jj,jjj
        integer pointr(*)
        logical debug
c
c       local
c
        double precision norm,test(6)
        integer i,j,k,nsel3
        double precision dmass1(maxpt),dmass2(maxpt)
        double precision x1(maxpt),y1(maxpt),z1(maxpt)

        nsel3=3*nselec
c
c The zero order values of the constraints (which are not zero)
c are stored in sigma(1-6)
c sigma(1-3) - center of mass positions for x y z
c sigma(4-6) - Orientation value for x y & z ( zero before
c           orthonormalization)
c
        do 11 i=1,6
                sigma(i)=0.d0
11      continue
c
c calculate double precision mass and picking the
c selected atoms to a separate vector
c
        do 1 j=1,nselec
                i=pointr(j)
                dmass1(j)=dmass(i)
                dmass2(j)=1.d0/dmass(i)
                x1(j)=r0(1,i)
                y1(j)=r0(2,i)
                z1(j)=r0(3,i)
C               write(*,*), 'r0: ', k, r0(1,i), r0(2,i), r0(3,i)
1       continue
c
c calculate the gradient of the translation vector.
c Note the strange normalization: If a1(i) and a2(i) are
c vectors of constraint gradients and m(i) is the mass,
c the scalar product is defined by
c SUM a1(i)*a2(i)/m(i) 
c
        norm=0.d0
C$DOIT IVDEP
        do 2 j=1,nselec

                jj=j+nselec
                jjj=jj+nselec

                grdcmx(j)=dmass1(j)
                grdcmx(jj)=0.d0
                grdcmx(jjj)=0.d0

                grdcmy(j)=0.d0
                grdcmy(jj)=dmass1(j)
                grdcmy(jjj)=0.d0

                grdcmz(j)=0.d0
                grdcmz(jj)=0.d0
                grdcmz(jjj)=dmass1(j)

                grdlx(j)=0.d0
                grdlx(jj)=dmass1(j)*z1(j)
                grdlx(jjj)=-dmass1(j)*y1(j)

                grdly(j)=-grdlx(jj)
                grdly(jj)=0.d0
                grdly(jjj)=dmass1(j)*x1(j)

                grdlz(j)=-grdlx(jjj)
                grdlz(jj)=-grdly(jjj)
                grdlz(jjj)=0.d0
2       continue
        do 25 j=1,nselec
c
c actually it is: norm=norm+grdcm[xyz](i)*grdcm[xyz](i)/dmass(i)
c
                norm=norm+grdcmx(j)
c
c opportunity for constraints values calculation
c
                sigma(1)=sigma(1)+dmass1(j)*x1(j)
                sigma(2)=sigma(2)+dmass1(j)*y1(j)
                sigma(3)=sigma(3)+dmass1(j)*z1(j)
25      continue
        norm=1.d0/dsqrt(norm)
c
c normalize grdcmx grdcmy grdcmz and the constraint functions
c
        sigma(1)=sigma(1)*norm
        sigma(2)=sigma(2)*norm
        sigma(3)=sigma(3)*norm
        do 3 i=1,nsel3
                grdcmx(i)=grdcmx(i)*norm
                grdcmy(i)=grdcmy(i)*norm
                grdcmz(i)=grdcmz(i)*norm
3       continue
        if (debug) then
                do 31 i=1,6
                        test(i)=0.d0
31              continue
                do 32 i=1,nsel3
                        j=i-((i-1)/nselec)*nselec
                        test(1)=test(1)+grdcmx(i)*grdcmx(i)*dmass2(j)
                        test(2)=test(2)+grdcmy(i)*grdcmy(i)*dmass2(j)
                        test(3)=test(3)+grdcmz(i)*grdcmz(i)*dmass2(j)
                        test(4)=test(4)+grdcmx(i)*grdlx(i)*dmass2(j)
                        test(5)=test(5)+grdcmx(i)*grdly(i)*dmass2(j)
                        test(6)=test(6)+grdcmx(i)*grdlz(i)*dmass2(j)
32              continue
                write(stdo,*)' before orthog. norms of grdcm[xyz]'
                write(stdo,*)' scalar prod. of grdcmx & grdl[x-z] '
                write(stdo,*)(test(i),i=1,6)
        end if
c
c orthogonalize the vectors with mass weighting. Obviously after
c orthonormaliztion, the meaning of rotation in a given direction
c may be lost.
c
c ** grdlx
        call orthg(grdcmy,grdlx,nsel3,dmass2,sigma,2,4)
        call orthg(grdcmz,grdlx,nsel3,dmass2,sigma,3,4)
c ** grdly
        call orthg(grdcmx,grdly,nsel3,dmass2,sigma,1,5)
        call orthg(grdcmz,grdly,nsel3,dmass2,sigma,3,5)
        call orthg(grdlx ,grdly,nsel3,dmass2,sigma,4,5)
c ** grdlz
        call orthg(grdcmx,grdlz,nsel3,dmass2,sigma,1,6)
        call orthg(grdcmy,grdlz,nsel3,dmass2,sigma,2,6)
        call orthg(grdlx ,grdlz,nsel3,dmass2,sigma,4,6)
        call orthg(grdly ,grdlz,nsel3,dmass2,sigma,5,6)

c
c check that everything is orthonormal
c
        if (debug) then
        do 61 i=1,6
                test(i)=0.d0
61      continue
        do 7 i=1,nsel3
                j=i-((i-1)/nselec)*nselec
                test(1)=test(1)+grdcmx(i)*grdcmx(i)*dmass2(j)
                test(2)=test(2)+grdcmy(i)*grdcmy(i)*dmass2(j)
                test(3)=test(3)+grdcmz(i)*grdcmz(i)*dmass2(j)
                test(4)=test(4)+grdlx(i)*grdlx(i)*dmass2(j)
                test(5)=test(5)+grdly(i)*grdly(i)*dmass2(j)
                test(6)=test(6)+grdlz(i)*grdlz(i)*dmass2(j)
7       continue
        write(*,*)' normalization of pgrd grdcm[x,y,z] grdl[x,y,z] '
        write(*,*)(test(i),i=1,6)
        do 8 i=1,6
                test(i)=0.d0
8       continue
        do 9 i=1,nsel3
                j=i-((i-1)/nselec)*nselec
                test(1)=test(1)+grdlx(i)*grdcmx(i)*dmass2(j)
                test(2)=test(2)+grdlx(i)*grdcmy(i)*dmass2(j)
                test(3)=test(3)+grdlx(i)*grdcmz(i)*dmass2(j)
9       continue
        write(*,*)'scalar products (should be zero) of:'
        write(*,*)'grdlx-grdcmx grdlx-grdcmy grdlx-grdcmz '
        write(*,*)(test(i),i=1,3)
        end if
c --- end of debug statement
        return
        end
