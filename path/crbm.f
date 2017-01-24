                Subroutine crbm(r,scalar,sigma,grdcmx,grdcmy,grdcmz,
     1           grdlx,grdly,grdlz,divms,npt,nselec,pointr,debug,
     2           udata)
c
c Factor out rigid body motion from the current coordinate system
c
        double precision r(*),scalar(6),sigma(6),grdcmx(*),grdcmy(*)
        double precision grdcmz(*),grdlx(*),grdly(*),grdlz(*)
        double precision divms(*)
        integer pointr(*)
        integer npt,nselec,udata
        logical debug
c
c r          - coordinates
c scalar     - work vector used to store scalar products
c sigma      - fixed value of constraints (i.e. ideally scalar should be
c               equal to sigma)
c grdcm[x-z] - gradient of center of mass constraints
c grdl[x-z]  - gradient of rigid infitesimal rotation
c divms      - 1/m vector
c pointr     - a pointr to atoms which required the constraints
c               for a given structure
c npt      - number of atoms
c nselec     - number of selected atoms
c debug      - if .true. print debugging info.
c
c
c Note the coordinate constraints are of the form
c
c constraint = sum a(i)*r(i) + constant
c
c the lagrange multiplier (which we call here scalar) is
c scalar = sum a(i)*r(i) + constant
c
c LOCAL
        integer i,j,jj,jjj,k,npt3
        npt3 = 3 * npt
c
c calculate the scalar product of the constraints gradient and
c the "free" positions.
c
                do 1 i=1,6
                 scalar(i)=0.d0
1               continue

                do 2 j=1,nselec
                 divms(j)=1.0d0/divms(j)
                 i=pointr(j)
                 k = 3*(i-1)
                 jj=j+nselec
                 jjj=jj+nselec
C                write(6,*) 'grdlx(j),grdlx(jj),grdlx(jjj)',
C     1         grdlx(j),grdlx(jj),grdlx(jjj)
C                 write(6,*) 'grdly(j),grdly(jj),grdly(jjj)',
C     1          grdly(j),grdly(jj),grdly(jjj)
C                 write(6,*) 'grdlz(j),grdlz(jj),grdlz(jjj)',
C     1          grdlz(j),grdlz(jj),grdlz(jjj)
                 scalar(1)=scalar(1)+grdcmx(j)*r(k+1)
                 scalar(2)=scalar(2)+grdcmy(jj)*r(k+2)
                 scalar(3)=scalar(3)+grdcmz(jjj)*r(k+3)
                 scalar(4)=scalar(4)+grdlx(j)*r(k+1)+grdlx(jj)*r(k+2)
     1                   +grdlx(jjj)*r(k+3)
                 scalar(5)=scalar(5)+grdly(j)*r(k+1)+grdly(jj)*r(k+2)
     1                   +grdly(jjj)*r(k+3)
                 scalar(6)=scalar(6)+grdlz(j)*r(k+1)+grdlz(jj)*r(k+2)
     1                   +grdlz(jjj)*r(k+3)
2               continue
c
c add the constants to the scalar products
c
c NOTE: the constraints over the path direction and the rigid body
c motions underwent orthonormalization with respect to mass
c weighting. 
c
                do 3 i=1,6
C                       write(6,*) 'Scalar: ', scalar(i)
                        scalar(i)=scalar(i)-sigma(i)
3               continue
c
c calculate the "corrected" coordinates (finally)
c
c$DOIT IVDEP
                do 4 j=1,nselec
                 i=pointr(j)
                 jj=j+nselec
                 jjj=jj+nselec

                 k = 3*(i-1)
C         write(6,*) k,r(k+1),r(k+2),r(k+3)
                 r(k+1)=r(k+1)-divms(i)*(scalar(1)*grdcmx(j)
     1                  +scalar(4)*grdlx(j)
     2                  +scalar(5)*grdly(j)+scalar(6)*grdlz(j))

                 r(k+2)=r(k+2)-divms(i)*(
     1                  scalar(2)*grdcmy(jj)+scalar(4)*grdlx(jj)+
     2                  scalar(5)*grdly(jj)+scalar(6)*grdlz(jj))

                 r(k+3)=r(k+3)-divms(i)*(
     1                  scalar(3)*grdcmz(jjj)+scalar(4)*grdlx(jjj)+
     2                  scalar(5)*grdly(jjj)+scalar(6)*grdlz(jjj))
C         write(6,*) k,r(k+1),r(k+2),r(k+3)

4               continue
c
c test that the new coordinates satisfies the constraints
c
        if (.false.) then
                write(*,*)' npt nselec ',npt,nselec
                write(*,*)'divms = ',(divms(i),i=1,npt)

          do 5 i=1,6
                scalar(i)=0.d0
5         continue
          do 6 j=1,nselec
                i=pointr(j)
                jj=j+nselec
                jjj=jj+nselec

                k = 3*(i-1)

                scalar(1)=scalar(1)+r(k+1)*grdcmx(j)
                scalar(2)=scalar(2)+r(k+2)*grdcmy(jj)
                scalar(3)=scalar(3)+r(k+3)*grdcmz(jjj)
                scalar(4)=scalar(4)+grdlx(j)*r(k+1)+grdlx(jj)*r(k+2)
     1                   +grdlx(jjj)*r(k+3)
                scalar(5)=scalar(5)+grdly(j)*r(k+1)+grdly(jj)*r(k+2)
     1                   +grdly(jjj)*r(k+3)
                scalar(6)=scalar(6)+grdlz(j)*r(k+1)+grdlz(jj)*r(k+2)
     1                   +grdlz(jjj)*r(k+3)

6               continue
          do 7 i=1,6
                scalar(i)=scalar(i)-sigma(i)
7         continue
          write(udata,8)(scalar(i),i=1,6)
8         format(//,1x,' errors in coordinates constraints ',/,1x,
     1          4(e14.5),/,3(e14.5),//)
        end if
100     continue
         do i=1,npt
         divms(i)=1.0d0/divms(i)
         enddo
c        stop
        return
        end
