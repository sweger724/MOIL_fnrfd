        subroutine compare_structures()
        implicit none
c
c  Compute rms of two structures 
c
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/COUPLED_DYNA.BLOCK'
        include 'mpif.h'

C       Local variables
        integer status(MPI_STATUS_SIZE)
        double precision  coor3(3,maxpt)
        integer dest
        double precision rms, tmp_mass(maxpt)

        integer j, k
C       integer CAindex(maxmono)
        double precision CA(3,maxmono), CA3(3,maxmono)
C
C  Comunicate the data
C
C     Process 0 part
C
        write(*,*) 'ProcID, maxpt:   ', procID, maxpt
        write(*,*) 'ProcID, # of particles:   ', procID, npt

        if (procID.eq.0) then
          call MPI_RECV(coor3,3*npt,MPI_DOUBLE_PRECISION,
     &      MPI_ANY_SOURCE,DYNA_TAG,MPI_COMM_WORLD,status,MPIerr)

C       Send to process 0
          dest = 1 
          call MPI_SEND(coor,3*npt,MPI_DOUBLE_PRECISION,
     &              dest,DYNA_TAG,MPI_COMM_WORLD,MPIerr)
        end if
C
C     Process 1 part
C
        if (procID.eq.1) then
          dest = 0
C       Send to process 1
          call MPI_SEND(coor,3*npt,MPI_DOUBLE_PRECISION,
     &              dest,DYNA_TAG,MPI_COMM_WORLD,MPIerr)

          call MPI_RECV(coor3 ,3*npt,MPI_DOUBLE_PRECISION,
     &      MPI_ANY_SOURCE,DYNA_TAG,MPI_COMM_WORLD,status,MPIerr)
        end if
C
C Compute the Rmsd      
C
        if (.false.) then
           j=0
           do k=1,npt
                 if (ptnm(k).eq.'CA  ')  then
                  j=j+1
                  CAindex(j) = k
                 end if
           end do
           
           do k=1, (totmon-2)
                   j=CAindex(k)
                   CA(1,k)=Coor(1,j)
                   CA(2,k)=Coor(2,j)
                   CA(3,k)=Coor(3,j)

                   CA3(1,k)=Coor3(1,j)
                   CA3(2,k)=Coor3(2,j)
                   CA3(3,k)=Coor3(3,j)
                   tmp_mass(k) = ptms(k)
                   ptms(k) = 1.d0
                   rms_pick(k) = k
           end do
        
        
           call rmsd_weight(totmon-2,CA,CA3,rms,.false.,ptms)

           do k=1, (totmon-2)
             ptms(k) = tmp_mass(k)
           end do
           
        else
        if (procID.lt.3) then
          call rmsd_weight(npt,coor,coor3,rms,.false.,ptms)
          write(*,*) 'Weighted Rmsd: ', rms
        end if
        end if

        end
