        subroutine homology_force()
        implicit none
c
c  Add a penalty to the energy term for different 3D structure
c  first coordinates of alligned CA atoms are comunicated
c  then rmsd of these CA atoms is computed and if rmsd is bigger
C  than R_tresh an aditional force_const.(1-R_tresh/RMSD) gradient(RMSD^2)  
C  is added to CA atoms

C  when computing gradient of E_H the rotation matrix used in
c  RMSD computation is considered as a constant matrix... we
c  checked this in the experiments and it looks as a fair assumption

        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/COORD.BLOCK'
        include 'COMMON/DYNA.BLOCK'
        include 'COMMON/COUPLED_DYNA.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'
        include 'COMMON/OVERLAP.BLOCK'
        include 'mpif.h'

C       Local variables
        integer status(MPI_STATUS_SIZE)
        double precision  coor2CA(3,maxmono)
        integer j, i, l
        integer dest
        double precision tmp, tmp2(3), tmp_mass(maxpt)
C
C  Prepare the data
C
C       write(*,*),'******** AlignN: ', alignN
        do j=1,alignN
              coorCA(1,j)=coor(1,CAindex(j))
              coorCA(2,j)=coor(2,CAindex(j))
              coorCA(3,j)=coor(3,CAindex(j))
        end do
        write(*,*) 'coorCA(3,15): ', coorCA(3,15)
C
C  Comunicate the data
C
C     Process 0 part
C
        if (procID.eq.0) then
          call MPI_RECV(coor2CA ,3*alignN,MPI_DOUBLE_PRECISION,
     &      MPI_ANY_SOURCE,DYNA_TAG,MPI_COMM_WORLD,status,MPIerr)

C       Send to process 0
          dest = 1 
          call MPI_SEND(coorCA,3*alignN,MPI_DOUBLE_PRECISION,
     &              dest,DYNA_TAG,MPI_COMM_WORLD,MPIerr)
        end if
C
C     Process 1 part
C
        if (procID.eq.1) then
          dest = 0
C       Send to process 1
          call MPI_SEND(coorCA,3*alignN,MPI_DOUBLE_PRECISION,
     &              dest,DYNA_TAG,MPI_COMM_WORLD,MPIerr)

          call MPI_RECV(coor2CA ,3*alignN,MPI_DOUBLE_PRECISION,
     &      MPI_ANY_SOURCE,DYNA_TAG,MPI_COMM_WORLD,status,MPIerr)
        end if

C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
C       write(*,*) '============================='
C       write(*,*) 'Communication: '
C       write(*,*) 'procID: ', procID
C       write(*,*) 'coorCA(3,2): ', coorCA(3,2)
C       write(*,*) 'coor2CA(3,2): ', coor2CA(3,2)
C       write(*,*) '============================='
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

C
C Compute the homology force    
C
        do j=1,alignN
          tmp_mass(j) = ptms(j)
          ptms(j) = 1.d0
          rms_pick(j) = j
        end do
        
        call rmsd_weight(alignN,coorCA,coor2CA,rms,.false.,ptms)
        
        do j=1,alignN
          ptms(j) = tmp_mass(j)
        end do
        
        if( rms.gt.R_tresh) then
          tmp = -2.0d0 * E_homology_const *(1- R_tresh/rms) / alignN
          write(*,*) 'Homology force: ', tmp
          do i=1,alignN
             tmp2(1) = (coor2CA(1,i)-coorCA(1,i)) * tmp
             tmp2(2) = (coor2CA(2,i)-coorCA(2,i)) * tmp
             tmp2(3) = (coor2CA(3,i)-coorCA(3,i)) * tmp
C
C   Aply force to  atoms of l-th residues
C            
             l=align(i)
             j=CAindex(i)
          
             dpot(1,j)= dpot(1,j) + tmp2(1) *ptms(j) 
             dpot(2,j)= dpot(2,j) + tmp2(2) *ptms(j)
             dpot(3,j)= dpot(3,j) + tmp2(3) *ptms(j)
             
          end do
          
          write(*,*) 'Force: ', tmp2(1)*ptms(CAindex(2))
        end if

        end
