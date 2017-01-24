        subroutine broad_conn()
c
c boradcast connectivity data,( MPI VERSION)
c  Modified by V. Zaloj on August 10, 1999
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/CONVERT.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/ACTPARA.BLOCK'
c
      include 'mpif.h'

C@
      integer zero,lenst,comm,ierr

      lenst = 4
C@
c zero is the index of the PE with all the data on
c

      zero = 0

      write(stdo,*) ' in broadcast '
      write(stdo,*) ' before broadcasting totmon = ',totmon
      CALL MPI_BCAST(totmon,1,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      write(stdo,*) ' in broadcast after sending totmon'
      write(stdo,*) ' totmon = ',totmon
      CALL MPI_BCAST(npt,1,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(nb,1,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(nangl,1,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(ntors,1,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(nimp,1,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(lestyp,1,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(nbulk,1,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(nbeta,1,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(totex,1,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(totspe,1,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(nchg,1,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(totex3,1,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(totspe3,1,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(lesflag,1,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(hydroflag,1,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(arith,1,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)

      write(stdo,*) ' Fisrt conn'

      CALL MPI_BCAST(poimon,npt,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(lesid,npt,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(cplbl,npt,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(flagchr,npt,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)

      write(stdo,*) ' Second conn'

c      CALL MPI_BCAST(ptnm,(lenst*npt),MPI_CHARACTER,zero,
c     1  MPI_COMM_WORLD, ierr)

      write(stdo,*) ' Third conn'

      CALL MPI_BCAST(ptms,npt,MPI_DOUBLE_PRECISION,zero,
     1  MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(ptchg,npt,MPI_DOUBLE_PRECISION,zero, 
     1  MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(epsgm6,npt,MPI_DOUBLE_PRECISION,zero,
     1    MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(epsgm12,npt,MPI_DOUBLE_PRECISION,zero,
     1    MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(ptwei,npt,MPI_DOUBLE_PRECISION,zero, 
     1   MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(invms,npt,MPI_DOUBLE_PRECISION,zero, 
     1   MPI_COMM_WORLD, ierr)
c      CALL MPI_BCAST(bulk,(lenst*nbulk),MPI_CHARACTER,zero, 
c     1   MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(pbulk,nbulk,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
c      CALL MPI_BCAST(moname,(lenst*totmon),MPI_CHARACTER,zero, 
c     1   MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(poipt,(totmon+1),MPI_INTEGER,zero, 
     1   MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(betap,totmon,MPI_INTEGER,zero, 
     1   MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(betap,totmon,MPI_INTEGER,zero,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(cbeta,totmon,MPI_DOUBLE_PRECISION,zero, 
     1   MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(avg_hydro,1,MPI_DOUBLE_PRECISION,zero, 
     1   MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(hydro_scale,1,MPI_DOUBLE_PRECISION,zero, 
     1   MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(ib1,nb,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(ib2,nb,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(imb1,nmb,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(imb2,nmb,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(kbond,nb,MPI_DOUBLE_PRECISION,zero, 
     1   MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(req,nb,MPI_DOUBLE_PRECISION,zero, 
     1   MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(rmeq,nmb,MPI_DOUBLE_PRECISION,zero, 
     1   MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(d,nmb,MPI_DOUBLE_PRECISION,zero, 
     1   MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(alpha,nmb,MPI_DOUBLE_PRECISION,zero, 
     1   MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(arep,nmb,MPI_DOUBLE_PRECISION,zero, 
     1   MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(beta1,nmb,MPI_DOUBLE_PRECISION,zero, 
     1   MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(brep,nmb,MPI_DOUBLE_PRECISION,zero, 
     1   MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(iangl1,nangl,MPI_INTEGER,zero,MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(iangl2,nangl,MPI_INTEGER,zero,MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(iangl3,nangl,MPI_INTEGER,zero,MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(kangl,nangl,MPI_INTEGER,zero,MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(angleq,nangl,MPI_INTEGER,zero,MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(itor1,ntors,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(itor2,ntors,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(itor3,ntors,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(itor4,ntors,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(period,ntors,MPI_DOUBLE_PRECISION,zero, 
     1   MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(ktors1,ntors,MPI_DOUBLE_PRECISION,zero, 
     1   MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(ktors2,ntors,MPI_DOUBLE_PRECISION,zero, 
     1   MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(ktors3,ntors,MPI_DOUBLE_PRECISION,zero, 
     1   MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(phase,ntors,MPI_DOUBLE_PRECISION,zero, 
     1   MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(iimp1,ntors,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(iimp2,ntors,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(iimp3,ntors,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(iimp4,ntors,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(kimp,ntors,MPI_DOUBLE_PRECISION,zero, 
     1   MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(impeq,ntors,MPI_DOUBLE_PRECISION,zero, 
     1   MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(exc1,(npt+1),MPI_INTEGER,zero,MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(exc2,totex,MPI_INTEGER,zero,MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(spec1,(npt+1),MPI_INTEGER,zero,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(spec2,totspe,MPI_INTEGER,zero,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(p14,(3*totspe),MPI_DOUBLE_PRECISION,zero, 
     1   MPI_COMM_WORLD, ierr)


c energy parameters
c
      CALL MPI_BCAST(ctrue,1,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(ebyes,1,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(ethyes,1,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(etoyes,1,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(eimyes,1,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(evdyes,1,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(eelyes,1,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(ecnyes,1,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(esymyes,1,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(emyes,nmb,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(eteth_yes,1,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(lcent,1,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(nocut,1,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(shift,1,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(efyes,1,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(ehyes,1,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(specl,1,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(repyes,nmb,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(repyes0,1,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(emyes0,1,MPI_INTEGER,zero, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(eps,1,MPI_DOUBLE_PRECISION,zero,
     1   MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(e0fin,1,MPI_DOUBLE_PRECISION,zero, 
     1  MPI_COMM_WORLD, ierr)

      return
      end



