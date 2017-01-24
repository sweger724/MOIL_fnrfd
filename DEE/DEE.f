      subroutine DEE(niterations)
c
c Eliminates rotamers using the Dead-End-Elimination algorithm. 
c In this routine we use the inequality for rotamers i and j. 
c    (Goldstein, Biophysical J. 66,1335(1994)
c    E(g_i)-E(h_i)+SUM_{j}{ min_{f'} [ E(g_i,f'_j)-E(h_i,f'_j)]}
c
c In this subroutine we use, as an improvement on the above equation,
c that dead-end-pairs (g_i,f'_j) do not contribute to the sum. 
c In practice we used (as proposed by Desmet et all, (1995)) that the 
c energy in this case becomes very large. 
c 
c 
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/UNITS.BLOCK'
      include 'COMMON/CONNECT_SHORT.BLOCK'
      include 'COMMON/CONSPECL4.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ROTAINT.BLOCK'
      include 'COMMON/ROTAGEN.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/ENERGY.BLOCK'
      include 'COMMON/ENERGY_DEE.BLOCK'

      
      character*3 name
      integer namel,level
c
      double precision dtotal1,dtotal2,ediffmin1,ediffmin2
      double precision Eij1,Eij2
      integer i,j,removed,removed_loc,ind1,ind2,niterations,rotlastreal
      integer rot1,rot2,rotfirst,rotlast,rot1real,rot2real
      integer rotT,rotfirstT,rotlastT,rotTreal
      integer pointEij1,pointEij2,pointEijT
      integer nmonenhstart
c
      name = 'DEE'
      namel = 3
      level=1
c
c
      removed=999
      niterations=0
      nmonenhstart=nmonenhleft
c
      write (6,*) 
      write (6,*) 'DEE procedure:'
c
c
      do while (removed.gt.0)
c
         removed=0
         niterations=niterations+1
c
c........loop over all positions
         do i=1,nposenh
c     
            rotfirst=poirotenhaux(i-1)+1
            rotlast=poirotenhaux(i)
            removed_loc=0
c...........check if there is more than one left
            if (rotlast.gt.rotfirst) then
c
c
               do rot1=rotfirst,rotlast-1
c
                  rot1real=rotaux(rot1)
                  pointEij1=pointEij(rot1real-1)-poirotenh(i)
c
c.................check if it was not removed in a previous iteraction 
c.................of this loop.
                  if (kept(rot1real)) then
c
                     do rot2=rot1+1,rotlast
c
                        rot2real=rotaux(rot2)
                        pointEij2=pointEij(rot2real-1)-poirotenh(i)
c
c.......................check if it was not removed in a previous 
c.......................iteraction of this loop (check rot1 again).
                        if (kept(rot2real).and.kept(rot1real)) then
c
                           dtotal1=Eiback(rot1real)-Eiback(rot2real)
                           dtotal2=Eiback(rot2real)-Eiback(rot1real)
c     
c
                           do j=1,i-1
c     
                              ediffmin1=1.0d30
                              ediffmin2=1.0d30
c     
                              rotfirstT=poirotenhaux(j-1)+1
                              rotlastT=poirotenhaux(j)
c
                              do rotT=rotfirstT,rotlastT
                                 rotTreal=rotaux(rotT)
                                 pointEijT=pointEij(rotTreal-1)-
     &                                poirotenh(j)
c
                                 ind1=pointEijT+rot1real
                                 ind2=pointEijT+rot2real
c
                                 Eij1=Eij(ind1)
                                 Eij2=Eij(ind2)
c
                                 if(DEPij(ind1))then
                                    ediffmin1=dmin1(ediffmin1,
     &                                   hugeE-Eij2)
                                 else
                                    ediffmin1=dmin1(ediffmin1,Eij1-Eij2)
                                 end if
c
                                 if(DEPij(ind2))then
                                    ediffmin2=dmin1(ediffmin2,
     &                                   hugeE-Eij1)
                                 else
                                    ediffmin2=dmin1(ediffmin2,Eij2-Eij1)
                                 end if
c
                              end do
c     
                              dtotal1=dtotal1+ediffmin1
                              dtotal2=dtotal2+ediffmin2
c                     
                           end do
c     
                           do j=i+1,nposenh
c     
                              ediffmin1=1.0d30
                              ediffmin2=1.0d30
c     
                              rotfirstT=poirotenhaux(j-1)+1
                              rotlastT=poirotenhaux(j)
c     
                              do rotT=rotfirstT,rotlastT
                                 rotTreal=rotaux(rotT)
                                 ind1=pointEij1+rotTreal
                                 ind2=pointEij2+rotTreal
c
                                 Eij1=Eij(ind1)
                                 Eij2=Eij(ind2)
c
                                 if(DEPij(ind1))then
                                    ediffmin1=dmin1(ediffmin1,
     &                                   hugeE-Eij2)
                                 else
                                    ediffmin1=dmin1(ediffmin1,Eij1-Eij2)
                                 end if
c
                                 if(DEPij(ind2))then
                                    ediffmin2=dmin1(ediffmin2,
     &                                   hugeE-Eij1)
                                 else
                                    ediffmin2=dmin1(ediffmin2,Eij2-Eij1)
                                 end if
c
                              end do
c
                              dtotal1=dtotal1+ediffmin1
                              dtotal2=dtotal2+ediffmin2     
c     
                           end do
c
                           if (dtotal1.gt.0) then
                              kept(rot1real)=.false.
                              removed = removed+1
                              removed_loc=removed_loc+1
                           else if (dtotal2.gt.0) then
                              kept(rot2real)=.false.
                              removed = removed+1
                              removed_loc=removed_loc+1
                           end if
c
c
                        end if
c
                     end do
c
                  end if
               
               end do
c
               if (removed_loc.eq.rotlast-rotfirst+1) then
                  write (6,*) 'Problems with clashes with the ',
     &                 'backbone at position', i
                  call alert(name,namel,'all rotamers discarded ',
     &                 23,level)
               end if
c
c..............ACTUALLY ELIMINATE THE ROTAMERS
c
               if (removed_loc.gt.0) then
                  call removerot(i,rotfirst,rotlast)
               end if
c
c
            end if
c
         end do
c
         write (6,*) '  finished iteration',niterations
c
      end do
c
      write (6,*) '  #iteration',niterations
      write (6,*) '  #rotamers start',nmonenhstart
      write (6,*) '  #rotamers end  ',nmonenhleft
      write (6,*) '  final residues:'
      do j=1,nposenh
         rotfirst=poirotenhaux(j-1)+1
         rotlast=poirotenhaux(j)
         write (6,40) j,rotlast-rotfirst+1
         write (6,50)(intindrotenh(rotaux(i)),
     &        monames4(typerotenh(rotaux(i))),
     &        i=rotfirst,rotlast)
      end do
CDEB
CDEB      write (6,*) 'all residues flags'
CDEB      write (6,*)(kept(i),i=1,nmonenh)
c
 40   format(2x,'position #',i4,' (',i4,' rotamers)')
 50   format(8(1x,i4,1x,a4))
c
      return
      end
c     
c







