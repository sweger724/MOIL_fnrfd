      subroutine buildmolec()
c
c Subroutine to built a molecule from the rotamer already included. 
c All it does is to move the atoms from the end to the correct position,
c fixing the other indeces and names. 
c 
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/CONNECT_SHORT.BLOCK'
      include 'COMMON/CONSPECL4.BLOCK'
      include 'COMMON/COORD.BLOCK'
      include 'COMMON/ROTAGEN.BLOCK'
      include 'COMMON/DEBUG.BLOCK'
      include 'COMMON/ROTAINT.BLOCK'
      include 'COMMON/ENERGY_DEE.BLOCK'
      
      character*12 name
      integer namel,level

      double precision coor_temp(3,maxpt)
      integer newindex(maxpt)
      integer poipt_new(0:maxmono),iposenh,imono
      integer i,j,k,irot,typerot
      integer nptnew
c
      integer poimon_temp(maxpt),ptid_temp(maxpt),poimonold
      logical flagchr_temp(maxpt)
      character*4 ptnm_temp(maxpt),mononewname
      double precision ptms_temp(maxpt),ptchg_temp(maxpt)
      double precision epsgm6_temp(maxpt),epsgm12_temp(maxpt)
c
      name = 'buildmolec'
      namel = 10
c
      nptnew=0
c
c
      do i=1,npt
         do j=1,3
            coor_temp(j,i)=coor(j,i)
         end do
         newindex(i)=0
         poimon_temp(i)=poimon(i)
         ptid_temp(i)=ptid(i)
         ptnm_temp(i)=ptnm(i)
         ptms_temp(i)=ptms(i)
         ptchg_temp(i)=ptchg(i)
         epsgm6_temp(i)= epsgm6(i)
         epsgm12_temp(i)=epsgm12(i)
         flagchr_temp(i)=flagchr(i)

      end do
c
c
c
      do i=1,totmon
         if (ipickm(i).eq.1)then
c
            iposenh=iposenh+1
c...........copy the first part         
            do j=poipt(i-1)+1,poipt(i)
               nptnew=nptnew+1
               newindex(j)=nptnew
               do k=1,3
                  coor(k,nptnew)=coor_temp(k,j)
               end do
c
               poimon(nptnew)=poimon_temp(j)
               ptid(nptnew)=ptid_temp(j)
               ptnm(nptnew)=ptnm_temp(j)
               ptms(nptnew)=ptms_temp(j)
               ptchg(nptnew)=ptchg_temp(j)
               epsgm6(nptnew)= epsgm6_temp(j)
               epsgm12(nptnew)=epsgm12_temp(j)
               flagchr(nptnew)=flagchr_temp(j)
c
            end do
c
            poimonold=poimon_temp(poipt(i))
c
c...........copy the second part            
            irot=poirotenh(iposenh)
            typerot=typerotenh(irot)
            do j=poiatrotenh(irot-1)+1,poiatrotenh(irot)
               nptnew=nptnew+1
               newindex(j)=nptnew
               do k=1,3
                  coor(k,nptnew)=coor_temp(k,j)
               end do
c
               poimon(nptnew)=poimonold
               ptid(nptnew)=ptid_temp(j)
               ptnm(nptnew)=ptnm_temp(j)
               ptms(nptnew)=ptms_temp(j)
               ptchg(nptnew)=ptchg_temp(j)
               epsgm6(nptnew)= epsgm6_temp(j)
               epsgm12(nptnew)=epsgm12_temp(j)
               flagchr(nptnew)=flagchr_temp(j)
c
            end do
c

         else
c
c...........just copy            
            do j=poipt(i-1)+1,poipt(i)
               nptnew=nptnew+1
               newindex(j)=nptnew
               do k=1,3
                  coor(k,nptnew)=coor_temp(k,j)
               end do
c
               poimon(nptnew)=poimon_temp(j)
               ptid(nptnew)=ptid_temp(j)
               ptnm(nptnew)=ptnm_temp(j)
               ptms(nptnew)=ptms_temp(j)
               ptchg(nptnew)=ptchg_temp(j)
               epsgm6(nptnew)= epsgm6_temp(j)
               epsgm12(nptnew)=epsgm12_temp(j)
               flagchr(nptnew)=flagchr_temp(j)
c
            end do
c
         end if
c
         poipt_new(i)=nptnew
c
      end do
c
      if (nptnew.ne.npt)then
         level = 1
         call alert(name,namel,'wrong number of atoms',21,level)
      end if
c
      do i=1,totmon
         poipt(i)=poipt_new(i)
      end do
c
      do i=1,nposenh
c
c........find the monomer type used at the position. 
         irot=poirotenhaux(i)
         typerot=typerotenh(rotaux(irot))
c
         imono=indexposenh(i)
         mononewname=monames4(typerot)
         moname(imono)=mononewname(1:3)
c                
c
      end do
c
c              
      return
      end

c






