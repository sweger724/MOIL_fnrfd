        subroutine Andersen_thermostat(dr)

        implicit none

         include 'COMMON/LENGTH.BLOCK'
         include 'COMMON/CONNECT.BLOCK'
         include 'COMMON/VELOC.BLOCK'
         include 'COMMON/DYNA.BLOCK'
         include 'COMMON/NBLIST.BLOCK'
         
         real rrr(maxmono)
         double precision dr, vvv(6), mo, mh, S
         integer i,l, io, ih1, ih2
c@
         integer k

         call RANLUX(rrr,nwaters)
         k = 0

         io = dpoipt(idxtip3(1))-2
         ih1 = io + 1
         ih2 = io + 2
         mo = ptms(io)/(ptms(io)+2*ptms(ih1))
         mh = ptms(ih1)/(ptms(io)+2*ptms(ih1))

         do i=1,nwaters
           if (dr.gt.rrr(i)) then
            k = k+1
            io = dpoipt(idxtip3(i))-2
            ih1=io+1
            ih2=io+2

            do l = 1,3
              vvv(l) = velo(l,io)*mo + (velo(l,ih1) + velo(l,ih2))*mh
              
              velo(l,io)  = velo(l,io) - vvv(l)
              velo(l,ih1) = velo(l,ih1) - vvv(l)
              velo(l,ih2) = velo(l,ih2) - vvv(l)
            end do
           
            call draw_normal(vvv,3, tempi(1) / (ptms(io)+2*ptms(ih1)) ) 
            
            do l = 1,3
              velo(l,io)  = velo(l,io) + vvv(l)
              velo(l,ih1) = velo(l,ih1) + vvv(l)
              velo(l,ih2) = velo(l,ih2) + vvv(l)
            end do

           end if
         end do
         return
        end
