        subroutine Write2File(ufile,npt,rr,ee,massfac,MASSWEIGHT)
          
          implicit none

          integer j,npt,ufile,MASSWEIGHT
          double precision rr(3,*),ee,massfac(*)

          if (MASSWEIGHT.eq.1) then
            write(ufile) ee,
     1            (rr(1,j)*massfac(j),j=1,npt),
     2            (rr(2,j)*massfac(j),j=1,npt),
     3            (rr(3,j)*massfac(j),j=1,npt)
          else
            write(ufile) ee,
     1            (rr(1,j),j=1,npt),
     2            (rr(2,j),j=1,npt),
     3            (rr(3,j),j=1,npt)

          end if

        return
        end
