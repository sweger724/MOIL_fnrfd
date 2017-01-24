        integer function search(indx1,indx2,indx3,indx4,
     1              v1,v2,v3,v4,max,dim)
c
c search for the index i(=search on return) such that 
c       v[k=1,dim] = indx[k=1,dim]
c       max is the maximum length of the v's
c
        integer indx1,indx2,indx3,indx4,dim,max
        integer v1(*),v2(*),v3(*),v4(*)
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
c
c local
        integer i,level

        if (debug) then
         write(stdo,*)' in search dim max = ',dim,max
         do 11 i=1,max
          if (v1(i).ne.0) then
           write(stdo,*)' v1(',i,') ',v1(i),' v2(',i,') ',v2(i)
          end if
11       continue
        end if
        if (dim.eq.2) then
         do 1 i=1,max
          if (v1(i).eq.indx1 .and. v2(i).eq.indx2) then
           search = i
           return
          end if
1        continue
         write(stdo,100)indx1,indx2
100      format(1x,'cannot match bond for particle ids',2(1x,i7))
         level = 0
         call alert('search',6,'Bond not found',14,level)
        else if (dim.eq.3) then
         if (debug) write(stdo,*)' ANGLE SEARCH'
         do 2 i=1,max
          if (debug) then
           write(stdo,*)' v = ',v1(i),v2(i),v3(i)
          end if
          if ((v1(i).eq.indx1 .and. v2(i).eq.indx2 .and.
     1          v3(i).eq.indx3) .or. (v1(i).eq.indx3 .and.
     2          v2(i).eq.indx2 .and. v3(i).eq.indx1)) then
           search = i
           return
          end if
2        continue
         write(stdo,101)indx1,indx2,indx3
101      format(1x,'cannot match angle for particle ids',3(1x,i7))
         level = 1
         call alert('search',6,'Angle not found',14,level)
        else if (dim.eq.4) then
         do 3 i=1,max
           if (debug) then
            write(stdo,*)' indx1 indx2 indx3 indx4 '
            write(stdo,*)indx1,indx2,indx3,indx4
            write(stdo,*)'v1(i),v2(i),v3(i),v4(i)'
            write(stdo,*)v1(i),v2(i),v3(i),v4(i)
           end if
          if ((v1(i).eq.indx1 .and. v2(i).eq.indx2 .and.
     1          v3(i).eq.indx3 .and. v4(i).eq.indx4) .or.
     2          (v1(i).eq.indx4 .and. v2(i).eq.indx3 .and.
     3          v3(i).eq.indx2 .and. v4(i).eq.indx1)) then
           search = i
           return
          else if ((v1(i).eq.-999 .and. v2(i).eq.indx2 .and.
     1          v3(i).eq.indx3 .and. v4(i).eq.-999) .or.
     2          (v1(i).eq.-999 .and. v2(i).eq.indx3 .and.
     3          v3(i).eq.indx2 .and. v4(i).eq.-999)) then
           search = i
           return
          end if
3        continue
         write(stdo,102)indx1,indx2,indx3,indx4
102      format(1x,'cannot match torsion for particle ids',4(1x,i7))
         level = 0
         call alert('search',6,'Torsion not found',14,level)
         search = -1
        end if
        return
        end

        
        subroutine srch_imp(indx1,indx2,indx3,indx4,i1,i2,i3,i4,
     1              v1,v2,v3,v4,max,dim,id)
c
c search for the index i(=search on return) such that 
c done for improper torsions to get also the first and last
c indices correct
c       i1-i4 are the actual atomic indices which may be swapped
c               if needed (input and output)
c       v[k=1,dim] = indx[k=1,dim]
c       max is the maximum length of the v's
c
        integer indx1,indx2,indx3,indx4,dim,max
        integer i1,i2,i3,i4,id
        integer v1(*),v2(*),v3(*),v4(*)
        integer j
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
c
c local
        integer i,level

           if (debug) then
            write(stdo,*)' Improper torsion search '
            write(stdo,*)' indx1 indx2 indx3 indx4 '
            write(stdo,*)indx1,indx2,indx3,indx4
            write(stdo,*)'v1(i),v2(i),v3(i),v4(i)'
            do 10 i=1,max
            write(stdo,*)v1(i),v2(i),v3(i),v4(i)
10          continue
           end if
         do 1 i=1,max
          if (v1(i).eq.indx1) then
           if (v2(i).eq.indx2 .and. v3(i).eq.indx3
     1        .and. v4(i).eq.indx4) then
c             write(stdo,*) '1-2-3-4'
c no need to swap anything
                id = i
                return
           else if (v2(i).eq.indx2 .and. v3(i).eq.indx4
     1          .and. v4(i).eq.indx3) then
c             write(stdo,*)'1-2-4-3'
                j  = i4
                i4 = i3
                i3 = j
                id = i
                return
           else if (v2(i).eq.indx3 .and. v3(i).eq.indx4
     1          .and. v4(i).eq.indx2) then
c             write(stdo,*)'1-3-4-2'
                j  = i4
                i4 = i2
                i2 = i3
                i3 = j
                id = i
                return
           else if (v2(i).eq.indx3 .and. v3(i).eq.indx2
     1          .and. v4(i).eq.indx4) then
c             write(stdo,*)'1-3-2-4'
                j  = i2
                i2 = i3
                i3 = j
                id = i
                return
           else if (v2(i).eq.indx4 .and. v3(i).eq.indx2
     1          .and. v4(i).eq.indx3) then
c             write(stdo,*)'1-4-2-3'
                j  = i2
                i2 = i4
                i4 = i3
                i3 = j
                id = i
                return
           else if (v2(i).eq.indx4 .and. v3(i).eq.indx3
     1          .and. v4(i).eq.indx2) then
c             write(stdo,*)'1-4-3-2'
                j  = i2
                i2 = i4
                i4 = j
                id = i
                return
c       ileana -- adding the option of matching "X"atoms
          else if (v2(i).eq.-999.and.v3(i).eq.indx3
     1          .and. v4(i).eq.indx4) then
c               write(stdo,*)'1-X-3-4'
                id = i
                return
           else if (v2(i).eq.-999.and.v3(i).eq.indx4
     1          .and. v4(i).eq.indx3) then
c             write(stdo,*)'1-X-4-3'
                j  = i2
                i2 = i4
                i4 = j
                id = i
                return
         else if (v2(i).eq.-999.and.v3(i).eq.indx2
     1          .and. v4(i).eq.indx3) then
c           write(stdo,*)'1-X-2-3'
                j  = i2
                i2 = i4
                i4 = j
                id = i
                return
         else if (v2(i).eq.-999.and.v3(i).eq.indx2
     1          .and. v4(i).eq.indx4) then
c           write(stdo,*)'1-X-2-4'
                j  = i2
                i2 = i4
                i4 = j
                id = i
                return
         else if (v2(i).eq.-999.and.v3(i).eq.indx3
     1          .and. v4(i).eq.indx2) then
c           write(stdo,*)'1-X-3-2'
                j  = i2
                i2 = i4
                i4 = j
                id = i
                return
         else if (v2(i).eq.-999.and.v3(i).eq.indx4
     1          .and. v4(i).eq.indx2) then
c           write(stdo,*)'1-X-4-2'
                j  = i2
                i2 = i4
                i4 = j
                id = i
                return
          else if (v2(i).eq.-999.and.v3(i).eq.-999
     1          .and. v4(i).eq.indx4) then
c            write(stdo,*)'1-X-X-4'
                id = i
                return
         else if (v2(i).eq.-999.and.v3(i).eq.-999
     1          .and. v4(i).eq.indx3) then
c            write(stdo,*)'1-X-X-3'
                j  = i2
                i2 = i4
                i4 = j
                id = i
                return
         else if (v2(i).eq.-999.and.v3(i).eq.-999
     1          .and. v4(i).eq.indx2) then
c            write(stdo,*)'1-X-X-2'
                j  = i2
                i2 = i4
                i4 = j
                id = i
                return
         else if (v2(i).eq.indx2.and.v3(i).eq.-999
     1          .and. v4(i).eq.indx4) then
c          write(stdo,*)'1-2-X-4'
                id = i
                return
         else if (v2(i).eq.indx4.and.v3(i).eq.-999
     1          .and. v4(i).eq.indx2) then
c          write(stdo,*)'1-4-X-2'
                j  = i2
                i2 = i4
                i4 = j
                id = i
                return
         else if (v2(i).eq.indx3.and.v3(i).eq.-999
     1          .and. v4(i).eq.indx4) then
c            write(stdo,*)'1-3-X-4'
                j  = i2
                i2 = i4
                i4 = j
                id = i
                return
         else if (v2(i).eq.indx4.and.v3(i).eq.-999
     1          .and. v4(i).eq.indx3) then
c            write(stdo,*)'1-4-X-3'
                j  = i2
                i2 = i4
                i4 = j
                id = i
                return
         else if (v2(i).eq.indx2.and.v3(i).eq.-999
     1          .and. v4(i).eq.indx3) then
c            write(stdo,*)'1-2-X-3'
                j  = i2
                i2 = i4
                i4 = j
                id = i
                return
        else if (v2(i).eq.indx3.and.v3(i).eq.-999
     1          .and. v4(i).eq.indx2) then
c           write(stdo,*)'1-3-X-2'
                j  = i2
                i2 = i4
                i4 = j
                id = i
                return
        
            end if
C Alfredo. Include this to read the HEM information:
          else if (v3(i).eq.indx3.and.v1(i).eq.-999
     1            .and.v4(i).eq.-999.and.v2(i).eq.-999) then
c           write(stdo,*)'X-X-3-X'
                 id=i
                 return
          else if(v3(i).eq.indx2.and.v1(i).eq.-999
     1            .and.v4(i).eq.-999.and.v2(i).eq.-999) then
c           write(stdo,*)'X-3-X-X'
                 j=i3
                 i3=i2
                 i2=j
                 id=i
                 return
          else if(v3(i).eq.indx4.and.v1(i).eq.-999
     1            .and.v2(i).eq.-999.and.v4(i).eq.-999) then
c           write(stdo,*)'X-X-X-3'
                 j=i3
                 i3=i4
                 i4=j
                 id=i
                 return
c End of Alfredo modification                
          end if
1       continue

         id=max
         
         write(stdo,100)indx1,indx2,indx3,indx4
100      format(1x,'cannot match imptors for particle ids',4(1x,i7))
         level = 0
         call alert('search',6,'Imptors not found',14,level)
        
         return
        
        end
