        subroutine cmprs(vector,away,length)
c
c This routine check the element of the integer vector  "vector"
c if an element is equal "away", it is deleted and the vector is 
c compressed. length is modified in output
c
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        integer vector(*),away,length
        integer i,j

        if (length.eq.0) return

        if (debug) then
         write(stdo,*)' in VCMPS before action '
         write(stdo,*)' length away ',length,away
         write(stdo,*)' vector = ',(vector(i),i=1,length)
        end if
        if (vector(length).eq.away .and.length.eq.1) then
         length = 0
         return
        end if
        if (vector(length).eq.away) length = length-1
        i   = 0
1       continue
        i = i + 1
        if (vector(i).eq.away) then
         do 2 j=i,length-1
          vector(j) = vector(j+1)
2        continue
         length = length - 1
         i = i - 1
        end if
        if (i .eq. length) then
         if (debug) then
          write(stdo,*)' in VCMPS after action '
          write(stdo,*)' length away ',length,away
          write(stdo,*)' vector = ',(vector(i),i=1,length)
         end if
         return
        end if
        go to 1
        end

        subroutine cmprs2(vector1,vector2,away,length)
c
c This routine check the element of the integer vectors  "vector1/2"
c if an element is equal "away", it is deleted and the vectors are
c compressed. length is modified in output
c
        include 'COMMON/UNITS.BLOCK'
        include 'COMMON/DEBUG.BLOCK'
        integer vector1(*),vector2(*),away,length
        integer i,j

        if (length.eq.0) return
        i   = 0
        if (debug) then
         write(stdo,*) ' length away ',length,away
        end if

1       continue
        i = i + 1
        if (vector1(i).eq.away .or. vector2(i).eq.away) then
         do 2 j=i,length
          vector1(j) = vector1(j+1)
          vector2(j) = vector2(j+1)
2        continue
         length = length - 1
         i = i - 1
        end if
        if (i .eq. length) return
        go to 1
        end
