        function getCGid(name)
        implicit none
          character*(*) name
          integer getCGid

         getCGid = 6
         
         if (name .eq."HISZ"  ) getCGid = 1
         if (name .eq."ARGZ"  ) getCGid = 2
         if (name .eq."ASPZ"  ) getCGid = 3
         if (name .eq."ASNZ"  ) getCGid = 4
         if (name .eq."GLYZ"  ) getCGid = 5
         if (name .eq."ALAZ"  ) getCGid = 6
         if (name .eq."PROZ"  ) getCGid = 7
         if (name .eq."TYRZ"  ) getCGid = 8
         if (name .eq."TRPZ"  ) getCGid = 9
         if (name .eq."CYSZ"  ) getCGid = 10
         if (name .eq."LEUZ"  ) getCGid = 11
         if (name .eq."ILEZ"  ) getCGid = 12
         if (name .eq."METZ"  ) getCGid = 13
         if (name .eq."VALZ"  ) getCGid = 14
         if (name .eq."LYSZ"  ) getCGid = 15
         if (name .eq."GLUZ"  ) getCGid = 16
         if (name .eq."GLNZ"  ) getCGid = 17
         if (name .eq."SERZ"  ) getCGid = 18
         if (name .eq."THRZ"  ) getCGid = 19
         if (name .eq."PHEZ"  ) getCGid = 20
         if (name .eq."CGTR"  ) getCGid = 5

         if (name .eq. "H") getCGid = 1
         if (name .eq. "E") getCGid = 2
         if (name .eq. "T") getCGid = 3
         if (name .eq. "F") getCGid = 4
         if (name .eq. "P") getCGid = 5
         if (name .eq. "?") getCGid = 6

         return 

        end
