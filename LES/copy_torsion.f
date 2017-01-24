        subroutine copy_torsions()
        implicit none
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/TMP_CONNECT.BLOCK'

        integer i

        ntors = tntors

        do 1 i=1,tntors

                itor1(i)=titor1(i)
                itor2(i)=titor2(i)
                itor3(i)=titor3(i)
                itor4(i)=titor4(i)
                period(i)=tperiod(i)
                ktors1(i)=tktors1(i)
                ktors2(i)=tktors2(i)
                ktors3(i)=tktors3(i)
                phase1(i)=tphase1(i)
                phase2(i)=tphase2(i)
                phase3(i)=tphase3(i)
1       continue

        return
        end
