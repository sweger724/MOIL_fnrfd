        subroutine fill_torsion(k1,k2,k3,k4,i,weight)
        implicit none
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/TMP_CONNECT.BLOCK'

        integer k1,k2,k3,k4,i
        double precision weight

        tntors = tntors + 1

        titor1(tntors) = k1
        titor2(tntors) = k2
        titor3(tntors) = k3
        titor4(tntors) = k4
        tperiod(tntors) = period(i)
        tktors1(tntors) = ktors1(i)*weight
        tktors2(tntors) = ktors2(i)*weight
        tktors3(tntors) = ktors3(i)*weight
        tphase1 (tntors) = phase1(i)
        tphase2 (tntors) = phase2(i)
        tphase3 (tntors) = phase3(i)

        return
        end
