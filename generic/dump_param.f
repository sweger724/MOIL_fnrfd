        subroutine dump_param(stdo)
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONNECT.BLOCK'
        include 'COMMON/ENERGY.BLOCK'
        include 'COMMON/PARALLEL.BLOCK'

        integer stdo

        write(stdo,1)
1       format(////,1x, ' **** CONNECTIVITY PARAMETERS ')

        write(stdo,2)lesflag,hydroflag,arith,specl,prll_on_off
2       format(1x,' LES IS ON ? ',l1,/,1x,' Hydrophobicity ? ',l1,/
     1        ,1x,' Arith.aver? ',l1,/,1x,' Landau Zener   ? ',l1,/
     2        ,1x,' Parallel ?  ',l1,/)

        if (prll_on_off) write(stdo,3)my_pe
3     format(1x,'*The number of bonds angles torsions impropers 1-4',/
     1        ,1x,'* is for the current processor = ',i3,/
     2        ,1x,'* The rest are global variables.',/) 
     
        write(stdo,4)totmon,npt,nb,nangl,ntors,nimp,lestyp,
     1          nbulk,totex,totspe,nmb,nchg,nbeta
4       format( 1x,' total number of monomers  ',i6,/,
     1          1x,' total number of particles ',i6,/,
     2          1x,' total number of bonds     ',i6,/,
     3          1x,' total number of angles    ',i6,/,
     4          1x,' total number of torsions  ',i6,/,
     5          1x,' total number of impropers ',i6,/,
     6          1x,' number of les types       ',i6,/,
     7          1x,' number of Mol. systems    ',i6,/,
     8          1x,' number of non-bond excl.  ',i6,/,
     9          1x,' number of 1-4 interact.   ',i6,/,
     1          1x,' number of Morse bonds     ',i6,/,
     2          1x,' number of charged prtc    ',i6,/,
     3          1x,' number of Beta carbon     ',i6,/)

        if (specl) write(stdo,5)totex3,totspe3
5       format(1x,' Excited State (LZ model) ',/,
     1' number of exclusions at the exc. state ',i5,/,
     2' number of 1-4 inter. at the exc. state ',i5,/)

        write(stdo,6)
6       format(1x,' **** ENERGY PARAMETERS ')
        write(stdo,7)ebyes,ethyes,etoyes,eimyes,evdyes,eelyes,esymyes,
     1  ecnyes,emyes0,eteth_yes,efyes,ehyes,eenmyes,eCGyes,lcent
7     format(1x,'LOGICAL FLAGS FOR DIFFERENT ENERGIES (True/False)',/,
     1 1x,' bonds=',l1,/,' angles=',l1,/,' torsions=',l1,/,
     2 1x,' impropers=',l1,/,' van-der Waals=',l1,/,' electro =',l1,/,
     3 1x,' symmetry =',l1,/,' constraints  =',l1,/,' Morse   =',l1,/,
     4 1x,' tether =',l1,/,' finite vdw   =',l1,/,' Hydropho=',l1,/,
     5 1x,' ENM =',l1,/,' CG   =',l1,/,' center =',l1,//)
        if (nocut) write(stdo,8)
8      format(1x,' **NO BUFFER IS USED IN NON-BONDED CALCULATION',//) 
        return
        end
