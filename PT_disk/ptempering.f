      subroutine ptempering(e_rep,t_rep,accept)
        implicit none 
        
        include 'COMMON/LENGTH.BLOCK'
        include 'COMMON/CONVERT.BLOCK'
        include 'COMMON/PT.BLOCK'

        double precision e_rep(2),t_rep(2)   
        double precision prob,betas,deltas
        logical accept

        accept=.false.
                
        betas = ( 1./t_rep(2) - 1./t_rep(1) )/kboltzmann
        deltas = ( e_rep(2)-e_rep(1) )
        prob = exp( deltas*betas )

        if(prob.ge.RANDOM_ARRAY(random_index)) accept=.true.
        !if (replicaID.eq.1)
         write(6,*)"Exchange with probability", prob, 
     &     RANDOM_ARRAY(random_index) 
         write(6,*)"Energies to switch",e_rep(1),e_rep(2)
     
        return
      end
