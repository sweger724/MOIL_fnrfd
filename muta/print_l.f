      subroutine print_l (D_E_L,sD_E_L,lambda_max,n_lambda,lambda_step,
     &                    in_lambda,f_lambda,n_ens,BandF)

      integer n_lambda,lambda_max,ilambda,n_ens
      double precision D_E_L(lambda_max),lambda,lambda_step,
     &                 sD_E_L(lambda_max),df,in_lambda,f_lambda,du
      logical BandF,back

      open(83,FILE='lambda_res.log',STATUS='NEW',FORM='FORMATTED')

      df=0.0d0
      du=0.0d0

      back=.FALSE.

      write(83,*) lambda_step

      if (BandF) lambda_step=-lambda_step

      lambda=in_lambda-lambda_step
      
      DO ilambda=1,n_lambda

      lambda = lambda_step + lambda

      if (abs(lambda-f_lambda).lt.10.d0**(-6)) then

       lambda_step=-lambda_step

       back=.TRUE.

      endif

      D_E_L(ilambda) =D_E_L(ilambda)/dble(n_ens)

      sD_E_L(ilambda)=sD_E_L(ilambda)/dble(n_ens) - D_E_L(ilambda)**2

      sD_E_L(ilambda)=(sD_E_L(ilambda)/dble(n_ens))**(0.5)

      df = df + D_E_L(ilambda)*lambda_step

      if (.not.back) then

       du = du + D_E_L(ilambda)*lambda_step

      endif

      write(83,1000) lambda,D_E_L(ilambda),sD_E_L(ilambda)

      ENDDO

      write(83,1010) df,du,df/du

1000   FORMAT('<U1-U2>_',g20.10,'=',g20.10,'+/-',g20.10)
1010   FORMAT('df(O)=',g20.10,'du(->)=',g20.10,'df(O)/du(->)=',g20.10)

      return

      end
