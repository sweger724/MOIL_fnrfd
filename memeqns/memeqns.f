        program memeqns

        implicit none
        include 'MEMEQNS.BLOCK'
        integer mi,nmlst
        double precision omfpt, mfpt(maxmlst)
        double precision frac_p(maxmlst), b
        double precision fpt(maxmlst,maxtraj)
        double precision kp(maxmlst,maxbins)
        double precision km(maxmlst,maxbins)
        double precision p(maxmlst,maxqk)
        double precision q(maxmlst,maxqk)
        integer nfpt(maxmlst), nq, nbins
        character*128 allfpts

        call query_user ( allfpts, nmlst, mi, nq, b )

        call read_fpts ( allfpts, fpt, nmlst, nfpt )
        call fpts2frac ( fpt,nmlst,nfpt,mfpt,frac_p )
        call check_fpts ( frac_p, nmlst )

        if ( mode .eq. 'EQUIL' ) then
           call print_equilibrium ( frac_p, mfpt, nmlst )
        else
           call shallmfpt ( frac_p, mfpt, nmlst, mi, omfpt )
           write(*,1) omfpt
1          format(1x,' MFPT: ',f9.3)
           if ( mode .eq. 'MFPT2' ) then
              call rshallmfpt ( frac_p, mfpt, nmlst, mi, omfpt )
              write(*,2) omfpt
2             format(1x,' Reverse MFPT: ',f9.3)
           end if
        end if

        if ( nq .eq. 0 ) stop

        call build_histograms(fpt,nmlst,nfpt,b,kp,km,nbins)
        call evolveq ( kp, km, mi, nq, q, nmlst, nbins )
        call write_array ( q, nmlst, nq, 10 )
        call q2p ( q, kp, km, mi, nmlst, nq, nbins, p )
        call write_array ( p, nmlst, nq, 11 )

        stop
        end






        subroutine q2p ( q, kp, km, mi, nmlst, nq, nbins, p )

        implicit none
        include 'MEMEQNS.BLOCK'
        integer m, j, t, u, nmlst, nq, nbins, mi
        double precision q(maxmlst,*)
        double precision p(maxmlst,*)
        double precision kp(maxmlst,*)
        double precision km(maxmlst,*)
        double precision cpt(maxmlst,maxqk)

c       do 1000 m = 1, nmlst
c          p(m,1) = 0.d0
c 1000  continue


        do 1 m = 1, nmlst
           cpt(m,1) = 1.0
           do 2 j = 1, nq
              cpt(m,j+1) = cpt(m,j)
              if ( j .le. nbins ) then
                 cpt(m,j+1) = cpt(m,j+1) - ( kp(m,j) + km(m,j) )
              end if
 2         continue
 1      continue


        do 11 t = 1, nq
           do 12 m = 1, nmlst
              p(m,t) = 0.d0
              do 13 j = 1, t
                 p(m,t) = p(m,t) + ( q(m,j) * cpt(m,t+1-j) )
 13           continue
 12        continue
 11     continue

        return
        end






        subroutine write_array ( a, m, n, uout )

        implicit none
        include 'MEMEQNS.BLOCK'
        double precision a(maxmlst,*)
        integer i, j, m, n, uout

        do 1 j = 1, n
           write(uout,2) ( a(i,j), i=1,m )
 1      continue
2       format(1x,5(1x,e9.4))
        return
        end




        subroutine evolveq ( kp, km, mi, nq, q, nmlst, nbins )

        implicit none
        include 'MEMEQNS.BLOCK'
        double precision kp(maxmlst,*), km(maxmlst,*)
        double precision q(maxmlst,maxqk)
        integer mi, nq, t, u, nbins, m, j, nmlst

        do 10 m = 1, nmlst
           q(m,1) = 0.d0
 10     continue
        q(mi,1) = 1.d0

        do 1 t = 2, nq+1
           u = max ( 1, t - nbins )
           do 2 m = 1, nmlst
              q(m,t) = 0.d0
              do 3 j = u, t-1
                 if ( m .gt. 1 ) then
                    q(m,t) = q(m,t) + ( q(m-1,j) * kp(m-1,t-j) )
                 end if
                 if ( m .lt. nmlst ) then
                    q(m,t) = q(m,t) + ( q(m+1,j) * km(m+1,t-j) )
                 end if
 3            continue
 2         continue
 1      continue

        return
        end







        subroutine build_histograms(fpt,nmlst,nfpt,b,kp,km,nbins)

        implicit none
        include 'MEMEQNS.BLOCK'
        double precision fpt(maxmlst,*)
        double precision kp(maxmlst,*), km(maxmlst,*)
        double precision b, jreal
        integer m, i, j, nmlst, nfpt(*), nbins



        if ( mode .ne. 'EQUIL' ) then
           write(*,*)
     1       'ERROR: Reflecting BCs required for QK integration'
           stop
        end if

        nbins = 0

        do 1 m = 1, nmlst
           do 20 j = 1, maxbins
              kp(m,j) = 0.d0
              km(m,j) = 0.d0
 20        continue
           do 2 i = 1, nfpt(m)
              j = int((abs(fpt(m,i)) / b) + 0.9999999)
              if ( j .gt. nbins ) nbins = j
              if ( j .gt. maxbins ) then
                 write(*,*) 'ERROR: maxbins not large enough'
                 write(*,*) 'caused by fpt(',m,',',i,')=',fpt(m,i)
                 stop
              end if
              if ( fpt(m,i) .gt. 0 ) then
                 kp(m,j) = kp(m,j) + 1
              else
                 km(m,j) = km(m,j) + 1
              end if
 2         continue
           do 3 j = 1,nbins
              kp(m,j) = kp(m,j) / nfpt(m)
              km(m,j) = km(m,j) / nfpt(m)
 3         continue
 1      continue

        write(*,*) 'nbins = ', nbins

        return
        end






        subroutine check_fpts ( frac_p, nmlst )

        implicit none
        include 'MEMEQNS.BLOCK'
        integer i, nmlst
        double precision frac_p(*)

        if ( frac_p(1) .ne. 1.d0 ) then
           write(*,*) 'ERROR: first milestone leaks: frac_p =',
     1         frac_p(1)
           stop
        end if

        if ( (mode .ne. 'MFPT1') .and.
     1       (frac_p(nmlst) .ne. 0.d0) ) then
           write(*,*) 'ERROR: last milestone leaks: frac_p =',
     1          frac_p(nmlst)
           stop
        end if

        do 1 i = 2, nmlst - 1
           if ( frac_p(i) .eq. 0 ) then
              write(*,100) i
 100          format('ERROR: no forward flux at milestone ', i7)
              stop
           end if
           if ( (mode .ne. 'MFPT1') .and.
     1          (frac_p(i) .eq. 1) ) then
              write(*,200) i
 200          format('ERROR: no backward flux at milestone ', i7)
              stop
             end if
 1      continue

        return
        end






        subroutine shallmfpt ( frac_p, mfpt, nmlst, mi, omfpt )

        implicit none
        include 'MEMEQNS.BLOCK'
        integer i, mi, nmlst
        double precision frac_m(maxmlst), x(maxmlst)
        double precision bvec(maxmlst)
        double precision frac_p(*), mfpt(*)
        double precision omfpt

        if ( mi .eq. nmlst ) then
           omfpt = 0.d0
           return
        end if

        do 1 i = 1, nmlst-1
           frac_m(i) = 1.d0 - frac_p(i)
           bvec ( i ) = 0.d0
 1      continue

        bvec ( mi ) = 1.d0
        call solve_mat ( frac_p, frac_m, x, bvec, nmlst - 1 )

        omfpt = 0.d0
        do 2 i = 1, nmlst-1
           omfpt = omfpt - mfpt(i) * x(i)
 2      continue

        return
        end







        subroutine rshallmfpt ( frac_p, mfpt, nmlst, mi, omfpt )

        implicit none
        include 'MEMEQNS.BLOCK'
        integer i, nmlst, mi
        double precision frac_p(*), mfpt(*), omfpt
        double precision rmfpt(maxmlst), rfrac_p(maxmlst)

        do 1 i = 1, nmlst
           rmfpt(i) = mfpt ( nmlst-i+1 )
           rfrac_p(i) = 1.d0 - frac_p ( nmlst-i+1 )
 1      continue

        call shallmfpt ( rfrac_p, rmfpt,
     1                   nmlst, nmlst - mi + 1, omfpt )

        return
        end






        subroutine print_equilibrium ( frac_p, mfpt, nmlst )

        implicit none
        include 'MEMEQNS.BLOCK'
        integer i, nmlst
        double precision totp, peq(maxmlst)
        double precision frac_p(*), mfpt(*)

c       forward markov rate of mlst i-1
        double precision rp

c       backward markov rate of mlst i
        double precision rm

        totp = 1.d0

        peq(1) = 1.d0

        do 1 i = 2, nmlst
           rp = frac_p(i-1) / mfpt(i-1)
           rm = (1.d0 - frac_p(i)) / mfpt(i)
           peq(i) = peq(i-1) * rp / rm
           totp = totp + peq(i)
 1      continue

c       normalize
        do 2 i = 1, nmlst
           peq(i) = peq(i) / totp
 2      continue

        write(*,*) 'Equilibrium distribution:'
        write(*,3) (peq(i),i=1,nmlst)
3       format(1x,5(1x,e9.5))
        return
        end






        subroutine read_fpts ( allfpts, fpt, nmlst, nfpt )

        implicit none
        include 'MEMEQNS.BLOCK'
        character*128 fptname
        character*(*) allfpts
        double precision fpt(maxmlst,maxtraj)
        integer nmlst, fm, fi, nm, m, i, nfpt(maxmlst)
        
        fm = 10
        call check_file_exists ( allfpts )
        open(UNIT=fm,FILE=allfpts(1:128),status='old')

        if ( mode .eq. 'MFPT1' ) then
           nm = nmlst - 1
        else
           nm = nmlst
        end if

        do 1 m = 1, nm
           read(fm,*) fptname
           call check_file_exists ( fptname )
           fi = fm + m
           open ( UNIT=fi, FILE=fptname, status='old' )
           write(*,100) fptname
 100       format ( 'Reading ', A50 )

           i = 1
 10        read ( fi, *, end=20 ) fpt(m,i)
           i = i + 1
           go to 10
 20        close ( fi )
           nfpt(m) = i - 1
 1      continue

        close ( fm )
        return
        end




        subroutine fpts2frac(fpt,nmlst,nfpt,mfpt,frac_p)

        implicit none
        include 'MEMEQNS.BLOCK'
        double precision fpt(maxmlst,maxtraj)
        double precision mfpt(*), frac_p(*)
        integer nmlst, nm, m, i, nfpt(maxmlst)
        
        if ( mode .eq. 'MFPT1' ) then
           nm = nmlst - 1
        else
           nm = nmlst
        end if

        do 1 m = 1, nm
           mfpt(m) = 0.d0
           frac_p(m) = 0.d0
           do 2 i = 1, nfpt(m)
              mfpt(m) = mfpt(m) + abs(fpt(m,i))
              if ( fpt(m,i) .gt. 0.d0 ) then
                frac_p(m) = frac_p(m) + 1
              end if
 2         continue
           frac_p(m) = frac_p(m) / nfpt(m)
           mfpt(m) = mfpt(m) / nfpt(m)
 1      continue

        return
        end





        subroutine solve_mat ( fp, fm, x, b, n )

        implicit none
        include 'MEMEQNS.BLOCK'
        integer i, n
        double precision fp(*), fm(*), x(*), b(*)
        double precision d(maxmlst)

        d(n) = -1.d0

c       make lower triangular
        do 1 i = n-1, 1, -1
           d(i) = -1.d0 - fp(i) * fm(i+1) / d(i+1)
           b(i) = b(i) - fm(i+1) * b(i+1) / d(i+1)
 1      continue

c       solve
        x(1) = b(1) / d(1)
        do 2 i = 2,n
           x(i) = (b(i) - fp(i-1) * x(i-1)) / d(i)
 2      continue

        return
        end





        subroutine query_user(allfpts,nmlst,init_mlst,nq,b)

        implicit none
        include 'MEMEQNS.BLOCK'
        character*(*) allfpts
        character c
        integer init_mlst, nmlst, nq
        double precision b

        nq = 0

        write(*,*) 'Number of milestones: '
        read(5,*) nmlst

        write(*,*) 'Filename of fpt file list: '
        read(5,*) allfpts

        write(*,*) 'Equilibrium run? (y/n) '
        read(5,*) c
        
        if ( c .eq. 'y' ) then
           mode = 'EQUIL'

           write(*,*) 'Perform QK integration? (y/n) '
           read(5,*) c
           if ( c .eq. 'y' ) then
              write(*,*) 'Histogram bin width (in FPT units): '
              read(5,*) b
              write(*,*)'Number of QK integration steps (in bin units):'
              read(5,*) nq
              if ( nq + 1 .gt. maxqk ) then
                 write(*,*) 'ERROR: maxqk not large enough'
                 stop
              end if
           end if
        else
           write(*,*) 'Two-way MFPT? (y/n) '
           read(5,*) c
           if ( c .eq. 'y' ) then
              mode = 'MFPT2'
           else
              mode = 'MFPT1'
           end if
        end if


        if ( (mode(1:4) .eq. 'MFPT') .or. (nq .gt. 0) ) then
 10        write(*,*) 'Initial milestone: '
           read(5,*) init_mlst
           if ( (init_mlst .lt. 1) .or.
     1          (init_mlst .gt. nmlst) ) then
              write(*,*) 'ERROR: invalid initial milestone'
              goto 10
           end if
        end if


        if ( mode .eq. 'EQUIL' ) then
           write(*,*) 'Equilibrium run'
        else if ( mode .eq. 'MFPT1' ) then
           write(*,*) 'MFPT run'
        else if ( mode .eq. 'MFPT2' ) then
           write(*,*) 'Two-way MFPT run'
        end if 
        
        if ( nq .ne. 0 ) then
           write(*,1000) nq, b
        end if
 1000   format('QK integration for ', i5, ' steps with ',
     1         'bin width of ', f10.5)

        write(*,100) nmlst
        write(*,110) allfpts

        if ( (mode(1:4) .eq. 'MFPT') .or. (nq .gt. 0) ) then
           write(*,120) init_mlst
        end if

 100    format(/, i7, ' milestones ')
 110    format('fpt filenames are listed in ', A30 )
 120    format('All probability begins in milestone ', i3)

        return
        end






        subroutine check_file_exists ( name )

        implicit none
        include 'MEMEQNS.BLOCK'
        character*(*) name
        logical exist


        inquire (FILE=name,EXIST=exist)

        if ( .not. exist ) then
           write(*,1) name
 1         format ('ERROR: File ', A30, ' does not exist')
           stop
        end if
        end
