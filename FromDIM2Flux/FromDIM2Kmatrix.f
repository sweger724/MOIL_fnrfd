        program FromDIM2Kmatrix

C       When you run dim(run) the relevant transition data is written in
C       a line that starts with "dimdata", like this:
C       dimdata: 2 3  3 4  4 3   3240 580 990   38 720
C       Meaning of the eleven columns:
C       1,2: the first two integers are the ID of the anchors that were closer at the start of the 
C       forward trajectory
C       3,4: The second two pair of integers are the anchors that are closer at the end of the
C       trajectory
C       5,6: The third pair of integers are the anchors that were closer at the end of the 
C       backward trajectory (to check first hitting point)
C       7: the number of time steps for forward trajectory
C       8: last step the forward trajectory crossed the initial milestone
C       9: number of steps from last time the trajectory crossed the initial milestone
C       until trajectory reaches backward milestone  
C       10: index of the initial configuration for the trajectory
C       11: random number seed used
C
C       This program uses information in columns 1, 2, 3, 4 and 7 to construct the transition 
C       matrix and vector of milestone lifetime.
C
C       The data from dim(run) is input in file milesinput.dat
C       this file can contain the data from columns 1, 2, 3, 4 and 7 (four anchor style)
C       or data from columns 1, 2, 3, 4, 5, 6, and 7 (six anchor style)
C
C       The number of anchor style is specified in file milesinitial.dat
C       milesinitial.dat contains five rows.
C       In first row, four variables are readed:
C       1: (ncolu = 4 or 6) number of anchor style used in milesinput.dat
C       2: (nanc) number of anchors in your data (or maximum anchor index in your data)
C       3: (dt) time step in fs
C       4: (ndirec) it is 1, if going from anchor i to j defines a different
C           milestone than going from anchor j -> i
C          it is 0, if both anchor transitions correspond to the same milestone
C       From second to fifth row we input information needed to
C       construct transition matrix for flux calculation (with boundary
C       states provided by the user):
C       Second row: the number of anchors identified as "initial state"
C       Third row: the indexes of the corresponding initial state anchors
C       Fourth row: the number of anchors identified as "final state"
C       Fifth row: the indexes of the corresponding final state anchors
        double precision t,dt
        integer mxnc, mxtr,mxmi
C mxnc = maximum number of anchors
C mxtr = maximum number of matrix elements of K different than zero
C mxmi = maximum number of milestones
        parameter (mxnc =5000)
        parameter (mxmi = mxnc*10)
        parameter (mxtr = mxmi*10)
        integer n(mxnc,mxnc),nb(mxnc,mxnc),nmi(mxmi),nmj(mxmi)
        integer ntrann(mxtr),nfromi(mxmi)
        integer t1,t2,t3
        integer ntrani(mxtr),ntranj(mxtr),trant(mxtr)
        integer tfori(mxmi)
C the following are variable for histogram of Kij
C They are not tested yet
        double precision area(25),timeh,dth
        integer kh(25,500),irasave(25),trantime(25,1000)
        integer ntrantime(25)
C 
        integer ianch(1000),fanch(1000)
        integer nmbi(mxmi),nmbf(mxmi)
     
c       some defaults
        ndirec = 1
        ncolu = 6
        dt = 0.001d0
        nianch = 1
        nfanch = 1
c read how many columns, anchors and time step
        open(50,file='milesinitial.dat')
        read (50,*) ncolu, nanc, dt, ndirec
C open output files
        open(91,file='Kij_nonzero.dat',status='replace')
        open(92,file='milestone_lifetime.dat',status='replace')
        open(93,file='tij_nonzero.dat',status='replace')
        open(94,file='milestone_anchors.dat',status='replace')
        open(89,file='anchors_out.dat',status='replace')
        open(911,file='Kij_nonzero_for_q.dat',status='replace')
        open(922,file='milestone_lifetime_for_q.dat',status='replace')
C variation with time of a specific transition kernel
C not used currently
c        itra = 10

C zeroing matrix:
        do i = 1 , nanc
        do j = 1 , nanc
        n(i,j) = 0
        enddo
        enddo
C read milestoning data
       numdat = 0
        open(51,file='milesinput.dat')
       do while (.true.)
        if (ncolu.eq.4) then
C read anchors for initial and final interfaces, and time
        read(51,end=100,fmt=*) ip, jp, j1p, kp, t1
        else if (ncolu.eq.6) then
C read anchors for initial, final and fhp interfaces (this is the 
C current output of the more recent code)
        read(51,end=100,fmt=*) ip, jp, j1p, kp, j2p, k1p, t1, t2, t3
        else
        write(6,*) 'Something is wrong with milesinput.dat'
        endif
        numdat = numdat + 1
C need to consider a special case if
C milestones are not directional
C i.e. if milestones i->j and j->i are identical
        if (ndirec.eq.1) then
        i = ip
        j = jp
        j1 = j1p
        k = kp
        else if (ndirec.eq.0) then
C only consider one pair
        if (ip.lt.jp) then
        i = ip
        j = jp
        else
        i = jp
        j = ip
        endif
        if (j1p.lt.kp) then
        j1 = j1p
        k = kp
        else
        j1 = kp
        k = j1p
        endif
        endif
        
C first identify interfaces that are "on" in the data
        n(i,j) = 1
        nb(j1,k) = 1
        enddo
 100    continue
        close(51)

c create milestone indexes
        is = 1
        do i = 1 , nanc
        ncounti = 0
        do j = 1 , nanc
        if (i.ne.j.and.n(i,j).eq.1.and.nb(i,j).eq.1) then
        nmi(is) = i
        nmj(is) = j
        is = is + 1
        endif
        ncounti = ncounti + n(i,j) 
        enddo
C print anchors that are not contributing:
        if (ncounti.eq.0) write(89,*) i
        enddo

C number of active milestones
        isac = is - 1
        write(6,*) "number of milestones = ", isac

        do i = 1 , isac
        nfromi(i) = 0
        tfori(i) = 0
        write(94,*) i,nmi(i),nmj(i)
        enddo

        do i = 1, mxtr
        ntrann(i) = 0
        trant(i) = 0
        enddo

        iy = 0
        icktran = 0
        open(52,file='milesinput.dat')
        do i10 = 1, numdat
C read anchors for initial and final interfaces
        if (ncolu.eq.4) then
        read(52,*) ip, jp, j1p, kp, t1
        else if (ncolu.eq.6) then
        read(52,*) ip, jp, j1p, kp, j2p, k1p, t1, t2, t3
        endif
        if (ndirec.eq.1) then
        i = ip
        j = jp
        j1 = j1p
        k = kp
        else if (ndirec.eq.0) then
C only consider one pair
        if (ip.lt.jp) then
        i = ip
        j = jp
        else
        i = jp
        j = ip
        endif
        if (j1p.lt.kp) then
        j1 = j1p
        k = kp
        else
        j1 = kp
        k = j1p
        endif
        endif
        imiles = 0
        jmiles = 0
        do i2 = 1 ,isac
        if (nmi(i2).eq.i.and.nmj(i2).eq.j) imiles = i2
        if (nmi(i2).eq.j1.and.nmj(i2).eq.k) jmiles = i2
        enddo
        if (imiles.eq.0.or.jmiles.eq.0) go to 71

        if (iy.eq.0) then
        ntrani(1) = imiles
        ntranj(1) = jmiles
        ntrann(1) = ntrann(1) + 1
        trant(1) = trant(1) + t1
c        ntind = ntrann(1)
        iyy = iy + 1
        else

        do ix = 1 , iy
        i5 = ntrani(ix) 
        j5 = ntranj(ix) 
        if (imiles.eq.i5.and.jmiles.eq.j5) then
        ntrann(ix) = ntrann(ix) + 1
        trant(ix) = trant(ix) + t1
c        ntind = ntrann(ix)
        go to 71
        else
        if (ix.eq.iy) then
        ntrani(ix + 1) = imiles
        ntranj(ix + 1) = jmiles
        ntrann(ix + 1) = ntrann(ix + 1) + 1
        trant(ix + 1) = trant(ix +1) + t1
c        ntind = ntrann(ix + 1)
        iyy = iy + 1
        endif
        endif
        enddo
        endif
        iy = iyy

 71     continue

        enddo
        close(52)

C  Compute K and milestone and transition times
        
        do ix = 1 , iyy
        i = ntrani(ix)
        j = ntranj(ix)
        nfromi(i) = nfromi(i) + ntrann(ix)
        tfori(i) = tfori(i) + trant(ix)
        write(93,*) i,j,dble(trant(ix))/dble(ntrann(ix))
        enddo




        do ix = 1 , iyy
        i = ntrani(ix)
        j = ntranj(ix)
        write(91,*) i,j,dble(ntrann(ix))/dble(nfromi(i))
        enddo

        do i = 1, isac
        write(92,*) i,dble(tfori(i))*dt/dble(nfromi(i))
        enddo
        
C Write matrix with boundary conditions 
C First, need to read the number of anchors
C associated with the initial state
        read(50,*) nianch
C now read the initial anchor indexes
        read(50,*) (ianch(idat),idat=1,nianch)
C now the "final state"
        read(50,*) nfanch
C now read the final anchor indexes
        read(50,*) (fanch(idat),idat=1,nfanch)
        close(50)
C initialze some arrays:
        do i = 1 , isac
        nmbi(i) = 0
        nmbf(i) = 0
        enddo
        inumini = 0
C identify the final milestone:
        do i = 1 , isac
        do ifi = 1 , nfanch
        indfa = fanch(ifi)
        if (nmj(i).eq.indfa) nmbf(i)=1
        enddo
C and the initial milestones too:
        do ini = 1 , nianch
        india = ianch(ini)
        if (nmi(i).eq.india) then
        nmbi(i)=1
C I need an account of initial state milestones
        inumini = inumini + 1
        endif
        enddo
        enddo

C write in first line number of milestones
        write(911,*) isac
C Now the non-zero elements of the matrix
        do ix = 1 , iyy
        i = ntrani(ix)
        j = ntranj(ix)
        if (nmbf(i).eq.1) then
        else
        write(911,*) i,j,dble(ntrann(ix))/dble(nfromi(i))
        endif
        enddo
C write the transitions from final to initial
        do i = 1 , isac
        if (nmbf(i).eq.1) then
        do j = 1, isac
        if (nmbi(j).eq.1) then
        write(911,*) i,j,1.0d0/dble(inumini)
        endif
        enddo
        endif
        enddo

        do i = 1, isac
        if (nmbf(i).eq.0) then 
        write(922,*) i,dble(tfori(i))*dt/dble(nfromi(i))
        else if (nmbf(i).eq.1) then
        write(922,*) i,0.0d0
        endif
        enddo

        stop

cC The following is for histograms of K. No tested yet

c        do i = 1, icktran
c        do j = 1, 500
c        kh(i,j) = 0
c        enddo
c        area(i) = 0.0d0
c        enddo

c        dth = 2000.0
cC   compute K histogram for a specific milestone
c        do i = 1, icktran
c        iratmax = 1
c        numj = ntrantime(i)
c        if (numj.ne.0) then
c        do j = 1,numj
c        irat = int(float(trantime(i,j))/dth) + 1
c        Kh(i,irat) = Kh(i,irat) + 1
c        if (irat.gt.iratmax) iratmax = irat
c        enddo
c        endif
c        irasave(i) = iratmax
c        enddo

C normalize the histograms
c        do i = 1, icktran
c        numj = ntrantime(i)
c        if (numj.ne.0) then
c        iratmax = irasave(i)
c        do j= 1 , iratmax
c        timeh = float(j)*dth -dth*0.5
c        area(i) = area(i) + dth*Kh(i,j)
c        enddo
c        endif
c        enddo
cC print the histograms
c        do i = 1, icktran
c        numj = ntrantime(i)
c        if (numj.ne.0) then
c        iratmax = irasave(i)
c        do j= 1 , iratmax
c        iex = 10000*i
c        timeh = float(j)*dth -dth*0.5
c        write(iex,*) timeh, Kh(i,j)/area(i) 
c        enddo
c        endif
c        enddo
        
        end
