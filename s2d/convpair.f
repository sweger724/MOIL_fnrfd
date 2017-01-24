      subroutine convpair(d2v,ii,jj,kk)
c
      include 'COMMON/LENGTH.BLOCK'
      include 'COMMON/GETVEC.BLOCK'
c
      double precision d2v(3*maxpt2d,3*maxpt2d)
      integer ii,jj,kk
c
c.....xi xj ; xj xi
      d2v(ii,jj)     = d2v(ii,jj)     + offdiag(kk)
      d2v(jj,ii)     = d2v(jj,ii)     + offdiag(kk)
c.....yi xj ; xj yi ; yj xi ; xi yj
      d2v(ii+1,jj)   = d2v(ii+1,jj)   + offdiag(kk+1)
      d2v(jj,ii+1)   = d2v(jj,ii+1)   + offdiag(kk+1)
      d2v(ii,jj+1)   = d2v(ii,jj+1)   + offdiag(kk+1)
      d2v(jj+1,ii)   = d2v(jj+1,ii)   + offdiag(kk+1)
c.....zi xj ; xi zj ; xj zi ; zj xi
      d2v(ii+2,jj)   = d2v(ii+2,jj)   + offdiag(kk+2)
      d2v(ii,jj+2)   = d2v(ii,jj+2)   + offdiag(kk+2)
      d2v(jj+2,ii)   = d2v(jj+2,ii)   + offdiag(kk+2)
      d2v(jj,ii+2)   = d2v(jj,ii+2)   + offdiag(kk+2)
c.....yi yj ; yj yi
      d2v(ii+1,jj+1) = d2v(ii+1,jj+1) + offdiag(kk+3)
      d2v(jj+1,ii+1) = d2v(jj+1,ii+1) + offdiag(kk+3)
c.....zi yj ; yi zj ; zj yi ; yj zi
      d2v(ii+2,jj+1) = d2v(ii+2,jj+1) + offdiag(kk+4)
      d2v(ii+1,jj+2) = d2v(ii+1,jj+2) + offdiag(kk+4)
      d2v(jj+1,ii+2) = d2v(jj+1,ii+2) + offdiag(kk+4)
      d2v(jj+2,ii+1) = d2v(jj+2,ii+1) + offdiag(kk+4)
c.....zi zj ; zj zi
      d2v(ii+2,jj+2) = d2v(ii+2,jj+2) + offdiag(kk+5)
      d2v(jj+2,ii+2) = d2v(jj+2,ii+2) + offdiag(kk+5)
c
      return
      end
