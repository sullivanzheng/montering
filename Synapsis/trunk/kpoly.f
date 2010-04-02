c subroutine from Alex V. 
c calculate I(-1) and I(-2) Alexander Polynomial--put at ial(1) ial(2)
c
c example call kpoly(n,ial,ierr)

      subroutine kpoly(jr1,ial,ierr)
	implicit double precision (a-h,o-z)
      parameter (maxa=500,isi=203)
        common/matrix/x(maxa),y(maxa),z(maxa),bangle(maxa),rz(maxa)  
        common /ddx/ dx(maxa),dy(maxa),dz(maxa)
      dimension  ial(2)
      dimension  ix(isi),id(isi),ic1(isi),ic2(isi),cx(isi)
      dimension  da(isi,isi)
      data deps/1.d-7/
      jr2=jr1-2

c inventarization of intersections
c n4 - number of intersections
c cx(i) - x-value of i-th intersection
c ic1(i) - number of undergoing segment for i-th intersection
c ic2(i) - number of overgoing segment for i-th intersection

      n4=1
      t=-1.
      do 110 n1=1,jr2
      n11=n1+1
      jr3=jr1
      if(n1.lt.2)  jr3=jr3-1
      jr4=n1+2
      pdx1=dx(n1)
      pdy1=dy(n1)
      px1=x(n1)
      py1=y(n1)
      px11=x(n11)
      do 111 n2=jr4,jr3
      px2=x(n2)
      py2=y(n2)
      if(dabs(px2-px1).gt.2..or.dabs(py2-py1).gt.2.) goto 111 
      pdx2=dx(n2)
      pdy2=dy(n2)
      n21=n2+1
      px21=x(n21)
      d=pdx2*pdy1-pdx1*pdy2
      if(dabs(d).lt.deps) goto 111
      pdx12=pdx1*pdx2
      drx=(py2*pdx12-px2*pdx1*pdy2
     *    -py1*pdx12+px1*pdx2*pdy1)/d
       if(px1.le.px11) then
        if(drx.lt.px1.or.drx.ge.px11) goto 111
      else
        if(drx.le.px11.or.drx.gt.px1) goto 111
      endif
      if(px2.le.px21) then
        if(drx.lt.px2.or.drx.ge.px21) goto 111
      else
        if(drx.le.px21.or.drx.gt.px2) goto 111
      endif
      rx=drx
      if(n4.gt.isi) goto 1002
      rz1=(rx-x(n1))/dx(n1)*dz(n1)+z(n1)
      rz2=(rx-x(n2))/dx(n2)*dz(n2)+z(n2)
      cx(n4)=rx
      if(rz1-rz2.ge.0.) goto 401
      ic1(n4)=n1
      ic2(n4)=n2
      goto 402
 401  ic1(n4)=n2
      ic2(n4)=n1
 402  n4=n4+1
 111  continue
 110  continue
      n4=n4-1
      n4m=n4
c renumeration of intersections in the proper order

      if(n4.le.2) goto 31
      n41=n4-1
      do 403 n1=1,n41
      i1=ic1(n1)
      nv1=n1
      n11=n1+1
      do 404 n2=n11,n4
      if(ic1(n2)-i1) 405,406,404
 405  i1=ic1(n2)
      nv1=n2
      goto 404
 406  rl1=dabs(x(i1)-cx(nv1))
      rl2=dabs(x(i1)-cx(n2))
      if(rl1.lt.rl2) goto 404
      nv1=n2
 404  continue
      if(n1.ge.nv1) goto 403
      ir=ic1(n1)
      ic1(n1)=ic1(nv1)
      ic1(nv1)=ir
      ir=ic2(n1)
      ic2(n1)=ic2(nv1)
      ic2(nv1)=ir
      xr=cx(n1)
      cx(n1)=cx(nv1)
      cx(nv1)=xr
 403  continue

c Determination of the crossing type

 801  n4=n4m
      do 59 i=1,n4
      i1=ic1(i)
      i2=ic2(i)
      if((dx(i1)*dy(i2)-dx(i2)*dy(i1)).gt.0.) then
        id(i)=1
      else
        id(i)=-1
      endif
   59 continue


c determination of the number ix(i) of overgoing  part of line
c   for i-th intersection

      do 19 n1=1,n4
      nv=ic2(n1)
      n2=0
 72   n2=n2+1
      if(n2.ge.n4+1) goto 771
      if(ic1(n2)-nv) 72,73,71
 73   r1=dabs(x(nv)-cx(n1))
 74   r2=dabs(x(nv)-cx(n2))
      if(r1.lt.r2) goto 71
      goto 72
 771  ix(n1)=1
      goto 19
 71   ix(n1)=n2
 19   continue
      if(t.lt.-1.5) goto 802
c attempts for some detanglement

 230  mj=0
      i=0
 231  i=i+1
      if(i.ge.n4) goto 240
      if(ix(i).ne.i.and.ix(i).ne.i+1) goto 231
      mj=mj+1
      n4=n4-1
      if(n4.le.2) goto 31
      do 233 k=i,n4
      ix(k)=ix(k+1)
 233  continue
      do 234 k=1,n4
      if(ix(k).gt.i) ix(k)=ix(k)-1
 234  continue
      i=i-1
      goto 231

 240  if(ix(n4).ne.n4.and.ix(n4).ne.1) goto 232
      do 242 k=1,n4
      if(ix(k).ge.n4) ix(k)=1
 242  continue      
      mj=mj+1
      n4=n4-1
      if(n4.le.2) goto  31
 232  if(mj.gt.0) goto 230
      i=0
 235  i=i+1
      if(i.ge.n4-1) goto 244
      if(ix(i).ne.ix(i+1)) goto 235
      i1=i+1
      do 237 k=1,n4
      if(ix(k).eq.i1) goto 235
 237  continue
      mj=mj+1
      n4=n4-2
      if(n4.le.2) goto 31
      do 238 k=i,n4
      ix(k)=ix(k+2)
 238  continue
      do 239 k=1,n4
      if(ix(k).gt.i) ix(k)=ix(k)-2
 239  continue
      goto 235
 244  if(ix(n4-1).ne.ix(n4)) goto 236
      do 246 k=1,n4
      if(ix(k).eq.n4) goto 236
 246  continue
      mj=mj+1
      n4=n4-2
      if(n4.le.2) goto 31
      do 247 k=1,n4
      if(ix(k).eq.n4+1) ix(k)=1
 247  continue
 236  continue
      if(mj.gt.0) goto 230

c array assignments

 802  m4=n4-1
      do 500 ks=1,m4
      do 501 js=1,m4
      da(ks,js)=0
 501  continue
      js=ix(ks)
      if(js.eq.ks.or.js.eq.ks+1) then
        da(ks,ks)=-1.  
        da(ks,ks+1)=1.
      else
        if(id(ks).gt.0) then
          da(ks,ks)=1.  
          da(ks,ks+1)=-t
        else
          da(ks,ks)=-t  
          da(ks,ks+1)=1.
        endif
        da(ks,js)=t-1.
      endif
 500  continue
     
c calculation of determinant
 
      c=1.
      kmax=m4-1
      do 20 k=1,kmax
        if(dabs(da(k,k)).lt.deps) then
          jj=k+1
          c=-c
  50      if(jj.gt.m4) then 
            c=0.
            goto 90            
          endif
          if(dabs(da(jj,k)).gt.deps) then
            do 80 l=1,m4
            dr=da(k,l)
            da(k,l)=da(jj,l)
  80        da(jj,l)=dr
          else
            jj=jj+1
            goto 50
          endif
        endif

        jmin=k+1

        do 21 j=jmin,m4
          if(dabs(da(k,j)).gt.deps) then
            dr=da(k,j)/da(k,k)
            do 22 i=k,m4
              da(i,j)=da(i,j)-dr*da(i,k)
  22        continue
           endif
  21    continue
  20  continue

      do 30 i=1,m4
  30  c=c*da(i,i)
  90  c=dabs(c)

  790 if(c.gt.1.d7) then
        c=c/2.
        goto 790
      endif
      ikn=c+0.1
      if(t.eq.-2.) then
  789   if((ikn/2)*2.lt.ikn) goto 788
        ikn=ikn/2
        goto 789
  788   ial(2)=ikn
        return
      endif

      if(ikn.eq.1) then
        ial(1)=1
        ial(2)=1
        return
      elseif(ikn.eq.3) then
        ial(1)=3
        ial(2)=7
        return
      elseif(ikn.eq.5.or.ikn.eq.7) then
        ial(1)=ikn
        t=-2.
        goto 801
      elseif(ikn.gt.8) then
        ial(1)=ikn
        return
      endif
 
   31 ial(1)=1
      ial(2)=1
      return

 1002 ierr=2
      return
      end
 




