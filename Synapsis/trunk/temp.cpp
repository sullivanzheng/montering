int * CircularChain::kpoly(int in_flag[2]){

//  Translate from Alexander Vologodskii's program.
//  Subroutine KN calculates the values I(t) for t=-1 and t=-2. I(t) is related
//  to the Alexander Polynomial, D(t), by the next equation:

//				I(t)=abs(D(t)/t**n)

//  The value of n is the largest for which I(t) is the integer. The subroutine
//  puts I(-1) and I(-2) into array ial.

//  The subroutine parameters: 
//	maxnum - the number of straight segments of the closed chain;
//	ial - the array for I(-1) and I(-2) (output parameters);
//	ierr = 2 if the number of crossings on the chain
//  projection is too big. Normally its value is equal to 0.

//  The coordinates of the chain vertices are in the arrays x, y, z which have to
//  be in the common block /x/. Note that values of x(maxnum+1), y(maxnum+1)
//  Therefore we have to update the following statement:
	C[maxnum+1].x=C[0].x;
	C[maxnum+1].y=C[0].y;
	C[maxnum+1].z=C[0].z;
//	and
//  z(maxnum+1) should be defined and equal x(1), y(1), and z(1) correspondingly. Before
//  the subroutine call one should also define the vectors corresponding to the chain
//  segments: 

//  dx(i)=x(i+1)-x(i)
//  dy(i)=y(i+1)-y(i)
//  dz(i)=z(i+1)-z(i).

      const int isl=500,isi=270;
	  const double eps=0.001;
     
      //common /x/ x(maxa),y(maxa),z(maxa)
      //common /d/ dx(maxa),dy(maxa),dz(maxa)
      int  ial[2];
      int  ix[isl],id[isl],ic1[isl],ic2[isl],cx[isl];
      double pdx1,pdy1,pdx2,pdy2,d,pdx12,drx,deps,px1,py1,px2,py2,px11,px21;
      double da[isi][isi],c, dr,t;
      const double deps = 1.d-7;
	  
      int jr2=maxnum-2;

//  Search for intersections
//  n4 - total number of intersections
//  cx(i) - x-value of i-th intersection
//  ic1(i) - number of undergoing segment for i-th intersection
//  ic2(i) - number of overgoing segment for i-th intersection
	  int ierr = 0, ipr = 0;
	  int jr3,n11,jr4,n2,n21;
	  double d,rx,rz1,rz2;
g491: int n4 = 0;
      double t=-1.d0;
	  for (int n1 = 0; nl<=jr2; nl++){
		  n11=n1+1;
		  jr3=maxnum;
		  if(n1<2)  jr3=jr3-1;
		  jr4=n1+2;
		  pdx1=C[n1].dx;
		  pdy1=C[n1].dy;
		  px1=C[n1].x;
		  py1=C[n1].y;
		  px11=C[n11].x;
		  for (int n2=jr4;n2<=jr3;n2++){
			  px2=C[n2].x;
			  py2=C[n2].y;
			  if (dabs(px2-px1)>2.d0 || dabs(py2-py1)>2.d0) continue; 
			  pdx2=C[n2].dx;
			  pdy2=C[n2].dy;
			  n21=n2+1;
			  px21=C[n21].x;
			  d=pdx2*pdy1-pdx1*pdy2;
			  if(dabs(d)<deps) continue; 
			  pdx12=pdx1*pdx2;
			  drx=(py2*pdx12-px2*pdx1*pdy2 - py1*pdx12+px1*pdx2*pdy1)/d;
			  if(px1<=px11) {
				if(drx<px1 || drx>=px11) continue; 
			  }
			  else{
				if(drx<=px11 || drx>px1) continue; 
			  }
			  if(px2<=px21) {
				if(drx<px2 || drx>=px21) continue; 
			  }
			  else
			  {
				if(drx<=px21 || drx>px2) continue; 
			  }
			  rx=drx;
			  if(n4>isl) goto g1002
				  rz1=(rx-C[n1].x)/C[n1].dx*C[n1].dz+C[n1].z;
			  rz2=(rx-C[n2].x)/C[n2].dx*C[n2].dz+C[n2].z;
			  cx[n4]=rx;
			  if(rz1-rz2>=0.) goto g401;
			  ic1[n4]=n1;
			  ic2[n4]=n2;
			  goto g402;
g401:         ic1[n4]=n2;
			  ic2[n4]=n1;
g402:         n4=n4+1;
		  }
	  }

      n4=n4-1;
      n4i=n4;

// renumeration of intersections in the proper order

      if(n4<=2) goto 31;
      int n41=n4-1, nv1,n11,il,ir;
	  double xr;
	  for (int n1=0; n1<=n41;n1++){
		  i1=ic1[n1];
		  nv1=n1;
		  n11=n1+1;
		  for (int n2=n11;n2<=n4;n2++){
			  if(ic1[n2]-i1 < 0) {goto g405;}
			  else if (ic1(n2)-i1 == 0) {goto g406;}
			  else {goto 404;}
	   g405:  i1=ic1[n2];
			  nv1=n2;
			  goto g404;
       g406:  rl1=abs(C[i1].x-cx[nv1]);
			  rl2=abs(C[i1].x-cx[n2]);
			  if(rl1<rl2) goto g404;
			  nv1=n2;
		  }
		  if(n1>=nv1) continue;
		  ir=ic1[n1];
		  ic1[n1]=ic1[nv1];
		  ic1[nv1]=ir;
		  ir=ic2[n1];
		  ic2[n1]=ic2[nv1];
		  ic2[nv1]=ir;
		  xr=cx[n1];
		  cx[n1]=cx[nv1];
		  cx[nv1]=xr;
	  }

// Determination of the crossing type

	  for(int i=0;i<=n4;i++){
		  i1=ic1[i];
		  i2=ic2[i];
		  if((C[il].dx*C[i2].dy-C[i2].dx*C[i1].dy)>0.){
			id[i]=1;
		  else
			id[i]=-1;
		  }
	  }


//   Determination of the number of overpassing generator, ix(i),
//   for i-th intersection

	  for (n1=0;n1<=n4;n1++){
		  nv=ic2[n1];
		  n2=0;
g72:      n2=n2+1;
		  if(n2>=n4+1) goto g771;
		  if(ic1[n2]-nv < 0) goto g72;
		  else if (ic1[n2]-nv == 0) goto g73;
	      else goto g71;
g73:      r1=abs(C[nv].x-cx[n1]);
g74:      r2=abs(C[nv].x-cx[n2]);
		  if(r1<r2) goto g71;
		  goto g72;
g771:     ix[n1]=1;
		  goto g19;
g71:      ix[n1]=n2;
	  }
      
       
//  Attempts of detanglement

g230: int mj=0;
      int i=0;
g235: i=i+1;
      if(i>=n4-1) goto g244;
      if(ix[i]!=ix[i+1]) goto g235;
      if(id[i]*id[i+1]==1) goto g235;
      i1=i+1;
	  for (int k=0;k<=n4;k++){
		if(ix[k]==i1) goto g235;
	  }
      mj=mj+1;
      n4=n4-2;
      if(n4<=2) goto g31;
	  for (int k=i; k<=n4;k++){
		  ix[k]=ix[k+2];
		  id[k]=id[k+2];
g238: }
	  for (int k=0; k<=n4;k++){
	      if(ix[k]>i) ix[k]=ix[k]-2;
g239:  }
      goto g235;
g244: if(ix[n4-1]!=ix[n4]) goto g236;
      if(id[n4-1]*id[n4]==1) goto g236;
      for (int k=0;k<=n4;k++)  if(ix[k]==n4) goto g236;
      mj=mj+1;
      n4=n4-2;
      if ( n4 <= 2) goto g31;
      for (int k=0;k<=n4;k++) if(ix[k]==n4+1) ix[k]=1;
g236: if (mj > 0) goto g230;

//  Filling Alexander's matrix

      m4=n4-1;
      if(m4>isi) goto g1002;

g801: for(int ks=0; k<=m4;k++){
		  for (js=0;j<=m4;j++)  da[ks][js]=0;
		  js=ix[ks];
		  if(js==ks || js==ks+1) {
			da[ks][ks]=-1.; 
			da[ks][ks+1]=1.;
		  }
		  else{
			  if(id[ks]>0){
				  da[ks][ks]=1. ; 
				  da[ks][ks+1]=-t;
			  }
			  else{
				  da[ks][ks]=-t ;
				  da[ks][ks+1]=1.;
			  }
			  da(ks,js)=t-1.;
		  }
	  }
      
//  Calculation of the determinant
 
      c=1.;
      kmax=m4-1;

	  for(int k=0; k<=kmax;k++){
		  if(dabs(da[k][k])<deps){
			  int jj=k+1;
			  c=-c;
	          if(jj>m4) { 
				c=0.;
				goto g90;            
			  }
			  if(dabs(da(jj,k))>deps) then
				do 80 l=1,m4
				dr=da(k,l)
				da(k,l)=da(jj,l)
	  80        da(jj,l)=dr
			  else
				jj=jj+1
				goto 50
			  endif
		  }

        jmin=k+1

        do 21 j=jmin,m4
          if(dabs(da(k,j))>deps) then
            dr=da(k,j)/da(k,k)
            do 22 i=k,m4
              da(i,j)=da(i,j)-dr*da(i,k)
  22        continue
           endif
  21    continue
	  }

      do 30 i=1,m4
  30  c=c*da(i,i)
g90:  c=dabs(c)


  790 if(c>1.d7) then
        c=c/2.d0
        goto 790
      endif
      ikn=c+1.d-1
      if(t==-2.d0) then
  789   if((ikn/2)*2<ikn) goto 788
        ikn=ikn/2
        goto 789
  788   ial[1]=ikn
        return
      endif

      if(ikn==1) then
        ial[0]=1
        ial[1]=1
       return
      elseif(ikn==3) then
        ial[0]=3
        ial[1]=7
        return
      else
        ial[1]=ikn
        if(ikn>=82) then
          ial[2]=49000
          return
        else  
          t=-2.
          goto 801
        endif
      endif

   31 ial[0]=1
      ial[1]=1
      return

g1002:continue
      if(ipr==0) then
        do 345 i=1,maxnum
          dr=x(i)
          x(i)=z(i)
          z(i)=dr
          dr=dx(i)
          dx(i)=dz(i)
          dz(i)=dr
  345   continue
        ipr=1
        goto 491
      elseif(ipr==1) then
        do 346 i=1,maxnum
          dr=y(i)
          y(i)=z(i)
          z(i)=dr
          dr=dy(i)
          dy(i)=dz(i)
          dz(i)=dr
  346   continue
        ipr=2
        goto 491
      else
        ierr=2
        return
      endif
      end

}
