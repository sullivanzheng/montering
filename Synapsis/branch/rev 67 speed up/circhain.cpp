#include "chain.h"
using namespace std;

int CircularChain::updateAllBangle()		
{		
	C[0].bangle = calAngle(C[maxnum], C[0]);		
	for (int i = 1; i < maxnum + 1; i++)		
	{		
		C[i].bangle = calAngle(C[i - 1], C[i]);		
	}		
    cout <<"============Bangle Updated=============="<<endl;		
	for (int i = 1; i < maxnum + 1; i++)		
        cout <<"|"<< C[i].bangle << endl;		
    cout<<"------------------------------"<<endl;		
    return 0;	
}		

double CircularChain::updateBangle(int i)		
{		
	if (i > maxnum || i < 0)		
	{		
		printf("bangle update (num): refered to non-existing segment %d", i);		
		getchar();		
		exit(EXIT_FAILURE);		
	}		

	if (i == 0)		
		C[0].bangle = calAngle(C[maxnum], C[0]);		
	else		
		C[i].bangle = calAngle(C[i - 1], C[i]);		
	return C[i].bangle;		
}

void CircularChain::driftProof()
{
	for (int i = 1; i < maxnum + 1; i++)
	{
		C[i].x -= C[0].x;
		C[i].y -= C[0].y;
		C[i].z -= C[0].z;
	}
	C[0].x = C[0].y = C[0].z = 0.0;
}

CircularChain::CircularChain(char const *filename,const int length)
:Chain(filename,true,length){
	this->E_t_updateWrithe_E_t();
	this->updateKPoly();
	this->snapshot("ini.txt");
}

CircularChain::CircularChain(int length):Chain(true,length){
	this->E_t_updateWrithe_E_t();
	this->updateKPoly();
	this->snapshot("ini.txt");
}

double CircularChain::calG_bSum()
{
	double temp = 0;
	for (int i = 0; i < maxnum + 1; i++)
		temp += G_b(C[i].bangle);
	return temp;
}
int CircularChain::crankshaft(int m, int n, double a)
{
	//Crankshaft for an angle
	//The segment is from X(m) to X(n).
	//M is the rotation matrix.
	if ((m<0||n<0)||(n>maxnum || m>maxnum))
	{
		cout<<"Illegal values of m and n ("<<m<<','<<n<<')';
		exit(EXIT_FAILURE);
	}
	//how many moves are accepted.
	if (fabs(C[0].x) > (totsegnum) / 2 || fabs(C[0].y) > (totsegnum) / 2 || fabs(C[0].z) > (totsegnum) / 2)
	{
		cout << "Step #" << stats.auto_moves() << " Drift proof." << endl;
		driftProof();
	}
	double M[3][3];
	SetRotM_crankshaft(M, m, n, a);
	int i;
	i = m;
	while (i != n)
	{
		double v[3];
		double temp[3];
		v[0] = C[i].dx;
		v[1] = C[i].dy;
		v[2] = C[i].dz;
		mat33mulvec3(M, v, temp);
		C[i].dx = temp[0];
		C[i].dy = temp[1];
		C[i].dz = temp[2];
		int j;
		j = (i == maxnum ? 0 : i + 1);
		C[j].x = C[i].x + C[i].dx;
		C[j].y = C[i].y + C[i].dy;
		C[j].z = C[i].z + C[i].dz;
		i = j;
	}
	updateBangle(m);	
	updateBangle(n);
	return 0;
}

double CircularChain::dE_reptation(int m, int n, int move){
	//this operation will make the reptation movement between node m  to n-1
	//in the clockwise fashion, ie, the arc: m m+1 m+2 ...[maxnum 0 1]...n-2
	//n-1. Therefore, n is not included.

	//Mind that the program will wrap back to head when tail of the chain is hit
	//(i.e. m>n).

	//Test if m and n are illegal.
	if (m==n || m<0 || n<0 || m>maxnum || n>maxnum){
		cout<<"[In CircularChain::reptation]"
			"Illegal values of m and n for reptation movement"<<endl;
		exit(EXIT_FAILURE);
	}

	//Test if the segnum is within range.
	int rep_segnum=wrap(n-m,totsegnum);
	if (rep_segnum<reptation_minlen || rep_segnum>reptation_maxlen){
		cout<<"[In CircularChain::reptation]"
			"segment involved is either larger than repation_maxlen"
			"or smaller than reptation_minlen: rep_segnum="<<rep_segnum<<endl;
		exit(EXIT_FAILURE);
	}

	//Test if move is within range. move>0 means the reptation is forward, otherwise backward.
	if (move==0 || abs(move)>int(rep_segnum/2)+1){
		cout<<"[In CircularChain::reptation]"
			"move is not in range, its absolute value should be within (-rep_segnum/2,rep_segnum/2)"
			"and !=0: move="<<rep_segnum<<endl;
		getchar();
		exit(EXIT_FAILURE);
	}
	//rigid protection detect
	int nend=wrap(n-1,totsegnum);//exclude the last vertex. ie. n is now the last vector now.
	static segment temp[maxa];
	int i,j;

	//Lazy method: moving move(move < 0) is equivalent to moving rep_segnum+move;
	move=wrap(move,rep_segnum);

	//copy n-move ~ n-1 vector to temp;
	i=wrap(n-move,totsegnum);j=0;
	while (i!=n){
		temp[j]=this->C[i];
		++i;++j;
		i=wrap(i,totsegnum); //wrapping around;
	}

	//move m~n-move-1 to m+move ~ n - 1
	i=wrap(n-move-1,totsegnum);
	j=wrap(n-1,totsegnum);
	while (i!=wrap(m-1,totsegnum)){
		this->C[j]=this->C[i];
		--i;--j;
		i=wrap(i,totsegnum);j=wrap(j,totsegnum);
	}
	
	//move temp back to m~m+move-1
	i=0;j=m;
	while(i<=move-1){
		this->C[j]=temp[i];
		++i;++j;
		j=wrap(j,totsegnum);
	}

	//Update X.
	i=m;
	while(i!=n){
		int pre=wrap(i-1,totsegnum);
		C[i].x=C[pre].x+C[pre].dx;
		C[i].y=C[pre].y+C[pre].dy;
		C[i].z=C[pre].z+C[pre].dz;
		i=wrap(i+1,totsegnum);
	}
	
	//Return the delta E: the change of energy.
	//The bandangles have only been altered at three locations. m,n and m+move
	double _a[3];
	_a[0]=this->C[m].bangle;
	_a[1]=this->C[wrap(m+move,totsegnum)].bangle;
	_a[2]=this->C[n].bangle;
	double a[3];
	a[0]=this->updateBangle(m);
	a[1]=this->updateBangle(wrap(m+move,totsegnum));
	a[2]=this->updateBangle(n);
	double dE=0;
	for (i=0;i<3;i++){
		dE+=G_b(a[i])-G_b(_a[i]);
	}
	return dE;
}

double CircularChain::dE_TrialCrankshaft(int m, int n, double a)
{
	if ((m<0||n<0)||(n>maxnum || m>maxnum))
	{
		cout<<"Illegal values of m and n ("<<m<<','<<n<<')';
		exit(EXIT_FAILURE);
	}
	//how many steps of moves have been performed.
	double M[3][3];//M is the rotation matrix.
	SetRotM_crankshaft(M, m, n, a);
	
	double dE;
	double newA1;
	double newA2;
	double oldA1;
	double oldA2;
	segment C1;
	segment C2;
	C1 = C[m];
	C2 = C[n == 0 ? maxnum : n - 1];
	double v[3];
	double temp[3];
	v[0] = C1.dx;v[1] = C1.dy;v[2] = C1.dz;
	mat33mulvec3(M, v, temp);
	C1.dx = temp[0];C1.dy = temp[1];C1.dz = temp[2];
	newA1 = calAngle(C[m == 0 ? maxnum : m - 1], C1);
	v[0] = C2.dx;v[1] = C2.dy;v[2] = C2.dz;
	mat33mulvec3(M, v, temp);
	C2.dx = temp[0];C2.dy = temp[1];C2.dz = temp[2];
	newA2 = calAngle(C2, C[n]);
	oldA1 = C[m].bangle;
	//cout<<" ROTATE 1>"<<oldA1<<"\t"<<C[m].bangle;
	oldA2 = C[n].bangle;
	//cout<<"2>"<<oldA2<<"\t"<<C[n].bangle<<endl;

	double const C = 2.4e-28 / (1.3806503e-23 * 300) / 3.4e-10;
	//C=3e-19 erg.cm=3e-28 J.m. Change this unit to kT.basepairlength
	dE = (+G_b(newA1) + G_b(newA2) - G_b(oldA1) - G_b(oldA2));

	/*
	C_=2*pi*pi*C;
	bt_angle=bt_angle/360;bpo_angle=bpo_angle/360;
	temp=4.0*pi*pi*140/2/N+C_/N*(d_Lk+bt_angle*0+bpo_angle)^2;
	*/
	return dE;
}

void CircularChain::snapshot(char *filename)
{
	ofstream fh (filename);
	char buf[300];
	int i;
	if (!fh.good())
	{
		cout << "file not writable" << endl;
		getchar();
		exit(EXIT_FAILURE);
	}
	sprintf(buf, "%6d %20s", maxnum + 2, "[Snapshot of a circle]");
	fh << buf << endl;
	i=0;
	sprintf(buf, "%6d%4s%12.6f%12.6f%12.6f%6d%6d%6d", 
        i + 1, "F", C[i].x, C[i].y, C[i].z, 1, (i == 0 ? maxnum + 1 : i), (i + 1 == maxnum + 1 ? 1 : i + 2));
	fh << buf << endl;
	for ( i = 1; i <= maxnum; i++)
	{	
		int mark=protect_list[i];
		sprintf(buf, "%6d%4s%12.6f%12.6f%12.6f%6d%6d%6d", 
			i + 1, mark?"Cl":"C", C[i].x, C[i].y, C[i].z, 1, (i == 0 ? maxnum + 1 : i), (i + 1 == maxnum + 1 ? 1 : i + 2));
		fh << buf << endl;
	}

	//The first bead is repeated for convenient processing of convertJmol2For.py.
	i=0;
	int mark=protect_list[i];
	sprintf(buf, "%6d%4s%12.6f%12.6f%12.6f%6d%6d%6d", 
			totsegnum + 1, mark?"Cl":"C", C[i].x, C[i].y, C[i].z, 1, (i == 0 ? maxnum + 1 : i), (i + 1 == maxnum + 1 ? 1 : i + 2));
		fh << buf << endl;


    fh << endl << "Detailed Info" << endl;
	for ( i = 0; i <= maxnum-1; i++)
	{
		sprintf(buf, "seg %d |dX(%15.13f %15.13f %15.13f)| %15.10f X[i+1]-X %15.10f %15.10f", 
            i+1, C[i].dx,C[i].dy,C[i].dz,
            modu(C[i].dx, C[i].dy, C[i].dz), 
            modu(C[i + 1].x - C[i].x, 
                 C[i + 1].y - C[i].y, 
                 C[i + 1].z - C[i].z), 
            C[i].bangle);
		fh << buf << endl;
	}
	i=maxnum;
    sprintf(buf, "seg %d |dX(%15.13f %15.13f %15.13f)| %15.10f X[i+1]-X %15s %15.10f", 
            i+1, C[i].dx,C[i].dy,C[i].dz,
            modu(C[i].dx, C[i].dy, C[i].dz), 
            "--------",
            C[i].bangle);
	fh << buf << endl;
	fh.close();
}

int CircularChain::IEV( int in,  int ik){
// excluded volume effects
// iev=0 corresponds to intersection
// iev=1 corresponds to no intersection
// This subroutine is translated from A. Vologodskii Monte FORTRAN 77 program.
// ER AND ER2 are volume exclusion DIAMETER.

	if ((in<0||ik<0)||(in>maxnum || ik>maxnum))
	{
		cout<<"Illegal values of in and ik "
			"at int CircularChain::IEV( int in,  int ik) ("<<in<<','<<ik<<')';
		exit(EXIT_FAILURE);
	}
	
	int temp;
	if (in>ik) {temp=in;in=ik;ik=temp;}

	int iev=0,idiam=1;
	const double eps=1e-7;
	double xij,yij,zij,a2,ddd;
	double b,b2,rna,rnb,ak,bk,t;
	double er,er2;				//PAY ATTENTION, ER HERE IS THE VOLUME EXCLUSION DIAMETER.
	er=this->VolEx_R*2.0;		//That is why we need to time this->VolEx_R by 2.0

	for (int i=in;i<=ik;i++){				// do 2 i=in,ik
		for (int j=0;j<=maxnum;j++){		// do 3 j=1,jr1
			if (j >= in && j <= ik) continue;//if(j.ge.in.and.j.le.ik) goto 3
			if (abs(i-j) <= VEcutoff) continue;//(if(iabs(i-j).le.lll) goto 3
/*			if (protect_list[i]==1 || protect_list[j]==1) continue; */
			xij=this->C[j].x-this->C[i].x;    //xij=x(j)-x(i)
			yij=this->C[j].y-this->C[i].y;//yij=y(j)-y(i)
			zij=this->C[j].z-this->C[i].z;//zij=z(j)-z(i)
			a2=xij*xij+yij*yij+zij*zij;
			ddd=a2;//ddd=a2
			if (a2 >= (2+er)*(2+er)) continue;// if(a2.ge.4.+4.*er+er2) goto 3
			b=C[i].dx*C[j].dx+C[i].dy*C[j].dy+C[i].dz*C[j].dz;// b=dx(i)*dx(j)+dy(i)*dy(j)+dz(i)*dz(j)
			b2=1.0-b*b; //b2=1.-b*b
			rna=xij*C[i].dx+yij*C[i].dy+zij*C[i].dz;//rna=xij*dx(i)+yij*dy(i)+zij*dz(i)
			rnb=xij*C[j].dx+yij*C[j].dy+zij*C[j].dz;//rnb=xij*dx(j)+yij*dy(j)+zij*dz(j)
			if(fabs(b2) < eps) goto g1; // if dXi//dXj   //if(abs(b2).lt.eps) goto 1
			ak=(rna-rnb*b)/b2;
			bk=(-rnb+rna*b)/b2;
			if(( ak<0 || ak>1 ) || (bk<0 || bk >1) ) goto g1;//if((ak.lt.0..or.ak.gt.1.).or.(bk.lt.0..or.bk.gt.1.)) goto 1
			ddd=a2+bk*bk+ak*ak+2.*bk*rnb-2.*ak*rna-2.*ak*bk*b;
			goto g4;
g1:         if(rna<0 || rna>1) goto g5;
			ddd=a2-rna*rna;
g5:         if(rnb>0 || rnb<-1.) goto g4;
			t=a2-rnb*rnb;
			if(t<ddd) ddd=t;
g4:         if(ddd<er*er) idiam=0;
		}//3 continue
	if(idiam==0) return iev;//iev=0 here;
	}    //2 continue
	iev=1;
	return iev;
}

double CircularChain::Slow_E_t_updateWrithe_E_t(){
	double temp[3];double dX[3];double modX;
	double _writhe=0;
	for (int s=0;s<=maxnum;s++){
		for (int t=0;t<=maxnum;t++){
			if (s==t) continue;
			Xprod_exp(C[t].dx,C[t].dy,C[t].dz,
					  C[s].dx,C[s].dy,C[s].dz,temp);
			dX[0]=C[t].x - C[s].x;
			dX[1]=C[t].y - C[s].y;
			dX[2]=C[t].z - C[s].z;
			_writhe += dot_product(temp,dX)/pow(modu2V(dX),3.0/2.0)
				/modu(C[t].dx, C[t].dy, C[t].dz)
				/modu(C[s].dx, C[s].dy, C[s].dz);
		}
	}
	this->writhe = _writhe/4/PI;
	this->E_t= 2 * PI * PI * C_t / (totsegnum * bpperseg) *
				(dLk - this->writhe)*(dLk - this->writhe);
	return this->E_t;
}

double CircularChain::E_t_updateWrithe_E_t(){
	this->writhe = this->_fastWr();
	this->E_t= 2 * PI * PI * C_t / (totsegnum * bpperseg) *
				(dLk - this->writhe)*(dLk - this->writhe);
	return this->E_t;
}

double CircularChain::_fastWr(){
	/* FORTRAN origin:
	c calculation of knot type and part of Wr

      call KNDWR(jr1,topl,jwr,ierr)

	c calculation of the second part of Wr

      call BWR(jr1,2,jr1,beta)
      wr=jwr+beta
	*/	

	int kndwr,ierr=0;
	kndwr=this->_kndwr(this->topl,ierr);
	if (ierr!=0){
		cout<<"[In CircularChain::_fastWr()] ierr!=0"<<endl;
		exit(EXIT_FAILURE);
	}

	double bwr;
	bwr=this->_bwr(1,maxnum);

	return double(kndwr)+bwr;
}

#include "f2c.h"

int  CircularChain::kpoly(int ial[2], int &ierr)//TODO problematic. Knot formation still.
{
/* k3.f -- translated by f2c (version 20060506).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/
	
/*  Subroutine KN calculates the values I(t) for t=-1 and t=-2. I(t) is related */
/*  to the Alexander Polynomial, D(t), by the next equation: */
/* 				I(t)=abs(D(t)/t**n) */
/*  The value of n is the largest for which I(t) is the integer. The subroutine */
/*  puts I(-1) and I(-2) into array ial. */
/*  The subroutine parameters: */
/* 	jr1 - the number of straight segments of the closed chain; */
/* 	ial - the array for I(-1) and I(-2) (output parameters); */
/* 	ierr = 2 if the number of crossings on the chain */
/*  projection is too big. Normally its value is equal to 0. */
/*  The coordinates of the chain vertices are in the arrays x, y, z which have to */
/*  be in the common block /x/. Note that values of x(jr1+1), y(jr1+1), and */
/*  z(jr1+1) should be defined and equal x(1), y(1), and z(1) correspondingly. Before */
/*  the subroutine call one should also define the vectors corresponding to the chain */
/*  segments: */
/*  dx(i)=x(i+1)-x(i) */
/*  dy(i)=y(i+1)-y(i) */
/*  dz(i)=z(i+1)-z(i). */
    /* Initialized data */

    static doublereal deps = 1e-7;

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;
	const int isi=270;
    /* Local variables */
    static doublereal c__, d__;
    static integer i__, j, k, l;
    static doublereal t, x[maxa], y[maxa], z__[maxa];
    static integer i1, i2, n1, n2, m4, n4;
    static real r1, r2;
    static doublereal da[isi*isi];	
    static integer id[maxcross], n11, n21, n41, jj;
    static doublereal dr;
    static integer mj, cx[maxcross];
    static doublereal dx[maxa], dy[maxa], dz[maxa];
    static integer ir, ks, ix[maxcross], js, nv;
    static real rx, xr;
    static integer ic1[maxcross], ic2[maxcross], n4i, jr1, jr2, jr3, jr4;
    static real rl1, rl2;
    static integer nv1;
    static doublereal px1, py1, px2, py2;
    static real rz1, rz2;
    static integer ikn;
    static doublereal px11, px21;
    static integer ipr;
    static doublereal drx, pdx1, pdy1, pdx2, pdy2;
    static integer jmin;
    static doublereal pdx12;
    static integer kmax;

	//Interfacing my program data to this routine.
	//F2C automatically converted the indexing difference.
	//However mind that indexing variables such as jr1 etc still needs to +1 since the
	//program now index like x[index-1]. Therefore indexing variables still need to convert
	//from C convention to Fortran.
	for (int i=0; i<=maxnum; i++){
		x[i]=C[i].x;y[i]=C[i].y;z__[i]=C[i].z;
		dx[i]=C[i].dx;dy[i]=C[i].dy;dz[i]=C[i].dz;
	}
	x[maxnum+1]=C[0].x;y[maxnum+1]=C[0].y;z__[maxnum+1]=C[0].z;
	jr1=maxnum+1;

	//-------------------------------------------//
/* Function Parameter */
/* ParaOver */
    jr2 = jr1 - 2;
/*  Search for intersections */
/*  n4 - total number of intersections */
/*  cx(i) - x-value of i-th intersection */
/*  ic1(i) - number of undergoing segment for i-th intersection */
/*  ic2(i) - number of overgoing segment for i-th intersection */
    ierr = 0;
    ipr = 0;
L491:
    n4 = 1;
    t = -1.;
    i__1 = jr2;
    for (n1 = 1; n1 <= i__1; ++n1) {
	n11 = n1 + 1;
	jr3 = jr1;
	if (n1 < 2) {
	    --jr3;
	}
	jr4 = n1 + 2;
	pdx1 = dx[n1 - 1];
	pdy1 = dy[n1 - 1];
	px1 = x[n1 - 1];
	py1 = y[n1 - 1];
	px11 = x[n11 - 1];
	i__2 = jr3;
	for (n2 = jr4; n2 <= i__2; ++n2) {
	    px2 = x[n2 - 1];
	    py2 = y[n2 - 1];
	    if ((d__1 = px2 - px1, abs(d__1)) > 2. || (d__2 = py2 - py1, abs(
		    d__2)) > 2.) {
		goto L111;
	    }
	    pdx2 = dx[n2 - 1];
	    pdy2 = dy[n2 - 1];
	    n21 = n2 + 1;
	    px21 = x[n21 - 1];
	    d__ = pdx2 * pdy1 - pdx1 * pdy2;
	    if (abs(d__) < deps) {
		goto L111;
	    }
	    pdx12 = pdx1 * pdx2;
	    drx = (py2 * pdx12 - px2 * pdx1 * pdy2 - py1 * pdx12 + px1 * pdx2 
		    * pdy1) / d__;
	    if (px1 <= px11) {
		if (drx < px1 || drx >= px11) {
		    goto L111;
		}
	    } else {
		if (drx <= px11 || drx > px1) {
		    goto L111;
		}
	    }
	    if (px2 <= px21) {
		if (drx < px2 || drx >= px21) {
		    goto L111;
		}
	    } else {
		if (drx <= px21 || drx > px2) {
		    goto L111;
		}
	    }
	    rx = drx;
	    if (n4 > maxcross) {
		goto L1002;
	    }
	    rz1 = (rx - x[n1 - 1]) / dx[n1 - 1] * dz[n1 - 1] + z__[n1 - 1];
	    rz2 = (rx - x[n2 - 1]) / dx[n2 - 1] * dz[n2 - 1] + z__[n2 - 1];
	    cx[n4 - 1] = rx;
	    if (rz1 - rz2 >= 0.f) {
		goto L401;
	    }
	    ic1[n4 - 1] = n1;
	    ic2[n4 - 1] = n2;
	    goto L402;
L401:
	    ic1[n4 - 1] = n2;
	    ic2[n4 - 1] = n1;
L402:
	    ++n4;
L111:
	    ;
	}
/* L110: */
    }
    --n4;
    n4i = n4;
/* renumeration of intersections in the proper order */
    if (n4 <= 2) {
	goto L31;
    }
    n41 = n4 - 1;
    i__1 = n41;
    for (n1 = 1; n1 <= i__1; ++n1) {
	i1 = ic1[n1 - 1];
	nv1 = n1;
	n11 = n1 + 1;
	i__2 = n4;
	for (n2 = n11; n2 <= i__2; ++n2) {
	    if ((i__3 = ic1[n2 - 1] - i1) < 0) {
		goto L405;
	    } else if (i__3 == 0) {
		goto L406;
	    } else {
		goto L404;
	    }
L405:
	    i1 = ic1[n2 - 1];
	    nv1 = n2;
	    goto L404;
L406:
	    rl1 = (d__1 = x[i1 - 1] - cx[nv1 - 1], abs(d__1));
	    rl2 = (d__1 = x[i1 - 1] - cx[n2 - 1], abs(d__1));
	    if (rl1 < rl2) {
		goto L404;
	    }
	    nv1 = n2;
L404:
	    ;
	}
	if (n1 >= nv1) {
	    goto L403;
	}
	ir = ic1[n1 - 1];
	ic1[n1 - 1] = ic1[nv1 - 1];
	ic1[nv1 - 1] = ir;
	ir = ic2[n1 - 1];
	ic2[n1 - 1] = ic2[nv1 - 1];
	ic2[nv1 - 1] = ir;
	xr = (real) cx[n1 - 1];
	cx[n1 - 1] = cx[nv1 - 1];
	cx[nv1 - 1] = xr;
L403:
	;
    }
/* Determination of the crossing type */
    i__1 = n4;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i1 = ic1[i__ - 1];
	i2 = ic2[i__ - 1];
	if (dx[i1 - 1] * dy[i2 - 1] - dx[i2 - 1] * dy[i1 - 1] > 0.f) {
	    id[i__ - 1] = 1;
	} else {
	    id[i__ - 1] = -1;
	}
/* L59: */
    }
/*   Determination of the number of overpassing generator, ix(i), */
/*   for i-th intersection */
    i__1 = n4;
    for (n1 = 1; n1 <= i__1; ++n1) {
	nv = ic2[n1 - 1];
	n2 = 0;
L72:
	++n2;
	if (n2 >= n4 + 1) {
	    goto L771;
	}
	if ((i__2 = ic1[n2 - 1] - nv) < 0) {
	    goto L72;
	} else if (i__2 == 0) {
	    goto L73;
	} else {
	    goto L71;
	}
L73:
	r1 = (d__1 = x[nv - 1] - cx[n1 - 1], abs(d__1));
/* L74: */
	r2 = (d__1 = x[nv - 1] - cx[n2 - 1], abs(d__1));
	if (r1 < r2) {
	    goto L71;
	}
	goto L72;
L771:
	ix[n1 - 1] = 1;
	goto L19;
L71:
	ix[n1 - 1] = n2;
L19:
	;
    }
/*  Attempts of detanglement */
L230:
    mj = 0;
    i__ = 0;
L235:
    ++i__;
    if (i__ >= n4 - 1) {
	goto L244;
    }
    if (ix[i__ - 1] != ix[i__]) {
	goto L235;
    }
    if (id[i__ - 1] * id[i__] == 1) {
	goto L235;
    }
    i1 = i__ + 1;
    i__1 = n4;
    for (k = 1; k <= i__1; ++k) {
	if (ix[k - 1] == i1) {
	    goto L235;
	}
/* L237: */
    }
    ++mj;
    n4 += -2;
    if (n4 <= 2) {
	goto L31;
    }
    i__1 = n4;
    for (k = i__; k <= i__1; ++k) {
	ix[k - 1] = ix[k + 1];
	id[k - 1] = id[k + 1];
/* L238: */
    }
    i__1 = n4;
    for (k = 1; k <= i__1; ++k) {
	if (ix[k - 1] > i__) {
	    ix[k - 1] += -2;
	}
/* L239: */
    }
    goto L235;
L244:
    if (ix[n4 - 2] != ix[n4 - 1]) {
	goto L236;
    }
    if (id[n4 - 2] * id[n4 - 1] == 1) {
	goto L236;
    }
    i__1 = n4;
    for (k = 1; k <= i__1; ++k) {
	if (ix[k - 1] == n4) {
	    goto L236;
	}
/* L246: */
    }
    ++mj;
    n4 += -2;
    if (n4 <= 2) {
	goto L31;
    }
    i__1 = n4;
    for (k = 1; k <= i__1; ++k) {
	if (ix[k - 1] == n4 + 1) {
	    ix[k - 1] = 1;
	}
/* L247: */
    }
L236:
    if (mj > 0) {
	goto L230;
    }
/*  Filling Alexander's matrix */
    m4 = n4 - 1;
    if (m4 > 269) {
	goto L1002;
    }
L801:
    i__1 = m4;
    for (ks = 1; ks <= i__1; ++ks) {
	i__2 = m4;
	for (js = 1; js <= i__2; ++js) {
	    da[ks + js * isi - (isi+1)] = 0.;
/* L501: */
	}
	js = ix[ks - 1];
	if (js == ks || js == ks + 1) {
	    da[ks + ks * isi - (isi+1)] = -1.f;
	    da[ks + (ks + 1) * isi - (isi+1)] = 1.f;
	} else {
	    if (id[ks - 1] > 0) {
		da[ks + ks * isi - (isi+1)] = 1.f;
		da[ks + (ks + 1) * isi - (isi+1)] = -t;
	    } else {
		da[ks + ks * isi - (isi+1)] = -t;
		da[ks + (ks + 1) * isi - (isi+1)] = 1.f;
	    }
	    da[ks + js * isi - (isi+1)] = t - 1.f;
	}
/* L500: */
    }
/*  Calculation of the determinant */
    c__ = 1.f;
    kmax = m4 - 1;
    i__1 = kmax;
    for (k = 1; k <= i__1; ++k) {
	if ((d__1 = da[k + k * isi - (isi+1)], abs(d__1)) < deps) {
	    jj = k + 1;
	    c__ = -c__;
L50:
	    if (jj > m4) {
		c__ = 0.f;
		goto L90;
	    }
	    if ((d__1 = da[jj + k * isi - (isi+1)], abs(d__1)) > deps) {
		i__2 = m4;
		for (l = 1; l <= i__2; ++l) {
		    dr = da[k + l * isi - (isi+1)];
		    da[k + l * isi - (isi+1)] = da[jj + l * isi - (isi+1)];
/* L80: */
		    da[jj + l * isi - (isi+1)] = dr;
		}
	    } else {
		++jj;
		goto L50;
	    }
	}
	jmin = k + 1;
	i__2 = m4;
	for (j = jmin; j <= i__2; ++j) {
	    if ((d__1 = da[k + j * isi - (isi+1)], abs(d__1)) > deps) {
		dr = da[k + j * isi - (isi+1)] / da[k + k * isi - (isi+1)];
		i__3 = m4;
		for (i__ = k; i__ <= i__3; ++i__) {
		    da[i__ + j * isi - (isi+1)] -= dr * da[i__ + k * isi - (isi+1)];
/* L22: */
		}
		}
/* L21: */
	}
/* L20: */
    }
    i__1 = m4;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L30: */
	c__ *= da[i__ + i__ * isi - (isi+1)];
    }
L90:
    c__ = abs(c__);
L790:
    if (c__ > 1e7) {
	c__ /= 2.;
	goto L790;
    }
    ikn = (integer) (c__ + .1);
    if (t == -2.) {
L789:
	if (ikn / 2 << 1 < ikn) {
	    goto L788;
	}
	ikn /= 2;
	goto L789;
L788:
	ial[1] = ikn;
	return 0;
    }
    if (ikn == 1) {
	ial[0] = 1;
	ial[1] = 1;
	return 0;
    } else if (ikn == 3) {
	ial[0] = 3;
	ial[1] = 7;
	return 0;
    } else {
	ial[0] = ikn;
	if (ikn >= 82) {
	    ial[1] = 49000;
	    return 0;
	} else {
	    t = -2.f;
	    goto L801;
	}
    }
L31:
    ial[0] = 1;
    ial[1] = 1;
    return 0;
L1002:
    if (ipr == 0) {
	i__1 = jr1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dr = x[i__ - 1];
	    x[i__ - 1] = z__[i__ - 1];
	    z__[i__ - 1] = dr;
	    dr = dx[i__ - 1];
	    dx[i__ - 1] = dz[i__ - 1];
	    dz[i__ - 1] = dr;
/* L345: */
	}
	ipr = 1;
	goto L491;
    } else if (ipr == 1) {
	i__1 = jr1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dr = y[i__ - 1];
	    y[i__ - 1] = z__[i__ - 1];
	    z__[i__ - 1] = dr;
	    dr = dy[i__ - 1];
	    dy[i__ - 1] = dz[i__ - 1];
	    dz[i__ - 1] = dr;
/* L346: */
	}
	ipr = 2;
	goto L491;
    } else {
	ierr = 2;
	return 0;
	}
    return 0;
} /* kn_ */


int CircularChain::updateKPoly(){
	return 0; //TODO temporarily disabled.
	int ierr=0,ial[2];
	this->kpoly(ial,ierr);
	if (ierr!=0){
		cout<<"[kpoly error]"<<endl;
		exit(EXIT_FAILURE);
	}
	this->AlexPoly[0]=ial[0];
	this->AlexPoly[1]=ial[1];
	return 0;
}

int CircularChain::checkConsistancy(){
	static const double eps=1e-5;
	for (int i=0;i<=maxnum;i++){
		if (dabs(C[wrap(i+1,totsegnum)].x-C[i].x-C[i].dx) > eps ||
			dabs(C[wrap(i+1,totsegnum)].y-C[i].y-C[i].dy) > eps ||
			dabs(C[wrap(i+1,totsegnum)].z-C[i].z-C[i].dz) > eps ){
				cout<<"[Inconsistant]"<<i
					<<"("
					<<C[wrap(i+1,totsegnum)].x-C[i].x-C[i].dx<<","
					<<C[wrap(i+1,totsegnum)].y-C[i].y-C[i].dy<<","
					<<C[wrap(i+1,totsegnum)].z-C[i].z-C[i].dz<<")"<<endl;
		}
	}
	return 0;
}

allrigid::allrigid(char *configfile,CircularChain * target){
	using namespace std;
	ifstream f(configfile);
	if (!f.good()){
		cout<<"Rigidbody configfile does not exist!"<<std::endl;
		exit(EXIT_FAILURE);
	}

	//Clear global protect list to 0;
	for (int i=0;i<=maxnum;i++)
		protect_list[i]=0;
	
	vector<int> r_protect;
	vector<vector<double>>r_ref_v;
	int hard_eof=0;

	while (!f.eof() && !hard_eof){
		string tempbuf;
		getline(f,tempbuf);
		stringstream sline(tempbuf);
		string initial;
		sline>>initial;
		if (initial[0]=='#') continue;

		string token(initial);
				//lock points (vertices), indices from 0.
		if (token=="rigidstart"){
			r_protect.clear();
			r_ref_v.clear();
		}
		else if (token==string("vec")){
			vector<double> a(3,0);
		    sline>>a[0]>>a[1]>>a[2];
			r_ref_v.push_back(a);
		}
		else if (token==string("$")){ //protect
			while (!sline.eof()){
				int temp;
				sline>>temp;
				r_protect.push_back(temp);
			}
		}
		else if (token==string("rigidend")){
			cls_rigid temp_rigid(target, r_protect, r_ref_v);
			this->R.push_back(temp_rigid);
		}
		else if (token==string("[endoffile]")){
			hard_eof=1;
		}

	}
	for (int i=0;i<this->R.size();i++)
		for (int j=0;j<this->R[i].protect.size();j++)
			this->protect.push_back(R[i].protect[j]);
	this->update_allrigid_and_E();
}

double allrigid::update_allrigid_and_E(){
	//This section can be customized for different rigid body set.
	if (R.size()==0) {
		this->E = 0;
		return 0.;
	}
	for (int i=0;i<this->R.size();i++){
		R[i].update_ref_v_xyz();
	}
	this->E = 0;

	//Axis Orientation. ref_v_xyz[0]
	this->AxisBeta=betaVec1Vec2(this->R[0].ref_v_xyz[0],this->R[1].ref_v_xyz[0]);

	//Curvature radius direction orientation. ref_v_xyz[1]
	this->RadiusBeta =betaVec1Vec2(this->R[0].ref_v_xyz[1],this->R[1].ref_v_xyz[1]);

	//Center point of each rigid body. ref_v_xyz[2]
	double t0[3],t1[3];
	t0[0]=R[0].ref_v_xyz[2][0] + R[0].target->C[R[0].protect[0]].x;
	t0[1]=R[0].ref_v_xyz[2][1] + R[0].target->C[R[0].protect[0]].y;
	t0[2]=R[0].ref_v_xyz[2][2] + R[0].target->C[R[0].protect[0]].z;

	t1[0]=R[1].ref_v_xyz[2][0] + R[1].target->C[R[1].protect[0]].x;
	t1[1]=R[1].ref_v_xyz[2][1] + R[1].target->C[R[1].protect[0]].y;
	t1[2]=R[1].ref_v_xyz[2][2] + R[1].target->C[R[1].protect[0]].z;
	
	this->r=modu(t1[0]-t0[0],t1[1]-t0[1],t1[2]-t0[2]);
	
//	A potential from Quan Du and A. Vologodskii.
/*	double sigma2=1;
	double r0=0.1,q=0.8;
	double A=17;

	this->E=
	A*exp(-((AxisBeta-PI)*(AxisBeta-PI)+(RadiusBeta-PI)*(RadiusBeta-PI))/2/sigma2)
	*(pow(r0/r,q*2)-2*pow(r0/r,q));
*/

/*
	this->E += r * 50;
	this->E += AxisBeta*(-10);
	this->E += RadiusBeta *(-10);
*/

// My latest design.
	static double r0=1.5,r01=2.,a0=20.0/180.*PI,R0=20.0/180.*PI;
	static double A=12,B=10,C=5; //A+B=17
	static double re=2.0,ree=0.05,Ar=1;
	this->E=-B*exp(
		-(AxisBeta-PI)*(AxisBeta-PI)/(a0*a0)/2-r*r/(r0*r0)/2
		-(RadiusBeta-PI)*(RadiusBeta-PI)/(R0*R0)/2)
//		-C*exp(-(RadiusBeta-PI)*(RadiusBeta-PI)/(R0*R0)/2-r*r/(r0*r0)/2)
		-A*exp(-r*r/(r01*r01)/2);
//		-A*(r>re?exp(-r*r/(r01*r01)/2):exp(-re*re/(r01*r01)/2))
/*		+Ar*(r>re?
				log(2*PI*re*re):
				r>=ree?log(2*PI*r*r):
					log(2*PI*ree*ree)); */
    return this->E;
}


cls_rigid::cls_rigid(CircularChain* r_target, std::vector<int> r_protect,
	std::vector < vector<double> > r_ref_v):
target(r_target),protect(r_protect),ref_v(r_ref_v){
	//When cls_rigid is constructed, it will automatically update the global
	//variable "protect_list", which is the copy of .protect member (vector)
	// of cls_rigid's container class "all_rigid".
	//This copy is made for better global access of protect list.
	double Mv[3][3];
	Mv[0][0]=target->C[protect[0]].dx;
	Mv[1][0]=target->C[protect[0]].dy;
	Mv[2][0]=target->C[protect[0]].dz;

	Mv[0][1]=target->C[protect[1]].dx;
	Mv[1][1]=target->C[protect[1]].dy;
	Mv[2][1]=target->C[protect[1]].dz;

	Mv[0][2]=target->C[protect[2]].dx;
	Mv[1][2]=target->C[protect[2]].dy;
	Mv[2][2]=target->C[protect[2]].dz;

	for (int i=0;i<this->ref_v.size();i++){
		vector<double> a(3,0);

		//input vector
		double b[3];
		b[0]=ref_v[i][0];b[1]=ref_v[i][1];b[2]=ref_v[i][2];

		//output vector
		double o[3];
		mat33mulvec3(Mv,b,o);
		for (int j=0;j<3;j++) { a[j]=o[j]; }
		ref_v_xyz.push_back(a);
	}

	std::cout<<"Rigid body constructed, ref_v_xyz as follows:"<<std::endl;
	for (int i=0;i<this->ref_v_xyz.size();i++){
		std::cout<<'['<<i<<']';
		std::cout<<ref_v_xyz[i][0]<<" ";
		std::cout<<ref_v_xyz[i][1]<<" ";
		std::cout<<ref_v_xyz[i][2]<<" "<<std::endl;
	}

	//Initialize or append the global protection list for rigid bodies.
	for (int i=0;i<this->protect.size();i++){
		protect_list[this->protect[i]]=1;
	}
//	protect_list.push_back(this->protect.back + 1);

	std::cout<<"Constructing rigid body is over. The current status of global protect_list is:"<<endl;
	for (int i=0;i<maxa;i++)
		if (protect_list[i])
			std::cout<<'['<<i<<']'<<protect_list[i]<<std::endl;
}

void cls_rigid::update_ref_v_xyz(){
	double Mv[3][3];
	Mv[0][0]=target->C[protect[0]].dx;
	Mv[1][0]=target->C[protect[0]].dy;
	Mv[2][0]=target->C[protect[0]].dz;

	Mv[0][1]=target->C[protect[1]].dx;
	Mv[1][1]=target->C[protect[1]].dy;
	Mv[2][1]=target->C[protect[1]].dz;

	Mv[0][2]=target->C[protect[2]].dx;
	Mv[1][2]=target->C[protect[2]].dy;
	Mv[2][2]=target->C[protect[2]].dz;

	for (int i=0;i<this->ref_v.size();i++){
		vector<double> a(3,0);

		//input vector
		double b[3];
		b[0]=ref_v[i][0];b[1]=ref_v[i][1];b[2]=ref_v[i][2];

		//output vector
		double o[3];
		mat33mulvec3(Mv,b,o);
		for (int j=0;j<3;j++) { a[j]=o[j]; }
		ref_v_xyz[i]=a;
	}
}
