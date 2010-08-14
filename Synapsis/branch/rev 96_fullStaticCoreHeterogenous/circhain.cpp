#include "chain.h"
using namespace std;

long CircularChain::updateAllBangle()		
{		
	C[0].bangle = calAngle(C[maxnum], C[0]);		
	for (long i = 1; i < maxnum + 1; i++)		
	{		
		C[i].bangle = calAngle(C[i - 1], C[i]);		
	}		
    cout <<"============Bangle Updated=============="<<endl;		
	for (long i = 1; i < maxnum + 1; i++)		
        cout <<"|"<< C[i].bangle << endl;		
    cout<<"------------------------------"<<endl;		
    return 0;	
}		

double CircularChain::updateBangle(long i)		
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
	for (long i = 1; i < maxnum + 1; i++)
	{
		C[i].x -= C[0].x;
		C[i].y -= C[0].y;
		C[i].z -= C[0].z;
	}
	C[0].x = C[0].y = C[0].z = 0.0;
}

CircularChain::CircularChain(char const *filename,const long totsegnum)
:Chain(filename,true,totsegnum){
	this->snapshot("ini.txt");
}

double CircularChain::G_b(long n){
	return getg(n) * C[n].bangle * C[n].bangle ;
}

double CircularChain::getg(long n){
	//return rigidity constant of the starting point of segment C[n]
	static double a = 141./bpperunit; //persistence length in # of unit length.
	double g; //Rigidity constant s.t. E = g * bangle * bangle;
	long n_pre = (n==0) ? maxnum : (n-1);
	g = 1/( (C[n_pre].l+C[n].l) / a );
	return g;
}

double CircularChain::calG_bSum()
{
	double temp = 0;
	for (long i = 0; i < maxnum + 1; i++)
		temp += G_b(i);
	return temp;
}


long CircularChain::crankshaft(long m, long n, double a)
{
	//Crankshaft for an angle
	//The segment is from X(m) to X(n).
	//M is the rotation matrix.
	if ((m<0||n<0)||(n>maxnum || m>maxnum))
	{
		cout<<"Illegal values of m and n ("<<m<<','<<n<<')';
		exit(EXIT_FAILURE);
	}

	//Anti-drift.
	if (fabs(C[0].x) > (totsegnum) / 2 || fabs(C[0].y) > (totsegnum) / 2 || fabs(C[0].z) > (totsegnum) / 2)
	{
		cout << "Step #" << stats.auto_moves() << " Drift proof." << endl;
		driftProof();
	}

	double M[3][3];
	SetRotM_crankshaft(M, m, n, a);
	long i;
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
		long j;
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

double CircularChain::dE_reptation(long m, long n, long move){
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
	long rep_segnum=wrap(n-m,totsegnum);
	if (rep_segnum<reptation_minlen || rep_segnum>reptation_maxlen){
		cout<<"[In CircularChain::reptation]"
			"segment involved is either larger than repation_maxlen"
			"or smaller than reptation_minlen: rep_segnum="<<rep_segnum<<endl;
		exit(EXIT_FAILURE);
	}

	//Test if move is within range. move>0 means the reptation is forward, otherwise backward.
	if (move==0 || abs(move)>long(rep_segnum/2)+1){
		cout<<"[In CircularChain::reptation]"
			"move is not in range, its absolute value should be within (-rep_segnum/2,rep_segnum/2)"
			"and !=0: move="<<rep_segnum<<endl;
		getchar();
		exit(EXIT_FAILURE);
	}
	//rigid protection detect
	long nend=wrap(n-1,totsegnum);//exclude the last vertex. ie. n is now the last vector now.
	static segment temp[maxa];
	long i,j;

	//Lazy method: moving move(move < 0) is equivalent to moving rep_segnum+move;
	move=wrap(move,rep_segnum);



	//Illustration:
	//Before movement:
	//          |<-----n-move-m------> |<-----move----> |
	//-->-->--> |-->-->-->-->-->-->--> |-->-->-->-->--> |-->-->-->-->-->-->-->-->
	//          m                      n-move       n-1 n

	//After movement:
	//          |<-----move----> |<-----n-move-m------> |
	//-->-->--> |-->-->-->-->--> |-->-->-->-->-->-->--> |-->-->-->-->-->-->-->-->
	//          m                m+move             n-1 n


	//store the initial energy at cleaving points
	double initialG_b = G_b(m) + G_b(n) + G_b(wrap(n-move,totsegnum));

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
		long pre=wrap(i-1,totsegnum);
		C[i].x=C[pre].x+C[pre].dx;
		C[i].y=C[pre].y+C[pre].dy;
		C[i].z=C[pre].z+C[pre].dz;
		i=wrap(i+1,totsegnum);
	}
	
	//Return the delta E: the change of energy.
	this->updateBangle(m);
	this->updateBangle(wrap(m+move,totsegnum));
	this->updateBangle(n);
	double dE = G_b(m) + G_b(n) + G_b(wrap(m+move,totsegnum)) - initialG_b;
	return dE;
}

double CircularChain::dE_TrialCrankshaft(long m, long n, double a)
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

	//double const C = 2.4e-28 / (1.3806503e-23 * 300) / 3.4e-10;
	//C=3e-19 erg.cm=3e-28 J.m. Change this unit to kT.basepairlength
	double gm = getg(m), gn = getg(n);
	dE = gm * newA1 * newA1 + gn * newA2 * newA2
	   - gm * oldA1 * oldA1 - gn * oldA2 * oldA2;

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
	long i;
	if (!fh.good())
	{
		cout << filename << " file not writable" << endl;
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
		long mark=protect_list[i];
		sprintf(buf, "%6d%4s%12.6f%12.6f%12.6f%6d%6d%6d", 
			i + 1, mark?"Cl":"C", C[i].x, C[i].y, C[i].z, 1, (i == 0 ? maxnum + 1 : i), (i + 1 == maxnum + 1 ? 1 : i + 2));
		fh << buf << endl;
	}

	//The first bead is repeated for convenient processing of convertJmol2For.py.
	i=0;
	long mark=protect_list[i];
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

long CircularChain::IEV_closeboundary( long in,  long ik){
// closeboundary means from check it IEV from seg in's starting point to ik's starting point (ik-1's endpoint)
// excluded volume effects
// iev=0 corresponds to intersection
// iev=1 corresponds to no intersection
// This subroutine is adapted from:
//       http://softsurfer.com/Archive/algorithm_0106/algorithm_0106.htm.
// ER is volume exclusion DIAMETER.
	long temp;
	if (in>ik) {temp=in;in=ik;ik=temp;}

	if ((in<0||ik<1)||(in>maxnum-1 || ik>maxnum))
	{
		cout<<"Illegal values of in and ik "
			"at long CircularChain::IEV( long in,  long ik) ("<<in<<','<<ik<<')';
		exit(EXIT_FAILURE);
	}
	
	const float eps = 2e-7;
    
	//PAY ATTENTION, ER HERE IS THE VOLUME EXCLUSION DIAMETER.
	float er = this->VolEx_R*2.0;		//That is why we need to time this->VolEx_R by 2.0

	for (long i=in;i<=ik-1;i++){		
		for (long j=0;j<=maxnum;j++){

			if (j >= in && j <= ik-1) continue;

      		long tempi=i<j?i:j,tempj=i<j?j:i;
			if (tempj-tempi <= VEcutoff || tempi+totsegnum-tempj <= VEcutoff) continue;
			
			float wx,wy,wz,w2;
			//wji=Xi-Xj
			wx = C[i].x-C[j].x;
			wy = C[i].y-C[j].y;
			wz = C[i].z-C[j].z;
			//if far far away
			w2 = wx*wx + wy*wy + wz*wz;
			if (w2 >= (C[i].l + C[j].l + er)*(C[i].l + C[j].l + er)) continue;

			float a,b,c,d,e,D;
			a = C[i].dx * C[i].dx + C[i].dy * C[i].dy + C[i].dz * C[i].dz;
			b = C[i].dx * C[j].dx + C[i].dy * C[j].dy + C[i].dz * C[j].dz;
			c = C[j].dx * C[j].dx + C[j].dy * C[j].dy + C[j].dz * C[j].dz;

			d = C[i].dx * wx + C[i].dy * wy + C[i].dz * wz;
			e = C[j].dx * wx + C[j].dy * wy + C[j].dz * wz;

			D = a * c - b * b;
			
			float    sc, sN, sD = D;      // sc = sN / sD, default sD = D >= 0
			float    tc, tN, tD = D;      // tc = tN / tD, default tD = D >= 0

			// compute the line parameters of the two closest points
			if (D < eps) { // the lines are almost parallel
				sN = 0.0;        // force using point P0 on segment S1
				sD = 1.0;        // to prevent possible division by 0.0 later
				tN = e;
				tD = c;
			}
			else {                // get the closest points on the infinite lines
				sN = (b*e - c*d);
				tN = (a*e - b*d);
				if (sN < 0.0) {       // sc < 0 => the s=0 edge is visible
					sN = 0.0;
					tN = e;
					tD = c;
				}
				else if (sN > sD) {  // sc > 1 => the s=1 edge is visible
					sN = sD;
					tN = e + b;
					tD = c;
				}
			}

			if (tN < 0.0) {           // tc < 0 => the t=0 edge is visible
				tN = 0.0;
				// recompute sc for this edge
				if (-d < 0.0)
					sN = 0.0;
				else if (-d > a)
					sN = sD;
				else {
					sN = -d;
					sD = a;
				}
			}
			else if (tN > tD) {      // tc > 1 => the t=1 edge is visible
				tN = tD;
				// recompute sc for this edge
				if ((-d + b) < 0.0)
					sN = 0;
				else if ((-d + b) > a)
					sN = sD;
				else {
					sN = (-d + b);
					sD = a;
				}
			}
			// finally do the division to get sc and tc
			sc = (abs(sN) < eps ? 0.0 : sN / sD);
			tc = (abs(tN) < eps ? 0.0 : tN / tD);

			//get the difference of the two closest points
			//Vector   dP = w + (sc * u) - (tc * v); 
			float dPx,dPy,dPz;
			dPx = wx + (sc * C[i].dx) - (tc * C[j].dx);
			dPy = wy + (sc * C[i].dy) - (tc * C[j].dy);
			dPz = wz + (sc * C[i].dz) - (tc * C[j].dz);
			
			//Collision detected, subroutine returns with 0;
			if (dPx * dPx + dPy * dPy + dPz * dPz < er * er ) return 0;
		}
	}
	//No collision detected
	return 1;
}

long CircularChain::IEV_Alex_closeboundary( long in,  long ik, double info[3]){
// closeboundary means from check it IEV from seg in's starting point to ik's starting point (ik-1's endpoint)
// excluded volume effects
// iev=0 corresponds to intersection
// iev=1 corresponds to no intersection
// This subroutine is translated from A. Vologodskii Monte FORTRAN 77 program.
// ER AND ER2 are volume exclusion DIAMETER.

	long temp;
	if (in>ik) {temp=in;in=ik;ik=temp;}

	if ((in<0||ik<1)||(in>maxnum-1 || ik>maxnum))
	{
		cout<<"Illegal values of in and ik "
			"at long CircularChain::IEV( long in,  long ik) ("<<in<<','<<ik<<')';
		exit(EXIT_FAILURE);
	}

        long iev=0,idiam=1;
        const double eps=1e-7;
        double xij,yij,zij,a2,ddd;
        double b,b2,rna,rnb,ak,bk,t;
        double er,er2;                          //PAY ATTENTION, ER HERE IS THE VOLUME EXCLUSION DIAMETER.
        er=this->VolEx_R*2.0;           //That is why we need to time this->VolEx_R by 2.0

        for (long i=in;i<=ik-1;i++){                              // do 2 i=in,ik
                for (long j=0;j<=maxnum;j++){           // do 3 j=1,jr1
                        if (j >= in && j <= ik-1) continue;//if(j.ge.in.and.j.le.ik) goto 3

      					long tempi=i<j?i:j,tempj=i<j?j:i;
						if (tempj-tempi <= VEcutoff || tempi+totsegnum-tempj <= VEcutoff) continue;

/*                      if (protect_list[i]==1 || protect_list[j]==1) continue; */
                        xij=this->C[j].x-this->C[i].x;    //xij=x(j)-x(i)
                        yij=this->C[j].y-this->C[i].y;//yij=y(j)-y(i)
                        zij=this->C[j].z-this->C[i].z;//zij=z(j)-z(i)
                        a2=modu2(xij,yij,zij);//a2=xij*xij+yij*yij+zij*zij
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
g1:					if(rna<0 || rna>1) goto g5;
                        ddd=a2-rna*rna;
g5:			        if(rnb>0 || rnb<-1.) goto g4;
                        t=a2-rnb*rnb;
                        if(t<ddd) ddd=t;
g4:					if(ddd<er*er) {
						idiam=0;
						info[0]=i;
						info[1]=j;
						return iev;
					}
				if (idiam == 0) 
					return iev;
				}//3 continue
		}    //2 continue
        iev=1;
        return iev;
}

long CircularChain::IEV_with_rigidbody_closeboundary( long in,  long ik, double info[3]){
// closeboundary means from check it IEV from seg in's starting point to ik's starting point (ik-1's endpoint)
// excluded volume effects
// iev=0 corresponds to intersection
// iev=1 corresponds to no intersection
// This subroutine is adapted from:
//       http://softsurfer.com/Archive/algorithm_0106/algorithm_0106.htm.
// ER is volume exclusion DIAMETER.

// info will be used to passback the first i,j segment pair that collides.

	long temp;
	if (in>ik) {temp=in;in=ik;ik=temp;}

	if ((in<0||ik<1)||(in>maxnum-1 || ik>maxnum))
	{
		cout<<"Illegal values of in and ik "
			"at long CircularChain::IEV( long in,  long ik) ("<<in<<','<<ik<<')';
		exit(EXIT_FAILURE);
	}


	const float eps = 2e-7;
    
	//PAY ATTENTION, ER HERE IS THE VOLUME EXCLUSION DIAMETER.
	float er = this->VolEx_R*2.0;		//That is why we need to time this->VolEx_R by 2.0

	for (long i=in;i<=ik-1;i++){		
		for (long j=0;j<=maxnum;j++){
			if (j >= in && j <= ik-1) continue;

			long tempi,tempj;
			if (i<j){
				tempi=i;tempj=j;
			}
			else{
				tempi=j;tempj=i;
			}
			//TODO: test this new cutoff algorithm.
			if (tempj-tempi <= VEcutoff || tempi+totsegnum-tempj <= VEcutoff) continue;

			if ((protect_list[i-VolEx_cutoff_rigidbody]==1 ||
				protect_list[i+VolEx_cutoff_rigidbody+1]==1) &&
				protect_list[j]==1) continue;	

			float wx,wy,wz,w2;
			//wji=Xi-Xj
			wx = C[i].x-C[j].x;
			wy = C[i].y-C[j].y;
			wz = C[i].z-C[j].z;
			//if far far away
			w2 = wx*wx + wy*wy + wz*wz;
			if (w2 >= (C[i].l + C[j].l + er)*(C[i].l + C[j].l + er)) continue;

			float a,b,c,d,e,D;
			a = C[i].dx * C[i].dx + C[i].dy * C[i].dy + C[i].dz * C[i].dz;
			b = C[i].dx * C[j].dx + C[i].dy * C[j].dy + C[i].dz * C[j].dz;
			c = C[j].dx * C[j].dx + C[j].dy * C[j].dy + C[j].dz * C[j].dz;

			d = C[i].dx * wx + C[i].dy * wy + C[i].dz * wz;
			e = C[j].dx * wx + C[j].dy * wy + C[j].dz * wz;

			D = a * c - b * b;
			
			float    sc, sN, sD = D;      // sc = sN / sD, default sD = D >= 0
			float    tc, tN, tD = D;      // tc = tN / tD, default tD = D >= 0

			// compute the line parameters of the two closest points
			if (D < eps) { // the lines are almost parallel
				sN = 0.0;        // force using point P0 on segment S1
				sD = 1.0;        // to prevent possible division by 0.0 later
				tN = e;
				tD = c;
			}
			else {                // get the closest points on the infinite lines
				sN = (b*e - c*d);
				tN = (a*e - b*d);
				if (sN < 0.0) {       // sc < 0 => the s=0 edge is visible
					sN = 0.0;
					tN = e;
					tD = c;
				}
				else if (sN > sD) {  // sc > 1 => the s=1 edge is visible
					sN = sD;
					tN = e + b;
					tD = c;
				}
			}

			if (tN < 0.0) {           // tc < 0 => the t=0 edge is visible
				tN = 0.0;
				// recompute sc for this edge
				if (-d < 0.0)
					sN = 0.0;
				else if (-d > a)
					sN = sD;
				else {
					sN = -d;
					sD = a;
				}
			}
			else if (tN > tD) {      // tc > 1 => the t=1 edge is visible
				tN = tD;
				// recompute sc for this edge
				if ((-d + b) < 0.0)
					sN = 0;
				else if ((-d + b) > a)
					sN = sD;
				else {
					sN = (-d + b);
					sD = a;
				}
			}
			// finally do the division to get sc and tc
			sc = (abs(sN) < eps ? 0.0 : sN / sD);
			tc = (abs(tN) < eps ? 0.0 : tN / tD);

			//get the difference of the two closest points
			//Vector   dP = w + (sc * u) - (tc * v); 
			float dPx,dPy,dPz;
			dPx = wx + (sc * C[i].dx) - (tc * C[j].dx);
			dPy = wy + (sc * C[i].dy) - (tc * C[j].dy);
			dPz = wz + (sc * C[i].dz) - (tc * C[j].dz);
			
			//Collision detected, subroutine returns with 0;
			if (dPx * dPx + dPy * dPy + dPz * dPz < er * er ) {
				info[0]=i;
				info[1]=j;
				return 0;
			}
		}
	}
	//No collision detected
	return 1;
}

double CircularChain::E_t_updateWrithe_E_t(){
	return 0; //RESUME;
	this->writhe = this->_fastWr_topl_update();
	this->E_t= 2 * PI * PI * C_t / (this->contour_length * bpperunit) *
				(dLk - this->writhe)*(dLk - this->writhe);
	return this->E_t;
}

double CircularChain::_fastWr_topl_update(){
	/* FORTRAN origin:
	c calculation of knot type and part of Wr

      call KNDWR(jr1,topl,jwr,ierr)

	c calculation of the second part of Wr

      call BWR(jr1,2,jr1,beta)
      wr=jwr+beta
	*/	

	long kndwr,ierr=0;
	kndwr=this->_kndwr_topl_update(this->topl,ierr);
	if (ierr!=0){
		cout<<"[In CircularChain::_fastWr()] ierr!=0"<<endl;
		exit(EXIT_FAILURE);
	}

	double bwr;
	bwr=this->_bwr(1,maxnum);

	return double(kndwr)+bwr;
}

long CircularChain::checkConsistancy(){
	//return 1 if inconsistent.
	static const double eps=1e-5;
	long flag=0;
	for (long i=0;i<=maxnum;i++){
		if (fabs(C[wrap(i+1,totsegnum)].x-C[i].x-C[i].dx) > eps ||
			fabs(C[wrap(i+1,totsegnum)].y-C[i].y-C[i].dy) > eps ||
			fabs(C[wrap(i+1,totsegnum)].z-C[i].z-C[i].dz) > eps ){
				cout<<"[Inconsistant]"<<i
					<<"("
					<<C[wrap(i+1,totsegnum)].x-C[i].x-C[i].dx<<","
					<<C[wrap(i+1,totsegnum)].y-C[i].y-C[i].dy<<","
					<<C[wrap(i+1,totsegnum)].z-C[i].z-C[i].dz<<")"<<endl;
				flag=1;
		}
	}
	return flag;
}

allrigid::allrigid(char *configfile,CircularChain * target){
	using namespace std;
	ifstream f(configfile);
	if (!f.good()){
		cout<<"Rigidbody configfile does not exist!"<<std::endl;
		exit(EXIT_FAILURE);
	}

	//Clear global protect list to 0;
	for (long i=0;i<=maxnum;i++)
		protect_list[i]=0;
	
	vector<long> r_protect;
	vector< vector<double> >r_ref_v;
	long hard_eof=0;

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
				long temp;
				sline>>temp;
				r_protect.push_back(temp);
			}
		}
		else if (token==string("rigidend")){
			cls_rigid temp_rigid(target, r_protect, r_ref_v);
			this->R.push_back(temp_rigid);
		}
		else if (token==string("sphere")){
			sphere s;
			sline>>s.rg0>>s.x0>>s.rg1>>s.x1>>s.d;
			this->spheres.push_back(s);
		}
		else if (token==string("[endoffile]")){
			hard_eof=1;
		}

	}
	for (long i=0;i<this->R.size();i++)
		for (long j=0;j<this->R[i].protect.size();j++)
			this->protect.push_back(R[i].protect[j]);
	this->update_allrigid_and_E();
}

int allrigid::IEV_spheres(long m, long n){
	//return 0-intersetion 1-no intersection
	vector<sphere>::iterator it;
	for ( it=this->spheres.begin(); it < this->spheres.end() ; it++){

		CircularChain *t0=this->R[it->rg0].target;
		CircularChain *t1=this->R[it->rg1].target;
		long pt0=this->R[0].protect[it->x0];
		long pt1=this->R[1].protect[it->x1];
		double cx=(t0->C[pt0].x + t1->C[pt1].x)/2.;
		double cy=(t0->C[pt0].y + t1->C[pt1].y)/2.;
		double cz=(t0->C[pt0].z + t1->C[pt1].z)/2.;
			
		long p=m;
		for (p=m; wrap(p+1,totsegnum)!=n ; p=wrap(p+1,totsegnum) ){
			if (protect_list[p-VolEx_cutoff_rigidbody]==1 ||
				protect_list[p+VolEx_cutoff_rigidbody-2]==1) continue;

			if (modu(t0->C[p].x - cx,t0->C[p].y - cy,t0->C[p].z - cz) < it->d/2.)
						 return 0;

			if (modu(t0->C[p].x + 0.5 * t0->C[p].dx - cx,
				     t0->C[p].y + 0.5 * t0->C[p].dy - cy,
					 t0->C[p].z + 0.5 * t0->C[p].dz - cz) < it->d/2.)
						 return 0;
		}
	}
	return 1;
}

double allrigid::update_allrigid_and_E(){
	this->E = 0; //TODO: disabled.
	return 0;
	//This section can be customized for different rigid body set.
	if (R.size()==0) {
		this->E = 0;
		return 0.;
	}
	for (long i=0;i<this->R.size();i++){
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

	//Res Site I distance. 
	//
	//TODO: this may result in access violation if R[i].protect[0]<3!
	static const long dev=6;
	this->r_siteI=
		modu(R[0].target->C[R[0].protect[0]-dev].x-R[1].target->C[R[1].protect[0]-dev].x,
			 R[0].target->C[R[0].protect[0]-dev].y-R[1].target->C[R[1].protect[0]-dev].y,
			 R[0].target->C[R[0].protect[0]-dev].z-R[1].target->C[R[1].protect[0]-dev].z);
	
//	A potential from Quan Du and A. Vologodskii.
/*	double sigma2=1;
	double r0=0.1,q=0.8;
	double A=17;

	this->E=
	A*exp(-((AxisBeta-PI)*(AxisBeta-PI)+(RadiusBeta-PI)*(RadiusBeta-PI))/2/sigma2)
	*(pow(r0/r,q*2)-2*pow(r0/r,q));
*/


	this->E += r*50 ;
	this->E += AxisBeta * (-15);
	this->E += RadiusBeta * (-15);
	this->E += r_siteI * 50;

	//if the Complex I Site I -> Complex II Site II has the same distance as
	//       Complex II Site I -> Complex I Site II.
	//This will identify the symmetry of the complete synaptic complex (Complex I-II)
	double symm;
	for (long ii=0;ii<=2;ii++){
		symm+= abs(
				modu(R[0].target->C[R[0].protect[ii]].x-R[1].target->C[R[1].protect[ii]-dev].x,
					 R[0].target->C[R[0].protect[ii]].y-R[1].target->C[R[1].protect[ii]-dev].y,
					 R[0].target->C[R[0].protect[ii]].z-R[1].target->C[R[1].protect[ii]-dev].z)
			   -modu(R[1].target->C[R[1].protect[ii]].x-R[1].target->C[R[0].protect[ii]-dev].x,
					 R[1].target->C[R[1].protect[ii]].y-R[1].target->C[R[0].protect[ii]-dev].y,
					 R[1].target->C[R[1].protect[ii]].z-R[1].target->C[R[0].protect[ii]-dev].z));
	}
	this->E += symm*5;


	/*double bend=60.0;
	this->E += pow(R[0].target->C[R[0].protect[0]-3].bangle - bend/180*PI, 2)*50;
	this->E += pow(R[1].target->C[R[1].protect[1]-3].bangle - bend/180*PI, 2)*50;*/


/*
// My latest design.
	//A+B=17 suggested.
	static const double B=10, a0=20.0/180.*PI, r0=1.5, R0=20.0/180.*PI;
	static const double A=12,r0_siteII=2.0;
	static const double A_siteII_att=80.0,r0_siteII_att=30.0;//Additional SiteI site II attachment force

	static const double A_siteI=80.0, r0_siteI=6.0;

	this->E=-B*exp(
		-(AxisBeta-PI)*(AxisBeta-PI)/(a0*a0)/2-r*r/(r0*r0)/2
		-(RadiusBeta-PI)*(RadiusBeta-PI)/(R0*R0)/2)

		-A*exp(-r*r/(r0_siteII*r0_siteII)/2)
		-A_siteII_att*exp(-r*r/(r0_siteII_att*r0_siteII_att)/2)

		-A_siteI*exp(-r_siteI*r_siteI/(r0_siteI*r0_siteI)/2);
*/
    return this->E;
}


cls_rigid::cls_rigid(CircularChain* r_target, std::vector<long> r_protect,
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

	for (long i=0;i<this->ref_v.size();i++){
		vector<double> a(3,0);

		//input vector
		double b[3];
		b[0]=ref_v[i][0];b[1]=ref_v[i][1];b[2]=ref_v[i][2];

		//output vector
		double o[3];
		mat33mulvec3(Mv,b,o);
		for (long j=0;j<3;j++) { a[j]=o[j]; }
		ref_v_xyz.push_back(a);
	}

	std::cout<<"Rigid body constructed, ref_v_xyz as follows:"<<std::endl;
	for (long i=0;i<this->ref_v_xyz.size();i++){
		std::cout<<'['<<i<<']';
		std::cout<<ref_v_xyz[i][0]<<" ";
		std::cout<<ref_v_xyz[i][1]<<" ";
		std::cout<<ref_v_xyz[i][2]<<" "<<std::endl;
	}

	//Initialize or append the global protection list for rigid bodies.
	for (long i=0;i<this->protect.size();i++){
		protect_list[this->protect[i]]=1;
	}
//	protect_list.push_back(this->protect.back + 1);

	std::cout<<"Constructing rigid body is over. The current status of global protect_list is:"<<endl;
	for (long i=0;i<maxa;i++)
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

	for (long i=0;i<this->ref_v.size();i++){
		vector<double> a(3,0);

		//input vector
		double b[3];
		b[0]=ref_v[i][0];b[1]=ref_v[i][1];b[2]=ref_v[i][2];

		//output vector
		double o[3];
		mat33mulvec3(Mv,b,o);
		for (long j=0;j<3;j++) { a[j]=o[j]; }
		ref_v_xyz[i]=a;
	}
}
