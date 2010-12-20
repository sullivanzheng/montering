#include "chain.h"
using namespace std;

long CircularChain::normalizeAllBangle()		
{		
	
/*	
    cout <<"============Circular Chain All Bangle normalized =============="<<endl;		
	for (long i = 1; i < maxnum + 1; i++);		
        //cout <<"|"<< C[i].bangle << endl;		
    cout<<"------------------------------"<<endl;		
*/
	//deform last 4 segment to make perfect circle:
	static segment C2[maxa];
	int p=maxnum-3,flag;
	double V[3]={C[0].x-C[p].x,C[0].y-C[p].y,C[0].z-C[p].z};
	double length=moduV(V);
	double Lnow,L[3]={0};
	for (int i=p;i<=maxnum;i++){
		L[0]+=C[i].dx; L[1]+=C[i].dy; L[2]+=C[i].dz;
	}
	Lnow=moduV(L);

	double diff=modu(C[0].x-C[p].x-L[0],
			 C[0].y-C[p].y-L[1],
			 C[0].z-C[p].z-L[2]);

	//cout <<"------Normalize the chain end------"<<endl;
	//cout <<"Before: X[maxnum-3]--->X[0]: Real:"<<length
	//	<<"-sum(dx):"<<Lnow<<" =diff:"<<diff
	//	<<endl;
	

	if (diff<1e-7){
	//	cout <<"Difference ignorable."<<endl;
		goto normang;
	}
	flag=this->_deformReptSegments_updateInternalBangles_Normalize(p,3,length,Lnow,C2);

	if (flag==-1) {
	//	cout<<"No solution."<<endl;
	}
	else{

		L[0]=L[1]=L[2]=0;
		for (int i=p;i<=maxnum;i++){
			L[0]+=C2[i].dx; L[1]+=C2[i].dy; L[2]+=C2[i].dz;
		}
		//Rotate V1new to V1
		double rt[3], M[3][3],rotAngle;
		//rotate from L to V.
		Xprod(L,V,rt);
		norm(rt);
		rotAngle=betaArray12(L,V);
		this->SetRotM_halfchain(M,rt,rotAngle);	

		for (long i=p;i<=maxnum;i++){
			double vi[3]={C2[i].dx,C2[i].dy,C2[i].dz},vo[3];
			mat33mulvec3(M,vi,vo);
			C2[i].dx=vo[0];	C2[i].dy=vo[1];	C2[i].dz=vo[2];
		}

		//copy C2 to C1.
		for (int i=p;i<=maxnum;i++){
			C[i]=C2[i];
		}

		//Print the result
		length=modu(C[0].x-C[p].x,C[0].y-C[p].y,C[0].z-C[p].z);
		L[0]=L[1]=L[2]=0;
		for (int i=p;i<=maxnum;i++){
			L[0]+=C[i].dx; L[1]+=C[i].dy; L[2]+=C[i].dz;
		}
		Lnow=moduV(L);
		cout <<"Normalize chain end."<<endl
			<<"After : X[maxnum-3]--->X[0]: Real:"<<length
		<<"-sum(dx):"<<Lnow<<" =diff:"<<
		modu(C[0].x-C[p].x-L[0],
			 C[0].y-C[p].y-L[1],
			 C[0].z-C[p].z-L[2])
		<<endl;
	}


normang:	C[0].bangle = calAngle(C[maxnum], C[0]);		
	for (long i = 1; i < maxnum + 1; i++)		
	{		
		C[i].bangle = calAngle(C[i - 1], C[i]);		
	}	

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

long CircularChain::checkBangleConsistency(){
	long flag=0;
	for (int i=0;i<=maxnum;i++){
		double b = calAngle(C[wrap(i-1,totsegnum)], C[i]);
		if (fabs(b-C[i].bangle)>1e-12){
			flag=1;
			cout<< "[Bangle inconsistency] #"<<i<<": "
				<<"b_calculated="<<b<<" C[i].bangle="<<C[i].bangle
				<<" diff="<<b-C[i].bangle<<endl;
		}
	}
	if (flag==1) return 1;
	else return 0;
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
	if (specify_rigidity_list[n]<0){
		//return rigidity constant of the starting point of segment C[n]
		static double a = 141./bpperunit; //persistence length in # of unit length.
		double g; //Rigidity constant s.t. E = g * bangle * bangle;
		long n_pre = (n==0) ? maxnum : (n-1);
		g = 1/( (C[n_pre].l+C[n].l) / a );
		return g;
	}
	else{
		return specify_rigidity_list[n];
	}
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

// Bangle and energy changes for repatation movement:
//                                                                                                           
//        * *                             * * *            *-bangle updated (bangle+=deltabangle)            
//-----|->->->-------------------------|->->->->----------------------------                                 
//     |    |                          |      |                                                              
//     m1   m1+dm1                     m2     m2+dm2                                                         
//
//     +  * * * +                        +  * * +          +-bangle needed to be recalculated.               
//-----|->->->->-------------------------|->->->-----------------------------                                
//     |      |                          |    |                                                              
//     m1     m1+dm2                     |    m2+dm2                                                         
//                                       m2+dm2-dm1                                                          


double CircularChain::_adjustBangle(long m, long dm, double dBangle, 
									segment const C2[maxa], segment Ctemp[maxa]){
	//adjustBangle for m+1 m+2 ... m+dm, and return the length of final vector.
	//Ctemp is the array to hold the conformation after adjustion.

	//copy C2 to Ctemp;
	for (int i=0;i<=dm;i++){
		long p=wrap(m+i,totsegnum);
		Ctemp[p]=C2[p];
	}

	long p = wrap(m+ 1, totsegnum);
	long n = wrap(m+dm, totsegnum);
	do{
		long s=wrap(p-1,totsegnum);
		double v1[3]={	Ctemp[s].dx,Ctemp[s].dy,Ctemp[s].dz}, 
			v2[3]={ Ctemp[p].dx,Ctemp[p].dy,Ctemp[p].dz};
		double rt[3];
		Xprod(v1,v2,rt);

		//if v1 and v2 are co-linear, then take x-axis as rotation axis.
		if (moduV(rt)<1e-12){
			rt[0]=1.0;
			cout<<"CircularChain::_adjustBangle-- Adjacent repation segments v1 and v2 almost co-linear."
				" then take x-axis as rotation axis.";
			cout<<endl;
		}

		norm(rt);
		double M[3][3];
		this->SetRotM_halfchain(M,rt,dBangle);

		//repation segments length check: shold always close to 3.
		if (C[s].l<rept_min_seglength || C[p].l<rept_min_seglength) {
			cout<<"Repation segment selection mistake! The segment involved in a repation movement is shorter than rept_min_seglength";
			cout<<endl;
		}


		for (long ss=p;ss!=wrap(n+1,totsegnum);ss=wrap(ss+1,totsegnum)){
			double vv[3]={Ctemp[ss].dx,Ctemp[ss].dy,Ctemp[ss].dz},
				vvo[3];
			mat33mulvec3(M,vv,vvo);
			Ctemp[ss].dx=vvo[0];
			Ctemp[ss].dy=vvo[1];
			Ctemp[ss].dz=vvo[2];
		}

		if (p==n) break;
		else p=wrap(p+1,totsegnum);
	}while (1);
	
	double X[3]={0};
	for(p=m;;){
		X[0]+=Ctemp[p].dx;
		X[1]+=Ctemp[p].dy;
		X[2]+=Ctemp[p].dz;
		
		if (p==n) break;
		else p=wrap(p+1,totsegnum);
	}

	return moduV(X);
}

double CircularChain::_bisectBangle(long m, long dm, 
       double a1, double f1, double a2, double f2, double length,
	   segment const C2[maxa], segment Ctemp[maxa], double eps){
    f1-=length; f2-=length;
	double a0=a1,f0=f1;
	long count=0;
	while ( fabs(f0) > eps){
		count++;
		a0=(a1+a2)/2;
		f0=this->_adjustBangle(m,dm,a0,C2,Ctemp)-length;
		if (f0*f1<0){
			a2=a0;
			f2=f0;
		}
		else{
			a1=a0;
			f1=f0;
		}
	}
	return a0;
	
}

int CircularChain::_deformReptSegments_updateInternalBangles(long m, long dm, 
									  double length, double Lnow, segment C2[maxa]){
	//Solve reptation adaptation problem for vector m, m+1 ... m+dm.
    //Segments in C2 will be deformed to appropriate conformation for the V1, V2 exchange.

    const int ERROR=-1;
	const double eps=1e-8;
	double minA=10,maxA=0;

	//copy C to C2, the playfield.
	for (long i=0;i<=dm;i++){
		long p=wrap(m+i,totsegnum);	
		if (C[p].bangle<minA) minA=C[p].bangle;
		if (C[p].bangle>maxA) maxA=C[p].bangle;
		C2[p]=C[p];
	}
	
	double A0,f;
	static const int tryL=10;
	if (fabs(Lnow-length)<eps){
		A0=0;
	}
	else if (Lnow<length){
		double An=-maxA*0.9, Ap=0.9*PI-maxA>0?0.9*PI-maxA:0;
		double Atry[tryL]={An*0.2, An*0.4, An*0.6, An*0.8, An, Ap*0.2, Ap* 0.4, Ap*0.6, Ap*0.8, Ap};//{An,Ap,An/2,Ap/2};
		for (int t=0; t<tryL; t++){
			A0=Atry[t];
			static segment Ctemp[maxa];
			f=this->_adjustBangle(m,dm,A0,C2,Ctemp);
			if (f>length) goto bisect;
		}
		return ERROR;
	}else if(Lnow>length){
		double An=-maxA*0.9, Ap=0.9*PI-maxA>0?0.9*PI-maxA:0;
		double Atry[tryL]={Ap*0.2, Ap* 0.4, Ap*0.6, Ap*0.8, Ap, An*0.2, An*0.4, An*0.6, An*0.8, An};//{Ap,An,Ap/2,An/2};
		for (int t=0; t<tryL; t++){
			A0=Atry[t];
			static segment Ctemp[maxa];
			f=this->_adjustBangle(m,dm,A0,C2,Ctemp);
			if (f<length) goto bisect;
		}
		return ERROR;
	}else{
		cout<<"CircularChain::_deformReptSegments_updateInternalBangles Problem";
		exit(EXIT_FAILURE);
	}
	
bisect:
	static segment Ctemp[maxa];
	if (A0!=0){
		A0=this->_bisectBangle(m,dm,A0,f,0,Lnow,length,C2,Ctemp,eps);

		//Copy Ctemp to C2
		for (long i=0;i<=dm;i++){
			long p=wrap(m+i,totsegnum);
			C2[p]=Ctemp[p];
		}
		
		//update bangle within V1 and V2;
		for (long i=1;i<=dm;i++){
			long p=wrap(m+i,totsegnum);
			C2[p].bangle = fabs(C2[p].bangle + A0);
		}
	}
	return 0;
}


int CircularChain::_deformReptSegments_updateInternalBangles_Normalize(long m, long dm, 
									  double length, double Lnow, segment C2[maxa]){

	//SPEICIALLY CUSTOMIZED FOR END STITCHING OF A CIRCULAR CHAIN.
	//Solve reptation adaptation problem for vector m, m+1 ... m+dm.
    //Segments in C2 will be deformed to appropriate conformation for the V1, V2 exchange.

    const int ERROR=-1;
	const double eps=1e-9;
	double minA=10,maxA=0;

	//copy C to C2, the playfield.
	for (long i=0;i<=dm;i++){
		long p=wrap(m+i,totsegnum);	
		if (C[p].bangle<minA) minA=C[p].bangle;
		if (C[p].bangle>maxA) maxA=C[p].bangle;
		C2[p]=C[p];
	}
	
	double A0,f;
	static const int tryL=2;
	if (fabs(Lnow-length)<eps){
		A0=0;
	}
	else if (Lnow<length){
		double Atry[tryL]={-0.01,0.01};
		for (int t=0; t<tryL; t++){
			A0=Atry[t];
			static segment Ctemp[maxa];
			f=this->_adjustBangle(m,dm,A0,C2,Ctemp);
			if (f>length) goto bisect;
		}
		return ERROR;
	}else if(Lnow>length){
		double Atry[tryL]={0.01,-0.01};
		for (int t=0; t<tryL; t++){
			A0=Atry[t];
			static segment Ctemp[maxa];
			f=this->_adjustBangle(m,dm,A0,C2,Ctemp);
			if (f<length) goto bisect;
		}
		return ERROR;
	}else{
		cout<<"CircularChain::_deformReptSegments_updateInternalBangles Problem";
		exit(EXIT_FAILURE);
	}
	
bisect:
	static segment Ctemp[maxa];
	if (A0!=0){
		A0=this->_bisectBangle(m,dm,A0,f,0,Lnow,length,C2,Ctemp,eps);

		//Copy Ctemp to C2
		for (long i=0;i<=dm;i++){
			long p=wrap(m+i,totsegnum);
			C2[p]=Ctemp[p];
		}
		
		//update bangle within V1 and V2;
		for (long i=1;i<=dm;i++){
			long p=wrap(m+i,totsegnum);
			C2[p].bangle = fabs(C2[p].bangle + A0);
		}
	}
	return 0;
}


double CircularChain::dE_reptation_3_4(long m1, long dm1, 
									   long m2, long dm2, int& rejection_sign){
	//exchange vectors of m1~m1+dm1 with m2~m2+dm2, assuming m2 is at the positive direction of m1.
	//In this case, dm1, dm2 = 2,3 or 3,2 respectively.
    //rejection_sign=0  successfully exchanged segments.
	//rejection_sign=1	unsuccessful deformation.
	//rejection_sign=2	shortcut of unsuccess: the contour of the shorter is shorter than the head-tail vector of the longer.
	static segment C2[maxa];

	rejection_sign=0;
	long m1lastV=wrap(m1+dm1,totsegnum),m2lastV=wrap(m2+dm2,totsegnum);
	long m1end=wrap(m1+dm1+1,totsegnum),m2end=wrap(m2+dm2+1,totsegnum);

	//cal v1 and v2 length

	double V1[3]={
		C[m1end].x - C[m1].x,
		C[m1end].y - C[m1].y,
		C[m1end].z - C[m1].z
	};

	double V2[3]={
		C[m2end].x - C[m2].x,
		C[m2end].y - C[m2].y,
		C[m2end].z - C[m2].z
	};
	
	double L1=moduV(V1), L2=moduV(V2);
	int flag;
	//Try to fit shorter segment to longer segment
	if (dm1>dm2){

		if (C[m2].l*(dm2+1) < L1 ){
			rejection_sign=2;
			return 99999;
		}

		flag=this->_deformReptSegments_updateInternalBangles(m2, dm2, L1, L2, C2);
		if (flag==-1) {
			rejection_sign=1; 
			//char buf[200]; 
			//sprintf(buf,"tttmove_%05d_43exchange_rej3to4before.txt",this->stats.auto_moves.getTotCounts()); 
			//this->snapshotseg(buf,C,m1,wrap(m1+dm1,totsegnum));
			//sprintf(buf,"tttmove_%05d_43exchange_rej3to4trailends.txt",this->stats.auto_moves.getTotCounts()); 
			//this->snapshotseg(buf,C,m2,wrap(m2+dm2,totsegnum));
			return 99999;
		}

		flag=this->_deformReptSegments_updateInternalBangles(m1, dm1, L2, L1, C2);
		if (flag==-1) {
			rejection_sign=1; 
			//char buf[200]; 
			//sprintf(buf,"tttmove_%05d_43exchange_rej4to3before.txt",this->stats.auto_moves.getTotCounts()); 
			//this->snapshotseg(buf,C,m2,wrap(m2+dm2,totsegnum));
			//sprintf(buf,"tttmove_%05d_43exchange_rej4to3trailends.txt",this->stats.auto_moves.getTotCounts()); 
			//this->snapshotseg(buf,C,m1,wrap(m1+dm1,totsegnum));
			return 99999;
		}
	}else{

		if (C[m1].l*(dm1+1) < L2 ){
			rejection_sign=2;
			return 99999;
		}

		flag=this->_deformReptSegments_updateInternalBangles(m1, dm1, L2, L1, C2);
		if (flag==-1) {
			rejection_sign=1; 
			//char buf[200]; 
			//sprintf(buf,"tttmove_%05d_34exchange_rej3to4before.txt",this->stats.auto_moves.getTotCounts()); 
			//this->snapshotseg(buf,C,m1,wrap(m1+dm1,totsegnum));
			//sprintf(buf,"tttmove_%05d_34exchange_rej3to4trailends.txt",this->stats.auto_moves.getTotCounts()); 
			//this->snapshotseg(buf,C,m2,wrap(m2+dm2,totsegnum));
			return 99999;
		}
		flag=this->_deformReptSegments_updateInternalBangles(m2, dm2, L1, L2, C2);
		if (flag==-1) {
			rejection_sign=1; 
			//char buf[200]; 
			//sprintf(buf,"tttmove_%05d_34exchange_rej4to3before.txt",this->stats.auto_moves.getTotCounts()); 
			//this->snapshotseg(buf,C,m2,wrap(m2+dm2,totsegnum));
			//sprintf(buf,"tttmove_%05d_34exchange_rej4to3trailends.txt",this->stats.auto_moves.getTotCounts()); 
			//this->snapshotseg(buf,C,m1,wrap(m1+dm1,totsegnum));
			return 99999;
		}
	}

	//Rotate V1 to V2, V2 to V1.
	double V1new[3]={0},V2new[3]={0};
	for (int i=0;i<=dm1;i++){
		long p=wrap(m1+i,totsegnum);
		V1new[0]+=C2[p].dx;
		V1new[1]+=C2[p].dy;
		V1new[2]+=C2[p].dz;
	}
	for (int i=0;i<=dm2;i++){
		long p=wrap(m2+i,totsegnum);
		V2new[0]+=C2[p].dx;
		V2new[1]+=C2[p].dy;
		V2new[2]+=C2[p].dz;
	}

		double rt[3], M[3][3],rotAngle;
		//rotate v1 to v2's position
		Xprod(V1new,V2,rt);
		norm(rt);
		rotAngle=betaArray12(V1new,V2);
		this->SetRotM_halfchain(M,rt,rotAngle);	

		for (long i=0;i<=dm1;i++){
			long p=wrap(m1+i,totsegnum);
			double vi[3]={C2[p].dx,C2[p].dy,C2[p].dz},vo[3];
			mat33mulvec3(M,vi,vo);
			C2[p].dx=vo[0];	C2[p].dy=vo[1];	C2[p].dz=vo[2];
		}

		//rotate v2 to v1's position
		Xprod(V2new,V1,rt);
		norm(rt);
		rotAngle=betaArray12(V2new,V1);
		this->SetRotM_halfchain(M,rt,rotAngle);	

		for (long i=0;i<=dm2;i++){
			long p=wrap(m2+i,totsegnum);
			double vi[3]={C2[p].dx,C2[p].dy,C2[p].dz},vo[3];
			mat33mulvec3(M,vi,vo);
			C2[p].dx=vo[0];	C2[p].dy=vo[1];	C2[p].dz=vo[2];
		}
/*
		V1new[0]=V1new[1]=V1new[2]=0;
		V2new[0]=V2new[1]=V2new[2]=0;
		for (int i=0;i<=dm1;i++){
			long p=wrap(m1+i,totsegnum);
			V1new[0]+=C2[p].dx;
			V1new[1]+=C2[p].dy;
			V1new[2]+=C2[p].dz;
		}
		for (int i=0;i<=dm2;i++){
			long p=wrap(m2+i,totsegnum);
			V2new[0]+=C2[p].dx;
			V2new[1]+=C2[p].dy;
			V2new[2]+=C2[p].dz;
		}
*/
	//Calcualte initial energy:
	//Bangle will be changed will be quite complicated.
// Bangle and energy changes for repatation movement:
//                                                                                                           
//        * *                             * * *            *-bangle updated (bangle+=deltabangle)            
//-----|->->->-------------------------|->->->->----------------------------                                 
//     |    |                          |      |                                                              
//     m1   m1+dm1                     m2     m2+dm2                                                         
//
//     +  * * * +                        +  * * +          +-bangle needed to be recalculated.               
//-----|->->->->-------------------------|->->->-----------------------------                                
//     |      |                          |    |                                                              
//     m1     m1+dm2                     |    m2+dm2                                               
//                                       m2+dm2-dm1                                                          

		double E0=0;
		//* angle of V1
		for (long i=1;i<=dm1;i++) E0+=this->G_b(wrap(m1+i,totsegnum));
		//* angle of V2
		for (long i=1;i<=dm2;i++) E0+=this->G_b(wrap(m2+i,totsegnum));
		if (m1==m2end) E0+=G_b(m1); //in case m1 and m2end meets
		else E0+=(G_b(m1)+G_b(m2end));
		if (m1end==m2) E0+=G_b(m2); //in case m1end and m2 meets
		else E0+=(G_b(m1end)+G_b(m2));

		
	//Incorporate C2 into C.

		//shift segments between V1 and V2.
		if (wrap(m1+dm1+1,totsegnum)!=m2) {//V1 and V2 are not connected.
			long shift=dm2-dm1;
			if (shift>0){//shift rightward
				long p=wrap(m2-1,totsegnum);
				do{
					C[wrap(p+shift,totsegnum)]=C[p];

					if (p==wrap(m1+dm1+1,totsegnum)) break;
					p=wrap(p-1,totsegnum);
				}while(1);
			}else if (shift<0){
				shift=-shift;
				long p=wrap(m1+dm1+1,totsegnum);
				do{
					C[wrap(p-shift,totsegnum)]=C[p];

					if (p==wrap(m2-1,totsegnum)) break;
					p=wrap(p+1,totsegnum);
				}while(1);
			}else{
				cout<<"CircularChain::dE_reptation_3_4: dm1 should not be equal to dm2"<<endl;
				exit(EXIT_FAILURE);
			}
		}
		//Incorporation of V2 and V1.
		for (long i=0;i<=dm2;i++){
			long p=wrap(m1+i,totsegnum),q=wrap(m2+i,totsegnum);
			C[p]=C2[q];
		}

		for (long i=0;i<=dm1;i++){
			long p=wrap(m2+dm2-dm1+i,totsegnum),q=wrap(m1+i,totsegnum);
			C[p]=C2[q];
		}

	//recalculate X
		for (long p=m1;;){ //Note that C[m1].x is trash information, even though C[m1].X is not moved.
			long last=wrap(p-1,totsegnum);
			C[p].x=C[last].x+C[last].dx;
			C[p].y=C[last].y+C[last].dy;
			C[p].z=C[last].z+C[last].dz;

			if (p==wrap(m2+dm2,totsegnum)) break;
			p = wrap(p+1,totsegnum);
		}

	//recalculate angles


// Chart of angle updates
//     +  * * * +                        +  * * +          +-bangle needed to be recalculated.               
//-----|->->->->-------------------------|->->->-----------------------------                                
//     |      |                          |    |                                                              
//     m1     m1+dm2                     |    m2+dm2                                                         
//                                       m2+dm2-dm1       
// recalculate: m1 m1+dm2+1 m2+dm2-dm1 m2+dm2+1
	long m1_=m1,m2_=wrap(m2+dm2-dm1,totsegnum);
	long dm1_=dm2,dm2_=dm1,m1end_=wrap(m1+dm2+1,totsegnum),m2end_=m2end;
	this->updateBangle(m1_);
	this->updateBangle(m1end_);
	this->updateBangle(m2_);
	this->updateBangle(m2end_);

	//recalculate energy and find the difference:
	double E1=0;
	//* angle of V1
	for (long i=1;i<=dm1_;i++) E1+=this->G_b(wrap(m1_+i,totsegnum));
	//* angle of V2
	for (long i=1;i<=dm2_;i++) E1+=this->G_b(wrap(m2_+i,totsegnum));
	if (m1_==m2end_) E1+=G_b(m1_); //in case m1_ and m2end_ meets
	else E1+=(G_b(m1_)+G_b(m2end_));
	if (m1end_==m2_) E1+=G_b(m2_); //in case m1end_ and m2_ meets
	else E1+=(G_b(m1end_)+G_b(m2_));
	return E1-E0;
}

double CircularChain::dE_reptation_simple(long m, long n, long move){
	//this operation will make the reptation movement from vector m  to n-1
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


double CircularChain::dE_treadmill(double direction){
	//snapshot("before.txt");//DEBUG
	vector<int> updateBangleList;
	static int const stepforward=0;
	if (direction>0){ //positive shift.
		segment C0=this->C[maxnum];
		int nextp;
		for (int p=maxnum;;){
			nextp=p-1;
			if (nextp<0) break;
			if (this->C[wrap(nextp-stepforward,totsegnum)].l>rept_min_seglength && 
				this->C[wrap(nextp+stepforward,totsegnum)].l>rept_min_seglength ){
				this->C[p]=this->C[nextp];
				p--;
			}else{
				//            nextp                  p   p+1
				//            |                      |   |
				//--->--->#-->*-->==>==>==>==>....==>@-->--->--->
				//--->--->--->#-->==>==>==>==>....==>*-->@-->--->
				//                ^                  ^   ^
				//                |                  |   |
				//                -----BangleChange===----
				while(this->C[wrap(nextp+stepforward,totsegnum)].l<rept_min_seglength || 
					this->C[wrap(nextp-stepforward,totsegnum)].l<rept_min_seglength) nextp--;
				this->C[p]=this->C[nextp];
				updateBangleList.push_back(nextp+1);
				updateBangleList.push_back(p);
				updateBangleList.push_back(p+1);
				p=nextp;
			}
		}
		this->C[0]=C0;

	}else{//negative shift.
		segment C0=this->C[0];
		int nextp;
		for (int p=0;;){
			nextp=p+1;
			if (nextp>maxnum) break;
			if (this->C[wrap(nextp-stepforward,totsegnum)].l>rept_min_seglength && 
				this->C[wrap(nextp+stepforward,totsegnum)].l>rept_min_seglength ){
				this->C[p]=this->C[nextp];
				p++;
			}else{
				//            p                      nextp
				//            |                      |   
				//--->--->--->#-->==>==>==>==>....==>*-->@-->--->
				//--->--->#-->*-->==>==>==>==>....==>@-->--->--->
				//            ^   ^                  ^    
				//            |   |                  |    
				//            ---------BangleChange---    
				while(this->C[wrap(nextp+stepforward,totsegnum)].l<rept_min_seglength || 
					this->C[wrap(nextp-stepforward,totsegnum)].l<rept_min_seglength) nextp++;
				this->C[p]=this->C[nextp];
				updateBangleList.push_back(nextp);
				updateBangleList.push_back(p);
				updateBangleList.push_back(p+1);
				p=nextp;
			}
		}
		this->C[maxnum]=C0;
	}

	double oldE=0,newE=0,oldEt,newEt;
//	oldEt=this->calG_bSum();//DEBUG
	for (vector<int>::iterator it=updateBangleList.begin();it!=updateBangleList.end();it++){
		oldE+=this->G_b(*it);
		this->updateBangle(*it);
		newE+=this->G_b(*it);
	}
//	newEt=this->calG_bSum();//DEBUG


	//update x,y,z
	for (int i=1;i<=maxnum;i++){
		C[i].x=C[i-1].x+C[i-1].dx;
		C[i].y=C[i-1].y+C[i-1].dy;
		C[i].z=C[i-1].z+C[i-1].dz;
	}
#if 0
	if (newE-oldE>10) {
		char buf[200];
		sprintf(buf,"after[%05d]%5.2f-(%6.2f,%6.2f).txt",this->stats.auto_moves(),newE-oldE,oldEt,newEt);
		snapshot(buf);
	}//DEBUG
#endif
	return newE-oldE;
}
int CircularChain::snapshotseg(char *filename, segment const * Ct, int start, int end){
	//This subroutine will display segments (vectors) from start to end.
	ofstream fh (filename);
	char buf[300];
	long i;

	if (!fh.good())
	{
		cout << filename << " file not writable" << endl;
		getchar();
		exit(EXIT_FAILURE);
	}

	long dm=0;//count number of segments (vectors). Total nubmer of atoms displayed is dm+1.
	for (long p=start;;) {
		dm++;
		if (p==end) break;
		p=wrap(p+1,totsegnum);
	}
	sprintf(buf, "%6d %20s", dm + 1, "[Snapshot of segments from i_start to i_end]");
	fh << buf << endl;

	//Draw starting atom
	i=start;
	sprintf(buf, "%6d%4s%12.6f%12.6f%12.6f%6d%6d", 
        1, "F", Ct[i].x, Ct[i].y, Ct[i].z, 1, 2);
	fh << buf << endl;

	//Draw inside atom
	for ( long p = 1; p <= dm - 1 ; p++)
	{	
		i=wrap(start+p,totsegnum);
		long mark=protect_list[i];
		sprintf(buf, "%6d%4s%12.6f%12.6f%12.6f%6d%6d%6d", 
			p + 1, mark?"Cl":"C", Ct[i].x, Ct[i].y, Ct[i].z, 1, p , p + 2);
		fh << buf << endl;
	}

	//Draw tailing atom
	i=wrap(end+1,totsegnum);
	long mark=protect_list[i];
	sprintf(buf, "%6d%4s%12.6f%12.6f%12.6f%6d%6d", 
			dm + 1, mark?"Cl":"C", Ct[i].x, Ct[i].y, Ct[i].z, 1, dm);
		fh << buf << endl;


    fh << endl << "Detailed Info" << endl;
	for ( long p = 0; p <= dm - 1; p++)
	{
		i=wrap(start+p,totsegnum);
		sprintf(buf, "seg %d (realseg) %d |dX(%15.13f %15.13f %15.13f)| %15.10f X[i+1]-X %15.10f %15.10f", 
            p+1, i , Ct[i].dx,Ct[i].dy,Ct[i].dz,
            modu(Ct[i].dx, Ct[i].dy, Ct[i].dz), 
            modu(Ct[i + 1].x - Ct[i].x, 
                 Ct[i + 1].y - Ct[i].y, 
                 Ct[i + 1].z - Ct[i].z), 
            Ct[i].bangle);
		fh << buf << endl;
	}
	fh.close();

	return 0;
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

						//Volume exclusion waiver.
						//Make sure there is no collision check between a segment in the vinicity of rigid body core.
						if (
							  (
								(  
								   protect_list[i-VolEx_cutoff_rigidbody]==1 
								   ||
								   protect_list[i+VolEx_cutoff_rigidbody+1]==1
								) 
								&&
								protect_list[j]==1
							  )
							  ||
							  (
								(  
								   protect_list[j-VolEx_cutoff_rigidbody]==1 
								   ||
								   protect_list[j+VolEx_cutoff_rigidbody+1]==1
								) 
								&&
								protect_list[i]==1
							  )
						)continue;	

      					long tempi=i<j?i:j,tempj=i<j?j:i;
						if (tempj-tempi <= VEcutoff || tempi+totsegnum-tempj <= VEcutoff) continue;

						if (protect_list[i]==1 && protect_list[j]==1) {
							cout <<"Illegal segment collision detection appeared."
								"Should not detect collision inside core."<<endl;
							cout <<"This error checking step is not for two half sites simulation"<<endl;
							continue;
						}
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

	if ((in<0||ik<0)||(in>maxnum || ik>maxnum))
	{
		cout<<"Illegal values of in and ik "
			"at long CircularChain::IEV( long in,  long ik) ("<<in<<','<<ik<<')';
		exit(EXIT_FAILURE);
	}


	const float eps = 2e-7;
    
	//PAY ATTENTION, ER HERE IS THE VOLUME EXCLUSION DIAMETER.
	static float er = this->VolEx_R*2.0;		//That is why we need to time this->VolEx_R by 2.0
	static float cutoff2=(this->max_seglength+er)*(this->max_seglength+er);
	float ernow;

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
			//Volume exclusion waiver.
			//Make sure there is no collision check between a segment in the vinicity of rigid body core.

			//Strict check
			/*
			if (
				  (
					(  
					   protect_list[i-VolEx_cutoff_rigidbody]==1 
					   ||
					   protect_list[i+VolEx_cutoff_rigidbody+1]==1
					) 
					&&
					protect_list[j]==1
				  )
				  ||
				  (
					(  
					   protect_list[j-VolEx_cutoff_rigidbody]==1 
					   ||
					   protect_list[j+VolEx_cutoff_rigidbody+1]==1
					) 
					&&
					protect_list[i]==1
				  )
			)continue;	
			*/

			float wx,wy,wz,w2;
			//wji=Xi-Xj
			wx = C[i].x-C[j].x;
			wy = C[i].y-C[j].y;
			wz = C[i].z-C[j].z;
			//if far far away
			w2 = wx*wx + wy*wy + wz*wz;

//This slow cutoff:	if (w2 >= (C[i].l + C[j].l + er)*(C[i].l + C[j].l + er)) continue;
// has been replaced by a faster verson:
			if (w2 > cutoff2) continue;
			
			ernow=er;
			//if segments include rigidbody segments, er will be reducted.
			if (C[i].l<rept_min_seglength || C[j].l<rept_min_seglength) ernow=2.0/5.0*er;

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
			if (dPx * dPx + dPy * dPy + dPz * dPz < ernow * ernow ) {
				info[0]=i;
				info[1]=j;
				return 0;
			}
		}
	}
	//No collision detected
	return 1;
}

long CircularChain::IEV_with_rigidbody_closeboundary_fullChain(double info[3]){
// Check the entire chain for any collision.
// iev=0 corresponds to intersection
// iev=1 corresponds to no intersection
// This subroutine is adapted from:
//       http://softsurfer.com/Archive/algorithm_0106/algorithm_0106.htm.
// ER is volume exclusion DIAMETER.

// info will be used to passback the first i,j segment pair that collides.
	
	const float eps = 2e-7;
    
	//PAY ATTENTION, ER HERE IS THE VOLUME EXCLUSION DIAMETER.
	static float er = this->VolEx_R*2.0;		//That is why we need to time this->VolEx_R by 2.0
	static float cutoff2=(this->max_seglength+er)*(this->max_seglength+er);
	float ernow;

	for (long i=0;i<=maxnum-1;i++){		
		for (long j=i+1;j<=maxnum;j++){

			if (j-i <= VEcutoff || i+totsegnum-j <= VEcutoff) continue;

			
			float wx,wy,wz,w2;
			//wji=Xi-Xj
			wx = C[i].x-C[j].x;
			wy = C[i].y-C[j].y;
			wz = C[i].z-C[j].z;
			//if far far away
			w2 = wx*wx + wy*wy + wz*wz;

//This slow cutoff:	if (w2 >= (C[i].l + C[j].l + er)*(C[i].l + C[j].l + er)) continue;
// has been replaced by a faster verson:
			if (w2 > cutoff2) continue;

			ernow=er;
			//if segments include rigidbody segments, er will be reducted.
			if (C[i].l<rept_min_seglength || C[j].l<rept_min_seglength) ernow=2.0/5.0*er;

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
			if (dPx * dPx + dPy * dPy + dPz * dPz < ernow * ernow ) {
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
	kndwr=this->_kndwr_topl_update_stable(this->topl,ierr);
	if (ierr!=0){
		cout<<"[In CircularChain::_fastWr()] ierr!=0"<<endl;
		exit(EXIT_FAILURE);
	}

	double bwr;
	bwr=this->_bwr(1,maxnum);

	return double(kndwr)+bwr;
}

long CircularChain::checkConsistency(double eps){
	
	//return 1 if inconsistent.
	long flag=0;
	for (long i=0;i<=maxnum;i++){
		if (fabs(C[wrap(i+1,totsegnum)].x-C[i].x-C[i].dx) > eps ||
			fabs(C[wrap(i+1,totsegnum)].y-C[i].y-C[i].dy) > eps ||
			fabs(C[wrap(i+1,totsegnum)].z-C[i].z-C[i].dz) > eps ){
				double m=modu(
					C[wrap(i+1,totsegnum)].x-C[i].x,
					C[wrap(i+1,totsegnum)].y-C[i].y,
					C[wrap(i+1,totsegnum)].z-C[i].z);
				double dm=modu(C[i].dx,C[i].dy,C[i].dz);
				double diff=modu(
					C[wrap(i+1,totsegnum)].x - C[i].x - C[i].dx,
					C[wrap(i+1,totsegnum)].y - C[i].y - C[i].dy,
					C[wrap(i+1,totsegnum)].z - C[i].z - C[i].dz);
				cout<<"[Inconsistant]"<<i<<endl
					<<"X-X modu("
					<<C[wrap(i+1,totsegnum)].x-C[i].x<<","
					<<C[wrap(i+1,totsegnum)].y-C[i].y<<","
					<<C[wrap(i+1,totsegnum)].z-C[i].z<<")="
					<<m<<endl
					<<"dX  modu("
					<<C[i].dx<<","
					<<C[i].dy<<","
					<<C[i].dz<<")="
					<<dm<<endl
					<<"diff    ("
					<<C[wrap(i+1,totsegnum)].x-C[i].x-C[i].dx<<","
					<<C[wrap(i+1,totsegnum)].y-C[i].y-C[i].dy<<","
					<<C[wrap(i+1,totsegnum)].z-C[i].z-C[i].dz<<")="
					<<diff<<endl;
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
	//Global specify_rigidity_list list to -1.0;
	for (long i=0;i<=maxnum;i++){
		protect_list[i]=0;
		specify_rigidity_list[i]=-1.0;
	}
	
	vector<long> r_protect;
	vector<long> r_ref_vec_basis; //The vectors of segments as new basis.
	vector< vector<double> >r_ref_v;//The rigid body vector expressed using the new basis.
	long r_anchor; 

	long hard_eof=0;
	long specify_rigidity_flag=0;

	while (!f.eof() && !hard_eof){
		string tempbuf;
		getline(f,tempbuf);
		stringstream sline(tempbuf);
		string initial;
		sline>>initial;
		if (initial[0]=='#') continue;
		
		

		string token(initial);

		//RIGID BODY DEFINITION
		//lock points (vertices), indices from 0.
		if (token=="rigidstart"){
			r_protect.clear();
			r_ref_v.clear();
			r_ref_vec_basis.clear();
		}
		else if (token==string("vec")){
			vector<double> a(3,0);
		    sline>>a[0]>>a[1]>>a[2];
			r_ref_v.push_back(a);
		}
		else if (token==string("$")){ //protect
			int counter=0;
			for(;;){
				long temp;
				sline>>temp;
				counter++;
				if (temp<0) break;
				if (counter>300){
					std::cout<<"allrigid::allrigid problem. ref_vec_basis token is followed by too many "
						"basis. Usually, only 3 is necessary"<<endl;
					exit(EXIT_FAILURE);
				}
				r_protect.push_back(temp);
			}
		}
		else if (token==string("ref_vec_basis")){
			int counter=0;
			for(;;){	
				long temp;
				sline >> temp;
				counter++;
				if (temp<0) break;
				if (counter>300){
					std::cout<<"allrigid::allrigid problem. ref_vec_basis token is followed by too many "
						"basis. Usually, only 3 is necessary"<<endl;
					exit(EXIT_FAILURE);
				}
				r_ref_vec_basis.push_back(temp);
			}
		}
		else if (token==string("anchor")){
			long temp;
			sline>>temp;
			r_anchor=temp;
		}
		
		else if (token==string("rigidend")){
			cls_rigid temp_rigid(target, r_protect, r_ref_v, r_ref_vec_basis, r_anchor);
			this->R.push_back(temp_rigid);
		}

		//Artificial rigidity specification.
		else if (token==string("specify_rigidity")){
			specify_rigidity_flag=1;
		}
		else if (token==string("end_specify_rigidity")){
			if (specify_rigidity_flag==1) {
				specify_rigidity_flag=0;
			}
			else{
				cout<<"Encounter end_specify_rigidity before specify_rigidity"<<endl;
		        exit(EXIT_FAILURE);
			}
		}
		else if (token==string(">")){
			if (specify_rigidity_flag==1){
				int a; double b;
				sline >>a>>b;
				specify_rigidity_list[a]=b;
			}
			else{
				cout<<"Encouter '>' while not in a specify_rigidity - end_specify_rigidity block"<<endl;
				exit(EXIT_FAILURE);
			}
		}

		else if (token==string("sphere")){
			sphere s;
			sline>>s.rg0>>s.x0>>s.rg1>>s.x1>>s.d;
			this->spheres.push_back(s);
		}

		else if (token==string("sphere_simple")){
			sphere_simple s;
			sline>>s.x>>s.y>>s.z>>s.r;
		}

		else if (token==string("TURN_ON_IEV_WHEN_Q_SMALLER_THAN")){
			sphere_simple s;
			sline>>this->IEV_EFFECTIVE_THRESHOLD;
		}

		else if (token==string("[endoffile]")){
			hard_eof=1;
		}

	}

	//Global rigid body initialization.
	for (long i=0;i<this->R.size();i++)
		for (long j=0;j<this->R[i].protect.size();j++)
			this->protect.push_back(R[i].protect[j]);
	this->update_allrigid_and_E();
	this->IEV_spheres(0,0);
}

int allrigid::IEV_spheres(long m, long n){
	//TODO: add simple judgement for sphere_simple.

	//For complicated spheres that locate the center by relative value of 
	//return 0-intersetion 1-no intersection
	if (this->IEV_spheres_effective==0) return 1;

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

		// Volume exclusion waiver.
		// Diagram of how the volume exclusion waiver works.
		// e.g.VolEx_cutoff_rigidbody=3
		// protect list index:
		//    -4-3-2-1 0 1 2 3 4 5 6 7 8  9 101112
		//    [********************************]  
		//->->->->->|->*>*>*>*>*>*>*>*>*>|->->->->
		//   |1 2 3 |-------rigid--------|1	2 3 |
		//->: free vector; *>: vector with protected (fixed) starting point.
		for (p=m; wrap(p+1,totsegnum)!=n ; p=wrap(p+1,totsegnum) ){
			if (t0->C[p].l < rept_min_seglength) continue; 

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
	this->unbiasedE = 0;
	//return 0;
	//This section can be customized for different rigid body set.
	if (R.size()==0) {
		this->E = 0;
		return 0.;
	}
	for (long i=0;i<this->R.size();i++){
		R[i].update_ref_v_xyz();
	}
	this->E = 0;

	//Axis Orientation. ref_v_xyz[3]=z
	this->AxisBeta=3.14159-betaVec1Vec2(this->R[0].ref_v_xyz[3],this->R[1].ref_v_xyz[3]); //z

	//Curvature radius direction orientation. ref_v_xyz[1]=x
	this->RadiusBeta =3.14159-betaVec1Vec2(this->R[0].ref_v_xyz[1],this->R[1].ref_v_xyz[1]);//x

	//Center point of each rigid body. ref_v_xyz[0]+anchor segment coordinate
	double t0[3],t1[3];
	t0[0]=R[0].ref_v_xyz[0][0] + R[0].target->C[R[0].protect[R[0].anchor] ].x;
	t0[1]=R[0].ref_v_xyz[0][1] + R[0].target->C[R[0].protect[R[0].anchor] ].y;
	t0[2]=R[0].ref_v_xyz[0][2] + R[0].target->C[R[0].protect[R[0].anchor] ].z;

	t1[0]=R[1].ref_v_xyz[0][0] + R[1].target->C[R[1].protect[R[1].anchor] ].x;
	t1[1]=R[1].ref_v_xyz[0][1] + R[1].target->C[R[1].protect[R[1].anchor] ].y;
	t1[2]=R[1].ref_v_xyz[0][2] + R[1].target->C[R[1].protect[R[1].anchor] ].z;
	
	this->r=modu(t1[0]-t0[0],t1[1]-t0[1],t1[2]-t0[2]);

	//Res Site I distance. 
	//
	//TODO: this may result in access violation if R[i].protect[0]<3!
	
	this->r_siteI=
		modu(R[0].target->C[R[0].protect[0]].x-R[1].target->C[R[1].protect[0]].x,
			 R[0].target->C[R[0].protect[0]].y-R[1].target->C[R[1].protect[0]].y,
			 R[0].target->C[R[0].protect[0]].z-R[1].target->C[R[1].protect[0]].z);

	double tempa,tempb;
	{
	double X,Y,Z;
	X=R[0].target->C[R[0].protect[0]].x-(R[0].ref_v_xyz[4][0]+t0[0]);
	Y=R[0].target->C[R[0].protect[0]].y-(R[0].ref_v_xyz[4][1]+t0[1]);
	Z=R[0].target->C[R[0].protect[0]].z-(R[0].ref_v_xyz[4][2]+t0[2]);
	tempa=modu(X,Y,Z);
	}
	{
	 double X,Y,Z;
	 X=R[1].target->C[R[1].protect[0]].x-(R[1].ref_v_xyz[4][0]+t1[0]);
	 Y=R[1].target->C[R[1].protect[0]].y-(R[1].ref_v_xyz[4][1]+t1[1]);
	 Z=R[1].target->C[R[1].protect[0]].z-(R[1].ref_v_xyz[4][2]+t1[2]);
	 tempb=modu(X,Y,Z);
	}
	this->r_siteI_deviation= tempa + tempb;
	//Calcuate if site I's are all aligned to the -y direction.
	/*double yangle1=0,yangle2=0;
	{
		double Yx=R[0].ref_v_xyz[2][0],
			   Yy=R[0].ref_v_xyz[2][1],
			   Yz=R[0].ref_v_xyz[2][2];
		
		double x0=R[0].target->C[R[0].protect[0]].x,
			   y0=R[0].target->C[R[0].protect[0]].y,
			   z0=R[0].target->C[R[0].protect[0]].z;
		double dx0=x0-t0[0],
			   dy0=y0-t0[1],
			   dz0=z0-t0[2];
		yangle1=3.14159-acos((Yx*dx0 + Yy*dy0 + Yz*dz0)/modu(Yx,Yy,Yz)/modu(dx0,dy0,dz0));
	}
	
	{
		double Yx=R[1].ref_v_xyz[2][0],
			   Yy=R[1].ref_v_xyz[2][1],
			   Yz=R[1].ref_v_xyz[2][2];
		
		double x0=R[1].target->C[R[1].protect[0]].x,
			   y0=R[1].target->C[R[1].protect[0]].y,
			   z0=R[1].target->C[R[1].protect[0]].z;
		double dx0=x0-t1[0],
			   dy0=y0-t1[1],
			   dz0=z0-t1[2];
		yangle2=3.14159-acos((Yx*dx0 + Yy*dy0 + Yz*dz0)/modu(Yx,Yy,Yz)/modu(dx0,dy0,dz0));
	}
	this->siteI_direction=fabs(yangle1+yangle2-0.67071509866412082);
*/
	
	//---------------Artificial Reaction Coordinate------------
	double D1=digitalneg(r, 3.0, 2.3/1.);
	double D2=digitalneg(r, 1.0, 2.3/0.2)
		     *digitalneg(AxisBeta, 40./180*3.14, 2.3/(10./180*3.14))
		     *digitalneg(RadiusBeta, 40./180*3.14, 2.3/(10./180*3.14));

	Q = r;
	if (r<7.0){
		Q = Q -
			put(AxisBeta+RadiusBeta, 30./180*3.14 + 30./180*3.14, 5/(30./180*3.14 + 30./180*3.14))
			*D1;
		if (r<4.0 && AxisBeta<90/180.0*3.14159 && RadiusBeta<90/180.0*3.14159){
			Q = Q -
				(
			    put(r_siteI_deviation, 2. , 5./2.)
//			   +put(siteI_direction, 90./180*3.14 , 2.5/(90./180*3.14))
				)*D2; 
		}
	}

	if (Q<IEV_EFFECTIVE_THRESHOLD){
		this->IEV_spheres_effective = 1;
	}else{
		this->IEV_spheres_effective = 0;
	}

	//---------------Reshape Biasing potential-----------------
	this->unbiasedE = 0;
	double E11,E12,E21,E22,Er;

	E11=E12=E21=E22=0;

	Er=0; //(-0.0)*exp(-r*r/2/8.0/8.0);

	if (r<7.0){
		E11=-put(AxisBeta+RadiusBeta, 180./180*3.14 + 180./180*3.14, 5./(180./180*3.14 + 180./180*3.14));
		//E11= (-5)*exp(-AxisBeta*AxisBeta/(2*(20./180.*3.14)*(20./180.*3.14)));
		//E12= (-5)*exp(-RadiusBeta*RadiusBeta/(2*(20./180.*3.14)*(20./180.*3.14)));

		if (r<4.0 && AxisBeta<90/180.0*3.14159 && RadiusBeta<90/180.0*3.14159){
			E21=-put(r_siteI_deviation, 6. , 10.0/2.);
			E22=0;//-put(siteI_direction, 360./180*3.14 , 10.0/(360./180*3.14));
			//E21=(-8)*exp(-(siteI_direction-2.0)*(siteI_direction-2.0)/(2*(1.0*1.0)) );
			//E22=(-8)*exp(-r_siteI*r_siteI/(2*(20./180.*3.14)*(20./180.*3.14)));

		}
	} 
	
	this->unbiasedE = Er+(E11+E12)*D1 + (E21+E22)*D2;
	this->E = this->unbiasedE +  U.getBiasingE(Q);


	
//	A potential from Quan Du and A. Vologodskii.
/*	double sigma2=1;
	double r0=0.1,q=0.8;
	double A=17;

	this->E=
	A*exp(-((AxisBeta-PI)*(AxisBeta-PI)+(RadiusBeta-PI)*(RadiusBeta-PI))/2/sigma2)
	*(pow(r0/r,q*2)-2*pow(r0/r,q));



	double r0=0,r0std=11.0,a0=10;
	double E0=-a0*exp(-(r-r0)*(r-r0)/2/r0std);
	//double E11= fabs(AxisBeta-50/180.0*3.14159) * (0);
	//double E12= fabs(RadiusBeta-30/180.0*3.14159) * (0);

	double E21= fabs(AxisBeta-0/180.0*3.14159) * (5);
	double E22= fabs(RadiusBeta-0/180.0*3.14159) * (5);
	double E4= r_siteI * 0;
	double E5= siteI_direction * (0);

	this->E=E0
		//+(E11+E12)*digital(r,2.5,3)
		+(E21+E22)*(1-digital(r,2.0,																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																						3))
		+E4+E5;
*/


    return this->E;
}


cls_rigid::cls_rigid(CircularChain* r_target, std::vector<long> r_protect,
					 std::vector < vector<double> > r_ref_v, std::vector <long> r_ref_vec_basis, long r_anchor):
target(r_target),protect(r_protect),ref_v(r_ref_v), ref_vec_basis(r_ref_vec_basis), anchor(r_anchor){
	//When cls_rigid is constructed, it will automatically update the global
	//variable "protect_list", which is the copy of .protect member (vector)
	// of cls_rigid's container class "all_rigid".
	//This copy is made for better global access of protect list.
	double Mv[3][3];
	Mv[0][0]=target->C[protect[this->ref_vec_basis[0]]].dx;
	Mv[1][0]=target->C[protect[this->ref_vec_basis[0]]].dy;
	Mv[2][0]=target->C[protect[this->ref_vec_basis[0]]].dz;

	Mv[0][1]=target->C[protect[this->ref_vec_basis[1]]].dx;
	Mv[1][1]=target->C[protect[this->ref_vec_basis[1]]].dy;
	Mv[2][1]=target->C[protect[this->ref_vec_basis[1]]].dz;

	Mv[0][2]=target->C[protect[this->ref_vec_basis[2]]].dx;
	Mv[1][2]=target->C[protect[this->ref_vec_basis[2]]].dy;
	Mv[2][2]=target->C[protect[this->ref_vec_basis[2]]].dz;

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
	Mv[0][0]=target->C[protect[this->ref_vec_basis[0]]].dx;
	Mv[1][0]=target->C[protect[this->ref_vec_basis[0]]].dy;
	Mv[2][0]=target->C[protect[this->ref_vec_basis[0]]].dz;

	Mv[0][1]=target->C[protect[this->ref_vec_basis[1]]].dx;
	Mv[1][1]=target->C[protect[this->ref_vec_basis[1]]].dy;
	Mv[2][1]=target->C[protect[this->ref_vec_basis[1]]].dz;

	Mv[0][2]=target->C[protect[this->ref_vec_basis[2]]].dx;
	Mv[1][2]=target->C[protect[this->ref_vec_basis[2]]].dy;
	Mv[2][2]=target->C[protect[this->ref_vec_basis[2]]].dz;

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
