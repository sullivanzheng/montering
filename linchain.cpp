//////////////////////////////
#include "chain.h" 
#include <assert.h>
using namespace std;
int LinearChain::updateAllBangleKinkNum()
{
	stats.kink_num = 0;
	for (int i = 1; i < maxnum + 1; i++)
	{
		C[i].bangle = calAngle(C[i - 1], C[i]);
		if (C[i].bangle > KINKLOWERBOUND)
			stats.kink_num++;
	}/*
    cout <<"============Bangle and KinkNum Updated=============="<<endl;
	for (int i = 1; i < maxnum + 1; i++)
        cout <<"| "<< C[i].bangle << endl;
    cout<<"------------------------------"<<endl; */
    return 0;
}
int LinearChain::updateBangleKinkNum(int i)
{   
    if (i==0){
#ifdef DEBUG
        cout <<"bangle updata 0, ignored"<<endl;
#endif
        return 0;
    }
	if (i < 0 || i > maxnum)
	{
		printf("bangle update (num): refered to non-existing segment %d\n", i);
		getchar();
        exit(EXIT_FAILURE);
	}
	//double temp=C[i].bangle;
	if (C[i].bangle > KINKLOWERBOUND)
		stats.kink_num--;
	C[i].bangle = calAngle(C[i - 1], C[i]);
	if (C[i].bangle > KINKLOWERBOUND)
		stats.kink_num++;
	return 0;
}

LinearChain::LinearChain(int length)
:Chain(true,length)
{}

LinearChain::LinearChain(char const *filename,int length)
:Chain(filename,false,length){}

double LinearChain::calG_bSum()
{
	double temp = 0;
	for (int i = 1; i < maxnum + 1; i++)
		temp += G_b(C[i].bangle);
	return temp;
}
double LinearChain::deltaE_TrialCrankshaft(int m, int n, double a)
{   
	if (m<0 || n<1 || m>maxnum || n>maxnum + 1 || m>n)
	{
		cout<<"Illegal values of m and n ("<<m<<','<<n<<')';
		exit(EXIT_FAILURE);
	}
	//how many steps of moves have been performed.
	double M[3][3];
	SetRotM_crankshaft(M,m,n,a);
	//M is now the rotation matrix.
	double dE;
    double newA1,oldA1,newA2,oldA2;
    if (m>0){
        segment C1;oldA1 = C[m].bangle;
	    C1 = C[m];   
        double v[3];
	    double temp[3];
        v[0] = C1.dx;
        v[1] = C1.dy;
        v[2] = C1.dz;
        mat33mulvec3(M, v, temp);
        C1.dx = temp[0];
        C1.dy = temp[1];
        C1.dz = temp[2];
        newA1 = calAngle(C[m-1], C1);
    }
    else{
        newA1=oldA1=0;
    }
    
    if (n<maxnum+1){
	    segment C2;	oldA2 = C[n].bangle;
	    C2 = C[n-1];
	    double v[3];
	    double temp[3];
        v[0] = C2.dx;
	    v[1] = C2.dy;
	    v[2] = C2.dz;
	    mat33mulvec3(M, v, temp);
	    C2.dx = temp[0];
	    C2.dy = temp[1];
	    C2.dz = temp[2];
	    newA2 = calAngle(C2, C[n]);
    }
    else{    
        newA2=oldA2=0;
    }
    /*
	int new_kink_num = stats.kink_num + 
		(
		-(oldA1 > KINKLOWERBOUND ? 1 : 0) 
		- (oldA2 > KINKLOWERBOUND ? 1 : 0) 
		+ (newA1 > KINKLOWERBOUND ? 1 : 0) 
		+ (newA2 > KINKLOWERBOUND ? 1 : 0));
	double const C = 2.4e-28 / (1.3806503e-23 * 300) / 3.4e-10;
	//C=3e-19 erg.cm=3e-28 J.m. Change this unit to kT.basepairlength*/
    dE =+G_b(newA1)- G_b(oldA1) 
        +G_b(newA2)- G_b(oldA2);
		/* TORTIONAL STRESS TAKEN INTO CONSIDERATION
		+ 2 * PI * PI * C / (totsegnum * bpperseg) * 
		(
		+(Lk - totsegnum * bpperseg / 10.5 + DELTA_TW_K * new_kink_num) 
		* (Lk - totsegnum * bpperseg / 10.5 + DELTA_TW_K * new_kink_num) 
		- (Lk - totsegnum * bpperseg / 10.5 + DELTA_TW_K * stats.kink_num) 
		* (Lk - totsegnum * bpperseg / 10.5 + DELTA_TW_K * stats.kink_num));
		*/
#ifdef DEBUG
    cout<<"trialCrankshaft: "<<m<<"  "<<n<<endl;
#endif
	return dE;
}

int LinearChain::crankshaft(int m, int n, double a, 
                            bool trialMoveFlag)
{
	//Crankshaft for an angle
	//The segment is from X(m) to X(n).
	//M is the rotation matrix.
    //0<=m<=maxnum-1,1<=n<=maxnum
	if (m<0||n<1 || m>maxnum||n>maxnum+1 || m>n)
	{
		cout<<"Illegal values of m and n ("<<m<<','<<n<<')';
		exit(EXIT_FAILURE);
	}
    if (!trialMoveFlag)	{
        stats.accepts++;
        stats.crankSuccess++;
    }
	double M[3][3];
	SetRotM_crankshaft(M,m,n,a);
	int i,j;
	i = m;
	while (i != n)
	{
		double v[3];
		double temp[3];
		v[0] = C[i].dx;	v[1] = C[i].dy;	v[2] = C[i].dz;
		mat33mulvec3(M, v, temp);
		C[i].dx = temp[0];	C[i].dy = temp[1];	C[i].dz = temp[2];
		j=i+1;
        if (j<maxnum+1){
		    C[j].x = C[i].x + C[i].dx;
		    C[j].y = C[i].y + C[i].dy;
		    C[j].z = C[i].z + C[i].dz;
        }
		i = j;
	}
   	updateBangleKinkNum(m);
    if (n<maxnum+1)	updateBangleKinkNum(n);
#ifdef DEBUG
    cout <<"++AcceptedCrankshaft"<<endl;
#endif
	return 0;
}


double LinearChain::deltaE_TrialHalfChain(int m, double rv[3],double a)
{	
	if (m<0||m>maxnum)
	{
		cout<<"Illegal values of m ("<<m<<')';
		exit(EXIT_FAILURE);
	}

    if (m==0) return 0.0;
	//how many steps of moves have been performed.
	double M[3][3];
	SetRotM_halfchain(M,rv,a);
	//M is now the rotation matrix.
	double dE;
	double newA,oldA;
	segment C_;
	C_ = C[m];
	double v[3];
	double temp[3];
	v[0] = C_.dx;
	v[1] = C_.dy;
	v[2] = C_.dz;
	mat33mulvec3(M, v, temp);
	C_.dx = temp[0];
	C_.dy = temp[1];
	C_.dz = temp[2];
	newA = calAngle(C[m-1], C_);
	oldA = C[m].bangle;
	
	/*int new_kink_num = stats.kink_num + 
		(
		-(oldA > KINKLOWERBOUND ? 1 : 0) 
		+ (newA > KINKLOWERBOUND ? 1 : 0));
	double const C = 2.4e-28 / (1.3806503e-23 * 300) / 3.4e-10;
	//C=3e-19 erg.cm=3e-28 J.m. Change this unit to kT.basepairlength */
	dE = (+G_b(newA)- G_b(oldA));
		/* NO TORTIONAL STRESS TAKEN INTO CONSIDERATION 
		+ 2 * PI * PI * C / (totsegnum * bpperseg) * 
		(
		+(Lk - totsegnum * bpperseg / 10.5 + DELTA_TW_K * new_kink_num) 
		* (Lk - totsegnum * bpperseg / 10.5 + DELTA_TW_K * new_kink_num) 
		- (Lk - totsegnum * bpperseg / 10.5 + DELTA_TW_K * stats.kink_num) 
		* (Lk - totsegnum * bpperseg / 10.5 + DELTA_TW_K * stats.kink_num));*/
#ifdef DEBUG
    cout <<"trialHalfChain "<<m<<endl;
#endif
	return dE;
}
int LinearChain::halfChain(int m, double rv[3],double a,bool trialMoveFlag)
{
	if (m<0||m>maxnum)
	{
		cout<<"Illegal values of m ("<<m<<')';
		exit(EXIT_FAILURE);
	}

    if (!trialMoveFlag) {
        stats.accepts++;
        this->stats.halfChainSuccess++;
    }
    //if (m==0) return 0;
	double M[3][3];
	SetRotM_halfchain(M,rv,a);
	int i,j;
	for (i=m;i<maxnum;i++)
	{
		double v[3];
		double temp[3];
		v[0] = C[i].dx;	v[1] = C[i].dy;	v[2] = C[i].dz;
		mat33mulvec3(M, v, temp);
		C[i].dx = temp[0];	C[i].dy = temp[1];	C[i].dz = temp[2];
		j=i+1;
		C[j].x = C[i].x + C[i].dx;
		C[j].y = C[i].y + C[i].dy;
		C[j].z = C[i].z + C[i].dz;
	}
	double v[3];
	double temp[3];
	v[0] = C[maxnum].dx;v[1] = C[maxnum].dy;v[2] = C[maxnum].dz;
	mat33mulvec3(M, v, temp);
	C[maxnum].dx = temp[0];	C[maxnum].dy = temp[1];	C[maxnum].dz = temp[2];
	updateBangleKinkNum(m);
#ifdef DEBUG
    cout<<"++AcceptedHalfChain"<<endl;
#endif
	return 0;
}

bool LinearChain::trialLigateAfterHalfChainOK(int m, double rv[3],double a,
        double endToEndDistanceThreshold,
        double endToEndAngleThreshold) {
    if (m<0||m>maxnum)
	{
		cout<<"Illegal values of m ("<<m<<')';
		exit(EXIT_FAILURE);
	}
    //If both ends are far away from each other.
    if (m==0){
        if (this->getEndToEndDistance()>endToEndDistanceThreshold) {
            cout <<'F';
            return false;
        }
        if (this->getEndToEndAngle()>endToEndAngleThreshold)
            return false;
        else
            return true;
    }

    double M[3][3];
	SetRotM_halfchain(M,rv,a);

    double v[3]={
        C[maxnum].x+C[maxnum].dx-C[m].x,
        C[maxnum].y+C[maxnum].dy-C[m].y,
        C[maxnum].z+C[maxnum].dz-C[m].z
    };
    double temp[3];
    mat33mulvec3(M, v,temp);
    temp[0]+=C[m].x;
    temp[1]+=C[m].y;
    temp[2]+=C[m].z;
    if (    abs(temp[0]-C[0].x)>endToEndDistanceThreshold
        ||  abs(temp[1]-C[0].y)>endToEndDistanceThreshold
        ||  abs(temp[2]-C[0].z)>endToEndDistanceThreshold
        ||  modu(temp[0]-C[0].x,temp[1]-C[0].y,temp[2]-C[0].z)
            >endToEndDistanceThreshold){
            return false;
    }
    segment C_=C[maxnum];
    double v2[3]={C_.dx,C_.dy,C_.dz},temp2[3];
    mat33mulvec3(M,v2,temp2);
    C_.dx=temp2[0];C_.dy=temp2[1];C_.dz=temp2[2];
    if (calAngle(C[0],C_)>endToEndAngleThreshold)
        return false;
    else
//        cout<<"good move"<<endl;
        return true;
}

bool LinearChain::trialLigateAfterCrankshaftOK(int m, int n, double a, 
   double endToEndDistanceThreshold,
   double endToEndAngleThreshold){
    if (m<0||m>maxnum || n<1 || n>maxnum+1 || m>n )
    {
		cout<<"Crankshaft Illegal values of m and n("<<m<<','<<n<<')'<<endl;
		exit(EXIT_FAILURE);
	}

    //If both ends are far away from each other.
    if (this->getEndToEndDistance()>endToEndDistanceThreshold) {
        cout <<'F';
        return false;
    }
    
    //If the crankshaft rotation does not change the endToEndAngle at all.
    if ((m!=0 && n!=maxnum+1) 
        || (m==0 && n==maxnum+1)) {
        /* cout <<"MH"; */ 
        if (this->getEndToEndAngle()>endToEndAngleThreshold)
            return false;
        else
            return true;
    }

    //If the crankshaft rotation change the endToEndAngle.
    //cout<<'E';
    double M[3][3];
    SetRotM_crankshaft(M,m,n,a);
    if ( m==0 ) {
        double v[3]={C[0].dx,C[0].dy,C[0].dz};
        double temp[3];
        mat33mulvec3(M,v,temp);
        segment C_;
        C_.dx=temp[0];
        C_.dy=temp[1];
        C_.dz=temp[2];
        if (calAngle(C_,C[maxnum])>endToEndAngleThreshold) 
            return false;
        else
            return true;
    }
    else if(n==maxnum+1){
        double v[3]={C[maxnum].dx,C[maxnum].dy,C[maxnum].dz};
        double temp[3];
        mat33mulvec3(M,v,temp);
        segment C_;
        C_.dx=temp[0];
        C_.dy=temp[1];
        C_.dz=temp[2];
        if (calAngle(C_,C[0])>endToEndAngleThreshold) 
            return false;    
        else
            return true;
    }

    //if no proper return value, exception.
    assert(1==0);
}