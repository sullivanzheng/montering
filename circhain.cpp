//////////////////////////////////
#include "chain.h"
using namespace std;
int CircularChain::updateAllBangleKinkNum()
{
	stats.kink_num = 0;
	C[0].bangle = calAngle(C[maxnum], C[0]);
	if (C[0].bangle > KINKLOWERBOUND)
		stats.kink_num++;
	for (int i = 1; i < maxnum + 1; i++)
	{
		C[i].bangle = calAngle(C[i - 1], C[i]);
		if (C[i].bangle > KINKLOWERBOUND)
			stats.kink_num++;
	}
    cout <<"============Bangle and KinkNum Updated=============="<<endl;
	for (int i = 1; i < maxnum + 1; i++)
        cout <<"| "<< C[i].bangle << endl;
    cout<<"------------------------------"<<endl;
    return 0;
}
int CircularChain::updateBangleKinkNum(int i)
{
	if (i > maxnum || i < 0)
	{
		printf("bangle update (num): refered to non-existing segment %d", i);
		getchar();
		exit(EXIT_FAILURE);
	}
	//double temp=C[i].bangle;
	if (C[i].bangle > KINKLOWERBOUND)
		stats.kink_num--;
	if (i == 0)
		C[0].bangle = calAngle(C[maxnum], C[0]);
	else
		C[i].bangle = calAngle(C[i - 1], C[i]);
	if (C[i].bangle > KINKLOWERBOUND)
		stats.kink_num++;
	//if (C[i].bangle>15.*PI/180.) 
	//	cout <<"UPDT "<<i<<' '<<temp*180./PI<<"->"<<C[i].bangle*180./PI<<endl;
	return 0;
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

CircularChain::CircularChain(char const *filename,int length)
:Chain(filename,true,length)
{
	Lk = (int)(bpperseg * totsegnum / 10.5 + 0.5);
}

double CircularChain::calG_bSum()
{
	double temp = 0;
	for (int i = 0; i < maxnum + 1; i++)
		temp += G_b(C[i].bangle);
	return temp;
}
int CircularChain::crankshaft(int m, int n, double a, 
                              bool trialMoveFlag)
{
	//Crankshaft for an angle
	//The segment is from X(m) to X(n).
	//M is the rotation matrix.
	stats.accepts++;
	if ((m<0||n<0)||(n>maxnum || m>maxnum))
	{
		cout<<"Illegal values of m and n ("<<m<<','<<n<<')';
		exit(EXIT_FAILURE);
	}
	//how many moves are accepted.
	if (fabs(C[0].x) > (totsegnum) / 2 || fabs(C[0].y) > (totsegnum) / 2 || fabs(C[0].z) > (totsegnum) / 2)
	{
		cout << "Step #" << stats.moves() << " Drift proof." << endl;
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
	updateBangleKinkNum(m);
	updateBangleKinkNum(n);
	return 0;
}

double CircularChain::deltaE_TrialCrankshaft(int m, int n, double a)
{
	countmove();
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
	/*if (fabs(newA1)>22*PI/180.0 || fabs(newA2)>22*PI/180.0) {
	cout <<endl<<"trykink!"
	<<endl<<oldA1*180.0/PI
	<<"->"<<newA1*180.0/PI<<'\t'
	<<oldA2*180.0/PI<<"->"
	<<newA2*180.0/PI<<endl;
	}*/
	int new_kink_num = stats.kink_num + 
        (-(oldA1 > KINKLOWERBOUND ? 1 : 0) 
         -(oldA2 > KINKLOWERBOUND ? 1 : 0) 
         + (newA1 > KINKLOWERBOUND ? 1 : 0) 
         + (newA2 > KINKLOWERBOUND ? 1 : 0));
	double const C = 2.4e-28 / (1.3806503e-23 * 300) / 3.4e-10;
	//C=3e-19 erg.cm=3e-28 J.m. Change this unit to kT.basepairlength
	dE = (+G_b(newA1) + G_b(newA2) - G_b(oldA1) - G_b(oldA2)) 
        + 2 * PI * PI * C / (totsegnum * bpperseg) * 
        (+(Lk - totsegnum * bpperseg / 10.5 + DELTA_TW_K * new_kink_num) 
        *(Lk - totsegnum * bpperseg / 10.5 + DELTA_TW_K * new_kink_num) 
        -(Lk - totsegnum * bpperseg / 10.5 + DELTA_TW_K * stats.kink_num) 
        * (Lk - totsegnum * bpperseg / 10.5 + DELTA_TW_K * stats.kink_num));
	/*
	C_=2*pi*pi*C;
	bt_angle=bt_angle/360;bpo_angle=bpo_angle/360;
	temp=4.0*pi*pi*140/2/N+C_/N*(d_Lk+bt_angle*0+bpo_angle)^2;
	*/
	return dE;
}