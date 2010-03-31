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

CircularChain::CircularChain(char const *filename,const int length)
:Chain(filename,true,length)
{
	Lk = (int)(bpperseg * totsegnum / 10.5 + 0.5);
}

CircularChain::CircularChain(int length):Chain(true,length){
	Lk = (int)(bpperseg * totsegnum / 10.5 + 0.5);
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
	updateBangleKinkNum(m);
	updateBangleKinkNum(n);
	return 0;
}

double CircularChain::deltaE_TrialCrankshaft_countMove(int m, int n, double a)
{
	this->auto_updt_stats();
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

	int new_kink_num = stats.kink_num + 
        (-(oldA1 > KINKLOWERBOUND ? 1 : 0) 
         -(oldA2 > KINKLOWERBOUND ? 1 : 0) 
         + (newA1 > KINKLOWERBOUND ? 1 : 0) 
         + (newA2 > KINKLOWERBOUND ? 1 : 0));
	double const C = 2.4e-28 / (1.3806503e-23 * 300) / 3.4e-10;
	//C=3e-19 erg.cm=3e-28 J.m. Change this unit to kT.basepairlength
	dE = (+G_b(newA1) + G_b(newA2) - G_b(oldA1) - G_b(oldA2));/* 
        + 2 * PI * PI * C / (totsegnum * bpperseg) * 
        (+(Lk - totsegnum * bpperseg / 10.5 + DELTA_TW_K * new_kink_num) 
        *(Lk - totsegnum * bpperseg / 10.5 + DELTA_TW_K * new_kink_num) 
        -(Lk - totsegnum * bpperseg / 10.5 + DELTA_TW_K * stats.kink_num) 
        * (Lk - totsegnum * bpperseg / 10.5 + DELTA_TW_K * stats.kink_num));*/
//Torsional stress removed for now.

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
	sprintf(buf, "%6d %20s", maxnum + 1, "[Snapshot of a circle]");
	fh << buf << endl;
	i=0;
	sprintf(buf, "%6d%4s%12.6f%12.6f%12.6f%6d%6d%6d", 
        i + 1, "F", C[i].x, C[i].y, C[i].z, 1, (i == 0 ? maxnum + 1 : i), (i + 1 == maxnum + 1 ? 1 : i + 2));
	fh << buf << endl;
	for ( i = 1; i <= maxnum; i++)
	{
		sprintf(buf, "%6d%4s%12.6f%12.6f%12.6f%6d%6d%6d", 
            i + 1, "C", C[i].x, C[i].y, C[i].z, 1, (i == 0 ? maxnum + 1 : i), (i + 1 == maxnum + 1 ? 1 : i + 2));
		fh << buf << endl;
	}

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


allrigid::allrigid(char *configfile,CircularChain * target){
	using namespace std;
	ifstream f(configfile);
	if (!f.good()){
		cout<<"Rigidbody configfile does not exist!"<<std::endl;
		exit(EXIT_FAILURE);
	}
	
	vector<int> r_protect;
	vector<vector<double>>r_ref_v;

	while (!f.eof()){
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

	}
	for (int i=0;i<this->R.size();i++)
		for (int j=0;j<this->R[i].protect.size();j++)
			this->protect.push_back(R[i].protect[j]);
	this->update_allrigid_and_E();
}

double allrigid::update_allrigid_and_E(){
	//This section can be customized for different rigid body set.
	for (int i=0;i<this->R.size();i++){
		R[i].update_ref_v();
	}
	this->E=modu2(R[0].target->C[R[0].protect[1]].x - R[1].target->C[R[1].protect[1]].x,
		R[0].target->C[R[0].protect[1]].y - R[1].target->C[R[1].protect[1]].y,
		R[0].target->C[R[0].protect[1]].z - R[1].target->C[R[1].protect[1]].z);
	return this->E;
}


cls_rigid::cls_rigid(CircularChain* r_target, std::vector<int> r_protect,
	std::vector < vector<double> > r_ref_v):
target(r_target),protect(r_protect),ref_v(r_ref_v){
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

	for (int i=0;i<this->ref_v_xyz.size();i++){
		std::cout<<ref_v_xyz[i][0]<<" ";
		std::cout<<ref_v_xyz[i][1]<<" ";
		std::cout<<ref_v_xyz[i][2]<<" "<<std::endl;
	}
}

void cls_rigid::update_ref_v(){
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

	/* test printouts. */
	
	for (int i=0;i<this->ref_v_xyz.size();i++){
		std::cout<<">"<<this->target->stats.auto_moves()<<std::endl;
		
		std::cout<<this->target->C[protect[1]].dx<<" ";
		std::cout<<this->target->C[protect[1]].dy<<" ";
		std::cout<<this->target->C[protect[1]].dz<<" "<<std::endl;

		std::cout<<ref_v_xyz[i][0]<<" ";
		std::cout<<ref_v_xyz[i][1]<<" ";
		std::cout<<ref_v_xyz[i][2]<<" "<<std::endl;
	}
/*
rigid.txt is as follows:

#rigid body 1
rigidstart
$ 8 9 10
vec 0,1,0
rigidend


#rigid body 2
rigidstart
$ 108 109 110
vec 0,1,0
rigidend

Should see the same vector.
*/

}