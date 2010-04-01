
#include "chain.h"
using namespace std;
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
:Chain(filename,true,length){}

CircularChain::CircularChain(int length):Chain(true,length){}

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


int CircularChain::IEV( int in,  int ik){
// excluded volume effects
// iev=0 corresponds to intersection
// iev=1 corresponds to no intersection
// This subroutine is translated from A. Vologodskii Monte FORTRAN 77 program.

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
	double er,er2;
	er=this->VolEx_R;

	for (int i=in;i<=ik;i++){				// do 2 i=in,ik
		for (int j=0;j<=maxnum;j++){		// do 3 j=1,jr1
			if (j >= in && j <= ik) continue;//if(j.ge.in.and.j.le.ik) goto 3
			if (abs(i-j) <= VEcutoff) continue;//(if(iabs(i-j).le.lll) goto 3
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

allrigid::allrigid(char *configfile,CircularChain * target){
	using namespace std;
	ifstream f(configfile);
	if (!f.good()){
		cout<<"Rigidbody configfile does not exist!"<<std::endl;
		exit(EXIT_FAILURE);
	}
	
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
	for (int i=0;i<this->R.size();i++){
		R[i].update_ref_v_xyz();
	}

	this->E = 0;

	//Axis Orientation. ref_v_xyz[0]
	this->E += betaVec1Vec2(this->R[0].ref_v_xyz[0],this->R[1].ref_v_xyz[0])*(-10);

	//Curvature radius direction orientation. ref_v_xyz[1]
	this->E += betaVec1Vec2(this->R[0].ref_v_xyz[1],this->R[1].ref_v_xyz[1])*(-10);

	//Center point of each rigid body. ref_v_xyz[2]
	double t0[3],t1[3];
	t0[0]=R[0].ref_v_xyz[2][0] + R[0].target->C[R[0].protect[0]].x;
	t0[1]=R[0].ref_v_xyz[2][1] + R[0].target->C[R[0].protect[0]].y;
	t0[2]=R[0].ref_v_xyz[2][2] + R[0].target->C[R[0].protect[0]].z;

	t1[0]=R[1].ref_v_xyz[2][0] + R[1].target->C[R[1].protect[0]].x;
	t1[1]=R[1].ref_v_xyz[2][1] + R[1].target->C[R[1].protect[0]].y;
	t1[2]=R[1].ref_v_xyz[2][2] + R[1].target->C[R[1].protect[0]].z;
	
	this->E += modu(t1[0]-t0[0],t1[1]-t0[1],t1[2]-t0[2]) * 50;

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

	std::cout<<"Rigid body constructed"<<std::endl;
	for (int i=0;i<this->ref_v_xyz.size();i++){
		
		std::cout<<ref_v_xyz[i][0]<<" ";
		std::cout<<ref_v_xyz[i][1]<<" ";
		std::cout<<ref_v_xyz[i][2]<<" "<<std::endl;
	}
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

	/* test printouts. 
	
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