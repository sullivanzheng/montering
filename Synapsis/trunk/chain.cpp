#include "chain.h" 
using namespace std;


int Chain::readIniFile(char const *filename)
{
	ifstream file_st (filename);
    cout << endl<<"Chain::readIniFile: Reading Ini coordinate file: "<<filename<<endl
        <<"Number of atoms or segments: "<<maxnum+1<<endl;
    if (maxnum<10 || maxnum>maxa-1){
        cout <<"Too many or too few atoms"<<endl;
        getchar();
        exit(EXIT_FAILURE);
    }
	if (!file_st.good())
	{
		cout << "File do not exist!";
		getchar();
		exit(EXIT_FAILURE);
	}
	C[0].x = C[0].y = C[0].z = 0;
	file_st >> C[0].dx;
	for (int i = 1; i < maxnum + 1; i++)
	{
		file_st >> C[i].dx;
		C[i].x = C[i - 1].x + C[i - 1].dx;
	}
	file_st >> C[0].dy;
	for (int i = 1; i < maxnum + 1; i++)
	{
		file_st >> C[i].dy;
		C[i].y = C[i - 1].y + C[i - 1].dy;
	}
	file_st >> C[0].dz;
	for (int i = 1; i < maxnum + 1; i++)
	{
		file_st >> C[i].dz;
		C[i].z = C[i - 1].z + C[i - 1].dz;
	}

	gCal ins_gCal;

	for (int i=0; i <=maxnum ; i++){
		C[i].l=modu(C[i].dx,C[i].dy,C[i].dz);
		double g1=ins_gCal.FindGbyK(DNA_KUHN/(bpperunit*C[i].l));
		double g2=ins_gCal.FindGbyK(DNA_KUHN/(bpperunit*C[wrap(i-1,totsegnum)]);
		C[i].g=1.0/(1.0/g1+1.0/g2);
		cout<<"[Initializing rigidity] l="<<C[i].l<<" g="<<C[i].g<<endl;
	}

	this->total_length = 0;
	for (int i=0; i <= maxnum ; i++){
		this -> total_length += C[i].l;
	}
	file_st.close();
	return 0;
}

void Chain::SetRotM_halfchain(double M[3][3], double rv[3], double a)
{
	//rv 
	double cosa;
	double sina;
	cosa = cos(a);
	sina = sin(a);
	M[0][0]=1.0;	M[0][1]=0.0;	M[0][2]=0.0;
	M[1][0]=0.0;	M[1][1]=1.0;	M[1][2]=0.0;
	M[2][0]=0.0;	M[2][1]=0.0;	M[2][2]=1.0;
	mat33mulscalOW(M, cosa);
	//cos(a)*eye(3);
	double M2[3][3];
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			M2[i][j] = rv[i] * rv[j];
	mat33mulscalOW(M2, (1 - cosa));
	double M3[3][3] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
	/// <>/////////////</>
	/*0*/
	M3[0][1] = rv[2];
	M3[0][2] = -rv[1];
	M3[1][0] = -rv[2] /*0*/;
	M3[1][2] = rv[0];
	M3[2][0] = rv[1];
	M3[2][1] = -rv[0];//cross product matrix
	/*0*/
	/// <>/////////////</>
	mat33mulscalOW(M3, sina);
	mat33addOW(M, M2);
	mat33addOW(M, M3);
	//See http://en.wikipedia.org/wiki/Euler_Angles for algorithm
	//See also http://en.wikipedia.org/wiki/Rotation_around_a_fixed_axis
	//And of course the general idea is to decompose any vector x to x//e (e is rotation vector)
	//and x_perpendicular to e (x_). x_ Xprd with e you have a x_2 perp. to x_ and e.
	//Express the rotated x using combination of x// x_ and x_2 and convert this algorithm to 
	//a matrix.
}

void Chain::SetRotM_crankshaft(double M[3][3],int m, int n, double a)
{
	double rv_norm;
	double rv[3];
	//rotation vector;
    if (n<maxnum+1){
        rv[0] = C[n].x - C[m].x;
	    rv[1] = C[n].y - C[m].y;
	    rv[2] = C[n].z - C[m].z;
    }
    else if (n==maxnum+1)
    {
        rv[0] = C[maxnum].x + C[maxnum].dx - C[m].x;
	    rv[1] = C[maxnum].y + C[maxnum].dy - C[m].y;
        rv[2] = C[maxnum].z + C[maxnum].dz - C[m].z;
    }
    else
    {
        cout <<"SetRotM_crankshaft: Illegal n number:"<<n<<endl;
        exit(EXIT_FAILURE);
    }
	rv_norm = sqrt(dot_product(rv, rv));
	rv[0] = rv[0] / rv_norm;
	rv[1] = rv[1] / rv_norm;
	rv[2] = rv[2] / rv_norm;
	//normalize rotation vector
	SetRotM_halfchain(M, rv, a);
}

void Chain::normalize(){
    for (int i=0;i<=maxnum-1;i++){
        double temp=modu(C[i].dx,C[i].dy,C[i].dz);
        C[i].dx=C[i].dx/temp*C[i].l;
        C[i].dy=C[i].dy/temp*C[i].l;
        C[i].dz=C[i].dz/temp*C[i].l;
        C[i+1].x=C[i].x+C[i].dx;
        C[i+1].y=C[i].y+C[i].dy;
        C[i+1].z=C[i].z+C[i].dz;
    }
    double temp=modu(C[maxnum].dx,C[maxnum].dy,C[maxnum].dz);
    C[maxnum].dx=C[maxnum].dx/temp*C[maxnum].l;
    C[maxnum].dy=C[maxnum].dy/temp*C[maxnum].l;
    C[maxnum].dz=C[maxnum].dz/temp*C[maxnum].l;
}

double Chain::calAngle(segment &C1, segment &C2)
{
	double temp;
	double angle;
	temp = C1.dx * C2.dx + C1.dy * C2.dy + C1.dz * C2.dz;
    temp = temp/(modu2(C1.dx,C1.dy,C1.dz)*modu(C2.dx,C2.dy,C2.dz));
	if (temp > 1 || temp < -1)
	{
		printf("Step# %d: Warning, cos(bangle)=%5.3f is larger than 1\n", stats.auto_moves, temp);
		angle = 0.0;
	}
	else
		angle = acos(temp);
	return angle;
}

Chain::Chain(char const *filename, bool circular,int r_numofseg)
{   
    this->endToEndSampleCycle=this->defaultSampleCycle;
    if (r_numofseg>maxa||r_numofseg<5){
         cout<<"The chain is too long or too short. Simulation aborted."<<endl
             <<"Chain with 5~910 segments are allowed."<<endl;
         getchar();
         exit(EXIT_FAILURE);
    }
    maxnum=(this->_numofseg=r_numofseg)-1;
	readIniFile(filename);
    //updateAllBangleKinkNum(); Calling virtual function is dangerous.		
	updateAllBangle_Ini(circular);
    stats.resetStat();
}


//updateAllBangle_Ini is supposed to be used in constructor.		
//Therefore it contains a "curcular" parameter that make it not virtual.		
void Chain::updateAllBangle_Ini(bool circular = false)		
{		
	if (circular){		
		C[0].bangle = calAngle(C[maxnum], C[0]);		
		for (int i = 1; i < maxnum + 1; i++)		
		{		
			C[i].bangle = calAngle(C[i - 1], C[i]);		
		}		
	}		
	else		
	{		
		for (int i = 1; i < maxnum + 1; i++)		
		{		
			C[i].bangle = calAngle(C[i - 1], C[i]);		
		}		
		C[0].bangle=0.0;		
	}		
	cout <<"============BangleInitiated=============="<<endl;		
	for (int i = 1; i < maxnum + 1; i++)		
		cout <<'['<<i<<']'<< C[i].bangle;
	cout<<"[END]"<<endl;
}
void Chain::snapshot(char *filename)
{
	//A not-working well version
	ofstream fh (filename);
	char buf[300];
	if (!fh.good())
	{
		cout << "file not writable" << endl;
		getchar();
		exit(EXIT_FAILURE);
	}
	sprintf(buf, "%6d %20s", maxnum + 1, "[General Chain Snapshot, treat as Linear]");
	fh << buf << endl;
    int	i=0;
	sprintf(buf, "%6d%4s%12.6f%12.6f%12.6f%6d%6d", 
        i + 1, "F", C[i].x, C[i].y, C[i].z, 1, 1);
	fh << buf << endl;
	for (i = 1; i <= maxnum-1; i++)
	{
		sprintf(buf, "%6d%4s%12.6f%12.6f%12.6f%6d%6d%6d", 
            i + 1, "C", C[i].x, C[i].y, C[i].z, 1, (i == 0 ? maxnum + 1 : i), (i + 1 == maxnum + 1 ? 1 : i + 2));
		fh << buf << endl;
	}
    sprintf(buf, "%6d%4s%12.6f%12.6f%12.6f%6d%6d", 
        maxnum + 2, "N", 
        C[maxnum].x + C[maxnum].dx,
        C[maxnum].y + C[maxnum].dy,
        C[maxnum].z + C[maxnum].dz, 1, maxnum);
	fh << buf << endl;

	//Detailed information.
    fh << endl << "Detailed Info" << endl;
	for (i = 0; i < maxnum; i++)
	{
		sprintf(buf, "seg %d |dX(%15.13f %15.13f %15.13f)|=%15.10f X[i+1]-X=%15.10f l=%15.10f b_ang=%15.10f", 
            i, C[i].dx,C[i].dy,C[i].dz,
            modu(C[i].dx, C[i].dy, C[i].dz), 
            modu(C[i + 1].x - C[i].x, 
                 C[i + 1].y - C[i].y, 
                 C[i + 1].z - C[i].z),
		    C[i].l,
            C[i].bangle);
		fh << buf << endl;
	}
	i=maxnum;
	sprintf(buf, "seg %d |dX(%15.13f %15.13f %15.13f)|=%15.10f X[i+1]-X=%15.10f l=%15.10f b_ang=%15.10f", 
            i+1, C[i].dx,C[i].dy,C[i].dz,
            modu(C[i].dx, C[i].dy, C[i].dz), 
            "----------",
			C[i].l,
			C[i].bangle);
	fh << buf << endl;
	fh.close();
}