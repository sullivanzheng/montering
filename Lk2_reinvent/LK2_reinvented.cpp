#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <cassert>
#include "chain.h"

using namespace std;
const int maxintsec=200;


struct cross{
	int seg_over;
	int seg_under;
	int sign;
	double s_over;
	double s_under;
	double endpos;
};

struct Generator{
	int id;
	int cutby;
	int pre;
	int next;
	bool alive;
	cross cr;
};

int DET_Doolittle_LU_Decomp_Pivot(double da[maxintsec][maxintsec], int n){
	const double eps=1e-12;
	int pivot[maxintsec];
	double max;

	for(int k=0;k<n;k++){
		//find pivot row.      
		max=fabs(da[k][k]);
		pivot[k]=k;
		for(int i=k+1;i<n;i++){
			if(fabs(da[i][k]) > max){
			   max=fabs(da[i][k]);
			   pivot[k]=i;
			}
		}
	
          
		//pivoting: switch row k with pivot row.
		if (pivot[k]!=k){
			for (int i=0; i<n; i++){
				max=da[k][i];
				da[k][i]=da[pivot[k]][i];
				da[pivot[k]][i]=max;
			}
		}

		//if singularity occurs, return error.
		if (fabs(da[k][k])<eps){
			return -999999;
		}
          
		//lower triangular matrix elements of column k.
		for (int i=k+1;i<n;i++){
			da[i][k]=da[i][k]/da[k][k];
		}
         
		//upper triangular matrix elements          
		for(int i=k+1;i<n;i++){
			for (int j=k+1;j<n;j++){
				da[i][j]=da[i][j]-da[i][k]*da[k][j];
			}
		}
	}

	double c=1.0;
	for (int i=0; i<n; i++){
		c*=da[i][i];
	}
	return c;
}

bool cmpfun(cross c1, cross c2){
	return (c1.endpos<c2.endpos);
}

cross intsec(int i, int j, double x[], double y[], double z[]){
///algorithm ref:http://softsurfer.com/Archive/algorithm_0104/algorithm_0104B.htm
///if yes, return +1 if + crossing, -1 if - crossing; if no crossing, return 0.

	const double eps=1e-12;
	cross X={-999,0,0,0.0,0.0,0.0};

	double d,u1=x[i+1]-x[i],u2=y[i+1]-y[i],
		     v1=x[j+1]-x[j],v2=y[j+1]-y[j],
			 w1=x[i]  -x[j],w2=y[i]  -y[j];

    d=v1*u2-v2*u1; //vXu.
    if (fabs(d)<eps) return X;
	
	double s=(v2*w1-v1*w2)/d;
    if (s<0 || s>=1) return X;

	double t=(u1*w2-u2*w1)/(-d);
    if (t<0 || t>=1) return X;

	if ( z[i]+s*(z[i+1]-z[i]) < z[j]+t*(z[j+1]-z[j]) ){
		X.seg_under=i;
		X.seg_over=j;
		X.sign= d>0?1:-1;
		X.s_under=s;
		X.s_over=t;
		X.endpos=X.seg_under+X.s_under;
	}else{
		X.seg_under=j;
		X.seg_over=i;
		X.sign= d>0?-1:1;
		X.s_under=t;
		X.s_over=s;
		X.endpos=X.seg_under+X.s_under;
	}
    return X;
}

int CircularChain::AP(long vertM, long vertN,double s, double t){
/* L1 is the number of segments in the first contour.
   L2 is the number of segments in the second contour.
   s,t are parameters of AlexanderPolynomial(s,t) */
	

	//Store the intersections in different fashions.
	//Mind that x[0]=x[L1] x[L1+1]=x[L1+L2+1]
	const double pi=3.141592653589793;

	long L1=vertM+(maxnum-vertN+1);//Number of seg in the 1st circle
	long L2=vertN-vertM;           //Number of seg in the 2nd circle
	double x[maxa],y[maxa],z[maxa];

	cross cr[maxintsec],cr_temp;
	double id[maxintsec],ix[maxintsec];
	int L20=L1+1,L2end=L1+L2+1;

	long i,j;
	for (i=0,j=0;i<=vertM-1;i++,j++){
		x[j]=C[i].x;y[j]=C[i].y;z[j]=C[i].z;
	}

	for (i=vertN,j=vertM;i<=maxnum;i++,j++){
		x[j]=C[i].x;y[j]=C[i].y;z[j]=C[i].z;
	}
	
	x[l1]=C[0].x;y[l1]=C[0].y;z[l1]=C[0].z;

	for (i=vertM,j=l1+1;i<=vertN-1;i++,j++){
		x[j]=C[i].x;y[j]=C[i].y;z[j]=C[i].z;
	}
	x[l1+l2+1]=C[vertM].x;
	y[l1+l2+1]=C[vertM].y;
	z[l1+l2+1]=C[vertM].z;


	//Enumerate all intersections
	int N=0;
	for (int i=0;i<=L2end-3;i++){
		if (i==L1) continue;
		for (int j=i+2;j<=L2end-1;j++){
			if (j==L1) continue;
			if ((i==0 && j==L1-1) || (i==L20 && j==L2end-1)) continue;
			cr_temp=intsec(i,j,x,y,z);
			if (cr_temp.seg_over ==-999){
				;
			}
			else{
				cr[N] = cr_temp;
				++N;
			}
		}
	}

	if (N==0) return 0; //Special occasion: no intersection at all.

	//Sort intersections
	sort(cr,cr+N,cmpfun);

	//Enumerate generators
	int g[maxintsec]; //if generator i is cut by j then g[i]=j

	int M=0;//contour one has generator from 0~M.
	if (cr[0].endpos > L1 || cr[N-1].endpos < L1) return 0; //Special occasion: self-intersections only.
	else{
		while (!(cr[M].endpos<L1 && cr[M+1].endpos>L1) ) {
			M++;
			if (M==N-1) break; //Overflow prevention.
		}
	}
	assert (M < N-1 && M >= 0);

	for (int i=0; i<N; i++){
		double overpos=cr[i].s_over + cr[i].seg_over;
		assert (!(overpos>L1 && overpos<L1+1));
		if (overpos < cr[0].endpos) g[i]=0;
		else if (overpos > cr[N-1].endpos) g[i]=M+1;
		else if (cr[M].endpos < overpos && overpos<L1) g[i]=0;
		else if (L1 < overpos && overpos < cr[M+1].endpos) g[i]=M+1;
		else{
			assert ( (cr[0].endpos<overpos && overpos<cr[M].endpos)
				   ||(cr[M+1].endpos<overpos && overpos<cr[N-1].endpos));
			int mm=1;
			while (!(cr[mm-1].endpos<overpos && overpos<cr[mm].endpos)) mm++;
			g[i]=mm;
			assert(g[i]!=0 && g[i]!=M+1);
		}
	}

	//Disentanglement
	Generator G[maxintsec];

	for (int i=0; i<N; i++){
		G[i].id=i;
		G[i].cutby=g[i];
		G[i].pre=i-1;
		G[i].next=i+1;
		G[i].alive=true;
		G[i].cr=cr[i];
	}
	G[0].pre=M;G[M].next=0;
	G[M+1].pre=N-1;G[N-1].next=M+1;	

	int totdisentang=0;
	int L1num=M+1,L2num=N-1-(M+1)+1;
	bool exitflag=false;
	
	//goto DA;
	while (!exitflag){
		exitflag=true;

//		Type Ia: g[i]==i	 Type Ib: g[i]=i+1 i and i+1 on the same contour.
		bool exitflagI=false;
		while (!exitflagI){
			exitflagI=true;
			for (int i=0; i<N; i++){
				if (!G[i].alive) continue;
				if (G[i].cutby==G[i].id || G[i].cutby==G[i].next){
					totdisentang++;
					if (G[i].id<=M) L1num--; else L2num--;
					if (L1num==0 || L2num==0) return 0;
					G[i].alive=false;
					for (int j=0;j<N;j++){
						if (!G[j].alive) continue;
						if (G[j].cutby==G[i].id) G[j].cutby=G[G[i].next].id;
					}
					G[G[i].pre].next=G[G[i].next].id;					
					G[G[i].next].pre=G[G[i].pre].id;
					i--;
					exitflag=false;exitflagI=false;
				}
			}
		}
		;
		//Type II g[i]==j and g[i+1]==j
		for (int i=0; i<N; i++){
			if (!G[i].alive) continue;
			if (G[i].cutby==G[G[i].next].cutby){
				for (int r=0; r<N; r++){
					if (!G[r].alive) continue;
					if (G[r].cutby==G[i].next) 
						goto skip;
				}
				totdisentang++;

				if (G[i].id<=M) L1num-=2; else L2num-=2;
				if (L1num==0 || L2num==0) return 0;

				G[i].alive=false;
				G[G[i].next].alive=false;
				
				int next2=G[G[i].next].next;
				for (int j=0;j<N;j++){
					if (!G[j].alive) continue;
					if (G[j].cutby==G[i].next || G[j].cutby==G[i].id) {
						G[j].cutby=next2;
					}
				}

				G[next2].pre=G[i].pre;
				G[G[i].pre].next=next2;
				i--;
				exitflag=false;
			}
skip:		;
		}
	}

DA:	int newid=0;
	for (int i=0; i<N; i++){
		if (G[i].alive){
			G[i].id=newid;
			cr[newid]=cr[i];
			newid++;
		}
	}
	newid--;

	assert (newid+1==L1num+L2num);

	for (int i=0; i<N; i++){
		if (!G[i].alive) continue;
		g[G[i].id]=G[G[i].cutby].id;
	}
	M=L1num-1;
	N=L1num+L2num;
	
	//Fill up the matrix
	double da[maxintsec][maxintsec];
	for (int k=0; k<N; k++)
		for (int i=0; i<N; i++)
			da[k][i]=0.;

	for (int k=0; k<N; k++){
		int i=g[k];
		if(k<=M && M>0) {//---------------------------------------------------
			int k1=k+1;
			if(k==M) k1=0; //wrap around.
			if (i==k || i==k1){     //a. Trivial intersections TODO: k+1 or k1
				da[k][k]=-1;
				da[k][k1]=1;
			}
			else if(i<=M) {   		//b. Internal intersection within L1.
				if(cr[k].sign < 0) {      //Type I or (-) crossing
				  da[k][k]=1;
				  da[k][k1]=-s;
				  da[k][i]=s-1;
				}else{				    //Type II or (+) crossing
				  da[k][k]=-s;
				  da[k][k1]=1;
				  da[k][i]=s-1;
				}
			}
			else{ assert (i>M);		//c. L1 generator cut by L2
				if(cr[k].sign < 0){        //Type I or (-) crossing
				  da[k][k]=1;
				  da[k][k1]=-t;
				  da[k][i]=s-1;
				}else{				     //Type II or (+) crossing
				  da[k][k]=-t;
				  da[k][k1]=1;
				  da[k][i]=s-1;
				}
			}
		}
		else if(k==M && M==0){//---------------------------------------------------
			assert(i>M);
			da[k][k]=1-t;
			da[k][i]=s-1;
		}
		else if(k>M && N-1 > M+1){//---------------------------------------------------
			int k1=k+1;
			if(k==N-1) k1=M+1;
			if (i==k || i==k1){       //a. Trivial intersections TODO: k+1 or k1
				da[k][k]=-1;
				da[k][k1]=1;
			}
			else if(i>M) {			  //b. Internal intersection within L2.
				if(cr[k].sign < 0) {
					da[k][k]=1;
					da[k][k1]=-t;
					da[k][i]=t-1;
				}else{
					da[k][k]=-t;
					da[k][k1]=1;
					da[k][i]=t-1;
				}
			}else{				//c. L2 generator cut by L1.
				assert (i<=M);
				if(cr[k].sign < 0){
					da[k][k]=1;
					da[k][k1]=-s;
					da[k][i]=t-1;
				}else{
					da[k][k]=-s;
					da[k][k1]=1;
					da[k][i]=t-1;
				}
			}
		}
		else if(k==N-1 && N-1==M+1) {//----------------------------------------------
			assert (i<=M);
			da[k][k]=1-s;
			da[k][i]=t-1;
		}
	}

	double temp = DET_Doolittle_LU_Decomp_Pivot(da,N-1)/(1-t); //only feed N-1 x N-1 remainder
	long re=floor(temp+0.5);
	while (re%2 ==0) re=re/2;
	return re;
}

int main(){
	ofstream fpdbg;
	double x[maxa],y[maxa],z[maxa];
	fpdbg.open("dbg.log");
	ifstream fp("Cor_for.txt",ifstream::in);
	int i=0;
	while (fp.good()){
		fp>>x[i]>>y[i]>>z[i];
		fpdbg <<"Line "<<i<<" read. "<<x[i]<<' '<<y[i]<<' '<<z[i]<<endl;
		++i;
	}
	fp.close();
	double temp=AP(3,212,-1.0,-2.0,0.0,x,y,z);
	fpdbg <<"RESULT  "<<temp<<endl;
	fpdbg.close();
	return 0;
}


