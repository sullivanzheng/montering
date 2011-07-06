#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <cassert>

using namespace std;
const int maxa=900;
const int maxintsec=200;
ofstream fpdbg;

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

int Doolittle_LU_Decomposition_with_Pivoting(double *A, int pivot[], int n)
{
////////////////////////////////////////////////////////////////////////////////
//  int Doolittle_LU_Decomposition_with_Pivoting(double *A, int pivot[],      //
//                                                                    int n)  //
//                                                                            //
//  Description:                                                              //
//     This routine uses Doolittle's method to decompose the n x n matrix A   //
//     into a unit lower triangular matrix L and an upper triangular matrix U //
//     such that A = LU.                                                      //
//     The matrices L and U replace the matrix A so that the original matrix  //
//     A is destroyed.  Note!  In Doolittle's method the diagonal elements of //
//     L are 1 and are not stored.                                            //
//     The LU decomposition is convenient when one needs to solve the linear  //
//     equation Ax = B for the vector x while the matrix A is fixed and the   //
//     vector B is varied.  The routine for solving the linear system Ax = B  //
//     after performing the LU decomposition for A is                         //
//                      Doolittle_LU_with_Pivoting_Solve.                     //
//                                                                            //
//     The Doolittle method with partial pivoting is:  Determine the pivot    //
//     row and interchange the current row with the pivot row, then assuming  //
//     that row k is the current row, k = 0, ..., n - 1 evaluate in order the //
//     following pair of expressions                                          //
//       U[k][j] = A[k][j] - (L[k][0]*U[0][j] + ... + L[k][k-1]*U[k-1][j])    //
//                                 for j = k, k+1, ... , n-1                  //
//       L[i][k] = (A[i][k] - (L[i][0]*U[0][k] + . + L[i][k-1]*U[k-1][k]))    //
//                          / U[k][k]                                         //
//                                 for i = k+1, ... , n-1.                    //
//       The matrix U forms the upper triangular matrix, and the matrix L     //
//       forms the lower triangular matrix.                                   //
//                                                                            //
//  Arguments:                                                                //
//     double *A       Pointer to the first element of the matrix A[n][n].    //
//     int    pivot[]  The i-th element is the pivot row interchanged with    //
//                     row i.                                                 //
//     int     n       The number of rows or columns of the matrix A.         //
//                                                                            //
//  Return Values:                                                            //
//     0  Success                                                             //
//    -1  Failure - The matrix A is singular.                                 //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double A[N][N];                                                        //
//     int    pivot[N];                                                       //
//                                                                            //
//     (your code to create matrix A)                                         //
//     err = Doolittle_LU_Decomposition_with_Pivoting(&A[0][0], pivot, N);    //
//     if (err < 0) printf(" Matrix A is singular\n");                        //
//     else { printf(" The LU decomposition of A is \n");                     //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//      
   int i, j, k, p;
   double *p_k, *p_row, *p_col;
   double max;


//    For each row and column, k = 0, ..., n-1,
 
   for (k = 0, p_k = A; k < n; p_k += n, k++) {

//    find the pivot row

      pivot[k] = k;
      max = fabs( *(p_k + k) );
      for (j = k + 1, p_row = p_k + n; j < n; j++, p_row += n) {
         if ( max < fabs(*(p_row + k)) ) {
            max = fabs(*(p_row + k));
            pivot[k] = j;
            p_col = p_row;
         }
      }

//    and if the pivot row differs from the current row, then
//    interchange the two rows.
   
      if (pivot[k] != k)
         for (j = 0; j < n; j++) {
            max = *(p_k + j);
            *(p_k + j) = *(p_col + j);
            *(p_col + j) = max;
         }

//    and if the matrix is singular, return error


      if ( *(p_k + k) == 0.0 ) return -1;

//    otherwise find the lower triangular matrix elements for column k. 

      for (i = k+1, p_row = p_k + n; i < n; p_row += n, i++) {
         *(p_row + k) /= *(p_k + k);
      }  

//    update remaining matrix

      for (i = k+1, p_row = p_k + n; i < n; p_row += n, i++)
         for (j = k+1; j < n; j++)
            *(p_row + j) -= *(p_row + k) * *(p_k + j);
   }

   return 0;
}



void printGenerators(Generator G[], int N){
	//fpdbg <<"------------------------"<<endl;
	//for (int i=0; i<N; i++){
	//	if (!G[i].alive) continue;
	//	fpdbg << "G[" <<i<<"] " << (G[i].alive?'.':'X') 
	//		<<" cutby "	<< G[i].cutby << " <" << G[i].pre << ',' << G[i].next << "> " <<endl;
	//}
	//fpdbg <<endl;

	for (int i=0; i<N; i++)
		if (G[i].alive==true) assert(G[G[i].cutby].alive==true);

	Generator GG[maxintsec];
	int j=0;
	for (int i=0; i<N; i++){
		if (!G[i].alive) continue;
		GG[j]=G[i];
		j++;
	}

	//int count=0;
	//for (int i=1; i<j-1; i++){
	//	if (GG[i].next!=GG[i+1].id || GG[i].pre!=GG[i-1].id){
	//		fpdbg << "----GG[" <<i<<"] id "<<GG[i].id <<" cutby "	<< GG[i].cutby << " <" << GG[i].pre << ',' << GG[i].next << "> " <<endl;
	//		count++;
	//	}
	//}

	//if (count>2){
	//	fpdbg <<"PROBLEM!!!"<<endl;
	//	fpdbg <<endl;
	//}

	int newid=0;
	for (int i=0; i<N; i++){
		if (G[i].alive){
			GG[i]=G[i];
			GG[i].id=newid;
			newid++;
		}
	}
	newid--;

	int g[maxintsec];
	fpdbg <<"---Generators in reduced form---"<<endl;
	for (int i=0; i<N; i++){
		if (!G[i].alive) continue;
		g[GG[i].id]=GG[GG[i].cutby].id;
		fpdbg<<"ix["<<GG[i].id+1<<"] "<<g[GG[i].id]+1
			<<" cross: "<<GG[i].cr.sign<<" "<<GG[i].cr.seg_under<<"  "<<GG[i].cr.s_under<<endl;
	}
	fpdbg<<"###"<<endl<<endl;
}
int printGenerators_checkDisentangsSegment(Generator G[], int N, int ID){
	int newid=0;
	for (int i=0; i<N; i++){
		if (!G[i].alive) continue;
		newid++;
		if (G[i].id==ID) return newid;
	}
	return -1;
}

void saveMatrixToFile(double da[maxintsec][maxintsec], int N, char* filename){
	//da is the matrix to print. M is a NxN matrix indexed from (0,0) to (N-1, N-1).
	ofstream fp(filename);
	fp << "a=["<<endl;
	for (int i=0;i<N;i++){
		for (int j=0;j<N;j++)
			fp <<da[i][j]<<' ';
		fp << endl;
	}
	fp <<"];"<<endl;
	fp <<"a=a(1:end-1,1:end-1);"<<endl;
	fp <<"det(sym(a))"<<endl;
	fp <<"s=det(sym(a))/sym(3);"<<endl;
	fp <<"while mod(s,2)==0;s=s/2;end;"<<endl;
	fp <<"s"<<endl;
	fp.close();
}

double abs_det_lu(double *A, int n){
	int pivot[maxintsec];
	Doolittle_LU_Decomposition_with_Pivoting(A,pivot,n);
	double det=1;
	for (int i=0;i<n;i++) {
		double c=*(A+n*i+i);
		//fpdbg <<c<<endl;
		det=det*c;
	}
	//fpdbg <<"Determinant="<<det;
	
	//fpdbg <<det-long(det)<<endl;
	//while (fabs(det/2.0-long(det/2.0))<1e-13) det=det/2.0;
	return det;
}

int readGenerator(char* filename, int g[], int &M){
	int N=0;
	ifstream fp(filename);
	while (fp.good()){
		fp>>g[N];
		if (g[N]==-1) { M=N-1; continue;}
		N++;
	}
	return N;
}



bool cmpfun(cross c1, cross c2){
	return (c1.endpos<c2.endpos);
}

cross intsec(int i, int j, double x[maxa], double y[maxa], double z[maxa]){
///algorithm ref:http://softsurfer.com/Archive/algorithm_0104/algorithm_0104B.htm
///if yes, return +1 if + crossing, -1 if - crossing; if no crossing, return 0.

	const double eps=1e-12;
	cross X={-999,0,0,0.0,0.0,0.0};

	double d,u1=x[i+1]-x[i],u2=y[i+1]-y[i],
		     v1=x[j+1]-x[j],v2=y[j+1]-y[j],
			 w1=x[i]-x[j],  w2=y[i]-y[j];

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

int AP(int L1,int L2,double s, double t, double alpha, 
	   double x1[maxa], double y1[maxa], double z1[maxa]){
/* L1 is the number of segments in the first contour.
   L2 is the number of segments in the second contour.
   s,t are parameters of AlexanderPolynomial(s,t) */
	

/*	inventarization of intersections
	s1 s2(i) - s-value of i-th intersection
	ic1(i) - number of undergoing segment for i-th intersection
	ic2(i) - number of overgoing segment for i-th intersection
	id(i) - overpassing or underpassing
	ix(i) - generator number */

/*  Mind that x[0]=x[L1] x[L1+1]=x[L1+L2+1] */
	const double pi=3.141592653589793;
	double x[maxa],y[maxa],z[maxa];

	alpha=alpha/180.0*pi;
	for (int i=0; i<=L1+L2+1; i++){
		x[i]=x1[i];
		y[i]=cos(alpha)*y1[i]-sin(alpha)*z1[i];
		z[i]=sin(alpha)*y1[i]+cos(alpha)*z1[i];
		//fpdbg <<x[i]<<' '<<y[i]<<' '<<z[i]<<endl;
	}

	cross cr[maxintsec],cr_temp;
	double id[maxintsec],ix[maxintsec];
	int L20=L1+1,L2end=L1+L2+1;

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

	//for (int i=0; i<N; i++){
	//	fpdbg <<"X["<<i<<"]"<<cr[i].seg_under<<" cut by "
	//		<<cr[i].seg_over<<" Position "<<cr[i].endpos<<endl;
	//}

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
	
	//for (int i=0; i<N; i++){
	//	fpdbg <<"X[gen# "<<i<<"  cut by "<<g[i]<<"]   X-type:  "<<cr[i].sign/*<<cr[i].seg_under<<" cut by "
	//		<<cr[i].seg_over<<" Position "<<cr[i].endpos*/<<endl;
	//}
	//


	//Disentanglement
	Generator G[maxintsec];

	//N=readGenerator("gen.txt",g,M);
	//goto DA;
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

	printGenerators(G,N);
	int totdisentang=0;
	int L1num=M+1,L2num=N-1-(M+1)+1;
	bool exitflag=false;
	const int TER=0;
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
					fpdbg <<"["<<totdisentang<<"]Type I at "<<i<<endl;
					if (totdisentang > TER) {fpdbg<<totdisentang<<endl;goto t2;}
	fpdbg<<"IN FORTRAN EXPRESSION: "<<printGenerators_checkDisentangsSegment(G,N,i)<<endl;
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
	printGenerators(G,N);
					
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
				fpdbg <<"["<<totdisentang<<"]Type II at "<<i<<endl;
				if (totdisentang > TER) {fpdbg<<totdisentang<<endl;goto t2;}
				fpdbg<<"IN FORTRAN EXPRESSION: "<<printGenerators_checkDisentangsSegment(G,N,i)<<endl;

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
printGenerators(G,N);
				i--;
				exitflag=false;
				//TODO: restore this. exitflagII=false;
			}
skip:		;
		}
	}
t2:	fpdbg <<"<<<<Final Generators>>>>"<<endl;
	printGenerators(G,N);

	int newid=0;
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
DA:	double da[maxintsec][maxintsec];
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


	saveMatrixToFile(da,N,"matrix.txt");
	double lin[maxintsec*maxintsec];
	int p=0;
	for (int i=0;i<N-1;i++){
		for (int j=0;j<N-1;j++){
			lin[p]=da[i][j];
			p++;
		}
	}
	double temp = abs_det_lu(&lin[0],N-1)/(1-t); //only feed N-1 x N-1 remainder
	long re=floor(temp+0.5);
	while (re%2 ==0) re=re/2;
	return re;
}

int main(){
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
	double temp=AP(111,111,-1.0,-2.0,0.0,x,y,z);
	fpdbg <<"RESULT  "<<temp<<endl;
	fpdbg.close();
	return 0;
}


