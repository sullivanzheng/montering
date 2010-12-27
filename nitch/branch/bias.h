#ifndef BIAS_H
#define BIAS_H

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
extern double tempE(double);

class BiasingPotential{
public:
	static const int MAXBIN=500;
	
	double G[MAXBIN+2],K[MAXBIN+2];
	unsigned long HS[MAXBIN+2],f[MAXBIN+2];
	
	double binstart,binend;
	int binnum;
	unsigned long counter;

	long energy_update_cycle;

	double HHR;

private:
	double binsize;
	BiasingPotential();
	std::ofstream fG;
	std::ofstream fK;
	std::ofstream fHS;
	std::ofstream ff;
	std::ofstream fx;
	inline int index(double a){
		if(a<binstart)
			return 0;
		else if (a>=binend)
			return binnum+1;
		else{
			return int((a-binstart)/binsize)+1;
		}
	}

	int updateG(){
		double sumK=0;
		for (int i=0; i<=binnum+1; i++) sumK+=K[i];
		for (int i=0; i<=binnum+1; i++) G[i]=-log(K[i]*(binnum+2)/sumK);

		return 0;
	}

	double updateRepository_HHRc(){
		unsigned long HSmin=HS[1],HSmax=HS[1];
		for (int i=0; i<=binnum+1;i++){
			//ff<<f[i]<<' ';
		}
		//ff<<std::endl;

		for (int i=0; i<=binnum+1;i++){
			HS[i]+=f[i];
			K[i]+=double(f[i])*exp(-G[i]);
			f[i]=0;

			if (i>=2 && i<=binnum){
				if (HS[i]<HSmin) HSmin=HS[i];
				else if (HS[i]>HSmax) HSmax=HS[i];
			}
		}


		double HHRc = (HSmin==0)?1e5:double(HSmax)/double(HSmin);
		std::cout<<"HHRc="<<HHRc<<std::endl;
		return HHRc;
	}

	double updateall(){
		for (double x=0;x<=55;x+=0.1){
			//fG<<-getBiasingE(x)<<' ';
		}
		//fG<<std::endl;

		double HHRc=this->updateRepository_HHRc();
		if (HHRc>10 || HHRc > 1.10 * this->HHR){
			this->updateG();
			if (this->HHR>10) {
				for (int i=0; i<=binnum+1;i++) K[i]*=1.0;
			}
		}
		else{
			if (HHRc<this->HHR) this->HHR=HHRc;
		}
		

		for (int i=0; i<=binnum+1;i++){
			
			//fHS<<HS[i]<<' ';
			//fK<<K[i]<<' ';
		}

		//fHS<<std::endl;
		//fK<<std::endl;
		return this->HHR;
	}

public:
	BiasingPotential(double rbinstart, double rbinend, double rbinnum, long renergy_update_cycle):
	binstart(rbinstart),binend(rbinend),binnum(rbinnum),HHR(1e5),counter(0),energy_update_cycle(renergy_update_cycle),
	fG("fG.txt"), fK("fK.txt"), fHS("fHS.txt"), ff("ff.txt"),fx("fx.txt"){
		if (rbinnum>MAXBIN){
			std::cout<<"Exceed the maximum number of bins"<<std::endl;
			exit(3);
		}

		for (int i=0;i<=binnum+1;i++) {
			K[i]=0.0001;
			HS[i]=0;f[i]=0;
		}
		binsize=(binend-binstart)/binnum;
		this->updateG();
	}

	~BiasingPotential(){
		fG.close();fK.close();fHS.close();ff.close();fx.close();
	}

	inline double collect(double a){
		
		counter++;
		if (counter%1==0){
			f[index(a)]++;//fx<<a<<std::endl;
		}
		if (counter==energy_update_cycle){
			//fx<<std::endl;
			this->updateall();
			counter=0;
			return this->HHR;
		}
		return -999;
	}

	inline double getBiasingE(double a){
		int n=index(a);
		return -G[n];
		/*if (n==0 || n==binnum+1)
			return -G[n];
		double G11=G[n-1],G0=G[n],G1=G[n+1];
		double a0=a-(binstart+binsize*(n-1));
		double p=a0/binsize,out;
		if (p<0.5){
			out=(0.5-p)*G11+(0.5+p)*G0;
		}
		else{
			out=(p-0.5)*G1+(1.5-p)*G0;
		}
		return -out;*/
	}

	int pickle(char* filename){
		std::ofstream fp(filename);
		if (!fp.good()){
			std::cout<<"BiasingPotential::pickle could not write to file:"<<filename<<std::endl;
			exit(3);
		}
		char buf[400];
		sprintf(buf,"%4s %3s %8s %12s %12s %14s %10s","#","SYM","binstart","G","K","HS","f");
		fp<<buf<<std::endl;
		for (int i=0; i<=binnum+1; i++){
			sprintf(buf,"%4d %3s %8.3f %12g %12g %14d %10d",
				i,i==0?"NEG":(i==binnum+1?"POS":"---"),
				binstart+(i-1)*binsize,G[i],K[i],HS[i],f[i]);
			fp<<buf<<std::endl;
		}
		fp <<"---Column end---"<<std::endl;
		fp <<"binstart "<<binstart<<std::endl;
		fp <<"binend "<<binend<<std::endl;
		fp <<"binnum "<<binnum<<std::endl;
		fp <<"counter "<<counter<<std::endl;
		fp <<"HHR "<<HHR<<std::endl;

		fp.close();
		return 0;
	}

	int load(char* filename){
		std::ifstream fp(filename);
		if (!fp.good()){
			std::cout<<"BiasingPotential::load could not load file:"<<filename<<std::endl;
			exit(3);
		}
		char buf[400];
		fp.getline(buf,400);

		int i=0;
		while(!fp.eof()){
			fp.getline(buf,400);
			if (strcmp(buf,"---Column end---")==0) {
				break;
			}
			std::stringstream s(buf);
			s>>i>>buf>>buf;
			s>>G[i]>>K[i]>>HS[i]>>f[i];
			std::cout<<i<<" G= "<<G[i]<<std::endl;
		}

		
		fp.getline(buf,400);
		std::stringstream(buf)>>buf>>binstart;
		fp.getline(buf,400);
		std::stringstream(buf)>>buf>>binend;
		fp.getline(buf,400);
		std::stringstream(buf)>>buf>>binnum;
		fp.getline(buf,400);
		std::stringstream(buf)>>buf>>counter;
		fp.getline(buf,400);
		std::stringstream(buf)>>buf>>HHR;

		binsize=(binend-binstart)/binnum;
		fp.close();

		return 0;
	}
};


class BiasingPotential2D{
public:
	static const int MAXBIN_m=500,MAXBIN_n=100;
	
	double G[MAXBIN_m+2][MAXBIN_n+2],K[MAXBIN_m+2][MAXBIN_n+2];
	unsigned long HS[MAXBIN_m+2][MAXBIN_n+2],f[MAXBIN_m+2][MAXBIN_n+2];
	
	double binstart_m,binstart_n,binend_m,binend_n;
	int binnum_m,binnum_n;
	unsigned long counter;

	long energy_update_cycle;

	double HHR;

private:
	double binsize_m,binsize_n;
	BiasingPotential2D();
	std::ofstream fG;
	std::ofstream fK;
	std::ofstream fHS;
	std::ofstream ff;
	std::ofstream fx;
	inline int index_m(double m){
		if(m<binstart_m)
			return 0;
		else if (m>=binend_m)
			return binnum_m+1;
		else{
			return int((m-binstart_m)/binsize_m)+1;
		}
	}
	inline int index_n(double n){
		if(n<binstart_n)
			return 0;
		else if (n>=binend_n)
			return binnum_n+1;
		else{
			return int((n-binstart_n)/binsize_n)+1;
		}
	}

	int updateG(){
		double sumK=0;
		for (int i=0; i<=binnum_m+1; i++){
			for (int j=0;j<=binnum_n+1;j++){
				sumK+=K[i][j];
			}
		}
		
		for (int i=0; i<=binnum_m+1; i++){
			for (int j=0;j<=binnum_n+1;j++){
				G[i][j]=-log(K[i][j]*(binnum_m+2)*(binnum_n+2)/sumK);
			}
		}
		
		return 0;
	}

	double updateRepository_HHRc(){
		unsigned long HSmin=HS[1][1],HSmax=HS[1][1];
		for (int i=0; i<=binnum_m+1;i++){
				for (int j=0; j<=binnum_n+1; j++){
				//ff<<f[i]<<' ';
				}
		}
		//ff<<std::endl;

		for (int i=0; i<=binnum_m+1;i++){
				for (int j=0; j<=binnum_n+1; j++){
					HS[i][j]+=f[i][j];
					K[i][j]+=double(f[i][j])*exp(-G[i][j]);
					f[i][j]=0;

					if (i>=2 && i<=binnum_m && j>=2 && j<=binnum_n){
						if (HS[i][j]<HSmin) HSmin=HS[i][j];
						else if (HS[i][j]>HSmax) HSmax=HS[i][j];
				}
			}
		}


		double HHRc = (HSmin==0)?1e5:double(HSmax)/double(HSmin);
		std::cout<<"HHRc="<<HHRc<<std::endl;
		return HHRc;
	}

	double updateall(){
		for (double x=0;x<=55;x+=0.1){
			//fG<<-getBiasingE(x)<<' ';
		}
		//fG<<std::endl;

		double HHRc=this->updateRepository_HHRc();
		if (HHRc>10 || HHRc > 1.10 * this->HHR){
			this->updateG();
			if (this->HHR>10) {
				for (int i=0; i<=binnum_m+1;i++) for(int j=0; j<=binnum_n+1;j++) K[i][j]*=1.0;
			}
		}
		else{
			if (HHRc<this->HHR) this->HHR=HHRc;
		}
		

		for (int i=0; i<=binnum_m+1;i++){
			for (int j=0; j<=binnum_n+1;j++){
			
			//fHS<<HS[i]<<' ';
			//fK<<K[i]<<' ';
			}
		}

		//fHS<<std::endl;
		//fK<<std::endl;
		return this->HHR;
	}

public:
	BiasingPotential2D(double rbinstart_m, double rbinend_m, double rbinnum_m, 
					 double rbinstart_n, double rbinend_n, double rbinnum_n,
					 long renergy_update_cycle):
	binstart_m(rbinstart_m),binend_m(rbinend_m),binnum_m(rbinnum_m),
    binstart_n(rbinstart_n),binend_n(rbinend_n),binnum_n(rbinnum_n),
		HHR(1e5),counter(0),energy_update_cycle(renergy_update_cycle),
	fG("fG.txt"), fK("fK.txt"), fHS("fHS.txt"), ff("ff.txt"),fx("fx.txt"){
		if (rbinnum_m>MAXBIN_m || rbinnum_n>MAXBIN_n){
			std::cout<<"Exceed the maximum number of bins"<<std::endl;
			exit(3);
		}

		for (int i=0;i<=binnum_m+1;i++) {
			for (int j=0;j<=binnum_n+1;j++){
				K[i][j]=0.0001;
				HS[i][j]=0;f[i][j]=0;
			}
		}
		binsize_m=(binend_m-binstart_m)/binnum_m;
		binsize_n=(binend_n-binstart_n)/binnum_n;
		this->updateG();
	}

	~BiasingPotential2D(){
		fG.close();fK.close();fHS.close();ff.close();fx.close();
	}

	inline double collect(double a,double b){
		
		counter++;
		if (counter%1==0){
			f[index_m(a)][index_n(b)]++;//fx<<a<<std::endl;
		}
		if (counter==energy_update_cycle){
			//fx<<std::endl;
			this->updateall();
			counter=0;
			return this->HHR;
		}
		return -999;
	}

	inline double getBiasingE(double a, double b){
		int m=index_m(a);
		int n=index_n(b);

		return -G[m][n];
		/*if (n==0 || n==binnum+1)
			return -G[n];
		double G11=G[n-1],G0=G[n],G1=G[n+1];
		double a0=a-(binstart+binsize*(n-1));
		double p=a0/binsize,out;
		if (p<0.5){
			out=(0.5-p)*G11+(0.5+p)*G0;
		}
		else{
			out=(p-0.5)*G1+(1.5-p)*G0;
		}
		return -out;*/
	}

	int pickle(char* filename){
		std::ofstream fp(filename);
		if (!fp.good()){
			std::cout<<"BiasingPotential::pickle could not write to file:"<<filename<<std::endl;
			exit(3);
		}
		char buf[400];
		sprintf(buf,"%4s %4s %3s %8s %12s %12s %14s %10s","#m","#n","SYM","binstart","G","K","HS","f");
		fp<<buf<<std::endl;
		for (int j=0; j<=binnum_n+1;j++){
			for (int i=0; i<=binnum_m+1; i++){
				sprintf(buf,"%4d %4d %3s %8.3f %8.3f %12g %12g %14d %10d",
					i,j,i==0?"NEG":(i==binnum_m+1?"POS":"---"),
					binstart_m+(i-1)*binsize_m, binstart_n+(j-1)*binsize_n, 
					G[i][j],K[i][j],HS[i][j],f[i][j]);
				fp<<buf<<std::endl;
			}
			fp<<"---Column "<<j<<" end---"<<std::endl;
		}
		fp <<"===Table end==="<<std::endl;
		fp <<"binstart_m "<<binstart_m<<std::endl;
		fp <<"binend_m "<<binend_m<<std::endl;
		fp <<"binnum_m "<<binnum_m<<std::endl;

		fp <<"binstart_n "<<binstart_n<<std::endl;
		fp <<"binend_n "<<binend_n<<std::endl;
		fp <<"binnum_n "<<binnum_n<<std::endl;
		
		fp <<"counter "<<counter<<std::endl;
		fp <<"HHR "<<HHR<<std::endl;



		fp.close();
		return 0;
	}

	int load(char* filename){
		std::ifstream fp(filename);
		if (!fp.good()){
			std::cout<<"BiasingPotential::load could not load file:"<<filename<<std::endl;
			exit(3);
		}
		char buf[400];
		fp.getline(buf,400);

		int i=0,j=0;
		while(!fp.eof()){
			fp.getline(buf,400);
			if (strcmp(buf, "===Table end===")==0) break;
			std::stringstream s(buf);
			s>>i>>j>>buf>>buf;
			if (s.fail()) {
				std::cout<<"Failed to read line:"<<s.str()<<std::endl;
				continue;
			}
			s>>G[i][j]>>K[i][j]>>HS[i][j]>>f[i][j];
			std::cout<<i<<' '<<j<<" G= "<<G[i][j]<<std::endl;
		}

		
		fp.getline(buf,400);
		std::stringstream(buf)>>buf>>binstart_m;
		fp.getline(buf,400);
		std::stringstream(buf)>>buf>>binend_m;
		fp.getline(buf,400);
		std::stringstream(buf)>>buf>>binnum_m;
		
		fp.getline(buf,400);
		std::stringstream(buf)>>buf>>binstart_n;
		fp.getline(buf,400);
		std::stringstream(buf)>>buf>>binend_n;
		fp.getline(buf,400);
		std::stringstream(buf)>>buf>>binnum_n;

		fp.getline(buf,400);
		std::stringstream(buf)>>buf>>counter;
		fp.getline(buf,400);
		std::stringstream(buf)>>buf>>HHR;

		binsize_m=(binend_m-binstart_m)/binnum_m;
		binsize_n=(binend_n-binstart_n)/binnum_n;

		fp.close();

		return 0;
	}
};

#if 0
class BiasingPotentialE{
public:
	static const int MAXBIN=500;
	
	double E[MAXBIN+2],K[MAXBIN+2];
	unsigned long HS[MAXBIN+2],f[MAXBIN+2];
	
	double binstart,binend;
	int binnum;
	unsigned long counter;

	double HHR;

private:
	double binsize;
	BiasingPotentialE();
	std::ofstream fE;
	std::ofstream fK;
	std::ofstream fHS;
	std::ofstream ff;
	std::ofstream fx;
	inline int index(double a){
		if(a<binstart)
			return 0;
		else if (a>=binend)
			return binnum+1;
		else{
			return int((a-binstart)/binsize)+1;
		}
	}

	int updateE(){
		long sumHS=0;
		for (int i=0; i<=binnum+1; i++) sumHS+=HS[i];
		E[0]=E[binnum+1]=-9999;
		for (int i=1; i<=binnum; i++) {
			E[i]=(E[i]*(HS[i]-f[i])+K[i]*1)/(HS[i]+0*f[i]);
			//if (f[i]==0) E[i]*=0.95;
			K[i]=0;
		}
		return 0;
	}

	double updateRepository_HHRc(){
		unsigned long HSmin=HS[1],HSmax=HS[1];
		for (int i=0; i<=binnum+1;i++){
			//ff<<f[i]<<' ';
		}
		//ff<<std::endl;

		for (int i=0; i<=binnum+1;i++){
			HS[i]+=f[i];
			f[i]=0;
			if (i>=2 && i<=binnum){
				if (HS[i]<HSmin) HSmin=HS[i];
				else if (HS[i]>HSmax) HSmax=HS[i];
			}
		}

		double HHRc = (HSmin==0)?1e5:double(HSmax)/double(HSmin);
		std::cout<<"HHRc="<<HHRc<<std::endl;
		return HHRc;
	}

	double updateall(){
		for (double x=0;x<=55;x+=0.1){
			fE<<-getBiasingE(x)<<' ';
		}
		fE<<std::endl;

		for (int i=0; i<=binnum+1;i++){
			
			fHS<<HS[i]<<' ';
			fK<<K[i]<<' ';
		}

		fHS<<std::endl;
		fK<<std::endl;
		this->HHR = this->updateRepository_HHRc();
		this->updateE();


		return this->HHR;
	}

public:
	BiasingPotentialE(double rbinstart, double rbinend, double rbinnum):
	binstart(rbinstart),binend(rbinend),binnum(rbinnum),HHR(1e5),counter(0),
	fE("fE.txt"), fK("fK.txt"), fHS("fHS.txt"), ff("ff.txt"),fx("fx.txt"){
		if (rbinnum>MAXBIN){
			std::cout<<"Exceed the maximum number of bins"<<std::endl;
			exit(3);
		}

		for (int i=0;i<=binnum+1;i++) {
			E[i]=0;K[i]=0;
			HS[i]=1;f[i]=0;
		}
		binsize=(binend-binstart)/binnum;
	}

	~BiasingPotentialE(){
		fE.close();fK.close();fHS.close();ff.close();fx.close();
	}

	inline double collect(double a, double in_E){
		counter++;
		if (counter%1==0){
			int id=index(a);
			f[id]++;fx<<in_E<<' ';
			K[id]=K[id]+in_E;
		}
		if (counter==10000){
			fx<<std::endl;
			this->updateall();
			counter=0;
			return this->HHR;
		}
		return -999;
	}

	inline double getBiasingE(double a){
		int n=index(a);
		return -E[n];
		/*if (n==0 || n==binnum+1)
			return -E[n];
		double E11=E[n-1],E0=E[n],E1=E[n+1];
		double a0=a-(binstart+binsize*(n-1));
		double p=a0/binsize,out;
		if (p<0.5){
			out=(0.5-p)*E11+(0.5+p)*E0;
		}
		else{
			out=(p-0.5)*E1+(1.5-p)*E0;
		}
		return -out;*/
	}

	int pickle(char* filename){
		std::ofstream fp(filename);
		if (!fp.good()){
			std::cout<<"BiasingPotential::pickle could not write to file:"<<filename<<std::endl;
			exit(3);
		}
		char buf[400];
		sprintf(buf,"%4s %3s %8s %12s %12s %14s %10s","#","SYM","binstart","E","K","HS","f");
		fp<<buf<<std::endl;
		for (int i=0; i<=binnum+1; i++){
			sprintf(buf,"%4d %3s %8.3f %12g %12g %14d %10d",
				i,i==0?"NEG":(i==binnum+1?"POS":"---"),
				binstart+(i-1)*binsize,E[i],K[i],HS[i],f[i]);
			fp<<buf<<std::endl;
		}
		fp <<"binstart "<<binstart<<std::endl;
		fp <<"binend "<<binend<<std::endl;
		fp <<"binnum "<<binnum<<std::endl;
		fp <<"counter "<<counter<<std::endl;
		fp <<"HHR "<<HHR<<std::endl;

		fp.close();
		return 0;
	}

	int load(char* filename){
		std::ifstream fp(filename);
		if (!fp.good()){
			std::cout<<"BiasingPotential::load could not load file:"<<filename<<std::endl;
			exit(3);
		}
		char buf[400];
		fp.getline(buf,400);

		int i=0;
		while(!fp.eof()){
			fp.getline(buf,400);
			std::stringstream s(buf);
			s>>i>>buf>>buf;
			s>>E[i]>>K[i]>>HS[i]>>f[i];
			std::cout<<i<<" E= "<<E[i]<<std::endl;
		}

		
		fp.getline(buf,400);
		std::stringstream(buf)>>buf>>binstart;
		fp.getline(buf,400);
		std::stringstream(buf)>>buf>>binend;
		fp.getline(buf,400);
		std::stringstream(buf)>>buf>>binnum;
		fp.getline(buf,400);
		std::stringstream(buf)>>buf>>counter;
		fp.getline(buf,400);
		std::stringstream(buf)>>buf>>HHR;

		binsize=(binend-binstart)/binnum;
		fp.close();
	}
};
#endif


#endif /* BIAS_H */