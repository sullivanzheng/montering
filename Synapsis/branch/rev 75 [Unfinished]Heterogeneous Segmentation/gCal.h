//Simpson rule, adaptive
#ifndef CALCULUS_H
#define CALCULUS_H

#include <cmath>
#include <iostream>

class gCal{
public:
	double re_a_simp(double (gCal::*integrand)(double), double a, double m,double b,
					 double fa, double fm, double fb, double err){
		using namespace std;
		static int level=0;
		//cout<<"LEVEL:"<<level<<" ["<<a<<","<<b<<"] Err="<<err<<endl;
		double Q1=(b-a)/6*(fa+4*fm+fb);
		double ma=(m+a)/2,mb=(m+b)/2;
		double fma=(*this.*integrand)(ma),fmb=(*this.*integrand)(mb);
		double Q2=(b-a)/12*(fa+4*fma+2*fm+4*fmb+fb);
		//cout<<"Q1="<<Q1<<" Q2="<<Q2<<"fabs(Q2-Q1)/15="<<fabs(Q2-Q1)/15<<endl;
 		if (fabs(Q2-Q1)/15<err && level>3){
			return Q2;
		}
		else{
			++level;
			double temp=re_a_simp(integrand,a,ma,m,fa,fma,fm,err<1e-15?err:err/2)
				+re_a_simp(integrand,m,mb,b,fm,fmb,fb,err<1e-15?err:err/2);
			--level;
			return temp;
		}
	}

	double a_simp(double (gCal::*integrand)(double), double a, double b,double err){
		double m=(a+b)/2,ma=(m+a)/2,mb=(m+b)/2;
		double fa=(*this.*integrand)(a);
		double fb=(*this.*integrand)(b);
		double fm=(*this.*integrand)(m);
		double fma=(*this.*integrand)(ma);
		double fmb=(*this.*integrand)(mb);
		double Q1=(b-a)/6*(fa+4*fm+fb);
		double Q2=(b-a)/12*(fa+4*fma+2*fm+4*fmb+fb);
		if (fabs(Q2-Q1)/15<err){
			return Q2;
		}
		else{
			return re_a_simp(integrand,a,ma,m,fa,fma,fm,err<1e-15?err:err/2)
				+re_a_simp(integrand,m,mb,b,fm,fmb,fb,err<1e-15?err:err/2);
		}
	}

	double bisect(double (gCal::*f)(double),double x0, double x1,double err) {
		using namespace std;
		double f0=(*this.*f)(x0),f1=(*this.*f)(x1);
		double fm,xm;
		while(fabs(x0-x1)>err){
			xm=(x0+x1)/2;
			fm=(*this.*f)(xm);
			double t=fm*f0;
			if (t>0)
				x0=xm;
			else if(t==0)
				return xm;
			else 
				x1=xm;
		}
		return (x0+x1)/2;
	}

	double g;
	double m;

	double numerator(double x){
		return cos(x)*sin(x)*exp(-g*x*x);
	}

	double denominator(double x){
		return sin(x)*exp(-g*x*x);
	}

	double MeanCosMinusM(double _g){
		static const double pi=3.14159265358979;
		this->g=_g;
		double a=a_simp(&gCal::numerator,0,0.1,1e-13)+a_simp(&gCal::numerator,0.1,pi,1e-13);
		double b=a_simp(&gCal::denominator,0,0.1,1e-13)+a_simp(&gCal::denominator,0.1,pi,1e-13);
		return a/b-this->m;
	}

public:
	double FindGbyK(double k){ //k is number of seg per kuhn length.
		this->m=(k-1.)/(k+1.);
		return bisect(&gCal::MeanCosMinusM,0.0001,200,1e-7);
	}
/*
private:
	void _printMeanCos(){		
		//A test function
		static const double pi=3.14159265358979;
		using namespace std;
		ofstream fp("Meancos.txt");
		for (double s=1;s<=282;s+=1){
			cout<<s<<" "<<this->FindGbyK(s)<<endl;
//			this->g=s;
//			cout<<this->g<<endl;
//			double a=a_simp(&gCal::numerator,0,0.1,1e-13)+a_simp(&gCal::numerator,0.1,pi,1e-13);
//			double b=a_simp(&gCal::denominator,0,0.1,1e-13)+a_simp(&gCal::denominator,0.1,pi,1e-13);
//			fp<<s<<"  "<<a/b<<endl;
			
		}
		fp.close();
	}*/
};


#endif /*CALCULUS_H*/