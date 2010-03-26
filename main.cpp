
#include <string> //used in map<string, string> and config reading.
#include <sstream>//string to number convertion
#include <cstring>//using itoa.
#include "MCbox.h"
#include "file\configldr.h"
using namespace std;

void calBp(int segNumNow,ofstream &mainlog,
           double r_theta_k,
           double r_h,
           double minE2EDis,//in bp
           double minE2EAng,//in deg
           unsigned long tot_moves,
           int PNum,//Number of conditional probabilities.
           unsigned long seed)
{
    const int FILENUM=7;
    char prefix[FILENUM][100];
    char buf[80];
    for (int pre=0;pre<7;pre++)
        sprintf(prefix[pre],"%d-%d",segNumNow,pre);
    sprintf(buf,"P_bp_%d.txt",segNumNow);
    ofstream fp(buf);
    const int max_simu_segments=25;
    double paceAng=pow(minE2EAng/180.0,1.0/PNum),
        paceDist=pow(minE2EDis/(segNumNow*bpperseg),1.0/PNum);
    double P_r0_r1[max_simu_segments];
    double PkCondInV0[MAXKINKNUM]={0};
    unsigned long KinkCountInV0[MAXKINKNUM]={0};
    statqueue <double> Jfinal;

    for (int circleNum=0;circleNum<FILENUM;circleNum++){
        MCbox_linear
            metro(prefix[circleNum],segNumNow,r_theta_k,r_h,seed);
        metro.dnaChain.dispChainCord();
        for (int j=PNum-1;j>=0;j--){
            bool kinkLog=(j==PNum-1)?true:false;
            P_r0_r1[PNum-1-j]=metro.performSimpleCondP
                (
                (j==PNum-1)?tot_moves/PNum*10:tot_moves/PNum, //more sampling for ligation conformation
                segNumNow*pow(paceDist,j), 
                180.0*pow(paceAng,j),
                segNumNow*pow(paceDist,(j+1)),
                180.0*pow(paceAng,(j+1)),
                kinkLog
                );
            cout <<'.';
#if ALLOW_ADDITIONAL_RUNS
            //If counts is too small, then run more times to improve accuracy.
            double counts=P_r0_r1[simu_segments-1-j]*(tot_moves/(double)simu_segments);
            if (counts<3000){
                cout <<"Counts too small. Increase number of moves."<<endl;
                int more=3000/counts+1;
                more=(more>50?50:more);
                cout <<"Make "<<more<<" more runs."<<endl;
                double temp_p;
                for (int i_more=0;i_more<more;i_more++){
                    temp_p=metro.performMetropolisLinearLigationCondP
                        (
                        tot_moves/simu_segments,
                        segNumNow*pow(pace,j),
                        180.0*pow(pace,j),
                        segNumNow*pow(pace,(j+1)),
                        180.0*pow(pace,(j+1)),
                        kinkLog
                        );
                }
                P_r0_r1[simu_segments-1-j]=temp_p;
            }
#endif
            if (j==PNum-1){
                //metro.logAngleDist("_ligation_kinknum.txt");
                metro.getAngleDist(KinkCountInV0);
                cout<<"\n====LIGATION KINK COUNT AND PROBABILITY=====\n"<<endl;
                fp<<"\n====LIGATION KINK COUNT AND PROBABILITY=====\n"<<endl;
                for (int i=0;i<MAXKINKNUM;i++){
                    if (KinkCountInV0[i]>0){
                        PkCondInV0[i]=double(KinkCountInV0[i])
                            /double(metro.dnaChain.stats.ligateAccepts);
                        cout<<i <<"kink(s) - " <<KinkCountInV0[i]<<"\t"<<PkCondInV0[i]<<endl;
                        fp  <<i <<"kink(s) - " <<KinkCountInV0[i]<<"\t"<<PkCondInV0[i]<<endl;
                    }
                }
            }
            metro.logAcceptsEndToEndDistanceAngle();
            metro.dnaChain.stats.resetStat();
            metro.clearAngleStats();
            //sprintf(buf,"snap%d_%d.txt",circleNum,j);
        }

        double P_r0_rN=1.0;
        for (int i=0;i<PNum;i++){
            P_r0_rN=P_r0_rN*P_r0_r1[i];
        }
        double PGoodTorsionCondInV0=0;
        for (int i=0;i<MAXKINKNUM;i++){
            PGoodTorsionCondInV0+=PkCondInV0[i]*P_t_over_2t(totsegnum*bpperseg,i);
        }
        double coef=10.0/6.022141*1.5
            /(pow(0.34*minE2EDis,3)*PI)
            /(1-cos(PI*minE2EAng/180.0));
        fp<<"Minicircle"<<totsegnum*bpperseg<<endl
            <<" P(V0|V)(NoTorsional):"<<P_r0_rN<<endl
            <<" P(CorrectTorsion|V0):"<<PGoodTorsionCondInV0<<endl
            <<" Jfactor(no torsional):"<<P_r0_rN*coef<<endl
            <<" noTorsional coefficient:"<<coef<<endl
            <<"bp Jfactor:  "<<P_r0_rN*PGoodTorsionCondInV0*coef<<endl;
        cout<<"Minicircle"
            <<totsegnum*bpperseg<<"bp Jfactor:  "
            <<P_r0_rN*PGoodTorsionCondInV0*coef<<endl;
        Jfinal.push(P_r0_rN*PGoodTorsionCondInV0*coef);
    }
    fp<<"=========================="<<endl
        <<"Jfinal: "<<Jfinal.getMean()<<endl
        <<"StdevMean: "<<Jfinal.getStdevOfMean()
        <<'('<<Jfinal.getStdevOfMean()/Jfinal.getMean()<<')'<<endl;
    fp.close();
    sprintf(buf,"%d seg: Jfinal%12.6e +/- %12.6e(%4.2f)\n",segNumNow,Jfinal.getMean(),
        Jfinal.getStdevOfMean(),Jfinal.getStdevOfMean()/Jfinal.getMean());
    mainlog<<buf;
}
void main_do(){
	unsigned long seed,tot_moves;
	double r_theta_k,r_h,mine2edis,mine2eang;
    int startLength,endLength,segNum,bpInterval;
	dict config=readconfig("config.txt");
	stringstream(config["RANDSEED"]) >>seed;
	stringstream(config["THETA_K"])>>r_theta_k;
	stringstream(config["H"])>>r_h;
	stringstream(config["g"])>>g;//global
	stringstream(config["DELTA_TW_K"])>>DELTA_TW_K;
	stringstream(config["BPPERSEG"])>>bpperseg;//base pair number per rigid segment.
    stringstream(config["STARTSEGNUM"])>>startLength;
    stringstream(config["BPINTERVAL"])>>bpInterval;
    stringstream(config["ENDSEGNUM"])>>endLength;
    stringstream(config["TOTMOVES"])>>tot_moves;
    stringstream(config["MINE2EDIS"])>>mine2edis;
    stringstream(config["MINE2EANG"])>>mine2eang;
    stringstream(config["SEGNUM"])>>segNum;
    ofstream mfp("mainlog.txt");
    for (int chainLength=startLength;chainLength<=endLength;chainLength+=bpInterval){
        calBp(chainLength,mfp,r_theta_k,r_h,mine2edis,mine2eang,tot_moves,segNum,seed);
    }
    mfp.close();
}

void unittestEndToEnd() {
    //test sample
    g=4.909912581;
    ofstream fp("unittest.txt");
    calBp(20,fp,0,0,20,180,10000000,1,165);
}

void unittestP_t_over_2t(){
    ofstream fp("P.txt");
    for (int i=30;i<1000;i++){
        fp<<i<<"  "<<P_t_over_2t(i,0)<<endl;
    }
    fp.close();
}

void main(){
   main_do();
    //unittestP_t_over_2t();
    //unittestEndToEnd();
}
