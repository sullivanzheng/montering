#include "MCbox.h" 
MCbox_circular::MCbox_circular(
    char const *r_filePrefix, 
    int const length, 
    double const r_theta_k, 
    double const r_h, 
    unsigned long r_seeding)
    :fp_log(strcat_noOW(buf, strBufSize,const_cast<char*>(r_filePrefix), "_log.txt")), 
    dnaChain(strcat_noOW(buf, strBufSize, const_cast<char*>(r_filePrefix), ".vmc"),length), 
    seeding(r_seeding)
{
    strcpy(filePrefix, r_filePrefix);
    //Initialize the global variables related to 
    //the properties of chain.
    h = r_h;
    theta_k = r_theta_k * PI / 180.0;
    totsegnum = length;
    crank_max_length = length / 2;
    //Initialize statistical variables.
    for (int i = 0; i < 180; i++)
        anglenums[i] = 0;
    for (int i = 0; i < MAXKINKNUM + 1; i++)
        N_kink[i] = 0;
    //Log file handle.
    //dnaChain.snapshot(strcat_noOW(buf, strBufSize, filePrefix, "_ini.txt"));
    logParameters();
}
void MCbox_circular::logAngleDist(char *suffix)
{
    strcat_noOW(buf, strBufSize, filePrefix, "_angleDist.txt");
    ofstream fp (strcat_noOW(buf, strBufSize, buf, suffix));
    fp << "The number of n-kink conformations" << endl;
    fp << 'n' << '\t' << "Count" << endl;
    for (int i = 0; i < MAXKINKNUM; i++)
        fp << i << '\t' << N_kink[i] << endl;
    fp << "Angle distribution" << endl;
    for (int i = 0; i < 180; i++)
        fp << i << '\t' << anglenums[i] << endl;
    fp << endl;
    fp.close();
}
void MCbox_circular::clearAngleStats(void)
{
    for (int i = 0; i < 180; i++)
        anglenums[i] = 0;
    for (int i = 0; i < MAXKINKNUM + 1; i++)
        N_kink[i] = 0;
}
void MCbox_circular::appendAngleStats(void)
{
    unsigned short int kink_counter = 0;
    for (int i = 0; i < maxnum; i++)
    {
        anglenums[int(floor(dnaChain.C[i].bangle * 180 / PI))]++;
        //stat angle distribution
        if (dnaChain.C[i].bangle > KINKLOWERBOUND)
            kink_counter++;
    }
    if (kink_counter > MAXKINKNUM)
    {
        cout << endl << "In step: " << dnaChain.stats.moves() << "Number of kink is larger than " << MAXKINKNUM << endl;
        cout << endl << "Simulation is terminated" << endl;
        getchar();
        exit(EXIT_FAILURE);
    }
    else
        N_kink[kink_counter]++;
}
void MCbox_circular::logAccepts(void)
{
    sprintf(buf, "%12d %12d ", 
        dnaChain.stats.moves.getTotCounts(), 
        dnaChain.stats.accepts.getNumber());
    fp_log << buf;
}
void MCbox_circular::performMetropolisCircularCrankOnly(long monte_step)
{
    MTRand53 mt(seeding);
    for (int moves = 1; moves <= monte_step; moves++)
    {
        //MAKE MOVES
        double rotAng;
        double selection;
        selection = drand(1.0);
        if (selection > P_SMALLROTATION)
        {
            rotAng = drand(-maxRotAng, maxRotAng);
            if (rotAng <= 0)
                rotAng -= maxRotAng;
            else
                rotAng += maxRotAng;
        }
        else
            rotAng = drand(-maxRotAng, maxRotAng);
        int m;
        int n;
        m = irand(maxnum + 1);
        n = (m + (irand(crank_max_length - crank_min_length) + crank_min_length))%totsegnum;
        double dE;
        dE = dnaChain.deltaE_TrialCrankshaft(m, n, rotAng);
        if (dE < 0)
            dnaChain.crankshaft(m, n, rotAng);
        else
        {
            double exp_E;
            double r;
            exp_E = exp(-dE);
            r = mt();
            if (r < exp_E)
                dnaChain.crankshaft(m, n, rotAng);
            //if (r<exp_E && exp_E<0.006)
            //	cout<<"M "<<m<<" N "<<n<<endl;
        }
    }
}
void MCbox_circular::logParameters(void){
		fp_log 
			<<"========================PARAMETERS========================"<<endl
			<<">FilePrefix and path	= "<<filePrefix<<endl
			<<" maxa	= "<<	maxa	<<endl
			<<" MAXKINKNUM	= "<<	MAXKINKNUM	<<endl
			<<" KINKLOWERBOUND	="<<	KINKLOWERBOUND	<<endl
			<<">ENERGY CURVE: theta_k="<<theta_k*180.0/PI<<" h="<<	h	<<endl
			<<" crank_min_length	= "<<	crank_min_length	<<endl
			<<">maxRotAng	= "<<	maxRotAng	<<endl
			<<" maxnum	= "<<	maxnum	<<endl
			<<">totsegnum * bpperseg = total bp \t"<<totsegnum<<'*'<<bpperseg<<'='<<totsegnum*bpperseg<<endl
			<<" crank_max_length	= "<<	crank_max_length	<<endl
			<<" P_SMALLROTATION	="<< P_SMALLROTATION<<endl
			<<" DELTA_TW_K =	"<< DELTA_TW_K << "//In number of turns"<<endl
			<<" Chain Lk ="<<dnaChain.Lk<<endl
			<<"================ENERGY CURVE PARAS===================="<<endl
			<<" theta_k "<<theta_k<<" H "<<h<<" g "<<g<<endl
			<<endl;
	}




//////////////////////////////////////////////////////////////////////////////////////////


MCbox_linear::MCbox_linear(
    char const *r_filePrefix, 
    int const length, 
    double const r_theta_k, 
    double const r_h, 
    unsigned long r_seeding,
    double r_crankratio,
    double r_rotpercentage)
    :fp_log(strcat_noOW(buf, strBufSize,const_cast<char*>(r_filePrefix), "_log.txt")), 
    //dnaChain(strcat_noOW(buf, strBufSize, const_cast<char*>(r_filePrefix), ".vmc"),length), 
    dnaChain(length), 
    seeding(r_seeding),
    crankratio(r_crankratio),
    rotpercentage(r_rotpercentage)
{
    strcpy(filePrefix, r_filePrefix);
    //Initialize the global variables related to 
    //the properties of chain.
    h = r_h;
    theta_k = r_theta_k * PI / 180.0;
    maxnum = length - 1;
    totsegnum = length;
    crank_max_length = length / 2;
    //Initialize statistical variables.
    for (int i = 0; i < 180; i++)
        anglenums[i] = 0;
    for (int i = 0; i < MAXKINKNUM + 1; i++)
        N_kink[i] = 0;
    //Log file handle.
    //dnaChain.snapshot(strcat_noOW(buf, strBufSize, filePrefix, "_ini.txt"));
    logParameters();
}
void MCbox_linear::logAngleDist(char *suffix)
{
    strcat_noOW(buf, strBufSize, filePrefix, "_angleDist.txt");
    ofstream fp (strcat_noOW(buf, strBufSize, buf, suffix));
    fp << "The number of n-kink conformations" << endl;
    fp << 'n' << '\t' << "Count" << endl;
    for (int i = 0; i < MAXKINKNUM; i++)
        fp << i << '\t' << N_kink[i] << endl;
    fp << "Angle distribution" << endl;
    for (int i = 0; i < 180; i++)
        fp << i << '\t' << anglenums[i] << endl;
    fp << endl;
    fp.close();
}
void MCbox_linear::clearAngleStats(void)
{
    for (int i = 0; i < 180; i++)
        anglenums[i] = 0;
    for (int i = 0; i < MAXKINKNUM + 1; i++)
        N_kink[i] = 0;
}
void MCbox_linear::appendAngleStats(void)
{
    unsigned short int kink_counter = 0;
    for (int i = 0; i < maxnum; i++)
    {
        anglenums[int(floor(dnaChain.C[i].bangle * 180 / PI))]++;
        //stat angle distribution
        if (dnaChain.C[i].bangle > KINKLOWERBOUND)
            kink_counter++;
    }
    if (kink_counter > MAXKINKNUM)
    {
        cout << endl << "In step: " << dnaChain.stats.moves() << "Number of kink is larger than " << MAXKINKNUM << endl;
        cout << endl << "Simulation is terminated" << endl;
        getchar();
        exit(EXIT_FAILURE);
    }
    else
        N_kink[kink_counter]++;
}
void MCbox_linear::logAccepts(void)
{
    sprintf(buf, "%12d %12d \n", 
        dnaChain.stats.moves.getTotCounts(), 
        dnaChain.stats.accepts.getNumber());
    fp_log << buf;
}

void MCbox_linear::logAcceptsEndToEndDistance(void)
{
    sprintf(buf, "%12d %12d %12.6f+/-%12.6f\n", 
        dnaChain.stats.moves.getTotCounts(), 
        dnaChain.stats.accepts.getNumber(),
        dnaChain.stats.endToEndDistance2.getMean(),
        dnaChain.stats.endToEndDistance2.getStdev());
    fp_log << buf;
    fp_log.flush();
}
void MCbox_linear::logAcceptsEndToEndDistanceAngle(void)
{   
    sprintf(buf, "%12d ACC%12d LIG%12d P=%9.5f" //Probility
        "HC%12d CSH%12d " //Move stats
        "E2eDis:%12.6f+/-%12.6f[%12.6f,%12.6f] Ang:%12.6f+/-%12.6f[%12.6f,%12.6f]\n", //Move stats 2.
        dnaChain.stats.moves.getTotCounts(), 
        dnaChain.stats.accepts.getNumber(),
        dnaChain.stats.ligateAccepts,
        double(dnaChain.stats.ligateAccepts)/dnaChain.stats.moves.getNumber(),
        dnaChain.stats.halfChainSuccess,
        dnaChain.stats.crankSuccess,
        dnaChain.stats.endToEndDistance2.getMean(),
        dnaChain.stats.endToEndDistance2.getStdev(),
        dnaChain.stats.endToEndDistance2.getminItem(),
        dnaChain.stats.endToEndDistance2.getmaxItem(),
        dnaChain.stats.endToEndAngle.getMean()/PI*180.0,
        dnaChain.stats.endToEndAngle.getStdev()/PI*180.0,
        dnaChain.stats.endToEndAngle.getminItem()/PI*180.0,
        dnaChain.stats.endToEndAngle.getmaxItem()/PI*180.0);
    fp_log << buf;
    fp_log.flush();
}

void MCbox_linear::performMetropolisLinearCrankOnly(long monte_step)
{
    MTRand53 mt(seeding);
    for (int mt_moves = 1; mt_moves <= monte_step; mt_moves++)
    {
        //MAKE MOVES
        dnaChain.countmove();
        double rotAng;
        double selection;
        selection = drand(1.0);
        if (selection > P_SMALLROTATION)
        {
            rotAng = drand(-maxRotAng, maxRotAng);
            if (rotAng <= 0)
                rotAng -= maxRotAng;
            else
                rotAng += maxRotAng;
        }
        else
            rotAng = drand(-maxRotAng, maxRotAng);
        int m;
        int n;
        m = irand(maxnum + 1);
        n = (m + (irand(crank_max_length - crank_min_length) + crank_min_length))%totsegnum;
        double dE;
        dE = dnaChain.deltaE_TrialCrankshaft(m, n, rotAng);
        if (dE < 0)
            dnaChain.crankshaft(m, n, rotAng);
        else
        {
            double exp_E;
            double r;
            exp_E = exp(-dE);
            r = mt();
            if (r < exp_E)
                dnaChain.crankshaft(m, n, rotAng);
            //if (r<exp_E && exp_E<0.006)
            //	cout<<"M "<<m<<" N "<<n<<endl;
        }
    }
}
void MCbox_linear::logParameters(void){
		fp_log 
			<<"========================PARAMETERS========================"<<endl
            <<"RANDOM SEED \t"<<seeding<<endl
			<<">FilePrefix and path	= "<<filePrefix<<endl
			<<" maxa	= "<<	maxa	<<endl
			<<" MAXKINKNUM	= "<<	MAXKINKNUM	<<endl
			<<" KINKLOWERBOUND	="<<	KINKLOWERBOUND	<<endl
			<<">ENERGY CURVE: theta_k="<<theta_k*180.0/PI<<" h="<<	h	<<endl
			<<" crank_min_length	= "<<	crank_min_length	<<endl
			<<">maxRotAng	= "<<	maxRotAng	<<endl
			<<" maxnum	= "<<	maxnum	<<endl
			<<">totsegnum * bpperseg = total bp \t"<<totsegnum<<'*'<<bpperseg<<'='<<totsegnum*bpperseg<<endl
			<<" crank_max_length	= "<<	crank_max_length	<<endl
			<<" P_SMALLROTATION	="<< P_SMALLROTATION<<endl
			<<" DELTA_TW_K =	"<< DELTA_TW_K << "//In number of turns"<<endl
			<<" Chain Lk = 0 //Do not have this feature for a linear chain"<<endl
			<<"================ENERGY CURVE PARAS===================="<<endl
			<<" theta_k "<<theta_k<<" H "<<h<<" g "<<g<<endl
			<<endl;
	}

void MCbox_linear::performMetropolisLinearCrankHalfChain(long monte_step){
    MTRand53 mt(seeding);
//    double x2=0.0,x=0.0;int times=0;
    for (int moves = 1; moves <= monte_step; moves++)
    {
        //MAKE MOVES
        dnaChain.countmove();
        double rotAng;
        double moveselection;
        moveselection=drand(1.0);

#if 0       
        if (moves%100==0){
        double temp=dnaChain.getEndToEndDistance();
        temp=temp*temp;
        x2+=temp*temp;
        x+=temp;
        times++;
        this->fp_log <<temp<<endl;
        }
#endif
        if (moveselection < crankratio){
            cout <<"DAMNED IT!!CRANK!!"<<endl;
            double selection = drand(1.0);
            if (selection > P_SMALLROTATION)
            {
                rotAng = drand(-maxRotAng, maxRotAng);
                if (rotAng > 0)
                    rotAng += maxRotAng;
                else
                    rotAng -= maxRotAng;
            }
            else
                rotAng = drand(-maxRotAng, maxRotAng);
            int m; int n;
            m = irand(maxnum+1);
            while (abs((n=irand(maxnum+1))-m)<crank_min_length);
            //n=(m+irand(crank_max_length-crank_min_length)+crank_min_length)%totsegnum;
            //fp_log<<m<<'\t'<<n<<endl;
            if (n<m){
                int mid;
                mid=n;n=m;m=mid;
            }
            double dE;
            dE = dnaChain.deltaE_TrialCrankshaft(m, n, rotAng);
#if 0 //Energy correctness check
            double E1,E2,dE2;
            E1=dnaChain.calG_bSum();
            dnaChain.crankshaft(m,n,rotAng);
            E2=dnaChain.calG_bSum();
            dE2=E2-E1;
            dnaChain.crankshaft(m,n,-rotAng);
            if (dE!=0 && abs(dE2-dE)>1e-13){
                cout<<"Energy Inaccurate! Alert!"<<endl;
            }
#endif
            if (dE < 0)
                dnaChain.crankshaft(m, n, rotAng);
            else
            {
                double exp_E=exp(-dE);
                double r=mt();
                if (r < exp_E)
                    dnaChain.crankshaft(m, n, rotAng);
                //if (r<exp_E && exp_E<0.006)
                //	cout<<"M "<<m<<" N "<<n<<endl;
            }
        }
        else{
            double phi=drand(-1.0,1.0),theta=drand(0.0,2*PI);
            double rv[3]={
                sqrt(1-phi*phi)*cos(theta),
                sqrt(1-phi*phi)*sin(theta),
                phi};
            double chi=drand(-PI*rotpercentage,PI*rotpercentage);
            int m=irand(0,maxnum+1);
            double dE=dnaChain.deltaE_TrialHalfChain(m,rv,chi);

#if 0
            double E1,E2,dE2;
            E1=dnaChain.calG_bSum();
            dnaChain.halfChain(m,rv,chi);
            E2=dnaChain.calG_bSum();
            dE2=E2-E1;
            dnaChain.halfChain(m,rv,-chi);
            if (dE!=0 && abs(dE2-dE)>1e-13){
                cout<<"Energy Inaccurate! Alert!"<<endl;
            }
#endif
            if (dE < 0)
                dnaChain.halfChain(m,rv,chi);
            else
            {
                double exp_E=exp(-dE);
                double r=mt();
                if (r<exp_E)
                    dnaChain.halfChain(m,rv,chi);
            }
        }
    }
/*    cout <<"-------->e2eDist:  "<<x/times
        <<" stdev:"<<sqrt((x2-x*x/double(times))/times)
        <<" datanum: "<<times<<endl;*/
}

void MCbox_linear::performMetropolisLinearTryLigation
        (long monte_step,
        double endToEndDistanceThreshold,
        double r_endToEndAngleThresholdDeg)
{
    MTRand53 mt(seeding);
    double endToEndAngleThreshold;
    endToEndAngleThreshold=r_endToEndAngleThresholdDeg/180*PI;
    for (int moves = 1; moves <= monte_step; moves++)
    {
        dnaChain.countmove();
        double rotAng;
        double moveselection;
        moveselection=drand(1.0);
        if (moveselection < crankratio){
            double selection = drand(1.0);
            if (selection > P_SMALLROTATION)
            {
                rotAng = (mt()*2-1)*maxRotAng;
                if (rotAng > 0)
                    rotAng += maxRotAng;
                else
                    rotAng -= maxRotAng;
            }
            else
                rotAng = (mt()*2-1)*maxRotAng;
            int m; int n;
            m = irand(maxnum-1)+1;//////PAY ATTENTION
                                //////NOT IRAND(MAXNUM+1);
                                //////PIVOT 1~MAXNUM-1
            while (abs((n=irand(maxnum-1)+1)-m)<crank_min_length);
            if (n<m){
                int mid;
                mid=n;n=m;m=mid;
            }
            double dE;
            dE = dnaChain.deltaE_TrialCrankshaft(m, n, rotAng);
            if (dE < 0)
                dnaChain.crankshaft(m, n, rotAng);
            else
            {
                double exp_E=exp(-dE);
                double r=mt();
                if (r < exp_E)
                    dnaChain.crankshaft(m, n, rotAng);
            }
        }
        else{
            double phi=mt()*2-1,theta=mt()*2*PI;//(0.0,2*PI);
            double rv[3]={
                sqrt(1-phi*phi)*cos(theta),
                sqrt(1-phi*phi)*sin(theta),
                phi};
            double chi=(mt()*2-1)*PI*rotpercentage;
            int m=irand(0,maxnum+1);
            double dE=dnaChain.deltaE_TrialHalfChain(m,rv,chi);
            bool ligationOK=dnaChain.trialLigateAfterHalfChainOK(
                m,rv,chi,endToEndDistanceThreshold,
                endToEndAngleThreshold);
            if (!ligationOK){
                ;//do nothing.
            }
            else{
                if (dE < 0)
                    dnaChain.halfChain(m,rv,chi);
                else
                {
                    double exp_E=exp(-dE);
                    double r=mt();
                    if (r<exp_E)
                        dnaChain.halfChain(m,rv,chi);
                }
            }
        }
    }
}


double MCbox_linear::performMetropolisLinearLigationCondP
      (long monte_step,
      double endToEndDistanceThreshold,
      double r_endToEndAngleThresholdDeg,
      double endToEndDistanceThreshold_ligate,
      double r_endToEndAngleThresholdDeg_ligate,
      bool conformationBinningOn)
{
    MTRand53 mt(seeding);
    double endToEndAngleThreshold=
        r_endToEndAngleThresholdDeg/180*PI,
           endToEndAngleThreshold_ligate=
        r_endToEndAngleThresholdDeg_ligate/180*PI;
    dnaChain.stats.ligateAccepts=0;
    //ofstream mclog("big.txt");
    for (int moves = 1; moves <= monte_step; moves++)
    {
        dnaChain.countmove();

        char buf[20];
        sprintf(buf,"%10d)",dnaChain.stats.moves());
        //mclog<<buf;

        double moveselection;
        moveselection=drand(1.0);
        if (moveselection < crankratio){
            double rotAng;
            double selection = drand(1.0);
            if (selection > P_SMALLROTATION)
            {
                rotAng = (mt()*2)-1*maxRotAng;
                if (rotAng > 0)
                    rotAng += maxRotAng;
                else
                    rotAng -= maxRotAng;
            }
            else
                rotAng = (mt()*2-1)*maxRotAng;
            int m; int n;
            // 0<=m<=maxnum+1 
            m = irand(maxnum)+1;
            // 0<=n<=maxnum+1 && abs(n-m)>=crank_min_length
            while (abs((n=irand(maxnum)+1)-m)<crank_min_length);
            if (n<m){
                int mid;
                mid=n;n=m;m=mid;
            }
            //mclog <<'c'<<' '<<m<<','<<n<<','<<rotAng;
            if (dnaChain.trialLigateAfterCrankshaftOK(m,n,rotAng,
                endToEndDistanceThreshold,endToEndAngleThreshold)){
                    //mclog <<"L(+)"<<dnaChain.getEndToEndDistance()<<','<<dnaChain.getEndToEndAngle();
                double dE = dnaChain.deltaE_TrialCrankshaft(m, n, rotAng);
                //mclog <<'E'<<dE<<',';
                if (dE < 0){
                    //mclog <<'+';
                    dnaChain.crankshaft(m, n, rotAng);
                }
                else
                {
                    double exp_E=exp(-dE);
                    double r=mt();
                    if (r < exp_E){
                        //mclog <<'-';
                        dnaChain.crankshaft(m, n, rotAng);
                    }
                    else{
                        //mclog <<'x';
                    }
                }
                
            }
            else {
                //mclog <<"L(-),"<<dnaChain.getEndToEndDistance()<<','<<dnaChain.getEndToEndAngle()<<",X";/*do nothing */
            }
        }
        else{
            double phi=mt()*2-1,theta=mt()*2*PI;//(0.0,2*PI);
            double rv[3]={
                sqrt(1-phi*phi)*cos(theta),
                sqrt(1-phi*phi)*sin(theta),
                phi};
            double chi=(mt()*2-1)*PI*rotpercentage;
            int m=irand(0,maxnum+1);
            //mclog <<'h'<<' '<<m<<','<<theta<<','<<phi<<','<<chi<<',';
            if (dnaChain.trialLigateAfterHalfChainOK(
                m,rv,chi,endToEndDistanceThreshold,
                endToEndAngleThreshold)){
                //mclog <<"L(+),"<<dnaChain.getEndToEndDistance()<<','<<dnaChain.getEndToEndAngle()<<',';
                double dE=dnaChain.deltaE_TrialHalfChain(m,rv,chi);
                //mclog <<'E'<<dE<<',';
                if (dE < 0){
                    //mclog <<'+';
                    dnaChain.halfChain(m,rv,chi);
                }
                else
                {
                    double exp_E=exp(-dE);
                    double r=mt();
                    if (r<exp_E){
                        //mclog <<'-';
                        dnaChain.halfChain(m,rv,chi);
                    }
                    else{
                        //mclog <<'x';
                    }
                }
            }
            else {
                //mclog <<"L(-),"<<dnaChain.getEndToEndDistance()<<','<<dnaChain.getEndToEndAngle()<<",X";
                /* do nothing */
            }
        }

        //Statistics after one movement

        //Record how many conformation are in a smaller volume.
        if (dnaChain.getEndToEndDistance()>endToEndDistanceThreshold_ligate ||
            dnaChain.getEndToEndAngle()>endToEndAngleThreshold_ligate){
            //mclog <<'.';
        }
        else{
            //mclog <<'*';
            dnaChain.stats.ligateAccepts++;
            //Binning the configurations according to the number of disruptions in one conformation.
            if (conformationBinningOn){
                this->appendAngleStats();
            }
        }
        //cout <<dnaChain.stats.moves();
#if 0
        if (1) {
            buf[20];
            sprintf(buf,"s%06d.txt",dnaChain.stats.moves());
            dnaChain.snapshot(buf);
        }
#endif
    //mclog <<endl;


    }//end of movement loop.
    /*char buf[200];
    sprintf(buf,"\nDIS-dis:%6.3f-%6.3f ANG-ang:%6.3f-%6.3f",
        endToEndDistanceThreshold,
        endToEndDistanceThreshold_ligate,
        r_endToEndAngleThresholdDeg,        
        r_endToEndAngleThresholdDeg_ligate);
    this->fp_log<<buf<<endl;*/
    return (double)dnaChain.stats.ligateAccepts/(double)dnaChain.stats.moves();
}

double MCbox_linear::performSimpleCondP(long monte_step,
                              double endToEndDistanceThreshold,
                              double r_endToEndAngleThresholdDeg,
                              double endToEndDistanceThreshold_ligate,
                              double r_endToEndAngleThresholdDeg_ligate,
                              bool conformationBinningOn){
    const bool dbg=false;
    
    ofstream anglelog("anglelog.txt");
    MTRand53 mt(seeding);
    double endToEndAngleThreshold=
        r_endToEndAngleThresholdDeg/180*PI,
           endToEndAngleThreshold_ligate=
        r_endToEndAngleThresholdDeg_ligate/180*PI;
    dnaChain.stats.ligateAccepts=0;
    for (int moves = 1; moves <= monte_step; moves++)
    {   
        dnaChain.countmove();
        const double ERROR_TOL=1e-12;
        double moveselection;
        moveselection=drand(1.0);
        if (moveselection < crankratio){
            double rotAng;
            rotAng = (mt()*2-1)*maxRotAng;
            int m; int n;
            // 0<=m<=maxnum+1 
            m = irand(maxnum)+1;
            // 0<=n<=maxnum+1 && abs(n-m)>=crank_min_length
            while (abs((n=irand(maxnum)+1)-m)<crank_min_length);
            if (n<m){
                int mid;
                mid=n;n=m;m=mid;
            }
            double dE = dnaChain.deltaE_TrialCrankshaft(m,n,rotAng);
            dnaChain.crankshaft(m,n,rotAng,true);
         
            bool accept;
            accept=false;
            if (dnaChain.getEndToEndDistance()<endToEndDistanceThreshold
                && dnaChain.getEndToEndAngle()<endToEndAngleThreshold
                && mt()<exp(-dE)){accept=true;}
            if (!accept){
                dnaChain.crankshaft(m,n,-rotAng,true);
            }
            else{
                dnaChain.stats.accepts++;
                dnaChain.stats.crankSuccess++;
                /*if (dnaChain.stats.crankSuccess>380){
                    char tempbuf[30];
                    sprintf(tempbuf,"mv%d_cr_%d.txt",dnaChain.stats.moves(),dnaChain.stats.crankSuccess);
                    dnaChain.snapshot(tempbuf);
                }*/
            }
            if (dbg) anglelog<<'#'<<moves<<':'<<dnaChain.getEndToEndAngle()<<' '<<accept<<endl;
        }
        else{
            double phi=mt()*2-1,theta=mt()*2*PI;//(0.0,2*PI);
            double rv[3]={
                sqrt(1-phi*phi)*cos(theta),
                sqrt(1-phi*phi)*sin(theta),
                phi};
            double chi=(mt()*2-1)*PI*rotpercentage;
            int m=irand(0,maxnum+1);

            double dE = dnaChain.deltaE_TrialHalfChain(m,rv,chi);

            dnaChain.halfChain(m,rv,chi,true);

            bool accept;
            accept=false;
            if (dnaChain.getEndToEndDistance()<endToEndDistanceThreshold 
                && dnaChain.getEndToEndAngle()<endToEndAngleThreshold
                && mt()<exp(-dE)) {accept=true;}
            if (!accept){
                dnaChain.halfChain(m,rv,-chi,true);
            }
            else{
                dnaChain.stats.accepts++;
                dnaChain.stats.halfChainSuccess++;
            }
            if (dbg)anglelog<<'!'<<moves<<':'<<dnaChain.getEndToEndAngle()<<' '<<accept<<endl;
        }

        //Statistics after one movement

        //Record how many conformation are in a smaller volume.
        if (dnaChain.getEndToEndDistance()>endToEndDistanceThreshold_ligate ||
            dnaChain.getEndToEndAngle()>endToEndAngleThreshold_ligate){
        }
        else{
            dnaChain.stats.ligateAccepts++;
            //Binning the configurations according to the number of disruptions in one conformation.
            if (conformationBinningOn){
                this->appendAngleStats();
            }
        }
    }
    if (dbg)anglelog.close();
    return (double)dnaChain.stats.ligateAccepts/(double)dnaChain.stats.moves();
}

