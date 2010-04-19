#include "MCbox.h" 
MCbox_circular::MCbox_circular(
    char const *configFile){	
	using namespace std;
	//Read Config File;
	config=readconfig(configFile);

	stringstream(config[string("length")])>>totsegnum;
	maxnum= totsegnum -1;

	stringstream(config[string("crank_min_length")])>>crank_min_length;
	stringstream(config[string("crank_max_length")])>>crank_max_length;
	stringstream(config[string("VEcutoff")])>>VEcutoff;

	stringstream(config[string("g")])>>g;

	stringstream(config[string("bpperseg")])>>bpperseg;

	stringstream(config[string("maxRotAng")])>>maxRotAng;
	maxRotAng=maxRotAng/180*PI;

	stringstream(config[string("P_SMALLROTATION")])>>P_SMALLROTATION;
	
	strcpy(filePrefix,config[string("filePrefix")].c_str());

	this->fp_log = 
		new ofstream(strcat_noOW(buf, strBufSize,const_cast<char*>(filePrefix), "_log.txt"));
    this->dnaChain=
		new CircularChain(strcat_noOW(buf, strBufSize, const_cast<char*>(filePrefix), ".vmc"),totsegnum);
    
	stringstream(config[string("VolEx_R")])>>this->dnaChain->VolEx_R;

	stringstream(config[string("seeding")])>>this->seeding;

    //Initialize statistical variables.
    for (int i = 0; i < 180; i++)
        anglenums[i] = 0;
    //Log file handle.
    //dnaChain->snapshot(strcat_noOW(buf, strBufSize, filePrefix, "_ini.txt"));
    logParameters();
	dnaChain->snapshot("000.txt");
}

void MCbox_circular::logAngleDist(char *suffix)
{
    strcat_noOW(buf, strBufSize, filePrefix, "_angleDist.txt");
    ofstream fp (strcat_noOW(buf, strBufSize, buf, suffix));
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
}

void MCbox_circular::pushAngleStats(void)
{
    for (int i = 0; i < maxnum; i++)
    {
        anglenums[int(floor(dnaChain->C[i].bangle * 180 / PI))]++;
    }
}

void MCbox_circular::logAccepts(void)
{
    sprintf(buf, "%12d %12d ", 
        dnaChain->stats.auto_moves.getTotCounts(), 
        dnaChain->stats.accepts.getNumber());
    (*fp_log) << buf;
}

void MCbox_circular::performMetropolisCircularCrankOnly(long monte_step)
{
    MTRand53 mt(seeding);
	
	allrigid RG("rigid.txt", this->dnaChain);

	long SNAPSHOT_INTERVAL;
	std::stringstream(config["SNAPSHOT_INTERVAL"])>>SNAPSHOT_INTERVAL;

	long STAT_INTERVAL;
	std::stringstream(config["STAT_INTERVAL"])>>STAT_INTERVAL;

	for (int moves = 1; moves <= monte_step; moves++)
    {
		//MAKE MOVES

		//generate rotation axis
        int m,n;
		int flag=0;
		while(flag==0){
			m = irand(maxnum + 1);
			n = (m + (irand(crank_max_length - crank_min_length) + crank_min_length))%totsegnum;
			flag=1;int i;
			for(i=0;i<RG.protect.size();i++){
				if (m==RG.protect[i] || n==RG.protect[i]) {
					flag=0;
					break;
				}
			}
		}
		//NOTE ABOUT M AND N. PAY ATTENTION.
		//M<N is not necessary here since deltaE_TrialCrankshaft_countMove and crankshaft
		//support M>N and automatically wraps the iterator when it hits the tail to the 
		//beginning of the chain.
		//Therefore, all energy evaluation program should be careful with chain segment 
		//iteration due to wrapping problem.

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

        double dE, cacheRE;

		//bend energy change.
        dE=dnaChain->deltaE_TrialCrankshaft_countMove(m, n, rotAng);

		//old rigid body energy.
		cacheRE=RG.E;
		dnaChain->crankshaft(m,n,rotAng);
		RG.update_allrigid_and_E();
		//total energy change.
		dE+=(RG.E - cacheRE);
		
		int E_condition=0;
		int IEV_condition=0;
		int topo_condition=0;

		if (dE < 0){
			E_condition=1;
		}
        else
        {
            double exp_E;
            double r;
            exp_E = exp(-dE);
            r = mt();
			if (r < exp_E){
                 E_condition=1;
			}
			else{
				 E_condition=0;
			}
        }

		if (this->dnaChain->IEV(m,n)==1){
			IEV_condition=1;
		}
		else{
			IEV_condition=0;
		}

		int ial[2],ierr;
		this->dnaChain->kpoly(ial, ierr);
		if (ierr!=0){
			cout<<"The crossing on the chain is too many!"<<endl;
			exit(EXIT_FAILURE);
		}

		if (ial[0]==1 && ial[1]==1){
			topo_condition=1;
		}else{
			topo_condition=0;
		}

		if (E_condition==1 && IEV_condition==1 && topo_condition==1){
				this->dnaChain->stats.accepts++;		
		}
		else{
				dnaChain->crankshaft(m,n,-rotAng);
				RG.update_allrigid_and_E();
		}

		if (moves%SNAPSHOT_INTERVAL==0){
			sprintf(buf,"%s%09d.txt",filePrefix,moves);
			dnaChain->snapshot(buf);
		}

		if (moves%STAT_INTERVAL==0){
/*			(*fp_log)<<"accepted:"<<dnaChain->stats.accepts()
				<<" in moves "<<dnaChain->stats.auto_moves()
				<<'['<<float(dnaChain->stats.accepts())/dnaChain->stats.auto_moves()
				<<']'<<endl;
			dnaChain->stats.accepts.lap();
			dnaChain->stats.auto_moves.lap();
*/			(*fp_log)<<"["<<moves<<"] ";
			(*fp_log)<<"move_trial["<<m<<","<<n<<"] ";
			(*fp_log)<<"Flags(E,IEV,topo)"<<"["<<E_condition<<"(dE="<<dE<<"),"
				<<IEV_condition<<","<<topo_condition<<"]";
		    (*fp_log)<<"tp_trial("<<ial[0]<<','<<ial[1]<<")";//<<endl;
			(*fp_log)<<" r "<<RG.r<<" Ax "<<180-RG.AxisBeta/PI*180
				<<" Ra "<<180-RG.RadiusBeta/PI*180<<" E "<<RG.E<<endl;
//			double gyration_ratio=this->calcGyration();
//			dnaChain->stats.gyration_ratio.push(gyration_ratio);
//			for (int i=0;i<=maxnum;i++)
//				dnaChain->stats.anglelist[i].push(dnaChain->C[i].bangle);
		}
		
	}
	for (int i=0;i<=maxnum;i++){
		(*fp_log)<<i<<" "<<dnaChain->stats.anglelist[i].getMean()<<" "
			<<dnaChain->stats.anglelist[i].getStdev()<<" "<<endl;
	}
	(*fp_log)<<endl;
}

void MCbox_circular::logParameters(void){
		(*fp_log) 
			<<"========================CONSTANTS========================"<<endl
			<<" PI"<<PI
			<<" maxa	= "<<	maxa	<<endl
			<<" crank_min_length	= "<<	crank_min_length	<<endl
			<<"===================GLOBAL VARIABLES====================="<<endl
			<<" maxnum	= "<<	maxnum	<<endl
			<<" crank_max_length (auto_generated)= "<<	crank_max_length	<<endl
			<<" g "<<g<<endl
			<<">totsegnum * bpperseg = total bp \t"<<totsegnum<<'*'<<bpperseg<<'='<<totsegnum*bpperseg<<endl
			<<">maxRotAng	= "<<	maxRotAng	<<endl
			<<" P_SMALLROTATION	="<< P_SMALLROTATION<<endl
			<<">FilePrefix = "<<filePrefix<<endl
			<<endl;
}

double MCbox_circular::calcGyration(void){
	int i;
	double meanx,meany,meanz;
	meanx=meany=meanz=0;
	for(i=0;i<=maxnum;i++){
		meanx+=dnaChain->C[i].x;
		meany+=dnaChain->C[i].y;
		meanz+=dnaChain->C[i].z;
	}
    meanx/=(maxnum+1);
	meany/=(maxnum+1);
	meanz/=(maxnum+1);

	double G=0;
	for (i=0;i<=maxnum;i++){
		G=G+(dnaChain->C[i].x-meanx)*(dnaChain->C[i].x-meanx)
			+(dnaChain->C[i].y-meany)*(dnaChain->C[i].y-meany)
			+(dnaChain->C[i].z-meanz)*(dnaChain->C[i].z-meanz);
	}
	G=G/(maxnum+1);
	return G;
}


