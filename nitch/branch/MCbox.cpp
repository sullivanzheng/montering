#include "MCbox.h" 
MCbox_circular::MCbox_circular(char const *configFile){	
	using namespace std;
	//Read Config File;
	config=readconfig(configFile);
	for (long i = 0; i<maxa; i++){
		protect_list[i]=0;
	}
	//#############################GLOBAL##############################
	stringstream(config[string("totsegnum")])>>totsegnum;
	maxnum= totsegnum -1;

	stringstream(config[string("crank_min_length")])>>crank_min_length;
	stringstream(config[string("crank_max_length")])>>crank_max_length;
	stringstream(config[string("VEcutoff")])>>VEcutoff;

	stringstream(config[string("bpperunit")])>>bpperunit;

	stringstream(config[string("maxRotAng")])>>maxRotAng;
	maxRotAng=maxRotAng/180*PI;
	stringstream(config[string("P_SMALLROTATION")])>>P_SMALLROTATION;
	
	stringstream(config[string("P_REPT")])>>P_REPT;
	stringstream(config[string("P_REPT_SIMP")])>>P_REPT_SIMP;
	stringstream(config[string("reptation_maxlen")])>>reptation_maxlen;
	stringstream(config[string("reptation_minlen")])>>reptation_minlen;
	stringstream(config[string("rept_move_range")])>>rept_move_range;
	stringstream(config[string("rept_min_seglength")])>>rept_min_seglength;

	stringstream(config[string("P_TREADMILL")])>>P_TREADMILL;

	stringstream(config[string("P_HALFCHAIN")])>>P_HALFCHAIN;

	stringstream(config[string("RBAUS_LOAD_LAST")])>>RBAUS_LOAD_LAST;
	stringstream(config[string("RBAUS_COLLECT_ENABLED")])>>RBAUS_COLLECT_ENABLED;


	stringstream(config[string("VolEx_cutoff_rigidbody")])>>VolEx_cutoff_rigidbody;
	stringstream(config[string("VolEx_R_rigid")])>>VolEx_R_rigid;

	stringstream(config[string("initial_guess_siteII_umbrella_energy")])>>initial_guess_siteII_umbrella_energy;
	stringstream(config[string("initial_guess_siteI_umbrella_energy")])>>initial_guess_siteI_umbrella_energy;

	stringstream(config[string("BREAKANGLE")])>>BREAKANGLE;
	stringstream(config[string("SPECIAL_ANGLE")])>>SPECIAL_ANGLE;
	stringstream(config[string("G2")])>>G2;


	//#############################MCBox Variables#####################
	strcpy(filePrefix,config[string("filePrefix")].c_str());

	this->fp_log = 
		new ofstream(strcat_noOW(buf, strBufSize,const_cast<char*>(filePrefix), "_log.txt"));
    this->dnaChain=
		new CircularChain(strcat_noOW(buf, strBufSize, const_cast<char*>(filePrefix), ".vmc"),totsegnum);
    
	stringstream(config[string("seeding")])>>this->seeding;

	//#######################CircularChain Variables##################
	stringstream(config[string("VolEx_R")])>>this->dnaChain->VolEx_R;
	stringstream(config[string("dLk")])>>this->dnaChain->dLk;

	/*if (dnaChain->checkConsistency(5e-5)==1){
		cout<<"Chain inconsistency occurs. In most cases this is caused by inconsistency "
			"between alleged number of segment in _config file and the actual number of "
			"segments."<<endl;
		exit(EXIT_FAILURE);
	}*/
	dnaChain->normalize_X_bangle(false);

	//-------Initialize dnaChain writhe info--------
	dnaChain->E_t_updateWrithe_E_t();

    //Initialize statistical variables.
    for (long i = 0; i < 180; i++)
        anglenums[i] = 0;
    //Log file handle.
    //dnaChain->snapshot(strcat_noOW(buf, strBufSize, filePrefix, "_ini.txt"));
    logParameters();

}

void MCbox_circular::logAngleDist(char *suffix)
{
    strcat_noOW(buf, strBufSize, filePrefix, "_angleDist.txt");
    ofstream fp (strcat_noOW(buf, strBufSize, buf, suffix));
    fp << "Angle distribution" << endl;
    for (long i = 0; i < 180; i++)
        fp << i << '\t' << anglenums[i] << endl;
    fp << endl;
    fp.close();
}

void MCbox_circular::clearAngleStats(void)
{
    for (long i = 0; i < 180; i++)
        anglenums[i] = 0;
}

void MCbox_circular::pushAngleStats(void)
{
    for (long i = 0; i < maxnum; i++)
    {
        anglenums[long(floor(dnaChain->C[i].bangle * 180 / PI))]++;
    }
}

void MCbox_circular::logAccepts(void)
{
    sprintf(buf, "%12d %12d ", 
        dnaChain->stats.auto_moves.getTotCounts(), 
        dnaChain->stats.crk_accepts.getNumber());
    (*fp_log) << buf;
}

double MCbox_circular::ligationP(){
	double dis_f=100,ang_f=PI;
	double dis=0.1,ang=5.0/180.0*PI;
	int n=10;

	double dd=pow(dis_f/dis,1.0/n),da=pow(ang_f/ang,1.0/n);
	
	double totalP=1.0,condP=1.0;
	double d,a;
	int i;
	//cout << this->performMetropolisCircularCrankRept(1e6,100,PI,50,PI);
	
	for (d=dis*dd,a=ang*da,i=0;i<n;d*=dd,a*=da,i++){
		cout << "Simlulation: DIS ("<<d<<" -> "<<d/dd
		<<") ANG (" <<a<<" -> "<<a/da<<")" << i <<" before totalP=" <<totalP <<endl;

		(*fp_log) << "Simlulation: DIS ("<<d<<" -> "<<d/dd
		<<") ANG (" <<a<<" -> "<<a/da<<")" << i <<" before totalP=" <<totalP <<endl;
		condP=this->performMetropolisCircularCrankRept(1e8,d,a,d/dd,a/da);

		totalP*=condP;
		
		cout << ">>condP=" <<condP<<endl;

		(*fp_log)<< ">>condP=" <<condP<<endl;
	}

	cout << "=============== Ligation P:" <<totalP<<" ==============="<<endl;

	(*fp_log)<< "=============== Ligation P:" <<totalP<<" ==============="<<endl;

	double coef=1/6.022141 * (3.0/ (4.0 * PI * (pow(0.34*dis,3))) ) *  2.0/(1-cos(ang)) * 1e4; //Verify
	(*fp_log)<< ">>>>>>>>>>>>>>> J factor: "<<coef*totalP << endl; //TODO
	return totalP;
}

double MCbox_circular::performMetropolisCircularCrankRept(long monte_step,
									  double endToEndDistanceThreshold,
									  double endToEndAngleThreshold,
									  double endToEndDistanceThreshold_ligate,
									  double endToEndAngleThreshold_ligate)
{
    MTRand53 mt(seeding);
	allrigid RG("_rigid.cfg", this->dnaChain);

	long SNAPSHOT_INTERVAL;
	std::stringstream(config["SNAPSHOT_INTERVAL"])>>SNAPSHOT_INTERVAL;

	long STAT_INTERVAL;
	std::stringstream(config["STAT_INTERVAL"])>>STAT_INTERVAL;

    long RBAUS_PICKLE_INTERVAL;
	std::stringstream(config["RBAUS_PICKLE_INTERVAL"])>>RBAUS_PICKLE_INTERVAL;

	long EXOTIC_LK_SNAPSHOT;
	std::stringstream(config["EXOTIC_LK_SNAPSHOT"])>>EXOTIC_LK_SNAPSHOT;

	long tal[2]={0,0},ter=0;

	dnaChain->stats.resetStat();

	if (RBAUS_LOAD_LAST) {
		U.load("ArtificialPotential.txt");
		RG.update_allrigid_and_E();
	}



	double E=dnaChain->calG_bSum();

	int debugsignal=0;

	for (long moves = 1; moves <= monte_step; moves++)
    {
		//MAKE MOVES

		//NOTES ON M AND N. PAY ATTENTION.
		//M<N is not required here since dE_TrialCrankshaft and crankshaft
		//support M>N and automatically wraps the iterator to the 
		//beginning of the chain when it hits the tail.
		//Therefore, all energy evaluation program should be careful with chain segment 
		//iteration due to wrapping problem.

		double dE, cacheRE, cacheE_t, cacheWrithe, cacheTopl, WrChangeInTrialMove,E_tChangeInTrialMove;
		long m,n; double rotAng;
		long E_condition=0,IEV_condition=0,Boundary_condition=0;
		double info[3]={0,0,0},info_old[3]={0,0,0};
		if (RBAUS_COLLECT_ENABLED) {
			double flag;
			flag=U.collect(RG.Q);

			//Artificial potential U changes, and RE.E should be updated.
			if (flag>0) RG.update_allrigid_and_E(); 
		}
		double selection=drand(1.0);
		int movement=0;
		movement=(selection<P_REPT? 1 : (								//1: Reptation
						selection<P_REPT+P_REPT_SIMP? 2:(				//2: Simplified reptation
							selection<P_REPT+P_REPT_SIMP+P_TREADMILL? 3:(//3: Treadmill
								selection<P_REPT+P_REPT_SIMP+P_TREADMILL+P_HALFCHAIN?4: //4:Halfchain
							0)											//0: otherwise: crankshaft
				  )));
		static char const movement_symbol[][4]={"CRK","REP","SRP","TRD","HLF"};
		

		if (movement==0){//==================================================================
		//Crankshaft movement.
			//generate rotation axis, avoiding rigid body.
			long testp;long testflag;
			do{
				m=irand(1,maxnum+1);
				n=wrap(m+irand(crank_min_length,crank_max_length+1),totsegnum);
			}while(n==0);

			//(*this->fp_log)<<m<<' '<<n<<endl;

			//generate rotation angle.
			double selection;
			selection = drand(1.0);
			if (selection > P_SMALLROTATION)
			{
				rotAng = drand(-maxRotAng, maxRotAng)*4;
				//if (rotAng <= 0)
				//	rotAng -= maxRotAng;
				//else
				//	rotAng += maxRotAng;
			}
			else
				rotAng = drand(-maxRotAng, maxRotAng);

			//Stats.
			this->dnaChain->auto_updt_stats();
			this->dnaChain->stats.crk_counts++;

			//old Conformation
			static segment backC[maxa];
			for (int i=0;i<=maxnum;i++)	backC[i]=dnaChain->C[i];

			//bend energy change and movement.
			dE=dnaChain->dE_TrialCrankshaft(m, n, rotAng);
			dnaChain->crankshaft(m,n,rotAng);
			
			//Flags: 0 - uninitialized -1 - not satisfied +1 - satisfied
			E_condition=0;
			IEV_condition=0;
			Boundary_condition=0;

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
					 E_condition=-1;
				}
			}
			
			if (E_condition==1){//rigid_IEV_condition
				if (this->dnaChain->IEV_with_rigidbody_closeboundary(m,n,info)==1){
						IEV_condition=1;
					}
					else{
						IEV_condition=-1;
				}
			}
			if (IEV_condition==1){
				Boundary_condition=1;
			}

			if (E_condition==1 && IEV_condition==1 && Boundary_condition==1){
					this->dnaChain->stats.crk_accepts++;		
			}
			else{					
					for (int i=0;i<=maxnum;i++) dnaChain->C[i]=backC[i];
			}
		}//End Crankshaft movement.
		else if(movement==4){
			this->dnaChain->auto_updt_stats();
			this->dnaChain->stats.hlf_counts++;

            double phi=mt()*2-1,theta=mt()*2*PI;//(0.0,2*PI);
            double rv[3]={
                sqrt(1-phi*phi)*cos(theta),
                sqrt(1-phi*phi)*sin(theta),
                phi};
            double chi=(mt()*2-1) * maxRotAng;
            int m=irand(0,maxnum+1);

			double dE = dnaChain->deltaE_TrialHalfChain(m,rv,chi);
			dnaChain->halfChain(m,rv,chi);


			E_condition=0;
			IEV_condition=0;
			Boundary_condition=0;

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
					 E_condition=-1;
				}
			}
			

			if (E_condition==1){
				if (dnaChain->getEndToEndDistance()<endToEndDistanceThreshold 
					&& dnaChain->getEndToEndAngle()<endToEndAngleThreshold
					) {Boundary_condition=1;}
				else
					{Boundary_condition=-1;}
			}

			if (Boundary_condition==1){//rigid_IEV_condition
				if (this->dnaChain->IEV_with_rigidbody_closeboundary(m,maxnum+1,info)==1){
						IEV_condition= 1;
					}
					else{
						IEV_condition=-1;
				}
			}
			
            if (E_condition!=1 || IEV_condition!=1 || Boundary_condition!=1 ){
                dnaChain->halfChain(m,rv,-chi);
            }
            else{
                dnaChain->stats.hlf_accepts++;
            }
		}


goon:	if (E_condition==1 && IEV_condition==1 && Boundary_condition==1)
				E+=dE;
		if (dnaChain->getEndToEndDistance() < endToEndDistanceThreshold_ligate 
			&& dnaChain->getEndToEndAngle() < endToEndAngleThreshold_ligate)
			dnaChain->stats.ligation_count++;


		debugsignal=0;
		static unsigned long const CONSISTENCY_CHECK=7321;
		if (moves%CONSISTENCY_CHECK==0 && debugsignal==1){

			//Engergy tracking check.
			/*
			double recalE;
			recalE=dnaChain->calG_bSum()+RG.E+dnaChain->E_t;
			if (fabs(recalE-E)>1e-5){
				cout<<" Energy inconsistency at move: "<<moves
					<<" recalE="<<recalE<<" E = sum(dE) = "<<E;
				cout<<endl;
			}
			*/

			// length and bangle consistency check.
			if (dnaChain->checkConsistency()==1){
				cout<<"--"<<moves<<" dX and X[i+1]-X[i] inconsistency, probably caused by reptation";
				cout<<endl<<endl;
			}

			if (dnaChain->checkBangleConsistency()==1){
				cout<<"--"<<moves<<" Bangle inconsistency, probably caused by reptation";
				cout<<endl<<endl;
			}
			//IEV full-check.
			if(this->dnaChain->IEV_with_rigidbody_closeboundary_fullChain(info)!=1){
					cout<<"--"<<moves<<" IEV error,"<<info[0]<<','<<info[1];
					cout<<endl;
			}

			//Topology consistency
			long topo[2]={0},errorcode=0,topo2[2]={0},
	        		errorcode2=0,errorcodeKNDWR=0;

			this->dnaChain->kpoly(topo,errorcode);
			this->dnaChain->kpoly2(topo2,errorcode2);

			double dtopo;
			this->dnaChain->_kndwr_topl_update_stable(dtopo, errorcodeKNDWR);

			if (topo[0]!=int(dtopo) || topo2[0]!=topo[0]){
				cout<<moves<<": topology inconsistency. Deadlock triggerred for debug";
				cout<<endl;

				char filebuf[200];
				sprintf(buf,"topo%08d_%2d,%2d_ERR%2d _%2d,%2d_ERR%2d.txt",
					moves,topo[0],topo[1],errorcode,
					      topo2[0],topo[1],errorcode2);
				this->dnaChain->snapshot(buf);

				//Repeat error:
				for (;;){//DEBUG
					this->dnaChain->_kndwr_topl_update_stable(dtopo, errorcodeKNDWR);
					this->dnaChain->kpoly(topo,errorcode);
					this->dnaChain->kpoly2(topo2,errorcode2);
				}
			}
		}

		if (moves%SNAPSHOT_INTERVAL==0 ){
			sprintf(buf,"%s%09d.txt",filePrefix,moves);
			cout <<buf<<endl;
			dnaChain->snapshot(buf);
		}

		if (RBAUS_COLLECT_ENABLED && moves%RBAUS_PICKLE_INTERVAL==0){
			sprintf(buf,"[%09d]Potential profile updated.",moves);
			cout <<buf<<endl;
			U.pickle("ArtificialPotential.txt");
		}

		if (moves%STAT_INTERVAL==0){
			(*fp_log).precision(5);
			char buf[400];
			sprintf(buf,"crk_acpt:%3d/%3d[%5.2f%%] NOW:m,n,beta(%3d,%3d),%5.1f "
						"hlf_actp:%3d/%3d[%5.2f%%] in %4d moves || ",
					dnaChain->stats.crk_accepts(),dnaChain->stats.crk_counts(),
					float(dnaChain->stats.crk_accepts())/dnaChain->stats.crk_counts()*100,
					m,n,rotAng*180.0/PI,

					dnaChain->stats.hlf_accepts(),dnaChain->stats.hlf_counts(),
					float(dnaChain->stats.hlf_accepts())/dnaChain->stats.hlf_counts()*100,

					dnaChain->stats.auto_moves()
				);
			(*fp_log)<<buf;

			sprintf(buf,"DIS %8.3f ANG %8.3f ||",
				dnaChain->getEndToEndDistance(),dnaChain->getEndToEndAngle()/PI*180.0);
			(*fp_log)<<buf;

			sprintf(buf,"LIGATION: interval %d/%d[%5.2f%%] all %d/%d[%5.2f%%]",
				dnaChain->stats.ligation_count(),dnaChain->stats.auto_moves(),
				float(dnaChain->stats.ligation_count())/dnaChain->stats.auto_moves(),

				dnaChain->stats.ligation_count.getTotCounts(),
				dnaChain->stats.auto_moves.getTotCounts(),
				float(dnaChain->stats.ligation_count.getTotCounts())/
				dnaChain->stats.auto_moves.getTotCounts()
			);
			(*fp_log)<<buf;	

			dnaChain->stats.crk_accepts.lap();
			dnaChain->stats.crk_counts.lap();
			dnaChain->stats.rpt_accepts.lap();
			dnaChain->stats.rpt_counts.lap();
			dnaChain->stats.rpt_rejection_count.lap();
			dnaChain->stats.rpt_rejection_quickrej.lap();
			dnaChain->stats.rptsimp_accepts.lap();
			dnaChain->stats.rptsimp_counts.lap();
			dnaChain->stats.tdm_accepts.lap();
			dnaChain->stats.tdm_counts.lap();
			dnaChain->stats.hlf_accepts.lap();
			dnaChain->stats.hlf_counts.lap();
			dnaChain->stats.auto_moves.lap();
			dnaChain->stats.ligation_count.lap();

			(*fp_log)<<endl;
		}
	}
	double condP=float(dnaChain->stats.ligation_count.getTotCounts())/
			dnaChain->stats.auto_moves.getTotCounts();
	(*fp_log)<<"====Conditional probability for this is:"<<condP<<endl;
	return	condP;
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
			<<">totsegnum "<<totsegnum<<"  bpperunit"<<bpperunit<<endl
			<<">maxRotAng	= "<<	maxRotAng	<<endl
			<<" P_SMALLROTATION	="<< P_SMALLROTATION<<endl
			<<">FilePrefix = "<<filePrefix<<endl
			<<endl;
}

double MCbox_circular::calcGyration(void){
	long i;
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