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
	stringstream(config[string("reptation_maxlen")])>>reptation_maxlen;
	stringstream(config[string("reptation_minlen")])>>reptation_minlen;
	stringstream(config[string("rept_move_range")])>>rept_move_range;
	stringstream(config[string("rept_min_seglength")])>>rept_min_seglength;

	stringstream(config[string("RBAUS_LOAD_LAST")])>>RBAUS_LOAD_LAST;
	stringstream(config[string("RBAUS_COLLECT_ENABLED")])>>RBAUS_COLLECT_ENABLED;


	stringstream(config[string("VolEx_cutoff_rigidbody")])>>VolEx_cutoff_rigidbody;

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

	//-------Initialize dnaChain writhe info--------
	dnaChain->E_t_updateWrithe_E_t();

    //Initialize statistical variables.
    for (long i = 0; i < 180; i++)
        anglenums[i] = 0;
    //Log file handle.
    //dnaChain->snapshot(strcat_noOW(buf, strBufSize, filePrefix, "_ini.txt"));
    logParameters();

	if (dnaChain->checkConsistency()==1){
		cout<<"Chain inconsistency occurs. In most cases this is caused by inconsistency"
			"between alleged number of segment in _config file and the actual number of"
			" segments."<<endl;
		exit(EXIT_FAILURE);
	}
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

void MCbox_circular::performMetropolisCircularCrankRept(long monte_step)
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
	this->dnaChain->kpoly(tal,ter);
	*fp_log<<"Initial KPoly:"<<tal[0]<<','<<tal[1]<<' '<<ter<<endl;

	if (RBAUS_LOAD_LAST) {
		U.load("ArtificialPotential.txt");
		RG.update_allrigid_and_E();
	}

	double E=dnaChain->calG_bSum()+RG.E+dnaChain->E_t;

	for (long moves = 1; moves <= monte_step; moves++)
    {
		//MAKE MOVES

		//NOTES ON M AND N. PAY ATTENTION.
		//M<N is not required here since dE_TrialCrankshaft and crankshaft
		//support M>N and automatically wraps the iterator to the 
		//beginning of the chain when it hits the tail.
		//Therefore, all energy evaluation program should be careful with chain segment 
		//iteration due to wrapping problem.

		double dE, cacheRE, cacheE_t, cacheWrithe;
		long m,n;
		long E_condition=0,IEV_condition=0,topo_condition=0,rigid_IEV_condition=0;
		double info[3]={0,0,0},info_old[3]={0,0,0};
		if (RBAUS_COLLECT_ENABLED) {
			double flag;
			flag=U.collect(RG.Q);

			//Artificial potential U changes, and RE.E should be updated.
			if (flag>0) RG.update_allrigid_and_E(); 
		}
		
		if (drand(1.0)>=P_REPT){//==================================================================
		//Crankshaft movement.
			//generate rotation axis, avoiding rigid body.
			long testp;long testflag;
			do{
				m=irand(maxnum+1);
				n=wrap(m+irand(crank_min_length,crank_max_length+1),totsegnum);
			}while(protect_list[m]==1 || protect_list[n]==1);

			//(*this->fp_log)<<m<<' '<<n<<endl;

			//generate rotation angle.
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

			//Stats.
			this->dnaChain->auto_updt_stats();
			this->dnaChain->stats.crk_counts++;

			//old rigid body energy.
			cacheRE=RG.E;
			//old writhe energy;
			cacheE_t=dnaChain->E_t;
			cacheWrithe=dnaChain->writhe;

			//bend energy change and movement.
			dE=dnaChain->dE_TrialCrankshaft(m, n, rotAng);
			dnaChain->crankshaft(m,n,rotAng);

			//update energies.
			RG.update_allrigid_and_E();
			this->dnaChain->E_t_updateWrithe_E_t();
			//total energy change.
			dE= dE + (RG.E - cacheRE) + (dnaChain->E_t - cacheE_t);
			
			//Flags: 0 - uninitialized -1 - not satisfied +1 - satisfied
			E_condition=0;
			IEV_condition=0;
    		topo_condition=0;
			rigid_IEV_condition=0;

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
				if (RG.IEV_spheres(m,n)==1){
					rigid_IEV_condition=1;
				}
				else{
					rigid_IEV_condition=-1;
				}
			} 
			
			if (rigid_IEV_condition==1){//rigid_IEV_condition
				/*int IEVflag=this->dnaChain->IEV_with_rigidbody_closeboundary(m,n,info);
				int IEVflag_old=this->dnaChain->IEV_Alex_closeboundary(m,n,info_old);
				if (IEVflag!=IEVflag_old){
					char filebuf[100];
					sprintf(filebuf,"IEVerr_c%d,%d_%010d_(new %d[%3.0f,%3.0f],old %d[%3.0f,%3.0f]).txt",
						m,n,moves,IEVflag,info[0],info[1],IEVflag_old,info_old[0],info_old[1]);
					this->dnaChain->snapshot(filebuf);
				}*/
				if (this->dnaChain->IEV_with_rigidbody_closeboundary(m,n,info)==1){
						IEV_condition=1;
					}
					else{
						IEV_condition=-1;
				}
			}

			if (IEV_condition==1){
				if (this->dnaChain->topl<1.5){ 
					topo_condition=1;
				}
				else{
					topo_condition=-1;
				}
			}

			if (E_condition==1 && rigid_IEV_condition==1 
				&& IEV_condition==1  && topo_condition==1){
					this->dnaChain->stats.crk_accepts++;		
			}
			else{
					dnaChain->crankshaft(m,n,-rotAng);
					RG.update_allrigid_and_E();
					dnaChain->writhe=cacheWrithe;
					dnaChain->E_t=cacheE_t;
					dE=0;
			}
		}//End Crankshaft movement.
		else{//==================================================================================
		//Reptation movement.

/*			if (dnaChain->checkConsistency()==1){
				cout<<moves<<"inconsistency before a reptation movement"<<endl;
			}
*/
			long testp;long testflag1, testflag2;
			long m3,n3,m4,n4;
			do{
				m3=irand(maxnum+1); //3-segment
				n3=wrap(m3+2,totsegnum);


				//m3~m3+2 m3+3~~~~~~~~~~~~|-----4---|
				//|--3----|<--max+1-3-4-->|-----4---|
				//        +0             +max+1-3-4
				m4=wrap( (m3+3)+irand((maxnum+1)-3-4+1), totsegnum );//4-segment
				n4=wrap( m4+3, totsegnum);

				//Check if containting any rigid body segments or connecting segments whose length is smaller than 3.0.
				testflag1=0, testflag2=0;

				testp=m3;
				while (testp!=wrap(n4+1,totsegnum)){
					if (dnaChain->C[testp].l < rept_min_seglength){
						testflag1=1;
						break;
					}
					testp=wrap(testp+1,totsegnum);
				}
				if (testflag1==0){
					testflag1=2;
					break;
				}

				testp=m4;
				while (testp!=wrap(n3+1,totsegnum)){
					if (dnaChain->C[testp].l < rept_min_seglength){
						testflag2=1;
						break;
					}
					testp=wrap(testp+1,totsegnum);
				}
				if (testflag2==0) {
					testflag2=2;
					break;
				}

			}while(1);
			
			//Stats:
			this->dnaChain->auto_updt_stats();
			this->dnaChain->stats.rpt_counts++;

     		//old rigid body energy.
			cacheRE=RG.E;

			//old writhe energy;
			cacheE_t=dnaChain->E_t;
			cacheWrithe=dnaChain->writhe;

			//old conformation stored.
		    static segment backC[maxa];
			for (int i=0;i<=maxnum;i++)
				backC[i]=dnaChain->C[i];

			//m2 should be at the positive direction of m1.
			long m1,dm1,m2,dm2;
			if (testflag1==2) {
				m1=m3;dm1=2;
				m2=m4;dm2=3;
			}else if (testflag2==2){
				m1=m4;dm1=3;
				m2=m3;dm2=2;
			}else{
				cout<<"testflag error at performMetropolisCircularCrankRept:Reptation movement."<<endl;
				exit(EXIT_FAILURE);
			}

			int rejection_sign;
			dE=dnaChain->dE_reptation_3_4(m1,dm1,m2,dm2,rejection_sign);	


			if (rejection_sign!=0){
				dE=0;
				for (int i=0;i<=maxnum;i++)	dnaChain->C[i]=backC[i];
				if (rejection_sign==1) dnaChain->stats.rpt_rejection_count++;
				if (rejection_sign==2) dnaChain->stats.rpt_rejection_quickrej++;
				goto goon;
			}

			//update energies.
			RG.update_allrigid_and_E();
			this->dnaChain->E_t_updateWrithe_E_t();

			//total energy change.
			dE= dE + (RG.E - cacheRE) + (dnaChain->E_t - cacheE_t);
			
			//Flags: 0 - uninitialized -1 - not satisfied +1 - satisfied
			E_condition=0;
			IEV_condition=0;
    		topo_condition=0;
			rigid_IEV_condition=0;

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
				if (RG.IEV_spheres(m,n)==1){
					rigid_IEV_condition=1;
				}
				else{
					rigid_IEV_condition=-1;
				}
			} 
			
			if (rigid_IEV_condition==1){//rigid_IEV_condition
				if (this->dnaChain->IEV_with_rigidbody_closeboundary(m,n,info)==1){
						IEV_condition=1;
					}
					else{
						IEV_condition=-1;
				}
			}

			if (IEV_condition==1){
				if (this->dnaChain->topl<1.5){ 
					topo_condition=1;
				}
				else{
					topo_condition=-1;
				}
			}

			if (E_condition==1 && rigid_IEV_condition==1 
				&& IEV_condition==1  && topo_condition==1){
				dnaChain->stats.rpt_accepts++;
			}
			else{
				for (int i=0;i<=maxnum;i++)	dnaChain->C[i]=backC[i];
				RG.update_allrigid_and_E();
				dnaChain->writhe=cacheWrithe;
				dnaChain->E_t=cacheE_t;
				dE=0;
			}

		}//End reptation movement.
		

goon:	E+=dE;

		// length and bangle consistency check.
		/*
		if (dnaChain->checkConsistency()==1){
			cout<<moves<<"dX and X[i+1]-X[i] inconsistency after a reptation movement";
			cout<<endl;
		}

		if (dnaChain->checkBangleConsistency()==1){
			cout<<moves<<"Bangle inconsistency after a reptation movement";
			cout<<endl;
		}*/

		//Engergy tracking check.
		/*
		double recalE;
		recalE=dnaChain->calG_bSum()+RG.E+dnaChain->E_t;
		if (fabs(recalE-E)>1e-5){
			cout<<"Energy inconsistency at move: "<<moves
				<<" recalE="<<recalE<<" E = sum(dE) = "<<E;
			cout<<endl;
		}*/
		

		if (moves%SNAPSHOT_INTERVAL==0){
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
			(*fp_log)<<"crk_accepted:"<<dnaChain->stats.crk_accepts()<<"/"<<dnaChain->stats.crk_counts ()
				<<'['<<float(dnaChain->stats.crk_accepts())/dnaChain->stats.crk_counts()*100<<"%] "
				<<" rpt_accepted:"<<dnaChain->stats.rpt_accepts()<<"/"<<dnaChain->stats.rpt_counts()
				<<'['<<float(dnaChain->stats.rpt_accepts())/dnaChain->stats.rpt_counts()*100<<"%,"
				<<"NoSolRej: "<<float(dnaChain->stats.rpt_rejection_count())/dnaChain->stats.rpt_counts()*100<<"% "
				<<"QuickRej: "<<float(dnaChain->stats.rpt_rejection_quickrej())/dnaChain->stats.rpt_counts()*100<<"% "
				<<"] in "<<dnaChain->stats.auto_moves()<< " moves ";
			dnaChain->stats.crk_accepts.lap();
			dnaChain->stats.crk_counts.lap();
			dnaChain->stats.rpt_accepts.lap();
			dnaChain->stats.rpt_counts.lap();
			dnaChain->stats.rpt_rejection_count.lap();
			dnaChain->stats.rpt_rejection_quickrej.lap();
			dnaChain->stats.auto_moves.lap();
//			(*fp_log)<<endl;


//----------Log acceptance and rigid body statistics.------------
			

			(*fp_log)<<"["<<moves<<"]";
			(*fp_log)<<" Q "<<RG.Q;
			long overpass = dnaChain->overpassing(RG.R[0].protect[0]+1,RG.R[1].protect[0]+1);
			(*fp_log)<<" + "<< overpass;
			long Lk_recomb = dnaChain->productLk(RG.R[0].protect[0]+1,RG.R[1].protect[0]+1);
			(*fp_log)<<" Lk_re "<< Lk_recomb;
//			(*fp_log)<<" move_trial["<<m<<","<<n<<"]";
//			(*fp_log)<<" Branch="<<dnaChain->getBranchNumber();
//			(*fp_log)<<" Winding[Wr,E_t]"<<dnaChain->writhe<<","<<dnaChain->E_t;

//			(*fp_log)	<<" Flags(E,rigidIEV,IEV,topo)"<<"[";
//			(*fp_log)	<<E_condition<<"(dE="<<dE<<"),";
//			(*fp_log)	<<rigid_IEV_condition<<",";
//			(*fp_log)	<<IEV_condition<<'('<<info[0]<<','<<info[1]<<')';
//			(*fp_log)	<<","<<topo_condition<<".topl:"<<dnaChain->topl;

			long ial[2],ierr=0;
			this->dnaChain->kpoly(ial,ierr);
			(*fp_log)	<<" KPoly("<<ial[0]<<','<<ial[1]<<')'<<"]";

//			Log AlexPoly(s,t)~Linking Number of recombination products.
			if ( EXOTIC_LK_SNAPSHOT==1 && Lk_recomb != 1 && overpass==+1 && RG.Q < -5.0 ){
				char LkSnapBuf[100];
				sprintf(LkSnapBuf,"%s_%09d_Lk(%d).txt",this->filePrefix,moves,Lk_recomb);
				this->dnaChain->snapshot(LkSnapBuf);
			}

//			Log Rigidbody status
			(*fp_log)
/*				<<" "<<RG.r
				<<" "<<RG.AxisBeta/PI*180
				<<" "<<RG.RadiusBeta/PI*180
				<<" "<<RG.r_siteI
				<<" "<<RG.siteI_direction
				<<" "<<RG.E
				<<" "<<RG.Q;*/
				<<" r "<<RG.r
				<<" Ax "<<180-RG.AxisBeta/PI*180
				<<" Ra "<<180-RG.RadiusBeta/PI*180
				<<" r_siteI "<<RG.r_siteI
//				<<" SiteIDir "<<RG.siteI_direction
				<<" E "<<RG.E;

//			Log Gyration Radius
/*			double gyration_ratio=this->calcGyration();
			dnaChain->stats.gyration_ratio.push(gyration_ratio);
			(*fp_log)<<endl<<"$"<<moves<<"] ";
			(*fp_log)<<"Rg "<<gyration_ratio
				<<"<"<<dnaChain->stats.gyration_ratio.getMean()<<"+/-"
				<<dnaChain->stats.gyration_ratio.getStdev()<<"/meanstd"
				<<dnaChain->stats.gyration_ratio.getStdevOfMean()<<">"<<endl;
*/

//			Log Chain angle statistics
/*			for (long i=0;i<=maxnum;i++)
				dnaChain->stats.anglelist[i].push(dnaChain->C[i].bangle);
*/			(*fp_log)<<endl;
		}
	}

	for (long i=0;i<=maxnum;i++){
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