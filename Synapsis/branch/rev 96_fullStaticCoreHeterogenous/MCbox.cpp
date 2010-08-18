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

	if (dnaChain->checkConsistancy()==1){
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

	long tal[2]={0,0},ter=0;
	this->dnaChain->kpoly(tal,ter);
	*fp_log<<"Initial KPoly:"<<tal[0]<<','<<tal[1]<<' '<<ter<<endl;

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

		
		if (drand(1.0)>P_REPT){//==================================================================
		//Crankshaft movement.
			//generate rotation axis, avoiding rigid body.
			long testp;long testflag;
			do{
				m=irand(maxnum+1);
				n=wrap(m+irand(crank_min_length,crank_max_length+1),totsegnum);
				//Check if containting any rigid body segments.
				testp=m;testflag=0;
				//Check from vertex m to n, ensure they are not protected inside rigidbody.
				if (protect_list[testp]==1){
					testflag=1;
				}
				else{
					do{
						testp=wrap(testp+1,totsegnum);
						if (protect_list[testp]==1){
							testflag=1;
							break;
						}
					}while (testp!=n);
				}
				//if (testflag==0) break;

				//From seg n to m-1
				//This section disabled since crank_max_length < totsegnum/2.
				//If it can't pass the n~m-1 test, it won't pass the following one either
				//since m~n-1 is longer than totsegnum/2.
			   /* testp=n;testflag=0;
				while (testp!=wrap(m,totsegnum)){
					if (protect_list[testp]==1){
						testflag=1;
						break;
					}
					testp=wrap(testp+1,totsegnum);
				}
				if (testflag==0) break;*/
			}while(testflag==1);

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
			}
		}//End Crankshaft movement.
		else{//==================================================================================
		//Reptation movement.
			//generate reptation segment.
			//The segment between vertices m and n will be changed.
			//that is vectors m~n-1
			long testp;long testflag;
			do{
				m=irand(maxnum+1);
				n=wrap(m+irand(reptation_minlen,reptation_maxlen+1),totsegnum);
				//Check if containting any rigid body segments.
				testp=m;testflag=0;
				while (testp!=wrap(n+1,totsegnum)){
					if (protect_list[testp]==1){
						testflag=1;
						break;
					}
					testp=wrap(testp+1,totsegnum);
				}
			}while(testflag==1);
			
			long rept_move;
			rept_move=0;
			while(rept_move==0)
				//[-rept_move_range,rept_move_range] excluding 0
				rept_move=irand(-rept_move_range,rept_move_range+1); 
			
			//Stats:
			this->dnaChain->auto_updt_stats();

     		//old rigid body energy.
			cacheRE=RG.E;

			//old writhe energy;
			cacheE_t=dnaChain->E_t;
			cacheWrithe=dnaChain->writhe;

			//bend energy change and movement.
			dE=dnaChain->dE_reptation(m,n,rept_move);

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
					dnaChain->dE_reptation(m,n,-rept_move);
					RG.update_allrigid_and_E();
					dnaChain->writhe=cacheWrithe;
					dnaChain->E_t=cacheE_t;
			}

		}//End reptation movement.

		if (moves%SNAPSHOT_INTERVAL==0){
			sprintf(buf,"%s%09d.txt",filePrefix,moves);
			cout <<buf<<endl;
			dnaChain->snapshot(buf);
		}

		if (moves%STAT_INTERVAL==0){
			(*fp_log)<<"accepted:"<<dnaChain->stats.crk_accepts()
				<<" rpt_accepted:"<<dnaChain->stats.rpt_accepts()
				<<" in moves "<<dnaChain->stats.auto_moves()
				<<'['
				<<float(dnaChain->stats.crk_accepts())/dnaChain->stats.auto_moves()<<","
				<<float(dnaChain->stats.rpt_accepts())/dnaChain->stats.auto_moves()
				<<']';
			dnaChain->stats.crk_accepts.lap();
			dnaChain->stats.rpt_accepts.lap();
			dnaChain->stats.auto_moves.lap();
			(*fp_log)<<endl;
//----------Log acceptance and rigid body statistics.------------
			long ial[2],ierr=0;
			this->dnaChain->kpoly(ial,ierr);
			(*fp_log)<<"["<<moves<<"]";
			(*fp_log)<<" move_trial["<<m<<","<<n<<"]";
//			(*fp_log)<<" Branch="<<dnaChain->getBranchNumber();
			(*fp_log)<<" Winding[Wr,E_t]"<<dnaChain->writhe<<","<<dnaChain->E_t;
			(*fp_log)<<" Flags(E,rigidIEV,IEV,topo)"<<"["
				<<E_condition<<"(dE="<<dE<<"),"
				<<rigid_IEV_condition<<","
				<<IEV_condition<<'('<<info[0]<<','<<info[1]<<')'
				<<","<<topo_condition<<".topl:"<<dnaChain->topl
				<<"KPoly("<<ial[0]<<','<<ial[1]<<')'<<"]";

//			Log AlexPoly(s,t)~Linking Number of recombination products.
			(*fp_log)<<" Lk_recomb="<<dnaChain->productLk(RG.R[0].protect[1],RG.R[1].protect[1]);

//			Log Rigidbody status
/*		    (*fp_log)<<" r "<<RG.r<<" Ax "<<180-RG.AxisBeta/PI*180
				<<" Ra "<<180-RG.RadiusBeta/PI*180<<" E "<<RG.E<<endl; */

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