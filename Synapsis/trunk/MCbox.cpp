#include "MCbox.h" 
MCbox_circular::MCbox_circular(
    char const *configFile){	
	using namespace std;
	//Read Config File;
	config=readconfig(configFile);
	for (int i = 0; i<maxa; i++){
		protect_list[i]=0;
	}
	//#############################GLOBAL##############################
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
	
	stringstream(config[string("P_REPT")])>>P_REPT;
	stringstream(config[string("reptation_maxlen")])>>reptation_maxlen;
	stringstream(config[string("reptation_minlen")])>>reptation_minlen;
	stringstream(config[string("rept_move_range")])>>rept_move_range;

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

    //Initialize statistical variables.
    for (int i = 0; i < 180; i++)
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

void MCbox_circular::performMetropolisCircularCrankRept(long monte_step)
{
    MTRand53 mt(seeding);
	
	allrigid RG("_rigid.cfg", this->dnaChain);

	long SNAPSHOT_INTERVAL;
	std::stringstream(config["SNAPSHOT_INTERVAL"])>>SNAPSHOT_INTERVAL;

	long STAT_INTERVAL;
	std::stringstream(config["STAT_INTERVAL"])>>STAT_INTERVAL;
/*	for (int i=0;i<40;i++){
		dnaChain->C[i].x=i*(i-1)/2;dnaChain->C[i].y=0;dnaChain->C[i].z=0;
		dnaChain->C[i].dx=i;dnaChain->C[i].dy=0;dnaChain->C[i].dz=0;
	}
	dnaChain->snapshot("1.txt");
	dnaChain->dE_reptation(1,32,-1);
	dnaChain->snapshot("2.txt");
	dnaChain->dE_reptation(1,32,1);
	dnaChain->snapshot("3.txt");
	exit(0);*/
	for (int moves = 1; moves <= monte_step; moves++)
    {
		//MAKE MOVES

		//NOTES ON M AND N. PAY ATTENTION.
		//M<N is not required here since dE_TrialCrankshaft and crankshaft
		//support M>N and automatically wraps the iterator to the 
		//beginning of the chain when it hits the tail.
		//Therefore, all energy evaluation program should be careful with chain segment 
		//iteration due to wrapping problem.

        double dE, cacheRE, cacheE_t;
		int m,n;
		int E_condition=0,IEV_condition=0,topo_condition=0;
		
		//TODO: remove this after debug.
		dnaChain->checkConsistancy();
		
		if (drand(1.0)>P_REPT){
		//Crankshaft movement.
			//generate rotation axis, avoiding rigid body.
			do{
				m = irand(maxnum + 1);
				n = (m + (irand(crank_max_length - crank_min_length) + crank_min_length))%totsegnum;
			}while (protect_list[m]==1 || protect_list[n]==1);

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

			//bend energy change and movement.
			dE=dnaChain->dE_TrialCrankshaft(m, n, rotAng);
			dnaChain->crankshaft(m,n,rotAng);

			//update energies.
			RG.update_allrigid_and_E();
			this->dnaChain->E_t_updateWrithe_E_t();
			this->dnaChain->updateKPoly();
			//total energy change.
			dE= dE + (RG.E - cacheRE) + (dnaChain->E_t - cacheE_t);
			
			//Flags: 0 - uninitialized -1 - not satisfied +1 - satisfied
			E_condition=0;
			IEV_condition=0;
    		topo_condition=0;

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
			

			if (this->dnaChain->IEV(m,n)==1){
				IEV_condition=1;
			}
			else{
				IEV_condition=0;
			}

			
			if (this->dnaChain->AlexPoly[0]==1 && this->dnaChain->AlexPoly[1]==1){
				topo_condition=1;
			}
			else{
				topo_condition=0;
			}

			if (E_condition==1 && IEV_condition==1 && topo_condition==1){
					this->dnaChain->stats.rpt_accepts++;		
			}
			else{
					dnaChain->crankshaft(m,n,-rotAng);
					RG.update_allrigid_and_E();
					dnaChain->E_t_updateWrithe_E_t();
					dnaChain->updateKPoly();
			}
		}//End Crankshaft movement.
		else{
		//Reptation movement.
			//generate reptation segment.
			//The segment between vertices m and n will be changed.
			//that is vectors m~n-1
			int testp;int testflag;
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
			
			int rept_move;
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

			//bend energy change and movement.
			dE=dnaChain->dE_reptation(m,n,rept_move);

			//update energies.
			RG.update_allrigid_and_E();
			this->dnaChain->E_t_updateWrithe_E_t();
			this->dnaChain->updateKPoly();

			//total energy change.
			dE= dE + (RG.E - cacheRE) + (dnaChain->E_t - cacheE_t);
			
			E_condition=0;
			IEV_condition=0;
    		topo_condition=0;

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

			if (this->dnaChain->AlexPoly[0]==1 && this->dnaChain->AlexPoly[1]==1){
				topo_condition=1;
			}
			else{
				topo_condition=0;
			}

			if (E_condition==1 && IEV_condition==1 && topo_condition==1){
				dnaChain->stats.accepts++;
			}
			else{
					dnaChain->dE_reptation(m,n,-rept_move);
					RG.update_allrigid_and_E();
					dnaChain->E_t_updateWrithe_E_t();
					dnaChain->updateKPoly();
			}

		}//End reptation movement.

		if (moves%SNAPSHOT_INTERVAL==0){
			sprintf(buf,"%s%09d.txt",filePrefix,moves);
			dnaChain->snapshot(buf);
		}

		if (moves%STAT_INTERVAL==0){
			(*fp_log)<<"accepted:"<<dnaChain->stats.accepts()
				<<" rpt_accepted:"<<dnaChain->stats.rpt_accepts()
				<<" in moves "<<dnaChain->stats.auto_moves()
				<<'['
				<<float(dnaChain->stats.accepts())/dnaChain->stats.auto_moves()<<","
				<<float(dnaChain->stats.rpt_accepts())/dnaChain->stats.auto_moves()
				<<']'<<endl;
			dnaChain->stats.accepts.lap();
			dnaChain->stats.rpt_accepts.lap();
			dnaChain->stats.auto_moves.lap();
		
//			Log acceptance and rigid body statistics.
			(*fp_log)<<"["<<moves<<"] ";
			(*fp_log)<<"move_trial["<<m<<","<<n<<"] ";
			(*fp_log)<<"[Wr,E_t]"<<dnaChain->writhe<<","<<dnaChain->E_t;
			(*fp_log)<<"Flags(E,IEV,topo)"<<"["<<E_condition<<"(dE="<<dE<<"),"
				<<IEV_condition<<","<<topo_condition<<" Kpoly("<<dnaChain->AlexPoly[0]<<","<<dnaChain->AlexPoly[1]<<")]";
			(*fp_log)<<" r "<<RG.r<<" Ax "<<180-RG.AxisBeta/PI*180
				<<" Ra "<<180-RG.RadiusBeta/PI*180<<" E "<<RG.E<<endl;

//			Log Gyration Radius
/*			double gyration_ratio=this->calcGyration();
			dnaChain->stats.gyration_ratio.push(gyration_ratio);
			(*fp_log)<<"["<<moves<<"] ";
			(*fp_log)<<"Rg "<<gyration_ratio<<endl;
*/
//			Log Chain angle statistics
/*			for (int i=0;i<=maxnum;i++)
				dnaChain->stats.anglelist[i].push(dnaChain->C[i].bangle);
*/		}
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