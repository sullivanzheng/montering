#include "MCbox.h" 
MCbox_circular::MCbox_circular(
    char const *configFile){	
	using namespace std;
	//Read Config File;
	config=readconfig(configFile);
	
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

	//#############################MCBox Variables#####################
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

void MCbox_circular::performMetropolisCircularCrankRept(long monte_step)
{
    MTRand53 mt(seeding);
	
	allrigid RG("rigid.cfg", this->dnaChain);

	long SNAPSHOT_INTERVAL;
	std::stringstream(config["SNAPSHOT_INTERVAL"])>>SNAPSHOT_INTERVAL;

	long STAT_INTERVAL;
	std::stringstream(config["STAT_INTERVAL"])>>STAT_INTERVAL;

	for (int moves = 1; moves <= monte_step; moves++)
    {
		//MAKE MOVES

		//NOTES ON M AND N. PAY ATTENTION.
		//M<N is not required here since deltaE_TrialCrankshaft_countMove and crankshaft
		//support M>N and automatically wraps the iterator to the 
		//beginning of the chain when it hits the tail.
		//Therefore, all energy evaluation program should be careful with chain segment 
		//iteration due to wrapping problem.

        double dE, cacheRE;
		int ial[2],ierr;
		int m,n;

		if (1){//TODO (drand(1.0)>P_REPT){
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

			this->dnaChain->kpoly(ial, ierr);
			if (ierr!=0){
				cout<<"The crossing on the chain is too many!"<<endl;
				exit(EXIT_FAILURE);
			}

			if (ial[0]==1 && ial[1]==1){
				topo_condition=1;
			}
			else{
				topo_condition=0;
			}

			if (E_condition==1 && IEV_condition==1 && topo_condition==1){
					this->dnaChain->stats.accepts++;		
			}
			else{
					dnaChain->crankshaft(m,n,-rotAng);
					RG.update_allrigid_and_E();
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
				n=(m+irand(reptation_minlen,reptation_maxlen+1))%totsegnum;
				//Check if containting any rigid body segments.
				testp=m;testflag=0;
				while (testp!=n){
					if (protect_list[testp]==1){
						testflag=1;
						break;
						testp++;
						testp=testp % totsegnum;
					}
				}
			}while(testflag==1);
			
			int rept_move;
			rept_move=irand(1,4); //1~3;
	
     		//old rigid body energy.
			cacheRE=RG.E;
			//bend energy change.
			dE=dnaChain->dE_reptation(m,n,rept_move);
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

			this->dnaChain->kpoly(ial, ierr);
			if (ierr!=0){
				cout<<"The crossing on the chain is too many!"<<endl;
				exit(EXIT_FAILURE);
			}

			if (ial[0]==1 && ial[1]==1){
				topo_condition=1;
			}
			else{
				topo_condition=0;
			}

			if (E_condition==1 && IEV_condition==1 && topo_condition==1){
					this->dnaChain->stats.accepts++;		
			}
			else{
					dnaChain->dE_reptation(m,n,-rept_move);
					RG.update_allrigid_and_E();
			}
		}//End reptation movement.

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
*/			
//			Log acceptance and rigid body statistics.
			(*fp_log)<<"["<<moves<<"] ";
			(*fp_log)<<"move_trial["<<m<<","<<n<<"] ";
		    (*fp_log)<<"tp_trial("<<ial[0]<<','<<ial[1]<<")";//<<endl;
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
*/
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