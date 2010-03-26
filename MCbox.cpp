#include "MCbox.h" 
MCbox_circular::MCbox_circular(
    char const *r_filePrefix, 
    int const length, 
	double const r_g,
    unsigned long r_seeding)
    :fp_log(strcat_noOW(buf, strBufSize,const_cast<char*>(r_filePrefix), "_log.txt")), 
    dnaChain(strcat_noOW(buf, strBufSize, const_cast<char*>(r_filePrefix), ".vmc"),length), 
    seeding(r_seeding)
{

	strcpy(filePrefix, r_filePrefix);
    //Initialize the global variables related to 
    //the properties of chain.
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
	dnaChain.snapshot("inisnapshot.txt");
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
        int m,n;
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
        }
		if (moves%1000==0){
			sprintf(buf,"%s%09d.txt",filePrefix,moves);
			dnaChain.snapshot(buf);
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
			<<" crank_min_length	= "<<	crank_min_length	<<endl
			<<">maxRotAng	= "<<	maxRotAng	<<endl
			<<" maxnum	= "<<	maxnum	<<endl
			<<">totsegnum * bpperseg = total bp \t"<<totsegnum<<'*'<<bpperseg<<'='<<totsegnum*bpperseg<<endl
			<<" crank_max_length	= "<<	crank_max_length	<<endl
			<<" P_SMALLROTATION	="<< P_SMALLROTATION<<endl
			<<" DELTA_TW_K =	"<< DELTA_TW_K << "//In number of turns"<<endl
			<<" Chain Lk ="<<dnaChain.Lk<<endl
			<<"================ENERGY CURVE PARAS===================="<<endl
			<<" g "<<g<<endl
			<<endl;
	}




