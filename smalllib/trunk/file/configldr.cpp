#include "configldr.h"
using namespace std;

dict readconfig(char const *filename)
{
	std::cout<<"[readconfig] is called"<<std::endl;
	ifstream fp(filename);
	if (!fp){
		cout<<"could not open the config file!"<<endl;
		getchar();exit(EXIT_FAILURE);
	}
	string bufK,bufV;
	map<string,string> config;
	while (getline(fp,bufK)){
        if (bufK.length()==0) {
            std::cout<<"[readconifg]-------blankline in config file encountered."<<std::endl;
            continue;
        }
        stringstream linebuf(bufK);
        linebuf << bufK;
		linebuf >> bufK >> bufV;
        if (bufK[0]=='#'){
        }
        else{
            std::cout <<bufK<<" : "<<bufV<<std::endl;
		    config[bufK]=string(bufV);
        }
	}
	fp.close();
	std::cout<<"[readconfig]configs are loaded successfully."<<std::endl;
	return config;
}