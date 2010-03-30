#include "configldr.h"
using namespace std;

dict readconfig(char const *filename)
{
	ifstream fp(filename);
	if (!fp){
		cout<<"could not open the config file!"<<endl;
		getchar();exit(EXIT_FAILURE);
	}
	string bufK,bufV;
	map<string,string> config;
	while (getline(fp,bufK)){
        if (bufK.length()==0) {
            std::cout<<"-------"<<std::endl;
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
	return config;
}