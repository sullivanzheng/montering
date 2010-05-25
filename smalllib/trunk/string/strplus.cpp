#include "strplus.h"
char* strcat_noOW(char *tostr, long tostr_size,char *str1, char *str2){
	if (strlen(str1)+strlen(str2)+1>tostr_size){
		std::cout<<"ERROR! String is too long for *to!"<<std::endl;
        
		getchar();
		exit(EXIT_FAILURE);
	}
	else{
		strcpy(tostr,str1);
		strcat(tostr,str2);
		return tostr;
	}
}