#ifndef CONFIGLDR_H
#define CONFIGLDR_H
#include <iostream> 
#include <fstream>
#include <sstream>
#include <string>
#include <map>

typedef std::map <std::string, std::string> dict;

dict readconfig(const char* filename);

#endif /* CONFIGLDR_H */