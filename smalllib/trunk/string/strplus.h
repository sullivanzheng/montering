
#include <iostream>
#include <cstring>
#include <string>
#include <locale>


#ifndef STRPLUS_H
#define STRPLUS_H


char* strcat_noOW(char *tostr, int tostr_size,char *str1, char *str2);

/*
std::string stoupper(std::string s)
{ 
  std::locale loc;
  std::string::iterator i = s.begin();
  std::string::iterator end = s.end();

  while (i != end) {
	  *i = std::toupper((unsigned char)*i,loc);
    ++i;
  }
  return s;
}

std::string stolower(std::string s)
{
  std::locale loc;
  std::string::iterator i = s.begin();
  std::string::iterator end = s.end();

  while (i != end) {
	  *i = std::tolower((unsigned char)*i,loc);
    ++i;
  }
  return s;
}
*/

#endif /* STRPLUS_H */