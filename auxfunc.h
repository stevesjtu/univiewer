#ifndef AUXFUNC_H
#define AUXFUNC_H

#include<vector>
#include<string>
#include<iostream>
#include<fstream>

#ifdef _WIN32
#include <windows.h>
#include <ShlObj.h>
#include <CommDlg.h>
#include <stdlib.h>
#endif

namespace univiewer {

void ArgParser(std::vector<std::string> &,
                        std::vector<std::string> &simple_output_result,
                        std::vector<std::string> &modelFiles,
											  std::vector<std::string> &dispFiles,
											  std::vector<std::string> &contFiles,
											  std::vector<std::string> &nodeDatafiles);


bool OpenFileDlg(std::string&, std::string&);

template<typename T>
std::ostream & operator<<(std::ostream & s, const std::vector<T> &vec) {
  if (vec.empty()) return s;
  s << "(";
  for(unsigned i=0; i < vec.size()-1; ++i){
    s << vec[i] << ", ";
  }
  s << vec[vec.size()-1] << ")";
  return s;
}

// Mdfile means markdown file, 
// which use # ## ### ... symbol to represent titles
class Mdfile {
public:
  Mdfile();
  Mdfile(const std::string &filename);
  virtual~Mdfile();

  void Open(const std::string &filename );
	void Close();
  virtual std::vector<double> GetDoubleArrayFrom(const std::string &title);
  virtual std::vector<int> GetIntArrayFrom(const std::string &title);
  virtual std::vector<unsigned int> GetUintArrayFrom(const std::string &title);
  virtual std::vector<std::string> GetStringFrom(const std::string &title); 
  void tryit();

protected:
  std::ifstream infile;
  std::streampos pos1;

  inline void JumpTo(const std::string &title) {
    infile.seekg(0, std::ios::beg);
    std::string str("");
    while(std::getline(infile, str)) { 
      if (str[0] == '#') {
        str = str.substr(1);
        trim(str);
        if (str.compare(title) == 0) {
          //pos1 = infile.tellg();
          break;
        }
      }
    } 
  }

  inline void trim(std::string &s) {
    if (s.empty()) return;
    s.erase(0,s.find_first_not_of(" "));
    s.erase(s.find_last_not_of(" ") + 1);
  }

};

}

#endif