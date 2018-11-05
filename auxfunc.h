#ifndef AUXFUNC_H
#define AUXFUNC_H

#include<vector>
#include<string>
#include<iostream>

void argParser(const int& argc, char* argv[], 
                        std::vector<std::string> &simple_output_result,
                        std::vector<std::string> &modelFiles,
											  std::vector<std::string> &dispFiles,
											  std::vector<std::string> &contFiles,
											  std::vector<std::string> &nodeDatafiles);


#endif