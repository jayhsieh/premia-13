#ifndef MATHTOOLS_H
#define MATHTOOLS_H


#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <cstdlib>

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

using namespace std;

inline void myerror(const std::string& error_message)
{
  std::cerr << "Error: " << error_message << std::endl;
  exit(1);
}



vector<double> progonka(int N, const double low, const double diag, const double up, 
			const vector<double> & right_side);

#endif



