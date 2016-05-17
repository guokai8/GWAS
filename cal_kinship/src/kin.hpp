#ifndef CAL_KIN_HPP
#define CAL_KIN_HPP
#include "read_data.hpp"
#include <iostream>
#include "arma/include/armadillo"
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
using namespace std;
using namespace arma;
arma::mat cal_kin(arma::mat &B,int par)
{
  //int par=20;
  int bin_num=B.n_rows;
  double line_num=B.n_cols;
  arma::mat res(line_num,line_num);
  for(int i=0;i<bin_num;++i){
    arma::mat k(line_num,par);
    arma::rowvec a=B.row(i);
    int ss=a.size();
    for(int j=1;j<par+1;++j){
        for(int l=0;l<ss;++l){
          int x=2*(a(l)==j);
          k(l,j-1)=x;
    }
    }
    arma::mat m=k*k.t();
    res=res+m;
  }
  double diagn=trace(res);
  double norm_factor=diagn/line_num;
  arma::mat res_norm=res/norm_factor;
  return res;
}
#endif
