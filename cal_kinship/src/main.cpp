#include "kin.hpp"
#include <iostream>
#include "arma/include/armadillo"
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
using namespace std;
using namespace arma;
int main(int argc,char **argv)
{
  if(argc <=1){
    std::cerr<<"Usage: "<<argv[0]<<" [infile]"<<std::endl;
    return -1;
  }
  //char filename=argv[1];
  arma::mat B=read_data(argv[1],",");
  int par=atoi(argv[2]);
  arma::mat res_norm=cal_kin(B,par);
  arma::mat res2;
  res2.set_size(res_norm.n_rows, res_norm.n_cols+2);
  res2.submat(0, 2, res_norm.n_rows - 1, res2.n_cols - 1) = res_norm;
  res2.col(0).fill(1);
  int nn=res2.n_rows;
  for(int k=0;k<nn;k++){
      res2(k,1)=k+1;
  }
  res2.save("kinship_res.csv",csv_ascii);
  vec eigval;
  mat eigvec;
  arma::eig_sym(eigval, eigvec, res_norm);
  eigval.save("eigval.csv",csv_ascii);
  eigvec.save("eigvec.csv",csv_ascii);
  return 0;
}
