#ifndef READ_DATA_HPP
#define READ_DATA_HPP
#include <iostream>
#include "arma/include/armadillo"
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
using namespace std;
using namespace arma;
arma::mat read_data(char *argv,std::string sep){
  ifstream input(argv);
  string lineStr;
     string colStr;
     stringstream strStream;
     const char* colChar;
     int iRow=0;
     char _sep;
     const char* tsep;
     tsep=sep.c_str();
     if(tsep[0]=='\\'){
       _sep='\t';
     }else{
       _sep=tsep[0];
     }
     vector<string> colNames;
     //Read the first line into a string vector
     getline(input,lineStr);
     strStream.str(lineStr);
     while(getline(strStream,colStr,_sep)){
       colNames.push_back(colStr);
     }
     strStream.clear();
     vector<string> rowNames;
     vector<vector<double> > countData;
     vector<double> rowData;
     double countValue;
     while(getline(input,lineStr)){
       strStream.str(lineStr);
       getline(strStream,colStr,_sep);
       rowNames.push_back(colStr);
       while(getline(strStream,colStr,_sep)){
         colChar=colStr.c_str();
         countValue=atof(colChar);
         rowData.push_back(countValue);
       }
       countData.push_back(rowData);
       strStream.clear();
       rowData.clear();
       ++iRow;
     }
     arma::mat res(countData.size(),countData[0].size());
     int nn=countData.size();
     for(int i=0;i<nn;++i){
       vector<double> tt=countData[i];
       int mm=tt.size();
       for(int j=0;j<mm;++j){
         res(i,j)=tt[j];
       }
       tt.clear();
     }
    return res;
   }
#endif
