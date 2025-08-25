#include "nr.h"
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
using namespace std;
using std::cout;

int main() {

        std::ofstream outfile;
  char command0[4500];
  char command1[4500]; 
  char command2[4500]; 
  char command3[4500]; 
  char command4[4500]; 
  char command5[4500]; 
        string data0;
        string data1;
        string string1, string2;
        int N_node, m, n, N_bin, ii, jj, N_MPI, na, j, kk, num_b, num_d, kkk, na1, na2, na3, na4, na5;
    Mat_IO_INT check1(5000, 5000), check2(5000, 5000);
    Mat_IO_DP dis1(10000, 2), dis2(10000, 2), dis3(10000,2), dis4(10000,2), dis5(10000,2);
        Mat_IO_DP dis(5000, 5000 ), count1(5000,5000), count2(5000,5000);
        Vec_IO_DP PE(300000);
        DP rr;
    //    int system;
   //     int m,mm,mmm,mmm,mmmm;
  vector<string> names;
  vector<string> names2;
  N_bin = 10000;   



  
  ifstream inputFile3("output1.txt",ios::in);
  if (!inputFile3) {
    cerr << "inputFile1 could not be opened" << endl;
    exit(1);
  }

    
  ifstream inputFile33("output2.txt",ios::in);
  if (!inputFile33) {
    cerr << "inputFile2 could not be opened" << endl;
    exit(1);
  }


  kkk = 0;

    inputFile3 >> string1 >> string1 >> string1 >> string1 >> string1 >> string1;
    inputFile3 >> string1 >> string1 >> string1;
    
  while (!inputFile3.eof())
  {
    inputFile3 >> dis2[kkk][0] >> dis2[kkk][1];
    kkk = kkk + 1;
  }

   na2 = kkk-1;


  kkk = 0;
    
    inputFile33 >> string1 >> string1 >> string1 >> string1 >> string1 >> string1;
    inputFile33 >> string1 >> string1 >> string1;
    
  while (!inputFile33.eof())
  {
    inputFile33 >> dis1[kkk][0] >> dis1[kkk][1];
    kkk = kkk + 1;
  }

  na1 = kkk-1;
    
  ifstream inputFile331("output3.txt",ios::in);
  if (!inputFile331) {
    cerr << "inputFile3 could not be opened" << endl;
    exit(1);
  }


  kkk = 0;

    inputFile331 >> string1 >> string1 >> string1 >> string1 >> string1 >> string1;
    inputFile331 >> string1 >> string1 >> string1;

  while (!inputFile331.eof())
  {
    inputFile331 >> dis3[kkk][0] >> dis3[kkk][1];
    kkk = kkk + 1;
  }

  na3 = kkk-1;


  ifstream inputFile332("output4.txt",ios::in);
  if (!inputFile332) {
    cerr << "inputFile4 could not be opened" << endl;
    exit(1);
  }


  kkk = 0;

    inputFile332 >> string1 >> string1 >> string1 >> string1 >> string1 >> string1;
    inputFile332 >> string1 >> string1 >> string1;

  while (!inputFile332.eof())
  {
    inputFile332 >> dis4[kkk][0] >> dis4[kkk][1];
    kkk = kkk + 1;
  }

  na4 = kkk-1;

   // cout << na2 << "   " << na1 << endl;

  ofstream outputFile3222("output_all.txt",ios::out);

      outputFile3222 << dis2[na2-1][1] << "    " << dis1[na1-1][1] <<"    " << dis3[na3-1][1] <<"    " << dis4[na4-1][1] << "    " <<  endl;
                 
}
                                       
