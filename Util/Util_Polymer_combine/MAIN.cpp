#include "nr.h"
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
using namespace std;


const int aC = 0;
const int aSi = 1;
const int aO = 2;
const int aH = 3;
const int aZn = 4;
const int aMg = 5;
const int aN = 6;
const int aTi = 7;
const int aI = 8;
const int aPb = 9;
const int aLa = 10;
const int aNa = 11;
const int aP = 12;
const int aF = 13;

const int numC = 3;

int namax;
int na;

int main() {

        std::sprintf;
	        std::ofstream outfile;
		std::ofstream outfile2;
        char command0[4500]; 
  char command1[4500]; 
  char command2[4500]; 
  char command3[4500]; 
  char command4[4500]; 
  char command5[4500];
    namax = 500000;
        string data0;
        string data1;
	string string1, string2, system1, system2;
        DP N_node, rr;
        DP v, sin_a, sin_b, sin_c, cos_a, cos_b, cos_c;
    DP avg_x0, avg_x1, avg_x, avg_y0, avg_y1, avg_y, avg_z0, avg_z1, avg_z;
    DP sumx, sumy, sumz, del_x0, del_x1, del_y0, del_y1, del_z0, del_z1;
    DP len0, len1, dd, dd2;
    Mat_IO_DP x(namax, numC);
    Vec_IO_INT elem(namax), b_id0(namax), b_id1(namax), ad0(namax), ad1(namax), h0(namax), h1(namax);
    
    string C("C");
    string Si("Si");
    string O("O");
    string H("H");
    string Zn("Zn");
    string Mg("Mg");
    string N("N");
    string Ti("Ti");
    string I("I");
    string Pb("Pb");
    string La("La");
    string Na("Na");
    string F("F");
    string P("P");
    string str1 = ".";
    string str2 = "+";
    string str3 = "-";
    string str4 = "C1";
    char* ptr;
    
     int numC = 3;
     int na, n_N, ii, jj, iii, jjj, nb, nnn, kkk, ik, naa, ikik;
    int logic1;
    DP penalty, penalty2;
    
    elem = -1;
     DP pi;
      Vec_IO_DP len_c(numC), hlen(numC), cen(numC), angle(numC), comm(namax);
    pi = 3.141592;
    //    int system;
   //     int m,mm,mmm,mmm,mmmm;
   
    // Empty vector holding all names from file
         vector<string> names;
         
    //         // Read names from file LineUp.txt
                 ifstream in("name.txt");
                     if(!in.is_open())
                             cout << "Unable to open file\n";
    
    //                             // this is wrong, by the way: while(in.good()){
                                     string word;
                                        while(getline(in, word))
                                                     names.push_back(word);
    //
     //                                                    sort(names.begin(), names.end());
    //
    //                                                         // Loop to print names

  //  cout << com[5] << endl;
    
    
    ikik = 0;
    
    ofstream outputFile10077("name_list.txt", ios::out);
    
  for (size_t i = 0; i < names.size(); i++) {
    
      penalty = 0;
      penalty2 = 0;

      
      n_N = 0;
      
      ofstream outputFile55("test_back.smi", ios::out);
      outputFile55 << names[i].c_str() << endl;
      
      sprintf (command1, "cp test_back.smi test.smi");
      int systemRet1000 =  system(command1);
      
      
      sprintf (command1, "taskset --cpu-list 0 obabel -ismi test.smi -omol2 --gen3d --partialcharge eem -O monomer.mol2");
      int systemRet2000 =  system(command1);
      
      sprintf (command1,"../Util/Util_monomer_matching/exe");
      int systemRet2001 =  system(command1);
 
      sprintf (command1, "taskset --cpu-list 0 obabel -imol2 monomer_reorder.mol2 -osmi -O test2.smi");
      int systemRet2002 =  system(command1);
      
      ifstream inputFile505("test2.smi", ios::in);
      inputFile505 >> system1 >> system2;
      
      string result_smiles = system1;
      string diamine_substitute = "(NC8=CC=C(N)C=C8)";
      
      // Find and replace the single-bonded O in -COOH patterns
      bool cooh_replaced = false;
      
      // Pattern 1: C(=O)O -> C(=O)[diamine]
      size_t pos = result_smiles.find("C(=O)O");
      if (pos != string::npos) {
          result_smiles.replace(pos + 5, 1, diamine_substitute); // Replace the O at position +5
          cooh_replaced = true;
      }
      
      // Pattern 2: C(O)=O -> C([diamine])=O
      if (!cooh_replaced) {
          pos = result_smiles.find("C(O)=O");
          if (pos != string::npos) {
              result_smiles.replace(pos + 2, 1, diamine_substitute); // Replace the O at position +2
              cooh_replaced = true;
          }
      }
      
      // Pattern 3: OC(=O) -> [diamine]C(=O)
      if (!cooh_replaced) {
          pos = result_smiles.find("OC(=O)");
          if (pos != string::npos) {
              result_smiles.replace(pos, 1, diamine_substitute); // Replace the O at the beginning
              cooh_replaced = true;
          }
      }
      
      // Pattern 4: O=C(O) -> O=C([diamine])
      if (!cooh_replaced) {
          pos = result_smiles.find("O=C(O)");
          if (pos != string::npos) {
              result_smiles.replace(pos + 4, 1, diamine_substitute); // Replace the O at position +4
              cooh_replaced = true;
          }
      }
      
      // Pattern 5: C(=O)(O) -> C(=O)([diamine])
      if (!cooh_replaced) {
          pos = result_smiles.find("C(=O)(O)");
          if (pos != string::npos) {
              result_smiles.replace(pos + 6, 1, diamine_substitute); // Replace the O at position +6
              cooh_replaced = true;
          }
      }

      // Replace remaining single-bonded O in -COOH patterns with Cl
      bool cooh_cl_replaced = false;
      
      // Pattern 1: C(=O)O -> C(=O)Cl
      pos = result_smiles.find("C(=O)O");
      if (pos != string::npos) {
          result_smiles.replace(pos + 5, 1, "Cl"); // Replace the O at position +5
          cooh_cl_replaced = true;
      }
      
      // Pattern 2: C(O)=O -> C(Cl)=O
      if (!cooh_cl_replaced) {
          pos = result_smiles.find("C(O)=O");
          if (pos != string::npos) {
              result_smiles.replace(pos + 2, 1, "Cl"); // Replace the O at position +2
              cooh_cl_replaced = true;
          }
      }
      
      // Pattern 3: OC(=O) -> ClC(=O)
      if (!cooh_cl_replaced) {
          pos = result_smiles.find("OC(=O)");
          if (pos != string::npos) {
              result_smiles.replace(pos, 1, "Cl"); // Replace the O at the beginning
              cooh_cl_replaced = true;
          }
      }
      
      // Pattern 4: O=C(O) -> O=C(Cl)
      if (!cooh_cl_replaced) {
          pos = result_smiles.find("O=C(O)");
          if (pos != string::npos) {
              result_smiles.replace(pos + 4, 1, "Cl"); // Replace the O at position +4
              cooh_cl_replaced = true;
          }
      }

      // Pattern 5: C(=O)(O) -> C(=O)(Cl)
      if (!cooh_cl_replaced) {
          pos = result_smiles.find("C(=O)(O)");
          if (pos != string::npos) {
              result_smiles.replace(pos + 6, 1, "Cl"); // Replace the O at position +6
              cooh_cl_replaced = true;
          }
      }
      
      ofstream outputFile555("monomer.smi", ios::out);
      outputFile555 << result_smiles << endl;
         
      sprintf (command1, "cp monomer.smi structures/monomer_%d.smi",ikik);
      int systemRet1005 =  system(command1);
      
      outputFile10077 << "monomer_"<< ikik << endl;
      
      ikik++;
      
}

 

  return 0;


                                                         
}
                                       
